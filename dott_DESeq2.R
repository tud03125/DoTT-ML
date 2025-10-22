#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(DESeq2)
  # Optional helpers (auto-used if present; no required flags)
  ihw_ok    <- requireNamespace("IHW", quietly = TRUE)
  apeglm_ok <- requireNamespace("apeglm", quietly = TRUE)
  ev_ok     <- requireNamespace("EnhancedVolcano", quietly = TRUE)
  gg_ok     <- requireNamespace("ggplot2", quietly = TRUE)
  edgeR_ok  <- requireNamespace("edgeR", quietly = TRUE)
})

# -----------------------------
# Args (unchanged)
# -----------------------------
parser <- ArgumentParser()
parser$add_argument("--counts_file", required=TRUE, help="Cleaned counts file")
parser$add_argument("--output_dir", required=TRUE, help="Output directory")
parser$add_argument("--conditions", required=TRUE, help="Comma-separated list of condition labels")
parser$add_argument("--bootstrap", type="logical", default=FALSE, help="Whether to perform bootstrapping")
parser$add_argument("--n_boot", type="integer", default=100, help="Number of bootstrap iterations")
parser$add_argument("--consensus_threshold", type="double", default=0.5, help="Fraction threshold for consensus calls")
args <- parser$parse_args()

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Tunables
# -----------------------------
alpha   <- 0.05
lfc_thr <- 1.0
set.seed(1L)
prefix <- "3_UTR_extended"

# Profile selection: "robust" (default) or "classic"
profile <- tolower(Sys.getenv("DOTT_DESEQ2_PROFILE", "robust"))
if (!profile %in% c("robust","classic")) {
  warning("Unknown DOTT_DESEQ2_PROFILE='", profile, "'. Falling back to 'robust'.")
  profile <- "robust"
}
message("DESeq2 profile: ", profile)

# -----------------------------
# Load counts
# -----------------------------
infer_sep <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv","txt")) "\t" else if (ext == "csv") "," else "\t"
}
sep <- infer_sep(args$counts_file)
tab <- fread(args$counts_file, sep=sep, data.table=FALSE, check.names=FALSE)

if ("Geneid" %in% colnames(tab)) {
  gene_ids <- tab$Geneid
  mat <- as.matrix(tab[, setdiff(colnames(tab), "Geneid"), drop=FALSE])
  rownames(mat) <- gene_ids
} else {
  gene_ids <- tab[[1]]
  mat <- as.matrix(tab[, -1, drop=FALSE])
  rownames(mat) <- gene_ids
}
storage.mode(mat) <- "integer"

# unique colnames for safety
if (any(duplicated(colnames(mat)))) {
  colnames(mat) <- make.unique(colnames(mat), sep = ".dup")
}

# -----------------------------
# Conditions -> design
# control = first unique label, experimental = second (matches your runs)
# -----------------------------
conds <- strsplit(args$conditions, ",")[[1]]
stopifnot(length(conds) == ncol(mat))
u <- unique(conds); stopifnot(length(u) == 2)
control <- u[1]; experimental <- u[2]

coldata <- data.frame(
  sample = colnames(mat),
  group  = factor(conds, levels = c(control, experimental)),
  row.names = colnames(mat)
)

dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata, design = ~ group)

# -----------------------------
# Prefiltering
# -----------------------------
if (profile == "robust" && edgeR_ok) {
  y <- edgeR::DGEList(counts = counts(dds), group = coldata$group)
  keep <- edgeR::filterByExpr(y, group = coldata$group)
  if (sum(keep) < 200) keep <- edgeR::filterByExpr(y, group = coldata$group, min.count = 1)
} else {
  keep <- rep(TRUE, nrow(dds))
}
dds <- dds[keep, , drop=FALSE]

# -----------------------------
# Fit (with fallback)
# -----------------------------
fit_deseq_with_fallback <- function(dds, quiet=TRUE) {
  fit1 <- try(DESeq(dds, fitType = "parametric", quiet = quiet), silent = TRUE)
  if (!inherits(fit1, "try-error")) return(fit1)
  message("Parametric dispersion fit failed; retrying fitType='local' ...")
  fit2 <- try(DESeq(dds, fitType = "local", quiet = quiet), silent = TRUE)
  if (!inherits(fit2, "try-error")) return(fit2)
  message("Local fit also failed; using gene-wise dispersion estimates + Wald test ...")
  dds2 <- estimateSizeFactors(dds)
  dds2 <- estimateDispersionsGeneEst(dds2, quiet = quiet)
  dispersions(dds2) <- mcols(dds2)$dispGeneEst
  dds2 <- nbinomWaldTest(dds2, betaPrior = FALSE)
  dds2
}

if (profile == "robust") {
  dds <- fit_deseq_with_fallback(dds, quiet = TRUE)
} else {
  fit1 <- try(DESeq(dds, fitType = "parametric", quiet = TRUE), silent = TRUE)
  if (inherits(fit1, "try-error")) {
    message("Classic profile: parametric fit failed; trying 'local'.")
    dds <- DESeq(dds, fitType = "local", quiet = TRUE)
  } else dds <- fit1
}

# -----------------------------
# Results
# -----------------------------
if (profile == "robust") {
  res <- results(
    dds, contrast = c("group", experimental, control),
    alpha = alpha, lfcThreshold = lfc_thr, altHypothesis = "greaterAbs",
    independentFiltering = TRUE
  )
  if (ihw_ok) {
    library(IHW)
    df <- as.data.frame(res)
    df$baseMean[is.na(df$baseMean)] <- 0
    fit <- ihw(pvalue ~ baseMean, data = df, alpha = alpha)
    res$padj <- adj_pvalues(fit)
  }
} else {
  res <- results(
    dds, contrast = c("group", experimental, control),
    alpha = alpha, independentFiltering = TRUE
  )
}

# optional LFC shrinkage for plots
shrunk <- NULL
if (apeglm_ok) {
  library(apeglm)
  rn <- resultsNames(dds)
  coef_guess <- paste0("group_", make.names(experimental), "_vs_", make.names(control))
  coef_name <- if (coef_guess %in% rn) coef_guess else rn[grep("^group_.*_vs_.*", rn)][1]
  if (!is.na(coef_name)) {
    suppressWarnings({
      shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    })
  }
}

# -----------------------------
# Write normalized counts + results (legacy filenames)
# -----------------------------
norm_counts <- counts(dds, normalized = TRUE)
fwrite(data.table(gene = rownames(norm_counts), as.data.frame(norm_counts)),
       file.path(args$output_dir, "normalized_counts.csv"), sep = ",")

idx_ctrl <- which(coldata$group == control)
idx_exp  <- which(coldata$group == experimental)
mean_ctrl <- rowMeans(norm_counts[, idx_ctrl, drop=FALSE])
mean_exp  <- rowMeans(norm_counts[, idx_exp,  drop=FALSE])

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$mean_control <- mean_ctrl[res_df$gene]
res_df$mean_experimental <- mean_exp[res_df$gene]

res_out <- res_df[, c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj",
                      "mean_control","mean_experimental")]
fwrite(res_out,
       file.path(args$output_dir, paste0(prefix, "_differential_analysis_results.csv")),
       sep = ",")

sig_abs <- subset(res_out, !is.na(padj) & padj <= alpha & abs(log2FoldChange) >= lfc_thr)
fwrite(sig_abs, file.path(args$output_dir, "significant_extended_genes.csv"), sep = ",")
fwrite(sig_abs, file.path(args$output_dir, "absolute_significant_extended_genes_with_individual_means.csv"), sep = ",")

# -----------------------------
# Plots (MA & Volcano)
# -----------------------------
thr_p  <- alpha
thr_fc <- lfc_thr

res_df$direction <- "NS"
sel <- !is.na(res_df$padj) & (res_df$padj < thr_p) & (abs(res_df$log2FoldChange) > thr_fc)
res_df$direction[ sel & res_df$log2FoldChange >  thr_fc] <- paste0("Up in ", experimental)
res_df$direction[ sel & res_df$log2FoldChange < -thr_fc] <- paste0("Up in ", control)
res_df$direction <- factor(
  res_df$direction,
  levels = c("NS", paste0("Up in ", experimental), paste0("Up in ", control))
)

if (gg_ok) {
  library(ggplot2)
  p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = direction)) +
    geom_point(alpha = 0.55, size = 0.8, na.rm = TRUE) +
    geom_hline(yintercept = c(-thr_fc, thr_fc), linetype = "dashed") +
    labs(
      title = sprintf("%s MA plot (%s vs %s)", prefix, experimental, control),
      subtitle = "Significant: padj < 0.05 & |LFC| > 1",
      x = "log10(baseMean + 1)",
      y = sprintf("log2 fold change (%s / %s)", experimental, control),
      color = NULL
    ) +
    theme_bw() + theme(legend.position = "top")
  ggsave(file.path(args$output_dir, paste0(prefix, "_MA_plot.svg")), p_ma,
         width = 14, height = 10, units = "in")
  ggsave(file.path(args$output_dir, paste0(prefix, "_MA_plot.png")), p_ma,
         width = 14, height = 10, units = "in", dpi = 300)
} else {
  cols <- setNames(c("grey60","firebrick","steelblue"),
                   c("NS", paste0("Up in ", experimental), paste0("Up in ", control)))
  svg(file.path(args$output_dir, paste0(prefix, "_MA_plot.svg")), width=14, height=10)
  with(res_df, {
    colv <- cols[direction]
    plot(log10(baseMean + 1), log2FoldChange, pch=20, cex=0.6, col=colv,
         xlab="log10(baseMean+1)",
         ylab=sprintf("log2 fold change (%s / %s)", experimental, control),
         main=sprintf("%s MA plot (%s vs %s)", prefix, experimental, control))
    abline(h=c(-thr_fc, thr_fc), lty=2)
    legend("topright", legend=names(cols), col=cols, pch=16, cex=0.8, bty="n")
  })
  dev.off()
}

vol_df <- if (!is.null(shrunk)) as.data.frame(shrunk) else res_df
vol_df$gene <- rownames(vol_df)

if (ev_ok) {
  library(EnhancedVolcano)
  vol <- EnhancedVolcano(
    vol_df,
    lab = vol_df$gene,
    x   = "log2FoldChange",
    y   = "padj",
    pCutoff  = thr_p,
    FCcutoff = thr_fc,
    title    = sprintf("%s Volcano (%s vs %s)", prefix, experimental, control),
    subtitle = "Significant: padj < 0.05 & |LFC| > 1",
    cutoffLineType = "dashed",
    labSize  = 2.6,
    drawConnectors = FALSE
  )
  if (gg_ok) {
    ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), vol,
           width = 14, height = 10, units = "in")
    ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.png")), vol,
           width = 14, height = 10, units = "in", dpi = 300)
  } else {
    svg(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), width=14, height=10)
    print(vol)
    dev.off()
  }
} else if (gg_ok) {
  library(ggplot2)
  vol_df$padj[!is.na(vol_df$padj) & vol_df$padj == 0] <- .Machine$double.xmin
  vol_df$neglog10padj <- -log10(vol_df$padj)
  p_vol <- ggplot(vol_df, aes(x = log2FoldChange, y = neglog10padj, color = res_df$direction)) +
    geom_point(alpha = 0.55, size = 0.8, na.rm = TRUE) +
    geom_vline(xintercept = c(-thr_fc, thr_fc), linetype = "dashed") +
    geom_hline(yintercept = -log10(thr_p), linetype = "dashed") +
    labs(
      title = sprintf("%s Volcano (%s vs %s)", prefix, experimental, control),
      x = sprintf("log2 fold change (%s / %s)", experimental, control),
      y = "-log10(padj)", color = NULL
    ) +
    theme_bw() + theme(legend.position = "top")
  ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), p_vol,
         width = 14, height = 10, units = "in")
  ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.png")), p_vol,
         width = 14, height = 10, units = "in", dpi = 300)
} else {
  cols <- setNames(c("grey60","firebrick","steelblue"),
                   c("NS", paste0("Up in ", experimental), paste0("Up in ", control)))
  ypadj <- res_df$padj
  ypadj[!is.na(ypadj) & ypadj == 0] <- .Machine$double.xmin
  svg(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), width=14, height=10)
  with(res_df, {
    colv <- cols[direction]
    plot(log2FoldChange, -log10(ypadj), pch=20, cex=0.6, col=colv,
         xlab=sprintf("log2 fold change (%s / %s)", experimental, control),
         ylab="-log10(padj)",
         main=sprintf("%s Volcano (%s vs %s)", prefix, experimental, control))
    abline(v=c(-thr_fc, thr_fc), lty=2); abline(h=-log10(thr_p), lty=2)
    legend("topright", legend=names(cols), col=cols, pch=16, cex=0.8, bty="n")
  })
  dev.off()
}

# -----------------------------
# Violin (points-only) with auto batch inference for PLOTTING
# -----------------------------
outdir <- args$output_dir
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggbeeswarm)
})
have_limma   <- requireNamespace("limma",   quietly = TRUE)
have_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
have_svglite <- requireNamespace("svglite", quietly = TRUE)

clean_names <- function(x) {
  x <- basename(x)
  x <- stringr::str_remove(x, "_sorted\\.bam$|\\.markdup\\.sorted\\.bam$|\\.bam$")
  x
}

# ---- NEW: robust, generic batch inference from sample names ----
infer_batch_generic <- function(sample_ids, min_per_level = 2) {
  # 1) Extract a leading alphanumeric token (e.g., IA, I, RNAseq, IP, Batch1)
  token1 <- sub("^([A-Za-z0-9]+).*", "\\1", sample_ids)
  tab1 <- table(token1)
  ok1 <- names(tab1)[tab1 >= min_per_level]

  if (length(ok1) >= 2) {
    b <- ifelse(token1 %in% ok1, token1, NA_character_)
    return(b)
  }

  # 2) If not informative, split on common separators and use first field
  split_first <- function(s) {
    parts <- unlist(strsplit(s, "[._:-]+"))
    if (length(parts) == 0) return(NA_character_)
    parts[1]
  }
  token2 <- vapply(sample_ids, split_first, character(1))
  tab2 <- table(token2)
  ok2 <- names(tab2)[tab2 >= min_per_level]
  if (length(ok2) >= 2) {
    b <- ifelse(token2 %in% ok2, token2, NA_character_)
    return(b)
  }

  # 3) Otherwise, return NA (no batch detected)
  rep(NA_character_, length(sample_ids))
}

# Save coldata for debugging/repro
readr::write_tsv(
  data.frame(sample = rownames(coldata), coldata),
  file.path(outdir, "sample_info.tsv")
)

plot_dott_violin <- function(outdir,
                             counts_csv = file.path(outdir, "normalized_counts.csv"),
                             results_csv = file.path(outdir, "3_UTR_extended_differential_analysis_results.csv"),
                             colData,
                             experimental_condition,    # string, e.g. "HCD"
                             sig_only = TRUE, padj_thr = 0.05, lfc_thr = 1.0,
                             batch_correct = TRUE,
                             svg_file = NULL,
                             width = 14, height = 10) {

  stopifnot(file.exists(counts_csv))
  cts <- readr::read_csv(counts_csv, show_col_types = FALSE)
  names(cts)[1] <- "gene"

  res <- NULL
  if (file.exists(results_csv)) {
    res <- readr::read_csv(results_csv, show_col_types = FALSE)
    names(res)[1] <- "gene"
  }

  # Harmonize IDs
  counts_samples <- clean_names(names(cts)[-1])
  names(cts) <- c("gene", counts_samples)

  if (!("sample" %in% colnames(colData))) colData$sample <- rownames(colData)
  colData$sample_id_clean <- clean_names(colData$sample)
  rownames(colData) <- make.unique(colData$sample_id_clean, sep = ".dup")

  # map group -> condition if needed
  if (!("condition" %in% colnames(colData)) && ("group" %in% colnames(colData))) {
    colData$condition <- colData$group
  }

  # ---- auto batch inference IF user didn't supply batch ----
  if (!("batch" %in% colnames(colData))) {
    inferred <- infer_batch_generic(rownames(colData))
    if (any(!is.na(inferred))) {
      # only keep if at least 2 levels AND each has >=2 samples
      lev <- table(inferred)
      keep_levels <- names(lev)[lev >= 2]
      if (length(keep_levels) >= 2) {
        inferred[!inferred %in% keep_levels] <- NA_character_
        colData$batch <- inferred
        message("Inferred batch levels: ", paste(keep_levels, collapse=", "))
      }
    }
  }

  # Keep overlap
  keep <- intersect(rownames(colData), counts_samples)
  if (length(keep) < 2) stop("No overlapping samples between counts and colData.")
  cts     <- cts[, c("gene", keep), drop = FALSE]
  colData <- colData[keep, , drop = FALSE]

  # optional significance filter
  title_text <- "DoTT expression across conditions"
  cts_plot <- cts
  if (isTRUE(sig_only) && !is.null(res) && all(c("padj","log2FoldChange") %in% names(res))) {
    sig_res <- dplyr::filter(res, !is.na(padj), padj < padj_thr, abs(log2FoldChange) > lfc_thr)
    if (nrow(sig_res) > 0) {
      cts_plot <- dplyr::semi_join(cts, dplyr::select(sig_res, gene), by = "gene")
      title_text <- paste0(title_text, " (significant genes only)")
    } else {
      title_text <- paste0(title_text, " (all DoTT regions; no padj<", padj_thr, " & |L2FC|>", lfc_thr, ")")
    }
  }

  # build log matrix
  mat  <- as.matrix(cts_plot[, -1, drop = FALSE])
  rownames(mat) <- cts_plot$gene
  logm <- log2(mat + 1)

  # condition (control first, then experimental)
  cond <- factor(colData$condition)
  if (!(experimental_condition %in% levels(cond))) {
    stop("experimental_condition '", experimental_condition, "' is not in colData$condition")
  }
  cond <- stats::relevel(cond, ref = setdiff(levels(cond), experimental_condition)[1])
  colData$condition <- cond

  # optional batch correction for visualization (preserve condition effect in design)
  use_log <- logm
  if (isTRUE(batch_correct) && "batch" %in% colnames(colData) && have_limma) {
    if (length(na.omit(unique(colData$batch))) >= 2) {
      design <- model.matrix(~ condition, data = colData)
      use_log <- limma::removeBatchEffect(logm, batch = colData$batch, design = design) # visualization only
      title_text <- paste0(title_text, " (batch-corrected)")
    }
  }

  # long format
  long_df <- as.data.frame(use_log) |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expr") |>
    dplyr::mutate(condition = colData[sample, "condition", drop = TRUE]) |>
    dplyr::mutate(condition = factor(condition, levels = levels(colData$condition)))

  # palette
  ctrl <- levels(long_df$condition)[1]; exp <- levels(long_df$condition)[2]
  pal  <- c(setNames("#2CA9E1", ctrl), setNames("#E76F51", exp))

  # plot
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = condition, y = expr, colour = condition)) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(x = NULL, y = "log2(normalized DoTT counts + 1)",
                  title = title_text,
                  subtitle = paste0(ctrl, " vs ", exp)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   plot.title.position = "plot",
                   plot.margin = ggplot2::margin(t = 22, r = 12, b = 10, l = 10))

  pts <- ggbeeswarm::geom_quasirandom(width = 0.18, alpha = 0.65, size = 0.55, show.legend = FALSE)
  if (have_ggrastr) p <- p + ggrastr::rasterise(pts, dpi = 300) else p <- p + pts

  med <- long_df |>
    dplyr::group_by(condition) |>
    dplyr::summarise(y = stats::median(expr), .groups = "drop") |>
    dplyr::mutate(x = as.numeric(condition) - 0.35, xend = as.numeric(condition) + 0.35)

  p <- p +
    ggplot2::geom_segment(data = med,
                          ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                          inherit.aes = FALSE, linewidth = 1.05, color = "black", lineend = "butt") +
    ggplot2::geom_point(data = med,
                        ggplot2::aes(x = as.numeric(condition), y = y),
                        inherit.aes = FALSE, shape = 21, size = 3.1, fill = "white",
                        color = "black", stroke = 1)

  wt  <- stats::wilcox.test(expr ~ condition, data = long_df, exact = FALSE)
  lab <- sprintf("Wilcoxon, p = %.15g", as.numeric(wt$p.value))
  p   <- p + ggplot2::annotate("text", x = 1.5, y = max(long_df$expr)*0.98, label = lab, size = 4.1)

  if (is.null(svg_file)) svg_file <- file.path(outdir, paste0("DoTT_violin_", exp, "_", ctrl, ".svg"))
  if (!have_svglite) warning("svglite not found; ggsave() will still write SVG via built-in device.")
  ggplot2::ggsave(svg_file, p, device = "svg", width = width, height = height, dpi = 300)
  message("Saved: ", svg_file)
}

# Call it
try({
  plot_dott_violin(outdir                 = outdir,
                   counts_csv             = file.path(outdir, "normalized_counts.csv"),
                   results_csv            = file.path(outdir, "3_UTR_extended_differential_analysis_results.csv"),
                   colData                = coldata,
                   experimental_condition = experimental,
                   sig_only               = TRUE,
                   padj_thr               = 0.05,
                   lfc_thr                = 1.0,
                   batch_correct          = TRUE)
}, silent = FALSE)

# -----------------------------
# Bootstrap consensus (unchanged)
# -----------------------------
if (isTRUE(args$bootstrap)) {
  message(sprintf("Bootstrapping %d iterations for consensus >= %.2f",
                  args$n_boot, args$consensus_threshold))
  B <- args$n_boot
  genes <- rownames(mat)
  call_mat <- matrix(0L, nrow = nrow(mat), ncol = B,
                     dimnames = list(genes, paste0("b", seq_len(B))))
  idx_ctrl <- which(coldata$group == control)
  idx_exp  <- which(coldata$group == experimental)

  for (b in seq_len(B)) {
    rs_ctrl <- sample(idx_ctrl, length(idx_ctrl), replace = TRUE)
    rs_exp  <- sample(idx_exp,  length(idx_exp),  replace = TRUE)
    m_b <- cbind(mat[, rs_ctrl, drop=FALSE], mat[, rs_exp, drop=FALSE])
    colnames(m_b) <- make.unique(c(colnames(mat)[rs_ctrl], colnames(mat)[rs_exp]), sep = ".rep")
    cd_b <- data.frame(group = factor(c(rep(control, length(rs_ctrl)),
                                        rep(experimental, length(rs_exp))),
                                      levels = c(control, experimental)),
                       row.names = colnames(m_b), stringsAsFactors = FALSE)
    dds_b <- DESeqDataSetFromMatrix(countData = m_b, colData = cd_b, design = ~ group)

    if (profile == "robust" && edgeR_ok) {
      yb <- edgeR::DGEList(counts = counts(dds_b), group = cd_b$group)
      keep_b <- edgeR::filterByExpr(yb, group = cd_b$group)
      if (sum(keep_b) < 200) keep_b <- edgeR::filterByExpr(yb, group = cd_b$group, min.count = 1)
    } else {
      keep_b <- rep(TRUE, nrow(dds_b))
    }

    dds_b <- dds_b[keep_b, , drop=FALSE]
    if (profile == "robust") {
      dds_b <- fit_deseq_with_fallback(dds_b, quiet = TRUE)
      res_b <- results(dds_b,
                       contrast = c("group", experimental, control),
                       alpha = alpha,
                       lfcThreshold = lfc_thr,
                       altHypothesis = "greaterAbs",
                       independentFiltering = TRUE)
      if (ihw_ok) {
        dfb <- as.data.frame(res_b)
        dfb$baseMean[is.na(dfb$baseMean)] <- 0
        fitb <- IHW::ihw(pvalue ~ baseMean, data = dfb, alpha = alpha)
        res_b$padj <- IHW::adj_pvalues(fitb)
      }
    } else {
      fit1b <- try(DESeq(dds_b, fitType="parametric", quiet=TRUE), silent=TRUE)
      if (inherits(fit1b, "try-error")) dds_b <- DESeq(dds_b, fitType="local", quiet=TRUE) else dds_b <- fit1b
      res_b <- results(dds_b,
                       contrast = c("group", experimental, control),
                       alpha = alpha,
                       independentFiltering = TRUE)
    }

    sig_b <- (!is.na(res_b$padj)) & res_b$padj <= alpha & abs(res_b$log2FoldChange) >= lfc_thr
    call_mat[rownames(res_b), b] <- as.integer(sig_b)
  }

  freq <- rowMeans(call_mat, na.rm = TRUE)
  consensus <- as.integer(freq >= args$consensus_threshold)
  fwrite(data.table(gene = rownames(call_mat),
                    consensus_fraction = freq,
                    consensus_call = consensus),
         file.path(args$output_dir, "bootstrap_consensus.tsv"), sep = "\t")
  writeLines(rownames(call_mat)[consensus == 1L],
             con = file.path(args$output_dir, "consensus_deseq2_genes_bootstrap.txt"))
}

# -----------------------------
# Done
# -----------------------------
writeLines(capture.output(sessionInfo()), con = file.path(args$output_dir, "sessionInfo.txt"))
message("DESeq2 module finished.")
