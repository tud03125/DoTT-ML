#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(DESeq2)
  # Optional helpers (auto-used if present; no required flags)
  ihw_ok    <- requireNamespace("IHW", quietly = TRUE)      # IHW FDR (preferred if present)
  apeglm_ok <- requireNamespace("apeglm", quietly = TRUE)   # LFC shrinkage
  ev_ok     <- requireNamespace("EnhancedVolcano", quietly = TRUE)
  gg_ok     <- requireNamespace("ggplot2", quietly = TRUE)
  edgeR_ok  <- requireNamespace("edgeR", quietly = TRUE)    # prefilter
  limma_ok  <- requireNamespace("limma", quietly = TRUE)    # removeBatchEffect for plotting
  bees_ok   <- requireNamespace("ggbeeswarm", quietly = TRUE)
  rastr_ok  <- requireNamespace("ggrastr", quietly = TRUE)
  readr_ok  <- requireNamespace("readr", quietly = TRUE)
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
# Prefiltering (edgeR::filterByExpr when available)
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
# Results (guarantee BH/IHW padj)
# -----------------------------
if (profile == "robust") {
  res <- results(
    dds, contrast = c("group", experimental, control),
    alpha = alpha, lfcThreshold = lfc_thr, altHypothesis = "greaterAbs",
    independentFiltering = TRUE
  )
} else {
  res <- results(
    dds, contrast = c("group", experimental, control),
    alpha = alpha, independentFiltering = TRUE
  )
}

# Prefer IHW if installed; otherwise enforce BH on non-NA p-values
if (ihw_ok) {
  df <- as.data.frame(res)
  df$baseMean[is.na(df$baseMean)] <- 0
  fit <- IHW::ihw(pvalue ~ baseMean, data = df, alpha = alpha)
  res$padj <- IHW::adj_pvalues(fit)
} else {
  notna <- !is.na(res$pvalue)
  res$padj[notna] <- p.adjust(res$pvalue[notna], method = "BH")
}

# -----------------------------
# optional LFC shrinkage for plots (FIXED)
# -----------------------------
shrunk <- NULL
if (apeglm_ok) {  # only attempt if apeglm is installed
  rn <- resultsNames(dds)
  # Build the expected coefficient name (works for design ~ group)
  coef_guess <- paste0("group_", make.names(experimental), "_vs_", make.names(control))
  coef_name <- if (coef_guess %in% rn) coef_guess else {
    cand <- rn[grep("^group_.*_vs_.*", rn)]
    if (length(cand)) cand[1] else NA_character_
  }

  if (!is.na(coef_name)) {
    # IMPORTANT: lfcShrink() is from DESeq2, not apeglm
    suppressWarnings({
      shrunk <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
    })
  } else {
    message("Could not determine coef for lfcShrink; skipping shrinkage.")
  }
} else {
  message("apeglm not installed; skipping LFC shrinkage.")
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
# Volcano (counts and plot)
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

# counts you asked for
n_total       <- nrow(res_df)
n_sig_total   <- sum(!is.na(res_df$padj) & res_df$padj < thr_p & abs(res_df$log2FoldChange) > thr_fc)
n_sig_exp_up  <- sum(!is.na(res_df$padj) & res_df$padj < thr_p & res_df$log2FoldChange >  thr_fc)
n_sig_ctrl_up <- sum(!is.na(res_df$padj) & res_df$padj < thr_p & res_df$log2FoldChange < -thr_fc)

fwrite(data.frame(
  total_tested = n_total,
  sig_total = n_sig_total,
  sig_up_experimental = n_sig_exp_up,
  sig_up_control = n_sig_ctrl_up,
  experimental = experimental,
  control = control,
  padj_cutoff = thr_p,
  lfc_cutoff = thr_fc
), file = file.path(args$output_dir, "volcano_counts.tsv"), sep = "\t")

cap_txt <- sprintf("Total=%d | padj<%.2g & |L2FC|>%g: %d | %s up: %d | %s up: %d",
                   n_total, thr_p, thr_fc, n_sig_total, experimental, n_sig_exp_up, control, n_sig_ctrl_up)

if (ev_ok) {
  suppressPackageStartupMessages(
    library(EnhancedVolcano)
  )
  vol_df <- if (!is.null(shrunk)) as.data.frame(shrunk) else res_df
  vol_df$gene <- rownames(vol_df)
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
  ) + ggplot2::labs(caption = cap_txt)
  if (gg_ok) {
    ggplot2::ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), vol,
                    width = 14, height = 10, units = "in")
    ggplot2::ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.png")), vol,
                    width = 14, height = 10, units = "in", dpi = 300)
  } else {
    svg(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), width=14, height=10)
    print(vol)
    dev.off()
  }
} else if (gg_ok) {
  suppressPackageStartupMessages(
    library(ggplot2)
  )
  ypadj <- res_df$padj
  ypadj[!is.na(ypadj) & ypadj == 0] <- .Machine$double.xmin
  res_df$neglog10padj <- -log10(ypadj)
  p_vol <- ggplot2::ggplot(res_df, ggplot2::aes(x = log2FoldChange, y = neglog10padj, color = direction)) +
    ggplot2::geom_point(alpha = 0.55, size = 0.8, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(-thr_fc, thr_fc), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(thr_p), linetype = "dashed") +
    ggplot2::labs(
      title = sprintf("%s Volcano (%s vs %s)", prefix, experimental, control),
      subtitle = "Significant: padj < 0.05 & |LFC| > 1",
      caption = cap_txt,
      x = sprintf("log2 fold change (%s / %s)", experimental, control),
      y = "-log10(padj)", color = NULL
    ) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "top")
  ggplot2::ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), p_vol,
                  width = 14, height = 10, units = "in")
  ggplot2::ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.png")), p_vol,
                  width = 14, height = 10, units = "in", dpi = 300)
}

# -----------------------------
# MA plot (kept)
# -----------------------------
if (gg_ok) {
  suppressPackageStartupMessages(
    library(ggplot2)
  )
  p_ma <- ggplot2::ggplot(res_df, ggplot2::aes(x = log10(baseMean + 1), y = log2FoldChange, color = direction)) +
    ggplot2::geom_point(alpha = 0.55, size = 0.8, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = c(-thr_fc, thr_fc), linetype = "dashed") +
    ggplot2::labs(
      title = sprintf("%s MA plot (%s vs %s)", prefix, experimental, control),
      subtitle = "Significant: padj < 0.05 & |LFC| > 1",
      x = "log10(baseMean + 1)",
      y = sprintf("log2 fold change (%s / %s)", experimental, control),
      color = NULL
    ) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "top")
  ggplot2::ggsave(file.path(args$output_dir, paste0(prefix, "_MA_plot.svg")), p_ma,
                  width = 14, height = 10, units = "in")
  ggplot2::ggsave(file.path(args$output_dir, paste0(prefix, "_MA_plot.png")), p_ma,
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

# -----------------------------
# Violin (points-only) + counts + Wilcoxon label
# -----------------------------
if (readr_ok && gg_ok && bees_ok) {
  suppressPackageStartupMessages({
    library(readr); library(dplyr); library(tidyr); library(stringr)
    library(ggplot2); library(ggbeeswarm)
    if (rastr_ok) library(ggrastr)
  })

  clean_names <- function(x) {
    x <- basename(x)
    x <- stringr::str_remove(x, "_sorted\\.bam$|\\.markdup\\.sorted\\.bam$|\\.bam$")
    x
  }

  # Read normalized counts and results we just wrote
  counts_csv <- file.path(args$output_dir, "normalized_counts.csv")
  results_csv <- file.path(args$output_dir, paste0(prefix, "_differential_analysis_results.csv"))
  cts <- readr::read_csv(counts_csv, show_col_types = FALSE)
  names(cts)[1] <- "gene"
  res_tbl <- readr::read_csv(results_csv, show_col_types = FALSE)
  names(res_tbl)[1] <- "gene"

  # keep only significant genes (padj & |LFC|)
  sig_res <- dplyr::filter(res_tbl, !is.na(padj), padj < alpha, abs(log2FoldChange) > lfc_thr)
  n_sig      <- nrow(sig_res)
  n_up_exp   <- sum(sig_res$log2FoldChange >  lfc_thr, na.rm = TRUE)
  n_up_ctrl  <- sum(sig_res$log2FoldChange < -lfc_thr, na.rm = TRUE)

  # Harmonize IDs, order samples left->right by control, experimental
  counts_samples <- clean_names(names(cts)[-1])
  names(cts) <- c("gene", counts_samples)
  ord <- c(colnames(mat)[idx_ctrl], colnames(mat)[idx_exp])
  ord <- clean_names(ord)
  cts <- cts[, c("gene", ord), drop=FALSE]

  # matrix -> log2
  expr_mat <- as.matrix(dplyr::semi_join(cts, dplyr::select(sig_res, gene), by = "gene")[, -1, drop=FALSE])
  rownames(expr_mat) <- dplyr::semi_join(cts, dplyr::select(sig_res, gene), by = "gene")$gene
  log_mat <- log2(expr_mat + 1)

  # colData for plotting
  plot_df <- data.frame(
    sample    = ord,
    condition = factor(c(rep(control, length(idx_ctrl)), rep(experimental, length(idx_exp))),
                       levels = c(control, experimental)),
    stringsAsFactors = FALSE, row.names = ord
  )

  use_log <- log_mat
  # optional batch effect removal for VISUALIZATION
  if (limma_ok) {
    # infer batch from "I"/"IA" prefix if present
    b <- ifelse(grepl("^IA", rownames(plot_df)), "IA",
         ifelse(grepl("^I",  rownames(plot_df)), "I", NA))
    if (length(na.omit(unique(b))) >= 2) {
      design <- model.matrix(~ condition, data = plot_df) # preserve condition effect
      use_log <- limma::removeBatchEffect(log_mat, batch = b, design = design)
    }
  }

  long_df <- as.data.frame(use_log) |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expr") |>
    dplyr::mutate(condition = plot_df[sample, "condition", drop=TRUE],
                  condition = factor(condition, levels = c(control, experimental)))
  
  # ====================== IHW-FDR for multiple Wilcoxon tests ======================
  # REQUIREMENT: a "panel" column telling us which panel each gene belongs to.
  # If you already have a panel vector per gene, join it here.
  # Example stub: put every gene into the same panel if you don't yet have categories.
  # Replace this with your real categories (e.g., "3UTR_extended", "antisense", "intronic", ...).
  if (!"panel" %in% names(long_df)) {
    long_df$panel <- "3UTR_extended"
  }
  
  # Bring in baseMean for a covariate (median baseMean per panel is a reasonable IHW covariate)
  # Join DESeq2 results baseMean onto long_df by gene:
  bm <- res_tbl[, c("gene","baseMean")]   # res_tbl was read from your *_differential_analysis_results.csv
  long_df <- dplyr::left_join(long_df, bm, by = "gene")
  
  # Compute one Wilcoxon p-value per panel
  panels <- split(long_df, long_df$panel)
  p_vec  <- vapply(panels, function(df) {
    stats::wilcox.test(expr ~ condition, data = df, exact = FALSE)$p.value
  }, numeric(1))
  
  # Choose an IHW covariate per panel (must be independent of null p's, but informative for power).
  # Here we use panel-level median baseMean across genes:
  covar <- vapply(panels, function(df) median(df$baseMean, na.rm = TRUE), numeric(1))
  
  # IHW adjustment across panels
  if (requireNamespace("IHW", quietly = TRUE) && length(p_vec) > 1L) {
    fit <- IHW::ihw(p_vec ~ covar, alpha = 0.05)
    padj_panels <- IHW::adj_pvalues(fit)
    fdr_tbl <- data.frame(panel = names(p_vec),
                          p_raw = as.numeric(p_vec),
                          padj_ihw = as.numeric(padj_panels),
                          covariate = as.numeric(covar),
                          stringsAsFactors = FALSE)
    # Keep for annotation below
  } else {
    # Fallbacks: if only one panel, IHW isn't defined; keep raw p (or use BH on length-1 ? identical)
    fdr_tbl <- data.frame(panel = names(p_vec),
                          p_raw = as.numeric(p_vec),
                          padj_ihw = p.adjust(p_vec, method = "BH"),
                          covariate = as.numeric(covar),
                          stringsAsFactors = FALSE)
  }
  # ================================================================================
  
  # palette
  pal  <- c(setNames("#2CA9E1", control), setNames("#E76F51", experimental))

  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = condition, y = expr, colour = condition)) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(
      x = NULL,
      y = "log2(normalized DoTT counts + 1)",
      title = "DoTT expression across conditions (significant genes only)",
      subtitle = sprintf("%s vs %s | n(sig genes)=%d  [Up in %s: %d | Up in %s: %d]",
                         control, experimental, n_sig, experimental, n_up_exp, control, n_up_ctrl)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                 plot.title.position = "plot",
                 plot.margin = ggplot2::margin(t = 22, r = 12, b = 10, l = 10)) +
    ggplot2::facet_wrap(~ panel, scales = "free_y")   # <— add facets by panel

  pts <- ggbeeswarm::geom_quasirandom(width = 0.18, alpha = 0.65, size = 0.55, show.legend = FALSE)
  if (rastr_ok) p <- p + ggrastr::rasterise(pts, dpi = 300) else p <- p + pts

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

  # Panel-wise IHW-FDR annotation (fixed)
  ann_df <- aggregate(expr ~ panel, data = long_df, FUN = max, na.rm = TRUE)
  ann_df$y <- ann_df$expr * 0.98
  
  # keep y!  (either of these is fine)
  ann_df <- dplyr::left_join(ann_df, fdr_tbl, by = "panel")
  # or:
  # ann_df <- merge(ann_df[, c("panel","y")], fdr_tbl, by = "panel", all.x = TRUE)
  
  ann_df$label <- ifelse(is.finite(ann_df$padj_ihw),
                         sprintf("Wilcoxon, FDR (IHW) = %.3g", ann_df$padj_ihw),
                         sprintf("Wilcoxon, p (raw) = %.3g", ann_df$p_raw))
  
  p <- p + ggplot2::geom_text(data = ann_df,
                              ggplot2::aes(x = 1.5, y = y, label = label),
                              inherit.aes = FALSE, size = 3.6)

  ggplot2::ggsave(file.path(args$output_dir, sprintf("DoTT_violin_%s_%s.svg", experimental, control)),
                  p, device = "svg", width = 14, height = 10, dpi = 300)
}

# -----------------------------
# Bootstrap consensus (unchanged core)
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
    } else {
      fit1b <- try(DESeq(dds_b, fitType="parametric", quiet=TRUE), silent=TRUE)
      if (inherits(fit1b, "try-error")) dds_b <- DESeq(dds_b, fitType="local", quiet=TRUE) else dds_b <- fit1b
      res_b <- results(dds_b,
                       contrast = c("group", experimental, control),
                       alpha = alpha,
                       independentFiltering = TRUE)
    }

    # enforce BH/IHW on bootstrap too
    if (ihw_ok) {
      dfb <- as.data.frame(res_b); dfb$baseMean[is.na(dfb$baseMean)] <- 0
      fitb <- IHW::ihw(pvalue ~ baseMean, data = dfb, alpha = alpha)
      res_b$padj <- IHW::adj_pvalues(fitb)
    } else {
      notna <- !is.na(res_b$pvalue)
      res_b$padj[notna] <- p.adjust(res_b$pvalue[notna], method = "BH")
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
