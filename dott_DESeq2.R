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
# Args: unchanged from your script
# -----------------------------
parser <- ArgumentParser()
parser$add_argument("--counts_file", required=TRUE, help="Cleaned counts file")
parser$add_argument("--output_dir", required=TRUE, help="Output directory")
parser$add_argument("--conditions", required=TRUE, help="Comma-separated list of condition labels")
parser$add_argument("--bootstrap", type="logical", default=FALSE, help="Whether to perform bootstrapping")
parser$add_argument("--n_boot", type="integer", default=100, help="Number of bootstrap iterations")
parser$add_argument("--consensus_threshold", type="double", default=0.5, help="Fraction threshold for consensus calls")
# NOTE: no new required CLI flags; profile is chosen via env var
args <- parser$parse_args()

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Tunables (fixed; same thresholds as your original pipeline)
# -----------------------------
alpha   <- 0.05
lfc_thr <- 1.0        # for output gating and robust thresholded test
set.seed(1L)
prefix <- "3_UTR_extended"  # legacy filename prefix

# Profile selection: "robust" (default) or "classic"
profile <- tolower(Sys.getenv("DOTT_DESEQ2_PROFILE", "robust"))
if (!profile %in% c("robust","classic")) {
  warning("Unknown DOTT_DESEQ2_PROFILE='", profile, "'. Falling back to 'robust'.")
  profile <- "robust"
}
message("DESeq2 profile: ", profile)

# -----------------------------
# Load counts (first col = gene id or 'Geneid')
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

# Ensure unique sample names (important for bootstrapping)
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
# - robust: design-aware filtering (edgeR::filterByExpr) if available
# - classic: no prefilter (mimics your original behavior)
# -----------------------------
if (profile == "robust" && edgeR_ok) {
  y <- edgeR::DGEList(counts = counts(dds), group = coldata$group)
  keep <- edgeR::filterByExpr(y, group = coldata$group)  # recommended for RNA-seq
  if (sum(keep) < 200) keep <- edgeR::filterByExpr(y, group = coldata$group, min.count = 1)
} else {
  keep <- rep(TRUE, nrow(dds))  # preserve original behavior
}
dds <- dds[keep, , drop=FALSE]

# -----------------------------
# Fitting (with robust fallbacks)
# - robust: parametric -> local -> gene-wise + Wald
# - classic: parametric; if it *fails*, fall back to local (to avoid crashes)
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
  return(dds2)
}

if (profile == "robust") {
  dds <- fit_deseq_with_fallback(dds, quiet = TRUE)
} else {
  fit1 <- try(DESeq(dds, fitType = "parametric", quiet = TRUE), silent = TRUE)
  if (inherits(fit1, "try-error")) {
    message("Classic profile: parametric fit failed; trying 'local' to avoid error.")
    dds <- DESeq(dds, fitType = "local", quiet = TRUE)
  } else {
    dds <- fit1
  }
}

# -----------------------------
# Results
# - robust: thresholded Wald test aligned to |LFC| > lfc_thr (1.0)
# - classic: standard two-sided Wald test at LFC=0; apply |LFC|>=1 post hoc
# -----------------------------
if (profile == "robust") {
  res <- results(
    dds,
    contrast = c("group", experimental, control),
    alpha = alpha,
    lfcThreshold = lfc_thr,
    altHypothesis = "greaterAbs",
    independentFiltering = TRUE
  )
  # Adaptive FDR (IHW) if available
  if (ihw_ok) {
    library(IHW)
    df <- as.data.frame(res)
    df$baseMean[is.na(df$baseMean)] <- 0
    fit <- ihw(pvalue ~ baseMean, data = df, alpha = alpha)
    res$padj <- adj_pvalues(fit)
  }
} else {
  # classic
  res <- results(
    dds,
    contrast = c("group", experimental, control),
    alpha = alpha,
    independentFiltering = TRUE
  )
}

# -----------------------------
# Optional LFC shrinkage for plots (not used for calling)
# apeglm requires 'coef' (not 'contrast')
# -----------------------------
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
# Normalized counts + group means
# -----------------------------
norm_counts <- counts(dds, normalized = TRUE)
fwrite(data.table(gene = rownames(norm_counts), as.data.frame(norm_counts)),
       file.path(args$output_dir, "normalized_counts.csv"), sep = ",")

idx_ctrl <- which(coldata$group == control)
idx_exp  <- which(coldata$group == experimental)
mean_ctrl <- rowMeans(norm_counts[, idx_ctrl, drop=FALSE])
mean_exp  <- rowMeans(norm_counts[, idx_exp,  drop=FALSE])

# -----------------------------
# Write legacy result files
# -----------------------------
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$mean_control <- mean_ctrl[res_df$gene]
res_df$mean_experimental <- mean_exp[res_df$gene]

res_out <- res_df[, c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj",
                      "mean_control","mean_experimental")]
fwrite(res_out,
       file.path(args$output_dir, paste0(prefix, "_differential_analysis_results.csv")),
       sep = ",")

# Significant tables at your gates (padj<=alpha & |LFC|>=1)
sig_abs <- subset(res_out, !is.na(padj) & padj <= alpha & abs(log2FoldChange) >= lfc_thr)
fwrite(sig_abs, file.path(args$output_dir, "significant_extended_genes.csv"), sep = ",")
fwrite(sig_abs, file.path(args$output_dir, "absolute_significant_extended_genes_with_individual_means.csv"), sep = ",")

# -----------------------------
# Plots (SVG): MA and Volcano
# -----------------------------

# 1) classify once, reuse in both plots
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

# -----------------
# MA plot
# -----------------
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
         width = 7, height = 5, units = "in")
  ggsave(file.path(args$output_dir, paste0(prefix, "_MA_plot.png")), p_ma,
         width = 7, height = 5, units = "in", dpi = 300)
} else {
  # base-R fallback with colors/legend
  cols <- setNames(c("grey60","firebrick","steelblue"),
                   c("NS", paste0("Up in ", experimental), paste0("Up in ", control)))
  svg(file.path(args$output_dir, paste0(prefix, "_MA_plot.svg")), width=7, height=5)
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

# -----------------
# Volcano plot
# -----------------
# Prefer padj on y (better matches significance gating) and SAVE the ggplot object explicitly.
# Use shrunk LFCs if available (for visualization only).
vol_df <- if (!is.null(shrunk)) as.data.frame(shrunk) else res_df
vol_df$gene <- rownames(vol_df)

if (ev_ok) {
  # EnhancedVolcano returns a ggplot object — use ggsave to avoid blank SVGs.
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
           width = 7, height = 5, units = "in")
    ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.png")), vol,
           width = 7, height = 5, units = "in", dpi = 300)
  } else {
    svg(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), width=7, height=5)
    print(vol)  # critical so the device isn't blank
    dev.off()
  }
} else if (gg_ok) {
  # pure ggplot volcano with padj, colored like MA
  library(ggplot2)
  # guard against padj==0 -> Inf on -log10
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
         width = 7, height = 5, units = "in")
  ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.png")), p_vol,
         width = 7, height = 5, units = "in", dpi = 300)
} else {
  # base-R fallback using padj and colored points
  cols <- setNames(c("grey60","firebrick","steelblue"),
                   c("NS", paste0("Up in ", experimental), paste0("Up in ", control)))
  ypadj <- res_df$padj
  ypadj[!is.na(ypadj) & ypadj == 0] <- .Machine$double.xmin
  svg(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), width=7, height=5)
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
# Optional bootstrap consensus (unchanged outputs)
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
    # fit per profile
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
      # classic inside bootstrap
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
