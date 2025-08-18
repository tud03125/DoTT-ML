#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(DESeq2)
  ihw_ok    <- requireNamespace("IHW", quietly = TRUE)
  apeglm_ok <- requireNamespace("apeglm", quietly = TRUE)
  ev_ok     <- requireNamespace("EnhancedVolcano", quietly = TRUE)
  gg_ok     <- requireNamespace("ggplot2", quietly = TRUE)
  edgeR_ok  <- requireNamespace("edgeR", quietly = TRUE)
})

# -----------------------------
# CLI (unchanged)
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
# Tunables (no new CLI flags)
# -----------------------------
alpha <- 0.05
mode <- "classic"        # default: classic test vs 0; set "thresholded" to test |LFC| > lfc_test_thr
lfc_test_thr <- 0.5      # used only in thresholded mode
lfc_call_thr <- 1.0      # final reporting gate (padj<=alpha & |LFC|>=lfc_call_thr)
use_one_sided_up <- FALSE  # if TRUE and using thresholded mode, test "greater" (up in experimental)

set.seed(1L)
prefix <- "3_UTR_extended"

# -----------------------------
# Load counts (first column is gene id or "Geneid")
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

# ensure unique sample names
if (any(duplicated(colnames(mat)))) {
  colnames(mat) <- make.unique(colnames(mat), sep = ".dup")
}

# -----------------------------
# Design
# -----------------------------
conds <- strsplit(args$conditions, ",")[[1]]
stopifnot(length(conds) == ncol(mat))
u <- unique(conds); stopifnot(length(u) == 2)
control_label <- u[1]; experimental_label <- u[2]

coldata <- data.frame(sample = colnames(mat),
                      group  = factor(conds, levels = c(control_label, experimental_label)),
                      row.names = colnames(mat))

# >>> sanitize factor levels up front (maps e.g. "HSV-1" -> "HSV.1")
levels(coldata$group) <- make.names(levels(coldata$group))
# re-derive control/experimental AFTER sanitization (order preserved)
control      <- levels(coldata$group)[1]
experimental <- levels(coldata$group)[2]

dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldata, design = ~ group)

# -----------------------------
# Prefilter (prefer edgeR::filterByExpr; else mild row-sum)
# -----------------------------
if (edgeR_ok) {
  y <- edgeR::DGEList(counts = counts(dds), group = coldata$group)
  keep <- edgeR::filterByExpr(y, group = coldata$group)   # CPM-aware, design-aware
} else {
  keep <- rowSums(counts(dds)) >= 5
}
dds <- dds[keep, , drop=FALSE]
n_kept <- nrow(dds); n_total <- nrow(mat)

# -----------------------------
# Fit
# -----------------------------
dds <- DESeq(dds, quiet = TRUE)

# -----------------------------
# Results
# -----------------------------
if (mode == "thresholded") {
  altH <- if (use_one_sided_up) "greater" else "greaterAbs"
  res <- results(dds,
                 contrast = c("group", experimental, control),
                 alpha = alpha,
                 lfcThreshold = lfc_test_thr,
                 altHypothesis = altH,
                 independentFiltering = TRUE,
                 cooksCutoff = FALSE)
} else {
  res <- results(dds,
                 contrast = c("group", experimental, control),
                 alpha = alpha,
                 independentFiltering = TRUE,
                 cooksCutoff = FALSE)
}

# Save BH padj before optional IHW
padj_BH <- res$padj
nsig_bh_at_call <- sum(!is.na(padj_BH) & padj_BH <= alpha & abs(res$log2FoldChange) >= lfc_call_thr)

# Adaptive FDR via IHW with guard (fallback to BH if it slashes yield too much)
used_ihw <- FALSE
if (ihw_ok) {
  library(IHW)
  df <- as.data.frame(res); df$baseMean[is.na(df$baseMean)] <- 0
  fit <- ihw(pvalue ~ baseMean, data = df, alpha = alpha)
  padj_IHW <- adj_pvalues(fit)
  nsig_ihw_at_call <- sum(!is.na(padj_IHW) & padj_IHW <= alpha & abs(res$log2FoldChange) >= lfc_call_thr)
  if (nsig_ihw_at_call >= 0.8 * nsig_bh_at_call) {
    res$padj <- padj_IHW; used_ihw <- TRUE
  }
}

# -----------------------------
# Optional LFC shrinkage for plots/ranking (apeglm requires coef)
# -----------------------------
shrunk <- NULL
if (apeglm_ok) {
  library(apeglm)
  rn <- resultsNames(dds)  # inspect to see the sanitized coef names
  coef_name <- paste0("group_", experimental, "_vs_", control)
  if (coef_name %in% rn) {
    shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  } else {
    # ultra-rare fallback: try "normal" shrinkage with contrast *only* for plotting
    # (calls downstream still based on 'res')
    shrunk <- tryCatch(
      lfcShrink(dds, contrast = c("group", experimental, control), type = "normal"),
      error = function(e) NULL
    )
  }
}

# -----------------------------
# Normalized counts & group means
# -----------------------------
norm_counts <- counts(dds, normalized = TRUE)
fwrite(data.table(gene = rownames(norm_counts), as.data.frame(norm_counts)),
       file.path(args$output_dir, "normalized_counts.csv"), sep = ",")

idx_ctrl <- which(coldata$group == control)
idx_exp  <- which(coldata$group == experimental)
mean_ctrl <- rowMeans(norm_counts[, idx_ctrl, drop = FALSE])
mean_exp  <- rowMeans(norm_counts[, idx_exp,  drop = FALSE])

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

sig_abs <- subset(res_out, !is.na(padj) & padj <= alpha & abs(log2FoldChange) >= lfc_call_thr)
fwrite(sig_abs, file.path(args$output_dir, "significant_extended_genes.csv"), sep = ",")
fwrite(sig_abs, file.path(args$output_dir, "absolute_significant_extended_genes_with_individual_means.csv"), sep = ",")

# -----------------------------
# Plots
# -----------------------------
if (gg_ok) {
  library(ggplot2)
  p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
    geom_point(alpha = 0.5, size = 0.7) +
    geom_hline(yintercept = c(-lfc_call_thr, lfc_call_thr), linetype = "dashed") +
    labs(title = sprintf("%s MA plot (%s vs %s)", prefix, experimental, control),
         x = "log10(baseMean + 1)", y = "log2 fold change") +
    theme_bw()
  ggsave(file.path(args$output_dir, paste0(prefix, "_MA_plot.svg")), p_ma,
         width = 6, height = 4, dpi = 300, device = "svg")
} else {
  svg(file.path(args$output_dir, paste0(prefix, "_MA_plot.svg")), width=6, height=4)
  with(res_df, {
    plot(log10(baseMean+1), log2FoldChange, pch=20, cex=0.6,
         xlab="log10(baseMean+1)", ylab="log2 fold change",
         main=sprintf("%s MA plot (%s vs %s)", prefix, experimental, control))
    abline(h=c(-lfc_call_thr, lfc_call_thr), lty=2)
  })
  dev.off()
}

if (ev_ok) {
  library(EnhancedVolcano)
  volcano_df <- if (!is.null(shrunk)) as.data.frame(shrunk) else as.data.frame(res)
  volcano_df$gene <- rownames(volcano_df)
  svg(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), width=7, height=5)
  EnhancedVolcano(volcano_df,
                  lab = volcano_df$gene,
                  x = 'log2FoldChange', y = 'pvalue',
                  pCutoff = alpha, FCcutoff = lfc_call_thr,
                  title = sprintf("%s Volcano (%s vs %s)", prefix, experimental, control),
                  labSize = 3)
  dev.off()
} else if (gg_ok) {
  library(ggplot2)
  volcano_df <- if (!is.null(shrunk)) as.data.frame(shrunk) else as.data.frame(res)
  volcano_df$neglog10p <- -log10(volcano_df$pvalue)
  p_vol <- ggplot(volcano_df, aes(x = log2FoldChange, y = neglog10p)) +
    geom_point(alpha = 0.6, size = 0.7) +
    geom_vline(xintercept = c(-lfc_call_thr, lfc_call_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    labs(title = sprintf("%s Volcano (%s vs %s)", prefix, experimental, control),
         x = "log2 fold change", y = "-log10(pvalue)") +
    theme_bw()
  ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), p_vol,
         width = 7, height = 5, dpi = 300, device = "svg")
}

# -----------------------------
# Optional bootstrap consensus
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

    m_b <- cbind(mat[, rs_ctrl, drop = FALSE], mat[, rs_exp, drop = FALSE])

    drawn_names <- c(colnames(mat)[rs_ctrl], colnames(mat)[rs_exp])
    colnames(m_b) <- make.unique(drawn_names, sep = ".rep")

    cd_b <- data.frame(
      group = factor(c(rep(control, length(rs_ctrl)),
                       rep(experimental, length(rs_exp))),
                     levels = c(control, experimental)),
      row.names = colnames(m_b),
      stringsAsFactors = FALSE
    )

    dds_b <- DESeqDataSetFromMatrix(countData = m_b, colData = cd_b, design = ~ group)

    if (edgeR_ok) {
      yb <- edgeR::DGEList(counts = counts(dds_b), group = cd_b$group)
      keep_b <- edgeR::filterByExpr(yb, group = cd_b$group)
    } else {
      keep_b <- rowSums(counts(dds_b)) >= 5
    }
    dds_b <- dds_b[keep_b, , drop = FALSE]
    dds_b <- DESeq(dds_b, quiet = TRUE)

    if (mode == "thresholded") {
      altH_b <- if (use_one_sided_up) "greater" else "greaterAbs"
      res_b <- results(dds_b, contrast = c("group", experimental, control),
                       alpha = alpha, lfcThreshold = lfc_test_thr,
                       altHypothesis = altH_b, independentFiltering = TRUE,
                       cooksCutoff = FALSE)
    } else {
      res_b <- results(dds_b, contrast = c("group", experimental, control),
                       alpha = alpha, independentFiltering = TRUE,
                       cooksCutoff = FALSE)
    }

    padj_bh_b <- res_b$padj
    if (ihw_ok) {
      dfb <- as.data.frame(res_b); dfb$baseMean[is.na(dfb$baseMean)] <- 0
      fitb <- IHW::ihw(pvalue ~ baseMean, data = dfb, alpha = alpha)
      padj_ihw_b <- IHW::adj_pvalues(fitb)
      nsig_bh_b  <- sum(!is.na(padj_bh_b)  & padj_bh_b  <= alpha & abs(res_b$log2FoldChange) >= lfc_call_thr)
      nsig_ihw_b <- sum(!is.na(padj_ihw_b) & padj_ihw_b <= alpha & abs(res_b$log2FoldChange) >= lfc_call_thr)
      res_b$padj <- if (nsig_ihw_b >= 0.8 * nsig_bh_b) padj_ihw_b else padj_bh_b
    }

    sig_b <- (!is.na(res_b$padj)) & res_b$padj <= alpha & abs(res_b$log2FoldChange) >= lfc_call_thr
    call_mat[rownames(res_b), b] <- as.integer(sig_b)
  }

  freq <- rowMeans(call_mat)
  consensus <- as.integer(freq >= args$consensus_threshold)

  data.table(gene = rownames(call_mat),
             consensus_fraction = freq,
             consensus_call = consensus) |>
    fwrite(file.path(args$output_dir, "bootstrap_consensus.tsv"), sep = "\t")

  writeLines(rownames(call_mat)[consensus == 1L],
             con = file.path(args$output_dir, "consensus_deseq2_genes_bootstrap.txt"))
}

# -----------------------------
# Run summary (sanity checks)
# -----------------------------
nsig_final <- sum(!is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= lfc_call_thr)
summary_lines <- c(
  sprintf("mode=%s  lfc_test_thr=%.3f  call_|LFC|>=%.3f  one_sided_up=%s",
          mode, lfc_test_thr, lfc_call_thr, use_one_sided_up),
  sprintf("kept_genes=%d / %d (prefilter=%s)", n_kept, n_total,
          if (edgeR_ok) "edgeR::filterByExpr" else "rowSums>=5"),
  sprintf("independentFiltering=TRUE  cooksCutoff=FALSE"),
  sprintf("IHW_used=%s  nsig_at_call (BH)=%d  nsig_at_call (final)=%d",
          used_ihw, nsig_bh_at_call, nsig_final)
)
writeLines(summary_lines, con = file.path(args$output_dir, "deseq2_run_summary.txt"))
writeLines(capture.output(sessionInfo()), con = file.path(args$output_dir, "sessionInfo.txt"))
message("DESeq2 module finished.")
