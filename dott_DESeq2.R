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
if (gg_ok) {
  library(ggplot2)
  p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
    geom_point(alpha = 0.5, size = 0.7) +
    geom_hline(yintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
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
    abline(h=c(-lfc_thr, lfc_thr), lty=2)
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
                  pCutoff = alpha, FCcutoff = lfc_thr,
                  title = sprintf("%s Volcano (%s vs %s)", prefix, experimental, control),
                  labSize = 3)
  dev.off()
} else if (gg_ok) {
  library(ggplot2)
  volcano_df <- if (!is.null(shrunk)) as.data.frame(shrunk) else as.data.frame(res)
  volcano_df$neglog10p <- -log10(volcano_df$pvalue)
  p_vol <- ggplot(volcano_df, aes(x = log2FoldChange, y = neglog10p)) +
    geom_point(alpha = 0.6, size = 0.7) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    labs(title = sprintf("%s Volcano (%s vs %s)", prefix, experimental, control),
         x = "log2 fold change", y = "-log10(pvalue)") +
    theme_bw()
  ggsave(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), p_vol,
         width = 7, height = 5, dpi = 300, device = "svg")
} else {
  svg(file.path(args$output_dir, paste0(prefix, "_Volcano_plot.svg")), width=7, height=5)
  with(as.data.frame(res), {
    plot(log2FoldChange, -log10(pvalue), pch=20, cex=0.6,
         xlab="log2 fold change", ylab="-log10(pvalue)",
         main=sprintf("%s Volcano (%s vs %s)", prefix, experimental, control))
    abline(v=c(-lfc_thr, lfc_thr), lty=2); abline(h=-log10(alpha), lty=2)
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
