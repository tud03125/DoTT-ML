#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unsupervised analyses for the DoTT pipeline.

Outputs (in output_dir):
  dott_genes_with_conditions.csv
  dott_genes_with_conditions_filtered.csv
  enrichment_comparison.csv
  enrichment_summary.txt
  replicate_consistency_CV.csv
  unsup_PCA_samples.csv
  unsup_PCA_samples.png
  unsup_KMeans_results.json
"""

import os
import json
import warnings
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

import matplotlib
matplotlib.use("Agg")  # HPC-safe plotting
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=RuntimeWarning)


# ---------------------------
# Helpers
# ---------------------------

def _safe_makedirs(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def _write_text(path: str, text: str) -> None:
    with open(path, "w") as f:
        f.write(text)


# ---------------------------
# Enrichment
# ---------------------------

def compare_enrichment_across_conditions(predictions_csv: str, output_dir: str, experimental_condition: str) -> Optional[str]:

    """
    Fisher's exact on predicted_DE (1/0) by condition (experimental vs control).
    Robust to degenerate cases where one row/column is all zeros.
    """
    df = pd.read_csv(predictions_csv)
    if "condition" not in df.columns or "predicted_DE" not in df.columns:
        print("Missing 'condition' or 'predicted_DE'; cannot perform enrichment analysis.")
        return None

    # Build contingency and FORCE a 2x2 with rows [exp, control] and cols [0,1]
    contingency = pd.crosstab(df["condition"], df["predicted_DE"])
    other_conditions = [c for c in contingency.index if c != experimental_condition]
    control_condition = other_conditions[0] if other_conditions else "Control"

    contingency = contingency.reindex(index=[experimental_condition, control_condition],
                                      columns=[0, 1],
                                      fill_value=0)

    # Pull counts safely
    exp_nonDE = int(contingency.loc[experimental_condition, 0])
    exp_DE    = int(contingency.loc[experimental_condition, 1])
    ctrl_nonDE = int(contingency.loc[control_condition, 0])
    ctrl_DE    = int(contingency.loc[control_condition, 1])

    # Surface the one-sided case
    if contingency.loc[control_condition].sum() == 0:
        print("Note: no genes assigned to the control condition; contingency will contain zeros.")

    # Fisher’s exact (2x2). Zeros are allowed; OR may be 0/inf/NaN. 
    # (scipy.stats.fisher_exact expects a 2x2 array of non-negative ints.)
    table = [[exp_DE, exp_nonDE], [ctrl_DE, ctrl_nonDE]]
    try:
        oddsratio, pvalue = fisher_exact(table)
    except Exception as e:
        print(f"Fisher exact failed on table {table}: {e}. Setting OR=nan, p=1.0")
        oddsratio, pvalue = float("nan"), 1.0

    # (Optional) write the 2x2 for quick inspection
    contingency.to_csv(os.path.join(output_dir, "enrichment_contingency_table.csv"))

    out_csv = os.path.join(output_dir, "enrichment_comparison.csv")
    pd.DataFrame({
        "experimental_condition": [experimental_condition],
        "control_condition": [control_condition],
        "exp_DE": [exp_DE],
        "exp_nonDE": [exp_nonDE],
        "ctrl_DE": [ctrl_DE],
        "ctrl_nonDE": [ctrl_nonDE],
        "odds_ratio": [oddsratio],
        "p_value": [pvalue]
    }).to_csv(out_csv, index=False)

    _write_text(os.path.join(output_dir, "enrichment_summary.txt"),
                "Fisher's exact test p-value: {}\nOdds Ratio: {}\n".format(pvalue, oddsratio))

    print("Enrichment written:", out_csv)
    return out_csv

# ---------------------------
# Replicate consistency (CV on log1p counts)
# ---------------------------

def assess_replicate_consistency(norm_counts_df: pd.DataFrame, predicted_genes: List[str], experimental_samples: List[str], output_dir: str) -> str:

    """
    CV across experimental replicates only, using log1p counts for stability.
    """
    genes = [g for g in predicted_genes if g in norm_counts_df.index]
    samples = [s for s in experimental_samples if s in norm_counts_df.columns]
    if not genes or not samples:
        print("No overlap between genes/samples and normalized_counts; skipping CV.")
        return ""

    sub = norm_counts_df.loc[genes, samples].copy()
    sub_log = np.log1p(sub)

    means = sub_log.mean(axis=1)
    sds = sub_log.std(axis=1, ddof=1)
    cv = sds / means.replace(0, np.nan)

    out = pd.DataFrame({
        "Gene": means.index,
        "log1p_mean": means.values,
        "log1p_sd": sds.values,
        "CV_log1p": cv.values
    }).sort_values("CV_log1p")

    out_path = os.path.join(output_dir, "replicate_consistency_CV.csv")
    out.to_csv(out_path, index=False)
    print("Replicate CV written:", out_path)
    return out_path


# ---------------------------
# PCA + KMeans on samples
# ---------------------------

def _most_variable_genes(df: pd.DataFrame, top_n: int = 1000) -> List[str]:
    variances = df.var(axis=1)
    return variances.sort_values(ascending=False).head(min(top_n, len(variances))).index.tolist()

def run_pca_and_kmeans(norm_counts_df: pd.DataFrame, sample_conditions: List[str], output_dir: str, predicted_genes: Optional[List[str]] = None, top_var_genes: int = 1000, n_clusters: int = 2) -> Tuple[str, str, str]:

    """
    PCA on samples (columns = samples, rows = genes).
    KMeans on the standardized sample matrix; compute silhouette.
    """
    _safe_makedirs(output_dir)

    X = norm_counts_df.copy()
    if predicted_genes:
        keep = [g for g in predicted_genes if g in X.index]
        if len(keep) >= 3:
            X = X.loc[keep]
        else:
            print("Too few predicted genes found; using all genes for PCA/KMeans.")

    keep_genes = _most_variable_genes(X, top_n=top_var_genes)
    X = X.loc[keep_genes]

    # log1p then standardize features (genes) across samples
    X_log = np.log1p(X)
    scaler = StandardScaler(with_mean=True, with_std=True)
    X_scaled = scaler.fit_transform(X_log.T)  # shape: samples x genes

    pca = PCA(n_components=2, random_state=42)
    pcs = pca.fit_transform(X_scaled)
    explained = pca.explained_variance_ratio_.tolist()

    samples = list(X.columns)
    if len(samples) != len(sample_conditions):
        print("Warning: len(samples) != len(sample_conditions); colors may be misaligned.")
    conds_for_samples = sample_conditions[:len(samples)]

    pca_df = pd.DataFrame({
        "sample": samples,
        "condition": conds_for_samples,
        "PC1": pcs[:, 0],
        "PC2": pcs[:, 1]
    })
    pca_csv = os.path.join(output_dir, "unsup_PCA_samples.csv")
    pca_df.to_csv(pca_csv, index=False)

    # Plot
    pca_png = os.path.join(output_dir, "unsup_PCA_samples.png")
    plt.figure(figsize=(7, 6))
    for cond in sorted(set(conds_for_samples)):
        mask = [c == cond for c in conds_for_samples]
        plt.scatter(pcs[mask, 0], pcs[mask, 1], label=cond)
    plt.xlabel("PC1 ({:.1f}% var)".format(100.0 * explained[0]))
    plt.ylabel("PC2 ({:.1f}% var)".format(100.0 * explained[1]))
    plt.title("PCA (samples)")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(pca_png, dpi=150)
    plt.close()

    # KMeans + silhouette (use explicit n_init for compatibility across sklearn versions)
    km = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
    clusters = km.fit_predict(X_scaled)
    sil = float("nan")
    try:
        if len(set(clusters)) > 1:
            sil = float(silhouette_score(X_scaled, clusters))
    except Exception:
        pass

    kmeans_json = os.path.join(output_dir, "unsup_KMeans_results.json")
    with open(kmeans_json, "w") as f:
        json.dump({
            "n_samples": int(X_scaled.shape[0]),
            "n_features_genes": int(X_scaled.shape[1]),
            "n_clusters": n_clusters,
            "silhouette_score": sil,  # -1..1; higher is better
            "cluster_assignments": {s: int(c) for s, c in zip(samples, clusters)},
            "pca_variance_explained": explained
        }, f, indent=2)

    print("PCA CSV:", pca_csv)
    print("PCA PNG:", pca_png)
    print("KMeans JSON:", kmeans_json)
    return pca_csv, pca_png, kmeans_json


# ---------------------------
# Means & condition assignment
# ---------------------------

def _compute_group_means_from_norm_counts(norm_counts_df: pd.DataFrame, conditions_list: List[str], experimental_condition: str) -> pd.DataFrame:

    """
    Compute per-gene means from normalized counts using provided sample conditions.
    Returns DataFrame with columns ['exp_mean','ctrl_mean'] indexed by gene.
    """
    if len(conditions_list) != norm_counts_df.shape[1]:
        raise ValueError("Length of --conditions does not match columns in normalized_counts.csv")

    col_to_cond = dict(zip(norm_counts_df.columns, conditions_list))
    exp_cols = [c for c in norm_counts_df.columns if col_to_cond[c] == experimental_condition]
    ctrl_cols = [c for c in norm_counts_df.columns if col_to_cond[c] != experimental_condition]
    if not exp_cols or not ctrl_cols:
        raise ValueError("Could not split samples into experimental/control groups from --conditions.")

    exp_mean = norm_counts_df[exp_cols].mean(axis=1)
    ctrl_mean = norm_counts_df[ctrl_cols].mean(axis=1)
    return pd.DataFrame({"exp_mean": exp_mean, "ctrl_mean": ctrl_mean})


def _assign_condition(dott_df: pd.DataFrame,
                      means_df: pd.DataFrame,
                      experimental_condition: str) -> pd.DataFrame:
    """
    Join exp/ctrl means to dott_df (by index=gene) and set:
      condition = experimental_condition if exp_mean > ctrl_mean else 'Control'
      predicted_DE = 1 if condition==experimental_condition else 0
    """
    merged = dott_df.join(means_df, how="left")
    if merged[["exp_mean", "ctrl_mean"]].isna().all().all():
        raise ValueError("Failed to compute group means (all NaN).")

    merged["condition"] = np.where(merged["exp_mean"] > merged["ctrl_mean"],
                                   experimental_condition, "Control")
    merged["predicted_DE"] = (merged["condition"] == experimental_condition).astype(int)
    return merged


# ---------------------------
# Orchestrator
# ---------------------------

def run_unsupervised_ml(predictions_file: str, norm_counts_file: str, output_dir: str, experimental_condition: str, conditions_list: Optional[List[str]] = None):

    """
    Steps:
      1) Read absolute predictions (index=gene) and normalized counts.
      2) Prefer mean_experimental/mean_control if present; else compute from normalized_counts.csv + --conditions.
      3) Assign 'condition' + 'predicted_DE'.
      4) Enrichment (Fisher).
      5) Replicate CV (log1p counts) across experimental replicates.
      6) PCA + KMeans on samples; color by original --conditions.
    """
    _safe_makedirs(output_dir)

    # 1) Read data
    dott_df = pd.read_csv(predictions_file, index_col=0)
    norm_counts_df = pd.read_csv(norm_counts_file, index_col=0)  # first col "gene" used as index

    # 2) Means
    have_precomputed = {"mean_experimental", "mean_control"}.issubset(dott_df.columns)
    if have_precomputed:
        means_df = dott_df[["mean_experimental", "mean_control"]].rename(
            columns={"mean_experimental": "exp_mean", "mean_control": "ctrl_mean"}
        )
    else:
        if conditions_list is None:
            raise ValueError("conditions_list is required to compute group means from normalized_counts.csv.")
        means_df = _compute_group_means_from_norm_counts(norm_counts_df, conditions_list, experimental_condition)

    # 3) Assign condition/predicted_DE
    dott_df = _assign_condition(dott_df, means_df, experimental_condition)
    updated_csv = os.path.join(output_dir, "dott_genes_with_conditions.csv")
    dott_df.to_csv(updated_csv)
    print("Updated predictions written:", updated_csv)

    # 4) Enrichment
    compare_enrichment_across_conditions(updated_csv, output_dir, experimental_condition)

    # 5) Experimental/control sample names from conditions_list
    conds = list(conditions_list) if conditions_list is not None else ["NA"] * norm_counts_df.shape[1]
    exp_samples = [s for s, c in zip(norm_counts_df.columns, conds) if c == experimental_condition]
    ctrl_samples = [s for s, c in zip(norm_counts_df.columns, conds) if c != experimental_condition]
    print("Experimental samples:", exp_samples)
    print("Control samples:", ctrl_samples)

    # Filter to experimental-condition predicted genes (for CV & PCA subset)
    dott_df_filtered = dott_df[dott_df["condition"] == experimental_condition]
    filtered_csv = os.path.join(output_dir, "dott_genes_with_conditions_filtered.csv")
    dott_df_filtered.to_csv(filtered_csv)

    predicted_genes = dott_df_filtered.index.tolist()

    # CV (experimental only)
    assess_replicate_consistency(norm_counts_df, predicted_genes, exp_samples, output_dir)

    # PCA + KMeans
    run_pca_and_kmeans(norm_counts_df, conds, output_dir, predicted_genes=predicted_genes)

    return updated_csv, filtered_csv
