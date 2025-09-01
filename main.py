#!/usr/bin/env python3
import os
import argparse
import subprocess
import numpy as np

from annotation import generate_saf
from featureCounts import run_featurecounts, filter_counts_for_deseq2
from merging_genes import merge_gene_regions
from coordinate_extraction import extract_coordinates
from utils import parse_conditions

def main():
    parser = argparse.ArgumentParser(description="DoTT Bioinformatics Pipeline")
    parser.add_argument("--gtf-file", required=True)
    parser.add_argument("--bam-files", nargs="+", required=True)
    parser.add_argument("--species", choices=["mm39","hg38","hg19"], required=True)
    parser.add_argument("--id-type", choices=["Ensembl_ID","Symbol"], default="Symbol")
    parser.add_argument("--extension", type=int, default=10000)
    parser.add_argument("--dynamic", action="store_true")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--conditions", required=True)
    parser.add_argument("--kgx-file", default="kgXref.txt.gz")
    parser.add_argument("--run_gsea", action="store_true")
    parser.add_argument("--unsupervised_ml", action="store_true")
    parser.add_argument("--supervised_ml", action="store_true")
    parser.add_argument("--experimental_condition")
    parser.add_argument("--ground_truth")
    parser.add_argument("--bootstrap", action="store_true")
    parser.add_argument("--n_boot", type=int, default=100)
    parser.add_argument("--consensus_threshold", type=float, default=0.5)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse conditions as a list aligned to BAM files order
    conds = parse_conditions(args.bam_files, args.conditions)

    saf_file = generate_saf(
        gtf_file=args.gtf_file,
        extension=args.extension,
        output_dir=args.output_dir,
        dynamic=args.dynamic,
        bam_file_for_dynamic=(args.bam_files[0] if args.dynamic else None),
        species=args.species,
        id_type=args.id_type,
        kgx_file=args.kgx_file
    )

    counts_file = run_featurecounts(saf_file, args.extension, args.bam_files, args.output_dir)
    cleaned_counts = filter_counts_for_deseq2(counts_file, args.output_dir)
    
    # Run DESeq2 (optionally bootstrap consensus)
    r_cmd = ["Rscript", "dott_DESeq2.R",
             "--counts_file", cleaned_counts,
             "--output_dir", args.output_dir,
             "--conditions", args.conditions]
    if args.bootstrap:
        r_cmd += ["--bootstrap", "TRUE",
                  "--n_boot", str(args.n_boot),
                  "--consensus_threshold", str(args.consensus_threshold)]
    subprocess.run(r_cmd, check=True)
    
    # Expected outputs from DESeq2 stage
    deseq2_results_file = os.path.join(args.output_dir, "3_UTR_extended_differential_analysis_results.csv")
    absolute_sig_file   = os.path.join(args.output_dir, "absolute_significant_extended_genes_with_individual_means.csv")
    norm_counts_file    = os.path.join(args.output_dir, "normalized_counts.csv")
    
    # Downstream DoTT annotation products
    merge_gene_regions(absolute_sig_file, saf_file, os.path.join(args.output_dir, "merged_DoTT_annotation.bed"))
    extract_coordinates(absolute_sig_file, saf_file, args.output_dir)
    
    # Optional GSEA
    if args.run_gsea:
        subprocess.run([
            "Rscript", "dott_gsea_preranked_list.R",
            "--input_file", absolute_sig_file,
            "--output_dir", args.output_dir
        ], check=True)
        print("GSEA pre-ranked list created.")
    
    # Optional UNSUPERVISED ML
    if args.unsupervised_ml:
        if not args.experimental_condition:
            raise ValueError("For unsupervised ML, provide --experimental_condition.")
        from Unsupervised_ML import run_unsupervised_ml
        # Pass the parsed conditions list so the module can compute means from normalized_counts
        run_unsupervised_ml(absolute_sig_file, norm_counts_file,
                            args.output_dir, args.experimental_condition,
                            conditions_list=conds)
    
    # Optional SUPERVISED ML
    if args.supervised_ml:
        if not args.experimental_condition or not args.ground_truth:
            raise ValueError("For supervised ML, provide --experimental_condition and --ground_truth.")
        from Supervised_ML import (
            compare_to_ground_truth,
            train_ml_classifier_cv,
            evaluate_ml_performance
        )
        gt_merged, gt_cm, gt_auc = compare_to_ground_truth(
            args.ground_truth, deseq2_results_file,
            args.gtf_file, args.output_dir
        )
        print(f"DESeq2 full ROC AUC: {gt_auc:.3f}")
        print("See DESeq2_cv_performance_metrics.csv for CV results.")

        cv_scores, final_model, pred_df, holdout_results, cv_results = \
            train_ml_classifier_cv(
                args.ground_truth, deseq2_results_file,
                args.gtf_file, args.output_dir,
                clf_cutoff=0.5, cv=10, test_size=0.2
            )
        print("RF+SMOTE CV AUCs:", cv_scores)
        print(f"Mean RF+SMOTE CV AUC: {np.mean(cv_scores):.3f}")

        merged_ml, merged_holdout, merged_cv, cm_df, training_metrics_df = \
            evaluate_ml_performance(
                pred_df, holdout_results,
                args.ground_truth, args.gtf_file,
                deseq2_results_file, args.output_dir,
                cv_results=cv_results
            )
        print(f"RF training ROC AUC: {training_metrics_df['ROC_AUC'].iloc[0]:.3f}")

    print("Pipeline finished successfully.")


if __name__ == "__main__":
    main()
