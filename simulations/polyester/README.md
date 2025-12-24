**Title:** Polyester simulation scripts (mm39)
**Purpose:** Regenerate the simulated mouse liver RNA-seq datasets used to train/evaluate DoTT-ML’s supervised module.

**Prerequisites**

-R + Bioconductor packages used in the scripts (polyester, rtracklayer, GenomicRanges, Biostrings, Rsamtools).

-Reference files placed locally in pipeline/:

  -pipeline/mm39.fa
  
  -pipeline/mm39_RefSeq.gtf

**Download references (example: UCSC)**
(Reference directories shown here are maintained by UCSC)

```
# from repo root
mkdir -p pipeline
cd pipeline

# genome fasta
wget -O mm39.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
gunzip -f mm39.fa.gz

# RefSeq GTF (UCSC “ncbiRefSeq”)
wget -O mm39.ncbiRefSeq.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.ncbiRefSeq.gtf.gz
gunzip -f mm39.ncbiRefSeq.gtf.gz
mv mm39.ncbiRefSeq.gtf mm39_RefSeq.gtf
```

**Run**

From repo root:

```
Rscript simulations/polyester/simulate_polyester_biologically_realistic_and_150bp_library_5_percent_DoTT.R
Rscript simulations/polyester/simulate_polyester_biologically_realistic_and_150bp_library_10_percent_DoTT.R
Rscript simulations/polyester/simulate_polyester_biologically_realistic_and_150bp_library_20_percent_DoTT.R
Rscript simulations/polyester/simulate_polyester_biologically_realistic_and_150bp_library_37_percent_DoTT.R
```

**Outputs**

- Generated FASTA intermediates and simulated read folders will be written under pipeline/ (see each script’s outdir = "pipeline/simulated_reads_biologically_realistic_*").
