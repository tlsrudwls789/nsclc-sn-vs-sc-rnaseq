# Comparative analysis of single-cell and single-nucleus RNA-seq in NSCLC

This repository contains custom R scripts used to generate the analyses and
figures presented in the manuscript submitted to *Communications Biology*.

The study compares single-cell RNA sequencing (scRNA-seq) and single-nucleus
RNA sequencing (snRNA-seq) data from non-small cell lung cancer (NSCLC) tissues,
with a focus on epithelial, immune, and stromal cell populations, transcriptional
programs, and cell–cell communication.

---

## Repository structure

Scripts are organized by manuscript results and figures. Each script corresponds
to a specific main or supplementary figure in the manuscript and can be executed
independently to reproduce the associated analyses and visualizations.

---

## Software requirements

Analyses were performed using R (>= 4.3.0) and the following key packages:

- Seurat (v4.4.0)
- monocle2
- SCENIC
- clusterProfiler
- tidyverse and related packages

Exact package versions and key analysis parameters are described in the Methods
section of the manuscript.

---

## Usage

Each script is designed to reproduce the analyses associated with a specific
result or figure.

Input and output directories can be specified via command-line arguments or
environment variables (e.g., `DATA_DIR`, `OUT_DIR`), as documented within each
script. Scripts automatically create output directories if they do not exist.

---

## Data availability

Raw sequencing data are deposited in the Gene Expression Omnibus (GEO) under
the accession number(s) listed in the manuscript. This repository contains 
the analysis code and processed data; raw sequencing files are not included 
and should be obtained from GEO.

---

## Code availability

All custom code used to generate the results and figures in this study is
provided in this repository.
