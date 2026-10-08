<h1 align="center">Complementary single-nucleus and single-cell RNA sequencing reveal tumor state evolution and microenvironmental rewiring in NSCLC</h1>

<p align="center">
  <a href="https://doi.org/10.1038/s42003-026-11040-y"><img src="https://img.shields.io/badge/Communications%20Biology-2026-0b5394" alt="Communications Biology 2026"></a>
  <a href="https://doi.org/10.1038/s42003-026-11040-y"><img src="https://img.shields.io/badge/DOI-10.1038%2Fs42003--026--11040--y-blue" alt="DOI"></a>
  <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310916"><img src="https://img.shields.io/badge/GEO-GSE310916-2e7d32" alt="GEO GSE310916"></a>
  <img src="https://img.shields.io/badge/R-%E2%89%A54.3.0-276DC3?logo=r&logoColor=white" alt="R >= 4.3.0">
  <img src="https://img.shields.io/badge/Seurat-v4.4.0-6f42c1" alt="Seurat v4.4.0">
</p>

<p align="center">
  Analysis code for <b>Shin, G.-J. <i>et al.</i></b>, <i>Communications Biology</i> (2026), published online 2 October 2026.
</p>

---

## Overview

This repository contains the R scripts that generate the main and supplementary figures of the article.
The study compares **single-cell RNA-seq (scRNA-seq)** and **single-nucleus RNA-seq (snRNA-seq)** of
non-small cell lung cancer (NSCLC) tissues from an in-house cohort (SMC) and a public cohort (UL; Renaut *et al.*),
covering cell composition, platform-associated transcriptional differences, alveolar differentiation,
the AT2-to-malignant transition and ligand–receptor rewiring.

## Figures and scripts

Each folder contains one figure-level script that can be run on its own.

**Main figures**

| Figure | Content | Script |
|---|---|---|
| Fig. 1 | Integrated SMC atlas: UMAPs, cell type frequencies, paired sn vs sc proportions | [`01_Fig1_SupFig1cd`](scripts/01_Fig1_SupFig1cd) |
| Fig. 2 | Epithelial, endothelial and fibroblast subtype proportions (SMC, UL) | [`02_Fig2`](scripts/02_Fig2) |
| Fig. 3 | Platform-biased DEGs and GO terms in five shared subtypes | [`03_Fig3_SupFig4`](scripts/03_Fig3_SupFig4) |
| Fig. 4 | AT2-to-AT1 trajectories (Monocle2) and TF regulon activity (SCENIC) | [`04_Fig4`](scripts/04_Fig4) |
| Fig. 5 | AT2-to-malignant pseudotime modules M1–M4 (snRNA-seq) | [`05_Fig5_SupFig6`](scripts/05_Fig5_SupFig6) |
| Fig. 6 | Ligand–receptor interactions across epithelial Early/Late states (LIANA) | [`06_Fig6`](scripts/06_Fig6) |

**Supplementary figures**

| Figure | Content | Script |
|---|---|---|
| Supp. Fig. 1a–b | DecontX contamination and post-DecontX QC metrics | [`07_SupFig1ab_QC`](scripts/07_SupFig1ab_QC) |
| Supp. Fig. 1c–d | Lineage marker FeaturePlots and sub-cell type UMAPs | [`01_Fig1_SupFig1cd`](scripts/01_Fig1_SupFig1cd) |
| Supp. Fig. 2 | CopyKAT vs inferCNV vs WES (CNVkit) CNV concordance | [`08_SupFig2_CNV_concordance`](scripts/08_SupFig2_CNV_concordance) |
| Supp. Fig. 3 | Public (UL) dataset: UMAPs, CopyKAT, cell type proportions | [`09_SupFig3_public_dataset`](scripts/09_SupFig3_public_dataset) |
| Supp. Fig. 4 | Jaccard clustering of enriched GO terms | [`03_Fig3_SupFig4`](scripts/03_Fig3_SupFig4) |
| Supp. Fig. 5 | Lung cancer cell-line validation (GSE303762) <sup>†</sup> | [`10_SupFig5_cell_lines`](scripts/10_SupFig5_cell_lines) |
| Supp. Fig. 6 | Robustness of the AT2-to-malignant program (core DEGs) | [`05_Fig5_SupFig6`](scripts/05_Fig5_SupFig6) |
| Supp. Fig. 7 | Early/Late state assignment (Delta = M4 − M1) and transfer to scRNA-seq | [`11_SupFig7_early_late`](scripts/11_SupFig7_early_late) |
| Supp. Fig. 8 | Patient-level Monocle2 trajectories (T35, T44) | [`12_SupFig8_patient_trajectories`](scripts/12_SupFig8_patient_trajectories) |

<sup>†</sup> The Supp. Fig. 5 script was re-written after publication following the Methods of the article; the original
figure-generation script was not retained. Parameters not stated in the Methods are marked `assumed` in the script.

## Requirements

- **R** ≥ 4.3.0; **Seurat** v4.4.0 (Supp. Fig. 5: Seurat v5)
- **monocle** (v2), **SCENIC**, **liana**, **clusterProfiler**, **harmony**
- **ComplexHeatmap**, **ggplot2**, **patchwork**, **cowplot** and the tidyverse
- Supp. Fig. 2: **CopyKAT**, **inferCNV** and **CNVkit** (WES, shell script)

Exact package versions and analysis parameters are given in the Methods of the article.

## Usage

Inputs and outputs are passed as `--key=value` arguments, or through the `DATA_DIR` and `OUT_DIR`
environment variables, as documented in the header of each script. Output folders are created automatically.

```bash
# Figure 1 and Supplementary Figure 1c-d
DATA_DIR=/path/to/seurat_objects OUT_DIR=results/Fig1 \
  Rscript scripts/01_Fig1_SupFig1cd/Fig1_SupFig1cd.R

# Supplementary Figure 7 (also writes the Early/Late labels used by Figure 6)
Rscript scripts/11_SupFig7_early_late/SupFig7_early_late_states.R \
  --nuc_rds=nuc_epithelial.rds --cell_rds=cell_epithelial.rds \
  --module_genes=module_genes_core6.rds --outdir=results/SupFig7
```

Run order for the AT2-to-malignant analyses: Fig. 5 → Supp. Fig. 7 (module scores, Early/Late labels) → Fig. 6.

<details>
<summary><b>Repository layout</b></summary>

```text
scripts/
├── 01_Fig1_SupFig1cd/
├── 02_Fig2/
├── 03_Fig3_SupFig4/
├── 04_Fig4/
├── 05_Fig5_SupFig6/
├── 06_Fig6/
├── 07_SupFig1ab_QC/
├── 08_SupFig2_CNV_concordance/      (R script + WES_CNVkit.sh)
├── 09_SupFig3_public_dataset/
├── 10_SupFig5_cell_lines/
├── 11_SupFig7_early_late/
└── 12_SupFig8_patient_trajectories/
```
</details>

## Data availability

The scRNA-seq and snRNA-seq data generated in this study are available from GEO under
[**GSE310916**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310916).
The paired public lung dataset (Renaut *et al.*) is available from Zenodo
([10.5281/zenodo.11205626](https://doi.org/10.5281/zenodo.11205626)) and the cell-line dataset from GEO
([GSE303762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303762)).
This repository contains analysis code only; raw and processed sequencing data are not included.
Other external datasets are listed in Supplementary Data 1 of the [article](https://doi.org/10.1038/s42003-026-11040-y).

## Code availability

All custom code used to generate the figures of this study is provided in this repository.

## Citation

If you use this code, please cite:

> Shin, G.-J., Kang, H., Jo, A. *et al.* Complementary single-nucleus and single-cell RNA sequencing
> reveal tumor state evolution and microenvironmental rewiring in NSCLC.
> *Communications Biology* (2026). https://doi.org/10.1038/s42003-026-11040-y

```bibtex
@article{Shin2026NSCLC,
  author  = {Shin, Gyeong-Jin and Kang, Huiram and Jo, Areum and others},
  title   = {Complementary single-nucleus and single-cell RNA sequencing reveal tumor state evolution and microenvironmental rewiring in NSCLC},
  journal = {Communications Biology},
  year    = {2026},
  doi     = {10.1038/s42003-026-11040-y}
}
```
