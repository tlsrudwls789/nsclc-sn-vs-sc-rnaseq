#!/usr/bin/env Rscript

############################################################
# SupFig8_patient_trajectories.R
# ----------------------------------------------------------
# Supplementary Fig. 8: the AT2-to-malignant transition is
# reproduced in individual patients (snRNA-seq).
#
# Panels:
#   a : Monocle2 trajectory of AT2 + malignant epithelial cells, patient T35
#   b : Monocle2 trajectory of AT2 + malignant epithelial cells, patient T44
#   In each panel: left = colored by identity (orange AT2, red malignant),
#                  right = colored by Pseudotime.
#
# Inputs (via --key=value):
#   --epi_rds       : snRNA-seq epithelial Seurat object (CopyKAT-annotated)
#   --samples       : patient IDs, comma-separated (default: T35,T44)
#   --outdir        : output directory (default: SupFig8_out)
#   --sample_col    : sample column (default: sample)
#   --celltype_col  : identity column (default: celltype3)
#   --sample_prefix : sample name prefix (default: LUNG-)
#   --sample_suffix : sample name suffix (default: -NUC)
#   --cores         : cores for reduceDimension (default: 32)
#   --monocle_patch : (optional) path to a patched Monocle2 order_cells.R
#                     (needed only if monocle::orderCells fails on R >= 4)
#
# Required meta columns in epi_rds:
#   sample    : e.g. LUNG-N35-NUC / LUNG-T35-NUC (normal + tumor of one patient)
#   celltype3 : epithelial identity; "AT2" (AT2 + Dividing_AT2) and
#               "Malignant cells" (CopyKAT aneuploid cells)
#   RNA assay with raw counts.
#
# Outputs (in --outdir):
#   SupFig8a_T35_trajectory.pdf, SupFig8b_T44_trajectory.pdf
#   SupFig8_<ID>_pseudotime.tsv (per-cell Pseudotime, State, identity)
############################################################

Sys.setenv(R_LIFECYCLE_VERBOSITY = "quiet")
options(lifecycle_verbosity = "quiet")
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(monocle)   # Monocle2
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

# -----------------------------
# 0) Args
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NA_character_) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

EPI_RDS       <- get_arg("--epi_rds")
SAMPLES       <- strsplit(get_arg("--samples", default = "T35,T44"), ",")[[1]]
OUT_DIR       <- get_arg("--outdir", default = "SupFig8_out")
SAMPLE_COL    <- get_arg("--sample_col", default = "sample")
CELLTYPE_COL  <- get_arg("--celltype_col", default = "celltype3")
SAMPLE_PREFIX <- get_arg("--sample_prefix", default = "LUNG-")
SAMPLE_SUFFIX <- get_arg("--sample_suffix", default = "-NUC")
N_CORES       <- as.integer(get_arg("--cores", default = "32"))
MONOCLE_PATCH <- get_arg("--monocle_patch")

if (is.na(EPI_RDS) || !file.exists(EPI_RDS)) stop("Provide --epi_rds=<epithelial snRNA Seurat RDS>")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
if (!is.na(MONOCLE_PATCH)) source(MONOCLE_PATCH)  # overrides orderCells()/project2MST()

# -----------------------------
# 1) Parameters (identical to the analysis scripts)
# -----------------------------
AT2_LABEL       <- "AT2"
MAL_LABEL       <- "Malignant cells"
ORDER_MEAN_EXPR <- 0.001   # ordering genes: dispersion_empirical > fit & mean > 0.001
MAX_COMPONENTS  <- 3       # DDRTree dimensions
REVERSE         <- TRUE    # orderCells(reverse = TRUE): root at the AT2 end
CELLTYPE_COLORS <- c("Malignant cells" = "#D53E4F",  # red
                     "AT2"             = "#FDAE61")  # orange

theme_traj <- theme(
  axis.title = element_blank(), axis.text = element_blank(),
  axis.line = element_blank(), axis.ticks = element_blank(),
  panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), panel.background = element_blank())

# -----------------------------
# 2) Load epithelial object
# -----------------------------
epi <- readRDS(EPI_RDS)
DefaultAssay(epi) <- "RNA"
stopifnot(all(c(SAMPLE_COL, CELLTYPE_COL) %in% colnames(epi@meta.data)))

# -----------------------------
# 3) Monocle2 per patient
# -----------------------------
# Normal + tumor epithelial cells of one patient -> AT2 + malignant -> DDRTree
run_patient_monocle <- function(seu, patient_id) {
  pid <- sub("^[A-Za-z]+", "", patient_id)                 # "T35" -> "35"
  smp <- paste0(SAMPLE_PREFIX, c("N", "T"), pid, SAMPLE_SUFFIX)
  seu <- seu[, seu@meta.data[[SAMPLE_COL]] %in% smp]
  message(sprintf("[%s] %d epithelial cells (%s)", patient_id, ncol(seu), paste(smp, collapse = ", ")))

  # counts (all genes) + metadata
  exprs <- GetAssayData(seu, assay = "RNA", layer = "counts")
  exprs <- as(exprs, "dgCMatrix")
  pd <- seu@meta.data[colnames(exprs), , drop = FALSE]
  pd$Cell_Name <- rownames(pd)
  pd$celltype  <- trimws(as.character(pd[[CELLTYPE_COL]]))

  # subset AT2 + malignant cells
  keep <- pd$celltype %in% c(AT2_LABEL, MAL_LABEL)
  exprs <- exprs[, keep]
  pd <- pd[keep, , drop = FALSE]
  pd$celltype <- factor(pd$celltype, levels = c(MAL_LABEL, AT2_LABEL))
  print(table(pd$celltype, pd[[SAMPLE_COL]]))

  gene_ann <- data.frame(gene_short_name = rownames(exprs), row.names = rownames(exprs))
  cds <- monocle::newCellDataSet(
    cellData = exprs,
    phenoData = new("AnnotatedDataFrame", data = pd),
    featureData = new("AnnotatedDataFrame", data = gene_ann),
    expressionFamily = negbinomial.size(),
    lowerDetectionLimit = 0)

  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)

  # ordering genes by dispersion
  disp <- dispersionTable(cds)
  ordering_genes <- subset(disp, dispersion_empirical > dispersion_fit & mean_expression > ORDER_MEAN_EXPR)
  message(sprintf("[%s] %d ordering genes", patient_id, nrow(ordering_genes)))
  cds <- monocle::setOrderingFilter(cds, ordering_genes$gene_id)

  cds <- monocle::reduceDimension(
    cds,
    max_components = MAX_COMPONENTS,
    reduction_method = "DDRTree",
    verbose = TRUE,
    cores = N_CORES,
    norm_method = "vstExprs",
    relative_expr = TRUE,
    scaling = TRUE)

  cds <- orderCells(cds, reverse = REVERSE)

  # sanity check: pseudotime should increase from AT2 to malignant cells
  mpt <- tapply(pData(cds)$Pseudotime, pData(cds)$celltype, median)
  message(sprintf("[%s] median pseudotime: %s", patient_id,
                  paste(names(mpt), round(mpt, 2), sep = "=", collapse = ", ")))
  if (mpt[AT2_LABEL] > mpt[MAL_LABEL]) warning(sprintf("[%s] AT2 is not at the trajectory start", patient_id))

  cds
}

# -----------------------------
# 4) Plots: identity (left) + Pseudotime (right)
# -----------------------------
plot_patient_trajectory <- function(cds) {
  p_id <- plot_cell_trajectory(cds, color_by = "celltype", cell_size = 2,
                               show_tree = TRUE, show_branch_points = FALSE, theta = 0) +
    scale_color_manual(name = NULL, values = CELLTYPE_COLORS) +
    theme_traj + theme(legend.position = "bottom")

  p_pt <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.1) +
    theme_traj + theme(legend.position = "right")

  p_id | p_pt
}

# -----------------------------
# 5) Run
# -----------------------------
for (i in seq_along(SAMPLES)) {
  sid   <- SAMPLES[i]
  panel <- letters[i]
  cds   <- run_patient_monocle(epi, sid)

  out_pdf <- file.path(OUT_DIR, sprintf("SupFig8%s_%s_trajectory.pdf", panel, sid))
  pdf(out_pdf, width = 8.5, height = 3.4)
  print(plot_patient_trajectory(cds))
  dev.off()

  pt <- pData(cds)[, c("Cell_Name", SAMPLE_COL, "celltype", "Pseudotime", "State")]
  write.table(pt, file.path(OUT_DIR, sprintf("SupFig8_%s_pseudotime.tsv", sid)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved: ", out_pdf)
}

message("Done. Outputs in: ", normalizePath(OUT_DIR))
