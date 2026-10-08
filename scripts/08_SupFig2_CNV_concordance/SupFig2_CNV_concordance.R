#!/usr/bin/env Rscript

############################################################
# SupFig2_CNV_concordance.R
# ----------------------------------------------------------
# Supplementary Fig. 2: concordance of CopyKAT CNV calls with
# inferCNV and with matched WES (CNVkit), scRNA-seq (CELL) platform.
#
# Panels (stacked heatmaps; columns = 10-Mb bins, chr1-22):
#   (1) Single cell : per-cell CopyKAT CNA, aneuploid cells
#   (2) inferCNV    : pseudobulk of malignant cells
#   (3) CopyKAT     : pseudobulk = mean of (1)
#   (4) WES         : CNVkit bin-level log2 copy ratio
#   Pearson r of (2) vs (3) and (4) vs (3) per sample (reported in legend).
#
# Inputs (via --args):
#   --copykat_dir   : CopyKAT outputs per sample
#                     LUNG-{S}-{P}_copykat_CNA_results.txt
#                     LUNG-{S}-{P}_copykat_prediction.txt
#   --infercnv_dir  : inferCNV runs per sample
#                     LUNG-{S}-{P}/run.final.infercnv_obj
#   --callmal_dir   : malignant-cell calls from inferCNV
#                     inferCNV_SD_vs_COR_LUNG-{S}-{P}_classification.txt
#                     (columns Row.names, cell_index; "Malignant cells")
#   --cnvkit_dir    : CNVkit output (WES_CNVkit.sh), Lung_{S}.cnr
#   --samples       : comma-separated tumor IDs (default T35,T36,T39,T44)
#   --platform      : CELL (default) or NUC
#   --outdir        : output directory (default SupFig2_out)
#
# Outputs:
#   - SupFig2_CNV_concordance.pdf
#   - SupFig2_CNV_concordance_pearson_r.csv
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(infercnv)   # class definition for readRDS of infercnv objects
})

# -----------------------------
# 0) Args
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NA_character_) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

CK_DIR      <- get_arg("--copykat_dir")
ICNV_DIR    <- get_arg("--infercnv_dir")
CALLMAL_DIR <- get_arg("--callmal_dir")
WES_DIR     <- get_arg("--cnvkit_dir")
SAMPLES     <- strsplit(get_arg("--samples", default = "T35,T36,T39,T44"), ",")[[1]]
PLATFORM    <- get_arg("--platform", default = "CELL")
OUTDIR      <- get_arg("--outdir", default = "SupFig2_out")

for (k in c("CK_DIR", "ICNV_DIR", "CALLMAL_DIR", "WES_DIR")) {
  if (is.na(get(k))) stop("missing required argument for ", k)
}
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Parameters
# -----------------------------
WIN     <- 1e7                  # 10-Mb bins
CHROMS  <- paste0("chr", 1:22)
N_CELLS <- 300                  # per-sample downsampling, single-cell panel only
center  <- function(v) v - median(v, na.rm = TRUE)
set.seed(1)

# -----------------------------
# 2) Per-cell CopyKAT CNA (aneuploid cells), binned
# -----------------------------
read_ck_cells <- function(s) {
  f  <- file.path(CK_DIR, sprintf("LUNG-%s-%s_copykat_CNA_results.txt", s, PLATFORM))
  pf <- file.path(CK_DIR, sprintf("LUNG-%s-%s_copykat_prediction.txt", s, PLATFORM))
  pred <- fread(pf)
  aneu <- make.names(pred[copykat.pred == "aneuploid", cell.names])
  hdr  <- names(fread(f, nrows = 0)); aneu <- intersect(aneu, hdr)
  if (length(aneu) == 0) stop("no aneuploid cells for ", s)
  dt <- fread(f, select = c("chrom", "chrompos", aneu))
  dt <- dt[chrom %in% 1:22]
  dt[, chr := paste0("chr", chrom)][, winidx := floor(chrompos / WIN)]
  agg <- dt[, lapply(.SD, mean), by = .(chr, winidx), .SDcols = aneu]
  key <- paste(agg$chr, agg$winidx, sep = ":")
  m <- t(as.matrix(agg[, ..aneu])); colnames(m) <- key
  m
}
ckl <- setNames(lapply(SAMPLES, read_ck_cells), SAMPLES)

# -----------------------------
# 3) inferCNV pseudobulk (malignant cells), binned
# -----------------------------
read_icnv <- function(s) {
  obj  <- readRDS(file.path(ICNV_DIR, sprintf("LUNG-%s-%s", s, PLATFORM), "run.final.infercnv_obj"))
  expr <- obj@expr.data
  go   <- obj@gene_order
  obs  <- unlist(obj@observation_grouped_cell_indices)
  obs.names <- colnames(expr)[obs]
  cls <- read.table(file.path(CALLMAL_DIR,
           sprintf("inferCNV_SD_vs_COR_LUNG-%s-%s_classification.txt", s, PLATFORM)),
           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  mal <- cls$Row.names[cls$cell_index == "Malignant cells"]
  use <- intersect(obs.names, mal)
  if (length(use) < 10) { warning("few malignant cells for ", s, "; using all obs"); use <- obs.names }
  prof <- rowMeans(expr[, use, drop = FALSE]) - 1
  chr  <- as.character(go$chr)
  chr  <- ifelse(grepl("^chr", chr), chr, paste0("chr", chr))
  df   <- data.frame(chr = chr, start = go$start, val = prof)
  df   <- df[df$chr %in% CHROMS, ]
  df$winidx <- floor(df$start / WIN)
  agg  <- aggregate(val ~ chr + winidx, df, mean)
  setNames(agg$val, paste(agg$chr, agg$winidx, sep = ":"))
}
icnv <- setNames(lapply(SAMPLES, read_icnv), SAMPLES)

# -----------------------------
# 4) WES CNVkit log2 copy ratio (.cnr), binned
# -----------------------------
read_wes <- function(s) {
  d <- read.table(file.path(WES_DIR, paste0("Lung_", s, ".cnr")),
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  d <- d[d$chromosome %in% CHROMS, c("chromosome", "start", "log2")]
  d$winidx <- floor(d$start / WIN)
  agg <- aggregate(log2 ~ chromosome + winidx, data = d, FUN = mean)
  setNames(agg$log2, paste(agg$chromosome, agg$winidx, sep = ":"))
}
wes <- setNames(lapply(SAMPLES, read_wes), SAMPLES)

# -----------------------------
# 5) Common bins across all sources, genomic order
# -----------------------------
keys <- Reduce(intersect, c(lapply(ckl, colnames), lapply(icnv, names), lapply(wes, names)))
ord  <- data.frame(key = keys, stringsAsFactors = FALSE)
ord$chr <- factor(sub(":.*", "", ord$key), levels = CHROMS)
ord$win <- as.integer(sub(".*:", "", ord$key))
ord  <- ord[order(ord$chr, ord$win), ]
keys <- ord$key
chr_split <- ord$chr

# -----------------------------
# 6) Matrices
# -----------------------------
ckbulk <- lapply(ckl, function(m) colMeans(m[, keys, drop = FALSE]))   # CopyKAT pseudobulk

set.seed(1)
ckl.ds <- lapply(ckl, function(m) {
  idx <- if (nrow(m) > N_CELLS) sample(nrow(m), N_CELLS) else seq_len(nrow(m))
  m[idx, keys, drop = FALSE]
})
mat1 <- do.call(rbind, ckl.ds)
sample1 <- factor(rep(SAMPLES, times = sapply(ckl.ds, nrow)), levels = SAMPLES)

mat2 <- t(sapply(icnv,   function(v) center(v[keys]))); rownames(mat2) <- SAMPLES
mat3 <- t(sapply(ckbulk, function(v) center(v[keys]))); rownames(mat3) <- SAMPLES
mat4 <- t(sapply(wes,    function(v) center(v[keys]))); rownames(mat4) <- SAMPLES

# -----------------------------
# 7) Pearson r vs CopyKAT pseudobulk
# -----------------------------
r_ic_ck  <- sapply(SAMPLES, function(s) cor(mat2[s, ], mat3[s, ], use = "complete.obs"))
r_wes_ck <- sapply(SAMPLES, function(s) cor(mat4[s, ], mat3[s, ], use = "complete.obs"))

r_tab <- data.frame(sample = SAMPLES,
                    r_inferCNV_vs_CopyKAT = round(r_ic_ck, 3),
                    r_WES_vs_CopyKAT      = round(r_wes_ck, 3),
                    n_bins = length(keys))
print(r_tab, row.names = FALSE)
write.csv(r_tab, file.path(OUTDIR, "SupFig2_CNV_concordance_pearson_r.csv"), row.names = FALSE)

# -----------------------------
# 8) Clipping and colours (red = gain, blue = loss)
# -----------------------------
clip_at <- function(m, q) as.numeric(quantile(abs(m), q, na.rm = TRUE))
clip1 <- clip_at(mat1, 0.95)
clip2 <- clip_at(mat2, 0.90)
clip3 <- clip_at(mat3, 0.90)
clip4 <- clip_at(mat4, 0.90)
clamp <- function(m, cl) pmax(pmin(m, cl), -cl)
ramp  <- function(cl) colorRamp2(c(-cl, 0, cl), c("#2166AC", "white", "#B2182B"))

# -----------------------------
# 9) Heatmaps
# -----------------------------
ht1 <- Heatmap(clamp(mat1, clip1), name = "single_cell", col = ramp(clip1),
               height = unit(9, "cm"),
               row_split = sample1, cluster_row_slices = FALSE,
               cluster_rows = TRUE, show_row_dend = FALSE, show_row_names = FALSE,
               row_title_gp = gpar(fontsize = 9),
               cluster_columns = FALSE, show_column_names = FALSE,
               column_split = chr_split, column_gap = unit(0, "mm"),
               column_title_gp = gpar(fontsize = 7),
               use_raster = TRUE, border = TRUE,
               heatmap_legend_param = list(title = "CopyKAT CNA\n(single cell)"))

bulk_ht <- function(m, cl, name, title) {
  Heatmap(clamp(m, cl), name = name, col = ramp(cl), height = unit(2.6, "cm"),
          cluster_rows = FALSE, show_row_names = TRUE,
          row_names_side = "left", row_names_gp = gpar(fontsize = 9),
          cluster_columns = FALSE, show_column_names = FALSE,
          column_split = chr_split, column_gap = unit(0, "mm"), column_title = NULL,
          use_raster = TRUE, border = TRUE,
          heatmap_legend_param = list(title = title))
}
ht2 <- bulk_ht(mat2, clip2, "inferCNV",           "inferCNV")
ht3 <- bulk_ht(mat3, clip3, "CopyKAT_pseudobulk", "CopyKAT\npseudobulk")
ht4 <- bulk_ht(mat4, clip4, "wes",                "WES log2")

out.pdf <- file.path(OUTDIR, "SupFig2_CNV_concordance.pdf")
pdf(out.pdf, width = 13, height = 10.5)
draw(ht1 %v% ht2 %v% ht3 %v% ht4, ht_gap = unit(5, "mm"), merge_legends = TRUE)
dev.off()
cat("wrote:", out.pdf, "\n")
