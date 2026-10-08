#!/usr/bin/env Rscript

############################################################
# SupFig7_early_late_states.R
# ----------------------------------------------------------
# Supplementary Fig. 7: Early-to-late state transfer and
# cross-platform validation (snRNA-seq -> scRNA-seq).
#
# Panels:
#   a-c  snRNA-seq UMAPs: Pseudotime, cell type-tissue state, Early/Late
#   d    snRNA-seq Delta (M4 - M1) density by state, Youden cutoff (dashed)
#   e    snRNA-seq Early/Late composition by state
#   f    scRNA-seq sample-wise AT2 / malignant counts and proportions
#   g-h  scRNA-seq UMAPs: cell type-tissue state, Delta categories
#   i    scRNA-seq Delta density by state (dashed line at --cell_vline)
#   j    scRNA-seq Early/Mid/Late composition by state
#
# Method:
#   - Module gene sets M1-M4 (core6) are taken from the Fig. 5 script
#     output (module_genes_core6_<salt>.rds); they are NOT recomputed.
#   - Module score = mean over genes of per-gene z-scores (RNA "data" layer).
#   - Delta = Score_M4 - Score_M1.
#   - snRNA-seq: Youden-optimal cutoff from AT2 vs Malignant cells
#     (Malignant = positive, predict positive if Delta >= t);
#     Stage_EL = Early (Delta < cutoff) / Late (Delta >= cutoff).
#   - scRNA-seq: same module genes; quantiles of Delta in AT2 + Malignant
#     cells: <= 20% early_like, >= 80% late_like, else mid_like (Stage_like).
#
# Inputs (via --args):
#   --nuc_rds      : snRNA-seq AT2/malignant Seurat (Fig. 5 output
#                    seu_with_ModuleScores_<salt>.rds; Score_M1..M4 used if
#                    present, otherwise computed from --module_genes)
#   --cell_rds     : scRNA-seq AT2/malignant Seurat (seu_AT2Malig_umap_<salt>.rds)
#   --module_genes : named list RDS (M1..M4 -> SYMBOL), from Fig. 5
#   --pseudotime_csv : (optional) CSV with columns cell, Pseudotime for NUC
#                      cells, used when seu has no Pseudotime column
#   --nuc_umap     : NUC reduction name (default umap_h13)
#   --cell_umap    : CELL reduction name (default umap_h14)
#   --cell_vline   : x of dashed line in panel i (default 0)
#   --outdir       : output directory (default SupFig7_out)
#
# Required meta columns (both objects):
#   Sample, celltype4 (or celltype4_clean: AT2 / Malignant),
#   celltype_TissueDiseaseState (AT2_N, AT2_T, Malignant_N, Malignant_T)
#   NUC also: Pseudotime (or --pseudotime_csv)
#
# Outputs (outdir):
#   SupFig7a-c_NUC_UMAP.pdf, SupFig7d_NUC_delta_byTDS.pdf,
#   SupFig7e_NUC_EarlyLate_byTDS.pdf, SupFig7f_CELL_sample_AT2_Malignant.pdf,
#   SupFig7g-h_CELL_UMAP.pdf, SupFig7i_CELL_delta_byTDS.pdf,
#   SupFig7j_CELL_Stage_byTDS.pdf,
#   cutoffs_NUC_Youden.csv, cutoffs_CELL_quantile.csv,
#   early_late_labels_NUC.csv, early_late_labels_CELL.csv,
#   seu_NUC_withStage.rds, seu_CELL_withStage.rds  (inputs of Fig. 6)
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(Matrix)
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

NUC_RDS        <- get_arg("--nuc_rds")
CELL_RDS       <- get_arg("--cell_rds")
MODULE_GENES   <- get_arg("--module_genes")
PSEUDOTIME_CSV <- get_arg("--pseudotime_csv")
NUC_UMAP       <- get_arg("--nuc_umap",  default = "umap_h13")
CELL_UMAP      <- get_arg("--cell_umap", default = "umap_h14")
CELL_VLINE     <- as.numeric(get_arg("--cell_vline", default = "0"))
OUTDIR         <- get_arg("--outdir", default = "SupFig7_out")

for (f in c(NUC_RDS, CELL_RDS, MODULE_GENES)) {
  if (is.na(f) || !file.exists(f)) stop("Missing input: ", f)
}
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) Settings (identical to source)
# -----------------------------
SAMPLE_COL   <- "Sample"
CELLTYPE_COL <- "celltype4"
CT_CLEAN_COL <- "celltype4_clean"
TDS_COL      <- "celltype_TissueDiseaseState"
DELTA_COL    <- "Delta_M4_M1"
Q_EARLY      <- 0.20
Q_LATE       <- 0.80
BASE_SIZE    <- 13

modules    <- paste0("M", 1:4)
score_cols <- paste0("Score_", modules)
tds_levels <- c("AT2_N", "AT2_T", "Malignant_N", "Malignant_T")

epi.cell.colors2     <- c("Malignant" = "#D53E4F", "AT2" = "#FDAE61")
celltype_TDS.colors2 <- c("AT2_N" = "#FEE08B", "AT2_T" = "#F46D43",
                          "Malignant_N" = "#D53E4F", "Malignant_T" = "#9E0142")
stage_cols_NUC  <- c("Early" = "#E69F00", "Late" = "#88419D")
stage_cols_CELL <- c("early_like" = "#E69F00", "mid_like" = "#BDBDBD", "late_like" = "#88419D")

# -----------------------------
# 2) Helpers
# -----------------------------
pdf_safe <- function(filename, expr, width = 7, height = 5) {
  grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

theme_plain <- function(base = BASE_SIZE) {
  theme_classic(base_size = base) +
    theme(axis.title    = element_text(color = "black", face = "plain", size = base + 1),
          axis.text     = element_text(color = "black", face = "plain", size = base),
          legend.title  = element_text(color = "black", face = "plain", size = base + 1),
          legend.text   = element_text(color = "black", face = "plain", size = base),
          plot.title    = element_text(color = "black", face = "plain", size = base + 2),
          plot.subtitle = element_text(color = "black", face = "plain", size = base))
}

# AT2 / Malignant label ("Malignant cells" -> "Malignant")
ensure_ct_clean <- function(seu) {
  if (!CT_CLEAN_COL %in% colnames(seu@meta.data)) {
    seu@meta.data[[CT_CLEAN_COL]] <-
      dplyr::recode(as.character(seu@meta.data[[CELLTYPE_COL]]), "Malignant cells" = "Malignant")
  }
  seu
}

clean_tds <- function(x) {
  x <- trimws(as.character(x))
  x <- dplyr::recode(x, "Malignantcells_N" = "Malignant_N", "Malignantcells_T" = "Malignant_T",
                     "Malignant cells_N" = "Malignant_N", "Malignant cells_T" = "Malignant_T")
  factor(x, levels = tds_levels)
}

# z-mean module score: per-gene z across cells, averaged over module genes.
# sd_mode = "drop"  : genes with sd == 0 removed (Fig. 5 / NUC scoring)
# sd_mode = "floor" : sd floored at 1e-6 (CELL scoring in source)
score_zmean <- function(seu, module_genes, sd_mode = c("drop", "floor"),
                        min_genes = 10, sd_floor = 1e-6) {
  sd_mode <- match.arg(sd_mode)
  X <- tryCatch(GetAssayData(seu, assay = "RNA", layer = "data"), error = function(e) NULL)
  if (is.null(X) || nrow(X) == 0) X <- GetAssayData(seu, assay = "RNA", layer = "counts")
  if (!inherits(X, "dgCMatrix")) X <- as(X, "dgCMatrix")

  gs   <- lapply(module_genes[modules], function(g) intersect(unique(g[!is.na(g) & nzchar(g)]), rownames(X)))
  Xsub <- X[unique(unlist(gs)), , drop = FALSE]
  mu   <- Matrix::rowMeans(Xsub)
  sdv  <- sqrt(pmax(0, Matrix::rowMeans(Xsub^2) - mu^2))
  if (sd_mode == "drop") sdv[!is.finite(sdv) | sdv == 0] <- NA_real_
  else                   sdv[!is.finite(sdv) | sdv < sd_floor] <- sd_floor

  out <- sapply(modules, function(m) {
    g <- gs[[m]]
    g <- g[is.finite(mu[g]) & is.finite(sdv[g])]
    message(sprintf("[SCORE] %s usable genes: %d", m, length(g)))
    if (length(g) < min_genes) return(rep(NA_real_, ncol(X)))
    Z <- (as.matrix(Xsub[g, , drop = FALSE]) - mu[g]) / sdv[g]
    colMeans(Z)
  })
  colnames(out) <- score_cols
  rownames(out) <- colnames(X)
  as.data.frame(out)
}

# Youden-optimal threshold (Malignant = positive; positive if Delta >= t)
youden_cutoff <- function(x, y) {
  ok <- is.finite(x) & !is.na(y)
  x <- x[ok]; pos <- (y[ok] == "Malignant")
  thr <- sort(unique(x))
  if (length(thr) > 4000) thr <- thr[seq(1, length(thr), length.out = 4000)]
  st <- do.call(rbind, lapply(thr, function(t) {
    pred <- x >= t
    tp <- sum(pred & pos); fp <- sum(pred & !pos)
    tn <- sum(!pred & !pos); fn <- sum(!pred & pos)
    tpr <- ifelse(tp + fn > 0, tp / (tp + fn), NA_real_)
    fpr <- ifelse(fp + tn > 0, fp / (fp + tn), NA_real_)
    data.frame(thr = t, tpr = tpr, fpr = fpr, youden = tpr - fpr)
  }))
  st %>% filter(is.finite(youden)) %>% slice_max(youden, n = 1, with_ties = FALSE)
}

# Delta density by TDS (legend swatches forced via invisible geom_tile)
plot_delta_tds <- function(df, vline, title) {
  df_legend <- data.frame(x = median(df$delta, na.rm = TRUE), y = 0,
                          grp = factor(tds_levels, levels = tds_levels))
  ggplot(df, aes(x = delta, fill = tds)) +
    geom_density(alpha = 0.40, color = NA) +
    geom_tile(data = df_legend, inherit.aes = FALSE, aes(x = x, y = y, fill = grp),
              width = 0.02, height = 0.02, alpha = 0) +
    scale_fill_manual(values = celltype_TDS.colors2, breaks = tds_levels, drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(alpha = 1, color = "black"),
                               keywidth = grid::unit(0.55, "cm"), keyheight = grid::unit(0.55, "cm"))) +
    geom_vline(xintercept = vline, linetype = "dashed") +
    theme_plain() +
    labs(title = title, x = "Delta score (M4-M1)", y = "Density",
         fill = "celltype_TissueDiseaseState")
}

# Stage composition (% within TDS)
plot_comp_tds <- function(df, stage_col, cols, title, subtitle) {
  d <- df %>% count(tds, stage = .data[[stage_col]], name = "n") %>%
    group_by(tds) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
  ggplot(d, aes(x = tds, y = pct, fill = stage)) +
    geom_col(color = "black", width = 0.75) +
    scale_fill_manual(values = cols, drop = FALSE) +
    theme_plain(13) +
    theme(axis.title  = element_text(color = "black", face = "plain", size = 13),
          axis.text.x = element_text(color = "black", face = "plain", size = 11),
          axis.text.y = element_text(color = "black", face = "plain", size = 12),
          legend.title = element_text(color = "black", face = "plain", size = 13),
          plot.title  = element_text(color = "black", face = "plain", size = 13),
          plot.subtitle = element_text(color = "black", face = "plain", size = 12)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "% within group", fill = "Stage")
}

# Label table for downstream use (Fig. 6)
write_labels <- function(seu, stage_col, file) {
  seu@meta.data %>% rownames_to_column("cell") %>%
    dplyr::select(cell, any_of(c(SAMPLE_COL, CELLTYPE_COL, CT_CLEAN_COL, TDS_COL, "Pseudotime")),
                  all_of(score_cols), all_of(DELTA_COL), all_of(stage_col)) %>%
    write.csv(file, row.names = FALSE)
  message("[SAVED] ", file)
}

# -----------------------------
# 3) Load
# -----------------------------
module_genes <- readRDS(MODULE_GENES)
stopifnot(is.list(module_genes), all(modules %in% names(module_genes)))

seu_nuc  <- ensure_ct_clean(readRDS(NUC_RDS))
seu_cell <- ensure_ct_clean(readRDS(CELL_RDS))
message("[INFO] NUC cells: ", ncol(seu_nuc), " | CELL cells: ", ncol(seu_cell))

if (!"Pseudotime" %in% colnames(seu_nuc@meta.data) && !is.na(PSEUDOTIME_CSV)) {
  pt <- read.csv(PSEUDOTIME_CSV, stringsAsFactors = FALSE)
  seu_nuc$Pseudotime <- pt$Pseudotime[match(colnames(seu_nuc), pt$cell)]
}

# =========================================================
# 4) snRNA-seq (NUC): scores, Delta, Youden cutoff, Stage_EL
# =========================================================
if (!all(score_cols %in% colnames(seu_nuc@meta.data))) {
  sc <- score_zmean(seu_nuc, module_genes, sd_mode = "drop")
  for (cc in score_cols) seu_nuc@meta.data[[cc]] <- sc[[cc]]
}
seu_nuc@meta.data[[DELTA_COL]] <- seu_nuc$Score_M4 - seu_nuc$Score_M1

md_nuc <- seu_nuc@meta.data
is_epi <- as.character(md_nuc[[CT_CLEAN_COL]]) %in% c("AT2", "Malignant")
best   <- youden_cutoff(md_nuc[[DELTA_COL]][is_epi], as.character(md_nuc[[CT_CLEAN_COL]])[is_epi])
cut_youden <- best$thr[1]
message("[CUTOFF] NUC Youden-optimal threshold: ", signif(cut_youden, 4))
write.csv(data.frame(method = "Youden_opt", cutoff = cut_youden, tpr = best$tpr, fpr = best$fpr),
          file.path(OUTDIR, "cutoffs_NUC_Youden.csv"), row.names = FALSE)

d <- seu_nuc@meta.data[[DELTA_COL]]
seu_nuc$Stage_EL <- factor(ifelse(!is.finite(d), NA, ifelse(d < cut_youden, "Early", "Late")),
                           levels = c("Early", "Late"))
seu_nuc$celltype_TissueDiseaseState <- clean_tds(seu_nuc@meta.data[[TDS_COL]])

# ---- a-c: UMAPs ----
emb <- Embeddings(seu_nuc, reduction = NUC_UMAP)
dummy <- data.frame(x = max(emb[, 1], na.rm = TRUE) + 1, y = max(emb[, 2], na.rm = TRUE) + 1,
                    celltype_TissueDiseaseState = factor("Malignant_N", levels = tds_levels))

p_a <- FeaturePlot(seu_nuc, reduction = NUC_UMAP, features = "Pseudotime", pt.size = 0.1) +
  scale_color_viridis_c(option = "plasma", direction = 1) +
  ggtitle("FeaturePlot: Pseudotime (NUC)")
p_b <- DimPlot(seu_nuc, reduction = NUC_UMAP, group.by = TDS_COL, pt.size = 0.1) +
  geom_point(data = dummy, aes(x = x, y = y, color = celltype_TissueDiseaseState),
             inherit.aes = FALSE, size = 2, alpha = 1) +
  scale_color_manual(values = celltype_TDS.colors2, breaks = tds_levels,
                     limits = tds_levels, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_void() + ggtitle("UMAP: celltype (NUC)")
p_c <- DimPlot(seu_nuc, reduction = NUC_UMAP, group.by = "Stage_EL", pt.size = 0.1) +
  scale_color_manual(values = stage_cols_NUC, drop = FALSE) +
  theme_void() + ggtitle("UMAP: Stage_EL (NUC)")

pdf_safe(file.path(OUTDIR, "SupFig7a-c_NUC_UMAP.pdf"),
         { print(p_a); print(p_b); print(p_c) }, width = 5.1, height = 4.0)

# ---- d: Delta density by TDS (AT2/Malignant cells), Youden cutoff ----
df_nuc <- seu_nuc@meta.data[is_epi, , drop = FALSE] %>%
  transmute(delta = .data[[DELTA_COL]], tds = celltype_TissueDiseaseState, Stage_EL = Stage_EL)

p_d <- plot_delta_tds(df_nuc, cut_youden, paste0("Delta distribution: ", DELTA_COL))
ggsave(file.path(OUTDIR, "SupFig7d_NUC_delta_byTDS.pdf"), p_d,
       width = 8.3, height = 2.5, units = "in", device = cairo_pdf)

# ---- e: Early/Late composition by TDS ----
p_e <- plot_comp_tds(df_nuc %>% filter(!is.na(tds), !is.na(Stage_EL)), "Stage_EL", stage_cols_NUC,
                     "Early/Late composition by celltype_TissueDiseaseState",
                     paste0("Cutoff (Youden) = ", signif(cut_youden, 4), " | Delta = ", DELTA_COL))
pdf_safe(file.path(OUTDIR, "SupFig7e_NUC_EarlyLate_byTDS.pdf"), print(p_e), width = 4.6, height = 2.8)

write_labels(seu_nuc, "Stage_EL", file.path(OUTDIR, "early_late_labels_NUC.csv"))
saveRDS(seu_nuc, file.path(OUTDIR, "seu_NUC_withStage.rds"))

# =========================================================
# 5) scRNA-seq (CELL): NUC modules -> scores, Delta, quantile categories
# =========================================================
sc <- score_zmean(seu_cell, module_genes, sd_mode = "floor")
for (cc in score_cols) seu_cell@meta.data[[cc]] <- sc[[cc]]
seu_cell@meta.data[[DELTA_COL]] <- seu_cell$Score_M4 - seu_cell$Score_M1

d <- seu_cell@meta.data[[DELTA_COL]]
is_epi_c  <- as.character(seu_cell@meta.data[[CT_CLEAN_COL]]) %in% c("AT2", "Malignant")
delta_epi <- d[is_epi_c & is.finite(d)]
q_lo <- unname(quantile(delta_epi, probs = Q_EARLY))
q_hi <- unname(quantile(delta_epi, probs = Q_LATE))
message("[CUTOFF] CELL quantiles: 20%=", signif(q_lo, 4), " | 80%=", signif(q_hi, 4))
write.csv(data.frame(Q_EARLY = Q_EARLY, Q_LATE = Q_LATE, delta_q_lo = q_lo, delta_q_hi = q_hi,
                     n_epi = sum(is_epi_c), n_all = ncol(seu_cell)),
          file.path(OUTDIR, "cutoffs_CELL_quantile.csv"), row.names = FALSE)

st <- rep("mid_like", length(d))
st[is.finite(d) & d <= q_lo] <- "early_like"
st[is.finite(d) & d >= q_hi] <- "late_like"
seu_cell$Stage_like <- factor(st, levels = names(stage_cols_CELL))
seu_cell$celltype_TissueDiseaseState <- clean_tds(seu_cell@meta.data[[TDS_COL]])

# ---- f: sample-wise AT2 / malignant counts and proportions ----
df_f <- seu_cell@meta.data[is_epi_c, , drop = FALSE] %>%
  transmute(Sample = as.character(.data[[SAMPLE_COL]]),
            ct = factor(.data[[CT_CLEAN_COL]], levels = c("AT2", "Malignant"))) %>%
  count(Sample, ct, name = "n") %>%
  group_by(Sample) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()

p_f1 <- ggplot(df_f, aes(x = Sample, y = n, fill = ct)) +
  geom_col(color = "black", width = 0.75) +
  scale_fill_manual(values = epi.cell.colors2, drop = FALSE) +
  theme_plain(12) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "AT2 / malignant cells per sample (CELL)", x = NULL, y = "Cells", fill = NULL)
p_f2 <- p_f1 + aes(y = pct) + labs(y = "% within sample")
pdf_safe(file.path(OUTDIR, "SupFig7f_CELL_sample_AT2_Malignant.pdf"),
         { print(p_f1); print(p_f2) }, width = 6.0, height = 3.0)

# ---- g-h: UMAPs ----
p_g <- DimPlot(seu_cell, reduction = CELL_UMAP, group.by = TDS_COL, pt.size = 0.2) +
  scale_color_manual(values = celltype_TDS.colors2, drop = FALSE) +
  theme_void() + ggtitle("UMAP: celltype (CELL)")
p_h <- DimPlot(seu_cell, reduction = CELL_UMAP, group.by = "Stage_like", pt.size = 0.2) +
  scale_color_manual(values = stage_cols_CELL, drop = FALSE) +
  theme_void() + ggtitle("UMAP: Stage_like (CELL)")
pdf_safe(file.path(OUTDIR, "SupFig7g-h_CELL_UMAP.pdf"),
         { print(p_g); print(p_h) }, width = 5.3, height = 4.0)

# ---- i: Delta density by TDS ----
df_cell <- seu_cell@meta.data %>%
  transmute(delta = as.numeric(.data[[DELTA_COL]]), tds = celltype_TissueDiseaseState,
            Stage_like = Stage_like) %>%
  filter(is.finite(delta))

p_i <- plot_delta_tds(df_cell, CELL_VLINE, "Delta distribution: delta")
ggsave(file.path(OUTDIR, "SupFig7i_CELL_delta_byTDS.pdf"), p_i,
       width = 8.3, height = 2.5, units = "in", device = cairo_pdf)

# ---- j: Early/Mid/Late composition by TDS ----
p_j <- plot_comp_tds(df_cell %>% filter(!is.na(tds), !is.na(Stage_like)), "Stage_like", stage_cols_CELL,
                     "Early/Mid/Late composition by celltype_TissueDiseaseState",
                     paste0("Delta = ", DELTA_COL, "  |  Stage_like quantile cutoffs"))
pdf_safe(file.path(OUTDIR, "SupFig7j_CELL_Stage_byTDS.pdf"), print(p_j), width = 4.6, height = 2.8)

write_labels(seu_cell, "Stage_like", file.path(OUTDIR, "early_late_labels_CELL.csv"))
saveRDS(seu_cell, file.path(OUTDIR, "seu_CELL_withStage.rds"))

message("=== [DONE] Supplementary Fig. 7 ===")
