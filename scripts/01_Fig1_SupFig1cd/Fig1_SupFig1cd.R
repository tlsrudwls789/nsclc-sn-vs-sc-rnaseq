suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(purrr)
  library(RColorBrewer)})

## =============================================================================
## Figure 1b-d and Supplementary Figure 1c-d (SMC cohort, snRNA-seq + scRNA-seq)
##   Fig. 1b  UMAPs: cell type, platform x tissue state, CopyKAT CNV status
##   Fig. 1c  major cell type frequencies per platform-tissue group (Wilcoxon)
##   Fig. 1d  within-patient paired proportions, snRNA-seq vs scRNA-seq
##   Sup. Fig. 1c  FeaturePlots of lineage markers
##   Sup. Fig. 1d  sub-cell type UMAPs per platform
## Usage: DATA_DIR=<dir with the integrated Seurat RDS> OUT_DIR=<output dir> Rscript <this script>
##   required meta columns: orig.ident, platform, TissueDiseaseState, celltype,
##   celltype2, celltype3, CopyKAT
## =============================================================================

## =========================
## 0) Config
## =========================
pc.num  <- 21
res.num <- 0.6
label   <- "LUNG_Integration"

data.dir <- Sys.getenv("DATA_DIR", "data")
out.dir  <- Sys.getenv("OUT_DIR",  "Fig1_out")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
setwd(out.dir)

sr3.rds <- file.path(
  data.dir,
  sprintf("LUNG_Integration.qc500_pc%s_res%s_Seurat4.4.0_final_wMeta.rds", pc.num, res.num))
stopifnot(file.exists(sr3.rds))
sr3 <- readRDS(sr3.rds)

## =========================
## 1) Helpers
## =========================
save_pdf <- function(p, filename, w, h, use_cairo = TRUE) {
  # GitHub portability:
  # - use_cairo=TRUE: nicer font handling if cairo available
  # - fallback to default pdf if cairo not available
  if (use_cairo && capabilities("cairo")) {
    ggsave(filename = filename, plot = p, width = w, height = h, units = "in", device = cairo_pdf)
  } else {
    ggsave(filename = filename, plot = p, width = w, height = h, units = "in", device = "pdf")
  }
}

theme_umap <- function(base_size = 9) {
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.width = unit(0.22, "cm"),
    legend.key.size  = unit(0.2,  "cm"),
    legend.background = element_blank(),
    legend.title = element_text(size = base_size),
    legend.text  = element_text(size = base_size, colour = "black"),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.line  = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
}

## unified paired df builder (works for celltype or celltype4)
make_pairwise_df <- function(nuc_obj, cell_obj,
                             group_label,
                             celltype_col = "celltype",
                             tissue_col   = "TissueDiseaseState",
                             orig_col     = "orig.ident") {
  
  stopifnot(celltype_col %in% colnames(nuc_obj@meta.data))
  stopifnot(celltype_col %in% colnames(cell_obj@meta.data))
  stopifnot(tissue_col   %in% colnames(nuc_obj@meta.data))
  stopifnot(tissue_col   %in% colnames(cell_obj@meta.data))
  stopifnot(orig_col     %in% colnames(nuc_obj@meta.data))
  stopifnot(orig_col     %in% colnames(cell_obj@meta.data))
  
  calc_prop <- function(obj, modality_label) {
    md <- obj@meta.data
    md %>%
      count(.data[[orig_col]], .data[[celltype_col]], .data[[tissue_col]]) %>%
      group_by(.data[[orig_col]], .data[[tissue_col]]) %>%
      mutate(total = sum(n), proportion = n / total * 100) %>%
      ungroup() %>%
      transmute(
        orig.ident  = .data[[orig_col]],
        celltype    = .data[[celltype_col]],
        TissueDiseaseState = .data[[tissue_col]],
        modality    = modality_label,
        sample_id   = gsub("LUNG-|-NUC|-CELL", "", orig.ident),
        TissueGroup = ifelse(TissueDiseaseState == "N", "Normal", "Tumor"),
        proportion  = proportion
      )
  }
  
  bind_rows(
    calc_prop(nuc_obj,  "snRNA-seq"),
    calc_prop(cell_obj, "scRNA-seq")
  ) %>%
    mutate(Group = group_label)
}

plot_paired <- function(df, group_name,
                        celltype_order = NULL,
                        split_by_tissue = FALSE,
                        label_samples = FALSE,
                        scales = "fixed",
                        base_size = 9) {
  
  d <- df %>% filter(Group == group_name)
  d$modality <- factor(d$modality, levels = c("snRNA-seq", "scRNA-seq"))
  
  if (!is.null(celltype_order)) d$celltype <- factor(d$celltype, levels = celltype_order)
  
  p <- ggplot(d, aes(x = modality, y = proportion, group = sample_id)) +
    geom_point(aes(color = modality), size = 2.4, alpha = 0.9) +
    geom_line(color = "gray50", alpha = 0.5) +
    scale_color_manual(values = c("snRNA-seq" = "#1f77b4", "scRNA-seq" = "#ff7f0e")) +
    labs(x = NULL, y = "Proportion (%)") +
    theme_bw(base_size = base_size) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid.major.x = element_blank()
    )
  
  if (label_samples) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = sample_id),
      size = 2.6, max.overlaps = 70, force = 1
    )
  }
  
  if (split_by_tissue) {
    d$TissueGroup <- factor(d$TissueGroup, levels = c("Normal", "Tumor"))
    p <- p + facet_grid(TissueGroup ~ celltype, scales = scales)
  } else {
    p <- p + facet_wrap(~ celltype, scales = "free_y", nrow = 1)
  }
  p
}

## =========================
## 2) Meta / factors / colors
## =========================
sr3$TissueDiseaseState <- factor(sr3$TissueDiseaseState, levels = c("N","T"))
sr3$platform <- factor(sr3$platform, levels = c("NUC","CELL"))

sr3$platform_TissueDiseaseState <- paste(sr3$platform, sr3$TissueDiseaseState, sep = "_")
p.TissueDiseaseState.order  <- c("NUC_N","NUC_T","CELL_N","CELL_T")
p.TissueDiseaseState.colors <- c(
  "NUC_N"  = "#000957", "NUC_T"  = "#481E14",
  "CELL_N" = "#B1F0F7", "CELL_T" = "#C7B7A3")

celltype.order <- c("Epithelial", "Endothelial", "Fibroblast", "Myeloid", "Cycling", "B/Plasma", "T/NK", "Mast")
celltype.colors <- c(
  "Epithelial"  = brewer.pal(11, "RdBu")[2],
  "Endothelial" = brewer.pal(9,  "PuOr")[7],
  "Fibroblast"  = brewer.pal(9,  "RdPu")[6],
  "Myeloid"     = brewer.pal(9,  "Oranges")[8],
  "Cycling"     = brewer.pal(9,  "Greens")[9],
  "B/Plasma"    = brewer.pal(9,  "PuBuGn")[5],
  "T/NK"        = brewer.pal(9,  "Greens")[6],
  "Mast"        = brewer.pal(9,  "BuPu")[9])
sr3$celltype <- factor(sr3$celltype, levels = celltype.order)

nuc.obj  <- subset(sr3, platform == "NUC")
cell.obj <- subset(sr3, platform == "CELL")

## =========================
## 3) Figure 1b-d
## =========================

## Fig 1b (middle): UMAP by platform x tissue state
p_fig1b <- DimPlot(
  sr3, reduction = "umap",
  group.by = "platform_TissueDiseaseState",
  pt.size = 0.1, label = FALSE) +
  scale_color_manual(values = p.TissueDiseaseState.colors, breaks = p.TissueDiseaseState.order) +
  theme_umap(base_size = 8) +
  guides(colour = guide_legend(override.aes = list(size = 2.5)))

save_pdf(p_fig1b, sprintf("Fig1b_UMAP_platformTissue_%s_pc%s.pdf", label, pc.num), w = 6, h = 6)

## Fig 1b (left): UMAP by cell type
p_fig1_celltype <- DimPlot(sr3, reduction = "umap", group.by = "celltype", pt.size = 0.1) +
  scale_color_manual(values = celltype.colors) +
  theme_umap(base_size = 9) +
  theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 2.5)))

save_pdf(p_fig1_celltype, sprintf("Fig1b_UMAP_celltype_%s_pc%s.pdf", label, pc.num), w = 6.5, h = 6)

## Fig 1b (right): UMAP by CopyKAT CNV status (cells without a call are shown as not.defined)
copykat.colors <- c("aneuploid" = "red", "diploid" = "#1F4172", "not.defined" = "darkgrey")
sr3$CopyKAT_plot <- ifelse(sr3$CopyKAT %in% c("aneuploid", "diploid"), as.character(sr3$CopyKAT), "not.defined")
p_fig1b_cnv <- DimPlot(sr3, reduction = "umap", group.by = "CopyKAT_plot", pt.size = 0.1) +
  scale_color_manual(values = copykat.colors) +
  theme_umap(base_size = 9) +
  theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 2.5)))
save_pdf(p_fig1b_cnv, sprintf("Fig1b_UMAP_CopyKAT_%s_pc%s.pdf", label, pc.num), w = 6.5, h = 6)

## Fig 1d: paired proportions (celltype), split Normal/Tumor
pair_df1d <- make_pairwise_df(nuc.obj, cell.obj, group_label = "Total", celltype_col = "celltype")
p_fig1d   <- plot_paired(pair_df1d, "Total", celltype_order = celltype.order,
                         split_by_tissue = TRUE, label_samples = TRUE, scales = "fixed", base_size = 8.5)

save_pdf(p_fig1d, sprintf("Fig1d_PairedProp_NormalTumor_%s_pc%s.pdf", label, pc.num), w = 9.5, h = 3.5)

## Fig 1c: frequency boxplot per cell type with Wilcoxon labels (raw P)
data_freq <- sr3@meta.data %>%
  group_by(orig.ident, platform_TissueDiseaseState, celltype) %>%
  summarise(n_celltype = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(n_total = sum(n_celltype),
         frequency = n_celltype / n_total) %>%
  ungroup() %>%
  filter(!is.na(frequency))

data_freq$platform_TissueDiseaseState <- factor(data_freq$platform_TissueDiseaseState, levels = p.TissueDiseaseState.order)

my_comparisons <- list(
  c("NUC_N", "NUC_T"), c("NUC_N", "CELL_N"),
  c("NUC_T", "CELL_T"), c("CELL_N","CELL_T"))

stat_res <- map_dfr(unique(data_freq$celltype), function(ct) {
  df_ct <- data_freq %>% filter(celltype == ct)
  map_dfr(my_comparisons, function(comp) {
    df_pair <- df_ct %>% filter(platform_TissueDiseaseState %in% comp)
    if (n_distinct(df_pair$platform_TissueDiseaseState) == 2) {
      pval <- wilcox.test(frequency ~ platform_TissueDiseaseState, data = df_pair)$p.value
      tibble(celltype = ct, group1 = comp[1], group2 = comp[2], p.value = pval)
    }
  })
}) %>%
  mutate(
    p.value.label = case_when(
      p.value < 0.05  ~ "***",
      p.value < 0.115 ~ "**",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(p.value.label))

y_max <- data_freq %>%
  group_by(celltype) %>%
  summarise(y.max = max(frequency, na.rm = TRUE) + 0.05, .groups = "drop")

stat_res <- stat_res %>%
  left_join(y_max, by = "celltype") %>%
  group_by(celltype) %>%
  mutate(y.position = y.max + 0.05 * (row_number() - 1)) %>%
  ungroup()

p_fig1c <- ggplot(data_freq, aes(x = platform_TissueDiseaseState, y = frequency)) +
  geom_boxplot(aes(fill = platform_TissueDiseaseState), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = platform_TissueDiseaseState), width = 0.15, size = 1, alpha = 0.8) +
  scale_fill_manual(values = p.TissueDiseaseState.colors) +
  scale_color_manual(values = c(
    "NUC_N"="#0C0950", "NUC_T"="#664343", "CELL_N"="#4DA1A9", "CELL_T"="#914F1E")) +
  labs(x = NULL, y = "Frequency") +
  facet_wrap(~ celltype, scales = "free_y", nrow = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13, face = "bold"),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "none") +
  stat_pvalue_manual(
    stat_res,
    label = "p.value.label",
    tip.length = 0.03,
    bracket.size = 0.5,
    size = 4.0,
    xmin = "group1", xmax = "group2",
    y.position = "y.position")

save_pdf(p_fig1c, sprintf("Fig1c_FrequencyBoxplot_Wilcox_%s_pc%s.pdf", label, pc.num), w = 16.2, h = 2.95)

## =========================
## 4) Supplementary Figure 1c-d
## =========================
target.list <- list(
  Epithelial   = c("EPCAM","KRT18"),
  Endothelial  = c("PECAM1","CLDN5"),
  Fibroblast   = c("DCN","THY1","COL1A1","COL1A2"),
  Tcell        = c("CD3D","CD3E","CD3G","TRAC"),
  Cycling      = c("STMN1"),
  NK           = c("NCR1","NKG7"),
  Bcell        = c("MS4A1","CD79A","IGHM"),
  Plasma       = c("CREB3L2"),
  Myeloid      = c("LYZ","MARCO","CD68","FCGR3A"),
  CD4T         = c("CD4","CD40LG","FOXP3"),
  CD8T         = c("CD8A","CD8B"),
  Mast         = c("KIT"))

plot_feature_and_save <- function(obj, gene, prefix, pt.size = 0.1) {
  p <- FeaturePlot(obj, features = gene, reduction = "umap", pt.size = pt.size) +
    theme_umap(base_size = 8) +
    theme(legend.position = "none")
  save_pdf(p, sprintf("%s_%s_%s_pc%s.pdf", prefix, gene, label, pc.num), w = 3, h = 3.3)
}

for (ct in names(target.list)) {
  for (g in target.list[[ct]]) {
    tryCatch(
      plot_feature_and_save(sr3, g, prefix = paste0("SupFig1c_FeaturePlot_", ct)),
      error = function(e) message("ERROR: ", ct, " / ", g, " : ", conditionMessage(e))
    )
  }
}

p_sup1d_nuc_ct2  <- DimPlot(nuc.obj,  reduction = "umap", group.by = "celltype2", pt.size = 0.1) + theme_umap(8)
p_sup1d_nuc_ct3  <- DimPlot(nuc.obj,  reduction = "umap", group.by = "celltype3", pt.size = 0.1) + theme_umap(8)
p_sup1d_cell_ct2 <- DimPlot(cell.obj, reduction = "umap", group.by = "celltype2", pt.size = 0.1) + theme_umap(8)
p_sup1d_cell_ct3 <- DimPlot(cell.obj, reduction = "umap", group.by = "celltype3", pt.size = 0.1) + theme_umap(8)

save_pdf(p_sup1d_nuc_ct2,  "SupFig1d_UMAP_NUC_celltype2.pdf",  w = 5.3, h = 4.8)
save_pdf(p_sup1d_nuc_ct3,  "SupFig1d_UMAP_NUC_celltype3.pdf",  w = 5.3, h = 4.8)
save_pdf(p_sup1d_cell_ct2, "SupFig1d_UMAP_CELL_celltype2.pdf", w = 5.3, h = 4.8)
save_pdf(p_sup1d_cell_ct3, "SupFig1d_UMAP_CELL_celltype3.pdf", w = 5.3, h = 4.8)

message("Done. Outputs saved to: ", out.dir)
