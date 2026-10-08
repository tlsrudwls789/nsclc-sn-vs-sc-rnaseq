#!/usr/bin/env Rscript

############################################################
# SupFig3_public_dataset.R
# ----------------------------------------------------------
# Supplementary Fig. 3: comprehensive analysis of the public
# NSCLC dataset (Renaut et al., "UL"; snRNA-seq + scRNA-seq).
#
#   a : UMAP, seven major cell types
#   b : UMAP, platform x tissue state (NUC_N, NUC_T, CELL_N, CELL_T)
#   c : UMAP, CopyKAT CNV status (aneuploid / diploid / not.defined)
#   d : per-sample major-lineage frequency across platform x tissue
#       groups (pairwise Wilcoxon, BH-FDR within each cell type)
#   e : paired sn vs sc sample-level proportions, Normal / Tumor
#       (paired t-test on NUC - CELL, Cohen's d, BH-FDR within tissue)
#   f : UMAP, annotated subtypes
#
# Inputs (via --key=value):
#   --public_rds   : integrated Seurat object (UMAP + CopyKAT calls)
#                    (default: data/Public_UL_final_wMeta_wCNV.CopyKAT.rds)
#   --celltype_col : major cell type column      (default: celltype2)
#   --subtype_col  : subtype column              (default: celltype)
#   --platform_col : platform column, NUC / CELL (default: platform)
#   --tissue_col   : tissue state column, N / T  (default: TissueDiseaseState)
#   --sample_col   : sample column               (default: orig.ident)
#   --copykat_col  : CopyKAT call column         (default: CopyKAT)
#   --outdir       : output directory            (default: SupFig3_out)
#
# Outputs:
#   - SupFig3a_UMAP_celltype.pdf
#   - SupFig3b_UMAP_platform_tissue.pdf
#   - SupFig3c_UMAP_CopyKAT.pdf
#   - SupFig3d_frequency_boxplot.pdf  (+ SupFig3d_Wilcoxon_BH.csv)
#   - SupFig3e_paired_proportions.pdf (+ SupFig3e_paired_ttest_BH.csv)
#   - SupFig3f_UMAP_subtypes.pdf
############################################################

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(stringr); library(ggplot2); library(ggrepel); library(RColorBrewer)
})
options(stringsAsFactors = FALSE)

# -----------------------------
# 0) Args
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NA_character_) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", key, "="), "", hit[1])
}

PUBLIC_RDS   <- get_arg("--public_rds", file.path("data", "Public_UL_final_wMeta_wCNV.CopyKAT.rds"))
CELLTYPE_COL <- get_arg("--celltype_col", "celltype2")
SUBTYPE_COL  <- get_arg("--subtype_col",  "celltype")
PLATFORM_COL <- get_arg("--platform_col", "platform")
TISSUE_COL   <- get_arg("--tissue_col",   "TissueDiseaseState")
SAMPLE_COL   <- get_arg("--sample_col",   "orig.ident")
COPYKAT_COL  <- get_arg("--copykat_col",  "CopyKAT")
OUTDIR       <- get_arg("--outdir", "SupFig3_out")

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
stopifnot(file.exists(PUBLIC_RDS))
sr3 <- readRDS(PUBLIC_RDS)

req <- c(CELLTYPE_COL, SUBTYPE_COL, PLATFORM_COL, TISSUE_COL, SAMPLE_COL, COPYKAT_COL)
miss <- setdiff(req, colnames(sr3@meta.data))
if (length(miss) > 0) stop("Missing metadata columns: ", paste(miss, collapse = ", "))

# -----------------------------
# 1) Helpers (shared with the Fig. 1 script)
# -----------------------------
save_pdf <- function(p, filename, w, h, use_cairo = TRUE) {
  f <- file.path(OUTDIR, filename)
  if (use_cairo && capabilities("cairo")) {
    ggsave(filename = f, plot = p, width = w, height = h, units = "in", device = cairo_pdf)
  } else {
    ggsave(filename = f, plot = p, width = w, height = h, units = "in", device = "pdf")
  }
}

theme_umap <- function(base_size = 9) {
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.width = unit(0.22, "cm"),
    legend.key.size  = unit(0.2,  "cm"),
    legend.background = element_blank(),
    legend.title = element_blank(),
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

umap_by <- function(obj, group_col, cols) {
  DimPlot(obj, reduction = "umap", group.by = group_col,
          pt.size = 0.1, label = FALSE, raster = FALSE) +
    scale_color_manual(values = cols, drop = FALSE) +
    theme_umap(base_size = 9) +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
}

# -----------------------------
# 2) Meta / factors / colors
# -----------------------------
sr3$celltype2          <- sr3@meta.data[[CELLTYPE_COL]]
sr3$celltype2[sr3$celltype2 == "B"] <- "B/Plasma"   # the figure labels this lineage B/Plasma
sr3$subtype            <- sr3@meta.data[[SUBTYPE_COL]]
sr3$platform           <- sr3@meta.data[[PLATFORM_COL]]
sr3$TissueDiseaseState <- sr3@meta.data[[TISSUE_COL]]
sr3$sample             <- sr3@meta.data[[SAMPLE_COL]]
sr3$CopyKAT            <- sr3@meta.data[[COPYKAT_COL]]

sr3$TissueDiseaseState <- factor(sr3$TissueDiseaseState, levels = c("N", "T"))
sr3$platform           <- factor(sr3$platform, levels = c("NUC", "CELL"))

p.TissueDiseaseState.order  <- c("NUC_N", "NUC_T", "CELL_N", "CELL_T")
p.TissueDiseaseState.colors <- c(
  "NUC_N"  = "#000957", "NUC_T"  = "#481E14",
  "CELL_N" = "#B1F0F7", "CELL_T" = "#C7B7A3")
sr3$platform_TissueDiseaseState <- factor(
  paste(sr3$platform, sr3$TissueDiseaseState, sep = "_"),
  levels = p.TissueDiseaseState.order)

celltype2.order <- c("Epithelial", "Endothelial", "Fibroblast", "Myeloid", "B/Plasma", "T/NK", "Mast")
celltype2.colors <- c(
  "Epithelial"  = brewer.pal(11, "RdBu")[2],
  "Endothelial" = brewer.pal(9,  "PuOr")[7],
  "Fibroblast"  = brewer.pal(9,  "RdPu")[6],
  "Myeloid"     = brewer.pal(9,  "Oranges")[8],
  "B/Plasma"    = brewer.pal(9,  "PuBuGn")[5],
  "T/NK"        = brewer.pal(9,  "Greens")[6],
  "Mast"        = brewer.pal(9,  "BuPu")[9])
sr3$celltype2 <- factor(sr3$celltype2, levels = celltype2.order)
if (any(is.na(sr3$celltype2))) warning("Cells with a major cell type outside celltype2.order: ",
                                       sum(is.na(sr3$celltype2)))

copykat.colors <- c("aneuploid" = "red", "diploid" = "#1F4172", "not.defined" = "darkgrey")

# -----------------------------
# 3) Panels a-c: UMAPs
# -----------------------------
save_pdf(umap_by(sr3, "celltype2", celltype2.colors),
         "SupFig3a_UMAP_celltype.pdf", w = 5.3, h = 5.3)
save_pdf(umap_by(sr3, "platform_TissueDiseaseState", p.TissueDiseaseState.colors),
         "SupFig3b_UMAP_platform_tissue.pdf", w = 5.3, h = 5.3)
save_pdf(umap_by(sr3, "CopyKAT", copykat.colors),
         "SupFig3c_UMAP_CopyKAT.pdf", w = 4.3, h = 3.0)

# -----------------------------
# 4) Panel d: per-sample frequency boxplot (Wilcoxon, BH per cell type)
# -----------------------------
data_freq <- sr3@meta.data %>%
  select(sample, platform_TissueDiseaseState, celltype2) %>%
  filter(!is.na(sample), !is.na(platform_TissueDiseaseState), !is.na(celltype2)) %>%
  group_by(sample, platform_TissueDiseaseState, celltype2) %>%
  summarise(n_celltype = n(), .groups = "drop") %>%
  group_by(sample, platform_TissueDiseaseState) %>%
  mutate(n_total = sum(n_celltype), frequency = n_celltype / n_total) %>%
  ungroup() %>%
  filter(!is.na(frequency))
data_freq$platform_TissueDiseaseState <- factor(
  as.character(data_freq$platform_TissueDiseaseState), levels = p.TissueDiseaseState.order)

my_comparisons <- list(
  c("NUC_N", "NUC_T"), c("NUC_N", "CELL_N"),
  c("NUC_T", "CELL_T"), c("CELL_N", "CELL_T"))

stat_res <- map_dfr(unique(as.character(data_freq$celltype2)), function(ct) {
  df_ct <- data_freq %>% filter(celltype2 == ct)
  map_dfr(my_comparisons, function(comp) {
    df_pair <- df_ct %>% filter(platform_TissueDiseaseState %in% comp)
    grp <- droplevels(df_pair$platform_TissueDiseaseState)
    n1 <- sum(grp == comp[1]); n2 <- sum(grp == comp[2])
    if (n_distinct(grp) == 2 && n1 > 0 && n2 > 0) {
      pval <- wilcox.test(frequency ~ platform_TissueDiseaseState, data = df_pair, exact = FALSE)$p.value
      tibble(celltype2 = ct, group1 = comp[1], group2 = comp[2], n1 = n1, n2 = n2, p.value = pval)
    } else NULL
  })
}) %>%
  filter(!(is.na(p.value) | is.infinite(p.value))) %>%
  group_by(celltype2) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  mutate(p.label = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE ~ NA_character_))
write.csv(stat_res, file.path(OUTDIR, "SupFig3d_Wilcoxon_BH.csv"), row.names = FALSE)

y_info <- data_freq %>%
  group_by(celltype2) %>%
  summarise(y.max = max(frequency, na.rm = TRUE),
            y.min = min(frequency, na.rm = TRUE), .groups = "drop") %>%
  mutate(y.range = ifelse(y.max - y.min == 0, 0.05, y.max - y.min),
         base_y  = y.max + 0.08 * y.range,
         tick    = 0.03 * y.range,
         celltype2 = as.character(celltype2))

x_levels <- levels(data_freq$platform_TissueDiseaseState)
stat_draw <- stat_res %>%
  filter(!is.na(p.label)) %>%
  left_join(y_info, by = "celltype2") %>%
  group_by(celltype2) %>%
  arrange(group1, group2, .by_group = TRUE) %>%
  mutate(y1 = base_y + (row_number() - 1) * (0.18 * y.range),
         y0 = y1 - tick,
         x1 = match(group1, x_levels),
         x2 = match(group2, x_levels),
         xmid = (x1 + x2) / 2) %>%
  ungroup() %>%
  filter(!is.na(x1), !is.na(x2)) %>%
  mutate(celltype2 = factor(celltype2, levels = levels(data_freq$celltype2)))

p_supfig3d <- ggplot(data_freq, aes(x = platform_TissueDiseaseState, y = frequency)) +
  geom_boxplot(aes(fill = platform_TissueDiseaseState), outlier.shape = NA, alpha = 0.7, width = 0.7) +
  geom_jitter(aes(color = platform_TissueDiseaseState), width = 0.15, size = 1.0, alpha = 0.8) +
  scale_fill_manual(values = p.TissueDiseaseState.colors, drop = FALSE) +
  scale_color_manual(values = c(
    "NUC_N" = "#0C0950", "NUC_T" = "#664343", "CELL_N" = "#4DA1A9", "CELL_T" = "#914F1E"), drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = "Frequency (0\u20131)") +
  facet_wrap(~ celltype2, scales = "free_y", nrow = 1) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13, face = "bold"),
    axis.text.y  = element_text(size = 13, colour = "black"),
    axis.title.y = element_text(size = 13, colour = "black"),
    strip.text   = element_text(size = 13, face = "bold", colour = "black"),
    legend.position = "right")

if (nrow(stat_draw) > 0) {
  p_supfig3d <- p_supfig3d +
    geom_segment(data = stat_draw, aes(x = x1, xend = x2, y = y1, yend = y1), inherit.aes = FALSE, linewidth = 0.45) +
    geom_segment(data = stat_draw, aes(x = x1, xend = x1, y = y0, yend = y1), inherit.aes = FALSE, linewidth = 0.45) +
    geom_segment(data = stat_draw, aes(x = x2, xend = x2, y = y0, yend = y1), inherit.aes = FALSE, linewidth = 0.45) +
    geom_text(data = stat_draw, aes(x = xmid, y = y1 + tick * 0.25, label = p.label),
              inherit.aes = FALSE, size = 4.0, vjust = 0)
}
save_pdf(p_supfig3d, "SupFig3d_frequency_boxplot.pdf", w = 14.0, h = 2.95)

# -----------------------------
# 5) Panel e: paired sn vs sc proportions (paired t-test, BH per tissue)
# -----------------------------
# Sample number shared by matched NUC / CELL libraries (e.g. "p033t" -> "33")
get_sample_key <- function(x) {
  x2  <- gsub("(?i)([_-]?NUC|[_-]?CELL)", "", as.character(x), perl = TRUE)
  num <- str_extract(x2, "(?i)(?<=p)\\d+|\\d+")
  ifelse(is.na(num), NA_character_, as.character(as.integer(num)))
}

calc_prop <- function(obj, platform_label) {
  obj@meta.data %>%
    transmute(sample, celltype2, TissueDiseaseState) %>%
    group_by(sample, TissueDiseaseState, celltype2) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(sample, TissueDiseaseState) %>%
    mutate(proportion = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(sample_num  = get_sample_key(sample),
           platform    = platform_label,
           TissueGroup = ifelse(TissueDiseaseState == "N", "Normal", "Tumor"),
           sample_lab  = ifelse(TissueGroup == "Tumor", paste0("T", sample_num), paste0("N", sample_num)),
           pair_id     = paste0(sample_num, "||", TissueGroup))
}

nuc.obj  <- subset(sr3, platform == "NUC")
cell.obj <- subset(sr3, platform == "CELL")
pair_df <- bind_rows(calc_prop(nuc.obj, "NUC"), calc_prop(cell.obj, "CELL")) %>%
  filter(!is.na(sample_num))

ann <- pair_df %>%
  select(sample_num, TissueGroup, celltype2, platform, proportion) %>%
  pivot_wider(names_from = platform, values_from = proportion) %>%
  filter(!is.na(NUC) & !is.na(CELL)) %>%
  group_by(TissueGroup, celltype2) %>%
  summarise(
    n_pairs   = n(),
    mean_diff = mean(NUC - CELL, na.rm = TRUE),
    sd_diff   = sd(NUC - CELL, na.rm = TRUE),
    p_value   = {
      d <- (NUC - CELL); d <- d[!is.na(d)]
      if (length(d) < 2 || sd(d) == 0) NA_real_ else t.test(d, mu = 0)$p.value
    },
    .groups = "drop") %>%
  group_by(TissueGroup) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    cohens_d  = mean_diff / sd_diff,
    dir_label = case_when(is.na(mean_diff) ~ "",
                          mean_diff > 0 ~ "\u2191NUC",
                          mean_diff < 0 ~ "\u2193NUC",
                          TRUE ~ "0"),
    label = case_when(is.na(p_adj)  ~ "n.d.",
                      p_adj < 0.001 ~ "FDR<0.001",
                      TRUE          ~ sprintf("FDR=%.3f", p_adj)))
write.csv(ann, file.path(OUTDIR, "SupFig3e_paired_ttest_BH.csv"), row.names = FALSE)

pair_df$platform    <- factor(pair_df$platform, levels = c("NUC", "CELL"))
pair_df$TissueGroup <- factor(pair_df$TissueGroup, levels = c("Normal", "Tumor"))
pair_df$celltype2   <- factor(pair_df$celltype2)

ypos <- pair_df %>%
  group_by(TissueGroup, celltype2) %>%
  summarise(y_max = max(proportion, na.rm = TRUE), .groups = "drop") %>%
  mutate(y_pos = y_max + 0.08 * pmax(y_max, 1))
ann <- ann %>%
  mutate(TissueGroup = factor(TissueGroup, levels = levels(pair_df$TissueGroup)),
         celltype2   = factor(celltype2,   levels = levels(pair_df$celltype2))) %>%
  left_join(ypos, by = c("TissueGroup", "celltype2"))

p_supfig3e <- ggplot(pair_df, aes(x = platform, y = proportion, group = pair_id)) +
  geom_line(color = "gray50", alpha = 0.55, linewidth = 0.35) +
  geom_point(aes(color = platform), size = 2.4, alpha = 0.9) +
  geom_text_repel(aes(label = sample_lab), size = 2.5, max.overlaps = 70, force = 1.1,
                  min.segment.length = 0, segment.color = "gray70") +
  scale_color_manual(values = c("NUC" = "#1f77b4", "CELL" = "#ff7f0e")) +
  labs(x = NULL, y = "Proportion (%)") +
  facet_grid(TissueGroup ~ celltype2, scales = "fixed") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  theme_bw(base_size = 8.5) +
  theme(
    strip.text  = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(6, 16, 6, 6)) +
  coord_cartesian(clip = "off") +
  geom_label(data = ann, aes(x = "NUC", y = y_pos, label = label), inherit.aes = FALSE,
             hjust = 0, vjust = 0, size = 2.8, label.size = 0.25, fill = "white")
save_pdf(p_supfig3e, "SupFig3e_paired_proportions.pdf", w = 9.5, h = 3.5)

# -----------------------------
# 6) Panel f: UMAP of annotated subtypes
# -----------------------------
subtype.colors <- c(
  colorRampPalette(rev(brewer.pal(9, "Reds")[4:9]))(6),   # Epithelial
  colorRampPalette(brewer.pal(9, "BuPu")[4:8])(5),        # Endothelial
  colorRampPalette(brewer.pal(9, "PuRd")[6:8])(3),        # Fibroblast 1
  colorRampPalette(brewer.pal(9, "RdPu")[6:8])(3),        # Fibroblast 2
  colorRampPalette(brewer.pal(9, "Blues")[4:6])(3),       # B / Plasma
  colorRampPalette(brewer.pal(9, "Greens")[7:9])(3),      # T 1
  colorRampPalette(brewer.pal(11, "PRGn")[9:11])(3),      # T 2
  colorRampPalette(brewer.pal(9, "BuGn")[7])(1),          # NK
  colorRampPalette(brewer.pal(9, "Oranges")[4:9])(5),     # Myeloid
  colorRampPalette(brewer.pal(11, "BrBG")[1])(1),         # DC
  colorRampPalette(brewer.pal(8, "Purples")[8])(1),       # Mast
  colorRampPalette(brewer.pal(9, "Greys")[4:8])(1))       # Other

p_supfig3f <- DimPlot(sr3, group.by = "subtype", cols = subtype.colors, label = FALSE,
                      raster = FALSE, pt.size = 0.1, seed = 21, shuffle = TRUE) +
  theme_umap(base_size = 9) +
  guides(colour = guide_legend(override.aes = list(size = 2.5)))
save_pdf(p_supfig3f, "SupFig3f_UMAP_subtypes.pdf", w = 8.8, h = 4.2)

message("DONE: Supplementary Fig. 3 panels saved in ", normalizePath(OUTDIR))
