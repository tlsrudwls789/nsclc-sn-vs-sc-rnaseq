# =============================================================================
# 02_Analysis_Result2 (Fig2, Sup Fig 3) — GitHub-ready SINGLE script (functions-only)
# -----------------------------------------------------------------------------
# Policy:
# - NO auto-detection of celltype columns (you must pass which one to use)
# - NO setwd(), NO file I/O, NO directory dependence
# - Pure functions: (Seurat objects + explicit column names) -> (dataframes / ggplots)
#
# Core outputs supported:
# 1) Fig2: Dot heatmap (dot size = proportion) comparing NUC vs CELL by sample
# 2) Sup Fig 3f-like: frequency boxplot (per-sample) + Wilcoxon p-value labels
#    - p-value threshold updated to p < 0.005
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

# =============================================================================
# 0) Small utilities
# =============================================================================

.stopif_cols <- function(md, cols) {
  miss <- setdiff(cols, colnames(md))
  if (length(miss) > 0) stop("Missing meta.data columns: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

.make_sample_short <- function(x) {
  x <- as.character(x)
  x <- gsub("LUNG-", "", x)
  x <- gsub("-NUC|-CELL", "", x)
  x <- gsub("_Nuc|_Cell|_Normal|_Tumor", "", x)
  x
}

# Standardize tissue to "Normal"/"Tumor"
# - If input already is Normal/Tumor (or variants), keep
# - If input is N/T, convert
.standardize_tissue <- function(x) {
  x <- as.character(x)
  # N/T style
  if (all(na.omit(unique(x)) %in% c("N", "T"))) {
    return(ifelse(x == "N", "Normal", ifelse(x == "T", "Tumor", NA)))
  }
  # Normal/Tumor style (case-insensitive)
  xl <- tolower(x)
  out <- ifelse(xl %in% c("normal", "n"), "Normal",
                ifelse(xl %in% c("tumor", "tumour", "t"), "Tumor", x))
  out
}

.scale_size_prop <- function(min_pt = 2, max_pt = 8) {
  scale_radius(
    name   = "Proportion (%)",
    limits = c(0, 100),
    breaks = c(20, 40, 60, 80, 100),
    range  = c(min_pt, max_pt)
  )
}

# =============================================================================
# 1) Core: proportion table (explicit column mapping)
# =============================================================================
# Required explicit columns in seu@meta.data:
# - sample_col: sample/orig ident
# - tissue_col: tissue state (N/T or Normal/Tumor)
# - platform_col: platform label (NUC/CELL)  [if you filter by platform]
# - celltype_col: chosen celltype level to use (e.g., celltype4 for SMC, celltype3 for UL)
#
# Returns: complete grid of (sample_id × celltype) with n and proportion (%)
get_celltype_proportion <- function(seu,
                                    celltype_col,
                                    sample_col,
                                    tissue_col,
                                    platform_col = NULL,
                                    tissue_group = c("Normal", "Tumor"),
                                    platform_value = NULL,
                                    celltype_levels = NULL) {
  stopifnot(inherits(seu, "Seurat"))
  tissue_group <- match.arg(tissue_group)
  
  md <- seu@meta.data
  need <- c(celltype_col, sample_col, tissue_col)
  if (!is.null(platform_col)) need <- c(need, platform_col)
  .stopif_cols(md, need)
  
  df <- md %>%
    transmute(
      sample_id = as.character(.data[[sample_col]]),
      tissue_std = .standardize_tissue(.data[[tissue_col]]),
      platform   = if (!is.null(platform_col)) as.character(.data[[platform_col]]) else NA_character_,
      celltype   = as.character(.data[[celltype_col]])
    ) %>%
    filter(tissue_std == tissue_group)
  
  if (!is.null(platform_col) && !is.null(platform_value)) {
    df <- df %>% filter(platform == platform_value)
  }
  
  df_raw <- df %>%
    count(sample_id, celltype, name = "n") %>%
    group_by(sample_id) %>%
    mutate(proportion = ifelse(sum(n) == 0, 0, n / sum(n) * 100)) %>%
    ungroup()
  
  if (is.null(celltype_levels)) celltype_levels <- sort(unique(df$celltype))
  
  all_comb <- expand.grid(
    sample_id = unique(df_raw$sample_id),
    celltype  = celltype_levels,
    stringsAsFactors = FALSE
  )
  
  left_join(all_comb, df_raw, by = c("sample_id", "celltype")) %>%
    mutate(n = ifelse(is.na(n), 0, n),
           proportion = ifelse(is.na(proportion), 0, proportion))
}

# =============================================================================
# 2) Fig2: Dot heatmap comparing NUC vs CELL (explicit mapping)
# =============================================================================
# You must provide:
# - epi/endo/fibro objects for NUC and CELL
# - meta columns explicitly for each object (they can be different names across objects)
#
# alveolar_types / bronchial_types are explicit character vectors in the chosen celltype space.
plot_dot_heatmap_two_platforms <- function(seu_nuc, seu_cell,
                                           # explicit mappings for NUC object
                                           nuc_celltype_col, nuc_sample_col, nuc_tissue_col,
                                           nuc_platform_col = NULL, nuc_platform_value = NULL,
                                           # explicit mappings for CELL object
                                           cell_celltype_col, cell_sample_col, cell_tissue_col,
                                           cell_platform_col = NULL, cell_platform_value = NULL,
                                           # plot spec
                                           tissue_group = c("Normal", "Tumor"),
                                           alveolar_types,
                                           bronchial_types,
                                           min_pt = 2, max_pt = 8,
                                           base_size = 11) {
  tissue_group <- match.arg(tissue_group)
  celltype_levels_plot <- c(alveolar_types, bronchial_types)
  
  df_n <- get_celltype_proportion(
    seu_nuc,
    celltype_col = nuc_celltype_col,
    sample_col   = nuc_sample_col,
    tissue_col   = nuc_tissue_col,
    platform_col = nuc_platform_col,
    tissue_group = tissue_group,
    platform_value = nuc_platform_value,
    celltype_levels = celltype_levels_plot
  ) %>%
    transmute(
      sample_short = .make_sample_short(sample_id),
      ctype = factor(celltype, levels = celltype_levels_plot),
      anatomical_group = ifelse(celltype %in% alveolar_types, "Alveolar", "Bronchial"),
      proportion = pmin(pmax(proportion, 0), 100),
      platform2 = "NUC"
    )
  
  df_c <- get_celltype_proportion(
    seu_cell,
    celltype_col = cell_celltype_col,
    sample_col   = cell_sample_col,
    tissue_col   = cell_tissue_col,
    platform_col = cell_platform_col,
    tissue_group = tissue_group,
    platform_value = cell_platform_value,
    celltype_levels = celltype_levels_plot
  ) %>%
    transmute(
      sample_short = .make_sample_short(sample_id),
      ctype = factor(celltype, levels = celltype_levels_plot),
      anatomical_group = ifelse(celltype %in% alveolar_types, "Alveolar", "Bronchial"),
      proportion = pmin(pmax(proportion, 0), 100),
      platform2 = "CELL"
    )
  
  df_plot <- bind_rows(df_n, df_c) %>%
    mutate(platform2 = factor(platform2, levels = c("NUC", "CELL")))
  
  ggplot(df_plot, aes(x = sample_short, y = ctype)) +
    geom_point(aes(size = proportion, fill = anatomical_group),
               shape = 21, color = "black") +
    .scale_size_prop(min_pt, max_pt) +
    scale_fill_manual(
      values = c(Alveolar = "#CADEFC", Bronchial = "#556FB5"),
      guide = guide_legend(title = "Group", override.aes = list(size = 5))
    ) +
    labs(x = "Sample", y = "Cell Type") +
    facet_grid(cols = vars(platform2), drop = FALSE) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, colour = "black"),
      axis.text.y = element_text(size = 9.5, colour = "black"),
      legend.position = "right"
    )
}

# =============================================================================
# 3) Sup Fig 3f: per-sample frequency boxplot + Wilcoxon labels (explicit mapping)
# =============================================================================
# - No auto-detect: you pass sample_col, group_col, platform_col, tissue_col
# - Creates platform_TissueDiseaseState internally: e.g., NUC_N, CELL_T
# - p-value threshold updated: p < 0.005
freq_boxplot_stats_explicit <- function(seu,
                                        group_col,         # e.g., "celltype3"
                                        sample_col,        # e.g., "orig.ident"
                                        platform_col,      # e.g., "platform" (NUC/CELL)
                                        tissue_col,        # e.g., "TissueDiseaseState" or "tissue"
                                        comparisons = list(
                                          c("NUC_N", "NUC_T"),
                                          c("NUC_N", "CELL_N"),
                                          c("NUC_T", "CELL_T"),
                                          c("CELL_N", "CELL_T")
                                        ),
                                        p_cut = 0.005,
                                        p_cut2 = 0.05,
                                        p_adjust = NULL) {
  stopifnot(inherits(seu, "Seurat"))
  md <- seu@meta.data
  .stopif_cols(md, c(group_col, sample_col, platform_col, tissue_col))
  
  df <- md %>%
    transmute(
      sample_id = as.character(.data[[sample_col]]),
      group     = as.character(.data[[group_col]]),
      platform  = as.character(.data[[platform_col]]),
      tissue_raw = as.character(.data[[tissue_col]])
    ) %>%
    mutate(
      TissueDiseaseState = ifelse(tolower(tissue_raw) %in% c("normal", "n"), "N",
                                  ifelse(tolower(tissue_raw) %in% c("tumor", "tumour", "t"), "T", tissue_raw)),
      platform_TissueDiseaseState = paste0(platform, "_", TissueDiseaseState)
    ) %>%
    filter(!is.na(group), !is.na(sample_id), !is.na(platform_TissueDiseaseState))
  
  # per-sample frequency of "group"
  data_freq <- df %>%
    group_by(sample_id, platform_TissueDiseaseState, group) %>%
    summarise(n_group = n(), .groups = "drop") %>%
    group_by(sample_id) %>%
    mutate(n_total = sum(n_group),
           frequency = n_group / n_total) %>%
    ungroup() %>%
    filter(!is.na(frequency))
  
  # order if canonical levels exist
  canonical_levels <- c("NUC_N","NUC_T","CELL_N","CELL_T")
  if (all(canonical_levels %in% unique(as.character(data_freq$platform_TissueDiseaseState)))) {
    data_freq$platform_TissueDiseaseState <- factor(data_freq$platform_TissueDiseaseState, levels = canonical_levels)
  } else {
    data_freq$platform_TissueDiseaseState <- factor(data_freq$platform_TissueDiseaseState)
  }
  
  # wilcox per facet + comp
  stat_res <- purrr::map_dfr(unique(data_freq$group), function(g) {
    df_g <- data_freq %>% filter(group == g)
    purrr::map_dfr(comparisons, function(comp) {
      df_pair <- df_g %>% filter(platform_TissueDiseaseState %in% comp)
      if (n_distinct(df_pair$platform_TissueDiseaseState) == 2) {
        pval <- wilcox.test(frequency ~ platform_TissueDiseaseState, data = df_pair)$p.value
        tibble(group = g, group1 = comp[1], group2 = comp[2], p.value = pval)
      } else NULL
    })
  })
  
  if (!is.null(p_adjust)) {
    stat_res <- stat_res %>% mutate(p.value.use = p.adjust(p.value, method = p_adjust))
  } else {
    stat_res <- stat_res %>% mutate(p.value.use = p.value)
  }
  
  # *** Threshold: p < 0.005 ***
  stat_res <- stat_res %>%
    mutate(
      p.value.label = case_when(
        p.value.use < p_cut  ~ "**",
        p.value.use < p_cut2 ~ "***",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(p.value.label))
  
  # y positions per facet
  y_max <- data_freq %>%
    group_by(group) %>%
    summarise(y.max = max(frequency, na.rm = TRUE) + 0.05, .groups = "drop")
  
  stat_res <- stat_res %>%
    left_join(y_max, by = "group") %>%
    group_by(group) %>%
    mutate(y.position = y.max + 0.05 * (row_number() - 1)) %>%
    ungroup()
  
  list(data_freq = data_freq, stat_res = stat_res)
}

plot_freq_boxplot_explicit <- function(freq_out,
                                       colors_fill = c("NUC_N"="#000957","NUC_T"="#481E14","CELL_N"="#B1F0F7","CELL_T"="#C7B7A3"),
                                       colors_point= c("NUC_N"="#0C0950","NUC_T"="#664343","CELL_N"="#4DA1A9","CELL_T"="#914F1E"),
                                       facet_nrow = 1) {
  data_freq <- freq_out$data_freq
  stat_res  <- freq_out$stat_res
  
  p <- ggplot(data_freq, aes(x = platform_TissueDiseaseState, y = frequency)) +
    geom_boxplot(aes(fill = platform_TissueDiseaseState), outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = platform_TissueDiseaseState), width = 0.15, size = 1, alpha = 0.8) +
    scale_fill_manual(values = colors_fill) +
    scale_color_manual(values = colors_point) +
    labs(x = NULL, y = "Frequency") +
    facet_wrap(~ group, scales = "free_y", nrow = facet_nrow) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13, face = "bold"),
      axis.text.y = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      strip.text = element_text(size = 13, face = "bold"),
      legend.position = "none"
    )
  
  if (nrow(stat_res) > 0) {
    p <- p + stat_pvalue_manual(
      stat_res,
      label = "p.value.label",
      tip.length = 0.03,
      bracket.size = 0.5,
      size = 4.0,
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position"
    )
  }
  p
}

# =============================================================================
# 4) Usage examples (commented)
# =============================================================================
# ---- Fig2 (SMC example; choose celltype4 explicitly) ----
# p_epi_N <- plot_dot_heatmap_two_platforms(
#   seu_nuc = epi.nuc.obj, seu_cell = epi.cell.obj,
#   nuc_celltype_col="celltype4", nuc_sample_col="orig.ident", nuc_tissue_col="TissueDiseaseState",
#   nuc_platform_col=NULL, nuc_platform_value=NULL,
#   cell_celltype_col="celltype4", cell_sample_col="orig.ident", cell_tissue_col="TissueDiseaseState",
#   cell_platform_col=NULL, cell_platform_value=NULL,
#   tissue_group="Normal",
#   alveolar_types=c("AT2","AT1 "),
#   bronchial_types=c("Secretory_Club","Ciliated","Secretory_Goblet")
# )
#
# ---- Fig2 (UL example; choose celltype3 explicitly; tissue is "tissue") ----
# p_epi_UL_N <- plot_dot_heatmap_two_platforms(
#   seu_nuc = epi.nuc.obj, seu_cell = epi.cell.obj,
#   nuc_celltype_col="celltype3", nuc_sample_col="orig.ident", nuc_tissue_col="tissue",
#   nuc_platform_col="platform", nuc_platform_value="NUC",
#   cell_celltype_col="celltype3", cell_sample_col="orig.ident", cell_tissue_col="tissue",
#   cell_platform_col="platform", cell_platform_value="CELL",
#   tissue_group="Normal",
#   alveolar_types=c("AT2","AT1 "),
#   bronchial_types=c("Secretory_Club","Ciliated","Secretory_Goblet")
# )
#
# ---- Sup Fig 3f-like frequency boxplot (explicit) ----
# out_freq <- freq_boxplot_stats_explicit(
#   seu = combined,
#   group_col = "celltype3",
#   sample_col = "orig.ident",
#   platform_col = "platform",
#   tissue_col = "TissueDiseaseState",  # or "tissue"
#   p_cut = 0.005
# )
# p_freq <- plot_freq_boxplot_explicit(out_freq)
# =============================================================================
