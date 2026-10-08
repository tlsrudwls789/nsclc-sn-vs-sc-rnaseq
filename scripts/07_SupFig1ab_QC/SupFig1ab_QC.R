#!/usr/bin/env Rscript

############################################################
# SupFig1ab_QC.R
# ----------------------------------------------------------
# Supplementary Fig. 1a-b: DecontX contamination and per-cell QC
# metrics for snRNA-seq (NUC) and scRNA-seq (CELL).
#
#   a : median DecontX contamination per sample, by platform
#   b : violin plots (log scale) of nCount_RNA, nFeature_RNA and
#       percent.mt, by platform, with filtering thresholds
#
# Inputs (via --key=value; defaults resolve under DATA_DIR):
#   --nuc_meta         : tab-delimited per-cell metadata (snRNA-seq, before QC filter)
#   --cell_meta        : tab-delimited per-cell metadata (scRNA-seq, before QC filter)
#                        Required columns: orig.ident, nCount_RNA, nFeature_RNA, percent.mt
#   --nuc_decontx_dir  : folder with {sample}_DecontX_res.rds for snRNA-seq samples
#   --cell_decontx_dir : folder with {sample}_DecontX_res.rds for scRNA-seq samples
#                        Each RDS is a SingleCellExperiment with colData column `contamination`
#   --outdir           : output directory (default OUT_DIR, else "SupFig1_out")
#
# Outputs:
#   - SupFig1a_DecontX_contamination.pdf
#   - SupFig1b_QC_violin.pdf
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(SingleCellExperiment)
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

DATA_DIR <- Sys.getenv("DATA_DIR", unset = "data")
OUT_DIR  <- Sys.getenv("OUT_DIR",  unset = "SupFig1_out")

NUC_META         <- get_arg("--nuc_meta",
                            file.path(DATA_DIR, "LUNG_NUC.qc500_metadata_nUMI_nFeature_pMito.txt"))
CELL_META        <- get_arg("--cell_meta",
                            file.path(DATA_DIR, "LUNG_CELL.qc500_metadata_nUMI_nFeature_pMito.txt"))
NUC_DECONTX_DIR  <- get_arg("--nuc_decontx_dir",  file.path(DATA_DIR, "NUC_DecontX"))
CELL_DECONTX_DIR <- get_arg("--cell_decontx_dir", file.path(DATA_DIR, "CELL_DecontX"))
OUTDIR           <- get_arg("--outdir", OUT_DIR)

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# -----------------------------
# 1) Parameters
# -----------------------------
# Platform colors
col.nuc  <- "#1f77b4"   # snRNA-seq
col.cell <- "#ff7f0e"   # scRNA-seq
plat.cols <- c("snRNA-seq" = col.nuc, "scRNA-seq" = col.cell)

# QC thresholds
nuc.nFeature.min  <- 500
nuc.mito.max      <- 2
cell.nFeature.min <- 500
cell.mito.max     <- 20

# DecontX samples
nuc.sample.names  <- c("LUNG-N35-NUC", "LUNG-N36-NUC", "LUNG-N38-NUC", "LUNG-N39-NUC", "LUNG-N44-NUC",
                       "LUNG-T35-NUC", "LUNG-T36-NUC", "LUNG-T38-NUC", "LUNG-T39-NUC", "LUNG-T44-NUC")
cell.sample.names <- c("LUNG-N35-5PADT", "LUNG-N36-5PADT", "LUNG-N39-5PADT",
                       "LUNG-N43-5PADT", "LUNG-N44-5PADT",
                       "LUNG-T35-5PADT", "LUNG-T36-5PADT", "LUNG-T39-5PADT",
                       "LUNG-T41-5PADT", "LUNG-T43-5PADT", "LUNG-T44-5PADT")

save_pdf <- function(p, file, w, h) {
  dev <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(filename = file.path(OUTDIR, file), plot = p, width = w, height = h,
         units = "in", device = dev)
  message(">> Saved: ", file.path(OUTDIR, file))
}

# -----------------------------
# 2) Panel a: DecontX contamination per sample
# -----------------------------
load_decontx_summary <- function(sample.names, decontx.dir, platform.label) {
  res <- lapply(sample.names, function(s) {
    f <- file.path(decontx.dir, paste0(s, "_DecontX_res.rds"))
    if (!file.exists(f)) {
      message("  Not found: ", f)
      return(NULL)
    }
    contam <- readRDS(f)$contamination
    tibble(
      Sample               = gsub("-5PADT", "-CELL", s),
      Platform             = platform.label,
      n_cells_decontx      = length(contam),
      median_contamination = round(median(contam, na.rm = TRUE), 4),
      mean_contamination   = round(mean(contam, na.rm = TRUE), 4)
    )
  })
  bind_rows(Filter(Negate(is.null), res))
}

contam.table <- bind_rows(
  load_decontx_summary(nuc.sample.names,  NUC_DECONTX_DIR,  "snRNA-seq"),
  load_decontx_summary(cell.sample.names, CELL_DECONTX_DIR, "scRNA-seq")
)
stopifnot(nrow(contam.table) >= 2)

p_contam <- ggplot(contam.table, aes(x = Platform, y = median_contamination * 100,
                                     color = Platform)) +
  geom_jitter(width = 0.18, size = 3, alpha = 0.85) +
  geom_boxplot(width = 0.35, fill = NA, outlier.shape = NA, linewidth = 0.7) +
  scale_color_manual(values = plat.cols) +
  labs(x = NULL, y = "Median DecontX contamination per sample (%)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"))

save_pdf(p_contam, "SupFig1a_DecontX_contamination.pdf", w = 2.8, h = 3.5)

# -----------------------------
# 3) Panel b: per-cell QC violins (all cells before QC filter)
# -----------------------------
stopifnot(file.exists(NUC_META), file.exists(CELL_META))
nuc.meta  <- read.table(NUC_META,  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cell.meta <- read.table(CELL_META, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

qc.df <- bind_rows(
  nuc.meta  %>% select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt) %>%
    mutate(Platform = "snRNA-seq"),
  cell.meta %>% select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt) %>%
    mutate(Platform = "scRNA-seq")
) %>%
  mutate(Platform = factor(Platform, levels = c("snRNA-seq", "scRNA-seq")))

violin_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text       = element_text(color = "black", size = base_size - 1),
      axis.title.y    = element_text(size = base_size - 1),
      axis.title.x    = element_blank(),
      legend.position = "none",
      plot.title      = element_text(size = base_size, face = "bold", hjust = 0.5)
    )
}

# nCount_RNA (log10, +1 offset)
p_nCount <- ggplot(qc.df %>% mutate(val = nCount_RNA + 1),
                   aes(x = Platform, y = val, fill = Platform)) +
  geom_violin(width = 0.9, alpha = 0.85, trim = TRUE, color = "white") +
  geom_boxplot(width = 0.07, fill = "white", outlier.shape = NA, color = "gray30") +
  scale_y_log10(
    breaks = c(100, 500, 1000, 5000, 10000, 50000, 100000, 500000),
    labels = c("100", "500", "1K", "5K", "10K", "50K", "100K", "500K")
  ) +
  scale_fill_manual(values = plat.cols) +
  labs(y = "nCount_RNA (log scale)") +
  violin_theme()

# nFeature_RNA (log10, +1 offset), threshold at 500
p_nFeature <- ggplot(qc.df %>% mutate(val = nFeature_RNA + 1),
                     aes(x = Platform, y = val, fill = Platform)) +
  geom_violin(width = 0.9, alpha = 0.85, trim = TRUE, color = "white") +
  geom_boxplot(width = 0.07, fill = "white", outlier.shape = NA, color = "gray30") +
  geom_hline(yintercept = nuc.nFeature.min + 1, linetype = "dashed",
             color = "gray40", linewidth = 0.7) +
  scale_y_log10(
    breaks = c(100, 200, 500, 1000, 2000, 5000, 10000),
    labels = c("100", "200", "500", "1K", "2K", "5K", "10K")
  ) +
  scale_fill_manual(values = plat.cols) +
  labs(y = "nFeature_RNA (log scale)") +
  violin_theme() +
  annotate("text", x = 0.52, y = (nuc.nFeature.min + 1) * 1.22,
           label = "500", size = 3, color = "gray40", hjust = 0)

# percent.mt (log10, +0.01 offset), platform-specific thresholds
p_mito <- ggplot(qc.df %>% mutate(val = percent.mt + 0.01),
                 aes(x = Platform, y = val, fill = Platform)) +
  geom_violin(width = 0.9, alpha = 0.85, trim = TRUE, color = "white") +
  geom_boxplot(width = 0.07, fill = "white", outlier.shape = NA, color = "gray30") +
  geom_hline(yintercept = nuc.mito.max  + 0.01, linetype = "dashed",
             color = col.nuc,  linewidth = 0.7) +
  geom_hline(yintercept = cell.mito.max + 0.01, linetype = "dashed",
             color = col.cell, linewidth = 0.7) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 2, 5, 10, 20, 50, 100),
    labels = c("0", "0.1", "1", "2", "5", "10", "20", "50", "100")
  ) +
  scale_fill_manual(values = plat.cols) +
  labs(y = "Mitochondrial fraction (%, log scale)") +
  violin_theme() +
  annotate("text", x = 0.52, y = (nuc.mito.max  + 0.01) * 1.5,
           label = "2%",  size = 3, color = col.nuc,  hjust = 0) +
  annotate("text", x = 0.52, y = (cell.mito.max + 0.01) * 1.3,
           label = "20%", size = 3, color = col.cell, hjust = 0)

# Shared platform legend
legend_plot <- ggplot(qc.df, aes(x = Platform, y = nFeature_RNA, fill = Platform)) +
  geom_violin() +
  scale_fill_manual(values = plat.cols, name = "Platform") +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_text(size = 10, face = "bold"),
    legend.text       = element_text(size = 10),
    legend.key.size   = unit(0.4, "cm"),
    legend.background = element_blank()
  )
shared_legend <- get_legend(legend_plot)

panel_row <- plot_grid(p_nCount, p_nFeature, p_mito,
                       nrow = 1, ncol = 3, align = "h", axis = "tb",
                       rel_widths = c(1, 1, 1))
fig_b <- plot_grid(panel_row, shared_legend, ncol = 1, rel_heights = c(1, 0.08))

save_pdf(fig_b, "SupFig1b_QC_violin.pdf", w = 7.5, h = 3.8)

message("Done. Outputs in: ", OUTDIR)
