#!/usr/bin/env Rscript

############################################################
# SupFig5_cell_lines_GSE303762.R
# ----------------------------------------------------------
# Supplementary Fig. 5: platform-associated differences and LUAD-vs-SCLC
# discrimination in an independent lung cancer cell-line dataset
# (LongBench, GSE303762; paired snRNA-seq = NUC and scRNA-seq = CELL).
#
# NOTE: Re-implementation written after publication following the Methods of
# the article; the original figure-generation script was not retained.
# Parameters not stated in the Methods are marked 'assumed'. Results may
# differ slightly from the published figure.
#
# Inputs (via --key=value):
#   --nuc_counts   : NUC count matrix (10x directory or .h5)
#   --cell_counts  : CELL count matrix (10x directory or .h5)
#   --meta         : metadata table (.csv/.tsv) with one row per barcode and
#                    columns given by --barcode_col, --platform_col (NUC/CELL),
#                    --cellline_col (Vireo donor / cell line), --type_col (LUAD/SCLC)
#   --outdir       : output directory (default: SupFig5_out)
#   --top_n        : GOBP terms shown per side in c and e (default: 10; assumed)
#
# Barcodes in --meta must match the column names of each count matrix
# (unprefixed); the platform column is used to pair them.
#
# Outputs:
#   SupFig5a ... SupFig5f PDFs
#   SupData5_*.csv : DEG tables, GOBP tables, concordance summary
############################################################

Sys.setenv(R_LIFECYCLE_VERBOSITY = "quiet")
options(lifecycle_verbosity = "quiet")
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(VennDiagram)
  library(grid)
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

NUC_COUNTS   <- get_arg("--nuc_counts")
CELL_COUNTS  <- get_arg("--cell_counts")
META_FILE    <- get_arg("--meta")
OUTDIR       <- get_arg("--outdir", default = "SupFig5_out")
BARCODE_COL  <- get_arg("--barcode_col",  default = "barcode")
PLATFORM_COL <- get_arg("--platform_col", default = "platform")
CELLLINE_COL <- get_arg("--cellline_col", default = "cell_line")
TYPE_COL     <- get_arg("--type_col",     default = "cancer_type")
TOP_N        <- as.integer(get_arg("--top_n", default = "10"))   # assumed

for (p in c(NUC_COUNTS, CELL_COUNTS, META_FILE)) {
  if (is.na(p) || !file.exists(p)) stop("Missing required input: ", p)
}
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(1)  # assumed

COL_PLAT <- c(NUC = "#1F77B4", CELL = "#FF7F0E")
COL_TYPE <- c(LUAD = "#7B3294", SCLC = "#D62728")
LUAD_LINES <- c("HCC827", "NCI-H1975", "NCI-H2228")
SCLC_LINES <- c("NCI-H211", "NCI-H526", "NCI-H146", "NCI-H69", "SHP-77")

# Parameters
N_PCS     <- 30     # assumed
LOGFC_THR <- 0.25   # Methods
MIN_PCT   <- 0.1    # Methods
PADJ_THR  <- 0.05   # Methods (BH)
MIN_GS    <- 10     # Methods
FDR_CAP   <- 300    # legend (visualization only)

# -----------------------------
# 1) Load, filter Vireo-unassigned cells, merge
# -----------------------------
read_counts <- function(path) {
  m <- if (grepl("\\.h5$", path, ignore.case = TRUE)) Read10X_h5(path) else Read10X(path)
  if (is.list(m)) m <- m[["Gene Expression"]]  # multi-modal h5: keep RNA
  m
}

meta <- as.data.frame(fread(META_FILE))
meta <- data.frame(
  barcode   = as.character(meta[[BARCODE_COL]]),
  platform  = toupper(as.character(meta[[PLATFORM_COL]])),
  cell_line = as.character(meta[[CELLLINE_COL]]),
  cancer_type = toupper(as.character(meta[[TYPE_COL]])),
  stringsAsFactors = FALSE
)
# Cells lacking a confident Vireo assignment are excluded (Methods)
bad <- is.na(meta$cell_line) | meta$cell_line == "" |
  tolower(meta$cell_line) %in% c("unassigned", "doublet", "na")
meta <- meta[!bad & meta$platform %in% c("NUC", "CELL") &
               meta$cancer_type %in% c("LUAD", "SCLC"), ]
meta$cell_id <- paste(meta$platform, meta$barcode, sep = "_")
cat(sprintf("[meta] %d confidently assigned cells retained\n", nrow(meta)))

make_obj <- function(path, plat) {
  m <- read_counts(path)
  colnames(m) <- paste(plat, colnames(m), sep = "_")
  mm <- meta[meta$platform == plat, ]
  keep <- intersect(colnames(m), mm$cell_id)
  if (length(keep) == 0) stop("No barcodes matched metadata for ", plat)
  o <- CreateSeuratObject(m[, keep], project = plat)
  rownames(mm) <- mm$cell_id
  AddMetaData(o, mm[keep, c("platform", "cell_line", "cancer_type")])
}

obj <- merge(make_obj(NUC_COUNTS, "NUC"), make_obj(CELL_COUNTS, "CELL"))
obj <- JoinLayers(obj)  # assumed: single joined RNA layer (Seurat v5)
obj$platform    <- factor(obj$platform, levels = c("NUC", "CELL"))
obj$cancer_type <- factor(obj$cancer_type, levels = c("LUAD", "SCLC"))
obj$cell_line   <- factor(obj$cell_line,
                          levels = unique(c(LUAD_LINES, SCLC_LINES, sort(unique(obj$cell_line)))))
obj$cell_line   <- droplevels(obj$cell_line)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
print(table(obj$cell_line, obj$platform))

# -----------------------------
# 2) Normalize, PCA, Harmony, UMAP (standard Seurat defaults; assumed)
# -----------------------------
obj <- NormalizeData(obj, verbose = FALSE)                       # LogNormalize, 1e4 (assumed)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)  # vst, 2000 (assumed)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = N_PCS, verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "platform", verbose = FALSE)  # Harmony defaults (assumed)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)

# -----------------------------
# Helpers: DE and GOBP
# -----------------------------
run_de <- function(o, group_col, id1, id2) {
  Idents(o) <- group_col
  de <- FindMarkers(o, ident.1 = id1, ident.2 = id2, test.use = "wilcox",
                    logfc.threshold = LOGFC_THR, min.pct = MIN_PCT, verbose = FALSE)
  de$gene <- rownames(de)
  de$sig  <- de$p_val_adj < PADJ_THR & abs(de$avg_log2FC) >= LOGFC_THR
  de
}

run_go <- function(genes) {
  if (length(genes) < MIN_GS) return(data.frame())
  ego <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = PADJ_THR,
                  minGSSize = MIN_GS)  # universe, maxGSSize, qvalueCutoff: defaults (assumed)
  if (is.null(ego)) return(data.frame())
  as.data.frame(ego)
}

go_pair <- function(de, pos_label, neg_label) {
  up_pos <- de$gene[de$sig & de$avg_log2FC > 0]
  up_neg <- de$gene[de$sig & de$avg_log2FC < 0]
  bind_rows(mutate(run_go(up_pos), side = pos_label),
            mutate(run_go(up_neg), side = neg_label))
}

# Signed bar plot: top_n terms per side, value = sign * -log10(p.adjust)
signed_bar <- function(go, sign_pos, title, fills) {
  if (nrow(go) == 0) return(ggplot() + ggtitle(paste(title, "(no terms)")))
  d <- go %>% group_by(side) %>% slice_min(p.adjust, n = TOP_N, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(score = ifelse(side == sign_pos, 1, -1) * -log10(p.adjust),
           term = factor(Description, levels = unique(Description[order(score)])))
  ggplot(d, aes(score, term, fill = side)) + geom_col(width = 0.75) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = fills) +
    labs(x = "Signed -log10(adjusted p)", y = NULL, fill = NULL, title = title) +
    theme_classic(base_size = 9)
}

save_csv <- function(df, name) fwrite(df, file.path(OUTDIR, name))

# -----------------------------
# a) Per-cell QC metrics
# -----------------------------
qc <- FetchData(obj, c("platform", "cancer_type", "nFeature_RNA", "nCount_RNA", "percent.mt")) %>%
  pivot_longer(c(nFeature_RNA, nCount_RNA, percent.mt), names_to = "metric") %>%
  mutate(metric = factor(metric, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                         c("Detected genes", "UMI counts", "Mitochondrial %")),
         value = log10(value + 1))  # log-scaled; +1 pseudocount (assumed)
pa <- ggplot(qc, aes(platform, value, fill = platform)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  facet_grid(metric ~ cancer_type, scales = "free_y") +
  scale_fill_manual(values = COL_PLAT) +
  labs(x = NULL, y = "log10(value + 1)") + theme_classic(base_size = 9) +
  theme(legend.position = "none")
ggsave(file.path(OUTDIR, "SupFig5a_QC_violin.pdf"), pa, width = 4, height = 6)

# -----------------------------
# b) Harmony-integrated UMAPs
# -----------------------------
pb <- DimPlot(obj, group.by = "platform", cols = COL_PLAT, raster = FALSE) |
  DimPlot(obj, group.by = "cancer_type", cols = COL_TYPE, raster = FALSE) |
  DimPlot(obj, group.by = "cell_line", raster = FALSE)
ggsave(file.path(OUTDIR, "SupFig5b_UMAP_harmony.pdf"), pb, width = 15, height = 4.5)

# -----------------------------
# c) Platform-associated DEGs and GOBP, within LUAD and within SCLC
#    (CELL vs NUC; positive = CELL)
# -----------------------------
pc <- list()
for (ct in c("LUAD", "SCLC")) {
  sub <- subset(obj, subset = cancer_type == ct)
  de  <- run_de(sub, "platform", "CELL", "NUC")
  save_csv(de, sprintf("SupData5_platformDEG_%s_CELLvsNUC.csv", ct))
  go  <- go_pair(de, "CELL", "NUC")
  save_csv(go, sprintf("SupData5_platformGOBP_%s.csv", ct))
  cat(sprintf("[c] %s: CELL-up %d, NUC-up %d DEGs\n", ct,
              sum(de$sig & de$avg_log2FC > 0), sum(de$sig & de$avg_log2FC < 0)))
  pc[[ct]] <- signed_bar(go, "CELL", paste0(ct, ": CELL vs NUC"), COL_PLAT)
}
ggsave(file.path(OUTDIR, "SupFig5c_platform_GOBP.pdf"),
       pc$LUAD / pc$SCLC, width = 7, height = 8)

# -----------------------------
# d) LUAD vs SCLC DE within each platform (positive = LUAD)
# -----------------------------
de_ct <- list()
pd <- list()
for (pl in c("NUC", "CELL")) {
  sub <- subset(obj, subset = platform == pl)
  de  <- run_de(sub, "cancer_type", "LUAD", "SCLC")
  de_ct[[pl]] <- de
  save_csv(de, sprintf("SupData5_cancertypeDEG_%s_LUADvsSCLC.csv", pl))
  de <- de %>% mutate(
    nlp = pmin(-log10(p_val_adj), FDR_CAP),  # p_val_adj = 0 maps to the cap
    cls = case_when(sig & avg_log2FC > 0 ~ "LUAD", sig & avg_log2FC < 0 ~ "SCLC", TRUE ~ "NS"))
  pd[[pl]] <- ggplot(de, aes(avg_log2FC, nlp, color = cls)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c(COL_TYPE, NS = "grey75")) +
    geom_hline(yintercept = -log10(PADJ_THR), linetype = 2, linewidth = 0.3) +
    geom_vline(xintercept = c(-LOGFC_THR, LOGFC_THR), linetype = 2, linewidth = 0.3) +
    labs(x = "log2FC (LUAD vs SCLC)", y = "-log10(FDR), capped at 300", color = NULL,
         title = sprintf("%s: LUAD-up %d, SCLC-up %d", pl,
                         sum(de$cls == "LUAD"), sum(de$cls == "SCLC"))) +
    theme_classic(base_size = 9)
}
ggsave(file.path(OUTDIR, "SupFig5d_volcano_LUADvsSCLC.pdf"),
       pd$NUC | pd$CELL, width = 9, height = 4)

# -----------------------------
# e) GOBP of LUAD-vs-SCLC DEGs (positive = SCLC)
# -----------------------------
go_ct <- list()
pe <- list()
for (pl in c("NUC", "CELL")) {
  go <- go_pair(de_ct[[pl]], "LUAD", "SCLC")
  go_ct[[pl]] <- go
  save_csv(go, sprintf("SupData5_cancertypeGOBP_%s.csv", pl))
  pe[[pl]] <- signed_bar(go, "SCLC", paste0(pl, ": SCLC vs LUAD"), COL_TYPE)
}
ggsave(file.path(OUTDIR, "SupFig5e_cancertype_GOBP.pdf"),
       pe$NUC / pe$CELL, width = 7, height = 8)

# -----------------------------
# f) Cross-platform overlap and concordance
#    concordance = |shared same-direction features| / min(|NUC set|, |CELL set|)
# -----------------------------
sets <- list(
  "LUAD-up genes" = lapply(de_ct, function(d) d$gene[d$sig & d$avg_log2FC > 0]),
  "SCLC-up genes" = lapply(de_ct, function(d) d$gene[d$sig & d$avg_log2FC < 0]),
  "LUAD-up GOBP"  = lapply(go_ct, function(g) if (nrow(g)) g$ID[g$side == "LUAD"] else character()),
  "SCLC-up GOBP"  = lapply(go_ct, function(g) if (nrow(g)) g$ID[g$side == "SCLC"] else character())
)

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
conc <- list()
pdf(file.path(OUTDIR, "SupFig5f_venn_concordance.pdf"), width = 8, height = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
for (i in seq_along(sets)) {
  nm <- names(sets)[i]
  a <- unique(sets[[nm]]$NUC); b <- unique(sets[[nm]]$CELL)
  sh <- length(intersect(a, b))
  pct <- if (min(length(a), length(b)) > 0) 100 * sh / min(length(a), length(b)) else NA
  conc[[nm]] <- data.frame(set = nm, NUC_only = length(a) - sh, shared = sh,
                           CELL_only = length(b) - sh, NUC_total = length(a),
                           CELL_total = length(b), concordance_pct = round(pct, 1))
  pushViewport(viewport(layout.pos.row = (i - 1) %/% 2 + 1, layout.pos.col = (i - 1) %% 2 + 1))
  g <- draw.pairwise.venn(length(a), length(b), sh, category = c("NUC", "CELL"),
                          fill = unname(COL_PLAT), alpha = 0.5, lwd = 1, cex = 1,
                          cat.cex = 1, scaled = TRUE, ind = FALSE, margin = 0.1)
  grid.draw(g)
  grid.text(sprintf("%s\nconcordance = %.1f%%", nm, pct), y = 0.95,
            gp = gpar(fontsize = 10))
  popViewport()
}
dev.off()

conc <- bind_rows(conc)
save_csv(conc, "SupData5_concordance_summary.csv")
cat("\n[f] Cross-platform concordance (published: 58.7 / 91.4 / 62.4 / 77.6 %)\n")
print(conc, row.names = FALSE)

saveRDS(obj, file.path(OUTDIR, "SupFig5_GSE303762_harmony.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
cat("[done] outputs in ", OUTDIR, "\n")
