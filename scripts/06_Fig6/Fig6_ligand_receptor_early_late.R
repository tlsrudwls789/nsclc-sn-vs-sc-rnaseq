#!/usr/bin/env Rscript

############################################################
# Figure 6a-b: ligand-receptor analysis of epithelial Early/Late states
# (Early/Late labels come from scripts/11_SupFig7_early_late)
# ----------------------------------------------------------
# End-to-end LR analysis for epithelial Early/Late states
# across platforms (NUC vs CELL), using LIANA (Consensus).
#
# Inputs (via --args):
#   --lr_meta_rds        : Seurat object that contains merged cells + meta columns
#                          (or any Seurat containing the target cells you want)
#   --stage_nuc_rds      : Seurat/meta RDS containing Stage_EL for NUC cells
#   --stage_cell_rds     : Seurat/meta RDS containing Stage_like for CELL cells
#   --orig_pub_rds       : ORIGINAL Public Seurat to rebuild counts (Assay5 layers)
#   --orig_smc_rds       : (optional) ORIGINAL SMC Seurat to rebuild counts (may be missing)
#   --outdir             : output directory
#
# Required meta columns in lr_meta_rds:
#   orig.ident (for Platform inference if Platform missing),
#   celltype2 (receiver super-groups; e.g. B/T/NK/Myeloid/Fibroblast/Endothelial/Epithelial),
#   celltype  (optional; used for AT2)
#   CopyKAT   (optional; used for malignant proxy)
#
# Outputs:
#   - liana_rank_aggregate__NUC.tsv.gz / __CELL.tsv.gz
#   - delta_edges__NUC.tsv.gz / __CELL.tsv.gz  (Late - Early)
#   - Fig6-like summary heatmaps (PDF):
#       * heatmap_Fig6a_4combo.pdf  (NUC_Early/Late, CELL_Early/Late)
#       * heatmap_Fig6b_delta.pdf   (NUC: Late-Early, CELL: Late-Early)
#
# Notes:
#   - We compute LR for epithelial receiver states: EPI_{NUC/CELL}_{Early/Late}.
#   - Early/Late definition comes from Stage maps:
#       NUC: Stage_EL in {Early, Mid, Late} -> Stage2 uses only Early/Late
#       CELL: Stage_like {early_like, mid_like, late_like} -> Early/Mid/Late -> Stage2 Early/Late
#
############################################################

Sys.setenv(R_LIFECYCLE_VERBOSITY = "quiet")
options(lifecycle_verbosity = "quiet")
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
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

LR_META_RDS   <- get_arg("--lr_meta_rds")
STAGE_NUC_RDS <- get_arg("--stage_nuc_rds")
STAGE_CELL_RDS<- get_arg("--stage_cell_rds")
ORIG_PUB_RDS  <- get_arg("--orig_pub_rds")
ORIG_SMC_RDS  <- get_arg("--orig_smc_rds", default = NA_character_)
OUTDIR        <- get_arg("--outdir", default = "LR_out")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

STOP <- function(...) stop(paste0(...), call. = FALSE)

req_file <- function(x, nm) {
  if (is.na(x) || x == "" || !file.exists(x)) STOP("[missing] ", nm, ": ", x)
}
req_file(LR_META_RDS,   "lr_meta_rds")
req_file(STAGE_NUC_RDS, "stage_nuc_rds")
req_file(STAGE_CELL_RDS,"stage_cell_rds")
req_file(ORIG_PUB_RDS,  "orig_pub_rds")
if (!is.na(ORIG_SMC_RDS) && ORIG_SMC_RDS != "" && !file.exists(ORIG_SMC_RDS)) {
  message("[WARN] orig_smc_rds not found -> proceed with pub only: ", ORIG_SMC_RDS)
  ORIG_SMC_RDS <- NA_character_
}

# -----------------------------
# 1) Stage map utilities
# -----------------------------
.clean_stage_cell <- function(x) sub("^.*::", "", as.character(x))

make_stage_maps <- function(nuc_stage_rds, cell_stage_rds) {
  seu_nuc_stage  <- readRDS(nuc_stage_rds)
  seu_cell_stage <- readRDS(cell_stage_rds)
  
  if (!("Stage_EL" %in% colnames(seu_nuc_stage@meta.data)))  STOP("[stage] nuc: Stage_EL missing")
  if (!("Stage_like" %in% colnames(seu_cell_stage@meta.data))) STOP("[stage] cell: Stage_like missing")
  
  map_nuc <- seu_nuc_stage@meta.data %>%
    rownames_to_column("cell_raw") %>%
    mutate(cell = .clean_stage_cell(cell_raw)) %>%
    transmute(cell, Stage_NUC = as.character(Stage_EL)) %>%
    distinct(cell, .keep_all = TRUE)
  
  map_cell <- seu_cell_stage@meta.data %>%
    rownames_to_column("cell_raw") %>%
    mutate(cell = .clean_stage_cell(cell_raw)) %>%
    transmute(cell, Stage_like = as.character(Stage_like)) %>%
    mutate(Stage_CELL = dplyr::case_when(
      Stage_like == "early_like" ~ "Early",
      Stage_like == "mid_like"   ~ "Mid",
      Stage_like == "late_like"  ~ "Late",
      TRUE ~ NA_character_
    )) %>%
    select(cell, Stage_CELL) %>%
    distinct(cell, .keep_all = TRUE)
  
  list(map_nuc = map_nuc, map_cell = map_cell)
}

infer_platform <- function(md) {
  if ("Platform" %in% colnames(md)) return(md$Platform)
  if (!("orig.ident" %in% colnames(md))) return(rep(NA_character_, nrow(md)))
  dplyr::case_when(
    grepl("NUC",  md$orig.ident)  ~ "NUC",
    grepl("CELL", md$orig.ident)  ~ "CELL",
    TRUE ~ NA_character_
  )
}

attach_stage_group4 <- function(seu, map_nuc, map_cell) {
  md <- seu@meta.data %>% rownames_to_column("cell")
  md$Platform <- infer_platform(md)
  
  message(sprintf("[stage-join] overlap nuc=%d / cell=%d (total=%d)",
                  sum(md$cell %in% map_nuc$cell),
                  sum(md$cell %in% map_cell$cell),
                  nrow(md)))
  
  md <- md %>%
    left_join(map_nuc,  by="cell") %>%
    left_join(map_cell, by="cell") %>%
    mutate(
      Stage3 = dplyr::case_when(
        Platform=="NUC"  ~ Stage_NUC,
        Platform=="CELL" ~ Stage_CELL,
        TRUE ~ NA_character_),
      Stage2 = dplyr::case_when(
        Stage3 %in% c("Early","Late") ~ Stage3,
        TRUE ~ NA_character_),
      Group4 = dplyr::case_when(
        Platform=="NUC"  & Stage2=="Early" ~ "NUC_Early",
        Platform=="NUC"  & Stage2=="Late"  ~ "NUC_Late",
        Platform=="CELL" & Stage2=="Early" ~ "CELL_Early",
        Platform=="CELL" & Stage2=="Late"  ~ "CELL_Late",
        TRUE ~ NA_character_
      )
    )
  
  md$Group4 <- factor(md$Group4, levels=c("NUC_Early","NUC_Late","CELL_Early","CELL_Late"))
  
  rownames(md) <- md$cell
  md$cell <- NULL
  seu@meta.data <- md
  validObject(seu)
  seu
}

# -----------------------------
# 2) Assay/Assay5 counts extraction
# -----------------------------
.get_counts_layer <- function(seu, assay="RNA", layer_try=c("counts", "counts.NUC", "counts.CELL")) {
  DefaultAssay(seu) <- assay
  a <- seu[[assay]]
  
  # Seurat v5: Assay5 with layers
  if ("layers" %in% slotNames(a)) {
    avail <- Layers(a)
    for (ly in layer_try) {
      if (ly %in% avail) {
        M <- GetAssayData(seu, assay = assay, layer = ly)
        if (!inherits(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
        return(M)
      }
    }
    STOP("[counts] No requested layers in Assay5. Available: ", paste(avail, collapse=", "))
  }
  
  # Seurat v4: Assay@counts
  M <- GetAssayData(seu, assay=assay, slot="counts")
  if (!inherits(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
  M
}

# Assay5 dimnames restore helper (for NULL dimnames layers)
assay5_restore_layer_dimnames <- function(a, layer_name) {
  if (!inherits(a, "Assay5")) STOP("[Assay5] not Assay5")
  if (!(layer_name %in% names(a@layers))) STOP("[Assay5] missing layer: ", layer_name)
  
  M <- a@layers[[layer_name]]
  if (!inherits(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
  dn <- dimnames(M)
  if (!is.null(dn) && !is.null(dn[[1]]) && !is.null(dn[[2]])) return(M)
  
  feat_map <- a@features@.Data
  cell_map <- a@cells@.Data
  
  if (!(layer_name %in% colnames(feat_map))) STOP("[Assay5] features LogMap has no: ", layer_name)
  if (!(layer_name %in% colnames(cell_map))) STOP("[Assay5] cells LogMap has no: ", layer_name)
  
  feats_use <- rownames(feat_map)[which(feat_map[, layer_name])]
  cells_use <- rownames(cell_map)[which(cell_map[, layer_name])]
  
  if (nrow(M) != length(feats_use)) STOP("[Assay5] nrow mismatch: ", layer_name)
  if (ncol(M) != length(cells_use)) STOP("[Assay5] ncol mismatch: ", layer_name)
  
  dimnames(M) <- list(feats_use, cells_use)
  M
}

# gene union + cbind
.cbind_union_genes <- function(mats) {
  mats <- mats[!vapply(mats, is.null, logical(1))]
  if (length(mats) == 0) STOP("[cbind] no matrices")
  
  mats <- lapply(mats, function(M) {
    if (!inherits(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
    if (is.null(rownames(M)) || is.null(colnames(M))) STOP("[cbind] NULL dimnames")
    M
  })
  
  all_genes <- Reduce(base::union, lapply(mats, rownames))
  
  mats2 <- lapply(mats, function(M) {
    miss <- setdiff(all_genes, rownames(M))
    if (length(miss) > 0) {
      extra <- Matrix::Matrix(0, nrow=length(miss), ncol=ncol(M), sparse=TRUE)
      rownames(extra) <- miss
      colnames(extra) <- colnames(M)
      M <- rbind(M, extra)
    }
    M[all_genes, , drop=FALSE]
  })
  
  if (length(mats2) == 1) return(mats2[[1]])
  do.call(base::cbind, mats2)
}

.filter_genes_by_detect <- function(counts, min_pct=0.005) {
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  det <- Matrix::rowSums(counts > 0)
  keep <- det >= ceiling(ncol(counts) * min_pct)
  counts[keep, , drop=FALSE]
}

.filter_small_groups <- function(counts, md, idents_col="LR_CellGroup", min_cells=50) {
  tab <- table(md[[idents_col]], useNA="ifany")
  keep_groups <- names(tab)[tab >= min_cells]
  keep_cells  <- rownames(md)[as.character(md[[idents_col]]) %in% keep_groups]
  list(
    counts = counts[, keep_cells, drop=FALSE],
    md     = md[keep_cells, , drop=FALSE])
}

# -----------------------------
# 3) LIANA prep: define LR_CellGroup (EPI receivers staged)
# -----------------------------
.short_label <- function(x) {
  x <- as.character(x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) > 60,
         paste0(substr(x, 1, 35), "__", substr(x, nchar(x)-20, nchar(x))),
         x)
}

make_LR_CellGroup <- function(seu) {
  md <- seu@meta.data %>% rownames_to_column("cell")
  
  if (!("celltype2" %in% colnames(md))) STOP("[meta] need celltype2")
  sender_ct <- if ("celltype" %in% colnames(md)) "celltype" else "celltype2"
  
  md <- md %>%
    mutate(
      LR_CellGroup_raw = dplyr::case_when(
        as.character(celltype2) == "Epithelial" & !is.na(Group4) ~ paste0("EPI_", as.character(Group4)),
        as.character(celltype2) == "Epithelial" &  is.na(Group4) ~ "EPI_Unstaged",
        !is.na(.data[[sender_ct]]) ~ as.character(.data[[sender_ct]]),
        TRUE ~ "Other/NA"
      )
    )
  
  seu$LR_CellGroup_raw <- md$LR_CellGroup_raw[match(rownames(seu@meta.data), md$cell)]
  # shorten labels to avoid huge factor names
  lvl_raw  <- sort(unique(seu$LR_CellGroup_raw))
  lvl_safe <- .short_label(lvl_raw)
  lvl_code <- sprintf("G%04d_%s", seq_along(lvl_safe), lvl_safe)
  map_tbl  <- setNames(lvl_code, lvl_raw)
  
  seu$LR_CellGroup <- factor(unname(map_tbl[seu$LR_CellGroup_raw]))
  attr(seu@misc, "LR_CellGroup_map") <- data.frame(raw=lvl_raw, code=lvl_code)
  
  seu
}

# -----------------------------
# 4) Build LIANA-friendly Seurat (v4 assay) + Normalize
# -----------------------------
.build_liana_obj <- function(counts, md) {
  obj <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = md, assay = "RNA")
  obj[["RNA"]] <- SeuratObject::CreateAssayObject(counts = counts)  # enforce v4 assay
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, assay="RNA", verbose=FALSE)
  obj
}

# -----------------------------
# 5) Run LIANA
# -----------------------------
if (!requireNamespace("liana", quietly = TRUE)) {
  STOP(
    "[ERROR] liana not installed.\n",
    "Install:\n",
    "  install.packages('remotes')\n",
    "  remotes::install_github('saezlab/liana', upgrade='never')\n"
  )
}
suppressPackageStartupMessages(library(liana))

.run_liana <- function(obj, idents_col="LR_CellGroup", min_cells=50, expr_prop=0.10, mode=c("fast","full")) {
  mode <- match.arg(mode)
  methods_fast <- c("sca", "connectome", "logfc")
  methods_full <- c("cellphonedb", "sca", "connectome", "logfc")
  methods_use  <- if (mode == "full") methods_full else methods_fast
  
  liana_wrap(
    obj,
    method = methods_use, resource = "Consensus",
    idents_col = idents_col,min_cells = min_cells,
    expr_prop = expr_prop,verbose = TRUE
  )
}

# -----------------------------
# 6) Post-LIANA: Late–Early delta
# -----------------------------
.get_state <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    str_detect(x, "EPI_.*Early") ~ "Early",
    str_detect(x, "EPI_.*Late")  ~ "Late",
    TRUE ~ "Other"
  )
}
.is_epi_state <- function(x) {
  x <- as.character(x)
  str_detect(x, "EPI_") & str_detect(x, "Early|Late")
}

.compute_delta <- function(rank_df, score_col="sca.LRscore", rank_col="aggregate_rank") {
  stopifnot(all(c("source","target","ligand.complex","receptor.complex") %in% colnames(rank_df)))
  stopifnot(score_col %in% colnames(rank_df))
  stopifnot(rank_col  %in% colnames(rank_df))
  
  df2 <- rank_df %>%
    mutate(state = .get_state(target)) %>%
    filter(.is_epi_state(target)) %>%
    group_by(source, ligand.complex, receptor.complex, state) %>%
    summarise(
      score = mean(.data[[score_col]], na.rm = TRUE),
      rank  = mean(.data[[rank_col]],  na.rm = TRUE),
      .groups = "drop"
    )
  
  wide <- df2 %>%
    tidyr::pivot_wider(names_from=state, values_from=c(score, rank), values_fill=NA_real_)
  
  wide %>%
    mutate(
      delta_score = score_Late - score_Early,
      delta_rank  = rank_Late  - rank_Early,
      delta_rank_as_score = (rank_Early - rank_Late))
}

# -----------------------------
# 7) Fig6-like heatmaps (summary)
# -----------------------------
zscore_rows <- function(m) {
  m <- as.matrix(m)
  t(apply(m, 1, function(x) {
    x2 <- x; x2[!is.finite(x2)] <- NA_real_
    if (all(is.na(x2))) return(rep(NA_real_, length(x2)))
    mu <- mean(x2, na.rm=TRUE); sdv <- sd(x2, na.rm=TRUE)
    if (is.na(sdv) || sdv == 0) return(x2 - mu)
    (x2 - mu) / sdv
  }))
}

make_fig6_heatmaps <- function(delta_nuc, delta_cell, out_pdf_a, out_pdf_b,
                               topk_per_combo=10,
                               min_non_na=3,
                               celltype2_order=NULL) {
  
  # Build long edges: per platform, score_Early/score_Late -> Combo
  dnuc <- delta_nuc %>% mutate(platform="NUC")
  dcel <- delta_cell %>% mutate(platform="CELL")
  edges <- bind_rows(dnuc, dcel) %>%
    mutate(
      lr_pair = paste(ligand.complex, receptor.complex, sep="|"),
      # source is sender (stromal/immune etc.)
      sender = as.character(source)) %>%
    pivot_longer(cols=c(score_Early, score_Late), names_to="state", values_to="score") %>%
    mutate(
      Combo = paste0(platform, "_", ifelse(state=="score_Early","Early","Late")),
      Combo = factor(Combo, levels=c("NUC_Early","NUC_Late","CELL_Early","CELL_Late"))) %>%
    filter(is.finite(score))
  
  # choose celltype order by frequency unless provided
  if (is.null(celltype2_order)) {
    celltypes <- edges %>% count(sender, name="n") %>% arrange(desc(n)) %>% pull(sender)
  } else {
    celltypes <- celltype2_order
  }
  
  # select LR pairs per combo by max score
  lr_combo <- edges %>%
    group_by(lr_pair, Combo) %>%
    summarise(score_lr=max(score, na.rm=TRUE), .groups="drop")
  
  combos <- levels(edges$Combo)
  sel <- unique(unlist(lapply(combos, function(cb) {
    lr_combo %>%
      filter(Combo==cb, is.finite(score_lr)) %>%
      arrange(desc(score_lr)) %>%
      slice_head(n=topk_per_combo) %>%
      pull(lr_pair)
  })))
  
  # aggregate celltype-level max score
  agg <- edges %>%
    filter(lr_pair %in% sel) %>%
    group_by(lr_pair, Combo, sender) %>%
    summarise(max_score=max(score, na.rm=TRUE), .groups="drop")
  
  make_mat <- function(cb) {
    df <- agg %>% filter(Combo==cb) %>%
      pivot_wider(names_from=sender, values_from=max_score) %>% as.data.frame()
    rownames(df) <- df$lr_pair; df$lr_pair <- NULL
    for (ct in celltypes) if (!ct %in% colnames(df)) df[[ct]] <- NA_real_
    df <- as.matrix(df[, celltypes, drop=FALSE])
    keep <- rowSums(!is.na(df)) >= min_non_na
    df <- df[keep, , drop=FALSE]
    z <- zscore_rows(df)
    z[is.na(z)] <- 0
    list(raw=df, z=z)
  }
  
  mats <- lapply(combos, make_mat); names(mats) <- combos
  common_rows <- Reduce(intersect, lapply(mats, function(x) rownames(x$z)))
  mats <- lapply(mats, function(x) { x$z <- x$z[common_rows,,drop=FALSE]; x })
  
  col_fun <- circlize::colorRamp2(c(-2,0,2), c("#2C7BB6","#FFFFFF","#D7191C"))
  
  # Fig6a: 4 panels
  ht_list <- NULL
  for (i in seq_along(combos)) {
    cb <- combos[i]
    mz <- mats[[cb]]$z
    ht <- Heatmap(
      mz,
      name = if (i==1) "z" else paste0("z_",cb),
      col = col_fun,
      show_heatmap_legend = (i==1),
      column_title = cb,
      column_title_gp = gpar(fontface="bold", fontsize=12),
      show_column_names = TRUE, column_names_rot = 90,
      column_names_gp = gpar(fontsize=8),
      show_row_names = (i==1),
      row_labels = gsub("\\|"," - ", rownames(mz)),
      row_names_gp = gpar(fontsize=9, fontface="bold.italic"),
      cluster_rows = TRUE, cluster_columns = FALSE,
      use_raster = TRUE    )
    ht_list <- if (is.null(ht_list)) ht else (ht_list + ht)
  }
  
  pdf(out_pdf_a, width=14, height=max(7, 0.25*nrow(mats[[1]]$z)+3))
  draw(ht_list, heatmap_legend_side="right", ht_gap=unit(3,"mm"))
  dev.off()
  
  # Fig6b: delta heatmap (Late - Early) for NUC vs CELL
  # celltype axis: sender, value: delta_score (already edge-level; collapse by max per sender)
  delta2long <- function(df, platform) {
    df %>%
      mutate(lr_pair=paste(ligand.complex, receptor.complex, sep="|"),
             sender=as.character(source),
             platform=platform) %>%
      group_by(lr_pair, sender) %>%
      summarise(delta=max(delta_score, na.rm=TRUE), .groups="drop")
  }
  dn <- delta2long(delta_nuc, "NUC")
  dc <- delta2long(delta_cell,"CELL")
  
  mk_delta_mat <- function(d) {
    m <- d %>% filter(lr_pair %in% common_rows) %>%
      pivot_wider(names_from=sender, values_from=delta) %>% as.data.frame()
    rownames(m) <- m$lr_pair; m$lr_pair <- NULL
    for (ct in celltypes) if (!ct %in% colnames(m)) m[[ct]] <- NA_real_
    as.matrix(m[, celltypes, drop=FALSE])[common_rows,,drop=FALSE]
  }
  
  mN <- mk_delta_mat(dn)
  mC <- mk_delta_mat(dc)
  
  d_all <- c(as.numeric(mN), as.numeric(mC)); d_all <- d_all[is.finite(d_all)]
  d_lim <- as.numeric(quantile(abs(d_all), probs=0.95, na.rm=TRUE)); d_lim <- max(d_lim, 0.05)
  col_fun_d <- circlize::colorRamp2(c(-d_lim,0,d_lim), c("#2C7BB6","#FFFFFF","#D7191C"))
  
  htN <- Heatmap(mN, name="Δ", col=col_fun_d, column_title="NUC (Late - Early)",
                 show_row_names=FALSE, cluster_columns=FALSE, use_raster=TRUE)
  htC <- Heatmap(mC, name="Δ2", col=col_fun_d, column_title="CELL (Late - Early)",
                 show_heatmap_legend=FALSE,
                 row_labels=gsub("\\|"," - ", rownames(mC)),
                 row_names_side="right",
                 row_names_gp=gpar(fontsize=9, fontface="bold.italic"),
                 cluster_columns=FALSE, use_raster=TRUE)
  
  pdf(out_pdf_b, width=12, height=max(7, 0.25*nrow(mN)+3))
  draw(htN + htC, heatmap_legend_side="right", ht_gap=unit(4,"mm"))
  dev.off()
  
  invisible(list(common_rows=common_rows, celltypes=celltypes))
}

# -----------------------------
# 8) MAIN
# -----------------------------
message("=== Load LR meta object ===")
lr <- readRDS(LR_META_RDS)
if (!inherits(lr, "Seurat")) STOP("[lr_meta_rds] not a Seurat object")

need_cols <- c("celltype2","orig.ident")
miss <- setdiff(need_cols, colnames(lr@meta.data))
if (length(miss) > 0) STOP("[lr_meta_rds] missing meta columns: ", paste(miss, collapse=", "))

message("=== Build stage maps ===")
maps <- make_stage_maps(STAGE_NUC_RDS, STAGE_CELL_RDS)

message("=== Attach Stage2/Stage3 + Group4 ===")
lr <- attach_stage_group4(lr, maps$map_nuc, maps$map_cell)
message("[Group4]"); print(table(lr$Group4, useNA="ifany"))

message("=== Define LR_CellGroup (coded) ===")
lr <- make_LR_CellGroup(lr)

# Filter out Other/NA for LIANA
n0 <- ncol(lr)
lr <- subset(lr, subset = !is.na(LR_CellGroup_raw) & LR_CellGroup_raw != "Other/NA")
message("[filter] removed ", n0 - ncol(lr), " cells; remain=", ncol(lr))

# Split platform cells
cells_nuc  <- rownames(lr@meta.data)[as.character(lr$Platform)=="NUC"]
cells_cell <- rownames(lr@meta.data)[as.character(lr$Platform)=="CELL"]
message("[cells] NUC=", length(cells_nuc), " | CELL=", length(cells_cell))

message("=== Load ORIGINAL objects to rebuild counts ===")
pub <- readRDS(ORIG_PUB_RDS)
smc <- NULL
if (!is.na(ORIG_SMC_RDS)) smc <- readRDS(ORIG_SMC_RDS)

# If pub is Assay5 and layers have NULL dimnames, restore for counts.NUC/CELL
DefaultAssay(pub) <- "RNA"
if (inherits(pub[["RNA"]], "Assay5")) {
  a <- pub[["RNA"]]
  if ("counts.NUC" %in% names(a@layers)) a@layers[["counts.NUC"]]  <- assay5_restore_layer_dimnames(a, "counts.NUC")
  if ("counts.CELL"%in% names(a@layers)) a@layers[["counts.CELL"]] <- assay5_restore_layer_dimnames(a, "counts.CELL")
  pub@assays[["RNA"]] <- a
}

rebuild_counts_for_platform <- function(platform, cells_vec) {
  mats <- list()
  
  # PUB
  cells_pub <- intersect(cells_vec, colnames(pub))
  if (length(cells_pub) > 0) {
    ly <- if (platform=="NUC") c("counts.NUC","counts") else c("counts.CELL","counts")
    m <- .get_counts_layer(pub, assay="RNA", layer_try=ly)
    m <- m[, cells_pub, drop=FALSE]
    mats[[paste0("PUB_",platform)]] <- m
  }
  
  # SMC optional
  if (!is.null(smc)) {
    cells_smc <- intersect(cells_vec, colnames(smc))
    if (length(cells_smc) > 0) {
      ly <- if (platform=="NUC") c("counts.NUC","counts") else c("counts.CELL","counts")
      m <- .get_counts_layer(smc, assay="RNA", layer_try=ly)
      m <- m[, cells_smc, drop=FALSE]
      mats[[paste0("SMC_",platform)]] <- m
    }
  }
  
  counts <- .cbind_union_genes(mats)
  md <- lr@meta.data[colnames(counts), , drop=FALSE]
  list(counts=counts, md=md)
}

message("=== Rebuild counts: NUC ===")
xN <- rebuild_counts_for_platform("NUC", cells_nuc)
message("[counts NUC] ", paste(dim(xN$counts), collapse=" x "))

message("=== Rebuild counts: CELL ===")
xC <- rebuild_counts_for_platform("CELL", cells_cell)
message("[counts CELL] ", paste(dim(xC$counts), collapse=" x "))

# Filter tiny groups and low-detect genes
MIN_CELLS <- 50
EXPR_PROP <- 0.10
RUN_MODE  <- "fast"   # fast|full
DO_GENE_FILTER <- TRUE
MIN_PCT_DETECT <- 0.005

message("=== Filter small groups + genes ===")
fN <- .filter_small_groups(xN$counts, xN$md, idents_col="LR_CellGroup", min_cells=MIN_CELLS)
fC <- .filter_small_groups(xC$counts, xC$md, idents_col="LR_CellGroup", min_cells=MIN_CELLS)

if (DO_GENE_FILTER) {
  fN$counts <- .filter_genes_by_detect(fN$counts, min_pct=MIN_PCT_DETECT)
  fC$counts <- .filter_genes_by_detect(fC$counts, min_pct=MIN_PCT_DETECT)
}

message("=== Build LIANA input objects ===")
objN <- .build_liana_obj(fN$counts, fN$md)
objC <- .build_liana_obj(fC$counts, fC$md)

message("=== Run LIANA: NUC ===")
resN <- .run_liana(objN, idents_col="LR_CellGroup", min_cells=MIN_CELLS, expr_prop=EXPR_PROP, mode=RUN_MODE)

message("=== Run LIANA: CELL ===")
resC <- .run_liana(objC, idents_col="LR_CellGroup", min_cells=MIN_CELLS, expr_prop=EXPR_PROP, mode=RUN_MODE)

# Aggregate ranks
rankN <- liana::liana_aggregate(resN)
rankC <- liana::liana_aggregate(resC)

out_rankN <- file.path(OUTDIR, "liana_rank_aggregate__NUC.tsv.gz")
out_rankC <- file.path(OUTDIR, "liana_rank_aggregate__CELL.tsv.gz")
fwrite(rankN, out_rankN, sep="\t")
fwrite(rankC, out_rankC, sep="\t")
message("[saved] ", out_rankN)
message("[saved] ", out_rankC)

# Delta
message("=== Compute Late-Early deltas ===")
dN <- .compute_delta(rankN, score_col="sca.LRscore", rank_col="aggregate_rank")
dC <- .compute_delta(rankC, score_col="sca.LRscore", rank_col="aggregate_rank")

out_dN <- file.path(OUTDIR, "delta_edges__NUC.tsv.gz")
out_dC <- file.path(OUTDIR, "delta_edges__CELL.tsv.gz")
fwrite(dN, out_dN, sep="\t")
fwrite(dC, out_dC, sep="\t")
message("[saved] ", out_dN)
message("[saved] ", out_dC)

# Fig6-style heatmaps
message("=== Make Fig6-like heatmaps ===")
out_fig6a <- file.path(OUTDIR, "heatmap_Fig6a_4combo.pdf")
out_fig6b <- file.path(OUTDIR, "heatmap_Fig6b_delta.pdf")
make_fig6_heatmaps(dN, dC, out_fig6a, out_fig6b, topk_per_combo=10, min_non_na=3)

# Save mapping (coded LR_CellGroup -> raw)
map_df <- attr(lr@misc, "LR_CellGroup_map")
if (!is.null(map_df)) {
  saveRDS(map_df, file.path(OUTDIR, "LR_CellGroup_mapping_raw2code.rds"))
}

message("=== DONE ===")
message("Outputs written to: ", OUTDIR)
