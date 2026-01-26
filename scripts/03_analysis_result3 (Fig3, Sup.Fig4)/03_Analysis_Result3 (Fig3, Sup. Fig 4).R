################################################################################
## 02_Analysis_Result2_GitHub_single_script.R
## Purpose:
## - In-house integrated dataset only
## - DEG (single-cell) + pseudo-bulk (sn vs sc) across 5 shared cell types
## - GO(BP) for each shared cell type, separately for sn(NUC) and sc(CELL)
## - Jaccard heatmap of GO terms (within each platform + optional cross-platform)
##
## Shared5 (fixed):
##   Endothelia_vascular_arterial, Endothelia_vascular_venous, AT1, AT2, Fibro_adventitial
##
## Assumptions (DO NOT auto-detect; explicit columns):
##   - platform column: "platform" with values c("NUC","CELL")
##   - tissue state column: "TissueDiseaseState" with values c("N","T")
##   - celltype column (fine): "celltype4" containing the shared5 labels
##   - sample identifier column: "orig.ident"
##   - pseudobulk unit column: "orig.ident" (sample-level)
##
## Usage:
##   source("02_Analysis_Result2_GitHub_single_script.R")
##   res <- run_result2_shared5_pipeline(
##     seu = seu,
##     shared5 = SHARED5,
##     tissue_keep = "N"   # Epithelial-normal-only constraint does NOT apply here
##   )
##   res$plots$jaccard_go_nuc
##   res$plots$jaccard_go_cell
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  library(ComplexHeatmap)
  library(circlize)
})

# -----------------------------
# 0) Fixed configuration
# -----------------------------
PLATFORM_COL <- "platform"              # NUC / CELL
TISSUE_COL   <- "TissueDiseaseState"    # N / T
CELLTYPE_COL <- "celltype4"             # fine celltype
SAMPLE_COL   <- "orig.ident"            # per-sample id
ASSAY_USE    <- "RNA"

PLATFORM_LEVELS <- c("NUC", "CELL")

SHARED5 <- c(
  "Endothelia_vascular_arterial",
  "Endothelia_vascular_venous",
  "AT1",
  "AT2",
  "Fibro_adventitial"
)

# -----------------------------
# 1) Small helpers
# -----------------------------
.assert_cols <- function(seu, cols) {
  md <- seu@meta.data
  missing <- setdiff(cols, colnames(md))
  if (length(missing) > 0) {
    stop("Missing required meta.data columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

.trim_levels <- function(x) {
  # normalize minor whitespace inconsistencies WITHOUT guessing columns
  if (is.character(x)) trimws(x) else x
}

subset_by_tissue_platform <- function(seu, tissue_keep = c("N","T"), platform_keep = c("NUC","CELL")) {
  .assert_cols(seu, c(PLATFORM_COL, TISSUE_COL))
  tissue_keep   <- unique(tissue_keep)
  platform_keep <- unique(platform_keep)
  subset(seu, subset = .data[[TISSUE_COL]] %in% tissue_keep & .data[[PLATFORM_COL]] %in% platform_keep)
}

subset_to_shared5 <- function(seu, shared5 = SHARED5) {
  .assert_cols(seu, CELLTYPE_COL)
  md <- seu@meta.data
  md[[CELLTYPE_COL]] <- .trim_levels(as.character(md[[CELLTYPE_COL]]))
  seu@meta.data[[CELLTYPE_COL]] <- md[[CELLTYPE_COL]]
  subset(seu, subset = .data[[CELLTYPE_COL]] %in% shared5)
}

# -----------------------------
# 2) (A) Single-cell DEG: within platform (NUC or CELL)
#    celltype vs rest (within shared5 only)
# -----------------------------
deg_singlecell_celltype_vs_rest <- function(
    seu,
    shared5 = SHARED5,
    platform = c("NUC","CELL"),
    assay = ASSAY_USE,
    slot = "data",
    test_use = "wilcox",
    min_pct = 0.1,
    logfc_threshold = 0.25,
    min_cells_ct = 30,
    min_cells_rest = 30
) {
  platform <- match.arg(platform, choices = PLATFORM_LEVELS)
  .assert_cols(seu, c(PLATFORM_COL, CELLTYPE_COL))
  
  DefaultAssay(seu) <- assay
  
  seu_pf <- subset(seu, subset = .data[[PLATFORM_COL]] == platform)
  seu_pf <- subset_to_shared5(seu_pf, shared5 = shared5)
  if (ncol(seu_pf) == 0) return(list())
  
  Idents(seu_pf) <- factor(seu_pf@meta.data[[CELLTYPE_COL]], levels = shared5)
  
  out <- list()
  for (ct in shared5) {
    n_ct   <- sum(Idents(seu_pf) == ct)
    n_rest <- sum(Idents(seu_pf) != ct)
    if (n_ct < min_cells_ct || n_rest < min_cells_rest) next
    
    deg <- FindMarkers(
      seu_pf,
      ident.1 = ct,
      ident.2 = setdiff(shared5, ct),
      test.use = test_use,
      min.pct = min_pct,
      logfc.threshold = logfc_threshold,
      slot = slot
    ) |>
      tibble::rownames_to_column("gene") |>
      dplyr::as_tibble() |>
      dplyr::mutate(
        platform = platform,
        celltype = ct,
        contrast = "celltype_vs_rest (within shared5; single-cell)"
      )
    
    out[[ct]] <- deg
  }
  out
}

# -----------------------------
# 3) (B) Pseudo-bulk DEG: sn vs sc within each celltype
#    - Aggregate to sample-level (orig.ident) average expression per gene
#    - Then do per-gene Wilcoxon test across samples: NUC vs CELL
# -----------------------------
pseudobulk_matrix <- function(
    seu,
    assay = ASSAY_USE,
    slot = "data",
    celltype,
    shared5 = SHARED5,
    tissue_keep = "N"
) {
  .assert_cols(seu, c(SAMPLE_COL, PLATFORM_COL, TISSUE_COL, CELLTYPE_COL))
  DefaultAssay(seu) <- assay
  
  seu2 <- subset(seu, subset = .data[[TISSUE_COL]] %in% tissue_keep & .data[[CELLTYPE_COL]] == celltype)
  if (ncol(seu2) == 0) return(NULL)
  
  md <- seu2@meta.data
  md[[SAMPLE_COL]]   <- as.character(md[[SAMPLE_COL]])
  md[[PLATFORM_COL]] <- as.character(md[[PLATFORM_COL]])
  
  # group key = sample-platform
  group_key <- paste(md[[SAMPLE_COL]], md[[PLATFORM_COL]], sep = "||")
  groups <- unique(group_key)
  
  mat <- GetAssayData(seu2, assay = assay, slot = slot)
  # average expression per group
  avg_list <- lapply(groups, function(g) {
    cols <- which(group_key == g)
    if (length(cols) == 1) as.matrix(mat[, cols, drop = FALSE])
    else Matrix::rowMeans(mat[, cols, drop = FALSE])
  })
  avg_mat <- do.call(cbind, lapply(avg_list, function(x) if (is.matrix(x)) x[,1] else x))
  colnames(avg_mat) <- groups
  avg_mat
}

deg_pseudobulk_sn_vs_sc <- function(
    seu,
    celltype,
    tissue_keep = "N",
    assay = ASSAY_USE,
    slot = "data",
    p_adj_method = "BH",
    lfc_eps = 1e-9
) {
  avg_mat <- pseudobulk_matrix(
    seu = seu,
    assay = assay,
    slot = slot,
    celltype = celltype,
    tissue_keep = tissue_keep
  )
  if (is.null(avg_mat)) return(NULL)
  
  # parse meta from columns: sample||platform
  col_df <- tibble::tibble(
    group = colnames(avg_mat),
    sample = sub("\\|\\|.*$", "", colnames(avg_mat)),
    platform = sub("^.*\\|\\|", "", colnames(avg_mat))
  )
  
  # require both platforms with >=2 samples each for a robust Wilcoxon
  nuc_cols  <- which(col_df$platform == "NUC")
  cell_cols <- which(col_df$platform == "CELL")
  if (length(nuc_cols) < 2 || length(cell_cols) < 2) return(NULL)
  
  # per gene test
  genes <- rownames(avg_mat)
  pvals <- numeric(length(genes))
  lfc   <- numeric(length(genes))
  
  for (i in seq_along(genes)) {
    x <- as.numeric(avg_mat[i, nuc_cols])
    y <- as.numeric(avg_mat[i, cell_cols])
    
    # pseudo-bulk logFC: mean(NUC) - mean(CELL)
    lfc[i] <- (mean(x) + lfc_eps) - (mean(y) + lfc_eps)
    
    pvals[i] <- tryCatch(
      wilcox.test(x = x, y = y, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  }
  
  tibble::tibble(
    gene = genes,
    p_val = pvals,
    p_val_adj = p.adjust(pvals, method = p_adj_method),
    avg_diff = lfc,
    celltype = celltype,
    contrast = "NUC_vs_CELL (pseudo-bulk; sample-level)"
  ) |>
    dplyr::arrange(p_val_adj)
}

# -----------------------------
# 4) GO(BP) enrichment from DEG table
# -----------------------------
go_from_deg <- function(
    deg_df,
    gene_col = "gene",
    p_col = "p_val_adj",
    lfc_col = NULL,            # single-cell: "avg_log2FC", pseudobulk: "avg_diff"
    p_cut = 0.005,             # <-- your request: 0.115 -> 0.005
    lfc_cut = 0.25,
    direction = c("up","down","both","none"),
    ont = "BP",
    max_terms = 30
) {
  direction <- match.arg(direction)
  if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)
  if (!all(c(gene_col, p_col) %in% colnames(deg_df))) return(NULL)
  
  df2 <- deg_df |>
    dplyr::filter(!is.na(.data[[p_col]]), .data[[p_col]] <= p_cut)
  
  if (!is.null(lfc_col) && lfc_col %in% colnames(df2) && direction != "none") {
    df2 <- df2 |> dplyr::filter(!is.na(.data[[lfc_col]]))
    if (direction == "up")   df2 <- df2 |> dplyr::filter(.data[[lfc_col]] >=  lfc_cut)
    if (direction == "down") df2 <- df2 |> dplyr::filter(.data[[lfc_col]] <= -lfc_cut)
    if (direction == "both") df2 <- df2 |> dplyr::filter(abs(.data[[lfc_col]]) >= lfc_cut)
  }
  
  genes <- unique(df2[[gene_col]])
  if (length(genes) < 10) return(NULL)
  
  map <- suppressMessages(
    clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  if (is.null(map) || nrow(map) < 10) return(NULL)
  
  eg <- enrichGO(
    gene = unique(map$ENTREZID),
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(NULL)
  
  as.data.frame(eg) |>
    dplyr::as_tibble() |>
    dplyr::slice_head(n = max_terms)
}

# -----------------------------
# 5) Jaccard heatmap utilities (GO term overlap)
# -----------------------------
make_go_keyed_tbl <- function(go_tables_named) {
  # go_tables_named: named list of tibbles each containing Description
  out <- list()
  for (nm in names(go_tables_named)) {
    tt <- go_tables_named[[nm]]
    if (is.null(tt) || nrow(tt) == 0 || !"Description" %in% colnames(tt)) next
    out[[nm]] <- tt |> dplyr::transmute(key = nm, Description = Description)
  }
  dplyr::bind_rows(out)
}

jaccard_matrix <- function(go_keyed_tbl, min_terms = 5) {
  if (is.null(go_keyed_tbl) || nrow(go_keyed_tbl) == 0) return(NULL)
  
  sets <- go_keyed_tbl |>
    dplyr::group_by(key) |>
    dplyr::summarise(terms = list(unique(Description)), .groups = "drop") |>
    dplyr::filter(lengths(terms) >= min_terms)
  
  keys <- sets$key
  if (length(keys) < 2) return(NULL)
  
  m <- matrix(NA_real_, length(keys), length(keys), dimnames = list(keys, keys))
  for (i in seq_along(keys)) {
    a <- sets$terms[[i]]
    for (j in seq_along(keys)) {
      b <- sets$terms[[j]]
      inter <- length(intersect(a, b))
      uni   <- length(union(a, b))
      m[i, j] <- ifelse(uni == 0, NA_real_, inter / uni)
    }
  }
  m
}

plot_jaccard_heatmap <- function(jacc_mat, title) {
  if (is.null(jacc_mat)) return(NULL)
  ComplexHeatmap::Heatmap(
    jacc_mat,
    name = "Jaccard",
    col = circlize::colorRamp2(c(0, 0.5, 1), c("#FFFFFF", "#FDB863", "#B2182B")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_title = title
  )
}

# -----------------------------
# 6) Main pipeline wrapper
# -----------------------------
run_result2_shared5_pipeline <- function(
    seu,
    shared5 = SHARED5,
    tissue_keep = "N",             # set "N" or c("N","T") as needed
    assay = ASSAY_USE,
    slot = "data",
    p_cut_go = 0.005,              # requested cutoff
    lfc_cut_go = 0.25
) {
  .assert_cols(seu, c(PLATFORM_COL, TISSUE_COL, CELLTYPE_COL, SAMPLE_COL))
  
  # Normalize the celltype labels only (avoid "auto detect")
  seu@meta.data[[CELLTYPE_COL]] <- .trim_levels(as.character(seu@meta.data[[CELLTYPE_COL]]))
  
  # Keep only the requested tissue(s)
  seu_sub <- subset_by_tissue_platform(seu, tissue_keep = tissue_keep, platform_keep = PLATFORM_LEVELS)
  
  # -------------------
  # A) single-cell DEG within each platform: celltype vs rest (shared5 only)
  # -------------------
  deg_sc <- list(
    NUC  = deg_singlecell_celltype_vs_rest(seu_sub, shared5 = shared5, platform = "NUC",  assay = assay, slot = slot),
    CELL = deg_singlecell_celltype_vs_rest(seu_sub, shared5 = shared5, platform = "CELL", assay = assay, slot = slot)
  )
  
  # GO from single-cell markers (separate sn/sc)
  go_sc <- list(NUC = list(), CELL = list())
  for (pf in names(deg_sc)) {
    for (ct in names(deg_sc[[pf]])) {
      go_sc[[pf]][[ct]] <- go_from_deg(
        deg_sc[[pf]][[ct]],
        gene_col = "gene",
        p_col = "p_val_adj",
        lfc_col = "avg_log2FC",
        p_cut = p_cut_go,
        lfc_cut = lfc_cut_go,
        direction = "up",
        ont = "BP",
        max_terms = 30
      )
    }
  }
  
  # -------------------
  # B) pseudo-bulk DEG: NUC vs CELL within each celltype
  # -------------------
  deg_pb <- map(shared5, ~deg_pseudobulk_sn_vs_sc(
    seu = seu_sub,
    celltype = .x,
    tissue_keep = tissue_keep,
    assay = assay,
    slot = slot
  ))
  names(deg_pb) <- shared5
  
  # GO from pseudo-bulk DEG (direction: "up" = NUC higher; "down" = CELL higher)
  go_pb_up <- list()
  go_pb_dn <- list()
  for (ct in shared5) {
    if (is.null(deg_pb[[ct]])) next
    go_pb_up[[ct]] <- go_from_deg(
      deg_pb[[ct]],
      gene_col = "gene",
      p_col = "p_val_adj",
      lfc_col = "avg_diff",
      p_cut = p_cut_go,
      lfc_cut = lfc_cut_go,
      direction = "up",
      ont = "BP",
      max_terms = 30
    )
    go_pb_dn[[ct]] <- go_from_deg(
      deg_pb[[ct]],
      gene_col = "gene",
      p_col = "p_val_adj",
      lfc_col = "avg_diff",
      p_cut = p_cut_go,
      lfc_cut = lfc_cut_go,
      direction = "down",
      ont = "BP",
      max_terms = 30
    )
  }
  
  # -------------------
  # C) Jaccard heatmaps (GO overlap)
  #   - within NUC: 5 celltypes
  #   - within CELL: 5 celltypes
  #   - optional cross-platform: NUC_CT vs CELL_CT (same CT) by combining keys
  # -------------------
  # NUC Jaccard
  go_nuc_named <- setNames(go_sc$NUC, paste0("NUC__", names(go_sc$NUC)))
  go_cell_named <- setNames(go_sc$CELL, paste0("CELL__", names(go_sc$CELL)))
  
  go_tbl_nuc  <- make_go_keyed_tbl(go_nuc_named)
  go_tbl_cell <- make_go_keyed_tbl(go_cell_named)
  
  j_nuc  <- jaccard_matrix(go_tbl_nuc,  min_terms = 5)
  j_cell <- jaccard_matrix(go_tbl_cell, min_terms = 5)
  
  ht_nuc  <- plot_jaccard_heatmap(j_nuc,  title = paste0("GO(BP) Jaccard — NUC (tissue=", paste(tissue_keep, collapse=","), ")"))
  ht_cell <- plot_jaccard_heatmap(j_cell, title = paste0("GO(BP) Jaccard — CELL (tissue=", paste(tissue_keep, collapse=","), ")"))
  
  # Cross-platform Jaccard (NUC__CT & CELL__CT in one matrix)
  go_tbl_both <- dplyr::bind_rows(go_tbl_nuc, go_tbl_cell)
  j_both <- jaccard_matrix(go_tbl_both, min_terms = 5)
  ht_both <- plot_jaccard_heatmap(j_both, title = paste0("GO(BP) Jaccard — NUC vs CELL (tissue=", paste(tissue_keep, collapse=","), ")"))
  
  list(
    inputs = list(
      shared5 = shared5,
      tissue_keep = tissue_keep,
      columns = list(
        platform = PLATFORM_COL,
        tissue = TISSUE_COL,
        celltype = CELLTYPE_COL,
        sample = SAMPLE_COL
      )
    ),
    deg = list(
      single_cell = deg_sc,          # deg$single_cell$NUC[[ct]], deg$single_cell$CELL[[ct]]
      pseudo_bulk = deg_pb           # deg$pseudo_bulk[[ct]]
    ),
    go = list(
      single_cell = go_sc,           # go$single_cell$NUC[[ct]], go$single_cell$CELL[[ct]]
      pseudo_bulk = list(
        up_NUC = go_pb_up,           # go$pseudo_bulk$up_NUC[[ct]]
        down_CELL = go_pb_dn         # go$pseudo_bulk$down_CELL[[ct]]
      )
    ),
    jaccard = list(
      nuc = j_nuc,
      cell = j_cell,
      nuc_cell = j_both
    ),
    plots = list(
      jaccard_go_nuc = ht_nuc,
      jaccard_go_cell = ht_cell,
      jaccard_go_nuc_cell = ht_both
    )
  )
}

################################################################################
## END
################################################################################
