###############################################################################
# Figure 4a-j: AT2-to-AT1 trajectories (Monocle2) and TF regulon activity (SCENIC)
# Seurat (sn/sc) -> Monocle2 (AT1/AT2) -> SCENIC summary
# Outputs: UMAP, DotPlot, Monocle pseudotime/cluster/celltype, TF-marker cor dotplot, AUC heatmap
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  
  library(monocle)          # Monocle2
  library(Matrix)
  
  library(AUCell)           # SCENIC output reading
  library(SCENIC)           # (optional full)
  library(ComplexHeatmap)
  library(circlize)
  
  library(reshape2)
  library(forcats)
})

set.seed(123)

###############################################################################
# 0) USER CONFIG (ONLY PLACE YOU EDIT)
###############################################################################

RUN_SCENIC_FULL <- FALSE  # default: use existing SCENIC results and make plots

CFG <- list(
  sn = list(
    label = "LUNG_EP_NUC",
    seurat_rds = "PATH/TO/sn_seurat.rds",   # <- set
    meta = list(
      tissue_col   = "TissueDiseaseState",  # N/T
      normal_value = "N",
      celltype_col = "celltype2",
      at1_label    = "AT1 ",
      at2_label    = "AT2",
      batch_col    = "orig.ident",
      cluster_col  = "seurat_clusters"),
    seurat = list(
      nfeatures = 2000,
      npcs      = 50,
      pc_use    = 13,
      res       = 0.3,
      pt_size   = 0.2),
    markers = list(
      AT1 = c("AGER","CLIC5"),
      AT2 = c("SFTPB","SFTPC","SFTPD","ETV5")),
    monocle = list(
      ordering_mean_expr = 0.001,
      reverse = TRUE),
    scenic = list(
      # If already computed (recommended for GitHub script)
      regulonAUC_rds_at1 = "PATH/TO/int/3.4_regulonAUC_at1.Rds",
      regulonAUC_rds_at2 = "PATH/TO/int/3.4_regulonAUC_at2.Rds",
      regulonTargets_csv = "PATH/TO/output/Step2_regulonTargetsInfo.csv",
      metadata_csv       = "PATH/TO/AT1_AT2_metadata.csv",
      
      # (optional) for full run
      dbDir = "PATH/TO/cisTarget_databases",
      dbs = c(
        "500bp"="hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
        "10kb" ="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"),
      nCores = 8)),
  
  sc = list(
    label = "LUNG_EP_CELL",
    seurat_rds = "PATH/TO/sc_seurat.rds",   # <- set
    meta = list(
      tissue_col   = "TissueDiseaseState",
      normal_value = "N",
      celltype_col = "celltype2",
      at1_label    = "AT1 ",
      at2_label    = "AT2",
      batch_col    = "orig.ident",
      cluster_col  = "seurat_clusters"),
    seurat = list(
      nfeatures = 2000,
      npcs      = 50,
      pc_use    = 13,
      res       = 0.3,
      pt_size   = 0.2 ),
    markers = list(
      AT1 = c("AGER","CLIC5"),
      AT2 = c("SFTPB","SFTPC","SFTPD","ETV5")
    ),
    monocle = list(
      ordering_mean_expr = 0.001,
      reverse = TRUE),
    scenic = list(
      regulonAUC_rds_at1 = "PATH/TO/int/3.4_regulonAUC_at1.Rds",
      regulonAUC_rds_at2 = "PATH/TO/int/3.4_regulonAUC_at2.Rds",
      regulonTargets_csv = "PATH/TO/output/Step2_regulonTargetsInfo.csv",
      metadata_csv       = "PATH/TO/AT1_AT2_metadata.csv",
      dbDir = "PATH/TO/cisTarget_databases",
      dbs = c(
        "500bp"="hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
        "10kb" ="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"),
      nCores = 8
    )
  )
)

OUT_DIR <- "Figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# 1) SMALL HELPERS
###############################################################################

.save_plot <- function(p, filename, w=4, h=3, dpi=600){
  ggsave(file.path(OUT_DIR, filename), p, width=w, height=h, dpi=dpi)
}
.save_pdf <- function(p, filename, w=4, h=3){
  ggsave(file.path(OUT_DIR, filename), p, width=w, height=h, units="in")
}

.clean_regulon_names <- function(x){
  x <- gsub(" \\(\\d+g\\)", "", x)
  x <- gsub("_extended", "", x)
  toupper(trimws(x))
}

###############################################################################
# 2) SEURAT: preprocess -> harmony(optional) -> umap/cluster -> plots
###############################################################################

run_seurat_core <- function(seu, cfg){
  m <- cfg$meta
  s <- cfg$seurat
  
  # subset normal
  stopifnot(m$tissue_col %in% colnames(seu@meta.data))
  seu <- subset(seu, seu@meta.data[[m$tissue_col]] == m$normal_value)
  
  DefaultAssay(seu) <- "RNA"
  
  # variable genes / cell cycle regress
  seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=s$nfeatures)
  if (exists("cc.genes")) {
    seu <- CellCycleScoring(seu, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=FALSE)
    seu <- ScaleData(seu, vars.to.regress=c("S.Score","G2M.Score"))
  } else {
    seu <- ScaleData(seu)
  }
  
  seu <- RunPCA(seu, features=VariableFeatures(seu), npcs=s$npcs)
  
  # harmony if available
  if (requireNamespace("harmony", quietly = TRUE)) {
    if (!m$batch_col %in% colnames(seu@meta.data)) stop("batch_col not found: ", m$batch_col)
    seu <- harmony::RunHarmony(seu, group.by.vars = m$batch_col)
    red <- "harmony"
  } else {
    red <- "pca"
  }
  
  seu <- RunUMAP(seu, reduction=red, dims=1:s$pc_use)
  seu <- FindNeighbors(seu, reduction=red, dims=1:s$pc_use)
  seu <- FindClusters(seu, resolution=s$res)
  
  seu
}

plot_seurat_umap_cluster <- function(seu, cfg, platform){
  s <- cfg$seurat
  p <- DimPlot(seu, reduction="umap", group.by=cfg$meta$cluster_col, label=TRUE, pt.size=s$pt_size) +
    theme_void() + theme(legend.position="none")
  .save_plot(p, sprintf("%s_%s_UMAP_cluster.png", cfg$label, platform), w=3.2, h=3.2)
  .save_pdf(p,  sprintf("%s_%s_UMAP_cluster.pdf", cfg$label, platform), w=3.2, h=3.2)
}

plot_seurat_dotplot_markers <- function(seu, cfg, platform){
  markers <- unlist(cfg$markers, use.names=FALSE)
  Idents(seu) <- cfg$meta$cluster_col
  
  p <- DotPlot(seu, features=markers, dot.min=0, dot.scale=5) +
    scale_x_discrete(position="top") +
    coord_flip() +
    theme_bw() +
    theme(
      axis.title=element_blank(),
      panel.grid=element_blank())
  
  .save_plot(p, sprintf("%s_%s_DotPlot_AT1AT2.png", cfg$label, platform), w=4.2, h=2.6)
  .save_pdf(p,  sprintf("%s_%s_DotPlot_AT1AT2.pdf", cfg$label, platform), w=4.2, h=2.6)
}

###############################################################################
# 3) MONOCLE2: make cds -> subset AT1/AT2 -> DDRTree -> orderCells -> plots
###############################################################################

make_monocle_cds_from_seurat <- function(seu, cfg){
  m <- cfg$meta
  # counts
  exprs <- GetAssayData(seu, assay="RNA", slot="counts")
  exprs <- as(exprs, "dgCMatrix")
  
  pd <- seu@meta.data
  pd$Cell_Name <- rownames(pd)
  
  # align cell order
  pd <- pd[colnames(exprs), , drop=FALSE]
  
  gene_ann <- data.frame(gene_short_name=rownames(exprs), row.names=rownames(exprs))
  pDataNew <- new("AnnotatedDataFrame", data=pd)
  fDataNew <- new("AnnotatedDataFrame", data=gene_ann)
  
  cds <- monocle::newCellDataSet(
    cellData = exprs,
    phenoData = pDataNew,
    featureData = fDataNew,
    expressionFamily = negbinomial.size(),
    lowerDetectionLimit = 0)
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  # ordering genes by dispersion
  disp <- dispersionTable(cds)
  ordering_genes <- subset(disp, dispersion_empirical > dispersion_fit & mean_expression > cfg$monocle$ordering_mean_expr)
  cds <- setOrderingFilter(cds, ordering_genes$gene_id)
  
  list(cds=cds, ordering_genes=ordering_genes)
}

run_monocle_at12 <- function(cds, cfg){
  m <- cfg$meta
  stopifnot(m$celltype_col %in% colnames(pData(cds)))
  
  # subset AT1/AT2 within normal epithelial object
  ct <- trimws(as.character(pData(cds)[[m$celltype_col]]))
  keep <- ct %in% c(trimws(m$at1_label), trimws(m$at2_label))
  cds_at <- cds[, keep]
  
  cds_at <- estimateSizeFactors(cds_at)
  cds_at <- estimateDispersions(cds_at)
  
  disp <- dispersionTable(cds_at)
  ordering_genes <- subset(disp, dispersion_empirical > dispersion_fit & mean_expression > cfg$monocle$ordering_mean_expr)
  cds_at <- setOrderingFilter(cds_at, ordering_genes$gene_id)
  
  cds_at <- reduceDimension(
    cds_at,
    max_components = 4,
    reduction_method = "DDRTree",
    verbose = TRUE,
    norm_method = "vstExprs",
    relative_expr = TRUE,
    scaling = TRUE)
  
  cds_at <- orderCells(cds_at, reverse = cfg$monocle$reverse)
  
  cds_at
}

plot_monocle_trajectories <- function(cds_at, cfg, platform){
  # Expect in phenoData: Pseudotime / State / seurat_clusters or cluster name
  pdat <- pData(cds_at)
  
  # robust cluster_name
  if (!"cluster_name" %in% colnames(pdat)) {
    if (cfg$meta$cluster_col %in% colnames(pdat)) {
      pdat$cluster_name <- paste0(platform, "_C", sprintf("%02d", as.integer(as.character(pdat[[cfg$meta$cluster_col]]))))
      pData(cds_at) <- pdat
    }
  }
  
  # pseudotime
  pdf(file.path(OUT_DIR, sprintf("%s_%s_Monocle_Pseudotime.pdf", cfg$label, platform)), width=4.4, height=2.4)
  print(
    plot_cell_trajectory(cds_at, color_by="Pseudotime", cell_size=0.2) +
      theme_void() + theme(legend.position="right") )
  dev.off()
  
  # cluster_name
  pdf(file.path(OUT_DIR, sprintf("%s_%s_Monocle_Cluster.pdf", cfg$label, platform)), width=4.4, height=2.8)
  print(
    plot_cell_trajectory(cds_at, color_by="cluster_name", cell_size=0.4) +
      theme_void() + theme(legend.position="bottom"))
  dev.off()
  
  # celltype
  pdf(file.path(OUT_DIR, sprintf("%s_%s_Monocle_Celltype.pdf", cfg$label, platform)), width=4.4, height=2.8)
  print(
    plot_cell_trajectory(cds_at, color_by=cfg$meta$celltype_col, cell_size=0.4) +
      theme_void() + theme(legend.position="bottom"))
  dev.off()
}

###############################################################################
# 4) SCENIC SUMMARY: (existing AUC) -> AT1/AT2 mean AUC -> heatmap
#                  + TF-marker correlation dotplot (expression-based)
###############################################################################

read_auc_mean_by_celltype <- function(cfg){
  sc <- cfg$scenic
  meta <- read.csv(sc$metadata_csv, row.names=1)
  meta[[cfg$meta$celltype_col]] <- trimws(as.character(meta[[cfg$meta$celltype_col]]))
  
  auc_at1 <- getAUC(readRDS(sc$regulonAUC_rds_at1))
  auc_at2 <- getAUC(readRDS(sc$regulonAUC_rds_at2))
  
  # align to meta rows
  cells_at1 <- intersect(rownames(meta)[meta[[cfg$meta$celltype_col]] == trimws(cfg$meta$at1_label)], colnames(auc_at1))
  cells_at2 <- intersect(rownames(meta)[meta[[cfg$meta$celltype_col]] == trimws(cfg$meta$at2_label)], colnames(auc_at2))
  
  at1_mean <- rowMeans(as.matrix(auc_at1[, cells_at1, drop=FALSE]), na.rm=TRUE)
  at2_mean <- rowMeans(as.matrix(auc_at2[, cells_at2, drop=FALSE]), na.rm=TRUE)
  
  # clean regulon names & aggregate by TF
  tf_names1 <- .clean_regulon_names(names(at1_mean))
  tf_names2 <- .clean_regulon_names(names(at2_mean))
  
  df <- data.frame(
    TF = tf_names1,
    AT1 = as.numeric(at1_mean),
    stringsAsFactors = FALSE)
  df2 <- data.frame(
    TF = tf_names2,
    AT2 = as.numeric(at2_mean),
    stringsAsFactors = FALSE)
  
  # merge by TF (some may be missing)
  mdf <- full_join(df, df2, by="TF")
  mdf$AT1[is.na(mdf$AT1)] <- 0
  mdf$AT2[is.na(mdf$AT2)] <- 0
  
  # aggregate duplicates (mean)
  mdf <- mdf %>%
    group_by(TF) %>%
    summarise(AT1=mean(AT1), AT2=mean(AT2), .groups="drop")
  
  mdf
}

plot_auc_heatmap <- function(mean_auc_df, cfg, platform, top_n=30){
  mean_auc_df <- mean_auc_df %>%
    mutate(delta = AT2 - AT1) %>%
    arrange(desc(abs(delta))) %>%
    slice_head(n=top_n)
  
  mat <- as.matrix(mean_auc_df[, c("AT1","AT2")])
  rownames(mat) <- mean_auc_df$TF
  
  col_fun <- colorRamp2(c(min(mat), median(mat), max(mat)), c("skyblue","white","firebrick"))
  
  png(file.path(OUT_DIR, sprintf("%s_%s_SCENIC_AUC_heatmap_top%d.png", cfg$label, platform, top_n)),
      width=1800, height=2800, res=600)
  Heatmap(mat, name="Regulon AUC",
          cluster_rows=TRUE, cluster_columns=FALSE,
          col=col_fun, show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=11, fontface="bold"))
  dev.off()
  
  pdf(file.path(OUT_DIR, sprintf("%s_%s_SCENIC_AUC_heatmap_top%d.pdf", cfg$label, platform, top_n)),
      width=4.0, height=6.5)
  Heatmap(mat, name="Regulon AUC",
          cluster_rows=TRUE, cluster_columns=FALSE,
          col=col_fun, show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=11, fontface="bold"))
  dev.off()
}

tf_marker_correlation_dotplot <- function(seu, tf_vec, markers_list, platform, out_prefix){
  expr <- GetAssayData(seu, assay="RNA", slot="data")  # log-normalized
  
  get_cor <- function(tfs, markers){
    tfs <- intersect(tfs, rownames(expr))
    markers <- intersect(markers, rownames(expr))
    if (length(tfs) == 0 || length(markers) == 0) return(NULL)
    cor(t(as.matrix(expr[tfs, , drop=FALSE])),
        t(as.matrix(expr[markers, , drop=FALSE])),
        method="pearson", use="pairwise.complete.obs")
  }
  
  cor_at1 <- get_cor(tf_vec, markers_list$AT1)
  cor_at2 <- get_cor(tf_vec, markers_list$AT2)
  
  to_df <- function(mat, celltype){
    if (is.null(mat)) return(NULL)
    df <- reshape2::melt(mat)
    colnames(df) <- c("TF","Marker","Pearson")
    df$CellType <- celltype
    df
  }
  
  df <- bind_rows(to_df(cor_at1,"AT1"), to_df(cor_at2,"AT2"))
  if (nrow(df) == 0) return(invisible(NULL))
  
  # rank TF order by max |r| in each facet
  tf_order <- df %>%
    group_by(CellType, TF) %>%
    summarise(score=max(abs(Pearson), na.rm=TRUE), .groups="drop") %>%
    arrange(CellType, desc(score), TF) %>%
    pull(TF) %>% unique()
  
  df$TF <- factor(df$TF, levels=tf_order)
  
  p <- ggplot(df, aes(x=Marker, y=TF, color=Pearson, size=abs(Pearson))) +
    geom_point(alpha=0.9) +
    facet_wrap(~CellType, scales="free_y") +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
    theme_linedraw() +
    theme(
      axis.text.x=element_text(angle=45, hjust=1),
      strip.text=element_text(face="bold", size=12)
    ) +
    labs(x="Marker", y="TF", color="Pearson r", size="|r|",
         title=sprintf("%s TF–Marker correlation (%s)", out_prefix, platform))
  
  ggsave(file.path(OUT_DIR, sprintf("%s_%s_TF_marker_correlation.pdf", out_prefix, platform)),
         p, width=5.2, height=3.8, units="in")
  ggsave(file.path(OUT_DIR, sprintf("%s_%s_TF_marker_correlation.png", out_prefix, platform)),
         p, width=5.2, height=3.8, dpi=600)
  
  p
}

###############################################################################
# 5) OPTIONAL: SCENIC FULL PIPELINE (off by default)
###############################################################################
run_scenic_full_optional <- function(exprMat, cellInfo, cfg, datasetTitle){
  sc <- cfg$scenic
  stopifnot(all(file.exists(file.path(sc$dbDir, sc$dbs))))
  
  scenicOptions <- initializeScenic(
    org="hgnc",
    dbDir=sc$dbDir,
    dbs=sc$dbs,
    datasetTitle=datasetTitle,
    nCores=sc$nCores)
  scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
  
  exprMat_f <- SCENIC::filterGeneExprMatrix(exprMat, scenicOptions)
  runCorrelation(exprMat_f, scenicOptions)
  runGenie3(exprMat_f, scenicOptions)
  runSCENIC_1_coexNetwork2modules(scenicOptions)
  runSCENIC_2_createRegulons(scenicOptions)
  runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_f)
  runSCENIC_4_aucell_binarize(scenicOptions, exprMat_f)
  
  scenicOptions
}

###############################################################################
# 6) MAIN: run for sn + sc
###############################################################################

run_one_platform <- function(platform, cfg){
  message("=== Running: ", platform, " / ", cfg$label, " ===")
  
  # load seurat
  seu <- readRDS(cfg$seurat_rds)
  
  # seurat core
  seu <- run_seurat_core(seu, cfg)
  
  # plots: umap + dotplot
  plot_seurat_umap_cluster(seu, cfg, platform)
  plot_seurat_dotplot_markers(seu, cfg, platform)
  
  # monocle
  monocle_obj <- make_monocle_cds_from_seurat(seu, cfg)
  cds_at <- run_monocle_at12(monocle_obj$cds, cfg)
  plot_monocle_trajectories(cds_at, cfg, platform)
  
  # SCENIC summary (existing AUC)
  mean_auc_df <- read_auc_mean_by_celltype(cfg)
  plot_auc_heatmap(mean_auc_df, cfg, platform, top_n=30)
  
  # TF list from regulonTargetsInfo (highConf) if available
  tf_vec <- NULL
  if (file.exists(cfg$scenic$regulonTargets_csv)) {
    rti <- read.csv(cfg$scenic$regulonTargets_csv)
    if ("highConfAnnot" %in% colnames(rti) && "TF" %in% colnames(rti)) {
      tf_vec <- unique(toupper(trimws(rti$TF[rti$highConfAnnot %in% c(TRUE, "TRUE")])))
    } else if ("TF" %in% colnames(rti)) {
      tf_vec <- unique(toupper(trimws(rti$TF)))
    }
  }
  if (is.null(tf_vec) || length(tf_vec) == 0) {
    # fallback: use top differential TFs by |delta|
    tf_vec <- mean_auc_df %>%
      mutate(delta=AT2-AT1) %>%
      arrange(desc(abs(delta))) %>%
      slice_head(n=50) %>%
      pull(TF)
  }
  
  tf_marker_correlation_dotplot(
    seu=seu,
    tf_vec=tf_vec,
    markers_list=cfg$markers,
    platform=platform,
    out_prefix=cfg$label)
  
  invisible(list(seurat=seu, cds_at=cds_at, mean_auc=mean_auc_df))
}

res_sn <- run_one_platform("sn", CFG$sn)
res_sc <- run_one_platform("sc", CFG$sc)

message("DONE. Figures saved under: ", OUT_DIR)
