############################################################
## Figure 5a-h and Supplementary Figure 6a-d
## Core6-only: Module definition -> enrichment -> corr heatmaps
##            -> attach M1–M4 scores + dominant module to Seurat
##            -> plots (violin/topgenes/fraction/umap)
##
## - Single final script (copy & paste ready)
## - Requires: cds_traj (monocle2 CellDataSet with Pseudotime),
##             seu (Seurat object),
##             de_all (named list of DEG tables per sample),
##             core_genes_6 (gene IDs present in exprs(cds_traj) rownames)
##             core6_genes (SYMBOL vector) for DEG-corr section
##
## NOTE:
##  - core5 is NOT used anywhere.
##  - aux5 / core5_genes not needed.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(Matrix)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(msigdbr)
  library(monocle)      # monocle2
  library(Seurat)
})

## =========================================================
## [A] CONFIG (edit here)
## =========================================================
## ---- REQUIRED (must exist BEFORE running): OUTDIR, CACHE_SALT, cds_traj, seu, de_all, core_genes_6, core6_genes
stopifnot(exists("OUTDIR"), exists("CACHE_SALT"))
stopifnot(exists("cds_traj"), inherits(cds_traj, "CellDataSet"))
stopifnot("Pseudotime" %in% colnames(pData(cds_traj)))
stopifnot(exists("seu"), inherits(seu, "Seurat"))
stopifnot(exists("de_all"), is.list(de_all), length(de_all) > 0)
stopifnot(exists("core_genes_6"), is.character(core_genes_6), length(core_genes_6) > 0)
stopifnot(exists("core6_genes"), is.character(core6_genes), length(core6_genes) > 0)

outdir <- OUTDIR
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## ---- run switches
RUN_MODULE_DEFINE   <- TRUE
RUN_ENRICHMENT      <- TRUE
RUN_DEG_CORR        <- TRUE
RUN_SEURAT_SCORE    <- TRUE
RUN_PLOTS_VIOLIN    <- TRUE
RUN_PLOTS_TOPGENES  <- TRUE
RUN_FRACTION_PLOTS  <- TRUE
RUN_UMAP_PLOTS      <- TRUE

## ---- module/binning/kmeans
N_BIN         <- 50
MIN_BIN_CELLS <- 20
NUM_MODULES   <- 4
KMEANS_SEED   <- 1
KMEANS_NSTART <- 50
KMEANS_ITMAX  <- 200

modules <- paste0("M", seq_len(NUM_MODULES))
module_colors <- c("M1"="#4E79A7","M2"="#F28E2B","M3"="#59A14F","M4"="#E15759")[modules]

## ---- enrichment settings
source_levels <- c("C6","MSigDB","KEGG","GOBP")
source_cols   <- c(GOBP="#7B2CBF", KEGG="#3A86FF", MSigDB="#2EC4B6", C6="#FF595E")

## ---- DEG correlation settings
CTTDS_COL <- "celltype_TissueDiseaseState"
GENE_COL  <- "SYMBOL"
LFC_COL   <- "log2FC_Malig_vs_AT2"
Q_COL     <- "p.adjust"

Q_CUT_WHOLE             <- 0.05
MIN_SUPPORT_SIDE_WHOLE  <- 1
MIN_SUPPORT_SIDE_CORE   <- 1

SCALE_MIN <- 0.2
SCALE_MAX <- 1.0
N_BREAKS  <- 201

## ---- plots settings
TOPGENE_N <- 20
VIOLIN_USE_GGH4X <- TRUE  # if ggh4x installed, use colored strips

message("=== [CONFIG] ===")
message("OUTDIR=", outdir, " | CACHE_SALT=", CACHE_SALT)
message("N_BIN=", N_BIN, " MIN_BIN_CELLS=", MIN_BIN_CELLS, " NUM_MODULES=", NUM_MODULES)

## =========================================================
## [B] Small helpers (IO)
## =========================================================
pdf_safe <- function(pdf_path, expr, width = 6, height = 4.5) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

save_rds_safe <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, path)
  message("[SAVE] ", path)
}

as_df_safe <- function(x){
  if (is.null(x)) return(NULL)
  df <- tryCatch(as.data.frame(x), error=function(e) NULL)
  if (is.null(df) || nrow(df)==0) return(NULL)
  df
}

ratio_to_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- as.character(x)
  sapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_real_)
    sp <- stringr::str_split(s, "/", simplify = TRUE)
    if (ncol(sp) < 2) return(NA_real_)
    suppressWarnings(as.numeric(sp[1]) / as.numeric(sp[2]))
  })
}

sym2entrez <- function(sym){
  sym <- unique(sym[!is.na(sym) & sym != ""])
  if (!length(sym)) return(character())
  AnnotationDbi::select(org.Hs.eg.db, keys = sym, keytype="SYMBOL", columns="ENTREZID") %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID) %>%
    pull(ENTREZID) %>%
    as.character() %>%
    unique()
}

## =========================================================
## [C] Module definition helpers
## =========================================================
.make_gene_bin_profile <- function(cds, genes, ord, bin_id, N_BIN, min_bin_cells=20) {
  genes <- intersect(genes, rownames(exprs(cds)))
  if (length(genes) == 0) stop("[.make_gene_bin_profile] No genes found in exprs(cds).")
  
  E_sp <- exprs(cds)[genes, ord, drop = FALSE]
  if (!inherits(E_sp, "dgCMatrix")) E_sp <- as(E_sp, "dgCMatrix")
  
  bin_n <- tabulate(bin_id, nbins = N_BIN)
  
  E_bin <- sapply(seq_len(N_BIN), function(b) {
    idx <- which(bin_id == b)
    if (length(idx) < min_bin_cells) return(rep(NA_real_, length(genes)))
    Matrix::rowMeans(E_sp[, idx, drop = FALSE])
  })
  if (is.vector(E_bin)) E_bin <- matrix(E_bin, nrow = length(genes), ncol = N_BIN)
  rownames(E_bin) <- genes
  colnames(E_bin) <- paste0("bin", seq_len(N_BIN))
  
  E_bin_z <- t(apply(E_bin, 1, function(x) {
    ok <- is.finite(x)
    if (sum(ok) < 3) return(rep(NA_real_, length(x)))
    mu <- mean(x[ok]); sdv <- stats::sd(x[ok])
    if (!is.finite(sdv) || sdv == 0) return(rep(NA_real_, length(x)))
    as.numeric((x - mu) / sdv)
  }))
  rownames(E_bin_z) <- rownames(E_bin)
  colnames(E_bin_z) <- colnames(E_bin)
  
  list(E_bin = E_bin, E_bin_z = E_bin_z, bin_n = bin_n)
}

define_modules_from_core <- function(
    cds_all, core_genes, tag,
    outdir, CACHE_SALT,
    N_BIN = 50,
    min_bin_cells = 20,
    K_MODULES = 4,
    seed = 1,
    nstart = 50,
    itmax = 200,
    force = TRUE,
    BIN_COR_CUTREE = 4,
    BIN_COR_BREAKS = c("quantile", "symmetric"),
    BIN_COR_NCOL   = 100,
    BIN_COR_QCLIP  = 0.98,
    BIN_COR_BORDER = NA,
    module_colors = NULL
) {
  stopifnot(is.character(tag), length(tag) == 1)
  stopifnot(is.character(outdir), length(outdir) == 1)
  stopifnot(is.character(CACHE_SALT), length(CACHE_SALT) == 1)
  stopifnot(is.character(core_genes), length(core_genes) > 0)
  
  if (!inherits(cds_all, "CellDataSet")) stop("cds_all must be a monocle2 CellDataSet.")
  if (!("Pseudotime" %in% colnames(pData(cds_all)))) stop("Pseudotime not found. Run orderCells() first.")
  BIN_COR_BREAKS <- match.arg(BIN_COR_BREAKS)
  
  subdir <- file.path(outdir, tag)
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  module_tbl_rds <- file.path(subdir, sprintf("module_tbl_fixed_%s_%s.rds", tag, CACHE_SALT))
  module_aux_rds <- file.path(subdir, sprintf("module_aux_fixed_%s_%s.rds", tag, CACHE_SALT))
  support_pdf    <- file.path(subdir, sprintf("module_support_plots_1_%s_%s.pdf", tag, CACHE_SALT))
  
  if (!isTRUE(force) && file.exists(module_tbl_rds) && file.exists(module_aux_rds)) {
    message("[LOAD modules] ", tag)
    module_tbl <- readRDS(module_tbl_rds)
    aux <- readRDS(module_aux_rds)
    aux$module_tbl <- module_tbl
    return(aux)
  }
  
  message("=== [DEFINE modules] tag=", tag, " | core genes=", length(core_genes))
  
  modules_local <- paste0("M", seq_len(K_MODULES))
  if (is.null(module_colors)) {
    base4 <- c("M1"="#4E79A7","M2"="#F28E2B","M3"="#59A14F","M4"="#E15759")
    if (K_MODULES <= 4) {
      module_colors <- base4[modules_local]
    } else {
      cols <- grDevices::colorRampPalette(unname(base4))(K_MODULES)
      module_colors <- setNames(cols, modules_local)
    }
  } else {
    if (is.null(names(module_colors))) stop("module_colors must be named (M1..).")
    miss <- setdiff(modules_local, names(module_colors))
    if (length(miss) > 0) stop("module_colors missing: ", paste(miss, collapse = ", "))
    module_colors <- module_colors[modules_local]
  }
  
  ## pseudotime order + binning
  pt_raw <- pData(cds_all)$Pseudotime
  pt <- suppressWarnings(as.numeric(as.character(pt_raw)))
  ok <- which(is.finite(pt))
  if (length(ok) < 500) stop("[define_modules_from_core] Too few finite Pseudotime cells: ", length(ok))
  
  ord <- ok[order(pt[ok])]
  pt_ord <- pt[ord]
  
  rng <- range(pt_ord, finite = TRUE)
  if (!all(is.finite(rng)) || diff(rng) <= 0) stop("Invalid Pseudotime range.")
  pt_norm <- (pt_ord - rng[1]) / diff(rng)
  
  bin_id <- as.integer(cut(pt_norm, breaks = N_BIN, labels = FALSE, include.lowest = TRUE))
  bin_id[is.na(bin_id)] <- N_BIN
  
  ## profiles
  prof <- .make_gene_bin_profile(
    cds = cds_all,
    genes = core_genes,
    ord = ord,
    bin_id = bin_id,
    N_BIN = N_BIN,
    min_bin_cells = min_bin_cells
  )
  
  E_bin_z <- as.matrix(prof$E_bin_z)
  bin_n   <- as.integer(prof$bin_n)
  
  if (ncol(E_bin_z) != N_BIN && nrow(E_bin_z) == N_BIN) E_bin_z <- t(E_bin_z)
  stopifnot(ncol(E_bin_z) == N_BIN)
  
  keep_gene <- apply(E_bin_z, 1, function(x) any(is.finite(x)))
  E_bin_z <- E_bin_z[keep_gene, , drop = FALSE]
  if (nrow(E_bin_z) < K_MODULES * 5) stop("Too few usable genes after NA-filtering: ", nrow(E_bin_z))
  
  ## kmeans on good bins (finite only)
  good_bins <- which(bin_n >= min_bin_cells)
  if (length(good_bins) < 3) {
    stop("[define_modules_from_core] Too few bins with enough cells: ",
         length(good_bins), " (need >= 3). Lower MIN_BIN_CELLS or N_BIN.")
  }
  
  E_km <- E_bin_z[, good_bins, drop = FALSE]
  keep_km_gene <- apply(E_km, 1, function(x) all(is.finite(x)))
  E_km <- E_km[keep_km_gene, , drop = FALSE]
  if (nrow(E_km) < K_MODULES * 5) {
    stop("[define_modules_from_core] Too few usable genes for kmeans after filtering: ", nrow(E_km))
  }
  
  set.seed(seed)
  km <- stats::kmeans(E_km, centers = K_MODULES, nstart = nstart, iter.max = itmax)
  
  gene_cluster_raw <- km$cluster
  names(gene_cluster_raw) <- rownames(E_km)
  
  centers <- km$centers
  cluster_ids <- sort(unique(gene_cluster_raw))
  cluster_peak <- apply(centers, 1, which.max)
  cluster_order <- order(cluster_peak)
  
  cluster2mod <- setNames(paste0("M", seq_along(cluster_order)),
                          as.character(cluster_ids[cluster_order]))
  
  gene_module <- cluster2mod[as.character(gene_cluster_raw)]
  names(gene_module) <- names(gene_cluster_raw)
  gene_module <- factor(gene_module, levels = modules_local)
  
  ## module table
  sym <- NA_character_
  if ("gene_short_name" %in% colnames(fData(cds_all))) {
    sym <- as.character(fData(cds_all)[names(gene_module), "gene_short_name"])
  }
  module_tbl <- tibble(
    gene_id = names(gene_module),
    SYMBOL  = sym,
    Module  = gene_module
  )
  saveRDS(module_tbl, module_tbl_rds)
  
  E_bin_z_use <- E_bin_z[names(gene_cluster_raw), , drop = FALSE]
  
  ## module mean curves + peaks
  df_mean <- map_dfr(modules_local, function(m) {
    gs <- names(gene_module)[gene_module == m]
    mat <- E_bin_z_use[gs, , drop = FALSE]
    tibble(Module = m, bin = seq_len(N_BIN), z = colMeans(mat, na.rm = TRUE))
  }) %>% mutate(Module = factor(Module, levels = modules_local))
  
  df_peak <- df_mean %>%
    group_by(Module) %>%
    summarise(peak_bin = bin[which.max(z)], .groups = "drop")
  
  ## bin-bin correlation
  bin_cor <- stats::cor(E_bin_z, method = "pearson", use = "pairwise.complete.obs")
  bin_cor[!is.finite(bin_cor)] <- 0
  
  ## gene peak bin distribution
  gene_peak_bin <- apply(E_bin_z, 1, function(x) if (all(!is.finite(x))) NA_integer_ else which.max(x))
  df_gene_peak <- tibble(
    gene_id  = names(gene_peak_bin),
    Module   = gene_module[names(gene_peak_bin)],
    peak_bin = as.integer(gene_peak_bin)
  ) %>% filter(!is.na(peak_bin)) %>% mutate(Module = factor(Module, levels = modules_local))
  
  ## plots
  p_curve <- ggplot(df_mean, aes(x = bin, y = z, color = Module)) +
    geom_line(linewidth = 1.1) +
    geom_vline(data = df_peak, aes(xintercept = peak_bin, color = Module),
               linetype = "dashed", linewidth = 0.7, alpha = 0.8, show.legend = FALSE) +
    scale_color_manual(values = module_colors, drop = FALSE) +
    theme_bw(base_size = 17.5) +
    theme(panel.grid = element_blank()) +
    labs(title = sprintf("[%s] Module mean (z)", tag),
         x = "Pseudotime bin", y = "Mean z-score")
  
  mat_cor <- bin_cor
  if (BIN_COR_BREAKS == "quantile") {
    vv <- as.numeric(mat_cor); vv <- vv[is.finite(vv)]
    if (length(vv) == 0) vv <- c(-1, 1)
    qs <- stats::quantile(vv, probs = seq(0, 1, length.out = BIN_COR_NCOL), na.rm = TRUE)
    breaks <- unique(as.numeric(qs))
    if (length(breaks) < 5) breaks <- seq(min(vv), max(vv), length.out = BIN_COR_NCOL)
  } else {
    vmax <- stats::quantile(abs(as.numeric(mat_cor)), probs = BIN_COR_QCLIP, na.rm = TRUE)
    if (!is.finite(vmax) || vmax <= 0) vmax <- 1
    breaks <- seq(-vmax, vmax, length.out = BIN_COR_NCOL)
  }
  cor_cols <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(BIN_COR_NCOL - 1)
  
  ph_cor <- pheatmap::pheatmap(
    mat_cor,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    cutree_rows = BIN_COR_CUTREE,
    cutree_cols = BIN_COR_CUTREE,
    breaks = breaks,
    color = cor_cols,
    show_rownames = FALSE,
    show_colnames = FALSE,
    border_color = BIN_COR_BORDER,
    legend = TRUE,
    silent = TRUE
  )
  
  p_peak_hist <- ggplot(df_gene_peak, aes(x = peak_bin, fill = Module)) +
    geom_histogram(binwidth = 1, boundary = 0.5, color = NA, alpha = 0.9) +
    scale_fill_manual(values = module_colors, drop = FALSE) +
    theme_bw(base_size = 17.5) +
    theme(panel.grid = element_blank()) +
    labs(title = sprintf("[%s] Gene peak-bin distribution", tag),
         x = "Pseudotime bin", y = "nGenes")
  
  pdf_safe(support_pdf, {
    print(p_curve)
    grid::grid.newpage(); grid::grid.draw(ph_cor$gtable)
    print(p_peak_hist)
  }, width = 5.4, height = 4.4)
  message("[DONE] module support PDF: ", support_pdf)
  
  aux <- list(
    tag = tag,
    subdir = subdir,
    module_tbl = module_tbl,
    module_colors = module_colors,
    modules = modules_local,
    E_bin_z = E_bin_z,
    bin_id = bin_id,
    ord = ord,
    bin_n = bin_n,
    km = km,
    gene_cluster_raw = gene_cluster_raw,
    gene_module = gene_module,
    df_mean = df_mean,
    df_peak = df_peak,
    bin_cor = bin_cor,
    df_gene_peak = df_gene_peak,
    support_pdf = support_pdf
  )
  saveRDS(aux, module_aux_rds)
  message("[SAVE] module aux RDS: ", module_aux_rds)
  
  aux
}

## =========================================================
## [D] Enrichment helpers + dotplot (aligned labels)
## =========================================================
run_gobp_all <- function(genes_sym, universe_sym){
  enrichGO(
    gene          = genes_sym,
    universe      = universe_sym,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
}

run_kegg_all <- function(genes_sym, universe_sym){
  eg  <- sym2entrez(genes_sym)
  ueg <- sym2entrez(universe_sym)
  if (length(eg) < 5 || length(ueg) < 10) return(NULL)
  enrichKEGG(
    gene          = eg,
    universe      = ueg,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
}

run_term2gene_all <- function(genes_sym, universe_sym, term2gene){
  genes_sym <- unique(genes_sym[!is.na(genes_sym) & genes_sym != ""])
  if (length(genes_sym) < 5) return(NULL)
  enricher(
    gene          = genes_sym,
    universe      = universe_sym,
    TERM2GENE     = term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
}

## MSigDB term2gene set (H + C2:CP + lung/cancer CGP subset)
ms_h <- msigdbr(species="Homo sapiens", category="H") %>%
  transmute(term = paste0("H:", gs_name), gene = gene_symbol) %>% distinct()

ms_cp <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP") %>%
  transmute(term = paste0("CP:", gs_name), gene = gene_symbol) %>% distinct()

kw_cgp_lung <- "lung|nsclc|sclc|luad|lusc|adenocarcinoma|carcinoma|egfr|kras|alk|met|erbb|tp53|myc"
ms_cgp_all <- msigdbr(species="Homo sapiens", category="C2", subcategory="CGP") %>%
  distinct(gs_name, gene_symbol)

TOPN_CGP <- 150
cgp_hit_sets <- ms_cgp_all %>%
  count(gs_name, name="n_genes") %>%
  filter(str_detect(gs_name, regex(kw_cgp_lung, ignore_case=TRUE))) %>%
  arrange(desc(n_genes)) %>%
  slice_head(n = TOPN_CGP)

ms_cgp_lung <- ms_cgp_all %>%
  semi_join(cgp_hit_sets, by="gs_name") %>%
  transmute(term = paste0("CGP:", gs_name), gene = gene_symbol) %>%
  distinct()

msig_term2gene <- bind_rows(ms_h, ms_cp, ms_cgp_lung) %>%
  filter(!is.na(gene), gene != "") %>% distinct()

msig_c6 <- msigdbr(species="Homo sapiens", category="C6") %>%
  transmute(term = paste0("C6:", gs_name), gene = gene_symbol) %>%
  distinct()

make_colored_label_vec <- function(term, src, source_cols) {
  term <- as.character(term)
  src  <- as.character(src)
  cols <- source_cols[src]
  cols[is.na(cols)] <- "black"
  sprintf("<span style='color:%s'>%s</span>", cols, term)
}

make_module_enrichment_dotplot_aligned <- function(
    outdir, tag, CACHE_SALT,
    modules,
    source_levels, source_cols,
    q_try = c(0.05,0.10,0.20,0.30,0.50,1.0),
    topK_per_module = 10,
    term_wrap = 55,
    width = 9.0,
    height = 9.3,
    ring_stroke  = 1.2,
    inner_stroke = 0.20,
    inner_scale  = 0.70
) {
  suppressPackageStartupMessages({
    library(ggtext)
    library(viridisLite)
  })
  
  subdir <- file.path(outdir, tag)
  rds_in <- file.path(subdir, sprintf("module_enrich_df_ALL_%s_%s.rds", tag, CACHE_SALT))
  stopifnot(file.exists(rds_in))
  
  enrich_df <- readRDS(rds_in)
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    warning("[", tag, "] enrich_df is empty -> skip dotplot.")
    return(invisible(NULL))
  }
  
  df0 <- enrich_df %>%
    filter(!is.na(Term), Term != "", is.finite(p.adjust)) %>%
    mutate(
      Source_raw = as.character(Source),
      Source = case_when(
        Source_raw %in% c("C6","c6") ~ "C6",
        Source_raw %in% c("MSigDB","MSIGDB","MsigDB","msigdb") ~ "MSigDB",
        Source_raw %in% c("KEGG","kegg") ~ "KEGG",
        Source_raw %in% c("GOBP","GO:BP","GO_BP","gobp","GO.BP","GOBP_BP") ~ "GOBP",
        TRUE ~ Source_raw
      ),
      Module   = factor(as.character(Module), levels = modules),
      Source   = factor(Source, levels = source_levels),
      Term0    = as.character(Term),
      TermW    = str_wrap(Term0, width = term_wrap),
      log10FDR = -log10(pmax(p.adjust, 1e-300))
    )
  
  ## adaptive pick
  df_pick <- NULL
  for (q in q_try) {
    tmp <- df0 %>%
      filter(is.finite(p.adjust), p.adjust <= q) %>%
      group_by(Module) %>%
      arrange(p.adjust, desc(log10FDR)) %>%
      slice_head(n = topK_per_module) %>%
      ungroup()
    if (nrow(tmp) > 0) { df_pick <- tmp; break }
  }
  if (is.null(df_pick) || nrow(df_pick) == 0) {
    warning("[", tag, "] no terms selected (even with q_try ladder).")
    return(invisible(NULL))
  }
  
  df2 <- df_pick %>%
    mutate(
      TermW_chr = as.character(TermW),
      TermLabPlain = TermW_chr,
      TermLab = make_colored_label_vec(TermW_chr, Source, source_cols)
    )
  
  ## y ordering by earliest module appearance
  term_order_tbl <- df2 %>%
    mutate(Module_num = match(as.character(Module), modules)) %>%
    group_by(TermLabPlain) %>%
    summarise(
      first_mod    = min(Module_num, na.rm = TRUE),
      best_q       = min(p.adjust, na.rm = TRUE),
      best_log10   = max(log10FDR, na.rm = TRUE),
      TermLab_any  = TermLab[1],
      .groups = "drop"
    ) %>%
    arrange(first_mod, best_q, desc(best_log10), TermLabPlain)
  
  term_levels <- term_order_tbl$TermLab_any
  df2 <- df2 %>% mutate(TermLab = factor(TermLab, levels = term_levels))
  
  ## size mapping
  size_var <- if ("GeneRatio_num" %in% names(df2) && any(is.finite(df2$GeneRatio_num))) "GeneRatio_num" else "Count"
  v <- df2[[size_var]]
  v <- ifelse(is.finite(v), v, NA_real_)
  if (all(is.na(v))) v <- rep(1, nrow(df2))
  rng <- range(v, na.rm = TRUE)
  
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    df2 <- df2 %>% mutate(size_outer = 7.0, size_inner = 7.0 * inner_scale)
  } else {
    df2 <- df2 %>%
      mutate(
        size_outer = scales::rescale(v, to = c(3.2, 11.0), from = rng),
        size_inner = size_outer * inner_scale
      )
  }
  
  cap_fill <- stats::quantile(df2$log10FDR, 0.90, na.rm = TRUE)
  cap_fill <- max(10, cap_fill)
  
  p <- ggplot(df2, aes(x = Module, y = TermLab)) +
    geom_point(
      aes(size = size_outer, colour = Source),
      shape = 21, fill = "white",
      stroke = ring_stroke, alpha = 1, na.rm = TRUE) +
    geom_point(
      aes(size = size_inner, fill = log10FDR),
      shape = 21, colour = "black",
      stroke = inner_stroke, alpha = 0.98, na.rm = TRUE) +
    scale_colour_manual(values = source_cols, drop = FALSE, name = "Source") +
    scale_fill_gradientn(
      colours = viridisLite::inferno(256),
      limits  = c(0, cap_fill),
      oob     = scales::squish,
      name    = expression(-log[10]~FDR)) +
    scale_size_identity(guide = "none") +
    scale_x_discrete(drop = FALSE) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid   = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(size = 15, face = "bold"),
      axis.text.y  = ggtext::element_markdown(size = 13, face = "bold"),
      legend.position = "right") +
    labs(
      title = sprintf("[%s] Module enrichment \n(top %d/; adaptive q-cut)", tag, topK_per_module),
      x = NULL, y = NULL)
  
  out_pdf <- file.path(subdir, sprintf("DOTPLOT_%s_topTerms_top%d_%s.pdf", tag, topK_per_module, CACHE_SALT))
  ggsave(out_pdf, p, width = width, height = height, useDingbats = FALSE)
  message("[DONE] dotplot: ", out_pdf)
  p
}

run_downstream_for_core <- function(aux, tag) {
  stopifnot(is.list(aux), !is.null(aux$module_tbl))
  subdir <- file.path(outdir, tag)
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  module_tbl <- aux$module_tbl
  modules_local <- aux$modules
  
  ## universe symbols
  if ("gene_short_name" %in% colnames(fData(cds_traj))) {
    universe_sym <- as.character(fData(cds_traj)$gene_short_name)
  } else {
    universe_sym <- as.character(module_tbl$SYMBOL)
  }
  universe_sym <- unique(universe_sym[!is.na(universe_sym) & universe_sym != ""])
  
  module_tbl_enrich <- module_tbl %>%
    mutate(SYMBOL = as.character(SYMBOL)) %>%
    filter(!is.na(SYMBOL), SYMBOL != "")
  
  enrich_raw <- list()
  for (m in modules_local) {
    genes_m <- module_tbl_enrich %>% filter(Module == m) %>% pull(SYMBOL) %>% unique()
    genes_m <- intersect(genes_m, universe_sym)
    
    if (length(genes_m) < 5) {
      enrich_raw[[m]] <- list(GOBP=NULL, KEGG=NULL, MSigDB=NULL, C6=NULL)
      next
    }
    enrich_raw[[m]] <- list(
      GOBP   = tryCatch(run_gobp_all(genes_m, universe_sym), error=function(e) NULL),
      KEGG   = tryCatch(run_kegg_all(genes_m, universe_sym), error=function(e) NULL),
      MSigDB = tryCatch(run_term2gene_all(genes_m, universe_sym, msig_term2gene), error=function(e) NULL),
      C6     = tryCatch(run_term2gene_all(genes_m, universe_sym, msig_c6), error=function(e) NULL)
    )
  }
  
  flatten_enrich <- function(enrich_raw){
    purrr::imap_dfr(enrich_raw, function(lst, mod){
      purrr::imap_dfr(lst, function(obj, src){
        df <- as_df_safe(obj)
        if (is.null(df)) return(NULL)
        df <- tibble::as_tibble(df, rownames = "RowKey")
        
        tibble::tibble(
          Core     = tag,
          Module   = mod,
          Source   = src,
          TermID   = if ("ID" %in% names(df)) df$ID else NA_character_,
          Term     = if ("Description" %in% names(df)) df$Description else NA_character_,
          p.adjust = as.numeric(df$p.adjust),
          Count    = as.numeric(df$Count),
          GeneRatio= df$GeneRatio,
          BgRatio  = df$BgRatio,
          geneID   = df$geneID
        ) %>%
          mutate(
            Module   = factor(Module, levels = modules_local),
            Source   = factor(Source, levels = source_levels),
            neglog10 = -log10(pmax(p.adjust, 1e-300)),
            GeneRatio_num = ratio_to_num(GeneRatio)
          )
      })
    })
  }
  
  enrich_df <- flatten_enrich(enrich_raw)
  saveRDS(enrich_df, file.path(subdir, sprintf("module_enrich_df_ALL_%s_%s.rds", tag, CACHE_SALT)))
  message("[SAVE] enrich_df: ", file.path(subdir, sprintf("module_enrich_df_ALL_%s_%s.rds", tag, CACHE_SALT)))
  
  p_dot <- make_module_enrichment_dotplot_aligned(
    outdir = outdir,
    tag = tag,
    CACHE_SALT = CACHE_SALT,
    modules = modules_local,
    source_levels = source_levels,
    source_cols = source_cols,
    topK_per_module = 10,
    term_wrap = 55,
    width = 9.0,
    height = 9.3,
    ring_stroke = 1.2,
    inner_stroke = 0.20,
    inner_scale = 0.70
  )
  
  list(enrich_df = enrich_df, p_dot = p_dot)
}

## =========================================================
## [E] DEG correlation heatmap helpers
## =========================================================
.assert_de_has_cols <- function(de, gene_col, lfc_col, q_col) {
  miss <- setdiff(c(gene_col, lfc_col, q_col), colnames(de))
  if (length(miss) > 0) stop("DE table missing columns: ", paste(miss, collapse=", "))
}

get_whole_deg_genes <- function(de_all, samples, q_cut,
                                gene_col="SYMBOL", q_col="p.adjust",
                                lfc_col="log2FC_Malig_vs_AT2") {
  genes <- unique(unlist(lapply(samples, function(sid) {
    de <- de_all[[sid]]
    if (is.null(de) || nrow(de)==0) return(character(0))
    de <- as.data.frame(de)
    .assert_de_has_cols(de, gene_col, lfc_col, q_col)
    de %>%
      filter(is.finite(.data[[q_col]]), .data[[q_col]] <= q_cut) %>%
      pull(.data[[gene_col]]) %>% as.character()
  }), use.names = FALSE))
  genes <- genes[!is.na(genes) & genes != ""]
  unique(genes)
}

make_lfc_matrix <- function(de_all, samples, genes,
                            gene_col="SYMBOL", lfc_col="log2FC_Malig_vs_AT2",
                            q_col="p.adjust",
                            agg_fun=c("first","mean")) {
  agg_fun <- match.arg(agg_fun)
  
  mat <- sapply(samples, function(sid) {
    de <- de_all[[sid]]
    if (is.null(de) || nrow(de)==0) return(rep(NA_real_, length(genes)))
    de <- as.data.frame(de)
    .assert_de_has_cols(de, gene_col, lfc_col, q_col)
    
    g <- as.character(de[[gene_col]])
    x <- suppressWarnings(as.numeric(de[[lfc_col]]))
    
    if (agg_fun == "mean") {
      tmp <- tapply(x, g, function(v) mean(v[is.finite(v)], na.rm=TRUE))
      v <- unname(tmp[genes])
      v[is.nan(v)] <- NA_real_
      return(as.numeric(v))
    } else {
      m <- match(genes, g)
      v <- x[m]
      return(as.numeric(v))
    }
  })
  
  rownames(mat) <- genes
  colnames(mat) <- samples
  mat
}

select_side_genes <- function(lfc_mat, side=c("AT2","Malig"), min_support=1) {
  side <- match.arg(side)
  if (side == "AT2") {
    keep <- rowSums(is.finite(lfc_mat) & lfc_mat < 0) >= min_support
  } else {
    keep <- rowSums(is.finite(lfc_mat) & lfc_mat > 0) >= min_support
  }
  rownames(lfc_mat)[keep]
}

compute_cor_mat <- function(lfc_mat, min_non_na_per_gene = 2) {
  keep_gene <- rowSums(is.finite(lfc_mat)) >= min_non_na_per_gene
  lfc2 <- lfc_mat[keep_gene, , drop=FALSE]
  if (nrow(lfc2) < 30) warning("[WARN] very few genes after NA-filter: n=", nrow(lfc2))
  cor(lfc2, use = "pairwise.complete.obs", method = "pearson")
}

.make_dist_from_cor <- function(cor_mat) {
  d <- 1 - cor_mat
  d[!is.finite(d)] <- 1
  diag(d) <- 0
  as.dist(d)
}

DEFAULT_COLFUN <- function(n) grDevices::colorRampPalette(c("#2C7BB6", "#FFFFFF", "#D7191C"))(n)
SEQ_WR_COLFUN  <- function(n) grDevices::colorRampPalette(c("#FFFFFF", "#D7191C"))(n)

plot_corr_heatmap_fixedrange <- function(cor_mat, main, out_prefix,
                                         vmin = 0.2, vmax = 1.0,
                                         n_breaks = 201,
                                         colors = NULL,
                                         palette = c("bwr", "wr"),
                                         cluster_method = "average",
                                         border_color = "grey75",
                                         show_rownames = TRUE,
                                         show_colnames = TRUE,
                                         fontsize_axis = 13) {
  palette <- match.arg(palette)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  
  cor_mat <- as.matrix(cor_mat)
  cor_clip <- cor_mat
  cor_clip[!is.finite(cor_clip)] <- NA_real_
  cor_clip <- pmin(pmax(cor_clip, vmin), vmax)
  
  breaks <- seq(vmin, vmax, length.out = n_breaks)
  
  if (is.null(colors)) {
    if (palette == "bwr") colors <- DEFAULT_COLFUN(length(breaks) - 1)
    else colors <- SEQ_WR_COLFUN(length(breaks) - 1)
  } else {
    if (length(colors) != (length(breaks) - 1)) {
      stop(sprintf("[ERROR] colors length (%d) must equal length(breaks)-1 (%d).",
                   length(colors), length(breaks) - 1))
    }
  }
  
  dist_mat <- .make_dist_from_cor(cor_mat)
  main_txt <- sprintf("%s (clip %.2f–%.2f)", main, vmin, vmax)
  
  grDevices::pdf(paste0(out_prefix, ".pdf"), width = 4.3, height = 3.9, useDingbats = FALSE)
  p <- pheatmap(
    cor_clip,
    color  = colors,
    breaks = breaks,
    clustering_distance_rows = dist_mat,
    clustering_distance_cols = dist_mat,
    clustering_method = cluster_method,
    border_color = border_color,
    main = main_txt,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    fontsize_row = fontsize_axis,
    fontsize_col = fontsize_axis,
    fontsize     = fontsize_axis
  )
  print(p)
  grDevices::dev.off()
  invisible(TRUE)
}

## =========================================================
## [F] Seurat module scores (z-score across cells per gene)
## =========================================================
.get_assay_mat <- function(seu, assay="RNA") {
  m <- tryCatch(Seurat::GetAssayData(seu, assay=assay, layer="data"),
                error = function(e) NULL)
  if (is.null(m) || nrow(m) == 0) {
    m <- Seurat::GetAssayData(seu, assay=assay, layer="counts")
  }
  m
}

compute_module_scores_z <- function(
    seu,
    module_gene_sets,          # named list: M1..M4 -> SYMBOL vector
    modules = paste0("M",1:4),
    assay_use = "RNA",
    min_genes_per_module = 10
) {
  stopifnot(all(modules %in% names(module_gene_sets)))
  
  X <- .get_assay_mat(seu, assay = assay_use)
  if (!inherits(X, "dgCMatrix")) X <- as(X, "dgCMatrix")
  
  module_gene_sets2 <- lapply(module_gene_sets[modules], function(g) intersect(unique(g), rownames(X)))
  
  n_per_mod <- sapply(module_gene_sets2, length)
  if (any(n_per_mod < min_genes_per_module)) {
    warning("[WARN] Some modules have < min_genes_per_module:\n",
            paste(names(n_per_mod), n_per_mod, sep=": ", collapse="\n"))
  }
  
  all_genes <- unique(unlist(module_gene_sets2, use.names = FALSE))
  Xsub <- X[all_genes, , drop = FALSE]
  
  gmean <- Matrix::rowMeans(Xsub)
  ex2   <- Matrix::rowMeans(Xsub^2)
  gvar  <- ex2 - gmean^2
  gvar[gvar < 0] <- 0
  gsd <- sqrt(gvar)
  gsd[!is.finite(gsd) | gsd == 0] <- NA_real_
  
  score_list <- list()
  
  for (m in modules) {
    genes_m <- module_gene_sets2[[m]]
    if (length(genes_m) < min_genes_per_module) {
      score_list[[paste0("Score_", m)]] <- rep(NA_real_, ncol(X))
      next
    }
    
    mu <- gmean[genes_m]
    sd <- gsd[genes_m]
    ok <- is.finite(sd) & is.finite(mu)
    
    genes_ok <- genes_m[ok]
    mu <- mu[ok]
    sd <- sd[ok]
    
    if (length(genes_ok) < min_genes_per_module) {
      score_list[[paste0("Score_", m)]] <- rep(NA_real_, ncol(X))
      next
    }
    
    Xm <- Xsub[genes_ok, , drop = FALSE]
    
    zsum <- rep(0, ncol(Xm))
    ng   <- 0L
    for (i in seq_len(nrow(Xm))) {
      xi <- as.numeric(Xm[i, ])
      zi <- (xi - mu[i]) / sd[i]
      zsum <- zsum + zi
      ng <- ng + 1L
    }
    score_list[[paste0("Score_", m)]] <- zsum / max(1L, ng)
  }
  
  as.data.frame(score_list, row.names = colnames(X))
}

get_dom_module <- function(df, score_cols, modules, tie_method = c("first","random")) {
  tie_method <- match.arg(tie_method)
  m <- as.matrix(df[, score_cols, drop=FALSE])
  m <- apply(m, 2, function(x) suppressWarnings(as.numeric(as.character(x))))
  
  dom_idx <- apply(m, 1, function(x) {
    if (all(!is.finite(x))) return(NA_integer_)
    mx <- max(x, na.rm = TRUE)
    wh <- which(x == mx)
    if (length(wh) == 1) return(wh)
    if (tie_method == "random") return(sample(wh, 1))
    wh[1]
  })
  out <- rep(NA_character_, nrow(df))
  ok <- !is.na(dom_idx)
  out[ok] <- modules[dom_idx[ok]]
  out
}

plot_fraction_dominant_module <- function(
    seu,
    group_col = "celltype_TissueDiseaseState",
    dom_col   = "DominantModule",
    modules   = paste0("M",1:4),
    module_colors = c("M1"="#4E79A7","M2"="#F28E2B","M3"="#59A14F","M4"="#E15759"),
    out_pdf = NULL,
    width = 6.8, height = 4.2
) {
  md <- seu@meta.data
  stopifnot(group_col %in% colnames(md), dom_col %in% colnames(md))
  
  df <- md %>%
    filter(!is.na(.data[[group_col]]), !is.na(.data[[dom_col]])) %>%
    mutate(
      Group = factor(.data[[group_col]]),
      Dom   = factor(as.character(.data[[dom_col]]), levels = modules)
    ) %>%
    count(Group, Dom, name="n") %>%
    group_by(Group) %>%
    mutate(frac = n / sum(n),
           n_group = sum(n)) %>%
    ungroup()
  
  p <- ggplot(df, aes(x = Group, y = frac, fill = Dom)) +
    geom_col(color="black", linewidth=0.25) +
    scale_fill_manual(values = module_colors, drop=FALSE) +
    theme_bw(base_size=14) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle=30, hjust=1)) +
    labs(x=NULL, y="Fraction", fill="Dominant module",
         title="Dominant module fraction by celltype_TissueDiseaseState") +
    geom_text(
      data = df %>% distinct(Group, n_group, .keep_all = TRUE),
      aes(x = Group, y = 1.02, label = paste0("n=", n_group)),
      inherit.aes = FALSE, size = 3.6) +
    coord_cartesian(ylim = c(0, 1.08))
  
  if (!is.null(out_pdf)) {
    grDevices::pdf(out_pdf, width=width, height=height, useDingbats=FALSE)
    print(p); grDevices::dev.off()
    message("[DONE] ", out_pdf)
  }
  p
}

## =========================================================
## [G] RUN: core6 only
## =========================================================
TAG <- "core6"
SUBDIR <- file.path(outdir, TAG)
dir.create(SUBDIR, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------
## 1) Module definition (core6)
## ---------------------------------------------------------
if (RUN_MODULE_DEFINE) {
  aux6 <- define_modules_from_core(
    cds_all = cds_traj,
    core_genes = core_genes_6,
    tag = TAG,
    outdir = outdir,
    CACHE_SALT = CACHE_SALT,
    N_BIN = N_BIN,
    min_bin_cells = MIN_BIN_CELLS,
    K_MODULES = NUM_MODULES,
    seed = KMEANS_SEED,
    nstart = KMEANS_NSTART,
    itmax = KMEANS_ITMAX,
    force = TRUE,
    module_colors = module_colors)
} else {
  aux_rds <- file.path(SUBDIR, sprintf("module_aux_fixed_%s_%s.rds", TAG, CACHE_SALT))
  stopifnot(file.exists(aux_rds))
  aux6 <- readRDS(aux_rds)
  if (is.null(aux6$module_tbl)) {
    tbl_rds <- file.path(SUBDIR, sprintf("module_tbl_fixed_%s_%s.rds", TAG, CACHE_SALT))
    aux6$module_tbl <- readRDS(tbl_rds)
  }
  message("[LOAD] aux6 from RDS: ", aux_rds)
}

## module genes list (SYMBOL) save
module_genes_core6 <- split(aux6$module_tbl$SYMBOL, aux6$module_tbl$Module)
module_genes_core6 <- lapply(module_genes_core6, function(x) unique(x[!is.na(x) & x != ""]))

save_rds_safe(module_genes_core6, file.path(outdir, TAG, sprintf("module_genes_%s_%s.rds", TAG, CACHE_SALT)))

for (m in names(module_genes_core6)) {
  write.table(
    module_genes_core6[[m]],
    file = file.path(outdir, TAG, sprintf("genes_%s_%s_%s.txt", m, TAG, CACHE_SALT)),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
}

## ---------------------------------------------------------
## 2) Enrichment (core6)
## ---------------------------------------------------------
if (RUN_ENRICHMENT) {
  res6 <- run_downstream_for_core(aux6, TAG)
}

## ---------------------------------------------------------
## 3) DEG correlation heatmaps (ALL DEG vs core6 DEG)
## ---------------------------------------------------------
if (RUN_DEG_CORR) {
  SAMPLES_CORR <- names(de_all)
  stopifnot(length(SAMPLES_CORR) >= 3)
  
  whole_genes <- get_whole_deg_genes(
    de_all, SAMPLES_CORR, q_cut = Q_CUT_WHOLE,
    gene_col = GENE_COL, q_col = Q_COL, lfc_col = LFC_COL)
  message("[INFO] whole DEG genes (union, q<=", Q_CUT_WHOLE, "): ", length(whole_genes))
  
  lfc_whole <- make_lfc_matrix(
    de_all, SAMPLES_CORR, whole_genes,
    gene_col = GENE_COL, lfc_col = LFC_COL, q_col = Q_COL, agg_fun = "first")
  
  genes_whole_AT2   <- select_side_genes(lfc_whole, side="AT2",   min_support=MIN_SUPPORT_SIDE_WHOLE)
  genes_whole_Malig <- select_side_genes(lfc_whole, side="Malig", min_support=MIN_SUPPORT_SIDE_WHOLE)
  
  cor_whole_AT2   <- compute_cor_mat(lfc_whole[genes_whole_AT2, , drop=FALSE])
  cor_whole_Malig <- compute_cor_mat(lfc_whole[genes_whole_Malig, , drop=FALSE])
  
  plot_corr_heatmap_fixedrange(
    cor_mat = cor_whole_AT2,
    main = "[ALL DEG] AT2",
    out_prefix = file.path(outdir, paste0("corr_at2_all_", CACHE_SALT)),
    vmin = SCALE_MIN, vmax = SCALE_MAX, n_breaks = N_BREAKS,
    palette = "bwr")
  
  plot_corr_heatmap_fixedrange(
    cor_mat = cor_whole_Malig,
    main = "[ALL DEG] Malignant",
    out_prefix = file.path(outdir, paste0("corr_malignant_all_", CACHE_SALT)),
    vmin = SCALE_MIN, vmax = SCALE_MAX, n_breaks = N_BREAKS,
    palette = "bwr")
  
  core6_symbol <- unique(as.character(core6_genes))
  core6_symbol <- core6_symbol[core6_symbol %in% whole_genes]
  message("[INFO] core6 genes used (core6 ∩ wholeDEG): ", length(core6_symbol))
  
  lfc_core6 <- make_lfc_matrix(
    de_all, SAMPLES_CORR, core6_symbol,
    gene_col = GENE_COL, lfc_col = LFC_COL, q_col = Q_COL, agg_fun = "first")
  
  genes_core6_AT2   <- select_side_genes(lfc_core6, side="AT2",   min_support=MIN_SUPPORT_SIDE_CORE)
  genes_core6_Malig <- select_side_genes(lfc_core6, side="Malig", min_support=MIN_SUPPORT_SIDE_CORE)
  
  cor_core6_AT2   <- compute_cor_mat(lfc_core6[genes_core6_AT2, , drop=FALSE])
  cor_core6_Malig <- compute_cor_mat(lfc_core6[genes_core6_Malig, , drop=FALSE])
  
  plot_corr_heatmap_fixedrange(
    cor_mat = cor_core6_AT2,
    main = "[Core6 DEG] AT2",
    out_prefix = file.path(outdir, paste0("corr_at2_core_", CACHE_SALT)),
    vmin = SCALE_MIN, vmax = SCALE_MAX, n_breaks = N_BREAKS,
    palette = "bwr")
  
  plot_corr_heatmap_fixedrange(
    cor_mat = cor_core6_Malig,
    main = "[Core6 DEG] Malignant",
    out_prefix = file.path(outdir, paste0("corr_malignant_core_", CACHE_SALT)),
    vmin = SCALE_MIN, vmax = SCALE_MAX, n_breaks = N_BREAKS,
    palette = "bwr")
  
  message("=== [DONE] 4 correlation heatmaps saved to OUTDIR ===")
}

## ---------------------------------------------------------
## 4) Attach M1–M4 scores + dominant module to Seurat
## ---------------------------------------------------------
if (RUN_SEURAT_SCORE) {
  stopifnot(CTTDS_COL %in% colnames(seu@meta.data))
  
  score_df <- compute_module_scores_z(
    seu = seu,
    module_gene_sets = module_genes_core6,
    modules = modules,
    assay_use = "RNA",
    min_genes_per_module = 10)
  
  for (sc in colnames(score_df)) seu[[sc]] <- score_df[[sc]]
  
  score_cols <- paste0("Score_", modules)
  seu$DominantModule <- get_dom_module(seu@meta.data, score_cols, modules, tie_method="first")
  seu$DominantModule <- factor(seu$DominantModule, levels = modules)
  
  save_rds_safe(seu, file.path(outdir, sprintf("seu_with_ModuleScores_%s.rds", CACHE_SALT)))
}

## ---------------------------------------------------------
## 5) Plot-1: Violin module score by AT2 vs Malignant
## ---------------------------------------------------------
if (RUN_PLOTS_VIOLIN) {
  stopifnot("celltype4" %in% colnames(seu@meta.data))
  score_cols <- paste0("Score_", modules)
  
  seu$ct2_simple <- as.character(seu$celltype4)
  seu$ct2_simple <- dplyr::case_when(
    grepl("^AT2\\b", seu$ct2_simple) ~ "AT2",
    grepl("(?i)^malignant|malignant", seu$ct2_simple, perl = TRUE) ~ "Malignant",
    TRUE ~ NA_character_
  )
  seu$ct2_simple <- factor(seu$ct2_simple, levels = c("AT2","Malignant"))
  
  df_long <- seu@meta.data %>%
    select(ct2_simple, all_of(score_cols)) %>%
    filter(!is.na(ct2_simple)) %>%
    pivot_longer(cols = all_of(score_cols),
                 names_to = "Module", values_to = "Score") %>%
    mutate(Module = gsub("^Score_", "", Module),
           Module = factor(Module, levels = modules))
  
  strip_bg <- NULL
  if (VIOLIN_USE_GGH4X && requireNamespace("ggh4x", quietly = TRUE)) {
    strip_bg <- ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = unname(module_colors),
        color = "black",
        linewidth = 0.3 ))
  }
  
  p_violin <- ggplot(df_long, aes(x = ct2_simple, y = Score)) +
    geom_violin(aes(fill = Module), trim = TRUE, color = "black", linewidth = 0.25) +
    geom_boxplot(width = 0.13, outlier.size = 0.2, linewidth = 0.25) +
    { if (!is.null(strip_bg)) ggh4x::facet_wrap2(~ Module, nrow = 2, strip = strip_bg)
      else facet_wrap(~ Module, nrow = 2) } +
    scale_fill_manual(values = module_colors, drop = FALSE, guide = "none") +
    theme_bw(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      axis.text.x = element_text(face = "plain", color = "black", size = 14),
      axis.text.y = element_text(face = "plain", color = "black", size = 14.7),
      strip.background = element_blank()) +
    labs(x = "Cell type", y = "Module score", title = "Module score by cell type")
  
  out_pdf1 <- file.path(outdir, sprintf("Violin_ModuleScore_AT2vsMalig_%s.pdf", CACHE_SALT))
  ggsave(out_pdf1, p_violin, width = 4.2, height = 4.2, useDingbats = FALSE)
  message("[DONE] ", out_pdf1)
}

## ---------------------------------------------------------
## 6) Plot-2: Top genes per module (specificity-weighted)
## ---------------------------------------------------------
if (RUN_PLOTS_TOPGENES) {
  stopifnot("ct2_simple" %in% colnames(seu@meta.data))
  X <- tryCatch(Seurat::GetAssayData(seu, assay = "RNA", layer = "data"),
                error = function(e) NULL)
  if (is.null(X) || nrow(X) == 0) X <- Seurat::GetAssayData(seu, assay = "RNA", layer = "counts")
  if (!inherits(X, "dgCMatrix")) X <- as(X, "dgCMatrix")
  
  mt <- aux6$module_tbl %>%
    mutate(SYMBOL = as.character(SYMBOL),
           Module = as.character(Module)) %>%
    filter(!is.na(SYMBOL), SYMBOL != "", Module %in% modules) %>%
    distinct(SYMBOL, Module)
  
  genes_use <- intersect(mt$SYMBOL, rownames(X))
  mt2 <- mt %>% filter(SYMBOL %in% genes_use)
  
  cells_at2   <- rownames(seu@meta.data)[seu$ct2_simple == "AT2"]
  cells_malig <- rownames(seu@meta.data)[seu$ct2_simple == "Malignant"]
  stopifnot(length(cells_at2) > 20, length(cells_malig) > 20)
  
  eps <- 1e-6
  mu_at2   <- Matrix::rowMeans(X[genes_use, cells_at2, drop = FALSE])
  mu_malig <- Matrix::rowMeans(X[genes_use, cells_malig, drop = FALSE])
  
  df_gene <- data.frame(
    SYMBOL   = genes_use,
    mu_AT2   = as.numeric(mu_at2),
    mu_Malig = as.numeric(mu_malig),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      max_mu = pmax(mu_AT2, mu_Malig),
      logFC  = log2((mu_Malig + eps) / (mu_AT2 + eps)),
      specificity = abs(mu_Malig - mu_AT2) / (mu_Malig + mu_AT2 + eps),
      score_spec  = abs(logFC) * max_mu * (1 + specificity)
    ) %>%
    inner_join(mt2, by = "SYMBOL") %>%
    mutate(Module = factor(Module, levels = modules))
  
  df_top <- df_gene %>%
    group_by(Module) %>%
    arrange(desc(score_spec)) %>%
    slice_head(n = TOPGENE_N) %>%
    ungroup() %>%
    group_by(Module) %>%
    mutate(SYMBOL_f = factor(SYMBOL, levels = rev(SYMBOL))) %>%
    ungroup()
  
  strip_bg2 <- NULL
  if (requireNamespace("ggh4x", quietly = TRUE)) {
    strip_bg2 <- ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = unname(module_colors),
        color = "black",
        linewidth = 0.3))
  }
  
  p_top <- ggplot(df_top, aes(x = score_spec, y = SYMBOL_f)) +
    geom_col(aes(fill = Module), color = "black", linewidth = 0.25) +
    { if (!is.null(strip_bg2)) ggh4x::facet_wrap2(~ Module, scales = "free_y", ncol = 2, strip = strip_bg2)
      else facet_wrap(~ Module, scales = "free_y", ncol = 2) } +
    scale_fill_manual(values = module_colors, drop = FALSE, guide = "none") +
    theme_bw(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      axis.text.x = element_text(face = "bold", color = "black", size = 14.7),
      axis.text.y = element_text(face = "italic", color = "black", size = 14.7),
      strip.background = element_blank()) +
    labs(
      x = "Specificity-weighted contribution score",
      y = NULL,
      title = "Top genes per module")
  
  out_pdf2 <- file.path(outdir, sprintf("TopGenesPerModule_contribution_AT2vsMalig_%s.pdf", CACHE_SALT))
  ggsave(out_pdf2, p_top, width = 7.8, height = 9.6, useDingbats = FALSE)
  message("[DONE] ", out_pdf2)
}

## ---------------------------------------------------------
## 7) Fraction plot: DominantModule by celltype_TissueDiseaseState
## ---------------------------------------------------------
if (RUN_FRACTION_PLOTS) {
  out_pdf_frac <- file.path(outdir, sprintf("Fraction_DominantModule_by_%s_%s.pdf", "celltypeTDS", CACHE_SALT))
  plot_fraction_dominant_module(
    seu,
    group_col = CTTDS_COL,
    dom_col   = "DominantModule",
    modules   = modules,
    module_colors = module_colors,
    out_pdf = out_pdf_frac)
}

## ---------------------------------------------------------
## 8) UMAP plots (DominantModule + Score_Mi)
## ---------------------------------------------------------
if (RUN_UMAP_PLOTS) {
  if (!("umap" %in% names(seu@reductions))) {
    warning("[WARN] seu has no UMAP reduction named 'umap'. Skipping UMAP plots.")
  } else {
    p1 <- DimPlot(seu, group.by = "DominantModule", cols = module_colors, reduction = "umap") +
      ggtitle("Dominant module (M1–M4)")
    out1 <- file.path(outdir, sprintf("UMAP_DominantModule_%s.pdf", CACHE_SALT))
    ggsave(out1, p1, width=5.8, height=4.8, useDingbats=FALSE)
    message("[DONE] ", out1)
    
    for (sc in paste0("Score_", modules)) {
      if (!sc %in% colnames(seu@meta.data)) next
      p <- FeaturePlot(seu, features = sc, reduction = "umap") + ggtitle(sc)
      out <- file.path(outdir, sprintf("UMAP_%s_%s.pdf", sc, CACHE_SALT))
      ggsave(out, p, width=5.8, height=4.8, useDingbats=FALSE)
      message("[DONE] ", out)
    }
  }
}

message("=== [ALL DONE] core6-only pipeline completed ===")
