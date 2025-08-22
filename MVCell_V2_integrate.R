#!/usr/bin/env Rscript

# MVCell standalone run script — Part 2 (annotation)
# Runs SingleR, ScType, and scATOMIC, then exports per-cell annotations and scores.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
})

set.seed(1)

# ----- CLI / paths / parallel -----
all_args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("--file=", all_args, value = TRUE)
script.dir <- if (length(file.arg) > 0) dirname(normalizePath(sub("--file=", "", file.arg))) else getwd()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript script_part2.R /path/to/file.h5ad [num_cores]")

h5ad_path <- normalizePath(args[1])
h5ad_dir  <- dirname(h5ad_path)
cat("H5AD file path:", h5ad_path, "\n")

benchmark.num = if (length(args) == 3) as.integer(args[3]) else NA

# setting parameters. adding benchmark numbers for control
#paths
out_dir  <- dirname(h5ad_path)

#filename base
out.base <- tools::file_path_sans_ext(basename(h5ad_path))
if (!is.na(benchmark.num)) {
  out.base <- paste0(out.base, "_", benchmark.num)
} else {
  out.base <- out.base
}

#names
out_rds  <- file.path(out_dir, paste0(out.base, "_sce.rds"))
out_mvcell = file.path(out_dir, "mvcell_results/cell_type_predictions.csv")
out_raw_annot = file.path(out_dir, "celltype_annotate/merged_with_pipeline.tsv")
out_raw_annlevel = file.path(out_dir, "celltype_annotate/merged_annlevel_transfered.tsv")


# ---- setting cores ----
available_cores <- parallel::detectCores(logical = TRUE)
ncore <- if (length(args) >= 2) as.integer(args[2]) else 0
if (is.na(ncore) || ncore < 0) ncore <- 0

if (ncore == 0) {
  # default is 4, but don't exceed available_cores - 2
  ncore <- min(4L, max(1L, available_cores - 2L))
}

# final safety clamp
ncore <- min(max(1L, ncore), max(1L, available_cores - 2L))

if (.Platform$OS.type == "windows") {
  BPPARAM <- BiocParallel::SnowParam(workers = ncore)
} else {
  BPPARAM <- BiocParallel::MulticoreParam(workers = ncore)
}
BiocParallel::register(BPPARAM)
message(sprintf("Using %d cores (of %d available)", ncore, available_cores))

# ----- Load MVCell results -----

if (!file.exists(out_mvcell)) stop("Missing MVCell results: ", out_mvcell)
mvcell.df = read.csv(out_mvcell)

if (!file.exists(out_rds)) stop("Missing preprocessed RDS: ", out_rds)
sce <- readRDS(out_rds)

meta = colData(sce)
if (!"barcode" %in% colnames(meta)) meta$barcode <- rownames(meta)

meta.m = merge(meta, mvcell.df, all.x = T, by = "barcode")
rownames(meta.m) = meta.m$barcode
meta.m = meta.m[colnames(sce),]
colData(sce) = meta.m

if (!"pca_logcounts" %in% assayNames(sce)) {
  dec.sce <- modelGeneVar(sce)
  hvgs = getTopHVGs(dec.sce, n=2000)
  
  sce = runPCA(sce, exprs_values = "logcounts", ncomponents = 50, subset_row = hvgs, scale = TRUE, name = "pca_logcounts")
  sce = runUMAP(sce, dimred = "pca_logcounts", ncomponents = 2, name = "umap_logcounts")
  
}


plotReducedDim(sce, dimred = "pca_logcounts", colour_by = "ann_level_1_prediction")
plotReducedDim(sce, dimred = "umap_logcounts", colour_by = "ann_level_1_prediction")
plotReducedDim(sce, dimred = "pca_logcounts", colour_by = "ann_level_2_prediction")

plotReducedDim(sce, dimred = "umap_logcounts", colour_by = "celltypist_ann1")


plotReducedDim(sce, dimred = "umap_logcounts", colour_by = "cell_type")
plotReducedDim(sce, dimred = "pca_logcounts", colour_by = "cell_type")

#testing more data
test.data = as.data.frame(data.table::fread(out_raw_annlevel))
meta.m = merge(meta.m, test.data, all.x = T, by = "barcode")
rownames(meta.m) = meta.m$barcode
meta.m = meta.m[colnames(sce),]
colData(sce) = meta.m0
table(meta.m$cell_type, meta.m$celltypist_ann1)
table(meta.m$cell_type, meta.m$scimilarity_ann1)
table(meta.m$cell_type, meta.m$scATOMIC_ann1)
table(meta.m$cell_type, meta.m$singleR_ann1)
table(meta.m$cell_type, meta.m$sctype_ann1)
table(meta.m$cell_type, meta.m$ann_level_1_prediction)

plotReducedDim(sce, dimred = "pca_logcounts", colour_by = "celltypist_annot")


