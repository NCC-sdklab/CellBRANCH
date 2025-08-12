#!/usr/bin/env Rscript

# MVCell standalone run script — Part 2 (annotation)
# Runs SingleR, ScType, and scATOMIC, then exports per-cell annotations and scores.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(BiocParallel)
  library(SingleR)
  library(celldex)
  library(openxlsx)
  library(HGNChelper)
  library(reticulate)
})

set.seed(1)

# ----- Helpers -----
log_step <- function(fmt, ...) message(sprintf(fmt, ...))

score.to.celltype <- function(celltype.score, top.rank = 3){
  # Safeguard when fewer rows than requested top.rank
  kmax <- pmin(top.rank, nrow(celltype.score))
  query.celltype <- matrix(NA_character_, ncol = ncol(celltype.score), nrow = kmax)
  colnames(query.celltype) <- colnames(celltype.score)
  
  query.celltype.score <- matrix(0, ncol = ncol(celltype.score), nrow = kmax)
  colnames(query.celltype.score) <- colnames(celltype.score)
  
  for (k in seq_len(ncol(query.celltype))) {
    query.table <- sort(celltype.score[, k], decreasing = TRUE)
    top_idx <- seq_len(kmax)
    query.celltype[, k] <- names(query.table)[top_idx]
    query.celltype.score[, k] <- query.table[top_idx]
  }
  list(top.celltype = query.celltype, top.score = query.celltype.score)
}

score.rank <- function(score.df, log.n = FALSE){
  x <- as.matrix(score.df)
  if (log.n) {
    x <- log(x - min(x) + 0.01)
  }
  z <- pnorm(scale(as.vector(x)))
  out <- matrix(z, nrow = nrow(score.df), ncol = ncol(score.df))
  colnames(out) <- colnames(score.df)
  rownames(out) <- rownames(score.df)
  out
}

# ----- CLI / paths / parallel -----
all_args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("--file=", all_args, value = TRUE)
script.dir <- if (length(file.arg) > 0) dirname(normalizePath(sub("--file=", "", file.arg))) else getwd()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript script_part2.R /path/to/file.h5ad [num_cores]")

h5ad_path <- normalizePath(args[1])
h5ad_dir  <- dirname(h5ad_path)
cat("H5AD file path:", h5ad_path, "\n")

out_base <- tools::file_path_sans_ext(basename(h5ad_path))
out_rds  <- file.path(h5ad_dir, paste0(out_base, "_sce.rds"))

available_cores <- parallel::detectCores(logical = TRUE)
ncore <- suppressWarnings(as.integer(if (length(args) >= 2) args[2] else NA))
if (is.na(ncore) || ncore < 1) ncore <- available_cores - 2
ncore <- max(1, min(ncore, max(1, available_cores - 2)))

if (.Platform$OS.type == "windows") {
  BPPARAM <- BiocParallel::SnowParam(workers = ncore)
} else {
  BPPARAM <- BiocParallel::MulticoreParam(workers = ncore)
}
BiocParallel::register(BPPARAM)
log_step("Using %d cores (of %d available)", ncore, available_cores)

# ----- Load preprocessed SCE from Part 1 -----
if (!file.exists(sce_rds)) stop("Missing preprocessed RDS: ", out_rds)
sce <- readRDS(out_rds)

# ----- ScType sources / DB -----
sctype_base <- Sys.getenv("MVCell_SCTYPE_DIR", unset = script.dir)
sc_src1 <- file.path(sctype_base, "sctype", "sctype_score_.R")
sc_src2 <- file.path(sctype_base, "sctype", "gene_sets_prepare.R")
sc_db   <- file.path(sctype_base, "sctype", "ScTypeDB_full.xlsx")
if (!file.exists(sc_src1) || !file.exists(sc_src2)) stop("Missing ScType source scripts under ", sctype_base)
if (!file.exists(sc_db)) stop("Missing ScType DB at ", sc_db)

source(sc_src1)
source(sc_src2)

# ----- Reticulate / scATOMIC env (best-effort) -----
ok_env <- FALSE
try({
  cb <- reticulate::conda_binary()
  cb.base <- sub("/bin/(conda|mamba)$", "", cb)
  ret.path <- file.path(cb.base, "envs", "r-reticulate-mvcell")
  if (dir.exists(ret.path)) { 
    reticulate::use_condaenv(ret.path, required = TRUE)
    ok_env <- TRUE
    message("Successfully activated conda environment: r-reticulate-mvcell")
  }
}, silent = TRUE)

if (!ok_env) {
  log_step("Note: couldn't auto-select conda env for reticulate. If scATOMIC needs Python deps, set RETICULATE_PYTHON or use_condaenv manually.")
}

library(scATOMIC)

# =====================
# 1) SingleR (BlueprintEncode)
# =====================
log_step("Running SingleR with BlueprintEncodeData()...")
ref.data <- celldex::BlueprintEncodeData()

singleR.res <- SingleR(
  test   = sce,
  ref    = ref.data,
  labels = ref.data[["label.main"]],
  BPPARAM = BPPARAM
)
score.output <- singleR.res$scores
rownames(score.output) <- rownames(singleR.res)
singleR.df <- t(score.output)
log_step("SingleR finished.")

# =====================
# 2) ScType
# =====================
log_step("Running ScType with default reference DB...")
sc_db_df <- openxlsx::read.xlsx(sc_db)
# Restrict to lung cell types by default; adjust if needed
gs.list <- gene_sets_prepare(path_to_db_file = sc_db, cell_type = "Lung")

sctype.df <- sctype_score(
  scRNAseqData = assay(sce, "logcounts"),
  scaled = FALSE,
  gs  = gs.list$gs_positive,
  gs2 = gs.list$gs_negative
)
log_step("ScType finished.")

# =====================
# 3) scATOMIC
# =====================
log_step("Running scATOMIC (this may take a while)...")
count.mat <- assay(sce, "counts")
# filter zero rows with DelayedMatrix-aware function if available
if (inherits(count.mat, "DelayedMatrix") && requireNamespace("DelayedMatrixStats", quietly = TRUE)) {
  keep <- DelayedMatrixStats::rowSums2(count.mat) > 0
} else {
  keep <- rowSums(as.matrix(count.mat)) > 0
}
count.mat <- count.mat[keep, , drop = FALSE]

cell.limit <- 30000
n.core.scatomic <- 1  # keep at 1 to avoid freezing; consider exposing as a flag

if (ncol(count.mat) > cell.limit) {
  log_step("Dataset too large. Splitting into chunks of ~%d cells.", cell.limit)
  num_groups <- ceiling(ncol(count.mat) / cell.limit)
  barcodes <- colnames(count.mat)
  set.seed(1)
  barcodes <- sample(barcodes)
  group_indices <- cut(seq_along(barcodes), breaks = num_groups, labels = FALSE)
  barcode_groups <- split(barcodes, group_indices)
  
  query.barc <- barcode_groups[[1]]
  query.count <- count.mat[, query.barc]
  cell.pred <- run_scATOMIC(rna_counts = query.count, mc.cores = n.core.scatomic)
  scatomic.df <- create_summary_matrix(
    prediction_list = cell.pred, use_CNVs = FALSE, modify_results = TRUE,
    mc.cores = n.core.scatomic, raw_counts = query.count, min_prop = 0.5
  )
  
  for (i in 2:length(barcode_groups)) {
    log_step("scATOMIC chunk %d / %d", i, length(barcode_groups))
    query.barc <- barcode_groups[[i]]
    query.count <- count.mat[, query.barc]
    cell.pred <- run_scATOMIC(rna_counts = query.count, mc.cores = n.core.scatomic)
    query_results_sca <- create_summary_matrix(
      prediction_list = cell.pred, use_CNVs = FALSE, modify_results = TRUE,
      mc.cores = n.core.scatomic, raw_counts = query.count, min_prop = 0.5
    )
    scatomic.df <- rbind(scatomic.df, query_results_sca)
  }
} else {
  cell.pred <- run_scATOMIC(rna_counts = count.mat, mc.cores = n.core.scatomic)
  scatomic.df <- create_summary_matrix(
    prediction_list = cell.pred, use_CNVs = FALSE, modify_results = TRUE,
    mc.cores = n.core.scatomic, raw_counts = count.mat, min_prop = 0.5
  )
}

# =====================
# Combine results
# =====================
# SingleR
singleR.main.score <- score.to.celltype(singleR.df)
# ScType
sctype.score <- score.to.celltype(sctype.df)

# Check shared columns (cells)
if (!identical(colnames(sctype.df), colnames(singleR.df))) {
  log_step("ERROR: barcode mismatch between ScType and SingleR.")
  stop("Barcode mismatch between ScType and SingleR.")
}

celltype.list  <- list(singleR_main = singleR.main.score$top.celltype,
                       sctype       = sctype.score$top.celltype)
celltype.score <- list(singleR_main = singleR.main.score$top.score,
                       sctype       = sctype.score$top.score)

# scATOMIC: assemble prediction + max scores and align columns
scatomic.pred <- matrix(scatomic.df$scATOMIC_pred, nrow = 1)
colnames(scatomic.pred) <- rownames(scatomic.df)

scatomic.raw.s <- scatomic.df[, grep("median_score_class", colnames(scatomic.df)), drop = FALSE]
scatomic.raw.s[is.na(scatomic.raw.s)] <- 0
scatomic.max <- apply(scatomic.raw.s, 1, max)
scatomic.score <- matrix(scatomic.max, nrow = 1)
colnames(scatomic.score) <- rownames(scatomic.df)

# Force both to the same cell order
target_cols <- colnames(sctype.df)
scatomic.pred  <- scatomic.pred[,  target_cols, drop = FALSE]
scatomic.score <- scatomic.score[, target_cols, drop = FALSE]

if (!identical(colnames(sctype.df), colnames(scatomic.pred))) {
  log_step("ERROR: barcode mismatch between ScType and scATOMIC.")
  stop("Barcode mismatch between ScType and scATOMIC.")
}

celltype.list  <- append(celltype.list,  list(scATOMIC = scatomic.pred))
celltype.score <- append(celltype.score, list(scATOMIC = scatomic.score))

# =====================
# Export per-cell annotations (top1) and scores
# =====================
output.df <- as.data.frame(matrix(0, nrow = ncol(singleR.df), ncol = length(celltype.list) * 2))
annot_score_vector <- as.vector(sapply(names(celltype.list), function(x) c(paste0(x, "_annot"), paste0(x, "_score"))))
colnames(output.df) <- annot_score_vector

for (j in seq_along(celltype.list)) {
  output.df[, j * 2 - 1] <- celltype.list[[j]][1, ]
  output.df[, j * 2]     <- as.numeric(celltype.score[[j]][1, ])
}

output.df$barcode <- colnames(sctype.df)

celltype.dir <- file.path(h5ad_dir, "celltype_annotate")
if (!dir.exists(celltype.dir)) dir.create(celltype.dir)

out_csv <- file.path(celltype.dir, "MVCell_annotate_per_cell_v1.csv")
write.table(output.df, file = out_csv, quote = FALSE, row.names = FALSE, sep = ",")

log_step("Saved per-cell annotations to: %s", out_csv)
log_step("Summary: %d cells, methods = %s", ncol(sce), paste(names(celltype.list), collapse = ", "))
