#testing to see if it is more efficient
#MVCELL standalone run script

#First define script path.
all_args <- commandArgs(trailingOnly = FALSE)

# Extract the --file argument safely
file.arg <- grep("--file=", all_args, value = TRUE)

if (length(file.arg) > 0) {
  script.path <- sub("--file=", "", file.arg)
  script.dir <- dirname(normalizePath(script.path))
} else {
  # Fallback if --file= is missing (e.g., RStudio or source'd interactively)
  script.dir <- getwd()
  warning("Could not determine script location. Using current working directory instead.")
}

#second, take input
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript script.R /path/to/file.h5ad [num_cores]")
}

h5ad_path <- normalizePath(args[1])
h5ad_dir = dirname(h5ad_path)
cat("H5AD file path:", h5ad_path, "\n")

available_cores <- parallel::detectCores()
ncore <- if (length(args) >= 2) as.integer(args[2]) else NA
if (is.na(ncore) || ncore < 1) {
  ncore <- 1
} else if (ncore >= available_cores) {
  ncore <- max(1, available_cores - 1)
}
#setting multicore here
if (.Platform$OS.type == "windows") {
  BPPARAM <- BiocParallel::SnowParam(workers = ncore)
} else {
  BPPARAM <- BiocParallel::MulticoreParam(workers = ncore)
}

# #testing paths. remove later
# script.path = "/mnt/gluster_server/jhlee/Projects/python-tools/MVCell/MVCell"
# args = "/mnt/dks_nas2/data/external/Lung/CCCA/Kim_et_al_2020/Kim2020.h5ad"
# ncore = 8

library(zellkonverter)
library(SingleCellExperiment)
library(scran)

sce = readH5AD(h5ad_path)
#checking
idx.count = which(assayNames(sce) == "counts")
if (length(idx.count) == 0) {
  message("counts matrix not found. Looking for X (default)...")
  if (sum(assayNames(sce) == "X") == 1) {
    message("Default matrix found. Will treat this as a counts matrix...")
    assayNames(sce)[which(assayNames(sce) == "X")] = "counts"
  } else{
    stop("There are no counts matrices identified. Please check your .h5ad.file")
  }
} else {
  message("Count matrix found.")
}


idx.log = which(assayNames(sce) == "logcounts")
if (length(idx.log) == 0) {
  message("logcounts matrix not found. Will generate...")
  
  min.clust = min(round(dim(sce)[2] * 0.1), 100)
  clust <- quickCluster(sce, min.size = min.clust, BPPARAM = BPPARAM) 
  sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1, BPPARAM = BPPARAM)
  if (min(sce$sizeFactor) <= 0) {
    message("Size factor of 0 or less detected. This is pretty bad sample.")
    fil.val = min(sce$sizeFactor[sce$sizeFactor > 0])
    idx.neg = sce$sizeFactor <= 0
    sce$sizeFactor[idx.neg] = fil.val
  }
  
  sce <- logNormCounts(sce)
  message("logcounts generated.")
} else {
  message("logcounts found.")
}

message("saving temporary sce_temp.rds")
saveRDS(sce, file = file.path(h5ad_dir, "sce_temp.rds"))

