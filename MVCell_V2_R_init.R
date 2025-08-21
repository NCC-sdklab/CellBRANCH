#!/usr/bin/env Rscript

# MVCell standalone run script

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(scran)
  library(scater)            # needed for logNormCounts() and sizeFactors()
  library(BiocParallel)
})

set.seed(1)


# ---- CLI args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript script.R /path/to/file.h5ad [num_cores]\n")
}
h5ad_path <- normalizePath(args[1])
if (!file.exists(h5ad_path)) stop("File not found: ", h5ad_path)

benchmark.num = if (length(args) == 3) as.integer(args[3]) else NA

# setting parameters. adding benchmark numbers for control
out_dir  <- dirname(h5ad_path)
out_base <- tools::file_path_sans_ext(basename(h5ad_path))

if (is.na(benchmark.num)) {
  out_rds  <- file.path(out_dir, paste0(out_base, "_sce.rds"))
} else {
  out_rds  <- file.path(out_dir, paste0(out_base, "_", benchmark.num, "_sce.rds"))
}

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


# ---- Parallel backend ----
if (.Platform$OS.type == "windows") {
  BPPARAM <- SnowParam(workers = ncore)
} else {
  BPPARAM <- MulticoreParam(workers = ncore)
}
BiocParallel::register(BPPARAM)

cat("H5AD file:", h5ad_path, "\n")
cat("Using", ncore, "cores (of", available_cores, "available)\n")


# ---- Read H5AD ----
# use_hdf5=TRUE keeps large matrices on disk; adjust to your environment
sce <- readH5AD(h5ad_path)


#checking first
if (file.exists(out_rds)) {
  message(out_rds, " already exists! Checking to make sure...")
  
  sce.out <- tryCatch(readRDS(out_rds), error = function(e) {
    message("Failed to read existing RDS: ", conditionMessage(e))
    NULL
  })
  
  if (!is.null(sce.out)) {
    count.check <- "counts"    %in% SummarizedExperiment::assayNames(sce.out)
    log.check   <- "logcounts" %in% SummarizedExperiment::assayNames(sce.out)
    dim.check   <- identical(dim(sce), dim(sce.out))  # or: all(dim(sce) == dim(sce.out))
    
    if (isTRUE(count.check) && isTRUE(log.check) && isTRUE(dim.check)) {
      message("All checks out. Skipping pre-processing...")
      quit(save = "no", status = 0)
    } else {
      message("Existing RDS failed checks. Rebuilding fresh from .h5ad...")
      # fall through
    }
  } else {
    message("Existing RDS unreadable. Rebuilding fresh from .h5ad...")
    # fall through
  }
} else {
  message("No RDS found. Making it fresh from .h5ad...")
  # fall through
}


# ---- Ensure 'counts' assay ----
if (!("counts" %in% assayNames(sce))) {
  message("'counts' assay not found. Checking for 'X'...")
  if ("X" %in% assayNames(sce)) {
    message("Found 'X'. Treating it as counts. NOTE: ensure 'X' is raw counts.")
    # Safer than renaming: explicitly copy to 'counts'
    assay(sce, "counts") <- assay(sce, "X")
  } else {
    stop("No 'counts' or 'X' assay found. Check the .h5ad file.")
  }
} else {
  message("'counts' assay found.")
}

# Optional: warn if counts don't look integer-ish
counts_mat <- assay(sce, "counts")
if (is.numeric(counts_mat)) {
  # handle DelayedArray/HDF5Array without full realization
  rng <- range(as.numeric(counts_mat[1:min(1e5, nrow(sce)), 1:min(1e5, ncol(sce))]))
  if (any(abs(rng - round(rng)) > 1e-8)) {
    warning("Counts assay appears non-integer; are you sure these are raw counts?")
  }
}

# ---- Ensure logcounts (library-size normalize + log1p) ----
if (!("logcounts" %in% assayNames(sce))) {
  message("'logcounts' not found. Computing size factors and log-normalized counts...")
  ncells <- ncol(sce)
  min.clust <- max(10L, min(round(ncells * 0.1), 100L))
  
  clust <- quickCluster(sce, min.size = min.clust, BPPARAM = BPPARAM)
  sce <- computeSumFactors(sce, cluster = clust, min.mean = 0.1, BPPARAM = BPPARAM)
  
  sf <- sizeFactors(sce)
  if (any(!is.finite(sf)) || any(sf <= 0)) {
    warning("Non-positive/invalid size factors detected; replacing with smallest positive value.")
    fill <- min(sf[sf > 0 & is.finite(sf)], na.rm = TRUE)
    sf[!is.finite(sf) | sf <= 0] <- fill
    sizeFactors(sce) <- sf
  }
  
  sce <- logNormCounts(sce)
  message("Generated 'logcounts'.")
} else {
  message("'logcounts' assay already present.")
}


saveRDS(sce, file = out_rds)
cat("Saved:", out_rds, "\n")

# Brief summary
cat(sprintf("Summary: %d genes × %d cells; assays = %s\n",
            nrow(sce), ncol(sce), paste(assayNames(sce), collapse = ", ")))
