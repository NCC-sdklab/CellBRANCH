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

#testing paths. remove later
script.path = "/mnt/gluster_server/jhlee/Projects/python-tools/MVCell/MVCell"
args = "/mnt/dks_nas2/data/external/Lung/CCCA/Kim_et_al_2020/Kim2020.h5ad"
ncore = 8





#score function
#porting nested functions
################################################################################
score.to.celltype = function(celltype.score, top.rank = 3){
  query.celltype = matrix(data = NA, ncol = ncol(celltype.score), nrow = top.rank)
  colnames(query.celltype) = colnames(celltype.score)
  
  query.celltype.score = matrix(data = 0, ncol = ncol(celltype.score), nrow = top.rank)
  colnames(query.celltype.score) = colnames(celltype.score)
  
  for (k in 1:ncol(query.celltype)) {
    query.table = sort(celltype.score[,k], decreasing = T)
    query.celltype[,k] = names(query.table)[1:top.rank]
    query.celltype.score[,k] = query.table[1:top.rank]
  }
  
  output.list = list(top.celltype = query.celltype, top.score = query.celltype.score)
  
  return(output.list)
} #consider calculating scores for this

#deprecated for now
score.normalize = function(score.df, log.norm = F){
  score.arr = as.vector(as.matrix(score.df))
  
  if (log.norm) {
    score.arr = score.arr - min(score.arr) + 0.01 #just to above 0
    score.arr = log(score.arr)
  }
  
  score.mean = mean(score.arr)
  score.sd = sd(score.arr)
  
  score.z = (score.arr - score.mean) / score.sd
  
  #capping
  score.z[score.z > 3] = 3
  score.z[score.z < -3] = -3
  
  score.zn = score.z / 3
  score.zmat <- matrix(score.zn, nrow = nrow(score.df), ncol = ncol(score.df))
  
  colnames(score.zmat) = colnames(score.df)
  rownames(score.zmat) = rownames(score.df)
  
  return(score.zmat)
}

#currently using. simpler to deal with
score.rank = function(score.df, log.n = F){
  
  if (log.n) {
    score.df = log(score.df - min(score.df) + 0.01)
  }
  
  # #normalization per cell types. this is going to be work
  # score.rmat = (apply(as.matrix(score.df), 2, rank)) 
  # score.rmat = 1 - score.rmat / nrow(score.rmat) #ranking and flipping rank so the high rank (1) becomes high value when normalized
  
  score.arr = as.vector(as.matrix(score.df))
  score.scale = pnorm(scale(score.arr)) #precentile using CDF
  
  score.rmat = matrix(score.scale, nrow = nrow(score.df), ncol = ncol(score.df))
  
  # score.rank = rank(score.arr)
  # score.rank = score.rank / max(score.rank)
  # 
  # score.rmat = matrix(score.rank, nrow = nrow(score.df), ncol = ncol(score.df))
  # 
  colnames(score.rmat) = colnames(score.df)
  rownames(score.rmat) = rownames(score.df)
  
  return(score.rmat)
}
################################################################################


#req min lib
library(zellkonverter)
library(SingleCellExperiment)
library(scran)

#for parallel
library(BiocParallel)

#singleR
source(file.path(script.path, "sctype/sctype_score_.R"))
source(file.path(script.path, "sctype/gene_sets_prepare.R"))

library(SingleR)
library(celldex)

#sc-type
library(openxlsx)
library(HGNChelper)

#scATOMIC
library(reticulate)
#extracting and specying path
cb <- conda_binary()
cb.base = sub("/bin/(conda|mamba)$", "", cb)
ret.path = file.path(cb.base, "envs", "r-reticulate")
use_condaenv(ret.path, required = TRUE)
reticulate::import("magic")
library(scATOMIC)


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
  message("logcounts found. Proceeding...")
}

message("Running singleR with Blueprint Encode Data.")

#ref.data = HumanPrimaryCellAtlasData() #this for second model
ref.data = BlueprintEncodeData() #this is the trained.

singleR.res = SingleR(test = sce,
                      ref = ref.data,
                      labels = ref.data[["label.main"]],
                      BPPARAM = BPPARAM)
score.output = singleR.res$scores
rownames(score.output) = rownames(singleR.res)
singleR.df = t(score.output)
message("SingleR finished.")

message("Running sctype with default reference file.")
db.full = file.path(script.path, "sctype/ScTypeDB_full.xlsx")
db.df = read.xlsx(db.full)
gs.list = gene_sets_prepare(path_to_db_file = db.full, cell_type = "Lung")
sctype.df = sctype_score(scRNAseqData = assay(sce, "logcounts"), 
                        scaled = FALSE, 
                        gs = gs.list$gs_positive, 
                        gs2 = gs.list$gs_negative)


message("Running sctype_score. This will take a while...")
count.mat = assay(sce, "counts")
count.mat = count.mat[rowSums(count.mat) > 0,]

cell.limit = 30000
n.core = 1 #fixing it to 1 because it keeps freezing.
if (ncol(count.mat) > cell.limit) { #might have to lower it
  message("dataset too large. separating to keep cells under ", cell.limit)
  
  num_groups <- ceiling(ncol(count.mat) / cell.limit)
  barcodes <- sample(colnames(count.mat)) #suffling
  
  # Split the barcodes into 'num_groups' nearly equal groups.
  group_indices <- cut(seq_along(barcodes), breaks = num_groups, labels = FALSE)
  barcode_groups <- split(barcodes, group_indices)
  
  message("Divided barcodes into ", length(barcode_groups), ". ")
  
  query.barc = barcode_groups[[1]]
  query.count = count.mat[,query.barc]
  
  cell.pred = run_scATOMIC(rna_counts = query.count, mc.cores = n.core)
  scatomic.df <- create_summary_matrix(prediction_list = cell.pred, use_CNVs = F, 
                                       modify_results = T, mc.cores = n.core, 
                                       raw_counts = query.count, min_prop = 0.5 )
  
  for (i in 2:length(barcode_groups)) {
    message("Performing chunk ", i, " out of ", length(barcode_groups))
    query.barc = barcode_groups[[i]]
    query.count = count.mat[,query.barc]
    
    cell.pred = run_scATOMIC(rna_counts = query.count, mc.cores = n.core)
    query_results_sca <- create_summary_matrix(prediction_list = cell.pred, use_CNVs = F, 
                                               modify_results = T, mc.cores = n.core, 
                                               raw_counts = query.count, min_prop = 0.5 )
    
    scatomic.df = rbind(scatomic.df, query_results_sca)
  }
  
} else {
  cell.pred = run_scATOMIC(rna_counts = count.mat, mc.cores = n.core)
  scatomic.df <- create_summary_matrix(prediction_list = cell.pred, use_CNVs = F, 
                                       modify_results = T, mc.cores = n.core, 
                                       raw_counts = count.mat, min_prop = 0.5 )
}


#combining
#checking, aggregating
singleR.main.score = score.to.celltype(singleR.df)
sctype.score = score.to.celltype(sctype.df)

#basic version
if (identical(colnames(sctype.df), colnames(singleR.df))) {
  celltype.list = list(singleR_main = singleR.main.score$top.celltype, sctype = sctype.score$top.celltype)
  celltype.score = list(singleR_main = singleR.main.score$top.score, sctype = sctype.score$top.score)
} else {
  log_step("barcode mismatch between sctype and singleR.")
  stop()
}

scatomic.raw.s = scatomic.df[,grep("median_score_class", colnames(scatomic.df))]
#NA score will be treated as 0
scatomic.raw.s[is.na(scatomic.raw.s)] = 0 #generally this applies to Any Cell, which means failed prediction
scatomic.score = matrix(apply(scatomic.raw.s, 1, max), nrow = 1)
scatomic.pred = matrix(data = scatomic.df$scATOMIC_pred, nrow = 1)
colnames(scatomic.pred) = rownames(scatomic.df)

#forcing match
scatomic.pred = scatomic.pred[,colnames(sctype.df), drop = F]

if (identical(colnames(sctype.df), colnames(scatomic.pred))) {
  celltype.list = append(celltype.list, list(scATOMIC = scatomic.pred))
  celltype.score = append(celltype.score, list(scATOMIC = scatomic.score))
} else {
  log_step("barcode mismatch between sctype and scATOMIC")
  stop()
}

#exporting the cell-level data
output.df = matrix(data = 0, nrow = ncol(singleR.df), ncol = length(celltype.list) *2)
output.df = as.data.frame(output.df)
annot_score_vector <- as.vector(sapply(names(celltype.list), function(x) c(paste0(x, "_annot"), paste0(x, "_score"))))
colnames(output.df) = annot_score_vector

for (j in 1:length(celltype.list)) {
  output.df[,j*2-1] = celltype.list[[j]][1, ]
  output.df[,j*2] = as.numeric(celltype.score[[j]][1, ])
}

output.df$barcode = colnames(sctype.df)
celltype.dir = file.path(h5ad_dir, "celltype_annotate")
if (!dir.exists(celltype.dir)) {
  dir.create(celltype.dir)
}
write.table(output.df, file = file.path(celltype.dir, "MVCell_annotate_per_cell_extended.csv"), quote = F, row.names = F, sep = ",")
