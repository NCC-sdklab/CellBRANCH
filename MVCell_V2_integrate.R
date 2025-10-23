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
# all_args <- commandArgs(trailingOnly = FALSE)
# file.arg <- grep("--file=", all_args, value = TRUE)
# script.dir <- if (length(file.arg) > 0) dirname(normalizePath(sub("--file=", "", file.arg))) else getwd()


#testing args
args = "/mnt/"

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

out_performance = file.path(out_dir, "performance_report")
if (!dir.exists(out_performance)) {
  dir.create(out_performance)
}

# ----- Load MVCell results -----

if (!file.exists(out_mvcell)) stop("Missing MVCell results: ", out_mvcell)
mvcell.df = read.csv(out_mvcell)

if (!file.exists(out_rds)) stop("Missing preprocessed RDS: ", out_rds)
sce <- readRDS(out_rds)

meta = as.data.frame(colData(sce))
if (!"barcode" %in% colnames(meta)) meta$barcode <- rownames(meta)

if (identical(meta$barcode, mvcell.df$barcode)) {
  message("barcodes consistent. Proceeding...")
} else {
  stop("barcode mismatch")
}

meta.m = merge(meta, mvcell.df, all.x = T, by = "barcode", sorf = F)
rownames(meta.m) = meta.m$barcode
meta.m = meta.m[colnames(sce),]

test.data = as.data.frame(data.table::fread(out_raw_annlevel))
meta.m = merge(meta.m, test.data, all.x = T, by = "barcode")
rownames(meta.m) = meta.m$barcode



# -------- Multiclass MCC (no extra packages) --------
mcc_multiclass <- function(truth, pred) {
  # truth and pred should be factors on the same level set
  truth <- factor(truth)
  pred  <- factor(pred, levels = levels(truth))
  M <- table(truth, pred)              # K x K
  
  t_k <- rowSums(M)                    # true per class
  p_k <- colSums(M)                    # pred per class
  c   <- sum(diag(M))                  # trace
  s   <- sum(M)                        # total
  
  num <- (c * s) - sum(t_k * p_k)
  den <- sqrt( (s^2 - sum(p_k^2)) * (s^2 - sum(t_k^2)) )
  
  if (den == 0) return(NA_real_)
  num / den
}

# -------- Per-class F1 and Macro-F1 (no MLmetrics) --------
per_class_f1 <- function(truth, pred) {
  lv <- levels(truth)
  f1 <- setNames(numeric(length(lv)), lv)
  
  for (k in lv) {
    tp <- sum(truth == k & pred == k, na.rm = TRUE)
    fp <- sum(truth != k & pred == k, na.rm = TRUE)
    fn <- sum(truth == k & pred != k, na.rm = TRUE)
    
    if ((tp + fp) == 0 || (tp + fn) == 0) {
      f1[k] <- NA_real_
    } else {
      precision <- tp / (tp + fp)
      recall    <- tp / (tp + fn)
      if ((precision + recall) == 0) {
        f1[k] <- NA_real_
      } else {
        f1[k] <- 2 * precision * recall / (precision + recall)
      }
    }
  }
  f1
}

# -------- Assemble results for your methods --------
library(ggplot2)
library(dplyr)
library(tidyr)

methods <- c("ann_level_1_prediction", "celltypist_ann1", "scimilarity_ann1",
             "scATOMIC_ann1", "singleR_ann1","sctype_ann1")

truth <- factor(meta.m$ann_level_1)

methods.name = c("MVCell", "celltypist", "SCimiliary", "scATOMIC", "singleR", "scType")


score_rows <- lapply(methods, function(m) {
  pred <- factor(meta.m[[m]], levels = levels(truth))
  
  # per-class F1 and macro-F1
  f1_k <- per_class_f1(truth, pred)
  f1_macro <- mean(f1_k, na.rm = TRUE)
  
  # overall MCC (multiclass)
  mcc_overall <- mcc_multiclass(truth, pred)
  
  # per-class MCC (one-vs-rest)
  lv <- levels(truth)
  mcc_k <- sapply(lv, function(ct) {
    bin_truth <- factor(ifelse(truth == ct, ct, "other"), levels = c(ct, "other"))
    bin_pred  <- factor(ifelse(pred  == ct, ct, "other"), levels = c(ct, "other"))
    mcc_multiclass(bin_truth, bin_pred)  # works for binary (2x2) as well
  }, USE.NAMES = TRUE)
  
  rbind(
    data.frame(method = m, metric = "F1_macro",  celltype = NA,       score = f1_macro),
    data.frame(method = m, metric = "MCC",       celltype = NA,       score = mcc_overall),
    data.frame(method = m, metric = "F1",        celltype = names(f1_k), score = unname(f1_k)),
    data.frame(method = m, metric = "MCC_cell",  celltype = names(mcc_k), score = unname(mcc_k))
  )
}) %>% dplyr::bind_rows()

results_ann1 <- bind_rows(score_rows)
results_ann1$algorithm <- methods.name[match(results_ann1$method, methods)]
results_ann1$algorithm = factor(results_ann1$algorithm, levels = methods.name)

write.csv(results_ann1, file = file.path(out_performance, "performance_report_ann1_all.csv"))



# -------- Overall (F1_macro + MCC) with stars --------
overall <- results_ann1 %>% 
  filter(metric %in% c("F1_macro", "MCC")) %>%
  group_by(metric) %>%
  mutate(best = score == max(score, na.rm = TRUE)) %>%
  ungroup()

overall.mcc.f1 = ggplot(overall, aes(x = algorithm, y = score, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    data = overall %>% filter(best),
    aes(label = "★", y = score),
    position = position_dodge(width = 0.7),
    vjust = -0.3, size = 6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Overall performance", x = "Method", y = "Score") +
  theme_minimal()

ggsave(filename = "Level1_performance_overall.png", plot = overall.mcc.f1, path = out_performance,
       width = 7, height = 5, bg = "white")


# build per-class F1 (you already have `results`)
perclass <- results_ann1 %>%
  filter(metric == "F1") %>%
  mutate(method = factor(method, levels = methods),
         celltype = factor(celltype, levels = unique(celltype))) %>%
  group_by(celltype) %>%
  mutate(best = score == max(score, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(celltype, method) %>%            # keep a stable order
  mutate(.rowid = row_number())

pos <- position_dodge(width = 0.7)

# 1) draw bars (but don't print yet)
p0 <- ggplot(perclass, aes(x = celltype, y = score, fill = algorithm)) +
  geom_col(position = pos, width = 0.6)

# 2) extract exact x positions of bars
gb <- ggplot_build(p0)
barpos <- gb$data[[1]][, c("x","y","PANEL")]
perclass_pos <- perclass
perclass_pos$x_star <- barpos$x
perclass_pos$y_star <- barpos$y + diff(range(barpos$y, na.rm = TRUE))*0.04  # a little above bar

# 3) add stars at those x positions
p_f1 <- p0 +
  geom_text(
    data = perclass_pos %>% filter(best),
    aes(x = x_star, y = y_star, label = "★"),
    inherit.aes = FALSE, size = 6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "F1 per Level‑1 cell type", x = "Cell type", y = "F1") +
  theme_minimal()

ggsave(filename = "Level1_performance_celltype_F1.png", plot = p_f1, path = out_performance,
       width = 7, height = 5, bg = "white")


#MCC
perclass_mcc <- results_ann1 %>%
  filter(metric == "MCC_cell") %>%
  mutate(method = factor(method, levels = methods),
         celltype = factor(celltype, levels = unique(celltype))) %>%
  group_by(celltype) %>%
  mutate(best = score == max(score, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(celltype, method) %>%            # keep a stable order
  mutate(.rowid = row_number())

perclass_mcc <- perclass_mcc %>%
  mutate(method = factor(method, levels = methods),
         celltype = factor(celltype, levels = unique(celltype))) %>%
  group_by(celltype) %>%
  mutate(best = score == max(score, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(celltype, method) %>%
  mutate(.rowid = row_number())

pos2 <- position_dodge(width = 0.7)

p0_mcc <- ggplot(perclass_mcc, aes(x = celltype, y = score, fill = algorithm)) +
  geom_col(position = pos2, width = 0.6)

gb2 <- ggplot_build(p0_mcc)
barpos2 <- gb2$data[[1]][, c("x","y","PANEL")]
perclass_mcc$x_star <- barpos2$x
perclass_mcc$y_star <- barpos2$y + diff(range(barpos2$y, na.rm = TRUE))*0.04

p_mcc <- p0_mcc +
  geom_text(
    data = perclass_mcc %>% filter(best),
    aes(x = x_star, y = y_star, label = "★"),
    inherit.aes = FALSE, size = 6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "MCC per Level‑1 cell type", x = "Cell type", y = "MCC") +
  theme_minimal()

ggsave(filename = "Level1_performance_celltype_MCC.png", plot = p_mcc, path = out_performance,
       width = 7, height = 5, bg = "white")


#level 2

methods <- c("ann_level_2_prediction", "celltypist_ann2", "scimilarity_ann2",
             "scATOMIC_ann2", "singleR_ann2",
             "sctype_ann2")
truth <- factor(meta.m$ann_level_2)


score_rows <- lapply(methods, function(m) {
  pred <- factor(meta.m[[m]], levels = levels(truth))
  
  # per-class F1 and macro-F1
  f1_k <- per_class_f1(truth, pred)
  f1_macro <- mean(f1_k, na.rm = TRUE)
  
  # overall MCC (multiclass)
  mcc_overall <- mcc_multiclass(truth, pred)
  
  # per-class MCC (one-vs-rest)
  lv <- levels(truth)
  mcc_k <- sapply(lv, function(ct) {
    bin_truth <- factor(ifelse(truth == ct, ct, "other"), levels = c(ct, "other"))
    bin_pred  <- factor(ifelse(pred  == ct, ct, "other"), levels = c(ct, "other"))
    mcc_multiclass(bin_truth, bin_pred)  # works for binary (2x2) as well
  }, USE.NAMES = TRUE)
  
  rbind(
    data.frame(method = m, metric = "F1_macro",  celltype = NA,       score = f1_macro),
    data.frame(method = m, metric = "MCC",       celltype = NA,       score = mcc_overall),
    data.frame(method = m, metric = "F1",        celltype = names(f1_k), score = unname(f1_k)),
    data.frame(method = m, metric = "MCC_cell",  celltype = names(mcc_k), score = unname(mcc_k))
  )
}) %>% dplyr::bind_rows()

results_ann2 <- bind_rows(score_rows)
results_ann2$algorithm <- methods.name[match(results_ann2$method, methods)]
results_ann2$algorithm = factor(results_ann2$algorithm, levels = methods.name)

write.csv(results_ann2, file = file.path(out_performance, "performance_report_ann2_all.csv"))


# -------- Overall (F1_macro + MCC) with stars --------
overall <- results_ann2 %>% 
  filter(metric %in% c("F1_macro", "MCC")) %>%
  group_by(metric) %>%
  mutate(best = score == max(score, na.rm = TRUE)) %>%
  ungroup()

overall.mcc.f1 = ggplot(overall, aes(x = algorithm, y = score, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    data = overall %>% filter(best),
    aes(label = "★", y = score),
    position = position_dodge(width = 0.7),
    vjust = -0.3, size = 6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Overall performance: Level 2", x = "Method", y = "Score") +
  theme_minimal()

ggsave(filename = "Level2_performance_overall.png", plot = overall.mcc.f1, path = out_performance,
       width = 7, height = 5, bg = "white")


# build per-class F1 (you already have `results`)
perclass <- results_ann2 %>%
  filter(metric == "F1") %>%
  mutate(method = factor(method, levels = methods),
         celltype = factor(celltype, levels = unique(celltype))) %>%
  group_by(celltype) %>%
  mutate(best = score == max(score, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(celltype, method) %>%            # keep a stable order
  mutate(.rowid = row_number())

pos <- position_dodge(width = 0.7)

# 1) draw bars (but don't print yet)
p0 <- ggplot(perclass, aes(x = celltype, y = score, fill = algorithm)) +
  geom_col(position = pos, width = 0.6)

# 2) extract exact x positions of bars
gb <- ggplot_build(p0)
barpos <- gb$data[[1]][, c("x","y","PANEL")]
perclass_pos <- perclass
perclass_pos$x_star <- barpos$x
perclass_pos$y_star <- barpos$y + diff(range(barpos$y, na.rm = TRUE))*0.04  # a little above bar

# 3) add stars at those x positions
p_f1 <- p0 +
  geom_text(
    data = perclass_pos %>% filter(best),
    aes(x = x_star, y = y_star, label = "★"),
    inherit.aes = FALSE, size = 6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "F1 per Level-2 cell type", x = "Cell type", y = "F1") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggsave(filename = "Level2_performance_celltype_F1.png", plot = p_f1, path = out_performance,
       width = 7, height = 5, bg = "white")


#MCC
perclass_mcc <- results_ann2 %>%
  filter(metric == "MCC_cell") %>%
  mutate(method = factor(method, levels = methods),
         celltype = factor(celltype, levels = unique(celltype))) %>%
  group_by(celltype) %>%
  mutate(best = score == max(score, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(celltype, method) %>%            # keep a stable order
  mutate(.rowid = row_number())

perclass_mcc <- perclass_mcc %>%
  mutate(method = factor(method, levels = methods),
         celltype = factor(celltype, levels = unique(celltype))) %>%
  group_by(celltype) %>%
  mutate(best = score == max(score, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(celltype, method) %>%
  mutate(.rowid = row_number())

pos2 <- position_dodge(width = 0.7)

p0_mcc <- ggplot(perclass_mcc, aes(x = celltype, y = score, fill = algorithm)) +
  geom_col(position = pos2, width = 0.6)

gb2 <- ggplot_build(p0_mcc)
barpos2 <- gb2$data[[1]][, c("x","y","PANEL")]
perclass_mcc$x_star <- barpos2$x
perclass_mcc$y_star <- barpos2$y + diff(range(barpos2$y, na.rm = TRUE))*0.04

p_mcc <- p0_mcc +
  geom_text(
    data = perclass_mcc %>% filter(best),
    aes(x = x_star, y = y_star, label = "★"),
    inherit.aes = FALSE, size = 6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "MCC per Level‑2 cell type", x = "Cell type", y = "MCC") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )


ggsave(filename = "Level2_performance_celltype_MCC.png", plot = p_mcc, path = out_performance,
       width = 7, height = 5, bg = "white")


