############################################################
# Script: 01_remove_unwanted_variation.R
# Purpose: Remove unwanted variation from ATAC-seq peak
#          counts using RUVSeq.
# Input: samplesheet.tsv, consensus peak count matrices
# Output: 1-RUVseq_output/
# Author: Bruna Toledo
############################################################

############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(RUVSeq)
  library(EDASeq)
  library(edgeR)
  library(RColorBrewer)
  library(readr)
  library(dplyr)
  library(stringr)
})

############################################################
# 2. Define input/output paths
############################################################

out_dir <- "1-RUVseq_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples_tsv <- "samplesheet.tsv"

peaks_paths <- list(
  Narrow = "NarrowPeaks/consensus_peaks.mRp.clN.featureCounts.txt",
  Broad  = "BroadPeaks/consensus_peaks.mRp.clN.featureCounts.txt"
)

out_base   <- "nPeaks"
out_prefix <- file.path(out_dir, out_base)

############################################################
# 3. Load sample metadata
############################################################

metadata <- read_tsv(samples_tsv, show_col_types = FALSE) %>%
  transmute(
    file      = as.character(file),
    sample    = as.character(sample),
    condition = factor(condition),
    subject   = as.character(subject),
    batch     = as.character(batch)
  ) %>%
  mutate(
    condition = factor(condition, levels = c("Pre-Abx", "Post-Abx")),
    subject   = factor(subject)
  )

############################################################
# 4. Define helper functions
############################################################

read_and_rename_counts <- function(path, samplesheet) {
  fc <- read_tsv(path, comment = "#", show_col_types = FALSE)
  
  stopifnot("Geneid" %in% names(fc))
  
  fc <- as.data.frame(fc)
  rownames(fc) <- fc$Geneid
  fc$Geneid <- NULL
  
  for (i in seq_len(nrow(samplesheet))) {
    old <- samplesheet$file[i]
    new <- samplesheet$sample[i]
    if (old %in% colnames(fc)) {
      colnames(fc)[colnames(fc) == old] <- new
    }
  }
  
  sample_order <- samplesheet$sample
  sample_cols_present <- intersect(sample_order, colnames(fc))
  other_cols <- setdiff(colnames(fc), sample_cols_present)
  
  fc <- fc[, c(other_cols, sample_cols_present), drop = FALSE]
  fc
}

filter_low_peaks_cpm <- function(count_matrix, min_cpm = 1, min_samples = 2) {
  message(
    "Filtering low-abundance peaks using CPM (min_cpm=",
    min_cpm, ", min_samples=", min_samples, ")"
  )
  
  cpm_mat <- edgeR::cpm(count_matrix)
  keep <- rowSums(cpm_mat > min_cpm) >= min_samples
  
  message(sum(keep), " peaks kept out of ", nrow(count_matrix))
  
  count_matrix[keep, , drop = FALSE]
}

run_ruvr_pipeline <- function(filtered_counts, samplesheet, out_dir, out_base, k = 1) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_prefix <- file.path(out_dir, out_base)
  
  pd <- samplesheet %>%
    filter(sample %in% colnames(filtered_counts)) %>%
    arrange(match(sample, colnames(filtered_counts)))
  
  stopifnot(identical(pd$sample, colnames(filtered_counts)))
  
  mat <- as.matrix(filtered_counts)
  storage.mode(mat) <- "integer"
  
  set_raw <- newSeqExpressionSet(
    counts = mat,
    phenoData = data.frame(
      row.names = pd$sample,
      condition = pd$condition,
      subject   = pd$subject,
      batch     = pd$batch,
      sample    = pd$sample
    )
  )
  
  n_cond   <- nlevels(pData(set_raw)$condition)
  cols_all <- brewer.pal(max(2, min(8, n_cond)), "Set2")[seq_len(n_cond)]
  cols     <- cols_all[pData(set_raw)$condition]
  
  png(paste0(out_prefix, "-1cpm_PCA_raw.png"), width = 1200, height = 1000, res = 150)
  plotPCA(set_raw, col = cols, cex = 1.2)
  dev.off()
  
  set_uq <- betweenLaneNormalization(set_raw, which = "upper")
  
  png(paste0(out_prefix, "-1cpm_PCA_uq.png"), width = 1200, height = 1000, res = 150)
  plotPCA(set_uq, col = cols, cex = 1.2)
  dev.off()
  
  design_ruvr <- model.matrix(~ condition, data = pData(set_uq))
  
  y <- DGEList(counts = EDASeq::counts(set_uq), group = pData(set_uq)$condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateDisp(y, design_ruvr)
  
  fit <- glmFit(y, design_ruvr)
  resids <- residuals(fit, type = "deviance")
  
  if (k > 0) {
    c_idx <- rownames(EDASeq::counts(set_uq))
    set_ruvr <- RUVr(set_uq, c_idx, k = k, res = resids)
  } else {
    set_ruvr <- set_uq
  }
  
  png(
    paste0(out_prefix, sprintf("-1cpm_PCA_ruvr_k%d.png", k)),
    width = 1200, height = 1000, res = 150
  )
  plotPCA(set_ruvr, col = cols, cex = 1.2)
  dev.off()
  
  norm_counts <- set_ruvr@assayData[["normalizedCounts"]]
  if (is.null(norm_counts)) norm_counts <- EDASeq::counts(set_ruvr)
  
  norm_dt <- data.frame(peak_id = rownames(norm_counts), norm_counts, check.names = FALSE)
  out_tsv <- paste0(out_prefix, sprintf("-1cpm_RUV_k%d.tsv", k))
  write.table(norm_dt, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  
  obj_name   <- sprintf("set_ruvr_%s_k%d", out_base, k)
  rdata_file <- file.path(out_dir, sprintf("%s_set_ruvr_k%d.RData", out_base, k))
  
  assign(obj_name, set_ruvr, envir = .GlobalEnv)
  save(list = obj_name, file = rdata_file)
  
  message("Saved RUVSeq object as '", obj_name, "' -> ", rdata_file)
  
  load(rdata_file, envir = .GlobalEnv)
  message("Loaded RUVSeq object into workspace as '", obj_name, "'")
  
  invisible(list(
    set_raw  = set_raw,
    set_uq   = set_uq,
    set_ruvr = set_ruvr,
    obj_name = obj_name,
    outfile  = out_tsv,
    rdata    = rdata_file
  ))
}

############################################################
# 5. Load consensus peak count matrices
############################################################

n_peaks <- read_and_rename_counts(peaks_paths$Narrow, metadata)
b_peaks <- read_and_rename_counts(peaks_paths$Broad, metadata)

n_peaks$peak_id <- rownames(n_peaks)

############################################################
# 6. Extract peak count matrices
############################################################

np_counts <- as.matrix(n_peaks[, metadata$sample, drop = FALSE])
bp_counts <- as.matrix(b_peaks[, metadata$sample, drop = FALSE])

############################################################
# 7. Filter low-count peaks
############################################################

lib_sizes_np <- colSums(np_counts)
cv_lib_np <- sd(lib_sizes_np) / mean(lib_sizes_np)

if (cv_lib_np > 0.3) {
  fd_np_counts <- filter_low_peaks_cpm(np_counts, min_cpm = 1, min_samples = 2)
} else {
  keep <- rowSums(np_counts > 25) >= 3
  fd_np_counts <- np_counts[keep, , drop = FALSE]
  message(sum(keep), " peaks kept after raw count filter.")
}

lib_sizes_bp <- colSums(bp_counts)
cv_lib_bp <- sd(lib_sizes_bp) / mean(lib_sizes_bp)

if (cv_lib_bp > 0.3) {
  fd_bp_counts <- filter_low_peaks_cpm(bp_counts, min_cpm = 1, min_samples = 2)
} else {
  keep <- rowSums(bp_counts > 25) >= 3
  fd_bp_counts <- bp_counts[keep, , drop = FALSE]
  message(sum(keep), " peaks kept after raw count filter.")
}

############################################################
# 8. Run RUVSeq normalization
############################################################

run_ruvr_pipeline(fd_np_counts, metadata, out_dir, out_base = "nPeaks", k = 3)

############################################################
# 9. Save workspace
############################################################

save.image(file = "RUVseq_workspace.RData")



############################################################
# End of script
############################################################