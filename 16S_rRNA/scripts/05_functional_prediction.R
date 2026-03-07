############################################################
# Script: 05_functional_prediction.R
# Purpose: Analyze PICRUSt2-predicted functional profiles
#          and identify differentially abundant MetaCyc
#          pathways and KEGG ortholog/pathway features
#          between Pre-Abx and Post-Abx samples.
# Input: PICRUSt2 output tables and curated sample metadata
# Output: annotated ALDEx2 result tables, summary plots,
#         and session information
# Author: Bruna Toledo
############################################################

## =========================================================
## 1. Load libraries
## =========================================================

suppressPackageStartupMessages({
  library(ALDEx2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(openxlsx)
  library(readr)
  library(ggpicrust2)
})

## =========================================================
## 2. Define input and output paths
## =========================================================

metadata_file <- "16S_rRNA/results/intermediate/mpse_taxonomy.rds"

picrust_dir <- "16S_rRNA/pipeline/Nextflow/picrust/all_output"
metacyc_annotation_file <- "16S_rRNA/pipeline/Nextflow/picrust/METACYC_path_abun_unstrat_descrip.tsv"

out_dir_tables  <- "16S_rRNA/results/tables/functional_prediction"
out_dir_figures <- "16S_rRNA/results/figures/functional_prediction"
out_dir_intermediate <- "16S_rRNA/results/intermediate"
out_dir_env <- "environment"

if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)
if (!dir.exists(out_dir_figures)) dir.create(out_dir_figures, recursive = TRUE)
if (!dir.exists(out_dir_intermediate)) dir.create(out_dir_intermediate, recursive = TRUE)
if (!dir.exists(out_dir_env)) dir.create(out_dir_env, recursive = TRUE)

## =========================================================
## 3. Load metadata from MPSE object
## =========================================================

mpse <- readRDS(metadata_file)
sample_metadata <- mp_extract_sample(mpse) %>% as.data.frame()

sample_metadata$condition <- factor(
  sample_metadata$condition,
  levels = c("Pre-Abx", "Post-Abx")
)

metadata <- sample_metadata[, c("Sample", "subject", "condition")]
rownames(metadata) <- metadata$Sample
metadata_samples <- rownames(metadata)

message("Sample metadata loaded successfully.")

## =========================================================
## 4. Define sample renaming map
## =========================================================

sample_labels <- c(
  "M061"  = "S6_T1",  "M071"  = "S7_T1",  "M081"  = "S8_T1",  "M091"  = "S9_T1",
  "M101"  = "S10_T1", "M111"  = "S11_T1", "M121"  = "S12_T1", "M131"  = "S13_T1",
  "S15T1" = "S15_T1", "S16T1" = "S16_T1", "S2T1"  = "S2_T1",  "S3T1"  = "S3_T2",
  "S4T1"  = "S4_T1",  "S5T1"  = "S5_T1",
  "M062"  = "S6_T2",  "M072"  = "S7_T2",  "M082"  = "S8_T2",  "M092"  = "S9_T2",
  "M102"  = "S10_T2", "M112"  = "S11_T2", "M122"  = "S12_T2", "M132"  = "S13_T2",
  "S15T2" = "S15_T2", "S16T2" = "S16_T2", "S2T2"  = "S2_T2",  "S3T2"  = "S3_T1",
  "S4T2"  = "S4_T2",  "S5T2"  = "S5_T2"
)

## =========================================================
## 5. Load PICRUSt2 output tables
## =========================================================

kegg_data <- ko2kegg_abundance(
  file.path(picrust_dir, "KO_metagenome_out", "pred_metagenome_unstrat.tsv")
)

metacyc_data <- read.delim(
  file.path(picrust_dir, "pathways_out", "path_abun_unstrat.tsv"),
  row.names = 1,
  check.names = FALSE
)

message("PICRUSt2 output tables loaded successfully.")

## =========================================================
## 6. Rename sample columns and align to metadata
## =========================================================

rename_and_subset_columns <- function(df, sample_map, target_samples) {
  new_names <- sample_map[colnames(df)]
  
  if (any(is.na(new_names))) {
    missing_cols <- colnames(df)[is.na(new_names)]
    stop("Missing sample label mapping for columns: ", paste(missing_cols, collapse = ", "))
  }
  
  colnames(df) <- unname(new_names)
  
  missing_target <- setdiff(target_samples, colnames(df))
  if (length(missing_target) > 0) {
    stop("The following metadata samples are missing from abundance table: ",
         paste(missing_target, collapse = ", "))
  }
  
  df[, target_samples, drop = FALSE]
}

kegg_abundance <- rename_and_subset_columns(kegg_data, sample_labels, metadata_samples)
metacyc_abundance <- rename_and_subset_columns(metacyc_data, sample_labels, metadata_samples)

## =========================================================
## 7. Convert abundance tables to integer matrices
## =========================================================

kegg_abundance_int <- round(kegg_abundance)
stopifnot(all(kegg_abundance_int == floor(kegg_abundance_int)))

metacyc_abundance_int <- round(metacyc_abundance)
stopifnot(all(metacyc_abundance_int == floor(metacyc_abundance_int)))

conds <- factor(metadata$condition, levels = c("Pre-Abx", "Post-Abx"))

message("Functional abundance tables aligned and converted to integer matrices.")

## =========================================================
## 8. Define helper functions
## =========================================================

run_aldex <- function(count_matrix, conditions) {
  aldex(
    reads = count_matrix,
    conditions = conditions,
    mc.samples = 128,
    test = "t",
    paired.test = TRUE,
    effect = TRUE,
    denom = "all"
  )
}

prepare_aldex_df <- function(aldex_res, annotation_df = NULL, annotation_col = NULL) {
  df <- as.data.frame(aldex_res)
  df$feature <- rownames(df)
  
  if (!is.null(annotation_df) && !is.null(annotation_col)) {
    df$description <- annotation_df[rownames(df), annotation_col]
  }
  
  df <- df %>%
    arrange(desc(effect)) %>%
    mutate(condition_enriched = ifelse(effect > 0, "Post-Abx", "Pre-Abx"))
  
  df
}

filter_significant_aldex <- function(df, p_thresh = 0.05) {
  df %>%
    filter(wi.eBH < p_thresh)
}

make_effect_plot <- function(df, feature_col, fill_col = "condition_enriched",
                             xlab = "", title = NULL, top_n = NULL) {
  plot_df <- df
  
  if (!is.null(top_n)) {
    plot_df <- plot_df %>%
      slice_max(order_by = abs(effect), n = top_n)
  }
  
  plot_df[[fill_col]] <- factor(
    plot_df[[fill_col]],
    levels = c("Pre-Abx", "Post-Abx")
  )
  
  plot_df[[feature_col]] <- reorder(plot_df[[feature_col]], plot_df$effect)
  
  ggplot(plot_df, aes_string(x = feature_col, y = "effect", fill = fill_col)) +
    geom_col(color = "black", width = 0.6) +
    coord_flip() +
    scale_fill_manual(values = c("Pre-Abx" = "#3C548866", "Post-Abx" = "#87255B66")) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.major = element_line(color = "grey96", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey96", linewidth = 0.2)
    ) +
    labs(
      title = title,
      x = xlab,
      y = "Effect size"
    )
}

## =========================================================
## 9. Export raw abundance tables
## =========================================================

kegg_export <- kegg_abundance %>%
  as.data.frame() %>%
  rownames_to_column("feature")

metacyc_export <- metacyc_abundance %>%
  as.data.frame() %>%
  rownames_to_column("feature")

write.xlsx(
  kegg_export,
  file.path(out_dir_tables, "PICRUSt2_KEGG_abundance.xlsx"),
  overwrite = TRUE
)

write.xlsx(
  metacyc_export,
  file.path(out_dir_tables, "PICRUSt2_MetaCyc_abundance.xlsx"),
  overwrite = TRUE
)

message("Raw functional abundance tables exported successfully.")

## =========================================================
## 10. MetaCyc differential abundance analysis
## =========================================================

message("Running ALDEx2 for MetaCyc pathways...")
aldex_picrust2_metacyc <- run_aldex(metacyc_abundance_int, conds)

metacyc_meta <- read.delim(
  metacyc_annotation_file,
  row.names = 1,
  check.names = FALSE
)

aldex_picrust2_metacyc_all <- prepare_aldex_df(
  aldex_picrust2_metacyc,
  annotation_df = metacyc_meta,
  annotation_col = "description"
)

aldex_picrust2_metacyc_sig <- filter_significant_aldex(aldex_picrust2_metacyc_all)

write.xlsx(
  aldex_picrust2_metacyc_all,
  file.path(out_dir_tables, "PICRUSt2_MetaCyc_ALDEx2_all.xlsx"),
  overwrite = TRUE
)

write.xlsx(
  aldex_picrust2_metacyc_sig,
  file.path(out_dir_tables, "PICRUSt2_MetaCyc_ALDEx2_significant.xlsx"),
  overwrite = TRUE
)

message("MetaCyc ALDEx2 results exported successfully.")

## =========================================================
## 11. Plot MetaCyc results
## =========================================================

if (nrow(aldex_picrust2_metacyc_sig) > 0) {
  p_metacyc <- make_effect_plot(
    df = aldex_picrust2_metacyc_sig,
    feature_col = "description",
    xlab = "MetaCyc pathway",
    title = "Differentially abundant MetaCyc pathways",
    top_n = min(20, nrow(aldex_picrust2_metacyc_sig))
  )
  
  ggsave(
    filename = file.path(out_dir_figures, "PICRUSt2_MetaCyc_effect_plot.pdf"),
    plot = p_metacyc,
    width = 8,
    height = 6
  )
}

## =========================================================
## 12. KEGG differential abundance analysis
## =========================================================

message("Running ALDEx2 for KEGG/KO features...")
aldex_picrust2_kegg <- run_aldex(kegg_abundance_int, conds)

aldex_picrust2_kegg_sig_raw <- subset(aldex_picrust2_kegg, wi.eBH < 0.05)
aldex_picrust2_kegg_sig_raw$p_adjust <- aldex_picrust2_kegg_sig_raw$wi.eBH
aldex_picrust2_kegg_sig_raw$feature <- rownames(aldex_picrust2_kegg_sig_raw)

aldex_picrust2_kegg_sig <- pathway_annotation(
  pathway = "KO",
  daa_results_df = aldex_picrust2_kegg_sig_raw,
  ko_to_kegg = TRUE
)

aldex_picrust2_kegg_sig <- aldex_picrust2_kegg_sig %>%
  arrange(desc(effect)) %>%
  mutate(condition_enriched = ifelse(effect > 0, "Post-Abx", "Pre-Abx"))

# Export full raw ALDEx2 KO results
aldex_picrust2_kegg_all <- as.data.frame(aldex_picrust2_kegg)
aldex_picrust2_kegg_all$feature <- rownames(aldex_picrust2_kegg_all)
aldex_picrust2_kegg_all <- aldex_picrust2_kegg_all %>%
  arrange(desc(effect)) %>%
  mutate(condition_enriched = ifelse(effect > 0, "Post-Abx", "Pre-Abx"))

write.xlsx(
  aldex_picrust2_kegg_all,
  file.path(out_dir_tables, "PICRUSt2_KEGG_ALDEx2_all.xlsx"),
  overwrite = TRUE
)

write.xlsx(
  aldex_picrust2_kegg_sig,
  file.path(out_dir_tables, "PICRUSt2_KEGG_ALDEx2_significant.xlsx"),
  overwrite = TRUE
)

message("KEGG ALDEx2 results exported successfully.")

## =========================================================
## 13. Plot KEGG results
## =========================================================

if (nrow(aldex_picrust2_kegg_sig) > 0) {
  kegg_plot_df <- aldex_picrust2_kegg_sig %>%
    filter(!is.na(pathway_name))
  
  if (nrow(kegg_plot_df) > 0) {
    p_kegg <- make_effect_plot(
      df = kegg_plot_df,
      feature_col = "pathway_name",
      xlab = "KEGG pathway",
      title = "Differentially abundant KEGG pathways",
      top_n = min(20, nrow(kegg_plot_df))
    )
    
    ggsave(
      filename = file.path(out_dir_figures, "PICRUSt2_KEGG_effect_plot.pdf"),
      plot = p_kegg,
      width = 8,
      height = 6
    )
  }
}

## =========================================================
## 14. Save result objects
## =========================================================

saveRDS(
  aldex_picrust2_metacyc,
  file.path(out_dir_intermediate, "PICRUSt2_MetaCyc_ALDEx2.rds")
)

saveRDS(
  aldex_picrust2_kegg,
  file.path(out_dir_intermediate, "PICRUSt2_KEGG_ALDEx2.rds")
)

saveRDS(
  aldex_picrust2_metacyc_sig,
  file.path(out_dir_intermediate, "PICRUSt2_MetaCyc_ALDEx2_significant.rds")
)

saveRDS(
  aldex_picrust2_kegg_sig,
  file.path(out_dir_intermediate, "PICRUSt2_KEGG_ALDEx2_significant.rds")
)

message("PICRUSt2 result objects saved successfully.")

## =========================================================
## 15. Save session information
## =========================================================

session_file <- file.path(out_dir_env, "sessionInfo_16S.txt")

writeLines(
  capture.output(sessionInfo()),
  session_file
)

message("Session information saved to: ", session_file)