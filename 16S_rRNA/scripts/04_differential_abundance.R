############################################################
# Script: 04_differential_abundance.R
# Purpose: Identify differentially abundant taxa between
#          Pre-Abx and Post-Abx samples using ALDEx2 across
#          multiple taxonomic ranks.
# Input: processed MPSE object with curated sample metadata
# Output: ALDEx2 result tables, filtered significant taxa,
#         and effect size plots
# Author: Bruna Toledo
############################################################

## =========================================================
## 1. Load libraries
## =========================================================

suppressPackageStartupMessages({
  library(MicrobiotaProcess)
  library(SummarizedExperiment)
  library(ALDEx2)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(openxlsx)
})

## =========================================================
## 2. Define input and output paths
## =========================================================

input_file <- "16S_rRNA/results/intermediate/mpse_taxonomy.rds"

out_dir_tables  <- "16S_rRNA/results/tables/differential_abundance"
out_dir_figures <- "16S_rRNA/results/figures/differential_abundance"
out_dir_intermediate <- "16S_rRNA/results/intermediate"

if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)
if (!dir.exists(out_dir_figures)) dir.create(out_dir_figures, recursive = TRUE)
if (!dir.exists(out_dir_intermediate)) dir.create(out_dir_intermediate, recursive = TRUE)

## =========================================================
## 3. Load processed MPSE object
## =========================================================

mpse <- readRDS(input_file)

message("Processed MPSE object loaded successfully.")

## =========================================================
## 4. Extract count table and taxonomy table
## =========================================================

otu_table <- assay(mpse, "RareAbundance") %>%
  as.data.frame() %>%
  rownames_to_column("OTU")

tax_table <- mp_extract_taxonomy(mpse) %>%
  as.data.frame()

sample_metadata <- mp_extract_sample(mpse) %>%
  as.data.frame()

## =========================================================
## 5. Filter low-prevalence features
## =========================================================

# Keep taxa present in at least 40% of samples
min_prev <- ceiling(0.4 * (ncol(otu_table) - 1))
message("Minimum prevalence threshold: ", min_prev, " samples")

otu_counts <- otu_table %>%
  column_to_rownames("OTU")

keep <- rowSums(otu_counts > 0) >= min_prev
message("Keeping ", sum(keep), " ASVs out of ", nrow(otu_counts))

otu_counts_filt <- otu_counts[keep, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("OTU")

## =========================================================
## 6. Define helper function to aggregate counts by taxonomy
## =========================================================

aggregate_by_tax <- function(otu_df, tax_df, rank) {
  otu_df %>%
    left_join(tax_df %>% select(OTU, !!sym(rank)), by = "OTU") %>%
    filter(!is.na(.data[[rank]]), .data[[rank]] != "") %>%
    select(-OTU) %>%
    group_by(.data[[rank]]) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
    column_to_rownames(var = rank) %>%
    as.matrix()
}

## =========================================================
## 7. Build taxonomic count matrices
## =========================================================

matrix_list <- list(
  Phylum = aggregate_by_tax(otu_counts_filt, tax_table, "Phylum"),
  Class  = aggregate_by_tax(otu_counts_filt, tax_table, "Class"),
  Order  = aggregate_by_tax(otu_counts_filt, tax_table, "Order"),
  Family = aggregate_by_tax(otu_counts_filt, tax_table, "Family"),
  Genus  = aggregate_by_tax(otu_counts_filt, tax_table, "Genus")
)

# Ensure rownames are present
for (level in names(matrix_list)) {
  mat <- matrix_list[[level]]
  
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0(level, "_", seq_len(nrow(mat)))
  }
  
  matrix_list[[level]] <- mat
}

message("Taxonomic count matrices created successfully.")

## =========================================================
## 8. Define conditions from sample metadata
## =========================================================

# Ensure sample order matches matrix column order
sample_order <- colnames(matrix_list[[1]])
sample_metadata <- sample_metadata[sample_order, , drop = FALSE]

conds <- factor(
  sample_metadata$condition,
  levels = c("Pre-Abx", "Post-Abx")
)

message("Condition vector prepared successfully.")

## =========================================================
## 9. Run ALDEx2 across taxonomic levels
## =========================================================

aldex_results <- list()

for (level in names(matrix_list)) {
  message("Running ALDEx2 for level: ", level)
  
  counts <- as.data.frame(matrix_list[[level]])
  counts[] <- lapply(counts, function(col) as.numeric(as.character(col)))
  
  aldex_res <- aldex(
    reads = counts,
    conditions = conds,
    mc.samples = 128,
    test = "t",
    paired.test = TRUE,
    effect = TRUE,
    denom = "all"
  )
  
  aldex_results[[level]] <- aldex_res
}

message("ALDEx2 completed for all taxonomic levels.")

## =========================================================
## 10. Export complete ALDEx2 results
## =========================================================

export_aldex_to_excel <- function(result_list, file_name) {
  wb <- createWorkbook()
  
  for (level in names(result_list)) {
    addWorksheet(wb, sheetName = level)
    writeData(
      wb,
      sheet = level,
      x = as.data.frame(result_list[[level]]) %>%
        rownames_to_column("feature")
    )
  }
  
  saveWorkbook(wb, file = file_name, overwrite = TRUE)
}

all_results_file <- file.path(out_dir_tables, "aldex2_all_levels.xlsx")
export_aldex_to_excel(aldex_results, all_results_file)

message("Complete ALDEx2 results exported to: ", all_results_file)

## =========================================================
## 11. Filter significant ALDEx2 results
## =========================================================

filter_significant_aldex <- function(aldex_result,
                                     p_thresh = 0.05,
                                     effect_thresh = 1) {
  df <- as.data.frame(aldex_result)
  df$feature <- rownames(df)
  
  df %>%
    filter(wi.eBH < p_thresh, abs(effect) > effect_thresh) %>%
    arrange(wi.eBH)
}

filtered_aldex_results <- list()

for (level in names(aldex_results)) {
  filtered_df <- filter_significant_aldex(aldex_results[[level]])
  
  if (nrow(filtered_df) > 0) {
    filtered_aldex_results[[level]] <- filtered_df
    message(level, ": ", nrow(filtered_df), " significant features found.")
  } else {
    message(level, ": no significant features found.")
  }
}

## =========================================================
## 12. Export filtered significant results
## =========================================================

filtered_results_file <- file.path(out_dir_tables, "aldex2_significant_levels.xlsx")
export_aldex_to_excel(filtered_aldex_results, filtered_results_file)

message("Filtered ALDEx2 results exported to: ", filtered_results_file)

## =========================================================
## 13. Define plotting function
## =========================================================

plot_aldex_effects <- function(filtered_df, level_name) {
  if (nrow(filtered_df) == 0) {
    return(NULL)
  }
  
  filtered_df <- filtered_df %>%
    mutate(
      direction = ifelse(effect > 0, "Post-Abx", "Pre-Abx"),
      direction = factor(direction, levels = c("Pre-Abx", "Post-Abx"))
    )
  
  ggplot(filtered_df, aes(x = reorder(feature, effect), y = effect, fill = direction)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("Differential abundance at", level_name, "level"),
      x = "",
      y = "Effect size",
      fill = ""
    ) +
    scale_fill_manual(values = c("Pre-Abx" = "#3C548866", "Post-Abx" = "#87255B99")) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

## =========================================================
## 14. Generate and save plots
## =========================================================

for (level in names(filtered_aldex_results)) {
  p <- plot_aldex_effects(filtered_aldex_results[[level]], level)
  
  if (!is.null(p)) {
    ggsave(
      filename = file.path(out_dir_figures, paste0("aldex2_", tolower(level), "_effect_plot.pdf")),
      plot = p,
      width = 7,
      height = 5
    )
  }
}

message("ALDEx2 effect size plots exported successfully.")

## =========================================================
## 15. Save results object for downstream use
## =========================================================

saveRDS(
  aldex_results,
  file.path(out_dir_intermediate, "aldex2_results.rds")
)

saveRDS(
  filtered_aldex_results,
  file.path(out_dir_intermediate, "aldex2_significant_results.rds")
)

message("ALDEx2 result objects saved successfully.")