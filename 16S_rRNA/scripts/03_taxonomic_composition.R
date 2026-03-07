############################################################
# Script: 03_taxonomic_composition.R
# Purpose: Calculate and visualize taxonomic composition
#          across samples and experimental groups using
#          relative abundance profiles from rarefied data.
# Input: processed MPSE object from diversity workflow
# Output: taxonomic composition plots and relative abundance
#         tables exported to .xlsx
# Author: Bruna Toledo
############################################################

## =========================================================
## 1. Load libraries
## =========================================================

suppressPackageStartupMessages({
  library(MicrobiotaProcess)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(openxlsx)
})

## =========================================================
## 2. Define input and output paths
## =========================================================

input_file <- "16S_rRNA/results/intermediate/mpse_diversity.rds"

out_dir_tables  <- "16S_rRNA/results/tables/taxonomic_composition"
out_dir_figures <- "16S_rRNA/results/figures/taxonomic_composition"

if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)
if (!dir.exists(out_dir_figures)) dir.create(out_dir_figures, recursive = TRUE)

## =========================================================
## 3. Load processed MPSE object
## =========================================================

mpse <- readRDS(input_file)

message("Processed MPSE object loaded successfully.")

## =========================================================
## 4. Calculate relative abundance
## =========================================================

mpse <- mpse %>%
  mp_cal_abundance(
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance(
    .abundance = RareAbundance,
    .group = condition
  )

message("Relative abundance calculated by sample and by condition.")

## =========================================================
## 5. Define helper function to extract relative abundance
## =========================================================

extract_relabund_tables <- function(mpse_obj,
                                    tax_level = "Phylum",
                                    abundance_assay = "RareAbundance",
                                    group_var = "condition") {
  
  relabund_df <- tryCatch({
    mp_extract_abundance(
      mpse_obj,
      .abundance = abundance_assay,
      taxa.class = !!rlang::sym(tax_level)
    )
  }, error = function(e) {
    stop("Failed to extract abundance data for taxonomic level: ", tax_level)
  })
  
  sample_col <- paste0(abundance_assay, "BySample")
  group_col  <- paste0(abundance_assay, "By", group_var)
  
  if (!sample_col %in% colnames(relabund_df)) {
    stop("Column not found: ", sample_col)
  }
  if (!group_col %in% colnames(relabund_df)) {
    stop("Column not found: ", group_col)
  }
  
  relabund_sample <- relabund_df %>%
    transmute(
      Taxon = label,
      data = map(.data[[sample_col]], ~ select(.x, Sample, starts_with("Rel")))
    ) %>%
    unnest(data) %>%
    pivot_wider(
      names_from = Sample,
      values_from = starts_with("Rel")
    )
  
  relabund_group <- relabund_df %>%
    transmute(
      Taxon = label,
      data = map(.data[[group_col]], ~ select(.x, !!sym(group_var), starts_with("Rel")))
    ) %>%
    unnest(data) %>%
    pivot_wider(
      names_from = !!sym(group_var),
      values_from = starts_with("Rel")
    )
  
  relabund_sample[] <- lapply(relabund_sample, function(x) {
    if (is.numeric(x)) format(x, scientific = FALSE, digits = 6) else x
  })
  
  relabund_group[] <- lapply(relabund_group, function(x) {
    if (is.numeric(x)) format(x, scientific = FALSE, digits = 6) else x
  })
  
  list(
    per_sample = relabund_sample,
    per_group  = relabund_group
  )
}

## =========================================================
## 6. Extract relative abundance tables
## =========================================================

phylum_tables <- extract_relabund_tables(mpse, tax_level = "Phylum")
class_tables  <- extract_relabund_tables(mpse, tax_level = "Class")
order_tables  <- extract_relabund_tables(mpse, tax_level = "Order")
family_tables <- extract_relabund_tables(mpse, tax_level = "Family")
genus_tables  <- extract_relabund_tables(mpse, tax_level = "Genus")

taxa_tables <- list(
  Phylum = phylum_tables,
  Class  = class_tables,
  Order  = order_tables,
  Family = family_tables,
  Genus  = genus_tables
)

message("Relative abundance tables extracted successfully.")

## =========================================================
## 7. Export relative abundance tables to Excel
## =========================================================

for (taxon in names(taxa_tables)) {
  wb <- createWorkbook()
  
  addWorksheet(wb, "per_sample")
  writeData(wb, "per_sample", taxa_tables[[taxon]]$per_sample)
  
  addWorksheet(wb, "per_group")
  writeData(wb, "per_group", taxa_tables[[taxon]]$per_group)
  
  outfile <- file.path(out_dir_tables, paste0("relative_abundance_", taxon, ".xlsx"))
  saveWorkbook(wb, file = outfile, overwrite = TRUE)
}

message("Relative abundance Excel files exported successfully.")

## =========================================================
## 8. Define plotting order and labels
## =========================================================

plot_order <- c(
  "S2_T1","S3_T1","S4_T1","S5_T1","S6_T1","S7_T1","S8_T1","S9_T1",
  "S10_T1","S11_T1","S12_T1","S13_T1","S15_T1","S16_T1",
  "S2_T2","S3_T2","S4_T2","S5_T2","S6_T2","S7_T2","S8_T2","S9_T2",
  "S10_T2","S11_T2","S12_T2","S13_T2","S15_T2","S16_T2"
)

plot_label <- c(
  "S2_T1"  = "S2",  "S3_T1"  = "S3",  "S4_T1"  = "S4",  "S5_T1"  = "S5",
  "S6_T1"  = "S6",  "S7_T1"  = "S7",  "S8_T1"  = "S8",  "S9_T1"  = "S9",
  "S10_T1" = "S10", "S11_T1" = "S11", "S12_T1" = "S12", "S13_T1" = "S13",
  "S15_T1" = "S15", "S16_T1" = "S16",
  "S2_T2"  = "S2",  "S3_T2"  = "S3",  "S4_T2"  = "S4",  "S5_T2"  = "S5",
  "S6_T2"  = "S6",  "S7_T2"  = "S7",  "S8_T2"  = "S8",  "S9_T2"  = "S9",
  "S10_T2" = "S10", "S11_T2" = "S11", "S12_T2" = "S12", "S13_T2" = "S13",
  "S15_T2" = "S15", "S16_T2" = "S16"
)

## =========================================================
## 9. Define helper function for abundance plots
## =========================================================

make_abundance_plots <- function(mpse_obj, tax_level, topn) {
  
  p_sample <- mpse_obj %>%
    mp_plot_abundance(
      .abundance = RareAbundance,
      .group = condition,
      taxa.class = !!rlang::sym(tax_level),
      topn = topn,
      relative = TRUE,
      order.by.feature = FALSE
    ) +
    scale_x_discrete(limits = plot_order, labels = plot_label)
  
  p_group <- mpse_obj %>%
    mp_plot_abundance(
      .abundance = RareAbundance,
      .group = condition,
      taxa.class = !!rlang::sym(tax_level),
      topn = topn,
      plot.group = TRUE
    )
  
  list(
    per_sample = p_sample,
    per_group  = p_group
  )
}

## =========================================================
## 10. Generate plots for each taxonomic level
## =========================================================

phylum_plots <- make_abundance_plots(mpse, "Phylum", topn = 7)
class_plots  <- make_abundance_plots(mpse, "Class",  topn = 10)
order_plots  <- make_abundance_plots(mpse, "Order",  topn = 10)
family_plots <- make_abundance_plots(mpse, "Family", topn = 10)
genus_plots  <- make_abundance_plots(mpse, "Genus",  topn = 20)

## =========================================================
## 11. Save plots
## =========================================================

ggsave(
  filename = file.path(out_dir_figures, "taxonomic_composition_phylum_per_sample.pdf"),
  plot = phylum_plots$per_sample,
  width = 10,
  height = 5
)

ggsave(
  filename = file.path(out_dir_figures, "taxonomic_composition_phylum_per_group.pdf"),
  plot = phylum_plots$per_group,
  width = 5,
  height = 4
)

ggsave(
  filename = file.path(out_dir_figures, "taxonomic_composition_class_per_sample.pdf"),
  plot = class_plots$per_sample,
  width = 10,
  height = 5
)

ggsave(
  filename = file.path(out_dir_figures, "taxonomic_composition_class_per_group.pdf"),
  plot = class_plots$per_group,
  width = 5,
  height = 4
)

ggsave(
  filename = file.path(out_dir_figures, "taxonomic_composition_order_per_sample.pdf"),
  plot = order_plots$per_sample,
  width = 10,
  height = 5
)

ggsave(
  filename = file.path(out_dir_figures, "taxonomic_composition_order_per_group.pdf"),
  plot = order_plots$per_group,
  width = 5,
  height = 4
)

ggsave(
  filename = file.path(out_dir_figures, "taxonomic_composition_family_per_sample.pdf"),
  plot = family_plots$per_sample,
  width = 10,
  height = 5
)

ggsave(
  filename = file.path(out_dir_figures, "taxonomic_composition_family_per_group.pdf"),
  plot = family_plots$per_group,
  width = 5,
  height = 4
)

ggsave(
  filename = file.path(out_dir_figures, "taxonomic_composition_genus_per_sample.pdf"),
  plot = genus_plots$per_sample,
  width = 10,
  height = 6
)

ggsave(
  filename = file.path(out_dir_figures, "taxonomic_composition_genus_per_group.pdf"),
  plot = genus_plots$per_group,
  width = 5,
  height = 5
)

message("Taxonomic composition plots exported successfully.")

## =========================================================
## 12. Save processed object for downstream analyses
## =========================================================

saveRDS(mpse, file.path("16S_rRNA/results/intermediate", "mpse_taxonomy.rds"))

message("Processed MPSE object saved to: 16S_rRNA/results/intermediate/mpse_taxonomy.rds")