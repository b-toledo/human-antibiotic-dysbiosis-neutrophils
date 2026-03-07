############################################################
# Script: 01_build_phyloseq.R
# Purpose: Import ASV tables, taxonomy assignments, and sample
#          metadata to construct a phyloseq object for downstream
#          microbiome analyses.
# Input: dada2_phyloseq.rds, tree.nwk
# Output: curated phyloseq and MPSE objects saved as .rds files;
#         optional ASV abundance tables exported to .xlsx
# Author: Bruna Toledo
############################################################


## =========================================================
## 1. Load libraries
## =========================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(MicrobiotaProcess)
  library(mia)
  library(qiime2R)
  library(dplyr)
  library(tibble)
  library(openxlsx)
})

## =========================================================
## 2. Define input and output paths
## =========================================================

phyloseq_path <- "16S_rRNA/pipeline/Nextflow/phyloseq/dada2_phyloseq.rds"
tree_path     <- "16S_rRNA/pipeline/Nextflow/qiime2/phylogenetic_tree/tree.nwk"

out_dir          <- "16S_rRNA/results/intermediate"
out_tables_dir   <- "16S_rRNA/results/tables"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(out_tables_dir)) dir.create(out_tables_dir, recursive = TRUE)

## =========================================================
## 3. Define sample naming and ordering
## =========================================================

# "S3T1" and "S3T2" samples were switched in the original input.
# The mapping below restores the correct sample labels.
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

sample_order <- c(
  "S2_T1",  "S2_T2",  "S3_T1",  "S3_T2",  "S4_T1",  "S4_T2",
  "S5_T1",  "S5_T2",  "S6_T1",  "S6_T2",  "S7_T1",  "S7_T2",
  "S8_T1",  "S8_T2",  "S9_T1",  "S9_T2",  "S10_T1", "S10_T2",
  "S11_T1", "S11_T2", "S12_T1", "S12_T2", "S13_T1", "S13_T2",
  "S15_T1", "S15_T2", "S16_T1", "S16_T2"
)

subject_order <- c(
  "S2", "S3", "S4", "S5", "S6", "S7", "S8",
  "S9", "S10", "S11", "S12", "S13", "S15", "S16"
)

## =========================================================
## 4. Import phyloseq object and phylogenetic tree
## =========================================================

ps_raw <- readRDS(phyloseq_path)
tree   <- read_tree(tree_path)

ps <- merge_phyloseq(ps_raw, tree)

## =========================================================
## 5. Curate taxonomy table
## =========================================================

ps <- ps %>%
  tax_mutate(
    Species_exact = NULL,
    confidence    = NULL,
    sequence      = NULL,
    Species = if_else(
      Species != "",
      paste(Genus, Species, sep = "_"),
      Species
    )
  )

message("Phyloseq object successfully imported and curated.")

## =========================================================
## 6. Convert phyloseq to MPSE object
## =========================================================

mpse <- as.MPSE(ps)

if (!exists("mpse")) {
  stop("MPSE object could not be created. Check the phyloseq object.")
}

message("MPSE object successfully created.")

## =========================================================
## 7. Recode condition labels
## =========================================================

condition_vec <- mpse@colData@listData[["condition"]]
condition_vec[condition_vec == "Control"] <- "Pre-Abx"
condition_vec[condition_vec == "Treated"] <- "Post-Abx"

mpse@colData@listData[["condition"]] <- factor(
  condition_vec,
  levels = c("Pre-Abx", "Post-Abx")
)

## =========================================================
## 8. Rename sample IDs
## =========================================================

old_sample_names <- mpse@colData@rownames
new_sample_names <- unname(sample_labels[old_sample_names])

if (any(is.na(new_sample_names))) {
  missing_ids <- old_sample_names[is.na(new_sample_names)]
  stop("Missing sample label mapping for: ", paste(missing_ids, collapse = ", "))
}

# Rename colData rownames
mpse@colData@rownames <- new_sample_names

# Rename assay column names
for (assay_name in names(mpse@assays@data@listData)) {
  assay <- mpse@assays@data@listData[[assay_name]]
  
  if (!is.null(colnames(assay))) {
    assay_old_names <- colnames(assay)
    assay_new_names <- unname(sample_labels[assay_old_names])
    
    if (any(is.na(assay_new_names))) {
      missing_ids <- assay_old_names[is.na(assay_new_names)]
      stop(
        "Missing assay sample label mapping in assay '",
        assay_name, "' for: ",
        paste(missing_ids, collapse = ", ")
      )
    }
    
    colnames(assay) <- assay_new_names
  }
  
  mpse@assays@data@listData[[assay_name]] <- assay
}

## =========================================================
## 9. Add subject variable and reorder samples
## =========================================================

subject_map <- setNames(
  gsub("_T[12]$", "", sample_order),
  sample_order
)

# Reorder colData
mpse@colData@rownames <- sample_order

# Reorder assay columns
for (assay_name in names(mpse@assays@data@listData)) {
  assay <- mpse@assays@data@listData[[assay_name]]
  
  if (!is.null(colnames(assay))) {
    assay <- assay[, sample_order, drop = FALSE]
  }
  
  mpse@assays@data@listData[[assay_name]] <- assay
}

# Add ordered subject factor
mpse@colData@listData[["subject"]] <- factor(
  subject_map[mpse@colData@rownames],
  levels = subject_order
)

message("Sample names, condition labels, and subject order successfully curated.")

## =========================================================
## 10. Export abundance tables
## =========================================================

ASV_Ab <- as.data.frame(mpse@assays@data@listData[["Abundance"]]) %>%
  rownames_to_column(var = "ASV")

wb <- createWorkbook()
addWorksheet(wb, "Unrarefied-Abundance")
writeData(wb, "Unrarefied-Abundance", ASV_Ab)

xlsx_out <- file.path(out_tables_dir, "ASV_tables.xlsx")
saveWorkbook(wb, xlsx_out, overwrite = TRUE)

message("ASV abundance table exported to: ", xlsx_out)

## =========================================================
## 11. Save curated objects
## =========================================================

saveRDS(ps,   file.path(out_dir, "phyloseq_curated.rds"))
saveRDS(mpse, file.path(out_dir, "mpse_curated.rds"))

message("Curated objects saved successfully:")
message("- ", file.path(out_dir, "phyloseq_curated.rds"))
message("- ", file.path(out_dir, "mpse_curated.rds"))