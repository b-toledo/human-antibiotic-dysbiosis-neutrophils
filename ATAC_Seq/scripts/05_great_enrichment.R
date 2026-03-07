############################################################
# Script: 05_great_enrichment.R
# Purpose: Perform GREAT enrichment analysis on
#          differentially accessible ATAC-seq regions.
# Input: annotated DARs and consensus peak coordinates
# Output: 5_GREAT/GREAT_enrichment.xlsx
# Author: Bruna Toledo
############################################################

############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(rGREAT)
  library(writexl)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

############################################################
# 2. Define input/output paths
############################################################

out_dir <- "5_GREAT"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

great_output_file <- file.path(out_dir, "GREAT_enrichment.xlsx")

############################################################
# 3. Check required objects
############################################################

required_objects <- c("dar_annot", "nPeaks")

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "The following objects are missing from the environment: ",
    paste(missing_objects, collapse = ", "),
    ". Please load the annotated DAR workspace before running this script."
  )
}

############################################################
# 4. Define helper functions
############################################################

dar_to_gr <- function(df) {
  GRanges(
    seqnames = df$seqnames,
    ranges   = IRanges(start = df$start, end = df$end)
  )
}

combine_tables_by_direction <- function(results_list) {
  dirs <- unique(sapply(results_list, function(x) x$direction))
  
  combined <- list()
  
  for (dir in dirs) {
    items <- results_list[sapply(results_list, function(x) x$direction == dir)]
    
    items <- items[!sapply(items, function(x) {
      is.null(x$table) || nrow(x$table) == 0
    })]
    
    if (length(items) == 0) next
    
    df <- bind_rows(lapply(items, function(x) {
      x$table %>% mutate(ontology = x$ontology)
    }))
    
    combined[[dir]] <- df
  }
  
  combined
}

############################################################
# 5. Define background peak set
############################################################

bk_peaks <- GRanges(
  seqnames = nPeaks$Chr,
  ranges   = IRanges(start = nPeaks$Start, end = nPeaks$End)
)

############################################################
# 6. Define GREAT settings
############################################################

ontologies <- c("GO:BP", "GO:CC", "GO:MP", "msigdb:H")
directions <- c("Up", "Down")

############################################################
# 7. Run GREAT enrichment
############################################################

great_results <- list()

for (dir in directions) {
  df_dir <- dar_annot %>% filter(direction == dir)
  
  if (nrow(df_dir) == 0) {
    warning("No peaks for direction ", dir)
    next
  }
  
  gr_dir <- dar_to_gr(df_dir)
  
  for (ont in ontologies) {
    res <- tryCatch(
      {
        great(
          gr_dir,
          ont,
          "TxDb.Hsapiens.UCSC.hg38.knownGene",
          basal_downstream = 2000,
          background = bk_peaks
        )
      },
      error = function(e) {
        message("GREAT failed for ", dir, " / ", ont, ": ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(res)) next
    
    tb <- getEnrichmentTable(res)
    
    great_results[[paste(dir, ont, sep = "_")]] <- list(
      direction    = dir,
      ontology     = ont,
      great_object = res,
      table        = tb
    )
  }
}

############################################################
# 8. Combine GREAT results by direction
############################################################

combined_direction_tables <- combine_tables_by_direction(great_results)

############################################################
# 9. Export results
############################################################

write_xlsx(
  combined_direction_tables,
  path = great_output_file
)

############################################################
# 10. Save workspace
############################################################

save.image(file = "GREAT.RData")



############################################################
# End of script
############################################################