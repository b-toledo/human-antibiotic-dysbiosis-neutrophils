############################################################
# Script: 08_atac_rna_correlation.R
# Purpose: Integrate ATAC-seq DARs with RNA-seq DEGs and
#          summarize concordant and discordant changes.
# Input: annotated DARs and RNA-seq DEG results
# Output: 8_ATACvsRNA/
# Author: Bruna Toledo
############################################################

############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggrepel)
  library(ggplot2)
  library(writexl)
  library(openxlsx)
  library(readxl)
  library(scales)
  library(tibble)
  library(rGREAT)
  library(GenomicRanges)
  library(IRanges)
})

############################################################
# 2. Define input/output paths
############################################################

out_dir <- "8_ATACvsRNA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

combined_table_file <- file.path(out_dir, "Combined_DAR-DEG_k2.xlsx")
plot_svg_file       <- file.path(out_dir, "ATACvsRNA.svg")
plot_tiff_file      <- file.path(out_dir, "ATACvsRNA.tiff")

great_dir <- file.path(out_dir, "GREAT_by_Category")
dir.create(great_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
# 3. Check required objects
############################################################

required_objects <- c("PostAbx_vs_PreAbx", "peak_annotation_df", "nPeaks")

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "The following objects are missing from the environment: ",
    paste(missing_objects, collapse = ", "),
    ". Please load the ATAC and RNA differential results before running this script."
  )
}

############################################################
# 4. Define helper functions
############################################################

run_great_for_categories <- function(region_sets_df,
                                     background_df,
                                     output_dir = great_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  stopifnot(all(c("Chr", "Start", "End") %in% colnames(background_df)))
  
  bk_peaks <- GRanges(
    seqnames = background_df$Chr,
    ranges   = IRanges(start = background_df$Start, end = background_df$End)
  )
  
  ontologies <- c("GO:BP", "GO:CC", "GO:MP", "msigdb:H")
  
  all_results <- list()
  
  categories <- unique(region_sets_df$group)
  categories <- categories[!is.na(categories)]
  
  for (cat_name in categories) {
    message("Running GREAT for: ", cat_name)
    
    df_cat <- region_sets_df %>% filter(group == cat_name)
    
    if (nrow(df_cat) < 1) {
      message("Skipping ", cat_name, " - no peaks")
      next
    }
    
    req_cols <- c("seqnames", "start", "end")
    if (!all(req_cols %in% colnames(df_cat))) {
      stop("region_sets_df must contain: ", paste(req_cols, collapse = ", "))
    }
    
    gr_cat <- GRanges(
      seqnames = df_cat$seqnames,
      ranges   = IRanges(start = df_cat$start, end = df_cat$end)
    )
    
    cat_results <- list()
    
    for (ont in ontologies) {
      message("  Ontology: ", ont)
      
      great_obj <- tryCatch(
        great(
          gr_cat,
          ont,
          "TxDb.Hsapiens.UCSC.hg38.knownGene",
          basal_downstream = 2000,
          background = bk_peaks
        ),
        error = function(e) {
          message("  GREAT failed for ", cat_name, " / ", ont, ": ", e$message)
          return(NULL)
        }
      )
      
      great_table <- tryCatch(
        {
          if (!is.null(great_obj)) getEnrichmentTable(great_obj) else NULL
        },
        error = function(e) {
          message("  Could not extract GREAT table for ", cat_name, " / ", ont)
          return(NULL)
        }
      )
      
      cat_results[[ont]] <- list(
        great_object = great_obj,
        table = great_table
      )
    }
    
    all_results[[cat_name]] <- cat_results
    
    sheets <- lapply(cat_results, function(x) {
      if (is.null(x$table) || nrow(x$table) == 0) {
        data.frame(Message = "No enriched terms")
      } else {
        x$table
      }
    })
    
    safe_name <- gsub("[^A-Za-z0-9_]+", "_", cat_name)
    
    write_xlsx(
      sheets,
      path = file.path(output_dir, paste0(safe_name, "_GREAT.xlsx"))
    )
  }
  
  save(all_results, file = file.path(output_dir, "all_GREAT_results.RData"))
  
  invisible(all_results)
}

############################################################
# 5. Prepare RNA-seq DEG table
############################################################

deg <- PostAbx_vs_PreAbx
gene_ids <- deg$gene_id

deg <- deg %>%
  rename(SYMBOL = gene_id)

############################################################
# 6. Identify overlap between DAR-associated genes and DEGs
############################################################

dar_ids <- peak_annotation_df$SYMBOL

overlap <- intersect(dar_ids, gene_ids)
overlap_genes <- c(dar_ids, gene_ids)

############################################################
# 7. Build integrated DAR/DEG table
############################################################

combined_df <- peak_annotation_df %>%
  filter(SYMBOL %in% overlap_genes) %>%
  left_join(
    deg %>% dplyr::select(SYMBOL, RNA_log2FC = log2FoldChange),
    by = "SYMBOL"
  ) %>%
  dplyr::select(SYMBOL, geneId, RNA_log2FC, ATAC_logFC = logFC, region_type, seqnames, start, end) %>%
  filter(!is.na(RNA_log2FC), !is.na(ATAC_logFC), !is.na(region_type))

combined_df <- combined_df %>%
  mutate(
    category = case_when(
      RNA_log2FC > 0 & ATAC_logFC > 0 ~ "Upregulated in both",
      RNA_log2FC < 0 & ATAC_logFC < 0 ~ "Downregulated in both",
      RNA_log2FC > 0 & ATAC_logFC < 0 ~ "Opposite",
      RNA_log2FC < 0 & ATAC_logFC > 0 ~ "Opposite",
      TRUE ~ "Other"
    )
  )

combined_df <- combined_df %>%
  mutate(
    group = case_when(
      RNA_log2FC > 0 & ATAC_logFC > 0 ~ "Upregulated in both",
      RNA_log2FC < 0 & ATAC_logFC < 0 ~ "Downregulated in both",
      RNA_log2FC > 0 & ATAC_logFC < 0 ~ "Opposite (Up-RNA/Down-ATAC)",
      RNA_log2FC < 0 & ATAC_logFC > 0 ~ "Opposite (Down-RNA/Up-ATAC)",
      TRUE ~ "Other"
    )
  )

combined_df$color <- case_when(
  combined_df$category == "Upregulated in both" ~ "#ff333399",
  combined_df$category == "Downregulated in both" ~ "#0066ff99",
  combined_df$category == "Opposite" ~ "#00000099",
  combined_df$category == "Other" ~ "gray"
)

write_xlsx(combined_df, path = combined_table_file)

############################################################
# 8. Generate ATAC vs RNA scatter plot
############################################################

genes_to_label <- c(
  "IL1B","SOD2","NR4A1","LIF","PTGS2","NCOR2",
  "NCF1C","HCG27","CRYBG1","ORMDL2","GRHPR","IP6K1",
  "MSRB1","CSK","PLK3","ETS2","PLD1","EDN1","CSRNP1",
  "NPIPB4","TCF4","TRIB1","GADD45B"
)

label_df <- combined_df[combined_df$SYMBOL %in% genes_to_label, ]

p_atac_rna <- ggplot(combined_df, aes(x = ATAC_logFC, y = RNA_log2FC)) +
  geom_point(
    aes(color = group, shape = region_type),
    size = 3,
    alpha = 0.6,
    show.legend = c(color = TRUE, shape = TRUE)
  ) +
  geom_label_repel(
    data = label_df,
    aes(label = SYMBOL),
    size = 1.8,
    max.overlaps = 15,
    box.padding = 0.3,
    point.padding = 0.3,
    segment.color = "gray50",
    label.size = 0.25,
    label.padding = 0.25,
    fill = "white"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  scale_color_manual(
    values = c(
      "Upregulated in both" = "red2",
      "Downregulated in both" = "#0066ff",
      "Opposite (Up-RNA/Down-ATAC)" = "black",
      "Opposite (Down-RNA/Up-ATAC)" = "gray30",
      "Other" = "gray"
    ),
    name = "Direction",
    guide = guide_legend(override.aes = list(shape = 16, size = 4, alpha = 0.6))
  ) +
  labs(
    x = "ATAC-seq log2 fold change (Accessibility)",
    y = "RNA-seq log2 fold change (Expression)"
  ) +
  theme_minimal(base_size = 6) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 6, hjust = 0.5, family = "Arial"),
    axis.title = element_text(size = 6, family = "Arial"),
    axis.text = element_text(size = 6, family = "Arial"),
    legend.title = element_text(size = 6, family = "Arial"),
    legend.text = element_text(size = 6, family = "Arial"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_atac_rna

############################################################
# 9. Save integration plot
############################################################

ggsave(
  plot_svg_file,
  p_atac_rna,
  width = 8,
  height = 6.5,
  units = "cm",
  dpi = 300,
  device = "svg"
)

ggsave(
  plot_tiff_file,
  p_atac_rna,
  width = 8,
  height = 6.5,
  units = "cm",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

############################################################
# 10. Summarize integrated categories
############################################################

summary_table <- combined_df %>%
  group_by(region_type, group) %>%
  summarise(n_genes = n(), .groups = "drop")

summary_table_cat <- combined_df %>%
  group_by(group) %>%
  summarise(n_genes = n(), .groups = "drop")

print(summary_table)
print(summary_table_cat)

############################################################
# 11. Prepare region sets for GREAT
############################################################

combined_regions_df <- peak_annotation_df %>%
  filter(SYMBOL %in% combined_df$SYMBOL) %>%
  left_join(
    combined_df %>% dplyr::select(SYMBOL, group),
    by = "SYMBOL"
  ) %>%
  filter(!is.na(group)) %>%
  dplyr::select(seqnames, start, end, SYMBOL, group, region_type)

head(combined_regions_df)

############################################################
# 12. Run GREAT enrichment by integration category
############################################################

great_results <- run_great_for_categories(
  region_sets_df = combined_regions_df,
  background_df = nPeaks,
  output_dir = great_dir
)

############################################################
# 13. Optional downstream summary plot from external table
############################################################

if (file.exists("GO-DARvsDEG.xlsx")) {
  GO_DARvsDEG <- read_excel("GO-DARvsDEG.xlsx")
  
  group_colors <- c(
    "Upregulated in both" = scales::alpha("red2", 0.6),
    "Downregulated in both" = scales::alpha("#0066ff", 0.55),
    "Opposite (Up-RNA/Down-ATAC)" = scales::alpha("black", 0.7),
    "Opposite (Down-RNA/Up-ATAC)" = scales::alpha("gray30", 0.55)
  )
  
  group_order <- c(
    "Upregulated in both",
    "Downregulated in both",
    "Opposite (Up-RNA/Down-ATAC)",
    "Opposite (Down-RNA/Up-ATAC)"
  )
  
  GO_DARvsDEG <- GO_DARvsDEG %>%
    mutate(
      Group = factor(Group, levels = group_order),
      Description = factor(
        Description,
        levels = rev(
          GO_DARvsDEG %>%
            arrange(Group, desc(zScore)) %>%
            pull(Description)
        )
      )
    )
  
  p_go_summary <- ggplot(GO_DARvsDEG, aes(x = Description, y = zScore, fill = Group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), hjust = -0.1, size = 3) +
    scale_fill_manual(values = group_colors, breaks = group_order) +
    coord_flip() +
    labs(
      x = "GO Biological Process",
      y = "Z-Score",
      fill = ""
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 20, 5, 5)
    )
  
  print(p_go_summary)
}

############################################################
# 14. Save workspace
############################################################

save.image(file = "ATACxRNA.RData")

############################################################
# 15. Save session information
############################################################

sink(paste0(out_prefix, "_sessionInfo.txt"))
print(sessionInfo())
sink()



############################################################
# End of script
############################################################