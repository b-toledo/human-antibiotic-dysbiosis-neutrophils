############################################################
# Script: 04_heatmap_accessibility.R
# Purpose: Generate heatmaps of normalized ATAC-seq
#          accessibility aggregated at gene level.
# Input: normalized counts, peak annotations, metadata
# Output: 4_Heatmap/Ht-ATAC.svg
# Author: Bruna Toledo
############################################################

############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
})

############################################################
# 2. Define input/output paths
############################################################

out_dir <- "4_Heatmap"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

heatmap_file <- file.path(out_dir, "Ht-ATAC.svg")

############################################################
# 3. Check required objects
############################################################

required_objects <- c("peak_annotation_df", "norm_counts", "metadata")

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "The following objects are missing from the environment: ",
    paste(missing_objects, collapse = ", "),
    ". Please load the annotated DAR workspace before running this script."
  )
}

############################################################
# 4. Prepare gene-level accessibility matrix
############################################################

anno_sub <- peak_annotation_df %>%
  dplyr::select(peak_id, SYMBOL) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  distinct()

counts_annot <- norm_counts[rownames(norm_counts) %in% anno_sub$peak_id, , drop = FALSE]

row_map <- anno_sub[match(rownames(counts_annot), anno_sub$peak_id), ]
stopifnot(nrow(row_map) == nrow(counts_annot))

counts_gene <- cbind(SYMBOL = row_map$SYMBOL, as.data.frame(counts_annot)) %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

counts_gene <- as.data.frame(counts_gene)
rownames(counts_gene) <- counts_gene$SYMBOL
counts_gene$SYMBOL <- NULL

############################################################
# 5. Define gene-level genomic region annotation
############################################################

region_by_gene <- peak_annotation_df %>%
  transmute(SYMBOL, region_type) %>%
  distinct() %>%
  group_by(SYMBOL) %>%
  summarise(
    region_type = {
      r <- unique(region_type)
      if ("Promoter" %in% r) {
        "Promoter"
      } else if ("Distal" %in% r) {
        "Distal"
      } else {
        "Trans"
      }
    },
    .groups = "drop"
  )

region_by_gene <- region_by_gene[region_by_gene$SYMBOL %in% rownames(counts_gene), ]
region_by_gene <- region_by_gene[match(rownames(counts_gene), region_by_gene$SYMBOL), , drop = FALSE]

stopifnot(identical(rownames(counts_gene), region_by_gene$SYMBOL))

############################################################
# 6. Align samples and z-score matrix
############################################################

if (all(metadata$sample %in% colnames(counts_gene))) {
  counts_gene <- counts_gene[, metadata$sample, drop = FALSE]
}

mat <- as.matrix(counts_gene)
mat_z <- t(scale(t(mat)))
mat_z[!is.finite(mat_z)] <- 0

############################################################
# 7. Build annotations
############################################################

anno_col_df <- data.frame(
  Condition = factor(metadata$condition, levels = c("Pre-Abx", "Post-Abx")),
  row.names = metadata$sample
)

anno_col_df <- anno_col_df[colnames(mat_z), , drop = FALSE]

anno_col <- HeatmapAnnotation(
  Condition = anno_col_df$Condition,
  col = list(
    Condition = c(
      "Pre-Abx"  = "#804D8066",
      "Post-Abx" = "#804D80cc"
    )
  ),
  annotation_legend_param = list(
    title_gp  = gpar(fontfamily = "Arial", fontsize = 6),
    labels_gp = gpar(fontfamily = "Arial", fontsize = 6)
  ),
  show_annotation_name = FALSE
)

region_fac <- factor(region_by_gene$region_type, levels = c("Promoter", "Distal", "Trans"))

reg_col <- c(
  Promoter = "#d97b29cc",
  Distal   = "#5ca972cc",
  Trans    = "#7a68b3cc"
)

anno_row <- rowAnnotation(
  Reg = region_fac,
  col = list(Reg = reg_col),
  annotation_legend_param = list(
    Reg = list(title = "Genome region"),
    title_gp  = gpar(fontfamily = "Arial", fontsize = 6),
    labels_gp = gpar(fontfamily = "Arial", fontsize = 6)
  )
)

############################################################
# 8. Define heatmap aesthetics
############################################################

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("royalblue2", "white", "firebrick2")
)

genes_to_label <- c(
  "IL1B", "SOD2", "CXCR2", "ELANE", "ITGB3", "IL10", "AKT3", "IRF5", "VAMP7",
  "SIGLEC9", "HCK", "GPR183", "LTA4H", "CD80", "NCOR2", "HDAC3",
  "RELB", "CXCR4", "ITGAM", "JAK1", "MAPK13", "IRF1", "NOD1", "CASP8", "EP300",
  "FOS", "JUNB", "FOSL2", "IRF2", "IRF1", "STAT2", "MAP3K8", "KAT7", "KAT8",
  "EP300", "SP1", "SP2"
)

############################################################
# 9. Generate heatmap
############################################################

ht <- Heatmap(
  mat_z,
  name = "Z-score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  top_annotation = anno_col,
  left_annotation = anno_row,
  show_row_names = TRUE,
  row_labels = ifelse(rownames(mat_z) %in% genes_to_label, rownames(mat_z), ""),
  row_names_gp = gpar(fontsize = 5.5),
  show_column_names = FALSE,
  clustering_distance_rows = "manhattan",
  clustering_distance_columns = "manhattan",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  heatmap_legend_param = list(
    legend_direction = "vertical",
    color_bar = "continuous",
    legend_width = unit(4, "cm"),
    title_gp = gpar(fontfamily = "Arial", fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontfamily = "Arial", fontsize = 6)
  )
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE
)

############################################################
# 10. Save heatmap
############################################################

svg(heatmap_file, width = 4, height = 4)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE
)

dev.off()

############################################################
# 11. Save workspace
############################################################

save.image(file = "Heatmap.RData")



############################################################
# End of script
############################################################