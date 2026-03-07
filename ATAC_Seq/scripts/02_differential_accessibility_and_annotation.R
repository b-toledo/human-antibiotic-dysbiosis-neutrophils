############################################################
# Script: 02_differential_accessibility_and_annotation.R
# Purpose: Identify differential accessible regions and
#          annotate ATAC-seq peaks by genomic region.
# Input: RUVseq ATAC workspace and consensus peak tables
# Output: 2_DAR-tables/
# Author: Bruna Toledo
############################################################

############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(ashr)
  library(readr)
  library(dplyr)
  library(DESeq2)
  library(writexl)
  library(stringr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(edgeR)
  library(ggplot2)
  library(scales)
})

############################################################
# 2. Define input/output paths
############################################################

out_dir <- "2_DAR-tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dar_table_file            <- file.path(out_dir, "DAR_cpm_k1.xlsx")
dar_annotated_file        <- file.path(out_dir, "DAR_cpm_k1_annotated.xlsx")
dar_summary_file          <- file.path(out_dir, "DAR_region_summary_k1.xlsx")
promoter_annotation_file  <- file.path(out_dir, "promoter_annotation.xlsx")
distal_annotation_file    <- file.path(out_dir, "distal_annotation.xlsx")
trans_annotation_file     <- file.path(out_dir, "trans_annotation.xlsx")
full_annotation_file      <- file.path(out_dir, "full_annotation.xlsx")
region_plot_file          <- file.path(out_dir, "DARs_by_region_and_regulation.svg")
distance_plot_file        <- file.path(out_dir, "DARs_distance_to_TSS.tiff")
volcano_region_file       <- file.path(out_dir, "DARs_volcano_by_region.tiff")
dar_rdata_file            <- "DAR_results_anno_k1.RData"
workspace_file            <- "Deseq2&Anno_k1.RData"

############################################################
# 3. Check required objects
############################################################

required_objects <- c("set_ruvr_nPeaks_k3", "set_ruvr_nPeaks_k2", "nPeaks", "metadata")

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "The following objects are missing from the environment: ",
    paste(missing_objects, collapse = ", "),
    ". Please load the ATAC RUVSeq workspace before running this script."
  )
}

############################################################
# 4. Define helper functions
############################################################

annotate_atac_peaks_regions <- function(
    dar_df,
    peak_df,
    promoter_output = promoter_annotation_file,
    distal_output   = distal_annotation_file,
    trans_output    = trans_annotation_file,
    full_output     = full_annotation_file
) {
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(dplyr)
  library(writexl)
  
  stopifnot(is.data.frame(dar_df), is.data.frame(peak_df))
  
  if (!"peak_id" %in% colnames(dar_df)) {
    dar_df$peak_id <- rownames(dar_df)
  }
  if (!"peak_id" %in% colnames(peak_df)) {
    peak_df$peak_id <- rownames(peak_df)
  }
  
  required_cols <- c("peak_id", "Chr", "Start", "End")
  missing_cols <- setdiff(required_cols, colnames(peak_df))
  if (length(missing_cols) > 0) {
    stop("peak_df is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  dar_peak <- peak_df[peak_df$peak_id %in% dar_df$peak_id, , drop = FALSE]
  if (nrow(dar_peak) == 0) {
    stop("No overlapping peak_id between dar_df and peak_df.")
  }
  
  dar_peak$Start <- as.integer(dar_peak$Start)
  dar_peak$End   <- as.integer(dar_peak$End)
  
  dar_gr <- GRanges(
    seqnames = dar_peak$Chr,
    ranges   = IRanges(start = dar_peak$Start, end = dar_peak$End),
    peak_id  = dar_peak$peak_id
  )
  
  peak_annotation_all <- annotatePeak(
    dar_gr,
    TxDb      = TxDb.Hsapiens.UCSC.hg38.knownGene,
    tssRegion = c(-1e6, 1e6),
    annoDb    = "org.Hs.eg.db"
  )
  
  peak_annotation_df <- as.data.frame(peak_annotation_all)
  
  in_promoter <- peak_annotation_df$distanceToTSS >= -2000 &
    peak_annotation_df$distanceToTSS <= 500
  abs_dist <- abs(peak_annotation_df$distanceToTSS)
  
  peak_annotation_df$region_type <- ifelse(
    in_promoter, "Promoter",
    ifelse(abs_dist <= 10000, "Distal", "Trans")
  )
  
  write_xlsx(peak_annotation_df, full_output)
  write_xlsx(peak_annotation_df %>% filter(region_type == "Promoter"), promoter_output)
  write_xlsx(peak_annotation_df %>% filter(region_type == "Distal"), distal_output)
  write_xlsx(peak_annotation_df %>% filter(region_type == "Trans"), trans_output)
  
  annotation_cols <- intersect(
    c(
      "seqnames", "start", "end", "width", "strand", "peak_id",
      "geneId", "SYMBOL", "GENENAME", "distanceToTSS", "region_type"
    ),
    colnames(peak_annotation_df)
  )
  
  dar_df_annotated <- dar_df %>%
    left_join(peak_annotation_df[, annotation_cols, drop = FALSE], by = "peak_id")
  
  list(
    dar_df_annotated    = dar_df_annotated,
    full_annotation     = peak_annotation_df,
    promoter_annotation = peak_annotation_df %>% filter(region_type == "Promoter"),
    distal_annotation   = peak_annotation_df %>% filter(region_type == "Distal"),
    trans_annotation    = peak_annotation_df %>% filter(region_type == "Trans")
  )
}

log2p1 <- function(m) log2(m + 1)

############################################################
# 5. Differential accessibility analysis
############################################################

design <- model.matrix(~ subject + W_1 + condition, data = pData(set_ruvr_nPeaks_k3))

y_rpf <- DGEList(
  counts = counts(set_ruvr_nPeaks_k3),
  group  = set_ruvr_nPeaks_k3$condition
)

y_rpf <- calcNormFactors(y_rpf, method = "TMM")
y_rpf <- estimateGLMCommonDisp(y_rpf, design)
y_rpf <- estimateGLMTagwiseDisp(y_rpf, design)

fit_rpf <- glmFit(y_rpf, design)
lrt_rpf <- glmLRT(fit_rpf)

summary(decideTests(lrt_rpf))

lrt_rpf_table <- as.data.frame(lrt_rpf[["table"]])

dar <- topTags(
  lrt_rpf,
  n = Inf,
  adjust.method = "BH",
  sort.by = "PValue",
  p.value = 0.05
)

dar_df <- as.data.frame(dar$table)
dar_df$peak_id <- rownames(dar_df)

write_xlsx(dar_df, path = dar_table_file)

############################################################
# 6. Annotate DARs by genomic region
############################################################

dar_output <- annotate_atac_peaks_regions(dar_df, nPeaks)

peak_annotation_df <- na.omit(as.data.frame(dar_output$dar_df_annotated))

peak_annotation_df <- peak_annotation_df %>%
  mutate(
    direction = case_when(
      logFC > 0 ~ "Up",
      logFC < 0 ~ "Down",
      TRUE ~ NA_character_
    )
  )

write_xlsx(peak_annotation_df, path = dar_annotated_file)

norm_counts <- set_ruvr_nPeaks_k2@assayData[["normalizedCounts"]]

############################################################
# 7. Summarize DAR classes
############################################################

dar_summary <- peak_annotation_df %>%
  filter(!is.na(direction)) %>%
  group_by(region_type, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(region_type) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  ungroup()

write_xlsx(dar_summary, path = dar_summary_file)

############################################################
# 8. Plot DAR composition by region and direction
############################################################

dar_annot <- peak_annotation_df

p_region <- dar_annot %>%
  filter(!is.na(direction)) %>%
  count(region_type, direction) %>%
  group_by(region_type) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ggplot(aes(x = region_type, y = n, fill = direction)) +
  geom_col(position = "fill", alpha = 0.3) +
  geom_text(
    aes(label = paste0(round(pct, 1), "%")),
    position = position_fill(vjust = 0.5),
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(
    values = c("Up" = "red3", "Down" = "royalblue3"),
    name = "Regulation"
  ) +
  labs(
    x = NULL,
    y = "% inside genomic region",
    fill = "Regulation",
    title = "DARs by genomic region and regulation"
  ) +
  theme_classic()

ggsave(
  filename = region_plot_file,
  plot = p_region,
  device = "svg",
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)

############################################################
# 9. Plot distance to TSS
############################################################

p_distance <- dar_annot %>%
  mutate(abs_dist_kb = abs(distanceToTSS) / 1000) %>%
  ggplot(aes(x = abs_dist_kb, color = direction)) +
  stat_ecdf(linewidth = 0.8) +
  coord_cartesian(xlim = c(0, 200)) +
  labs(
    x = "|Distance to TSS| (kb)",
    y = "CDF",
    title = "Distribution of distances to TSS"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = distance_plot_file,
  plot = p_distance,
  device = "tiff",
  width = 5,
  height = 4,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

############################################################
# 10. Volcano plot by genomic region
############################################################

p_volcano <- ggplot(dar_annot, aes(x = logFC, y = -log10(FDR), color = direction)) +
  geom_point(alpha = 0.6, size = 1) +
  facet_wrap(~ region_type, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(
    title = "Volcano by genomic region",
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = volcano_region_file,
  plot = p_volcano,
  device = "tiff",
  width = 6,
  height = 4,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

############################################################
# 11. PCA of normalized counts
############################################################

pca <- prcomp(t(log2p1(norm_counts)), center = TRUE, scale. = FALSE)

ev <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 2)

pca_df <- transform(
  as.data.frame(pca$x[, 1:2]),
  sample    = rownames(pca$x),
  subject   = factor(metadata$subject),
  condition = factor(metadata$condition, levels = c("Pre-Abx", "Post-Abx"))
)

group_colors <- c(
  "Pre-Abx"  = "#804D80CC",
  "Post-Abx" = "#804D80CC"
)

shape_map <- c(
  "Pre-Abx"  = 1,
  "Post-Abx" = 19
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = condition, label = subject)) +
  geom_point(size = 6) +
  geom_text(vjust = 2, size = 2.5) +
  scale_color_manual(values = group_colors, breaks = levels(pca_df$condition), name = NULL) +
  scale_shape_manual(values = shape_map, breaks = levels(pca_df$condition), name = NULL) +
  labs(
    title = "PCA of Normalized Counts",
    x = sprintf("PC1 (%.2f%%)", ev[1]),
    y = sprintf("PC2 (%.2f%%)", ev[2])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size = 10),
    axis.text        = element_text(size = 10)
  )

p_pca

############################################################
# 12. Save objects
############################################################

save(peak_annotation_df, file = dar_rdata_file)
save.image(file = workspace_file)



############################################################
# End of script
############################################################