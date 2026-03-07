############################################################
# Script: 07_motif_enrichment_deg_tfs_plot.R
# Purpose: Visualize motif enrichment restricted to
#          transcription factors that are DEGs in RNA-seq.
# Input: motif enrichment results and RNA-seq TF DEGs
# Output: 6_Motif-Enrichment/volcano_TF_DEGs_only.tiff
# Author: Bruna Toledo
############################################################

############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

############################################################
# 2. Define input/output paths
############################################################

out_dir <- "6_Motif-Enrichment"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

volcano_file <- file.path(out_dir, "volcano_TF_DEGs_only.tiff")

############################################################
# 3. Check required objects
############################################################

required_objects <- c("motif_df_filt", "TF_DEGs")

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "The following objects are missing from the environment: ",
    paste(missing_objects, collapse = ", "),
    ". Please load the motif enrichment and RNA-seq TF DEG objects before running this script."
  )
}

############################################################
# 4. Keep only TFs that are DE in RNA-seq
############################################################

tf_keep <- intersect(motif_df_filt$TF, TF_DEGs$gene_id)

df <- motif_df_filt %>%
  filter(TF %in% tf_keep) %>%
  mutate(
    log2OR_signed = ifelse(odds_ratio > 0, log2(odds_ratio), -log2(abs(odds_ratio))),
    neglog10FDR   = -log10(adj_p)
  )

############################################################
# 5. Select one motif per TF and direction
############################################################

best_tf <- df %>%
  group_by(DAR_direction, TF) %>%
  slice_min(adj_p, n = 1, with_ties = FALSE) %>%
  ungroup()

############################################################
# 6. Select top TFs for labeling
############################################################

top_tf <- best_tf %>%
  group_by(DAR_direction) %>%
  slice_min(adj_p, n = 13, with_ties = FALSE) %>%
  ungroup()

############################################################
# 7. Generate volcano plot
############################################################

p_volcano <- ggplot(best_tf, aes(x = log2OR_signed, y = neglog10FDR, color = DAR_direction)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
  geom_point(alpha = 0.7, size = 2) +
  ggrepel::geom_text_repel(
    data = top_tf,
    aes(label = TF),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    segment.size = 0.2
  ) +
  scale_color_manual(
    values = c("Open-Peaks" = "red", "Closed-Peaks" = "blue"),
    name = NULL
  ) +
  labs(
    x = "signed log2(odds ratio)",
    y = expression(-log[10] * "(FDR)"),
    title = "TF motifs enriched in DARs (DEG TFs only)"
  ) +
  theme_bw(base_size = 11)

p_volcano

############################################################
# 8. Save figure
############################################################

ggsave(
  volcano_file,
  p_volcano,
  width = 8,
  height = 4.5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

############################################################
# 9. Save workspace
############################################################

save.image(file = "Motif_plot.RData")



############################################################
# End of script
############################################################