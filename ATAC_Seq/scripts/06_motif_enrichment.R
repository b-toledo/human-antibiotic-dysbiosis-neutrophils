############################################################
# Script: 06_motif_enrichment.R
# Purpose: Perform motif enrichment analysis across
#          ATAC-seq DARs stratified by genomic region.
# Input: annotated DARs and motif enrichment resources
# Output: 6_Motif-Enrichment/
# Author: Bruna Toledo
############################################################

############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggrepel)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
  library(writexl)
})

############################################################
# 2. Define input/output paths
############################################################

out_dir <- "6_Motif-Enrichment"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sig_motif_file      <- file.path(out_dir, "TF_motif_FC1.5_Gen-reg.xlsx")
sig_motif_deg_file  <- file.path(out_dir, "TF_motif_FC1.5_Gen-reg_DEG.xlsx")
motif_plot_svg_file <- file.path(out_dir, "TF_VP.svg")
motif_plot_tiff_file <- file.path(out_dir, "TF_VP.tiff")

############################################################
# 3. Check required objects
############################################################

required_objects <- c(
  "peak_annotation_df",
  "background_peaks",
  "motifs",
  "Hs",
  "run_motif_enrichment",
  "motif_df_filt",
  "TF_DEGs"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "The following objects are missing from the environment: ",
    paste(missing_objects, collapse = ", "),
    ". Please load the motif enrichment workspace before running this script."
  )
}

############################################################
# 4. Define genomic regions to test
############################################################

region_list <- unique(na.omit(peak_annotation_df$region_type))

############################################################
# 5. Run motif enrichment by region and direction
############################################################

region_enrich <- lapply(region_list, function(reg) {
  peaks_reg <- peak_annotation_df %>% filter(region_type == reg)
  
  gr_reg_up <- GRanges(
    seqnames = peaks_reg$seqnames[peaks_reg$direction == "Up"],
    ranges = IRanges(
      peaks_reg$start[peaks_reg$direction == "Up"],
      peaks_reg$end[peaks_reg$direction == "Up"]
    )
  )
  
  gr_reg_down <- GRanges(
    seqnames = peaks_reg$seqnames[peaks_reg$direction == "Down"],
    ranges = IRanges(
      peaks_reg$start[peaks_reg$direction == "Down"],
      peaks_reg$end[peaks_reg$direction == "Down"]
    )
  )
  
  up_enrich <- run_motif_enrichment(gr_reg_up, background_peaks, motifs, Hs)
  down_enrich <- run_motif_enrichment(gr_reg_down, background_peaks, motifs, Hs)
  
  up_enrich <- up_enrich %>%
    mutate(
      region_type = reg,
      direction = "Up",
      log2OR_signed = log2(odds_ratio)
    )
  
  down_enrich <- down_enrich %>%
    mutate(
      region_type = reg,
      direction = "Down",
      log2OR_signed = -log2(odds_ratio)
    )
  
  bind_rows(up_enrich, down_enrich)
}) %>%
  bind_rows()

############################################################
# 6. Filter significant motifs
############################################################

region_colors <- c(
  "Promoter" = "#d97b29",
  "Distal"   = "#5ca972",
  "Trans"    = "#7a68b3"
)

sig_tf <- region_enrich %>%
  filter(adj_p < 0.05, odds_ratio > 1.5)

write_xlsx(sig_tf, path = sig_motif_file)

############################################################
# 7. Restrict to TFs that are also DEGs in RNA-seq
############################################################

tf_keep <- intersect(motif_df_filt$TF, TF_DEGs$gene_id)

df <- sig_tf %>%
  filter(TF %in% tf_keep)

write_xlsx(df, path = sig_motif_deg_file)

############################################################
# 8. Generate motif enrichment plot
############################################################

p_motif <- ggplot(df, aes(x = log2OR_signed, y = -log10(adj_p), color = region_type)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = df,
    aes(label = TF, color = region_type),
    size = 2,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.4,
    segment.size = 0.05
  ) +
  scale_color_manual(values = region_colors, name = "") +
  labs(
    x = "Mean log2(FC)",
    y = expression(-log[10] * "(FDR)"),
    title = "TF motif enrichment, DARs - Post-Abx vs Pre-Abx"
  ) +
  theme_bw(base_size = 6) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 6, hjust = 0.5, family = "Arial"),
    axis.title = element_text(size = 6, family = "Arial"),
    axis.text = element_text(size = 6, family = "Arial"),
    legend.title = element_text(size = 6, family = "Arial"),
    legend.text = element_text(size = 6, family = "Arial"),
    legend.position = "bottom"
  )

p_motif

############################################################
# 9. Save figures
############################################################

ggsave(
  motif_plot_svg_file,
  p_motif,
  width = 7,
  height = 6.5,
  units = "cm",
  dpi = 300,
  device = "svg"
)

ggsave(
  motif_plot_tiff_file,
  p_motif,
  width = 7,
  height = 6.5,
  units = "cm",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

############################################################
# 10. Save workspace
############################################################

save.image(file = "Motif.RData")



############################################################
# End of script
############################################################