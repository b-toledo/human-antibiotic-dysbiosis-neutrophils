############################################################
# Script: 03_pca_visualization.R
# Purpose: Generate PCA plots from normalized ATAC-seq
#          accessibility counts.
# Input: normalized counts and sample metadata
# Output: 3_PCA/PCA.svg
# Author: Bruna Toledo
############################################################

############################################################
# 1. Load libraries
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

############################################################
# 2. Define input/output paths
############################################################

out_dir <- "3_PCA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pca_file <- file.path(out_dir, "PCA2.svg")

############################################################
# 3. Check required objects
############################################################

required_objects <- c("norm_counts", "metadata")

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "The following objects are missing from the environment: ",
    paste(missing_objects, collapse = ", "),
    ". Please load the ATAC normalized count workspace before running this script."
  )
}

############################################################
# 4. Define helper functions
############################################################

log2p1 <- function(m) log2(m + 1)

############################################################
# 5. Run PCA
############################################################

pca <- prcomp(
  t(log2p1(norm_counts)),
  center = TRUE,
  scale. = FALSE
)

############################################################
# 6. Prepare plotting data
############################################################

ev <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 2)

pca_df <- transform(
  as.data.frame(pca$x[, 1:2]),
  sample    = rownames(pca$x),
  subject   = factor(metadata$subject),
  condition = factor(metadata$condition, levels = c("Pre-Abx", "Post-Abx"))
)

############################################################
# 7. Define plot aesthetics
############################################################

group_colors <- c(
  "Pre-Abx"  = "#804D80CC",
  "Post-Abx" = "#804D80CC"
)

shape_map <- c(
  "Pre-Abx"  = 1,
  "Post-Abx" = 19
)

############################################################
# 8. Generate PCA plot
############################################################

p_pca <- ggplot(
  pca_df,
  aes(PC1, PC2, color = condition, shape = condition, label = subject)
) +
  geom_point(size = 3) +
  geom_text(vjust = 2, size = 1) +
  scale_color_manual(
    values = group_colors,
    breaks = levels(pca_df$condition),
    name = NULL
  ) +
  scale_shape_manual(
    values = shape_map,
    breaks = levels(pca_df$condition),
    name = NULL
  ) +
  labs(
    title = "ATAC-seq - PCA of Normalized Counts",
    x = sprintf("PC1 (%.2f%%)", ev[1]),
    y = sprintf("PC2 (%.2f%%)", ev[2])
  ) +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 6),
    axis.text        = element_text(size = 6),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.2),
    axis.line.y.left   = element_line(color = "black", linewidth = 0.2),
    axis.ticks         = element_line(color = "black")
  )

p_pca

############################################################
# 9. Save figure
############################################################

ggsave(
  filename = pca_file,
  plot = p_pca,
  width = 6.5,
  height = 6,
  units = "cm"
)

############################################################
# 10. Save workspace
############################################################

save.image(file = "PCA.RData")



############################################################
# End of script
############################################################