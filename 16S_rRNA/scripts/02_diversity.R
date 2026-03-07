############################################################
# Script: 02_diversity.R
# Purpose: Perform rarefaction and calculate alpha and beta
#          diversity metrics to assess microbiota changes
#          following antibiotic treatment.
# Input: curated MPSE object
# Output: rarefaction plots, alpha and beta diversity plots,
#         and summary tables exported to .xlsx
# Author: Bruna Toledo
############################################################

## =========================================================
## 1. Load libraries
## =========================================================

suppressPackageStartupMessages({
  library(MicrobiotaProcess)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(patchwork)
  library(openxlsx)
})

## =========================================================
## 2. Define input and output paths
## =========================================================

input_file <- "16S_rRNA/results/intermediate/mpse_curated.rds"

out_dir_intermediate <- "16S_rRNA/results/intermediate"
out_dir_tables       <- "16S_rRNA/results/tables"
out_dir_figures      <- "16S_rRNA/results/figures"

if (!dir.exists(out_dir_intermediate)) dir.create(out_dir_intermediate, recursive = TRUE)
if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)
if (!dir.exists(out_dir_figures)) dir.create(out_dir_figures, recursive = TRUE)

## =========================================================
## 3. Load curated MPSE object
## =========================================================

mpse <- readRDS(input_file)

message("Curated MPSE object loaded successfully.")

## =========================================================
## 4. Rarefaction
## =========================================================

mpse <- mp_rrarefy(mpse)
message("Rarefaction completed successfully.")

mpse <- mp_cal_rarecurve(
  mpse,
  .abundance = RareAbundance,
  chunks = 1000
)

message("Rarefaction curves calculated successfully.")

## =========================================================
## 5. Plot rarefaction curves
## =========================================================

condition_colors_fill <- c(
  "Pre-Abx"  = "#3C548866",
  "Post-Abx" = "#87255B66"
)

condition_colors_line <- c(
  "Pre-Abx"  = "#3C548866",
  "Post-Abx" = "#87255B66"
)

p_rarefaction_samples <- mpse %>%
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve,
    .alpha = Observe
  ) +
  theme_minimal() +
  scale_fill_manual(values = condition_colors_fill) +
  scale_color_manual(values = condition_colors_line)

p_rarefaction_group <- mpse %>%
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve,
    .alpha = Observe,
    .group = condition
  ) +
  theme_minimal() +
  scale_fill_manual(values = condition_colors_fill) +
  scale_color_manual(values = condition_colors_line)

p_rarefaction_group_mean <- mpse %>%
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve,
    .alpha = "Observe",
    .group = condition,
    plot.group = TRUE
  ) +
  theme_minimal() +
  scale_fill_manual(values = condition_colors_fill) +
  scale_color_manual(values = condition_colors_line)

p_rarefaction_combined <- p_rarefaction_group + p_rarefaction_group_mean

ggsave(
  filename = file.path(out_dir_figures, "rarefaction_curves.pdf"),
  plot = p_rarefaction_combined,
  width = 10,
  height = 5
)

## =========================================================
## 6. Export abundance tables
## =========================================================

ASV_Ab <- as.data.frame(mpse@assays@data@listData[["Abundance"]]) %>%
  rownames_to_column(var = "ASV")

ASV_RAb <- as.data.frame(mpse@assays@data@listData[["RareAbundance"]]) %>%
  rownames_to_column(var = "ASV")

wb_asv <- createWorkbook()

addWorksheet(wb_asv, "Unrarefied-Abundance")
writeData(wb_asv, "Unrarefied-Abundance", ASV_Ab)

addWorksheet(wb_asv, "Rarefied-Abundance")
writeData(wb_asv, "Rarefied-Abundance", ASV_RAb)

asv_outfile <- file.path(out_dir_tables, "ASV_tables.xlsx")
saveWorkbook(wb_asv, asv_outfile, overwrite = TRUE)

message("ASV tables exported to: ", asv_outfile)

## =========================================================
## 7. Alpha diversity
## =========================================================

mpse <- mp_cal_alpha(mpse, .abundance = RareAbundance)
message("Alpha diversity metrics calculated successfully.")

p_alpha <- mpse %>%
  mp_plot_alpha(
    .group = condition,
    .alpha = c(Observe, Shannon, Pielou),
    test = "wilcox.test",
    test.args = list(paired = TRUE),
    textsize = 3.5,
    comparisons = list(c("Pre-Abx", "Post-Abx"))
  ) +
  scale_fill_manual(values = c("#3C548866", "#87255B66"), guide = "none") +
  scale_color_manual(values = c("#3C548866", "#87255B66"), guide = "none")

ggsave(
  filename = file.path(out_dir_figures, "alpha_diversity.pdf"),
  plot = p_alpha,
  width = 9,
  height = 4.5
)

alpha_df <- tibble(
  SampleID   = rownames(colData(mpse)),
  Subject    = mpse@colData@listData[["subject"]],
  condition  = mpse@colData@listData[["condition"]],
  Observe    = mpse@colData@listData[["Observe"]],
  Chao1      = mpse@colData@listData[["Chao1"]],
  ACE        = mpse@colData@listData[["ACE"]],
  Shannon    = mpse@colData@listData[["Shannon"]],
  Simpson    = mpse@colData@listData[["Simpson"]],
  Pielou     = mpse@colData@listData[["Pielou"]]
)

## =========================================================
## 8. Beta diversity
## =========================================================

mpse <- mp_decostand(mpse, .abundance = Abundance)
mpse <- mp_cal_dist(mpse, .abundance = hellinger, distmethod = "bray")
mpse <- mp_cal_pcoa(mpse, .abundance = hellinger, distmethod = "bray")

message("Beta diversity and PCoA calculated successfully.")

mpse <- mp_adonis(
  mpse,
  .abundance = hellinger,
  .formula = ~subject + condition,
  distmethod = "bray",
  permutations = 9999,
  action = "add"
)

adonis <- mp_extract_internal_attr(mpse, name = "adonis")

p_beta <- mpse %>%
  mp_plot_ord(
    .ord = pcoa,
    .group = condition,
    .color = condition,
    .size = 4.8,
    .alpha = 0.75,
    ellipse = TRUE,
    show.legend = FALSE,
    show.adonis = TRUE
  ) +
  scale_fill_manual(values = c("#3C548866", "#87255B66")) +
  scale_color_manual(values = c("#3C548866", "#87255B66"))

ggsave(
  filename = file.path(out_dir_figures, "beta_diversity_pcoa_bray.pdf"),
  plot = p_beta,
  width = 8,
  height = 6
)

## =========================================================
## 9. Extract Bray-Curtis distance matrix
## =========================================================

bray_list <- mpse@colData@listData[["bray"]]
sample_names <- mpse@colData@rownames
names(bray_list) <- sample_names

bray_df <- purrr::imap_dfr(bray_list, ~ mutate(.x, From = .y))
colnames(bray_df) <- c("Sample", "BrayCurtis", "From")

bray_df_clean <- bray_df %>%
  rename(Sample1 = From, Sample2 = Sample)

bray_df_sym <- bray_df_clean %>%
  bind_rows(
    bray_df_clean %>%
      rename(Sample1 = Sample2, Sample2 = Sample1)
  ) %>%
  distinct()

bray_matrix <- bray_df_sym %>%
  pivot_wider(names_from = Sample2, values_from = BrayCurtis) %>%
  column_to_rownames("Sample1") %>%
  as.data.frame()

bray_df_tb <- as.data.frame(bray_matrix)

## =========================================================
## 10. Save diversity metrics to workbook
## =========================================================

wb_div <- createWorkbook()

addWorksheet(wb_div, "Alpha-Diversity")
writeData(wb_div, "Alpha-Diversity", alpha_df)

addWorksheet(wb_div, "Bray-Curtis")
writeData(wb_div, "Bray-Curtis", bray_df_tb)

addWorksheet(wb_div, "Adonis_Permanova")
writeData(wb_div, "Adonis_Permanova", adonis)

div_outfile <- file.path(out_dir_tables, "microbiome_diversity.xlsx")
saveWorkbook(wb_div, div_outfile, overwrite = TRUE)

message("Diversity tables exported to: ", div_outfile)

## =========================================================
## 11. Save processed object for downstream analyses
## =========================================================

saveRDS(mpse, file.path(out_dir_intermediate, "mpse_diversity.rds"))

message("Processed MPSE object saved to: ",
        file.path(out_dir_intermediate, "mpse_diversity.rds"))