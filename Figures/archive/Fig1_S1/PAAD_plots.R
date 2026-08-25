library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsignif)

setwd("/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/PAAD")
data_table <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")
working_table <- bind_rows(data_table)

############ PAAD #########
PAAD_tis_df <- data_table$PAAD_TIS
PAAD_is_df <- data_table$PAAD_IS
PAAD_armcnv_df <- data_table$PAAD_ArmCNV

# Common patients
common_patients <- Reduce(intersect, list(
  rownames(PAAD_tis_df),
  rownames(PAAD_is_df),
  rownames(PAAD_armcnv_df)
))

# Subset
PAAD_tis_df <- PAAD_tis_df[common_patients, , drop = FALSE]
PAAD_is_df <- PAAD_is_df[common_patients, , drop = FALSE]
PAAD_armcnv_df <- PAAD_armcnv_df[common_patients, , drop = FALSE]

# Normalize IS scores
PAAD_is_df$IS_rank_norm <- (rank(PAAD_is_df$Ranked_Sum) - 1) / (length(PAAD_is_df$Ranked_Sum) - 1)
PAAD_is_df$IS_minmax <- (PAAD_is_df$Ranked_Sum - min(PAAD_is_df$Ranked_Sum)) /
                        (max(PAAD_is_df$Ranked_Sum) - min(PAAD_is_df$Ranked_Sum))

# Arms to plot
arms_to_plot <- c("Arm1q", "Arm6q", "Arm9p")

for (arm in arms_to_plot) {

  if (!(arm %in% colnames(PAAD_armcnv_df))) {
    message(paste("Skipping", arm, "- not found in PAAD_ArmCNV"))
    next
  }

  cnv_status <- case_when(
    PAAD_armcnv_df[[arm]] > 0.2  ~ "Gain",
    PAAD_armcnv_df[[arm]] < -0.2 ~ "Loss",
    TRUE                         ~ "Neutral"
  )

  comparisons <- list(c("Loss", "Neutral"), c("Neutral", "Gain"))

  ## TIS Plot
  plot_df_tis <- data.frame(
    TIS = PAAD_tis_df$TIS,
    CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
  )

  p_tis <- ggplot(plot_df_tis, aes(x = CNV_Status, y = TIS, fill = CNV_Status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, color = "black") +
    geom_signif(comparisons = comparisons, test = "wilcox.test", map_signif_level = TRUE, step_increase = 0.1) +
    scale_fill_manual(values = c("Loss"="blue", "Neutral"="gray", "Gain"="red")) +
    theme_minimal() +
    labs(title = paste("PAAD:", arm), x = "CNV Status", y = "TIS Score")
  ggsave(paste0("PAAD_", arm, "_TIS_Boxplot_pairwiseStars.png"), p_tis, width = 5, height = 6, dpi = 300)

  ## Rank-normalized IS Plot
  plot_df_rank <- data.frame(
    IS = PAAD_is_df$IS_rank_norm,
    CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
  )

  p_rank <- ggplot(plot_df_rank, aes(x = CNV_Status, y = IS, fill = CNV_Status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, color = "black") +
    geom_signif(comparisons = comparisons, test = "wilcox.test", map_signif_level = TRUE, step_increase = 0.1) +
    scale_fill_manual(values = c("Loss"="blue", "Neutral"="gray", "Gain"="red")) +
    theme_minimal() +
    labs(title = paste("PAAD:", arm), x = "CNV Status", y = "Rank-Normalized IS")
  ggsave(paste0("PAAD_", arm, "_ISrankNorm_Boxplot_pairwiseStars.png"), p_rank, width = 5, height = 6, dpi = 300)

  ## Min-max normalized IS Plot
  plot_df_minmax <- data.frame(
    IS = PAAD_is_df$IS_minmax,
    CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
  )

  p_minmax <- ggplot(plot_df_minmax, aes(x = CNV_Status, y = IS, fill = CNV_Status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, color = "black") +
    geom_signif(comparisons = comparisons, test = "wilcox.test", map_signif_level = TRUE, step_increase = 0.1) +
    scale_fill_manual(values = c("Loss"="blue", "Neutral"="gray", "Gain"="red")) +
    theme_minimal() +
    labs(title = paste("PAAD:", arm), x = "CNV Status", y = "IS (Min-Max Normalized)")
  ggsave(paste0("PAAD_", arm, "_ISminmaxNorm_Boxplot_pairwiseStars.png"), p_minmax, width = 5, height = 6, dpi = 300)
}