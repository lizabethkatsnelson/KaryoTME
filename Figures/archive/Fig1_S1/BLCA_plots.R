library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsignif)

setwd("/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/BLCA")
data_table <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")
working_table <- bind_rows(data_table)

############ BLCA #########
BLCA_tis_df <- data_table$BLCA_TIS
BLCA_is_df <- data_table$BLCA_IS
BLCA_armcnv_df <- data_table$BLCA_ArmCNV

# Common patients
common_patients <- Reduce(intersect, list(
  rownames(BLCA_tis_df),
  rownames(BLCA_is_df),
  rownames(BLCA_armcnv_df)
))

# Subset
BLCA_tis_df <- BLCA_tis_df[common_patients, , drop = FALSE]
BLCA_is_df <- BLCA_is_df[common_patients, , drop = FALSE]
BLCA_armcnv_df <- BLCA_armcnv_df[common_patients, , drop = FALSE]

# Normalize IS scores
BLCA_is_df$IS_rank_norm <- (rank(BLCA_is_df$Ranked_Sum) - 1) / (length(BLCA_is_df$Ranked_Sum) - 1)
BLCA_is_df$IS_minmax <- (BLCA_is_df$Ranked_Sum - min(BLCA_is_df$Ranked_Sum)) /
                        (max(BLCA_is_df$Ranked_Sum) - min(BLCA_is_df$Ranked_Sum))

# Arms to plot
arms_to_plot <- c("Arm9p", "Arm11p", "Arm17q")

for (arm in arms_to_plot) {

  if (!(arm %in% colnames(BLCA_armcnv_df))) {
    message(paste("Skipping", arm, "- not found in BLCA_armcnv_df"))
    next
  }

  cnv_status <- case_when(
    BLCA_armcnv_df[[arm]] > 0.2  ~ "Gain",
    BLCA_armcnv_df[[arm]] < -0.2 ~ "Loss",
    TRUE                         ~ "Neutral"
  )

  comparisons <- list(c("Loss", "Neutral"), c("Neutral", "Gain"))

  ## TIS Plot
  plot_df_tis <- data.frame(
    TIS = BLCA_tis_df$TIS,
    CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
  )

  p_tis <- ggplot(plot_df_tis, aes(x = CNV_Status, y = TIS, fill = CNV_Status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, color = "black") +
    geom_signif(comparisons = comparisons, test = "wilcox.test", map_signif_level = TRUE, step_increase = 0.1) +
    scale_fill_manual(values = c("Loss"="blue", "Neutral"="gray", "Gain"="red")) +
    theme_minimal() +
    labs(title = paste("BLCA:", arm, "- TIS Score by CNV Status"), x = "CNV Status", y = "TIS Score")
  ggsave(paste0("BLCA_", arm, "_TIS_Boxplot_pairwiseStars.png"), p_tis, width = 5, height = 6, dpi = 300)

  ## Rank-normalized IS Plot
  plot_df_rank <- data.frame(
    IS = BLCA_is_df$IS_rank_norm,
    CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
  )

  p_rank <- ggplot(plot_df_rank, aes(x = CNV_Status, y = IS, fill = CNV_Status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, color = "black") +
    geom_signif(comparisons = comparisons, test = "wilcox.test", map_signif_level = TRUE, step_increase = 0.1) +
    scale_fill_manual(values = c("Loss"="blue", "Neutral"="gray", "Gain"="red")) +
    theme_minimal() +
    labs(title = paste("BLCA:", arm, "- Rank-Normalized IS by CNV"), x = "CNV Status", y = "Rank-Normalized IS")
  ggsave(paste0("BLCA_", arm, "_ISrankNorm_Boxplot_pairwiseStars.png"), p_rank, width = 5, height = 6, dpi = 300)

  ## Min-max normalized IS Plot
  plot_df_minmax <- data.frame(
    IS = BLCA_is_df$IS_minmax,
    CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
  )

  p_minmax <- ggplot(plot_df_minmax, aes(x = CNV_Status, y = IS, fill = CNV_Status)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, color = "black") +
    geom_signif(comparisons = comparisons, test = "wilcox.test", map_signif_level = TRUE, step_increase = 0.1) +
    scale_fill_manual(values = c("Loss"="blue", "Neutral"="gray", "Gain"="red")) +
    theme_minimal() +
    labs(title = paste("BLCA:", arm, "- Min-Max Normalized IS by CNV"), x = "CNV Status", y = "IS (Min-Max Normalized)")
  ggsave(paste0("BLCA_", arm, "_ISminmaxNorm_Boxplot_pairwiseStars.png"), p_minmax, width = 5, height = 6, dpi = 300)
}