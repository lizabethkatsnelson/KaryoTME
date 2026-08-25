library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsignif)

setwd("/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/BLCA")
data_table <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")
working_table <- bind_rows(data_table)

############ BLCA #########
# Prepare data frame 
BLCA_is_df <- data_table$BLCA_IS
BLCA_armcnv_df <- data_table$BLCA_ArmCNV
common_patients <- intersect(rownames(BLCA_is_df), rownames(BLCA_armcnv_df))
BLCA_is_df <- BLCA_is_df[common_patients, , drop = FALSE]
BLCA_armcnv_df <- BLCA_armcnv_df[common_patients, , drop = FALSE]

# Arms to plot
arms_to_plot <- c("Arm9p", "Arm11p", "Arm17q")

for (arm in arms_to_plot) {

  if (!(arm %in% colnames(BLCA_armcnv_df))) {
    message(paste("Skipping", arm, "- not found in BLCA_ArmCNV"))
    next
  }

  # Assign CNV status
  cnv_status <- case_when(
    BLCA_armcnv_df[[arm]] > 0.2  ~ "Gain",
    BLCA_armcnv_df[[arm]] < -0.2 ~ "Loss",
    TRUE                         ~ "Neutral"
  )

  # Prepare plot data
  plot_df <- data.frame(
    IS = BLCA_is_df$Ranked_Sum,
    CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
  )

  # Define comparisons
  comparisons <- list(
    c("Loss", "Neutral"),
    c("Neutral", "Gain")
  )

  # Create the plot
  p <- ggplot(plot_df, aes(x = CNV_Status, y = IS, fill = CNV_Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
    geom_signif(
      comparisons = comparisons,
      test = "wilcox.test",
      map_signif_level = TRUE,
      step_increase = 0.1,
      textsize = 5
    ) +
    scale_fill_manual(values = c("Loss" = "#4F6DB8", "Neutral" = "#B8B8B8", "Gain" = "#EA6A6A")) +
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    labs(
      x = NULL,
      y = "Immune Score"
    )

  # Save plots
  file_base <- paste0("BLCA_", arm, "_IS_RankedSum")
  ggsave(paste0(file_base, ".png"), plot = p, width = 4, height = 5, dpi = 300)
  ggsave(paste0(file_base, ".pdf"), plot = p, width = 4, height = 5, device = "pdf")
}