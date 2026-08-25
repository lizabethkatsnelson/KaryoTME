library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsignif)

# Load data
data_table <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")
working_table <- bind_rows(data_table)

##### UCEC.MSI #####
UCEC.MSI_is_df <- data_table$UCEC.MSI_IS
UCEC.MSI_armcnv_df <- data_table$UCEC.MSI_ArmCNV
common_patients <- intersect(rownames(UCEC.MSI_is_df), rownames(UCEC.MSI_armcnv_df))
UCEC.MSI_is_df <- UCEC.MSI_is_df[common_patients, , drop = FALSE]
UCEC.MSI_armcnv_df <- UCEC.MSI_armcnv_df[common_patients, , drop = FALSE]

arm <- "Arm1q"

cnv_status <- case_when(
  UCEC.MSI_armcnv_df[[arm]] > 0.2  ~ "Gain",
  UCEC.MSI_armcnv_df[[arm]] < -0.2 ~ "Loss",
  TRUE                             ~ "Neutral"
)

plot_df <- data.frame(
  Ranked_Sum = UCEC.MSI_is_df$Ranked_Sum,
  CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
)

comparisons <- list(
  c("Loss", "Neutral"),
  c("Neutral", "Gain")
)

p <- ggplot(plot_df, aes(x = CNV_Status, y = Ranked_Sum, fill = CNV_Status)) +
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
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "none"
  ) +
  labs(
    title = "UCEC.MSI – Arm 1q",
    x = "CNV_Status",
    y = "Immune Score"
  )

ggsave("UCEC.MSI_Arm1q_IS.png", plot = p, width = 4, height = 5, dpi = 300, bg = "white")


##### UCEC.MSS #####
UCEC.MSS_is_df <- data_table$UCEC.MSS_IS
UCEC.MSS_armcnv_df <- data_table$UCEC.MSS_ArmCNV
common_patients <- intersect(rownames(UCEC.MSS_is_df), rownames(UCEC.MSS_armcnv_df))
UCEC.MSS_is_df <- UCEC.MSS_is_df[common_patients, , drop = FALSE]
UCEC.MSS_armcnv_df <- UCEC.MSS_armcnv_df[common_patients, , drop = FALSE]

arms_to_plot <- c("Arm1q")

for (arm in arms_to_plot) {

  cnv_status <- case_when(
    UCEC.MSS_armcnv_df[[arm]] > 0.2  ~ "Gain",
    UCEC.MSS_armcnv_df[[arm]] < -0.2 ~ "Loss",
    TRUE                             ~ "Neutral"
  )

  plot_df <- data.frame(
    Ranked_Sum = UCEC.MSS_is_df$Ranked_Sum,
    CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
  )

  comparisons <- list(
    c("Loss", "Neutral"),
    c("Neutral", "Gain")
  )

  p <- ggplot(plot_df, aes(x = CNV_Status, y = Ranked_Sum, fill = CNV_Status)) +
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
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = "none"
    ) +
    labs(
      title = "UCEC.MSS – Arm 1q",
      x = "CNV_Status",
      y = "Immune Score"
    )

  file_base <- paste0("UCEC.MSS_", arm, "_IS")
  ggsave(paste0(file_base, ".png"), plot = p, width = 4, height = 5, dpi = 300, bg = "white")
  # ggsave(paste0(file_base, ".pdf"), plot = p, width = 4, height = 5, device = "pdf", bg = "white")
}