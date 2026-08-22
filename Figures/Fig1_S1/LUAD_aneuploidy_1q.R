library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsignif)

# Load and preprocess
data_table <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")
LUAD_is_df <- data_table$LUAD_IS
LUAD_armcnv_df <- data_table$LUAD_ArmCNV
common_patients <- intersect(rownames(LUAD_is_df), rownames(LUAD_armcnv_df))
LUAD_is_df <- LUAD_is_df[common_patients, , drop = FALSE]
LUAD_armcnv_df <- LUAD_armcnv_df[common_patients, , drop = FALSE]

# Compute aneuploidy score
LUAD_is_df$Aneuploidy_Score_Abs <- rowSums(abs(LUAD_armcnv_df), na.rm = TRUE)

# Filter for bottom 75%
cutoff <- quantile(LUAD_is_df$Aneuploidy_Score_Abs, 0.75)
filtered_df <- LUAD_is_df[LUAD_is_df$Aneuploidy_Score_Abs <= cutoff, , drop = FALSE]
filtered_cnv <- LUAD_armcnv_df[rownames(filtered_df), , drop = FALSE]

# Define CNV status for Arm1q
arm <- "Arm1q"
cnv_status <- case_when(
  filtered_cnv[[arm]] > 0.2  ~ "Gain",
  filtered_cnv[[arm]] < -0.2 ~ "Loss",
  TRUE ~ "Neutral"
)

# Combine
plot_df <- data.frame(
  IS = filtered_df$Ranked_Sum,
  CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
)

# Manual Wilcoxon p-values
p_ng <- wilcox.test(IS ~ CNV_Status, data = plot_df %>% filter(CNV_Status %in% c("Neutral", "Gain")))$p.value
p_nl <- wilcox.test(IS ~ CNV_Status, data = plot_df %>% filter(CNV_Status %in% c("Neutral", "Loss")))$p.value

# Helper to convert to significance stars
format_p <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("NS.")
}

# Draw plot with line-based significance using geom_signif
p <- ggplot(plot_df, aes(x = CNV_Status, y = IS, fill = CNV_Status)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.85) +
  geom_signif(comparisons = list(c("Neutral", "Loss")),
              annotations = format_p(p_nl),
              y_position = max(plot_df$IS) * 1.05,
              tip_length = 0.01, textsize = 6, vjust = 0.3) +
  geom_signif(comparisons = list(c("Neutral", "Gain")),
              annotations = format_p(p_ng),
              y_position = max(plot_df$IS) * 0.95,
              tip_length = 0.01, textsize = 6, vjust = 0.3) +
  scale_fill_manual(values = c("Loss" = "#4F6DB8", "Neutral" = "#B8B8B8", "Gain" = "#EA6A6A")) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "none",
    plot.title = element_blank()
  ) +
  labs(
    y = "Immune Score"
  )

# Save
ggsave("LUAD_Arm1q_Top25Filtered_1q.pdf", plot = p, width = 4, height = 5, dpi = 300)