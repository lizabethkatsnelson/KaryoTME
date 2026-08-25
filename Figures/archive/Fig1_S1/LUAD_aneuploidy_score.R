library(tidyr)
library(dplyr)
library(ggplot2)

setwd("/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/LUAD")
data_table <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")
working_table <- bind_rows(data_table)

# Load data
LUAD_is_df <- data_table$LUAD_IS
LUAD_armcnv_df <- data_table$LUAD_ArmCNV
common_patients <- intersect(rownames(LUAD_is_df), rownames(LUAD_armcnv_df))
LUAD_is_df <- LUAD_is_df[common_patients, , drop = FALSE]
LUAD_armcnv_df <- LUAD_armcnv_df[common_patients, , drop = FALSE]

# Aneuploidy score
LUAD_is_df$Aneuploidy_Score_Abs <- rowSums(abs(LUAD_armcnv_df), na.rm = TRUE)

# Arms to define 
arms_to_plot <- c("Arm1q", "Arm5q", "Arm6p", "Arm9p", "Arm6q")

###  Linear Model with Aneuploidy as Covariate
# Prepare full dataset (no filtering)
cnv_status <- case_when(
  LUAD_armcnv_df[[arm]] > 0.2  ~ "Gain",
  LUAD_armcnv_df[[arm]] < -0.2 ~ "Loss",
  TRUE                         ~ "Neutral"
)

plot_df_lm <- data.frame(
  IS = LUAD_is_df$Ranked_Sum,
  CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain")),
  Aneuploidy = LUAD_is_df$Aneuploidy_Score_Abs
)

# Fit linear model to remove aneuploidy
plot_df_lm$CNV_Status <- factor(plot_df_lm$CNV_Status, levels = c("Neutral", "Loss", "Gain"))


model <- lm(IS ~ CNV_Status + Aneuploidy, data = plot_df_lm)
summary(model)  # View p-values

# Extract p-value for CNV_StatusGain and CNV_StatusNeutral
p_loss <- summary(model)$coefficients["CNV_StatusLoss", "Pr(>|t|)"]
p_gain <- summary(model)$coefficients["CNV_StatusGain", "Pr(>|t|)"]
# Annotated boxplot with adjusted p-values
ggplot(plot_df_lm, aes(x = CNV_Status, y = IS, fill = CNV_Status)) +
  geom_boxplot(outlier.shape = NA) +
  annotate("text", x = 2, y = max(plot_df_lm$IS), label = paste0("p = ", signif(p_neutral, 3)), size = 4) +
  annotate("text", x = 3, y = max(plot_df_lm$IS)*0.95, label = paste0("p = ", signif(p_gain, 3)), size = 4) +
  scale_fill_manual(values = c("Loss" = "#4F6DB8", "Neutral" = "#B8B8B8", "Gain" = "#EA6A6A")) +
  theme_minimal(base_size = 16) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = NULL, y = "Immune Score", title = "Linear Model (adjusted for Aneuploidy)")




# Filter top 25% most aneuploid
top25_cutoff <- quantile(LUAD_is_df$Aneuploidy_Score_Abs, probs = 0.75, na.rm = TRUE)
keep <- LUAD_is_df$Aneuploidy_Score_Abs <= top25_cutoff
filtered_df <- LUAD_is_df[keep, , drop = FALSE]
filtered_cnv <- LUAD_armcnv_df[rownames(filtered_df), , drop = FALSE]

# Define CNV status
arm <- "Arm1q"
cnv_status <- case_when(
  filtered_cnv[[arm]] > 0.2  ~ "Gain",
  filtered_cnv[[arm]] < -0.2 ~ "Loss",
  TRUE                       ~ "Neutral"
)

# Combine into plot_df
plot_df <- data.frame(
  IS = filtered_df$Ranked_Sum,
  CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
)

# Manually compute Wilcoxon p-values
p_neutral_vs_gain <- wilcox.test(
  IS ~ CNV_Status,
  data = plot_df %>% filter(CNV_Status %in% c("Neutral", "Gain"))
)$p.value

p_neutral_vs_loss <- wilcox.test(
  IS ~ CNV_Status,
  data = plot_df %>% filter(CNV_Status %in% c("Neutral", "Loss"))
)$p.value

# Format function
format_p <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("NS.")
}

# Plot with both annotations
# Plot with * or NS only (no numeric p-values)
p <- ggplot(plot_df, aes(x = CNV_Status, y = IS, fill = CNV_Status)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85) +
  annotate("text", x = 2.5, y = max(plot_df$IS) * 0.95,
           label = format_p(p_neutral_vs_gain),
           size = 6, fontface = "bold") +
  annotate("text", x = 1.5, y = max(plot_df$IS) * 1.05,
           label = format_p(p_neutral_vs_loss),
           size = 6, fontface = "bold") +
  scale_fill_manual(values = c("Loss" = "#4F6DB8", "Neutral" = "#B8B8B8", "Gain" = "#EA6A6A")) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Immune Score (Ranked Sum)"
  )
  
# Save
ggsave("LUAD_Arm1q_Top25Filtered_WilcoxonManual.pdf", plot = p, width = 4, height = 5, dpi = 300)