library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsignif)

setwd("/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/LUAD")
data_table <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")
working_table <- bind_rows(data_table)

############ LUAD #########
# Prepare data frame 
LUAD_is_df <- data_table$LUAD_IS
LUAD_armcnv_df <- data_table$LUAD_ArmCNV
common_patients <- intersect(rownames(LUAD_is_df), rownames(LUAD_armcnv_df))
LUAD_is_df <- LUAD_is_df[common_patients, , drop = FALSE]
LUAD_armcnv_df <- LUAD_armcnv_df[common_patients, , drop = FALSE]

# Now calculate absolute aneuploidy score
absolute_aneuploidy_score <- rowSums(abs(LUAD_armcnv_df), na.rm = TRUE)

# Assign to IS data frame (now they match)
LUAD_is_df$Aneuploidy_Score_Abs <- absolute_aneuploidy_score

arm <- "Arm1q"
cnv_status <- case_when(
  LUAD_armcnv_df[[arm]] > 0.2  ~ "Gain",
  LUAD_armcnv_df[[arm]] < -0.2 ~ "Loss",
  TRUE                         ~ "Neutral"
)

plot_df_lm <- data.frame(
  IS = LUAD_is_df$Ranked_Sum,
  CNV_Status = factor(cnv_status, levels = c("Neutral", "Loss", "Gain")),
  Aneuploidy = LUAD_is_df$Aneuploidy_Score_Abs
)

model <- lm(IS ~ CNV_Status + Aneuploidy, data = plot_df_lm)
summary(model)

######## Remove Top 15–20% Aneuploidy Samples
# Step 1: Define cutoff for top 15% (adjust to 0.80 for top 20%)
quantile_val <- quantile(LUAD_is_df$Aneuploidy_Score_Abs, probs = 0.75, na.rm = TRUE)

# Step 2: Filter out top 15%
keep <- LUAD_is_df$Aneuploidy_Score_Abs <= quantile_val
filtered_df <- LUAD_is_df[keep, , drop = FALSE]
filtered_cnv <- LUAD_armcnv_df[rownames(filtered_df), , drop = FALSE]

# Step 3: Define CNV status (e.g., 1q gain)
cnv_status <- case_when(
  filtered_cnv[["Arm1q"]] > 0.2  ~ "Gain",
  filtered_cnv[["Arm1q"]] < -0.2 ~ "Loss",
  TRUE                           ~ "Neutral"
)

# Step 4: Model
plot_df <- data.frame(
  IS = filtered_df$Ranked_Sum,
  CNV_Status = factor(cnv_status, levels = c("Neutral", "Loss", "Gain"))
)
model <- lm(IS ~ CNV_Status, data = plot_df)
summary(model)


#### remove top 15% + consdier aneuploidy as co-variate
# Step 1: Define cutoff for top 25% (adjust to 0.75 for top 25% kept)
quantile_val <- quantile(LUAD_is_df$Aneuploidy_Score_Abs, probs = 0.85, na.rm = TRUE)

# Step 2: Filter out top 25%
keep <- LUAD_is_df$Aneuploidy_Score_Abs <= quantile_val
filtered_df <- LUAD_is_df[keep, , drop = FALSE]
filtered_cnv <- LUAD_armcnv_df[rownames(filtered_df), , drop = FALSE]

# Step 3: Define CNV status (e.g., 1q gain)
cnv_status <- case_when(
  filtered_cnv[["Arm1q"]] > 0.2  ~ "Gain",
  filtered_cnv[["Arm1q"]] < -0.2 ~ "Loss",
  TRUE                           ~ "Neutral"
)

# Step 4: Model input dataframe (including aneuploidy as covariate!)
plot_df <- data.frame(
  IS = filtered_df$Ranked_Sum,
  CNV_Status = factor(cnv_status, levels = c("Neutral", "Loss", "Gain")),
  Aneuploidy = filtered_df$Aneuploidy_Score_Abs
)

# Step 5: Fit model with aneuploidy as covariate
model <- lm(IS ~ CNV_Status + Aneuploidy, data = plot_df)
summary(model)


#### Plot 

# Arms to plot
arms_to_plot <- c("Arm1q", "Arm5q", "Arm6p", "Arm9p", "Arm6q")

for (arm in arms_to_plot) {

  if (!(arm %in% colnames(LUAD_armcnv_df))) {
    message(paste("Skipping", arm, "- not found in LUAD_ArmCNV"))
    next
  }

  # Assign CNV status
  cnv_status <- case_when(
    LUAD_armcnv_df[[arm]] > 0.2  ~ "Gain",
    LUAD_armcnv_df[[arm]] < -0.2 ~ "Loss",
    TRUE                         ~ "Neutral"
  )

  # Prepare plot data
  plot_df <- data.frame(
    IS = LUAD_is_df$Ranked_Sum,
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
  file_base <- paste0("LUAD_", arm, "_IS_RankedSum")
  ggsave(paste0(file_base, ".png"), plot = p, width = 4, height = 5, dpi = 300)
  ggsave(paste0(file_base, ".pdf"), plot = p, width = 4, height = 5, device = "pdf")
}

