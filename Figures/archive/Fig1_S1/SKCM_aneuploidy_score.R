library(dplyr)

# Load data
data_table <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")

# Prepare SKCM data
SKCM_is_df <- data_table$SKCM_IS
SKCM_armcnv_df <- data_table$SKCM_ArmCNV

# Keep common patients
common_patients <- intersect(rownames(SKCM_is_df), rownames(SKCM_armcnv_df))
SKCM_is_df <- SKCM_is_df[common_patients, , drop = FALSE]
SKCM_armcnv_df <- SKCM_armcnv_df[common_patients, , drop = FALSE]

# Calculate absolute aneuploidy score
SKCM_is_df$Aneuploidy_Score_Abs <- rowSums(abs(SKCM_armcnv_df), na.rm = TRUE)

# Function to test filtering cutoff
check_pval_by_filter <- function(cutoff_percent) {
  quantile_cutoff <- quantile(SKCM_is_df$Aneuploidy_Score_Abs, probs = 1 - cutoff_percent, na.rm = TRUE)

  # Filter out top X% most aneuploid tumors
  keep <- SKCM_is_df$Aneuploidy_Score_Abs <= quantile_cutoff
  filtered_df <- SKCM_is_df[keep, , drop = FALSE]
  filtered_cnv <- SKCM_armcnv_df[rownames(filtered_df), , drop = FALSE]

  # Define CNV status for Arm1q
  cnv_status <- case_when(
    filtered_cnv$Arm1q > 0.2  ~ "Gain",
    filtered_cnv$Arm1q < -0.2 ~ "Loss",
    TRUE                      ~ "Neutral"
  )

  plot_df <- data.frame(
    IS = filtered_df$Ranked_Sum,
    CNV_Status = factor(cnv_status, levels = c("Loss", "Neutral", "Gain"))
  )

  # Compare Neutral vs Gain
  df_ng <- plot_df %>% filter(CNV_Status %in% c("Neutral", "Gain"))
  p_ng <- wilcox.test(IS ~ CNV_Status, data = df_ng)$p.value

  cat(sprintf("Top %2d%% aneuploidy removed → Wilcoxon p (Neutral vs Gain): %0.5f\n", cutoff_percent * 100, p_ng))
}

# Run for 15%, 20%, 25% removal
check_pval_by_filter(0.15)
check_pval_by_filter(0.20)
check_pval_by_filter(0.25)



#### CHeck 61 & 9p
library(dplyr)

# Step 1: Load SKCM data
SKCM_is_df <- data_table$SKCM_IS
SKCM_armcnv_df <- data_table$SKCM_ArmCNV
common_patients <- intersect(rownames(SKCM_is_df), rownames(SKCM_armcnv_df))
SKCM_is_df <- SKCM_is_df[common_patients, , drop = FALSE]
SKCM_armcnv_df <- SKCM_armcnv_df[common_patients, , drop = FALSE]

# Step 2: Add Aneuploidy Score
SKCM_is_df$Aneuploidy_Score_Abs <- rowSums(abs(SKCM_armcnv_df), na.rm = TRUE)

# Step 3: Loop through both arms
arms_to_check <- c("Arm6q", "Arm9p")

for (arm in arms_to_check) {
  
  # Step 4: Define CNV status for the arm
  cnv_status <- dplyr::case_when(
    SKCM_armcnv_df[[arm]] > 0.2  ~ "Gain",
    SKCM_armcnv_df[[arm]] < -0.2 ~ "Loss",
    TRUE                         ~ "Neutral"
  )

  # Step 5: Construct model data
  plot_df_lm <- data.frame(
    IS = SKCM_is_df$Ranked_Sum,
    Aneuploidy = SKCM_is_df$Aneuploidy_Score_Abs,
    CNV_Status = factor(cnv_status, levels = c("Neutral", "Loss", "Gain"))
  )

  # Step 6: Fit linear model
  model <- lm(IS ~ CNV_Status + Aneuploidy, data = plot_df_lm)

  # Step 7: Print summary
  cat(paste0("\n📊 Linear model for ", arm, ":\n"))
  print(summary(model))
}