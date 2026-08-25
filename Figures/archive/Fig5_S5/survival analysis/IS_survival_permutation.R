library(tidyverse)

# Get tumor types from data_table names
tumor_list <- unique(gsub("_(CytoCNV|IS).*", "", names(data_table)))

all_results <- list()

for (ct in tumor_list) {
  # Get CNV and IS data if they exist
  cnv_key <- paste0(ct, "_CytoCNV")
  is_key  <- paste0(ct, "_IS")

  if (!(cnv_key %in% names(data_table)) || !(is_key %in% names(data_table))) next

  cnv_df <- data_table[[cnv_key]]
  is_df  <- data_table[[is_key]]

  # Match sample IDs
  common <- intersect(rownames(cnv_df), rownames(is_df))
  if (length(common) < 10) next  # skip if not enough samples

  cnv_df <- cnv_df[common, , drop = FALSE]
  is_vec <- is_df[common, 1]  # assuming 1 column = TotalIS

  results_ct <- data.frame()

  for (cyto in grep("^[0-9XY]+[pq]", colnames(cnv_df), value = TRUE)) {
    cnv_status <- cnv_df[[cyto]]
    if (length(unique(cnv_status)) < 2) next

    fit <- lm(is_vec ~ cnv_status)
    est <- coef(summary(fit))["cnv_status", "Estimate"]
    p   <- coef(summary(fit))["cnv_status", "Pr(>|t|)"]
    direction <- ifelse(est > 0, "Gain", "Loss")

    results_ct <- rbind(results_ct, data.frame(
      OncoTree = ct,
      Cytoband = cyto,
      Comparison = direction,
      ImmuneEffect = est,
      P_immune = p
    ))
  }

  all_results[[ct]] <- results_ct
}

# Combine all results
immune_results_all <- bind_rows(all_results)

# Save to CSV
write.csv(immune_results_all, "results_immune_allTumors.csv", row.names = FALSE)















# CNV Immune Score vs Survival Association Comparison
# Comparing Figure 2 (immune associations) with survival data

library(tidyverse)
library(survival)
library(gtools)

# ===== 1. SURVIVAL DATA PROCESSING =====
# Load and process your clinical survival data

# Load input files
clin_pt  <- read.delim("/gpfs/data/davolilab/projects/MSK_IMPACT/data_clinical_patient.txt", skip = 4, check.names = FALSE)
clin_smp <- read.delim("/gpfs/data/davolilab/projects/MSK_IMPACT/data_clinical_sample.txt", skip = 4, check.names = FALSE)

# Load CNV table and clean SAMPLE_ID
cnv_raw <- read.delim("all_cnv_table.txt", check.names = FALSE, quote = "\"") %>%
  tibble::rownames_to_column("SAMPLE_ID")

clin_smp$SAMPLE_ID <- as.character(clin_smp$SAMPLE_ID)
cnv_raw$SAMPLE_ID  <- as.character(cnv_raw$SAMPLE_ID)

cytoband_cols <- grep("^[0-9XY]+[pq]", names(cnv_raw), value = TRUE)
arm_cols <- grep("^Arm", names(cnv_raw), value = TRUE)

cnv_cat <- cnv_raw %>%
  mutate(across(all_of(c(cytoband_cols, arm_cols)), ~ case_when(
    .x >  0.3 ~ "Gain",
    .x < -0.3 ~ "Loss",
    TRUE      ~ "None"
  )))

## MERGE
# Step 1: Prepare clean sample info
sample_info <- clin_smp %>%
  select(SAMPLE_ID, PATIENT_ID, CANCER_TYPE, ONCOTREE_CODE, SAMPLE_TYPE)

# Step 2: Join sample_info with CNV matrix
merged_sample <- sample_info %>%
  left_join(cnv_cat, by = "SAMPLE_ID")

# Step 3: Join patient-level clinical data
df <- clin_pt %>%
  select(PATIENT_ID, OS_MONTHS, OS_STATUS, SEX, DRUG_TYPE) %>%
  left_join(merged_sample, by = "PATIENT_ID") %>%
  mutate(OS = if_else(str_detect(OS_STATUS, "LIVING"), 0, 1))

# Generate survival results using ONCOTREE_CODE to match immune data
cytoband_cols_df <- grep("^[0-9XY]+[pq]", names(df), value = TRUE)
results_survival <- data.frame()

for (tumor in unique(df$ONCOTREE_CODE)) {  # Use ONCOTREE_CODE not CANCER_TYPE
  tumor_data <- df %>% filter(ONCOTREE_CODE == tumor, DRUG_TYPE == "PD-1/PDL-1")
  
  for (band in cytoband_cols_df) {
    for (comparison in c("Gain", "Loss")) {
      subset <- tumor_data %>%
        filter(!!sym(band) %in% c("None", comparison)) %>%
        mutate(CNV = factor(!!sym(band), levels = c("None", comparison)))
      
      if (nrow(subset) >= 10 && length(unique(subset$CNV)) == 2) {
        tryCatch({
          fit <- coxph(Surv(OS_MONTHS, OS) ~ CNV, data = subset)
          hr <- exp(coef(fit))
          p <- summary(fit)$coefficients[5]
          
          results_survival <- rbind(results_survival, data.frame(
            OncoTree = tumor,  # Match naming with immune data
            Cytoband = band,
            Comparison = comparison,
            HR = round(hr, 3),
            P_survival = p,
            log2HR = log2(hr),
            sign_logP_survival = -log10(p) * ifelse(log2(hr) > 0, 1, -1),
            N = nrow(subset)
          ))
        }, error = function(e){})
      }
    }
  }
}

# ===== 2. IMMUNE SCORE DATA =====
# (Your existing immune score code)
data_table <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")

# Get the cancer types that exist in your survival data
cancers <- unique(df$ONCOTREE_CODE)

immune_all <- data.frame()
for (ct in cancers) {
  cat("Processing", ct, "\n")
  is_df <- data_table[[paste0(ct, "_IS")]]
  cnv_df <- data_table[[paste0(ct, "_CytoCNV")]]
  
  if (is.null(is_df) || is.null(cnv_df)) next
  
  # Match samples
  common <- intersect(rownames(is_df), rownames(cnv_df))
  is_df <- is_df[common, , drop = FALSE]
  cnv_df <- cnv_df[common, , drop = FALSE]
  as_score <- rowSums(abs(cnv_df), na.rm = TRUE)
  
  for (cyto in grep("^[0-9XY]+[pq]", names(cnv_df), value = TRUE)) {
    cnv_status <- case_when(
      cnv_df[[cyto]] > 0.3 ~ "Gain",
      cnv_df[[cyto]] < -0.3 ~ "Loss",
      TRUE ~ "Neutral"
    )
    
    df_temp <- data.frame(
      IS = is_df$Ranked_Sum,
      CNV_Status = factor(cnv_status, levels = c("Neutral", "Loss", "Gain")),
      Aneuploidy = as_score
    )
    
    model <- lm(IS ~ CNV_Status + Aneuploidy, data = df_temp)
    coefs <- summary(model)$coefficients
    
    if ("CNV_StatusGain" %in% rownames(coefs)) {
      p <- coefs["CNV_StatusGain", "Pr(>|t|)"]
      eff <- coefs["CNV_StatusGain", "Estimate"]
      immune_all <- rbind(immune_all, data.frame(
        OncoTree = ct, 
        Cytoband = cyto, 
        Comparison = "Gain",
        P_immune = p,
        Effect_immune = eff,
        sign_logP_immune = -log10(p) * sign(eff)
      ))
    }
    
    if ("CNV_StatusLoss" %in% rownames(coefs)) {
      p <- coefs["CNV_StatusLoss", "Pr(>|t|)"]
      eff <- coefs["CNV_StatusLoss", "Estimate"]
      immune_all <- rbind(immune_all, data.frame(
        OncoTree = ct, 
        Cytoband = cyto, 
        Comparison = "Loss",
        P_immune = p,
        Effect_immune = eff,
        sign_logP_immune = -log10(p) * sign(eff)
      ))
    }
  }
}

# ===== 3. MERGE AND CORRELATE =====
# Merge immune and survival data
merged_data <- immune_all %>%
  inner_join(results_survival, 
             by = c("OncoTree", "Cytoband", "Comparison")) %>%
  filter(!is.na(sign_logP_immune) & !is.na(sign_logP_survival))

# Calculate correlation
correlation_result <- cor.test(merged_data$sign_logP_immune, 
                              merged_data$sign_logP_survival, 
                              method = "spearman")

cat("Spearman correlation between immune and survival associations:\n")
cat("Correlation:", round(correlation_result$estimate, 3), "\n")
cat("P-value:", signif(correlation_result$p.value, 3), "\n")

# ===== 4. PERMUTATION TEST =====
# Function to perform permutation test
perform_permutation_test <- function(merged_data, n_permutations = 1000) {
  # Observed correlation
  observed_cor <- cor(merged_data$sign_logP_immune, 
                     merged_data$sign_logP_survival, 
                     method = "spearman", use = "complete.obs")
  
  # Permutation test
  permuted_cors <- replicate(n_permutations, {
    # Scramble the survival associations while keeping immune fixed
    shuffled_survival <- sample(merged_data$sign_logP_survival)
    cor(merged_data$sign_logP_immune, shuffled_survival, 
        method = "spearman", use = "complete.obs")
  })
  
  # Calculate p-value
  p_value <- mean(abs(permuted_cors) >= abs(observed_cor))
  
  return(list(
    observed_correlation = observed_cor,
    permuted_correlations = permuted_cors,
    p_value = p_value
  ))
}

# Run permutation test
cat("Running permutation test...\n")
perm_result <- perform_permutation_test(merged_data, n_permutations = 1000)

cat("Permutation test results:\n")
cat("Observed correlation:", round(perm_result$observed_correlation, 3), "\n")
cat("Permutation p-value:", round(perm_result$p_value, 3), "\n")

# ===== 5. VISUALIZE RESULTS =====
# Scatter plot of immune vs survival associations
p1 <- ggplot(merged_data, aes(x = sign_logP_immune, y = sign_logP_survival)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    title = "CNV Associations: Immune Score vs Survival",
    subtitle = paste("Spearman r =", round(correlation_result$estimate, 3), 
                     ", p =", signif(correlation_result$p.value, 3)),
    x = "Signed -log10(p) - Immune Association",
    y = "Signed -log10(p) - Survival Association"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

# Histogram of permuted correlations
p2 <- ggplot(data.frame(cors = perm_result$permuted_correlations), aes(x = cors)) +
  geom_histogram(bins = 50, alpha = 0.7, fill = "lightblue", color = "black") +
  geom_vline(xintercept = perm_result$observed_correlation, 
             color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Permutation Test Results",
    subtitle = paste("Observed correlation:", round(perm_result$observed_correlation, 3),
                     "| Permutation p =", round(perm_result$p_value, 3)),
    x = "Permuted Correlation Coefficients",
    y = "Frequency"
  ) +
  theme_bw()

# Display plots
print(p1)
print(p2)

# ===== 6. FOCUS ON SPECIFIC REGIONS =====
# Look specifically at 9p and 1q regions as mentioned
focus_regions <- merged_data %>%
  filter(grepl("^9p|^1q", Cytoband)) %>%
  arrange(desc(abs(sign_logP_survival)))

cat("\nTop associations in 9p and 1q regions:\n")
print(focus_regions %>% 
      select(OncoTree, Cytoband, Comparison, sign_logP_immune, sign_logP_survival) %>%
      head(10))

# Save results
write.csv(merged_data, "immune_survival_comparison.csv", row.names = FALSE)
write.csv(focus_regions, "9p_1q_immune_survival_associations.csv", row.names = FALSE)

# Summary statistics
cat("\nSummary:\n")
cat("Total CNV-cancer type combinations analyzed:", nrow(merged_data), "\n")
cat("Cancer types included:", length(unique(merged_data$OncoTree)), "\n")
cat("Unique cytobands analyzed:", length(unique(merged_data$Cytoband)), "\n")