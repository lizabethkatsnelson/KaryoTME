library(tidyverse)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)

# Load input files
clin_pt  <- read.delim("data_clinical_patient.txt", skip = 4, check.names = FALSE)
clin_smp <- read.delim("data_clinical_sample.txt", skip = 4, check.names = FALSE)

# Load CNV table and clean SAMPLE_ID
cnv_raw <- read.delim("all_cnv_table.txt", check.names = FALSE, quote = "\"") %>%
  tibble::rownames_to_column("SAMPLE_ID")

# Convert SAMPLE_ID columns to character
clin_smp$SAMPLE_ID <- as.character(clin_smp$SAMPLE_ID)
cnv_raw$SAMPLE_ID  <- as.character(cnv_raw$SAMPLE_ID)

# Identify cytoband and arm-level CNV columns
cytoband_cols <- grep("^[0-9XY]+[pq]", names(cnv_raw), value = TRUE)
arm_cols <- grep("^Arm", names(cnv_raw), value = TRUE)

# Categorize CNVs into Gain / Loss / None
cnv_cat <- cnv_raw %>%
  mutate(across(all_of(c(cytoband_cols, arm_cols)), ~ case_when(
    .x >  0.2 ~ "Gain",
    .x < -0.2 ~ "Loss",
    TRUE      ~ "None"
  )))

### ADJUSTED FOR ANEUPLOIDY SCORE + FILTER TO ONLY PD1 PATIENTS ### 
# STEP 1: Calculate aneuploidy score
aneuploidy_cytoband <- cnv_raw %>%
  select(SAMPLE_ID, all_of(cytoband_cols)) %>%
  mutate(AneuploidyScore = rowSums(abs(across(all_of(cytoband_cols))), na.rm = TRUE)) %>%
  select(SAMPLE_ID, AneuploidyScore)

# STEP 2: Merge everything
sample_info_extended <- clin_smp %>%
  select(SAMPLE_ID, PATIENT_ID, CANCER_TYPE, ONCOTREE_CODE, TUMOR_PURITY) %>%
  left_join(clin_pt[, c("PATIENT_ID", "DRUG_TYPE")], by = "PATIENT_ID") %>%
  left_join(aneuploidy_cytoband, by = "SAMPLE_ID")

# Merge CNV calls
merged_sample_cytoband <- merge(sample_info_extended, cnv_cat, by = "SAMPLE_ID")

# Final data frame with survival
df_cyto <- merge(
  clin_pt[, c("PATIENT_ID", "OS_MONTHS", "OS_STATUS")],
  merged_sample_cytoband,
  by = "PATIENT_ID"
) %>%
  mutate(OS = if_else(str_detect(OS_STATUS, "LIVING"), 0, 1))

# STEP 3: Cox models (adjust for aneuploidy), only PD-1 patients
# HR = hazard of gain or loss/hazard of none 
# both a Gain and a Loss result for the same cytoband and same tumor
results_cytoband <- data.frame()
for (tumor in unique(df_cyto$CANCER_TYPE)) {
  tumor_data <- df_cyto %>%
    filter(CANCER_TYPE == tumor, DRUG_TYPE == "PD-1/PDL-1")  
  
  for (band in cytoband_cols) {
    for (comparison in c("Gain", "Loss")) {
      subset <- tumor_data %>%
        filter(!!sym(band) %in% c("None", comparison)) %>%
        mutate(CNV = factor(!!sym(band), levels = c("None", comparison))) %>%
        filter(!is.na(OS), !is.na(OS_MONTHS), !is.na(AneuploidyScore))
      
      if (nrow(subset) >= 10 && length(unique(subset$CNV)) == 2) {
        tryCatch({
          fit <- coxph(Surv(OS_MONTHS, OS) ~ CNV + AneuploidyScore, data = subset)
          cnv_term <- grep("^CNV", names(coef(fit)), value = TRUE)
          if (length(cnv_term) != 1) next
          
          hr <- exp(coef(fit)[cnv_term])
          ci <- round(exp(confint(fit)[cnv_term, ]), 3)
          p <- summary(fit)$coefficients[cnv_term, 5]
          
          results_cytoband <- rbind(results_cytoband, data.frame(
            CANCER_TYPE = tumor,
            Cytoband = band,
            Comparison = comparison,
            HR = round(hr, 3),
            log2HR = log2(hr),
            CI = paste(ci[1], ci[2], sep = " - "),
            P = signif(p, 3),
            N = nrow(subset)
          ))
        }, error = function(e) {
          message(paste("Skipped:", tumor, band, comparison, "-", e$message))
        })
      }
    }
  }
}

# STEP 4: Format & Plot — only significant, label arm at first cytoband
results_cytoband <- results_cytoband %>%
  filter(P < 0.05) %>%  # only significant
  mutate(
    Event = ifelse(Comparison == "Gain", "Gain", "Loss"),
    Chrom = str_extract(Cytoband, "^[0-9XY]+"),
    ArmLabel = str_extract(Cytoband, "^[0-9XY]+[pq]")
  )

cyto_levels <- results_cytoband %>%
  mutate(
    chr_num = as.numeric(str_extract(Cytoband, "^[0-9]+")),
    chr_num = ifelse(is.na(chr_num), ifelse(str_detect(Cytoband, "^X"), 23, 24), chr_num),
    arm_suffix = str_extract(Cytoband, "[pq]$")
  ) %>%
  arrange(chr_num, arm_suffix) %>%
  pull(Cytoband) %>%
  unique()

results_cytoband$Cytoband <- factor(results_cytoband$Cytoband, levels = cyto_levels)

# Label arm at FIRST cytoband (NOT middle)
axis_labels <- results_cytoband %>%
  group_by(ArmLabel) %>%
  summarise(Cytoband = Cytoband[1], .groups = "drop")

# Plotting (labeled by arm)
dir.create("logHR_plots_filtered", showWarnings = FALSE)

for (tumor in unique(results_cytoband$CANCER_TYPE)) {
  plot_data <- results_cytoband %>%
    filter(CANCER_TYPE == tumor)
  
  if (nrow(plot_data) == 0) next
  
  p <- ggplot(plot_data, aes(x = Cytoband, y = log2HR, fill = Event, group = Event)) +
    geom_bar(stat = "identity", width = 0.9, position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("Gain" = "#EA6A6A", "Loss" = "#4F6DB8")) +  # pink and blue
    scale_x_discrete(
      breaks = axis_labels$Cytoband,
      labels = axis_labels$ArmLabel
    ) +
    labs(
      title = paste0("Significant CNV Survival Effects in ", tumor),
      x = "Chromosome Arm",
      y = "log2(Hazard Ratio)",
      fill = "CNV Event"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, size = 6),
      plot.title = element_text(face = "bold", size = 12)
    )
  
  ggsave(
    filename = paste0("logHR_plots_filtered/", gsub(" ", "_", tumor), ".pdf"),
    plot = p, width = 10, height = 4
  )
}


## Plot by cytobands
# Optional: remove or keep P filter as you want
results_cytoband <- results_cytoband %>%
  filter(P < 0.05) %>%  # keep only significant if desired
  mutate(
    Event = ifelse(Comparison == "Gain", "Gain", "Loss")
  )

# Order cytobands along the x-axis
cyto_levels <- results_cytoband %>%
  mutate(
    chr_num = as.numeric(str_extract(Cytoband, "^[0-9]+")),
    chr_num = ifelse(is.na(chr_num), ifelse(str_detect(Cytoband, "^X"), 23, 24), chr_num),
    arm_suffix = str_extract(Cytoband, "[pq]$")
  ) %>%
  arrange(chr_num, arm_suffix) %>%
  pull(Cytoband) %>%
  unique()

results_cytoband$Cytoband <- factor(results_cytoband$Cytoband, levels = cyto_levels)

# Plot directory
dir.create("logHR_filtered_cytoband", showWarnings = FALSE)

# Loop by tumor
for (tumor in unique(results_cytoband$CANCER_TYPE)) {
  plot_data <- results_cytoband %>%
    filter(CANCER_TYPE == tumor)
  
  if (nrow(plot_data) == 0) next
  
  p <- ggplot(plot_data, aes(x = Cytoband, y = log2HR, fill = Event, group = Event)) +
    geom_bar(stat = "identity", width = 0.9, position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("Gain" = "#EA6A6A", "Loss" = "#4F6DB8")) +
    labs(
      title = paste0("Significant CNV Survival Effects in ", tumor),
      x = "Cytoband",
      y = "log2(Hazard Ratio)",
      fill = "CNV Event"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5, hjust = 1),
      plot.title = element_text(face = "bold", size = 12)
    )
  
  ggsave(
    filename = paste0("logHR_filtered_cytoband/", gsub(" ", "_", tumor), "_cytoband.pdf"),
    plot = p, width = 12, height = 4
  )
}