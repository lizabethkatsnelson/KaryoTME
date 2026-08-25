# Setup
setwd("/gpfs/data/proteomics/projects/")
library(tidyverse)

tilepath <- "/gpfs/data/proteomics/projects/Lisa/Image_Models/tiles/"

# =========================
# Load & process CPTAC
# =========================
cptac <- read.csv("Wenke/theme5_img/idx_files/pancan_imaging_all_tiles_index_TN.csv", check.names = FALSE)
colnames(cptac) # show column names 

# Subset tumor types
cptac_list <- list(
  UCEC = cptac[cptac$Tumor == "UCEC" & cptac$Tumor_Normal == "tumor", ],
  LUAD = cptac[cptac$Tumor == "LUAD" & cptac$Tumor_Normal == "tumor", ],
  HNSC = cptac[cptac$Tumor == "HNSCC" & cptac$Tumor_Normal == "tumor", ]
)

# Fix tile paths
cptac_list <- lapply(cptac_list, function(df) {
  df[, c("L1path", "L2path", "L3path")] <- apply(df[, c("L1path", "L2path", "L3path")], 2,
    function(x) gsub("../tiles/.+?/", tilepath, x))
  df[, !colnames(df) %in% c("Tumor", "Tumor_Normal", "tumor_code")]
})

# =========================
# Load & process TCGA
# =========================
tcga <- read.csv("/gpfs/data/proteomics/projects/Josh_imaging_TCGA/all_tiles.csv", check.names = FALSE)

tcga_ucec <- tcga[tcga$cancerType == "UCEC" & tcga$NT == "tumor", ]
tcga_luad <- tcga[tcga$cancerType == "LUAD" & tcga$NT == "tumor", ]
tcga_hnsc <- tcga[tcga$cancerType == "HNSC" & tcga$NT == "tumor", ]

# Unify formats
for (tcga_df in list(tcga_ucec, tcga_luad, tcga_hnsc)) {
  tcga_df <- tcga_df[, !colnames(tcga_df) %in% c("cancerType", "NT", "max")]
  colnames(tcga_df) <- colnames(cptac_ucec)
  tcga_df[, c("L1path", "L2path", "L3path")] <- apply(tcga_df[, c("L1path", "L2path", "L3path")], 2, function(x) gsub("/gpfs/data/proteomics/projects/Josh_imaging_TCGA/tiles/.+?/", tilepath, x))
  tcga_df$Slide_ID <- paste0(tcga_df$Patient_ID, "-", tcga_df$Slide_ID)
}

# =========================
# Merge and write tiles
# =========================
save_tiles <- function(cptac_df, tcga_df, tumor_name, out_file) {
  all_df <- rbind(cptac_df, tcga_df)
  all_df <- cbind(all_df[,1:2], Tumor = tumor_name, all_df[,3:5])
  write.csv(all_df, out_file, row.names = FALSE)
}

save_tiles(cptac_ucec, tcga_ucec, "UCEC", "/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/CPTAC/UCEC_all_tiles.csv")
save_tiles(cptac_luad, tcga_luad, "LUAD", "/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/CPTAC/LUAD_all_tiles.csv")
save_tiles(cptac_hnsc, tcga_hnsc, "HNSC", "/gpfs/data/davolilab/projects/melanoma_scrna/bulk_analysis/CPTAC/HNSC_all_tiles.csv")

# =========================
# CNV Labels
# =========================
# CPTAC labels
cptac_labels <- read.csv("/gpfs/scratch/lk2513/cptac_pancan/CPTAC_PanCan_Arm_Level_CopyNumber_NORMALIZED.csv", check.names = FALSE, row.names = 1)

ucec_lbl <- cptac_labels[cptac_labels$TumorType == "UCEC", c("MSS_status", "Arm_1q")] %>%
  filter(MSS_status == "MSS_High") %>%
  mutate(Arm_1q_Gain = ifelse(Arm_1q > 0.2, 1, 0)) %>%
  rownames_to_column("Patient_ID") %>%
  select(Patient_ID, Tumor = "UCEC", Arm_1q_Gain)

luad_lbl <- cptac_labels[cptac_labels$TumorType == "LUAD", "Arm_9p", drop=FALSE] %>%
  mutate(Arm_9p_Loss = ifelse(Arm_9p < -0.2, 1, 0)) %>%
  rownames_to_column("Patient_ID") %>%
  mutate(Tumor = "LUAD") %>%
  select(Patient_ID, Tumor, Arm_9p_Loss)

hnsc_lbl <- cptac_labels[cptac_labels$TumorType == "HNSCC", "Arm_9p", drop=FALSE] %>%
  mutate(Arm_9p_Loss = ifelse(Arm_9p < -0.2, 1, 0)) %>%
  rownames_to_column("Patient_ID") %>%
  mutate(Tumor = "HNSC") %>%
  select(Patient_ID, Tumor, Arm_9p_Loss)

# TCGA labels
tcga_rds <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")

tcga_ucec <- data.frame(Patient_ID = rownames(tcga_rds$UCEC.MSS_High_ArmCNV),
                        Arm_1q = tcga_rds$UCEC.MSS_High_ArmCNV$Arm1q) %>%
  mutate(Arm_1q_Gain = ifelse(Arm_1q > 0.2, 1, 0), Tumor = "UCEC") %>%
  select(Patient_ID, Tumor, Arm_1q_Gain)

tcga_luad <- data.frame(Patient_ID = rownames(tcga_rds$LUAD.ArmCNV),
                        Arm_9p = tcga_rds$LUAD.ArmCNV$Arm9p) %>%
  mutate(Arm_9p_Loss = ifelse(Arm_9p < -0.2, 1, 0), Tumor = "LUAD") %>%
  select(Patient_ID, Tumor, Arm_9p_Loss)

tcga_hnsc <- data.frame(Patient_ID = rownames(tcga_rds$HNSC.HPVneg_ArmCNV),
                        Arm_9p = tcga_rds$HNSC.HPVneg_ArmCNV$Arm9p) %>%
  mutate(Arm_9p_Loss = ifelse(Arm_9p < -0.2, 1, 0), Tumor = "HNSC") %>%
  mutate(Patient_ID = paste0(Patient_ID, "-tumor")) %>%
  select(Patient_ID, Tumor, Arm_9p_Loss)

# Save combined labels
write.csv(rbind(ucec_lbl, tcga_ucec), "Lisa/Image_Models/UCEC/arm1q_gain_labels.csv", row.names = FALSE)
write.csv(rbind(luad_lbl, tcga_luad), "Lisa/Image_Models/LUAD/arm9p_loss_labels.csv", row.names = FALSE)
write.csv(rbind(hnsc_lbl, tcga_hnsc), "Lisa/Image_Models/HNSC/arm9p_loss_labels.csv", row.names = FALSE)