setwd("/gpfs/data/proteomics/projects/")
library(tidyverse)

tilepath <- "/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/tiles/"


cptac <- read.csv("Wenke/theme5_img/idx_files/pancan_imaging_all_tiles_index_TN.csv", check.names = F)
cptac_hnsc <- cptac[cptac$Tumor %in% "HNSCC" & cptac$Tumor_Normal %in% "tumor", ]
cptac_hnsc <- cptac_hnsc[, !colnames(cptac_hnsc) %in% c("Tumor", "tumor_code", "Tumor_Normal")]
cptac_hnsc[,c("L1path", "L2path", "L3path")] <- apply(cptac_hnsc[,c("L1path", "L2path", "L3path")], 2, function(x) gsub("../tiles/HNSCC/", tilepath, x))


tcga <- read.csv("/gpfs/data/proteomics/projects/Josh_imaging_TCGA/all_tiles.csv", check.names = F)
tcga_hnsc <- tcga[tcga$cancerType %in% "HNSC" & tcga$NT %in% "tumor", ]
tcga_hnsc <- tcga_hnsc[, !colnames(tcga_hnsc) %in% c("cancerType", "NT", "max")]
colnames(tcga_hnsc) <- colnames(cptac_hnsc)
tcga_hnsc[,c("L1path", "L2path", "L3path")] <- apply(tcga_hnsc[,c("L1path", "L2path", "L3path")], 2, function(x) gsub("/gpfs/data/proteomics/projects/Josh_imaging_TCGA/tiles/HNSC/", tilepath, x))
tcga_hnsc$Slide_ID <- paste0(tcga_hnsc$Patient_ID, "-", tcga_hnsc$Slide_ID)


allhnsc <- rbind(cptac_hnsc, tcga_hnsc)
allhnsc <- cbind(allhnsc[,1:2], Tumor = rep("HNSC"), allhnsc[,3:5])
write.csv(allhnsc, "Lisa/Image_Models/HNSC/all_tiles.csv", row.names = F)
# unique(allhnsc$Patient_ID)


cptac_labels <- read.csv("/gpfs/scratch/lk2513/cptac_pancan/CPTAC_PanCan_Arm_Level_CopyNumber_NORMALIZED.csv", check.names = F, row.names = 1)
cptac_labels <- cptac_labels[cptac_labels$TumorType %in% "HNSCC", "Arm_9p", drop=F] %>% rownames_to_column("Patient_ID")
cptac_labels$Arm_9p_Loss <- ifelse(cptac_labels$Arm_9p < -0.2, 1, 0)
# table(cptac_labels$Arm_9p_Loss)/nrow(cptac_labels)

tcgaworkingtables <- readRDS("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/data/TCGA_working_tables.rds")
tcga_labels <- data.frame(Patient_ID = rownames(tcgaworkingtables$HNSC.HPVneg_ArmCNV), Arm_9p = tcgaworkingtables$HNSC.HPVneg_ArmCNV$Arm9p)
tcga_labels$Arm_9p_Loss <- ifelse(tcga_labels$Arm_9p < -0.2, 1, 0)
tcga_labels$Patient_ID <- paste0(tcga_labels$Patient_ID, "-tumor")
# table(tcga_labels$Arm_9p_Loss)/nrow(tcga_labels)


alllabels <- rbind(cptac_labels, tcga_labels)
# table(alllabels$Arm_9p_Loss)/nrow(alllabels)
alllabels <- cbind(alllabels[,1,drop=F], Tumor = rep("HNSC"), alllabels[,3,drop=F])
write.csv(alllabels, "Lisa/Image_Models/HNSC/arm9p_loss_labels.csv", row.names = F)



##### 9p loss model eval (single res 10x) - accidentally only ran cptac samples in this run
setwd("/gpfs/data/proteomics/projects/Lisa/Image_Models/HNSC/results/I6/")
library(tidyverse)
library(ggpubr)

traindf <- read.csv("data/trn_slide_idx.csv")
nrow(traindf)
valdf <- read.csv("data/val_slide_idx.csv")
nrow(valdf)
testdf <- read.csv("data/tst_slide_idx.csv")
nrow(testdf)


tsnedf <- read.csv("pred/tSNE_P_N.csv")
tsnedf$label <- ifelse(tsnedf$label == 1, "Arm9p_Loss", "Arm9p_Neutral_or_Gain")
ggscatter(tsnedf, x = "tsne_0", y = "tsne_1", color = "label")

# tsnedf
test_tile_pred <- read.csv("pred/tst_tile_pred.csv")
test_tile_pred_filt <- test_tile_pred[, !colnames(test_tile_pred) %in% c("L2path", "L3path")]



