######################################################## HEATMAP TOP ARM HITS ########################################################
setwd("/gpfs/data/davolilab/projects/PanCan_ImmuneEvasion/") # working dir
source("scripts/get_top_regression_hits.R")
library(ggpubr)
View(arms_reg_out$ACC)
workingtables <- readRDS("data/TCGA_working_tables.rds")
arms_reg_out <- readRDS("IS_CNV_Regression/TCGA/Arm/Regression_Outputs_is0.3_ArmCNV_continuous0.2.rds")
cyto_reg_out <- readRDS("IS_CNV_Regression/TCGA/Cytoband/Regression_Outputs_is0.3_cytobandcnv_continuous0.2.rds")
# remove = c("BRCA","CHOL","COAD", "COAD.MSI", "COAD.MSS", "COADREAD", "DLBC", "ESCA", "HNSC", 
#            "HNSC.TCGAdef.HPVneg", "HNSC.TCGAdef.HPVpos", "HNSC.HPVother", "KIPAN", "READ", "READ.MSI", "READ.MSS", "STES", "UCEC",
#            "ACC","LAML", "MESO", "PCPG", "SARC", "TGCT", "THYM", "UCS")
# cancers = unique(gsub("_.*", "", names(workingtables)))
# cancers = cancers[!cancers %in% remove]


tumor_order = c("UCEC.MSI", "UCEC.MSS", "OV", "CESC", "LUSC", "HNSC.HPVneg", "ESCA.SC","BLCA",
                "BRCA.pos", "BRCA.neg", "PRAD", "LUAD", "LIHC", "COADREAD.MSI", "COADREAD.MSS", "ESCA.AD", "STAD",
                "PAAD", "KICH", "KIRC", "KIRP", "LGG", "GBM", "SKCM", "UVM", "ACC", "MESO", "PCPG", "SARC", "TGCT", "UCS")
tumor_groups = data.frame(row.names=tumor_order, 
                          Group=c("Gyn", "Gyn", "Gyn", "Squamous", "Squamous", "Squamous", "Squamous", "Squamous", 
                                  "Adeno", "Adeno", "Adeno", "Adeno", "Adeno", "GI", "GI", "GI", "GI", "GI",
                                  "Kidney", "Kidney", "Kidney", "NC-Derived", "NC-Derived", "NC-Derived", "NC-Derived", 
                                  "Other", "Other","Other","Other","Other","Other"))
tumor_group_cols <- c("Gyn" = "#6C2DC7", "Squamous" = "#028A0F", "Adeno" = "#6E0B14", 
                      "GI" = "#FFC30B", "Kidney" = "#2B65EC", "NC-Derived" = "#9B111E", "Other" = "#008080")

cancers = unique(gsub("_.*", "", names(workingtables))) # all tumor types
remove = setdiff(cancers, tumor_order) # tumors to remove from analysis


### arms
armAS <- data.frame() #### dataframe with arm level aneuploidy scores (sum of all absolute values of arm level cnv's)
for (tt in tumor_order) {
  armdf <- workingtables[[paste0(tt, "_ArmCNV")]]
  armAS <- rbind(armAS, data.frame(TumorType = rep(tt), Arm_AS = rowSums(abs(armdf))))
}
armAS_mean <- armAS %>% group_by(TumorType) %>% summarize(MeanAS = mean(Arm_AS))

get_top_reg(reg_out = arms_reg_out, AS_df = armAS_mean, remove = remove, cnvlevel = "Arm", continuous = T, 
            tumor_order = tumor_order, tumor_groups = tumor_groups, tumor_group_cols = tumor_group_cols, savetable = T,
            signif_thresh = 0.05, plotvals = "logP", savepath = "IS_CNV_Regression/TCGA/Arm/Top_Regression_Hits/Continuous_Top_Arms_")

armsout <- read.csv("IS_CNV_Regression/TCGA/Arm/Top_Regression_Hits/Continuous_Top_Arms_table.csv")
armscount <- as.data.frame(table(armsout$AllHits))
write.csv(armscount, "IS_CNV_Regression/TCGA/Arm/Top_Regression_Hits/Continuous_Top_Arms_GainsLosses_Counts.csv", row.names = F)

armscount$Event <- ifelse(grepl("Gain",armscount$Var1), "Gain", "Loss")
armscount$Hit <- ifelse(grepl("Cold", armscount$Var1), "Cold", "Hot")
p1 <- ggbarplot(armscount, x = "Hit", y = "Freq", fill = "Event", position = position_dodge(0.8), palette = c("darkred", "darkblue"), title="Arm Hits") + 
  stat_compare_means(aes(group = Event), label =  "p.format")





### cytobands
cytoAS <- data.frame() #### dataframe with cytoband level aneuploidy scores (sum of all absolute values of cytoband level cnv's)
for (tt in cancers) {
  cytodf <- workingtables[[paste0(tt, "_CytoCNV")]]
  cytoAS <- rbind(cytoAS, data.frame(TumorType = rep(tt), Cyto_AS = rowSums(abs(cytodf))))
}
cytoAS_mean <- cytoAS %>% group_by(TumorType) %>% summarize(MeanAS = mean(Cyto_AS))

get_top_reg(reg_out = cyto_reg_out, AS_df = cytoAS_mean, remove = remove, cnvlevel = "Cytoband", continuous = T, 
            tumor_order = tumor_order,tumor_groups = tumor_groups,tumor_group_cols = tumor_group_cols, savetable = T,
            signif_thresh = 0.05, plotvals = "logP", savepath = "IS_CNV_Regression/TCGA/Cytoband/Top_Regression_Hits/Continuous_Top_Cytobands_")

cytosout <- read.csv("IS_CNV_Regression/TCGA/Cytoband/Top_Regression_Hits/Continuous_Top_Cytobands_table.csv")
cytoscount <- as.data.frame(table(cytosout$AllHits))
write.csv(cytoscount, "IS_CNV_Regression/TCGA/Cytoband/Top_Regression_Hits/Continuous_Top_Cytobands_GainsLosses_Counts.csv", row.names = F)

cytoscount$Event <- ifelse(grepl("Gain",cytoscount$Var1), "Gain", "Loss")
cytoscount$Hit <- ifelse(grepl("Cold", cytoscount$Var1), "Cold", "Hot")
p2 <- ggbarplot(cytoscount, x = "Hit", y = "Freq", fill = "Event", position = position_dodge(0.8), palette = c("darkred", "darkblue"), title="Cytoband Hits") + 
  stat_compare_means(aes(group = Event), label =  "p.format")

pdf("IS_CNV_Regression/TCGA/arm_cytobands_gains_losses_hits_barplots.pdf", width = 7, height = 4)
ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = T, legend = "right")
dev.off()
