library(readxl)
# library(openxlsx)
library(dplyr)
library(tidyverse)
library(ggpubr)

### chr arm order
chrarm_order <- expand.grid(c(paste0("Chr", 1:22), "ChrX"), c("p", "q"))
chrarm_order <- chrarm_order[order(chrarm_order$Var1),]
chrarm_order <- subset(chrarm_order, !(Var1 %in% c("Chr13", "Chr14", "Chr15", "Chr21", "Chr22") & Var2 == "p"))
chrarm_order <- gsub("Chr", "", paste(chrarm_order$Var1, chrarm_order$Var2, sep = ""))


#### read in genes from outputs
control_genes <- read.csv("ML_Models/control_genes.csv")
control_genesdf <- data.frame(Gene = control_genes$Gene, Control= control_genes$Control, Prediction = rep(NA), Random_Forest=rep(NA),
                              Chromosome = control_genes$Chromosome, Arm=control_genes$Arm, Cytoband = control_genes$Cytoband, 
                              Immune_Genes=rep(NA), Amplification_Frequency=rep(NA), Deletion_Frequency=rep(NA))
control_genesdf$Prediction <- ifelse(control_genesdf$Control == "P", "Positive", "Negative")


tumors <- excel_sheets("ML_Models/predicted_genes_RF_per_tumortype_mean_threshold.xlsx")

ml_outputs <- list()
for (tt in tumors) {
  ml_outputs[[tt]] <- read_xlsx("ML_Models/predicted_genes_RF_per_tumortype_mean_threshold.xlsx", sheet = tt)
}
names(ml_outputs)[3] <- "HNSC.HPVneg"
tumors=names(ml_outputs)
View(ml_outputs$BRCA)

### change threshold to top 75% and bottom 35%
# for (tt in tumors) {
#   df <- ml_outputs[[tt]]
  # posthresh <- unname(unlist(df %>% filter(Control %in% "P") %>% summarize(quantile(Random_Forest, probs = .7, rm.na=T)))) 
  # negthresh <- unname(unlist(df %>% filter(Control %in% "N") %>% summarize(quantile(Random_Forest, probs = .3, rm.na=T)))) 
  # df$Prediction_70_30 <- ifelse(df$Random_Forest > posthresh, "Positive", 
  #                               ifelse(df$Random_Forest < negthresh, "Loss", NA))
#   ml_outputs[[tt]] <- df
# }

### genes that are lost and predicted to be immune cold hits from ML models
predicted_genes <- list()  
for (tt in tumors) {
  df <- ml_outputs[[tt]]
  predicted_genes[[tt]] <- df[df$Prediction %in% "Positive" & df$Deletion_Frequency > 0.2, ] # pos predictions and del freq > 20%
}
View(predicted_genes$BRCA)


##### gene hits from regression models
reg_genes <- read.csv("IS_CNV_Regression/TCGA/Gene/Top_Gene_ImmuneCold_Hits.csv")
reg_genes$Chromosome <- factor(reg_genes$Chromosome, levels = c(paste0("Chr", 1:22), "ChrX", "ChrY"))
reg_genes$Arm <- factor(reg_genes$Arm, levels = c("p", "q"))
reg_genes <- reg_genes[reg_genes$TumorType %in% tumors, -c(1:2)] # filter for tumor types in ML models
reg_genes$TumorType <- factor(reg_genes$TumorType, levels = sort(unique(reg_genes$TumorType)))
reg_genes$ImmuneHot <- ifelse(reg_genes$Z_val > 0 & reg_genes$P_val < 0.05 & reg_genes$has_diff_signs == T, "Yes", "No") #### immune hot hits

### losses, immune cold
reg_genes_filt <- reg_genes[reg_genes$Event == "Loss" & reg_genes$Hits == "Yes",] # only losses and immune hits
reg_genes_filt <- reg_genes_filt[order(reg_genes_filt$TumorType, reg_genes_filt$Chromosome, reg_genes_filt$Arm),]

### losses, immune hot
reg_genes_hot <- reg_genes[reg_genes$Event == "Loss" & reg_genes$ImmuneHot == "Yes",] # only losses and immune hot genes


### intersect immune COLD hits from reg and ML models
all_intersect_genes <- list()
for (tt in tumors) {
  mlgenes <- rbind(control_genesdf, predicted_genes[[tt]])
  mlgenes$Cytoband <- paste0(gsub("Chr", "", mlgenes$Chromosome), mlgenes$Cytoband )
  
  reggenes <- reg_genes_filt[reg_genes_filt$TumorType %in% tt,]
  
  in_both <- intersect(mlgenes$Gene, reggenes$Gene) # intersection
  only_in_ML <- setdiff(mlgenes$Gene, reggenes$Gene) # only in ML output
  only_in_reg <- setdiff(reggenes$Gene, mlgenes$Gene) # only in reg output
  
  savedf <- data.frame(TumorType = rep(tt), Gene = c(in_both, only_in_ML, only_in_reg), 
                       Intersection = c(rep("Both", length(in_both)), rep("ML_Only", length(only_in_ML)), rep("Reg_Only", length(only_in_reg))))
  
  add_loc <- unique(rbind(data.frame(Gene = mlgenes$Gene, Chromosome = mlgenes$Chromosome, Arm = mlgenes$Arm, Cytoband = mlgenes$Cytoband),
                          data.frame(Gene = reggenes$Gene, Chromosome = reggenes$Chromosome, Arm = reggenes$Arm, Cytoband = reggenes$Cytoband)))
  
  savedf <- merge(savedf, add_loc, by="Gene")
  all_intersect_genes[[tt]] <- savedf
}
View(all_intersect_genes$BRCA)

### function to check for dup genes (multiple regions accidentally added)
check_dup_genes <- function(df) {
  checkgenes <- df[duplicated(df$Gene),"Gene"]
  test <- df[df$Gene %in% checkgenes,]
  return(test)
}
dupgenesdf <- data.frame()
for (tt in tumors) {
  dupgenesdf <- rbind(dupgenesdf, check_dup_genes(all_intersect_genes[[tt]]))
}
# dupgenesdf <- unique(dupgenesdf[,c("Gene", "Cytoband")])
fixgenes <- read.csv("testgenes.csv") ##### looked up in original tcga tables which gene locations are correct

for (tt in unique(dupgenesdf$TumorType)) { ### fix genes in the tables
  dupgenes <- dupgenesdf[dupgenesdf$TumorType == tt,] # which genes are dup per tt
  fixgenes_filt <- fixgenes[fixgenes$Gene %in% dupgenes$Gene,] # correct locations
  
  df <- all_intersect_genes[[tt]] # get tt df, fix genes
  df$keep <- ifelse(df$Gene %in% fixgenes_filt$Gene & df$Cytoband %in% fixgenes_filt$Cytoband, "keep", 
                    ifelse(df$Gene %in% fixgenes_filt$Gene & !df$Cytoband %in% fixgenes_filt$Cytoband, "drop", "keep"))
  all_intersect_genes[[tt]] <- df[df$keep == "keep",] %>% dplyr::select(-keep)
}

# lapply(all_intersect_genes, check_dup_genes) # check that all dups removed

# plotdfs <- plotlists <- list()
plotdf_all <- data.frame()
for (tt in tumors) {
  all_intersect_genes[[tt]]$Chromosome <- factor(all_intersect_genes[[tt]]$Chromosome, levels = c(paste0("Chr", 1:22), "ChrX", "ChrY") )
  all_intersect_genes[[tt]]$Arm <- factor(all_intersect_genes[[tt]]$Arm, levels = c("p", "q"))
  all_intersect_genes[[tt]]$ChrArm <- paste0(all_intersect_genes[[tt]]$Chromosome, all_intersect_genes[[tt]]$Arm)
  
  all_intersect_genes[[tt]]$Cytoband <- ifelse(all_intersect_genes[[tt]]$Cytoband %in% "NANA", NA, all_intersect_genes[[tt]]$Cytoband)
  all_intersect_genes[[tt]]$ChrArm <- ifelse(all_intersect_genes[[tt]]$ChrArm %in% "NANA", NA, all_intersect_genes[[tt]]$ChrArm)
  
  all_intersect_genes[[tt]] <- all_intersect_genes[[tt]][order(all_intersect_genes[[tt]]$Chromosome, all_intersect_genes[[tt]]$Arm),]
  
  all_intersect_genes[[tt]]$ChrArm <- factor(gsub("Chr", "", all_intersect_genes[[tt]]$ChrArm), levels = chrarm_order)
  all_intersect_genes[[tt]]$Intersection <- factor(all_intersect_genes[[tt]]$Intersection, levels = c("Both", "ML_Only", "Reg_Only"))
  
  ### create table to plot
  # plotdfs[[tt]] <- data.frame(table(all_intersect_genes[[tt]]$Intersection, all_intersect_genes[[tt]]$ChrArm))
  # plotlists[[tt]] <- ggbarplot(plotdfs[[tt]], x = "Var2", y = "Freq", fill = "Var1", xlab = "", ylab = "Count", 
  #                              palette = c("darkgreen", "yellow", "blue"), title = tt) 
  plotdf <- cbind(data.frame(table(all_intersect_genes[[tt]]$Intersection, all_intersect_genes[[tt]]$ChrArm)), TumorType=rep(tt))
  plotdf_all <- rbind(plotdf_all, plotdf)
}

# pdf("ML_Models/ML_0.7_0.3_predictionthresh_0.2delfreq_filt_Reg_Models_Genes_Intersect.pdf", width = 22, height = 12)
# ggarrange(plotlist = plotlists, common.legend = T, ncol=2, nrow=5)
# dev.off()

pdf("ML_Models/MLmodels_meanthresh_delfreq.2_RegModels_Genes_Intersect.pdf", width = 22, height = 12)
ggbarplot(plotdf_all, x = "Var2", y = "Freq", fill = "Var1", color = "Var1", 
          facet.by = "TumorType", xlab = "", ylab = "Count", 
          palette = c("darkgreen", "orange", "blue")) + theme_bw()
dev.off()





### fisher
# ML hits VS regression hits (losses associated with immune cold)
ftestout <- data.frame()
for (tt in tumors) {
  test <- all_intersect_genes[[tt]]
  cont <- table(factor(test$Intersection, levels = c("Both", "ML_Only", "Reg_Only") ))
  
  
  ct <-  matrix(c(cont[["Both"]], cont[["Reg_Only"]], cont[["ML_Only"]], length(setdiff(ml_outputs[[tt]]$Gene, test$Gene))),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(c("In_Reg", "Not_In_Reg"), c("In_ML", "Not_In_ML")))
  ftest <- fisher.test(ct, alternative = "greater")
  
  ftestout <- rbind(ftestout, data.frame(TumorType=tt, NumGenes_Both = ct[1,1], NumGenes_ML_only = ct[2,1], 
                                         NumGenes_Reg_only = ct[1,2], NumGenes_Neither = ct[2,2], 
                                         OddsRatio = unname(ftest$estimate), Pval = unname(ftest$p.value)))
}

write.csv(ftestout, "ML_Models/FisherTest_ML_RegLosses_ImmuneCOLD_Greater.csv")



### immune hot genes (losses) and ML model genes
ML_pos_Reg_Loss_Hot <- list()
for (tt in tumors) {
  mlgenes <- rbind(control_genesdf, predicted_genes[[tt]])
  reggenes <- reg_genes_hot[reg_genes_hot$TumorType %in% tt,]
  
  in_both <- intersect(mlgenes$Gene, reggenes$Gene) # intersection
  only_in_ML <- setdiff(mlgenes$Gene, reggenes$Gene) # only in ML output
  only_in_reg <- setdiff(reggenes$Gene, mlgenes$Gene) # only in reg output
  
  ML_pos_Reg_Loss_Hot[[tt]] <- data.frame(Both = length(in_both), ML_Only = length(only_in_ML), Reg_Only = length(only_in_reg))
}

ftestout_hot <- data.frame()
for (tt in tumors) {
  cont = ML_pos_Reg_Loss_Hot[[tt]]
  ct <-  matrix(c(cont[["Both"]], cont[["Reg_Only"]], cont[["ML_Only"]], length(setdiff(ml_outputs[[tt]]$Gene, test$Gene))),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(c("In_Reg", "Not_In_Reg"), c("In_ML", "Not_In_ML")))
  ftest <- fisher.test(ct, alternative = "greater")
  
  ftestout_hot <- rbind(ftestout_hot, data.frame(TumorType=tt, NumGenes_Both = ct[1,1], NumGenes_ML_only = ct[2,1], 
                                         NumGenes_Reg_only = ct[1,2], NumGenes_Neither = ct[2,2], 
                                         OddsRatio = unname(ftest$estimate), Pval = unname(ftest$p.value)))
}

write.csv(ftestout_hot, "ML_Models/FisherTest_ML_RegLosses_ImmuneHOT_Greater.csv")



#### compare all gene outputs (ML, reg, manual filt)
ML_genes <- predicted_genes 

Reg_genes <- list()
for (tt in unique(reg_genes_filt$TumorType)) {
  Reg_genes[[tt]] <- reg_genes_filt[reg_genes_filt$TumorType %in% tt,]
}

Manual_Filt_genes <- list()
for (tt in unique(reg_genes_filt$TumorType)) {
  if (file.exists(paste0("TCGA_ManualFiltering_Outputs/pval0.01/", tt, "_filtered_summary_table.csv"))) {
    test <- read.csv(paste0("TCGA_ManualFiltering_Outputs/pval0.01/", tt, "_filtered_summary_table.csv"))
    test <- test[test$Filtered_Hit %in% "Deletion",]
  } else {
    test <- data.frame()
  }
  Manual_Filt_genes[[tt]] <- test
}


compare_genes <- list()
for (tt in unique(reg_genes_filt$TumorType)) {
  mlgenes <- ML_genes[[tt]]$Gene
  reggenes <- Reg_genes[[tt]]$Gene
  manfiltgenes <- Manual_Filt_genes[[tt]]$Gene
  
  compare_genes[[tt]] <- list(ML=mlgenes, Reg=reggenes, ManFilt=manfiltgenes)
}


### compare all three
intersect3 <- list()
for (tt in names(compare_genes)) {
  print(tt)
  if (is.null(Reduce(intersect,compare_genes[[tt]])) | is_empty(Reduce(intersect,compare_genes[[tt]])) ) {
    intersect3[[tt]] <- c(NA)
  } else {
    intersect3[[tt]] <- Reduce(intersect,compare_genes[[tt]])
  }
  print(length(intersect3[[tt]]))
}


summary_tables_intersect3 <- data.frame()
for (tt in names(intersect3)) {
  if (length(intersect3[[tt]]) == 1) {
    if (is.na(intersect3[[tt]])) {
      test <- data.frame(Gene = NA, Cytoband = NA, TumorType = tt)
    } else {
      test <- Manual_Filt_genes[[tt]]
      test <- cbind(test[test$Gene %in% intersect3[[tt]] , c("Gene", "Cytoband")] %>% remove_rownames(), TumorType = rep(tt))
    }
  } else {
    test <- Manual_Filt_genes[[tt]]
    test <- cbind(test[test$Gene %in% intersect3[[tt]] , c("Gene", "Cytoband")] %>% remove_rownames(), TumorType = rep(tt))
  }
  summary_tables_intersect3 <- rbind(summary_tables_intersect3, test)
}
summary_tables_intersect3$Chromosome <- as.numeric(gsub("p.*|q.*", "", summary_tables_intersect3$Cytoband))
summary_tables_intersect3$Arm <- ifelse(grepl("p", summary_tables_intersect3$Cytoband), "p", "q")
summary_tables_intersect3 <- summary_tables_intersect3[order(summary_tables_intersect3$TumorType, 
                                                             summary_tables_intersect3$Chromosome, summary_tables_intersect3$Arm),]
summary_tables_intersect3$ChrArm <- factor(paste0(summary_tables_intersect3$Chromosome, summary_tables_intersect3$Arm), levels = chrarm_order)

nohits <- summary_tables_intersect3[is.na(summary_tables_intersect3$Gene),"TumorType"]
# nohits

pdf("Compare3_ML_REG_MANFILT.pdf", width = 8, height = 5)
ggplot(na.omit(summary_tables_intersect3), aes(x=ChrArm)) + geom_bar() + facet_wrap(~TumorType)
dev.off()

write.csv(summary_tables_intersect3, "Compare3_ML_REG_MANFILT.csv", row.names = F)



