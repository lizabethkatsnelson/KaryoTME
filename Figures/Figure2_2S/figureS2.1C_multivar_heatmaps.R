library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)

workingtables <- readRDS("data/TCGA_working_tables.rds")

tumor_order = c("UCEC.MSI", "UCEC.MSS", "OV", "CESC", "LUSC", "HNSC.HPVneg", "ESCA.SC","BLCA",
                "BRCA.pos", "BRCA.neg", "PRAD", "LUAD", "LIHC", "COADREAD.MSI", "COADREAD.MSS", "ESCA.AD", "STAD",
                "PAAD", "KICH", "KIRC", "KIRP", "LGG", "GBM", "SKCM", "UVM", "ACC", "MESO", "PCPG", "SARC", "TGCT", "UCS")

############### top events per tumor type regressions
linreg_arm <- read.csv("IS_CNV_Regression/TCGA/Linear_Regression_TopHits/ArmCNV_linreg_results.csv")
linreg_arm$ImmuneHit <- ifelse(linreg_arm$TumorType %in% "LUAD" & linreg_arm$Arm %in% "Arm1q", "Immune Cold", linreg_arm$ImmuneHit) ### add 1q gain cold
linreg_arm_cold <- linreg_arm %>% filter(ImmuneHit %in% c("Immune Cold")) # only cold hits
linreg_arm_cold$Event = ifelse(linreg_arm_cold$cnv_status < 0, "Loss", "Gain")

linreg_arm_cold_fdr_top <- linreg_arm_cold[linreg_arm_cold$FDR < 0.2, ]
linreg_arm_cold_fdr_top <- rbind(linreg_arm_cold_fdr_top, linreg_arm_cold[linreg_arm_cold$TumorType %in% "LUAD" & linreg_arm_cold$Arm %in% "Arm1q",])

linreg_arm_cold_fdr_top <- linreg_arm_cold_fdr_top[linreg_arm_cold_fdr_top$TumorType %in% tumor_order, ]
linreg_arm_cold_fdr_top <- linreg_arm_cold_fdr_top[order(linreg_arm_cold_fdr_top$TumorType, linreg_arm_cold_fdr_top$Arm),]

# get regression formulas per tumor type
reg_con_hits_tt <- as.data.frame(linreg_arm_cold_fdr_top %>% group_by(TumorType) %>% 
                                   summarize(TopHits_Formula = paste(paste0(Arm, "_", Event), collapse = " + ")) %>%
                                   mutate(NumHits = str_count(TopHits_Formula, "\\+") + 1))

cancers <- reg_con_hits_tt[reg_con_hits_tt$NumHits > 1, "TumorType"] # tumor types with more than one top hit (use for multi regression)

##### RUN MULTI VAR REG
multi_reg_out <- data.frame()
for (tt in cancers) {
  print(tt)
  
  armdf <- workingtables[[paste0(tt, "_ArmCNV")]] # arm cnv df
  isdf <- workingtables[[paste0(tt, "_IS")]] # immune score df
  
  events <- strsplit(reg_con_hits_tt[reg_con_hits_tt$TumorType == tt, "TopHits_Formula"], " \\+ ")[[1]] # get string top events per tt
  eventsdf <- as.data.frame(t(as.data.frame(strsplit(events, "_"))),row.names = NA) # df of top arms and events per tt
  
  ### convert arm df to binary for models
  regdf <- armdf[,eventsdf$V1] # filter arm df for arms of interest
  for (i in 1:ncol(regdf)) { # convert arms of interest to binary depending on gain or loss
    if (eventsdf[i,2] == "Loss") { # if top event is loss for this arm, change binary for loss 
      regdf[,i] <- ifelse(regdf[,i] < -0.2, 1, 0)
    } else { # if top event is gain for this arm, change binary for gain
      regdf[,i] <- ifelse(regdf[,i] > 0.2, 1, 0)
    }
  }
  colnames(regdf) <- paste0(colnames(regdf), "_", eventsdf$V2) # fix names
  
  print("Arm event names correct: ")
  print(all(colnames(regdf) == events))
  
  ### add IS to model df
  regdf <- merge(isdf[,"Ranked_Sum", drop=F], regdf, by=0) %>% column_to_rownames("Row.names") # add immune score
  regdf$IS <- c(scale(regdf$Ranked_Sum)) # scale IS
  regdf$IS <- ifelse(regdf$IS >= unname(quantile(regdf$IS, probs = 0.7)), 1,  # change to binary
                     ifelse(regdf$IS <= unname(quantile(regdf$IS, probs = 0.3)), 0, NA))
  
  ### run multi reg model
  # regform <- as.formula(paste0("IS ~ ", reg_con_hits_tt[reg_con_hits_tt$TumorType == tt, "TopHits_Formula"])) # regression formula
  regform <- as.formula(paste0("IS ~ ", paste(grep("Arm", names(regdf), value = T), collapse = "+")))
  print(regform)
  
  glm_out <- glm(regform, data = regdf, family = binomial(link = "logit")) # run logistic model
  saveglm <- as.data.frame(summary(glm_out)$coefficients)[-1,]
  if (any(is.na(coef(glm_out)))) { # if any coef NA, add back to table
    saveglm <- merge(saveglm, data.frame(Estimate = coef(glm_out)[-1]), by=0, all=T)[,-6] %>%
      column_to_rownames("Row.names") %>% rename("Estimate" = "Estimate.x")
    saveglm <- saveglm[grep("Arm", names(regdf), value = T),] # order arms
  } 
  saveglm <- cbind(TumorType=rep(tt), saveglm) %>% rownames_to_column("Event")
  multi_reg_out <- rbind(multi_reg_out, saveglm)
}


colnames(multi_reg_out) <- c("Event", "TumorType", "Coefficient", "StdError", "Zval", "Pval")
write.csv(multi_reg_out, "IS_CNV_Regression/TCGA/Arm/MultiVar_Regression/Top_LinReg_filttumors_tophits_MULTI_REGRESSION.csv")
# write.csv(multi_reg_out, "IS_CNV_Regression/TCGA/Arm/continuous_filttumors_tophits_CoOccur_MULTI_REGRESSION.csv")


##### plot results
# multi_reg_out <- read.csv("IS_CNV_Regression/TCGA/Arm/MultiVar_Regression/continuous_filttumors_tophits_MULTI_REGRESSION.csv", row.names = 1)
multi_reg_out <- read.csv("IS_CNV_Regression/TCGA/Arm/MultiVar_Regression/Top_LinReg_filttumors_tophits_MULTI_REGRESSION.csv", row.names = 1)
multi_reg_out <- multi_reg_out[!multi_reg_out$TumorType %in% "HNSC.HPVother",]

heatmaplist <- list()
# tt=unique(multi_reg_out$TumorType)[2]
for (tt in unique(multi_reg_out$TumorType)) {
  test <- multi_reg_out[multi_reg_out$TumorType %in% tt,]
  test_mat <- as.matrix(data.frame(row.names = gsub("_", " ", gsub("Arm", "", test$Event)), Z=test$Zval))
  
  hm <- Heatmap(test_mat, col = colorRamp2(c(-4, 0, 2), c("#521a63", "white", "#61631a")), 
                show_heatmap_legend = FALSE, rect_gp = gpar(col = "black", lwd = 1), 
                show_column_names = F, show_row_names = T, cluster_rows = F, cluster_columns = F, 
                row_names_side = "left", column_title = tt, 
                heatmap_width = unit(3.5, "cm"), heatmap_height = unit(2*nrow(test_mat), "cm"),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.1f", test_mat[i, j]), x, y, gp = gpar(fontsize = 10, col="white"))
                })
  heatmaplist[[tt]] <- grid.grabExpr(draw(hm))
}

names(heatmaplist)

pdf("IS_CNV_Regression/TCGA/Arm/MultiVar_Regression/Top_LinReg__filttumors_tophits_MULTI_REGRESSION_Heatmaps.pdf", width = 28, height = 12)
grid.arrange(heatmaplist$ACC, heatmaplist$BLCA, heatmaplist$BRCA.pos, heatmaplist$COADREAD.MSI, heatmaplist$GBM,
             heatmaplist$HNSC.HPVneg, heatmaplist$KIRP, heatmaplist$LGG, heatmaplist$LUAD, heatmaplist$LUSC, 
             heatmaplist$OV, heatmaplist$PAAD, heatmaplist$SARC, heatmaplist$STAD, heatmaplist$TGCT, ncol = 15, nrow=1) 
dev.off()
