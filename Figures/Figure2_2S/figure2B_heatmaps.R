source("../../R/get_top_regression_hits.R")
library(readxl)
library(openxlsx)

workingtables <- readRDS("data/TCGA_working_tables.rds")

cyto_linreg <- read.csv("IS_CNV_Regression/TCGA/linear_regression_Outputs_cytobandlvl.csv")
cyto_linreg$X <- NULL

tumortypes <- unique(cyto_linreg$TumorType)
all_cyto_linreg_results <- list()
for (tt in tumortypes) {
  d = cyto_linreg[cyto_linreg$TumorType %in% tt, ]
  d$FDR = p.adjust(d$P_val, method = "fdr") # add FDR to table
  d$pval_signif <- ifelse(d$P_val < 0.05, "Yes", "No") # pval less than 0.05?
  d$plot_hit <- d$`X.signedlog.P.adj.` * d$cnv_status # what value to plot in heatmap (negative for loss and positive for gain)
  d$plot_hit <- ifelse(d$P_val < 0.05, d$plot_hit, 0) # coerce -log(FDR) to 0 if pval not signif
  d <- d[,c("TumorType", "Chromosome", "Arm", "ChrArm", "Cytoband", "Coefficients", "T_val", "P_val", "pval_signif", "FDR", "X.signedlog.P.adj.", "cnv_status", "plot_hit")]
  d <- d %>% remove_rownames()
  all_cyto_linreg_results[[tt]] <- d
}


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



### cytobands Aneuploidy Score
cytoAS <- data.frame() #### dataframe with cytoband level aneuploidy scores (sum of all absolute values of cytoband level cnv's)
for (tt in tumor_order) {
  cytodf <- workingtables[[paste0(tt, "_CytoCNV")]]
  cytoAS <- rbind(cytoAS, data.frame(TumorType = rep(tt), Cytoband_AS = rowSums(abs(cytodf))))
}
cytoAS_mean <- cytoAS %>% group_by(TumorType) %>% summarize(MeanAS = mean(Cytoband_AS))

### top cytoband hits and heatmaps
get_top_linearreg(reg_out=all_cyto_linreg_results, AS_df = cytoAS_mean, remove = remove, cnvlevel = "Cytoband", 
                  continuous = T, tumor_order=tumor_order, tumor_groups=tumor_groups, tumor_group_cols=tumor_group_cols,
                  cluster_rows=F, cluster_columns=F, show_column_names=F, show_row_names=T, 
                  savepath = "IS_CNV_Regression/TCGA/Linear_Regression_TopHits/CytobandCNV_linreg_tophits", savetable = T) 


full_cyto_df <- data.frame()
for(tt in names(all_cyto_linreg_results)) {
  d <- all_cyto_linreg_results[[tt]]
  full_cyto_df <- rbind(full_cyto_df, d)
}
full_cyto_df$ImmuneHit <- ifelse(full_cyto_df$plot_hit < 0, "Immune Cold", ifelse(full_cyto_df$plot_hit > 0, "Immune Hot", ""))
write.csv(full_cyto_df, "IS_CNV_Regression/TCGA/Linear_Regression_TopHits/CytoCNV_linreg_results.csv", row.names = F)
