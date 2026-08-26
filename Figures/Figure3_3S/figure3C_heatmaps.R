source("../../R/get_top_regression_hits.R")
library(ggpubr)

workingtables <- readRDS("data/TCGA_working_tables.rds")

arms_linreg <- read.csv("IS_CNV_Regression/TCGA/linear_regression_Outputs_armlvl.csv")
arms_linreg$X <- NULL

tumortypes <- unique(arms_linreg$TumorType)
all_arm_linreg_results <- list()
for (tt in tumortypes) {
  d = arms_linreg[arms_linreg$TumorType %in% tt, ]
  d$FDR = p.adjust(d$P_val, method = "fdr") # add FDR to table
  d$pval_signif <- ifelse(d$P_val < 0.05, "Yes", "No") # pval less than 0.05?
  d$plot_hit <- d$`X.signedlog.P.adj.` * d$cnv_status # what value to plot in heatmap (negative for loss and positive for gain)
  d$plot_hit <- ifelse(d$P_val < 0.05, d$plot_hit, 0) # coerce -log(FDR) to 0 if pval not signif
  d <- d[,c("TumorType", "Arm", "Coefficients", "T_val", "P_val", "pval_signif", "FDR", "X.signedlog.P.adj.", "cnv_status", "plot_hit")]
  d <- d %>% remove_rownames()
  all_arm_linreg_results[[tt]] <- d
}
# View(all_arm_linreg_results$UCEC.MSI)

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

### arms Aneuploidy Score
armAS <- data.frame() #### dataframe with arm level aneuploidy scores (sum of all absolute values of arm level cnv's)
for (tt in tumor_order) {
  armdf <- workingtables[[paste0(tt, "_ArmCNV")]]
  armAS <- rbind(armAS, data.frame(TumorType = rep(tt), Arm_AS = rowSums(abs(armdf))))
}
armAS_mean <- armAS %>% group_by(TumorType) %>% summarize(MeanAS = mean(Arm_AS))


### xcell arm events linear regression
xcell_arm <- read_xlsx("Deconvolutions_CNV_Regressions/TCGA/TableS3_TCGA_xCell_LinRegressions_Tables.xlsx", sheet = 2)

tumortypes <- unique(xcell_arm$TumorType)
celltypes <- unique(xcell_arm$CellType)
all_xcell_arm_linreg_results <- list()
for (ct in celltypes) {
  d_cell = xcell_arm[xcell_arm$CellType %in% ct, ]
  for (tt in tumortypes) {
    d_cell_t <- d_cell[d_cell$TumorType %in% tt, ]
    # d_cell_t$FDR = p.adjust(d_cell_t$P_val, method = "fdr") # add FDR to table
    d_cell_t$pval_signif <- ifelse(d_cell_t$P_val < 0.05, "Yes", "No") # pval less than 0.05?
    d_cell_t$plot_hit <- d_cell_t$`-signedlog(P.adj)` * d_cell_t$cnv_status # what value to plot in heatmap (negative for loss and positive for gain)
    d_cell_t$plot_hit <- ifelse(d_cell_t$P_val < 0.05, d_cell_t$plot_hit, 0) # coerce -log(FDR) to 0 if pval not signif
    d_cell_t <- d_cell_t[,c("TumorType", "Arm", "Coefficients", "T_val", "P_val", "pval_signif", "P.adj", "-signedlog(P.adj)", "cnv_status", "plot_hit")]
    d_cell_t <- d_cell_t %>% remove_rownames()
    all_xcell_arm_linreg_results[[ct]][[tt]] <- d_cell_t
  }
}


for (ct in celltypes) {
  get_top_linearreg(reg_out=all_xcell_arm_linreg_results[[ct]], AS_df = armAS_mean, remove = remove, cnvlevel = "Arm", 
                    continuous = T, tumor_order=tumor_order, tumor_groups=tumor_groups, tumor_group_cols=tumor_group_cols,
                    cluster_rows=F, cluster_columns=F, show_column_names=F, show_row_names=T, 
                    savepath = paste0("Deconvolutions_CNV_Regressions/TCGA/Linear_Regression/xCell_", ct, "_ArmCNV_linreg_tophits"), savetable = T) 
}