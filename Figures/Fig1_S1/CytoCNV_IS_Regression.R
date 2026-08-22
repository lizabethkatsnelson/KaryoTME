source("R/immunescore_cnv_regression.R") # function to run regression models
source("R/plot_area_immunescore_cnv_regression.R") # function to plot outputs

################################ get cytoband tables ################################
### read in tcga dna, rna, IS tables
working_tables <- readRDS(file = "data/TCGA_working_tables.rds")

### subtypes
tcgasubtypes <- read.csv("../TCGA_Cox/ABSOLUTE.purity.txt", sep="\t") # table with cancer names and subtypes - all samples in each group
tcgasubtypes <- tcgasubtypes[!tcgasubtypes$cancer %in% "GBMLGG",] # remove gbm lgg
cnv_gene_files_path <- "../TCGA_ABSOLUTE_purity_ploidy_corrected/" # path to gene cnv tables
cancers <- unique(gsub("_.*", "", list.files(cnv_gene_files_path))) # list of cancer names

for (ct in cancers) {
  ### read in gene file - get cytoband table
  cnv_file <- paste0(cnv_gene_files_path, ct, "_2020-06-11_purity_rescale_all_data_by_genes.txt") # gene cnv file names
  cnvdf <- read.delim(cnv_file, check.names = F) # read in gene cnv table
  cytobands = unique(cnvdf$Cytoband) # cytoband order

  cytobanddf <- cnvdf[,-c(1:2)] %>% group_by(Cytoband) %>% summarise(across(everything(), mean)) # aggregate by mean cytoband per sample
  cytobanddf$Cytoband <- factor(cytobanddf$Cytoband, levels = cytobands) # factor cytoband column
  cytobanddf <- as.data.frame(t(cytobanddf[order(cytobanddf$Cytoband),] %>% column_to_rownames(var="Cytoband"))) # order by bands, cols=cyto, rows=samples

  working_tables[[paste0(ct, "_CytoCNV")]] <- cytobanddf # add cytoband cnv table to working tables object

  possible_subtypes <- unique(grep(paste0(ct, "."), tcgasubtypes$cancer, value = T, fixed = T)) # get subtypes
  if (ct=="READ") { possible_subtypes <- possible_subtypes[!grepl("COAD", possible_subtypes)] }
  if ( length(possible_subtypes) > 1 ) {
    for (sub in possible_subtypes) {
      cytobanddf_sub <- cytobanddf[rownames(cytobanddf) %in% tcgasubtypes[tcgasubtypes$cancer == sub, "sample"],]  # filter samples for subtype
      working_tables[[paste0(sub, "_CytoCNV")]] <- cytobanddf_sub # add to tables
    }
  }
}
working_tables <- working_tables[order(names(working_tables))] # order tables by cancer type
saveRDS(working_tables, file = "data/TCGA_working_tables.rds")



################################ REGRESSIONS ################################
## read in tables
working_tables <- readRDS(file = "data/TCGA_working_tables.rds")
cancers <- unique(gsub("_.*", "", names(working_tables)))

outpath = "IS_CNV_Regression/TCGA/Cytoband/"

all_outputs_bi <- list()
all_outputs_con <- list()

for (ct in cancers) {
  print(ct)

  cnv_df = working_tables[[paste0(ct, "_CytoCNV")]]
  cnv_df = cbind(AS = rowSums(abs(cnv_df)), cnv_df) # add Aneuploidy Score (AS)

  is_df = working_tables[[paste0(ct, "_IS")]][,"Ranked_Sum",drop=F]

  ### cnv binary
  reg_out_bi <- immunescore_cnv_regression(cnv_df = cnv_df, is_df = is_df, use_AS = T, cnvlevel = "Cytoband",
                                           cnv_thresh = 0.2, is_thresh = 0.3, CNV_binary = T, IS_binary = T,
                                           model = "log", savepath = paste0(outpath, "Outputs/", ct, "_is0.3_cytocnv_binary0.2_regression.csv"))
  all_outputs_bi[[ct]] <- reg_out_bi

  print("done with first regression")

  ### cnv continuous
  reg_out_con <- immunescore_cnv_regression(cnv_df = cnv_df, is_df = is_df, use_AS = T, cnvlevel = "Cytoband",
                                            cnv_thresh = 0.2, is_thresh = 0.3, CNV_binary = F, IS_binary = T,
                                            model = "log", savepath = paste0(outpath, "Outputs/", ct, "_is0.3_cytobandcnv_continuous0.2_regression.csv"))
  all_outputs_con[[ct]] <- reg_out_con

  print("done with second regression")
}

### save object
saveRDS(all_outputs_bi, file = "IS_CNV_Regression/TCGA/Cytoband/Regression_Outputs_is0.3_cytobandcnv_binary0.2.rds")
saveRDS(all_outputs_con, file = "IS_CNV_Regression/TCGA/Cytoband/Regression_Outputs_is0.3_cytobandcnv_continuous0.2.rds")

print("saved rds files")



############################## plotting ##############################

all_outputs_bi <- readRDS("IS_CNV_Regression/TCGA/Cytoband/Regression_Outputs_is0.3_cytobandcnv_binary0.2.rds")
all_outputs_con <- readRDS("IS_CNV_Regression/TCGA/Cytoband/Regression_Outputs_is0.3_cytobandcnv_continuous0.2.rds")

### plot binary outputs
print("Binary plots")
for (ct in names(all_outputs_bi)) {
  print(ct)
  
  # binary cnv, z val
  plot_cnv_IS_regression(df = all_outputs_bi[[ct]], cnvlevel = "Cytoband", scorename = "IS", stat = "Z_val", cnvtype = "bi", 
                         signif = 2, plottitle = "IS(bi0.3) ~ Cytoband_CNV(bi0.2) + AS", 
                         savepath = paste0("IS_CNV_Regression/TCGA/Cytoband/Plots/", ct, "_is0.3_cytobandcnv_binary0.2_reg_plot_Zval.pdf"))
  
  # binary cnv, p val
  plot_cnv_IS_regression(df = all_outputs_bi[[ct]], cnvlevel = "Cytoband", scorename = "IS", stat = "P_val", cnvtype = "bi", 
                         signif = 2, plottitle = "IS(bi0.3) ~ Cytoband_CNV(bi0.2) + AS", 
                         savepath = paste0("IS_CNV_Regression/TCGA/Cytoband/Plots/", ct, "_is0.3_cytobandcnv_binary0.2_reg_plot_logP.pdf"))
}
print("Done binary plots")



################################ PLOTTING ################################
### plot continuous outputs
print("Continuous plots")
for (ct in names(all_outputs_con)) {
  print(ct)
  
  # continuous cnv, z val
  plot_cnv_IS_regression(df = all_outputs_con[[ct]], cnvlevel = "Cytoband", scorename = "IS", stat = "Z_val", cnvtype = "con", 
                         signif = 2, plottitle = "IS(bi0.3) ~ Cytoband_CNV(con0.2) + AS", 
                         savepath = paste0("IS_CNV_Regression/TCGA/Cytoband/Plots/", ct, "_is0.3_cytobandcnv_continuous0.2_reg_plot_Zval.pdf"))
  # continuous cnv, p val
  plot_cnv_IS_regression(df = all_outputs_con[[ct]], cnvlevel = "Cytoband", scorename = "IS", stat = "P_val", cnvtype = "con", 
                         signif = 2, plottitle = "IS(bi0.3) ~ Cytoband_CNV(con0.2) + AS", 
                         savepath = paste0("IS_CNV_Regression/TCGA/Cytoband/Plots/", ct, "_is0.3_cytobandcnv_continuous0.2_reg_plot_logP.pdf"))
}
print("Done continuous plots")

