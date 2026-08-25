source("R/manual_omic_immunescore_corr_filtering.R") # filtering function
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
cytobandorder <- read.csv("cytobandorder.csv")

#### data
summary_tables <- readRDS("SummaryTables/TCGA_summary_tables.RDS") # RDS object with all summary tables in list
cancers=names(summary_tables) # all cancer types


#### filtering thresholds

filter_del = list("corr_DNA_RNA > 0 & corr_DNA_IS > 0 & corr_RNA_IS > 0", 
                  "corr_DNA_RNA_Pval < 0.01 & corr_DNA_IS_Pval < 0.01 & corr_RNA_IS_Pval < 0.01", # p val signif
                  # "corr_DNA_RNA_FDR < 0.5 & corr_DNA_IS_FDR < 0.5 & corr_RNA_IS_FDR < 0.5", # fdr thresh
                  "freq_del_0.2 > 0.2") # high del freq

filter_amp = list("corr_DNA_RNA > 0 & corr_DNA_IS < 0 & corr_RNA_IS < 0",
                  "corr_DNA_RNA_Pval < 0.01 & corr_DNA_IS_Pval < 0.01 & corr_RNA_IS_Pval < 0.01", # signif p val
                  # "corr_DNA_RNA_FDR < 0.5 & corr_DNA_IS_FDR < 0.5 & corr_RNA_IS_FDR < 0.5", # fdr thresh
                  "freq_amp_0.2 > 0.2") # high amp freq


### run for all cancers
for (tt in cancers) {
  print(tt)
  df=summary_tables[[tt]]
  manualfilt <- custom_filter(df = df, filter_del = filter_del, filter_amp = filter_amp,
                              plottitle = tt, savepath = paste0("TCGA_ManualFiltering_Outputs/pval0.01/", tt))
}


### compare numbers across tt and different filtering runs

# pval < 0.05
pvalonly_genecounts <- data.frame()
for (tt in cancers) {
  df <- read.csv(paste0("TCGA_ManualFiltering_Outputs/pval0.05/", tt, "_filtered_gene_counts.csv"))
  pvalonly_genecounts <- rbind(pvalonly_genecounts, cbind(TumorType=tt, df[nrow(df),c("Deletions_NumberGenes", "Amplifications_NumberGenes")]))
}


# pval < 0.01
pvalonly_low_genecounts <- data.frame()
for (tt in cancers) {
  df <- read.csv(paste0("TCGA_ManualFiltering_Outputs/pval0.01/", tt, "_filtered_gene_counts.csv"))
  pvalonly_low_genecounts <- rbind(pvalonly_low_genecounts, cbind(TumorType=tt, df[nrow(df),c("Deletions_NumberGenes", "Amplifications_NumberGenes")]))
}

# pval < 0.05, fdr < 0.5
pval_fdr_high_genecounts <- data.frame()
for (tt in cancers) {
  df <- read.csv(paste0("TCGA_ManualFiltering_Outputs/pval0.05_fdr0.5/", tt, "_filtered_gene_counts.csv"))
  pval_fdr_high_genecounts <- rbind(pval_fdr_high_genecounts, cbind(TumorType=tt, df[nrow(df),c("Deletions_NumberGenes", "Amplifications_NumberGenes")]))
}

# pval < 0.05, fdr < 0.2
pval_fdr_low_genecounts <- data.frame()
for (tt in cancers) {
  df <- read.csv(paste0("TCGA_ManualFiltering_Outputs/pval0.05_fdr0.2/", tt, "_filtered_gene_counts.csv"))
  pval_fdr_low_genecounts <- rbind(pval_fdr_low_genecounts, cbind(TumorType=tt, df[nrow(df),c("Deletions_NumberGenes", "Amplifications_NumberGenes")]))
}



##### plot/compare
pvalonly_genecounts_plot <- pvalonly_genecounts %>% 
  pivot_longer(cols = `Deletions_NumberGenes`:`Amplifications_NumberGenes`, names_to = "Event", values_to = "Count")
pvalonly_genecounts_plot$Event <- gsub("_NumberGenes", "", pvalonly_genecounts_plot$Event)
p1 <- ggbarplot(pvalonly_genecounts_plot, x="TumorType", y="Count", fill="Event", palette = c("red", "blue"), 
                title = "pval < 0.05", xlab = "", ylab = "") + coord_flip() + theme_bw()


pvalonly_low_genecounts_plot <- pvalonly_low_genecounts %>% 
  pivot_longer(cols = `Deletions_NumberGenes`:`Amplifications_NumberGenes`, names_to = "Event", values_to = "Count")
pvalonly_low_genecounts_plot$Event <- gsub("_NumberGenes", "", pvalonly_low_genecounts_plot$Event)
p2 <- ggbarplot(pvalonly_low_genecounts_plot, x="TumorType", y="Count", fill="Event", palette = c("red", "blue"), 
                title = "pval < 0.01", xlab = "", ylab = "") + coord_flip() + theme_bw()


pval_fdr_high_genecounts_plot <- pval_fdr_high_genecounts %>% 
  pivot_longer(cols = `Deletions_NumberGenes`:`Amplifications_NumberGenes`, names_to = "Event", values_to = "Count")
pval_fdr_high_genecounts_plot$Event <- gsub("_NumberGenes", "", pval_fdr_high_genecounts_plot$Event)
p3 <- ggbarplot(pval_fdr_high_genecounts_plot, x="TumorType", y="Count", fill="Event", palette = c("red", "blue"), 
                title = "pval < 0.05, FDR < 0.5", xlab = "", ylab = "") + coord_flip() + theme_bw()


pval_fdr_low_genecounts_plot <- pval_fdr_low_genecounts %>% 
  pivot_longer(cols = `Deletions_NumberGenes`:`Amplifications_NumberGenes`, names_to = "Event", values_to = "Count")
pval_fdr_low_genecounts_plot$Event <- gsub("_NumberGenes", "", pval_fdr_low_genecounts_plot$Event)
p4 <- ggbarplot(pval_fdr_low_genecounts_plot, x="TumorType", y="Count", fill="Event", palette = c("red", "blue"), 
                title = "pval < 0.05, FDR < 0.2", xlab = "", ylab = "") + coord_flip() + theme_bw()


pdf("gene_counts_filtering.pdf", width = 20, height = 10)
ggarrange(p1, p2, p3, p4, nrow=1, common.legend = T)
dev.off()



plotdf <- rbind(cbind(pvalonly_genecounts_plot, Filt=rep("pval < 0.05")),
                cbind(pvalonly_low_genecounts_plot, Filt=rep("pval < 0.01")),
                cbind(pval_fdr_high_genecounts_plot, Filt=rep("pval < 0.05, FDR < 0.5")),
                cbind(pval_fdr_low_genecounts_plot, Filt=rep("pval < 0.05, FDR < 0.1")))
plotdf$Filt <- factor(plotdf$Filt, levels = c("pval < 0.05", "pval < 0.01", "pval < 0.05, FDR < 0.5", "pval < 0.05, FDR < 0.1"))


p5 <- facet(ggbarplot(plotdf, x="TumorType", y="Count", fill="Event", palette = c("red", "blue"), xlab = "", ylab = ""),
            facet.by = "Filt", nrow = 1) + coord_flip() + theme_bw()

pdf("gene_counts_filtering.pdf", width = 14, height = 8)
p5
dev.off()




#### intersect top genes: p < 0.05 & fdr < 0.2
cancers <- gsub("_.*", "", list.files("TCGA_ManualFiltering_Outputs/pval0.05_fdr0.2/", pattern = "filtered_summary_table"))
remove <- c("BRCA","CHOL","COAD", "COAD.MSI", "COAD.MSS", "COADREAD", "DLBC", "ESCA", 
            "HNSC", "HNSC.TCGAdef.HPVneg", "HNSC.TCGAdef.HPVpos", "KIPAN", "READ", "READ.MSI", "READ.MSS", "STES", "UCEC")
cancers <- cancers[!cancers %in% remove]

gain_hits <- list()
loss_hits <- list()
genelocations <- data.frame()
for (tt in cancers) {
  df <- read.csv(paste0("TCGA_ManualFiltering_Outputs/pval0.05_fdr0.2/", tt, "_filtered_summary_table.csv"), check.names = F)
  
  genelocations <- unique(rbind(genelocations, df[,c("Gene", "Cytoband")]))
  
  gain_hits[[tt]] <- df[df$Filtered_Hit %in% "Amplification", "Gene"]
  loss_hits[[tt]] <- df[df$Filtered_Hit %in% "Deletion", "Gene"]
}
genelocations$Cytoband <- factor(genelocations$Cytoband, levels = cytobandorder$Cytoband)

### count how many tumor types with same genes
ampgenes <- enframe(gain_hits, name = "TumorTypes", value = "Gene") %>% unnest(Gene)
ampgene_counts <- ampgenes %>% group_by(Gene) %>% summarize(TumorTypes = paste(unique(TumorTypes), collapse = ", "), count = n())
ampgene_counts <- merge(genelocations, ampgene_counts, by="Gene", all.y=T)
ampgene_counts <- ampgene_counts[order(ampgene_counts$count, decreasing = T),]
write.csv(ampgene_counts, "TCGA_ManualFiltering_Outputs/top_amplified_genes_p<0.05_fdr<0.2.csv", row.names = F)


delgenes <- enframe(loss_hits, name = "TumorTypes", value = "Gene") %>% unnest(Gene)
delgenes_counts <- delgenes %>% group_by(Gene) %>% summarize(TumorTypes = paste(unique(TumorTypes), collapse = ", "), count = n())
delgenes_counts <- merge(genelocations, delgenes_counts, by="Gene", all.y=T)
delgenes_counts <- delgenes_counts[order(delgenes_counts$count, decreasing = T),]
write.csv(delgenes_counts, "TCGA_ManualFiltering_Outputs/top_deleted_genes_p<0.05_fdr<0.2.csv", row.names = F)




#### pval <0.05 plot
ampgenes <- read.csv("TCGA_ManualFiltering_Outputs/top_amplified_genes_p<0.05.csv", check.names = F)
delgenes <- read.csv("TCGA_ManualFiltering_Outputs/top_deleted_genes_p<0.05.csv", check.names = F)





