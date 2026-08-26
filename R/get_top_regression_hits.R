library(tidyr)
library(tidyverse)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)

cytoband.df <- read.csv("cytobandorder.csv")
cytobandlevels <- cytoband.df$Cytoband
armlevels <- c(paste0(rep(1:12, each=2), c("p","q")), 
               paste0(13:15, "q"), paste0(rep(16:20, each=2), c("p","q")), 
               paste0(21:22, "q"), "Xp", "Xq")



##### function to find top arm/cytoband hits and plot heatmap
get_top_linearreg <- function(reg_out, AS_df = NULL, remove = NULL, cnvlevel = c("Cytoband", "Arm"), 
                        continuous = T, tumor_order=NULL, tumor_groups=NULL, tumor_group_cols=NULL,
                        cluster_rows=F, cluster_columns=F, show_column_names=F, show_row_names=T, 
                        max_val=2, savepath = "", savetable = T) {
  
  ### FILTER TUMOR TYPES
  if (!is.null(remove)) { reg_out <- reg_out[!names(reg_out) %in% remove] } # filter out tumor types if provided
  cancers = names(reg_out) # get names of tumor types
  
  ### CONCAT INTO ONE DF
  allregoutdf <- data.frame()
  for (tt in cancers) {
    allregoutdf <- rbind(allregoutdf, 
                         cbind(TumorType=rep(tt), data.frame(reg_out[[tt]][,c(cnvlevel, "cnv_status", "plot_hit")], row.names = NULL)))
  }
  allregoutdf$Event = factor(ifelse(allregoutdf$cnv_status == -1, "Loss", "Gain"), levels = c("Gain", "Loss"))
  
  if (cnvlevel == "Cytoband") {
    allregoutdf$Cytoband <- factor(allregoutdf$Cytoband, levels = cytobandlevels)
    allregoutdf <- allregoutdf[order(allregoutdf$TumorType, allregoutdf$Cytoband, allregoutdf$cnv_status),]
  }
  if (cnvlevel == "Arm") {
    allregoutdf$Arm <- factor(gsub("Arm", "", allregoutdf$Arm), levels = armlevels)
    allregoutdf <- allregoutdf[order(allregoutdf$TumorType, allregoutdf$Arm, allregoutdf$cnv_status),]
  }
  
  ### associated with immune COLD or HOT?
  allregoutdf$ImmuneCold <- ifelse(allregoutdf$plot_hit < 0, T, F)
  allregoutdf$ImmuneHot <- ifelse(allregoutdf$plot_hit > 0, T, F)
  
  allregoutdf$plot_cold_hit <- ifelse(allregoutdf$ImmuneCold == T, allregoutdf$plot_hit, 0)
  allregoutdf$plot_cold_hit <- ifelse(allregoutdf$Event == "Gain", allregoutdf$plot_cold_hit * -1, allregoutdf$plot_cold_hit)
  
  allregoutdf$plot_hot_hit <- ifelse(allregoutdf$ImmuneHot == T, allregoutdf$plot_hit, 0)
  allregoutdf$plot_hot_hit <- ifelse(allregoutdf$Event == "Loss", allregoutdf$plot_hot_hit * -1, allregoutdf$plot_hot_hit)

  ### save top hits output
  if (savetable == T) {
    write.csv(allregoutdf, paste0(savepath, "table.csv"), row.names = F)
  }
  
  ### HM COLORS
  # colors and legends
  arm_colors <- c("p" = "#E0E0E1", "q" = "#B4B4B4")
  barplot_colors <- c("Loss" = "#8c86f0", "Gain" = "#f26b6b")
  barplot_legend <- Legend(labels = c("Loss", "Gain"), title = "Sum of Events", legend_gp = gpar(fill = barplot_colors))
  # hm_colors <- c("darkblue", "white", "darkred")
  hm_colors = colorRamp2(c(-max_val, 0, max_val), c("darkblue", "white", "darkred"))
  
  ### plots

  ###################### immune cold hits
  plot_cold <- unique(data.frame(TumorType = allregoutdf$TumorType, CNVlevel = allregoutdf[,cnvlevel], Value = allregoutdf$plot_cold_hit)) %>% 
    group_by(TumorType, CNVlevel) %>% filter(Value !=0 | n() == 1) %>% ungroup()
  plot_cold_hm <- spread(plot_cold, key = CNVlevel, value = Value) %>% column_to_rownames("TumorType") 
  plot_cold_hm[is.na(plot_cold_hm)] <- 0
    
  ### coerce extreme values to -max or max
  plot_cold_hm[,1:ncol(plot_cold_hm)] <- apply(plot_cold_hm[,1:ncol(plot_cold_hm)], 2, function(x) ifelse(x > max_val, max_val, ifelse(x < -max_val, -max_val, x)))
  
  cold_nohit_tumors <- rownames(plot_cold_hm[rowSums(abs(plot_cold_hm)) == 0,])
  print(paste0("Tumor types with no immune cold hits: ", paste(cold_nohit_tumors, collapse = ", ")))
  # plot_cold_hm <- plot_cold_hm[!rownames(plot_cold_hm) %in% cold_nohit_tumors,] ### remove tumor types with zero hits
    
  if (!is.null(tumor_order)) {
    tumor_order_cold <- tumor_order[!tumor_order %in% setdiff(tumor_order, rownames(plot_cold_hm))] # remove tumors that didnt have hits
    plot_cold_hm <- plot_cold_hm[tumor_order_cold,] # reorder heatmap
    
    tumor_groups_cold <- tumor_groups[rownames(plot_cold_hm),,drop=F]
  }
    
    
  ###################### immune hot hits
  plot_hot <- unique(data.frame(TumorType = allregoutdf$TumorType, CNVlevel = allregoutdf[,cnvlevel], Value = allregoutdf$plot_hot_hit)) %>% 
    group_by(TumorType, CNVlevel) %>% filter(Value !=0 | n() == 1) %>% ungroup()
  plot_hot_hm <- spread(plot_hot, key = CNVlevel, value = Value) %>% column_to_rownames("TumorType") 
  plot_hot_hm[is.na(plot_hot_hm)] <- 0
  
  ### coerce extreme values to -max or max
  plot_hot_hm[,1:ncol(plot_hot_hm)] <- apply(plot_hot_hm[,1:ncol(plot_hot_hm)], 2, function(x) ifelse(x > max_val, max_val, ifelse(x < -max_val, -max_val, x)))
  
  hot_nohit_tumors <- rownames(plot_hot_hm[rowSums(abs(plot_hot_hm)) == 0,])
  print(paste0("Tumor types with no immune hot hits: ", paste(hot_nohit_tumors, collapse = ", ")))
  # plot_hot_hm <- plot_hot_hm[!rownames(plot_hot_hm) %in% hot_nohit_tumors,] ### remove tumor types with zero hits
    
  if (!is.null(tumor_order)) {
    tumor_order_hot <- tumor_order[!tumor_order %in% setdiff(tumor_order, rownames(plot_hot_hm))] # remove tumors that didnt have hits
    plot_hot_hm <- plot_hot_hm[tumor_order_hot,] # reorder heatmap
    
    tumor_groups_hot <- tumor_groups[rownames(plot_hot_hm),,drop=F]
  }
    
  ################### top annotations
  topannotdf <- data.frame(row.names = colnames(plot_cold_hm), 
                           Chromosome = factor(gsub("([a-z]).*", "", colnames(plot_cold_hm)), levels = unique(gsub("([a-z]).*", "", colnames(plot_cold_hm)))),
                           Arm = factor(ifelse(grepl("p", colnames(plot_cold_hm)), "p", "q"), levels = c("p", "q")))
  top_annot = HeatmapAnnotation(df = topannotdf[,"Arm",drop=F], col = list(Arm = arm_colors), 
                                border = TRUE, annotation_name_side = "left")
    
  ################### side annotations
    
  ################################# COLD #######################################
  ## right barplots (tumor sums) +optional groups
  tumorsums_cold <- data.frame(Losses = apply(plot_cold_hm, 1, function(x) sum(x < 0)), 
                               Gains = apply(plot_cold_hm, 1, function(x) sum(x > 0)))
  right_annot_cold = rowAnnotation(x = anno_barplot(tumorsums_cold, gp = gpar(fill = barplot_colors),
                                                    axis_param = list(side = "bottom", labels_rot = 0)), show_annotation_name = FALSE)
  ## left barplots (AS)
  if (!is.null(AS_df)) {
    AS_df_cold <- AS_df[AS_df$TumorType %in% rownames(plot_cold_hm),] %>% column_to_rownames("TumorType")
    AS_df_cold <- AS_df_cold[rownames(plot_cold_hm),,drop=F] # reorder if there is a tumor order
    # left_annot_cold = rowAnnotation(Mean_AS = anno_barplot(AS_df_cold, axis_param = list(side = "bottom", labels_rot = 0, direction="reverse")))
    
  # } else if (!is.null(tumor_groups_cold)) {
    # left_annot_cold = rowAnnotation(Tumor_Group = tumor_groups_cold$Group, col = list(Tumor_Group = tumor_group_cols),
    #                                 show_annotation_name = c(Tumor_Group = FALSE), annotation_name_rot = 45)
    
  # } else if (!is.null(AS_df) & !is.null(tumor_groups_cold)) {
    left_annot_cold = rowAnnotation(Mean_AS = anno_barplot(AS_df_cold, axis_param = list(side = "bottom", labels_rot = 0, direction="reverse")), 
                                    Tumor_Group = tumor_groups_cold$Group, col = list(Tumor_Group = tumor_group_cols),
                                    show_annotation_name = c(Tumor_Group = FALSE), annotation_name_rot = 45)
    
  } else { left_annot_cold = NULL }
    
    
    
  ################################# HOT #######################################
  ## right barplots (tumor sums)
  tumorsums_hot <- data.frame(Losses = apply(plot_hot_hm, 1, function(x) sum(x < 0)),
                              Gains = apply(plot_hot_hm, 1, function(x) sum(x > 0)))
  right_annot_hot = rowAnnotation(x = anno_barplot(tumorsums_hot, gp = gpar(fill = barplot_colors),
                                                   axis_param = list(side = "bottom", labels_rot = 0)), show_annotation_name = FALSE)
  ## left barplots (AS)
  if (!is.null(AS_df)) {
    AS_df_hot <- AS_df[AS_df$TumorType %in% rownames(plot_hot_hm),] %>% column_to_rownames("TumorType")
    AS_df_hot <- AS_df_hot[rownames(plot_hot_hm),,drop=F] # reorder if there is a tumor order
    # left_annot_hot = rowAnnotation(Mean_AS = anno_barplot(AS_df_hot, axis_param = list(side = "bottom", labels_rot = 0, direction="reverse")))
  #   
  # } else if (!is.null(tumor_groups_hot)) {
  #   left_annot_hot = rowAnnotation(Tumor_Group = tumor_groups_hot$Group, col = list(Tumor_Group = tumor_group_cols),
  #                                  show_annotation_name = c(Tumor_Group = FALSE), annotation_name_rot = 45)
  #   
  # } else if (!is.null(AS_df) & !is.null(tumor_groups_hot)) {
    left_annot_hot = rowAnnotation(Mean_AS = anno_barplot(AS_df_hot, axis_param = list(side = "bottom", labels_rot = 0, direction="reverse")), 
                                   Tumor_Group = tumor_groups_hot$Group, col = list(Tumor_Group = tumor_group_cols),
                                   show_annotation_name = c(Tumor_Group = FALSE), annotation_name_rot = 45)
    
  } else { left_annot_hot = NULL }
    
    
  if (cnvlevel=="Arm") {
    ## bottom barplots (event sums) - COLD
    eventsums_cold <- data.frame(Losses = apply(plot_cold_hm, 2, function(x) sum(x < 0)), 
                                 Gains = apply(plot_cold_hm, 2, function(x) sum(x > 0)))
    bottom_annot_cold = HeatmapAnnotation(x = anno_barplot(eventsums_cold, gp = gpar(fill = barplot_colors),
                                                           axis_param = list(direction = "reverse")), show_annotation_name = FALSE)
    ## bottom barplots (event sums) - HOT
    eventsums_hot <- data.frame(Losses = apply(plot_hot_hm, 2, function(x) sum(x < 0)),
                                Gains = apply(plot_hot_hm, 2, function(x) sum(x > 0)))
    bottom_annot_hot = HeatmapAnnotation(x = anno_barplot(eventsums_hot, gp = gpar(fill = barplot_colors),
                                                          axis_param = list(direction = "reverse")), show_annotation_name = FALSE)
  } else {
    bottom_annot_cold=bottom_annot_hot=NULL
  }
    
  ################## HEATMAP 
  ######################################## IMMUNE COLD HEATMAP ######################################## 
  lgd_cold = Legend(col_fun = hm_colors, title = "Immune Cold Hits\n-log10(p.adj)",
                    legend_height = unit(2.5, "cm"), title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 10),
                    at = c(-max_val, (-max_val/2), 0, (max_val/2), max_val), # add breaks
                    labels = c(paste0("< -", max_val, " Loss"), paste0(-max_val/2), "0", paste0(max_val/2), paste0("> ", max_val, " Gain") ) ) # add labels
  
  cold_hm <- Heatmap(as.matrix(plot_cold_hm), show_heatmap_legend = F, top_annotation = top_annot, #name = "Immune Cold Hits\nSigned logP",
                     col = hm_colors,
                     cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = T, row_names_side = "left",
                     column_split = topannotdf$Chromosome, column_gap = unit(0, "mm"),
                     row_split = c(1:nrow(plot_cold_hm)), row_title = NULL, row_gap = unit(0, "mm"),
                     bottom_annotation = bottom_annot_cold, right_annotation = right_annot_cold, left_annotation = left_annot_cold,
                     border = T)
    
  ######################################## IMMUNE HOT HEATMAP ######################################## 
  lgd_hot = Legend(col_fun = hm_colors, title = "Immune Hot Hits\n-log10(p.adj)",
                   legend_height = unit(2.5, "cm"), title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 10),
                   at = c(-max_val, (-max_val/2), 0, (max_val/2), max_val), # add breaks
                   labels = c(paste0("< -", max_val, " Loss"), paste0(-max_val/2), "0", paste0(max_val/2), paste0("> ", max_val, " Gain") ) ) # add labels
  
  hot_hm <- Heatmap(as.matrix(plot_hot_hm), show_heatmap_legend = F, top_annotation = top_annot, #name = "Immune Hot Hits\nSigned logP",
                    col = hm_colors, 
                    cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = T, row_names_side = "left",
                    column_split = topannotdf$Chromosome, column_gap = unit(0, "mm"), 
                    row_split = c(1:nrow(plot_hot_hm)), row_title = NULL, row_gap = unit(0, "mm"),
                    bottom_annotation = bottom_annot_hot, right_annotation = right_annot_hot, left_annotation = left_annot_hot, 
                    border = T)
  
  ### save plots
  pdf(paste0(savepath, "COLD_hm.pdf"), width = 12, height = 6)
  draw(cold_hm, annotation_legend_list = list(lgd_cold, barplot_legend))
  junk <- dev.off()
  
  pdf(paste0(savepath, "HOT_hm.pdf"), width = 12, height = 6)
  draw(hot_hm, annotation_legend_list = list(lgd_hot, barplot_legend))
  junk <- dev.off()
  
  ### save top hits output
  if (savetable == T) {
    write.csv(allregoutdf, paste0(savepath, "table.csv"), row.names = F)
    write.csv(plot_cold_hm, paste0(savepath, "COLD_heatmap_matrix.csv"), row.names = T)
    write.csv(plot_hot_hm, paste0(savepath, "HOT_heatmap_matrix.csv"), row.names = T)
  }
  
  # return(list(allregoutdf = allregoutdf, hm_df = hmdf_filt))
}