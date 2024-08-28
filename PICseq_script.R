# generate artificial doublets (AD)
library(dplyr)
#library(e1071)
library(caret)
#library(doParallel)
library(doMC)  # For parallel backend

set.seed(1234)

# input marker genes
T_markers <- intersect(T_markers,row.names(T_cell.combined[["RNA"]]$counts))
B_markers <- intersect(B_markers,row.names(B_cell.combined[["RNA"]]$counts))

# subset the matrix by top DEGs for each cell type
DEGnumber <- 5
top_gene_per_cluster.T <- T_cell.markers.all %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = DEGnumber)
top_gene_per_cluster.B <- B_cell.markers.all %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = DEGnumber)
top_gene_per_cluster <- unique(c(top_gene_per_cluster.T$gene,top_gene_per_cluster.B$gene, # top genes
                              T_markers ,B_markers)) # marker genes


t_data.top <- T_cell.combined[["RNA"]]$counts[top_gene_per_cluster,]
b_data.top <- B_cell.combined[["RNA"]]$counts[top_gene_per_cluster,]

# calculate the total count per cell
T_UMI <- colSums(T_cell.combined[["RNA"]]$counts)
B_UMI <- colSums(B_cell.combined[["RNA"]]$counts)


###### select major cell types #####
major_T <- c("0_0","0_1","1_0","1_1","2","3","8")
major_B <- c("0","1","3_1", "9")

# annotate major cell types
T_cell.combined$major <- as.character(T_cell.combined$sub.cluster)
T_cell.combined$major[!T_cell.combined$major %in% major_T] <- "others"
B_cell.combined$major <- as.character(B_cell.combined$sub.cluster)
B_cell.combined$major[!B_cell.combined$major %in% major_B] <- "others"

# generate dataframe of artificial doublets (AD)
randSize <- 5000
combined_data.ls <- list()
for (t in major_T){
  print(t)
  t_ct <- names(T_cell.combined$major[T_cell.combined$major==t]) # select cell IDs of a T cell type
  for (b in major_B){
    print(b)
    AD_label <- paste0("T",t,"_","B",b)
    print(paste0("generate artificial doublets (AD): ",AD_label))
    b_ct <- names(B_cell.combined$major[B_cell.combined$major==b]) # select cell IDs of a B cell type
    # random sampling T and B cells
    t_rand <- sample(t_ct,randSize,replace = T)
    b_rand <- sample(b_ct,randSize,replace = T)
    # get gene expression matrix
    t_data <- t_data.top[,t_rand]
    b_data <- b_data.top[,b_rand]
    # get total counts per cell for selected cells
    T_counts_per_cell <- T_UMI[t_rand]
    B_counts_per_cell <- B_UMI[b_rand]
    # average the gene expression of artificial doublets
    if (identical(rownames(t_data),rownames(b_data))){
      comined_TB <- t_data + b_data
      comined_TB <- comined_TB/2 # downsampling
      comined_TB[comined_TB < 1] <- 0
      data_AD <- log1p(comined_TB/((T_counts_per_cell + B_counts_per_cell)/2)*1e6) # log CPM of AD
      data_AD <- t(data_AD)
      row.names(data_AD) <- paste0(AD_label,"_",seq(1,nrow(data_AD)))
      data_AD <- as.data.frame(data_AD)
      data_AD$AD_label <- AD_label # add the doublet label to the dataframe
    } else {print("Error: rownames are not matching"); break}
    combined_data.ls[[AD_label]] <- data_AD
  }
}

combined_data <- purrr::reduce(combined_data.ls, bind_rows)
combined_data$AD_label <- as.factor(combined_data$AD_label)

# export the matrix
write.csv(combined_data,"TB_combined_data.csv") # normalized without scale


# export PICs for prediction Major cell types
PICs_pred <- PICs.combined[["RNA"]]$counts[colnames(combined_data)[-ncol(combined_data)],]
# calculate the total count per cell
PICs_UMI <- colSums(PICs.combined[["RNA"]]$counts)
PICs_pred <- log1p(PICs_pred/PICs_UMI*1e6) # log of CPM
PICs_pred <- as.data.frame(t(PICs_pred))
# scale the data
PICs_pred_scaled <- scale(PICs_pred)

# export the matrix for model building
# faster but bigger size in csv
write.csv(PICs_pred,"PICs_pred.csv") # normalized without scale


### annotation tranlation ###
# # A translator of cell type ID to annotation results
# PICs.combined$predicted_PICs
# 
T_cell.combined$sub.cluster.char <- as.character(T_cell.combined$sub.cluster)

# update meatadata
T_cell.combined@meta.data$cell_type <- case_when(T_cell.combined@meta.data$sub.cluster.char ==
                                                   "0_0" ~ "TFH Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "0_1" ~ "Tregs",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "1_0" ~ "CD8+ Ex/Mem T Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "1_1" ~ "CD4+ Ex/Mem T Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "1_2" ~ "CCL3+ CCL4+ Ex/Mem T Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "2" ~ "CD4+ Naive T Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "3" ~ "CD8+ Naive T Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "4" ~ "CD8+ CM T Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "5" ~ "Lars2hi Malat1hi T cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "6" ~ "Trbv3+ Naive T Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "7_0" ~ "Myeloid-like T Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "7_1" ~ "Ki67+ T Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "8" ~ "CD4+ IFITM+ CISH+ T Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "9" ~ "NKT Cells",
                                                 T_cell.combined@meta.data$sub.cluster.char ==
                                                   "10" ~ "Trbv4+ T Cells")

B_cell.combined$sub.cluster.char <- as.character(B_cell.combined$sub.cluster)
B_cell.combined@meta.data$cell_type <- case_when(B_cell.combined@meta.data$sub.cluster.char ==
                                                   "0" ~ "FOBs",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "1" ~ "MZBs",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "2" ~ "Iglv1-3+ Immature B Cells",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "3_0" ~ "Zeb2lo B Cells",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "3_1" ~ "ABCs",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "4" ~ "B1 Cells",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "5" ~ "Igkv1-135+ B cells",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "6" ~ "Igkv10-96+ B Cells",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "7" ~ "Igkv1-110+ B Cells",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "8" ~ "Srmhi B Cells",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "9" ~ "GC B Cells",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "10" ~ "MBCs",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "11" ~ "Igkv1-117+ B Cells",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "12" ~ "PBs",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "13" ~ "Myeloid-Like B Cells",
                                                 B_cell.combined@meta.data$sub.cluster.char ==
                                                   "14" ~ "Igkv1-99+ B Cells")




# import PICs prediction result
PICs_predicted_v1 <- read.csv("PICs_predict_v1.csv")
PICs_anno <- PICs_predicted %>% select(c(Cell_ID,pred))

# set the coresponding cell type names after translation
translate_T <- unique(T_cell.combined$sub.cluster.char)
names(translate_T) <- unique(T_cell.combined$cell_type)
translate_B <- unique(B_cell.combined$sub.cluster.char)
names(translate_B) <- unique(B_cell.combined$cell_type)
# convert to annotated cell types
toTranslate_name <- PICs_anno$pred
for (t in major_T){
  for (b in major_B){
    T_name <- names(translate_T[translate_T==t])
    B_name <- names(translate_B[translate_B==b])
    toTranslate_name <- gsub(paste0("T",t,"_","B",b), paste0(T_name," and ",B_name), toTranslate_name)
  }
}
PICs_anno$pred_trans <- toTranslate_name

PICs_predicted_v1 <- PICs_predicted_v1 %>% select(c(Cell_ID,pred))
PICs_predicted_v2 <- PICs_anno
sum(PICs_predicted_v1$pred == PICs_predicted_v2$pred)/length(PICs_predicted_v2$pred)
sum(PICs_predicted_v1$Cell_ID == PICs_predicted_v2$Cell_ID)

# import the metadata of seurat object
identical(row.names(PICs.combined@meta.data),PICs_anno$Cell_ID) # double-check if the order of two vectors is exactly the same
PICs.combined <- AddMetaData(
  object = PICs.combined,
  metadata = PICs_anno$pred_trans,
  col.name = 'predicted_PICs'
)

# replace with the up-to-date annotations
predicted_PICs <- PICs.combined$predicted_PICs
predicted_PICs <- gsub("CD8\\+ EM T Cells", "CD8+ Ex/Mem T Cells", predicted_PICs)
predicted_PICs <- gsub("CX4CR1_hi CD4\\+ EM T Cells", "CD4+ Ex/Mem T Cells", predicted_PICs)
predicted_PICs <- gsub("Nlrp6\\+ CD4\\+ EM T Cells", "CD4+ IFITM+ CISH+ T Cells", predicted_PICs)
predicted_PICs <- gsub("Treg", "Tregs", predicted_PICs)
predicted_PICs_bk <- PICs.combined$predicted_PICs
PICs.combined$predicted_PICs <- predicted_PICs



# change of T-B interaction during aging
PICs.combined$predicted_PICs
PICs.combined$orig.ident
PICs.combined$PICs_age <- paste0(PICs.combined$predicted_PICs,"_",PICs.combined$orig.ident)

PIC4compare <- unique(PICs.combined$predicted_PICs)
PIC4compare <- PIC4compare[!PIC4compare=="unknown"]
Idents(PICs.combined) <- "PICs_age"
# for (p in PIC4compare){
#   DEG.age <- FindMarkers(PICs.combined, assay = "SCT", ident.1 =paste0(p,"_Aged_PICS"), ident.2 =paste0(p,"_Young_PICS"),
#                          min.cells.group = 3)
#   write.csv(DEG.age, paste0(p,"_AgedvsYoung.csv"))
# }
dir.create("DEGs_PICs_aged_vs_young")
for (p in PIC4compare) {
  tryCatch({
    DEG.age <- FindMarkers(PICs.combined, assay = "SCT", ident.1 = paste0(p, "_Aged_PICS"), ident.2 = paste0(p, "_Young_PICS"), min.cells.group = 3)
    p <- gsub("/","_or_",p)
    write.csv(DEG.age, paste0("DEGs_PICs_aged_vs_young/",p, "_AgedvsYoung.csv"))
  }, error = function(e) {
    # Error handling code
    cat("An error occurred in:", p, "\nError Message:", e$message, "\n")
  })
}

Idents(PICs.combined) <- "orig.ident"
DEG.age <- FindMarkers(PICs.combined, assay = "SCT", ident.1 = "Aged_PICS", ident.2 = "Young_PICS", min.cells.group = 3)
write.csv(DEG.age, paste0("DEGs_PICs_aged_vs_young/","Overall_AgedvsYoung.csv"))

# volcano plot
# volcano plot for each PICs
library(EnhancedVolcano)
vol.plot.ls <- list()
cell4vol <- list.files("/opt/home/buckcenter.org/fwu/PICseq/DEGs_PICs_aged_vs_young_v1")
dir.create("PICs_aged_vs_young_Volcano_v1")
for (d in cell4vol){
  markers4plot <- read.csv(paste0("DEGs_PICs_aged_vs_young_v1/",d))
  # sort genes
  topgenes_u <- markers4plot %>% filter(p_val_adj<0.05 & avg_log2FC>0) %>% arrange(desc(avg_log2FC)) %>% head(n=20)
  topgenes_u <- topgenes_u$X
  topgenes_d <- markers4plot %>% filter(p_val_adj<0.05 & avg_log2FC<0) %>% arrange(avg_log2FC) %>% head(n=20)
  topgenes_d <- topgenes_d$X
  topgenes <- c(topgenes_u,topgenes_d)
  # Create a label column, label only selected genes
  markers4plot$label <- ifelse(markers4plot$X %in% topgenes, markers4plot$X, "")
  vol.plot <- EnhancedVolcano(markers4plot,
                              lab = markers4plot$label,
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              #xlim =c(-1.5,2),
                              ylab = bquote(~-Log[10] ~ italic(adj.P)),
                              legendLabels = c("NS", expression(Log[2] ~ FC), "adj.p-value", expression(adj.p - value ~ and
                                                                                                    ~ log[2] ~ FC)),
                              FCcutoff = 0.25,
                              pCutoff = 0.05,
                              labSize = 5, 
                              drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                              title = gsub("_AgedvsYoung.csv","",d),
                              subtitle = bquote(italic("Aged vs Young")))
  ggsave(plot=vol.plot,paste0("PICs_aged_vs_young_Volcano_v1/",d,"_Volcano.pdf"), width=9,height=8)
  vol.plot.ls[[d]] <- vol.plot
}
library(patchwork)
combined_plot <- wrap_plots(vol.plot.ls, ncol = floor(length(vol.plot.ls)/4))
print(combined_plot)
ggsave("PICs_aged_vs_young_Volcano_v1/Volcano_PICs_AgedvsYoung.pdf",width=25,height=30)



### heatmaps of genes and UMIs

# select a major PIC each time; one B and T combination a time

T_markers_major <- unique(c("Cd4","Il21","Tox2","Pdcd1","Maf","Bcl6","Cd4","Foxp3","Areg","Il2ra","Ikzf2",
                            "Ctla4","Cd8a","Gzmk","Eomes","Cd38","Pdcd1","Ifng","Cd4","Eomes","Cx3cr1",
                            "Gzmk","Ifng","Cd38","Cd4","Ccr7","Sell","Cd8a",
                            "Ccr7","Sell",
                            "Cd4","Nlrp6","Ifitm1","Ifitm3","Ifitm2","Cish","Cd44"))
B_markers_major <- unique(c("Fcer2a","Fchsd2","Icosl","Ms4a4c","Zfp318",
                            "Cr2","Asb2","Myof","Cd1d1","S1pr3","Pik3r4",
                            "Zeb2","Adgre1","Cdkn2a","Itgax","Nt5e","Itgam","Tbx21",
                            "Gcsam","S1pr2","Mef2b","Jchain","Mki67","Mybl1","Aicda","Slc41a2","Bcl6","Zbtb20","Ighg2b","Ighg2c","Igha"
))

# input major marker genes
T_markers_major <- intersect(T_markers_major,row.names(T_cell.combined[["RNA"]]$counts))
B_markers_major <- intersect(B_markers_major,row.names(B_cell.combined[["RNA"]]$counts))

DefaultAssay(PICs.combined) <- "RNA"
PICs.combined <- NormalizeData(PICs.combined)
PICs.combined <- ScaleData(PICs.combined)


PIC_markerExp_cmb <- list()
PICs_label_cmb <- c()
for (t in major_T){
  for (b in major_B){
    T_name <- names(translate_T[translate_T==t])
    B_name <- names(translate_B[translate_B==b])
    PICs_select <- paste0(T_name," and ",B_name)
    print(PICs_select)
    PICs_select_IDs <- names(PICs.combined$predicted_PICs[PICs.combined$predicted_PICs == PICs_select])
    PICs_label <- rep(PICs_select,length(PICs_select_IDs))
    PICs_label_cmb <- c(PICs_label_cmb,PICs_label)
    
    PIC_markerExp_T <- PICs.combined@assays[["RNA"]]$data[T_markers_major,PICs_select_IDs]
    PIC_markerExp_B <- PICs.combined@assays[["RNA"]]$data[B_markers_major,PICs_select_IDs]
    PIC_markerExp <- rbind(PIC_markerExp_T,PIC_markerExp_B)
    PIC_markerExp_cmb[[PICs_select]] <- PIC_markerExp
  }
}
PIC_markerExp_cmb[["CD4+ IFITM+ CISH+ T Cells and ABCs"]] <- NULL
PIC_markerExp_cmb[["CD4+ IFITM+ CISH+ T Cells and GC B Cells"]] <- NULL
# PIC_markerExp_cmb[["Tregs and FOBs"]] <- NULL
# PIC_markerExp_cmb[["Tregs and MZBs"]] <- NULL
# PIC_markerExp_cmb[["Tregs and ABCs"]] <- NULL
# PIC_markerExp_cmb[["Tregs and GC B Cells"]] <- NULL

PICs_label_cmb <- PICs_label_cmb[!PICs_label_cmb %in% c("CD4+ IFITM+ CISH+ T Cells and ABCs","CD4+ IFITM+ CISH+ T Cells and GC B Cells")]
library(tidyverse)
PIC_markerExp_cmb.mx <- reduce(PIC_markerExp_cmb,cbind)

PIC_markerExp_cmb.mx <- as.matrix(PIC_markerExp_cmb.mx)

library(circlize)
library(viridis)
library(RColorBrewer)
# Define specific breakpoints
#breaks <- c(0, 0.5, 1, 1.5)
# Get the colors from the magma palette
colors <- magma(10,direction = -1)
# Create color function with colorRamp2
col_fun <- colorRamp2(breaks, colors)

# Define the range of your data for color mapping
data_range <- range(PIC_markerExp_cmb.mx, na.rm = TRUE)
# Create a color function using colorRamp2
color_function <- colorRamp2(c(data_range[1], data_range[2]), rev(magma(10)))


# annotate condition
age_anno <- PICs.combined$orig.ident[colnames(PIC_markerExp_cmb.mx)]
age_color <- brewer.pal(9,"YlGnBu")
age_color <- age_color[c(5,9)]
names(age_color) <- unique(age_anno)

ha.top = HeatmapAnnotation(
  Age = age_anno,
  col = list(Age = age_color)
)

# annotate cell types
T_cell_label <- gsub(" and .*","",PICs_label_cmb)
B_cell_label <- gsub(".* and ","",PICs_label_cmb)

# T cell colors
T_cells_colors <- brewer.pal(9, "YlOrRd")
T_cells_colors <- T_cells_colors[3:9]
names(T_cells_colors) <- unique(T_cell_label)
# B cell colors
B_cells_colors <- brewer.pal(9, "YlGn")
B_cells_colors <- B_cells_colors[c(4,6,8,9)]
names(B_cells_colors) <- unique(B_cell_label)

ha.bottom = HeatmapAnnotation(
  B_cells = B_cell_label,
  T_cells = T_cell_label,
  col = list(B_cells = B_cells_colors,
             T_cells = T_cells_colors
  )
)


pdf(file="heatmap_PICs_GeneExp_v1.pdf", width=14, height = 12)
ht <- Heatmap(PIC_markerExp_cmb.mx,
              col = rev(magma(10)),
              name = "log(normalized UMIs)",
              top_annotation = ha.top,
              bottom_annotation = ha.bottom,
              row_title_side = "left",
              row_names_side = "left",
              show_column_names = F,
              cluster_rows=F,
              cluster_columns=F,
              row_split = c(rep("T cells",length(T_markers_major)),rep("B cells",length(B_markers_major))),
              row_gap=unit(3, "mm"),
              column_split = PICs_label_cmb,
              column_gap = unit(0, "mm"), border = TRUE
)
ht<-draw(ht)
dev.off()

# UMIs: Create a scatter plot with smaller dots major cells only
#PICs_UMI_major <- PICs_UMI[PICs.combined$predicted_PICs %in% PIC4compare]
PICs_UMI_major <- PICs_UMI[colnames(PIC_markerExp_cmb.mx)]
pdf("PICs_UMI_major_v1.pdf", width = 10, height = 3)

plot(x = seq_along(PICs_UMI_major), y = log2(PICs_UMI_major), 
     main = paste0(length(PICs_UMI_major)," PICs"), 
     xlab = "", ylab = "UMIs(log)", 
     pch = 19,  # solid circles
     cex = 0.1,
     xaxt = "n",
     xaxs = "i")  # Do not draw the x-axis) # smaller dots
dev.off()


### Cell component changes ###
translated_major_T <- names(translate_T[translate_T %in% major_T])
prop_major_T <- T_cell.combined$cell_type[T_cell.combined$cell_type %in% translated_major_T]
table(prop_major_T)

prop_PICs4T <- PICs.combined$predicted_PICs[PICs.combined$predicted_PICs != "unknown"]
prop_PICs4T <- gsub(" and.*","",prop_PICs4T)
table(prop_PICs4T)

prop_PICs4B <- PICs.combined$predicted_PICs[PICs.combined$predicted_PICs != "unknown"]
prop_PICs4B <- gsub(".*and ","",prop_PICs4B)
table(prop_PICs4B)

# proportation of young, aged cells
prop_age <- function(prop_major,SeuratObj,Aged_ident){
  Major_ID <- names(prop_major)
  Y_ID <- names(SeuratObj$orig.ident[SeuratObj$orig.ident==Aged_ident])
  Y_Major_ID <- intersect(Major_ID,Y_ID)
  prop_Y_major <- prop_major[Y_Major_ID]
  return(prop_Y_major)
}

prop_Y_major_T <- prop_age(prop_major=prop_major_T, SeuratObj=T_cell.combined, Aged_ident="Young_T_Cells")
prop_A_major_T <- prop_age(prop_major=prop_major_T, SeuratObj=T_cell.combined, Aged_ident="Aged_T_Cells")

prop_Y_major_B <- prop_age(prop_major=prop_major_B, SeuratObj=B_cell.combined, Aged_ident="Young_B_Cells")
prop_A_major_B <- prop_age(prop_major=prop_major_B, SeuratObj=B_cell.combined, Aged_ident="Aged_B_Cells")


prop_Y_major_PIC4T <- prop_age(prop_major=prop_PICs4T, SeuratObj=PICs.combined, Aged_ident="Young_PICS")
prop_A_major_PIC4T <- prop_age(prop_major=prop_PICs4T, SeuratObj=PICs.combined, Aged_ident="Aged_PICS")

prop_Y_major_PIC4B <- prop_age(prop_major=prop_PICs4B, SeuratObj=PICs.combined, Aged_ident="Young_PICS")
prop_A_major_PIC4B <- prop_age(prop_major=prop_PICs4B, SeuratObj=PICs.combined, Aged_ident="Aged_PICS")


# T cell in PIC general
sub.propT_PIC_G <- cell_prop_bar_df(prop_major_T,prop_PICs4T,"T cells", "General")
# T cell in PIC young
sub.propT_PIC_Y <- cell_prop_bar_df(prop_Y_major_T,prop_Y_major_PIC4T,"T cells", "Young")
# T cell in PIC aged
sub.propT_PIC_A <- cell_prop_bar_df(prop_A_major_T,prop_A_major_PIC4T,"T cells", "Aged")

library(patchwork)
# Combine the plots and use `&` to apply a single guide to the combined plot
combined_plot <- (sub.propT_PIC_G / sub.propT_PIC_Y / sub.propT_PIC_A) & theme(legend.position = "bottom")
# Specify that the layout should keep the guide
combined_plot <- combined_plot + plot_layout(guides = 'collect', ncol = 1)
# Plot
combined_plot
ggsave("T_PICs_prop_v1.pdf",plot=combined_plot, width=9, height=10)


# B cell in PIC general
sub.propB_PIC_G <- cell_prop_bar_df(prop_major_B,prop_PICs4B,"B cells", "General")
# B cell in PIC young
sub.propB_PIC_Y <- cell_prop_bar_df(prop_Y_major_B,prop_Y_major_PIC4B,"B cells", "Young")
# B cell in PIC aged
sub.propB_PIC_A <- cell_prop_bar_df(prop_A_major_B,prop_A_major_PIC4B,"B cells", "Aged")

library(patchwork)
# Combine the plots and use `&` to apply a single guide to the combined plot
combined_plot <- (sub.propB_PIC_G / sub.propB_PIC_Y / sub.propB_PIC_A) & theme(legend.position = "bottom")
# Specify that the layout should keep the guide
combined_plot <- combined_plot + plot_layout(guides = 'collect', ncol = 1)
# Plot
combined_plot
ggsave("B_PICs_prop_v1.pdf",plot=combined_plot, width=9, height=10)


# calculate p-values of cell prop changes
T_P_value_prop_overall <- P_value_prop(prop_major_T,prop_PICs4T)
T_P_value_prop_Young <- P_value_prop(prop_Y_major_T,prop_Y_major_PIC4T)
T_P_value_prop_Aged <- P_value_prop(prop_A_major_T,prop_A_major_PIC4T)

B_P_value_prop_overall <- P_value_prop(prop_major_B,prop_PICs4B)
B_P_value_prop_Young <- P_value_prop(prop_Y_major_B,prop_Y_major_PIC4B)
B_P_value_prop_Aged <- P_value_prop(prop_A_major_B,prop_A_major_PIC4B)

# create a data frame to store the result
# T cells
T_P_value_prop_overall.df <- as.data.frame(T_P_value_prop_overall); T_P_value_prop_overall.df$Cell_Types <- names(T_P_value_prop_overall)
T_P_value_prop_Young.df <- as.data.frame(T_P_value_prop_Young); T_P_value_prop_Young.df$Cell_Types <- names(T_P_value_prop_Young)
T_P_value_prop_Aged.df <- as.data.frame(T_P_value_prop_Aged); T_P_value_prop_Aged.df$Cell_Types <- names(T_P_value_prop_Aged)
T_prop_P_values <- purrr::reduce(list(T_P_value_prop_overall.df,T_P_value_prop_Young.df,T_P_value_prop_Aged.df),
                                 ~dplyr::full_join(.x, .y, by = "Cell_Types"))
T_prop_P_values <- tibble::column_to_rownames(T_prop_P_values,"Cell_Types")
write.csv(T_prop_P_values,"T_prop_P_values.csv")
# B cells
B_P_value_prop_overall.df <- as.data.frame(B_P_value_prop_overall); B_P_value_prop_overall.df$Cell_Types <- names(B_P_value_prop_overall)
B_P_value_prop_Young.df <- as.data.frame(B_P_value_prop_Young); B_P_value_prop_Young.df$Cell_Types <- names(B_P_value_prop_Young)
B_P_value_prop_Aged.df <- as.data.frame(B_P_value_prop_Aged); B_P_value_prop_Aged.df$Cell_Types <- names(B_P_value_prop_Aged)
B_prop_P_values <- purrr::reduce(list(B_P_value_prop_overall.df,B_P_value_prop_Young.df,B_P_value_prop_Aged.df),
                                 ~dplyr::full_join(.x, .y, by = "Cell_Types"))
B_prop_P_values <- tibble::column_to_rownames(B_prop_P_values,"Cell_Types")
write.csv(B_prop_P_values,"B_prop_P_values.csv")



### generate expect gene expression pattern from artificial doublets (AD) - all combinations ###
gene2exam <- top_gene_per_cluster
t_data.exp <- T_cell.combined[["RNA"]]$counts[gene2exam,]
b_data.exp <- B_cell.combined[["RNA"]]$counts[gene2exam,]

randSize <- 5000
combined_data.exp.ls <- list()
for (t in unique(T_cell.combined$cell_type)){
  print(t)
  t_ct <- names(T_cell.combined$cell_type[T_cell.combined$cell_type==t]) # select cell IDs of a T cell type
  for (b in unique(B_cell.combined$cell_type)){
    print(b)
    AD_label <- paste0(t," and ",b)
    print(paste0("generate artificial doublets (AD): ",AD_label))
    b_ct <- names(B_cell.combined$cell_type[B_cell.combined$cell_type==b]) # select cell IDs of a B cell type
    # random sampling T and B cells
    t_rand <- sample(t_ct,randSize,replace = T)
    b_rand <- sample(b_ct,randSize,replace = T)
    # get gene expression matrix
    t_data <- t_data.exp[,t_rand]
    b_data <- b_data.exp[,b_rand]
    # get total counts per cell for selected cells
    T_counts_per_cell <- T_UMI[t_rand]
    B_counts_per_cell <- B_UMI[b_rand]
    # average the gene expression of artificial doublets
    #T_prop <- t_data/(T_counts_per_cell + B_counts_per_cell)
    #B_prop <- b_data/(T_counts_per_cell + B_counts_per_cell)
    #T_prop <- t_data/(t_data + b_data)
    
    # proportion of each gene in the doublet
    T_prop <- rowMeans(t_data)/(rowMeans(t_data) + rowMeans(b_data))
    B_prop <- 1-T_prop
    
    # get expected logCPM for comparison with observed value
    data_AD_exp <- log1p((t_data + b_data)/(T_counts_per_cell + B_counts_per_cell)*1e6) # logCPM of AD
    colnames(data_AD_exp) <- paste0(AD_label,"_",seq(1,ncol(data_AD_exp)))
    # get mean expression of logCPM for each of gene
    expected_expression <- rowMeans(data_AD_exp)
    
    # create a data frame to store the expected proportation and expression
    AD_exp_df <- data.frame(expected_expression=expected_expression,T_prop=T_prop, B_prop=B_prop)
    
    # add the doublet label to the dataframe
    AD_exp_df$AD_label <- AD_label
    combined_data.exp.ls[[AD_label]] <- AD_exp_df
  }
}
#combined_data_exp <- purrr::reduce(combined_data.exp.ls, bind_rows)
#combined_data_exp$AD_label <- as.factor(combined_data_exp$AD_label) # convert the label to factor


### generate expect gene expression pattern from artificial doublets (AD) - all combinations & all genes ###

# function for gene expression levels
obs_vs_exp_gene <- function(condition, Cell1data, Cell2data, PICdata,randSize){
  # random sampling T and B cells
  t_rand <- sample(Cell1data,randSize,replace = T)
  b_rand <- sample(Cell2data,randSize,replace = T)
  # get gene expression matrix
  t_data <- t_data.exp.all[,t_rand]
  b_data <- b_data.exp.all[,b_rand]
  # get total counts per cell for selected cells
  T_counts_per_cell <- T_UMI[t_rand]
  B_counts_per_cell <- B_UMI[b_rand]

  # proportion of each gene in the doublet
  T_prop <- rowMeans(t_data)/(rowMeans(t_data) + rowMeans(b_data))
  B_prop <- 1-T_prop
  
  # get expected logCPM for comparison with observed value
  data_AD_exp <- log1p((t_data + b_data)/(T_counts_per_cell + B_counts_per_cell)*1e6) # logCPM of AD
  colnames(data_AD_exp) <- paste0(AD_label,"_",seq(1,ncol(data_AD_exp)))
  # get mean expression of logCPM for each of gene
  expected_expression <- rowMeans(data_AD_exp)
  
  if (length(PICdata) > 0){
    # get observed PICs gene expression
    PICs_data <- PICs.combined[["RNA"]]$counts[,PICdata]
    PICs_counts_per_cell <- PICs_UMI[PICdata]
    data_PICs_obs <- log1p(PICs_data/PICs_counts_per_cell*1e6)
    observed_expression <- rowMeans(data_PICs_obs)
    #library(Matrix)
    # Step 1: Compute the standard deviation for each gene, considering only non-zero values
    # Convert the sparse matrix to a dense format temporarily for the non-zero values
    PICdense_data <- as.matrix(data_PICs_obs)
    # Replace zeros with NA to exclude them from the computation
    PICdense_data[PICdense_data == 0] <- NA
    gene_sds <- apply(PICdense_data, 1, sd, na.rm = TRUE)
    # Step 2: Count the number of non-zero observations for each gene
    non_zero_counts <- rowSums(data_PICs_obs != 0)
    # Step 3: Calculate the standard error for each gene
    gene_se <- gene_sds / sqrt(non_zero_counts)
    # gene_se now contains the standard error for each gene, considering only non-zero values
    ymin <- observed_expression - gene_se
    ymax <- observed_expression + gene_se
    
    # select overlapped genes
    genes2compare <- intersect(names(expected_expression),names(observed_expression))
    expected_expression <- expected_expression[genes2compare]
    observed_expression <- observed_expression[genes2compare]
    ymin <- ymin[genes2compare]
    ymax <- ymax[genes2compare]
    T_prop <- T_prop[genes2compare]
    B_prop <- B_prop[genes2compare]
    
    # create a data frame to store the expected proporation and expression
    AD_exp_df <- data.frame(gene=rep(names(expected_expression),each=2),expected_expression=rep(expected_expression,each=2),prop=c(rbind(T_prop, B_prop)), 
                            Expression_level=rep(expected_expression,each=2)*c(rbind(T_prop, B_prop)), 
                            Cell_Type=rep(c(paste0("Expected ",t),paste0("Expected ",b)),length(expected_expression)), 
                            Cell_Type2=rep(c("Expected T cells","Expected B cells"),length(expected_expression)),
                            PICs= paste0("Expected ",AD_label),AD_label=AD_label,
                            ymin=NA, ymax=NA,condition=condition)
    # add observed information
    Obs_df <- data.frame(gene=names(observed_expression),expected_expression=NA,prop=NA,
                         Expression_level=observed_expression,Cell_Type=AD_label,
                         PICs= paste0("Observed ",AD_label),AD_label=AD_label,
                         Cell_Type2="Observed",
                         ymin=ymin, ymax=ymax,condition=condition)
    Obs_Exp_df <- rbind(AD_exp_df,Obs_df)
    return(Obs_Exp_df)
  }
}


# normalized expression, proportionate, labels, conditions
t_data.exp.all <- T_cell.combined[["RNA"]]$counts
b_data.exp.all <- B_cell.combined[["RNA"]]$counts

randSize <- 5000
combined_data.exp.all.ls <- list()
for (t in unique(T_cell.combined$cell_type)){
  print(t)
  # select cell IDs of a T cell type
  t_ct_Young <- T_cell.combined@meta.data %>% filter(cell_type==t & orig.ident=="Young_T_Cells") %>% row.names() # young
  t_ct_Aged <- T_cell.combined@meta.data %>% filter(cell_type==t & orig.ident=="Aged_T_Cells") %>% row.names() # aged
  t_ct_overall <- names(T_cell.combined$cell_type[T_cell.combined$cell_type==t]) # overall
  for (b in unique(B_cell.combined$cell_type)){
    print(b)
    # select cell IDs of a B cell type
    b_ct_Young <- B_cell.combined@meta.data %>% filter(cell_type==b & orig.ident=="Young_B_Cells") %>% row.names() # young
    b_ct_Aged <- B_cell.combined@meta.data %>% filter(cell_type==b & orig.ident=="Aged_B_Cells") %>% row.names() # aged
    b_ct_overall <- names(B_cell.combined$cell_type[B_cell.combined$cell_type==b]) # overall
    # create name of doublet
    AD_label <- paste0(t," and ",b)
    print(paste0("generate artificial doublets (AD): ",AD_label))
    
    # select cell IDs of observed PICs
    PICs_ct_Young <- PICs.combined@meta.data %>% filter(predicted_PICs==AD_label & orig.ident=="Young_PICS") %>% row.names() # young
    PICs_ct_Aged <- PICs.combined@meta.data %>% filter(predicted_PICs==AD_label & orig.ident=="Aged_PICS") %>% row.names() # aged
    PICs_ct_overall <- names(PICs.combined$predicted_PICs[PICs.combined$predicted_PICs==AD_label]) # overall
    
      # Young
      if ((length(t_ct_Young) >= 3 | length(b_ct_Young) >= 3) & length(PICs_ct_Young) >= 3){
        print(paste0("generate artificial doublets (AD): ",AD_label," Young"))
        obs_vs_exp_gene_young <- obs_vs_exp_gene(condition="Young",Cell1data=t_ct_Young, Cell2data=b_ct_Young, PICdata=PICs_ct_Young,randSize=randSize)
      }
      # Aged
      if ((length(t_ct_Aged) >= 3 | length(b_ct_Aged) >= 3) & length(PICs_ct_Aged) >= 3){
        print(paste0("generate artificial doublets (AD): ",AD_label," Aged"))
        obs_vs_exp_gene_aged <- obs_vs_exp_gene(condition="Aged",Cell1data=t_ct_Aged, Cell2data=b_ct_Aged, PICdata=PICs_ct_Aged,randSize=randSize)
      }
      # overall
      if ((length(t_ct_overall) >= 3 | length(b_ct_overall) >= 3) & length(PICs_ct_overall) >= 3){
        print(paste0("generate artificial doublets (AD): ",AD_label," overall"))
        obs_vs_exp_gene_overall <- obs_vs_exp_gene(condition="Overall",Cell1data=t_ct_overall, Cell2data=b_ct_overall, PICdata=PICs_ct_overall,randSize=randSize)
      }
      # add to a list
      combined_data.exp.all.ls[[AD_label]] <- rbind(obs_vs_exp_gene_young,obs_vs_exp_gene_aged,obs_vs_exp_gene_overall)
  }
}

#combined_data.exp.all.df <- purrr::reduce(combined_data.exp.ls, bind_rows)
major_TB <- c()
for (t in major_T){
  t <- names(translate_T[translate_T==t])
  for (b in major_B){
    b <- names(translate_B[translate_B==b])
    AD_label <- paste0(t," and ",b)
    print(AD_label)
    major_TB <- c(major_TB,AD_label)
  }
}

# combined_data.exp.major.ls <- list()
# for (m in major_TB){
#   combined_data.exp.major.ls[[m]] <- combined_data.exp.all.ls[[m]]
# }

combined_data.exp.major.ls <- combined_data.exp.all.ls[names(combined_data.exp.all.ls) %in% major_TB]

combined_data.exp.major.df <- purrr::reduce(combined_data.exp.major.ls, bind_rows)


# Plotting with error bars added only for "PIC"
test.df <- combined_data.exp.major.df %>% filter(condition != "Overall") %>%
  distinct()

# ct_order <- c()
# for (n in names(combined_data.exp.all.ls)){
#   ct_order <- c(ct_order,c(paste0("Expected ",n),paste0("Observed ",n)))
# }
test.df <- test.df %>%
  mutate(
    condition = factor(condition, levels = c("Young", "Aged")),
    PICs = factor(PICs, levels = unique(test.df$PICs)),
    group = interaction(PICs, condition) # Create interaction for grouping
  )
#test.df$test <- paste0(test.df$PICs,"_",test.df$condition)
# Get the default color set
default_colors <- scales::hue_pal()(3)
# Change the third color to grey
default_colors[3] <- "darkgrey"

gene2look <- c("Ighv1−81","Ighv4−1","Pilra","Trbv3","C1qb","Igkv1−99")
gene2look <- c("Tbx21","Fcer2a","Il21","Zeb2","Ifng","Tnf","Hacd4","Stat5a")
gene2look <- c("Ccl3","Tnf","Hacd4","Stat5a")
#gene2look <- names(top10_absSum[1:6])


unique(test.df$gene)[unique(test.df$gene) %in% gene2look] # quick check if gene in the df

combined_plot.ObsvsExp.bar <- list()
for (g in gene2look){
  print(g)
  test.df.g <- test.df %>% filter(gene==g)
  # Plot
  if (nrow(test.df.g)>0){
    line_positions <- length(test.df.g$condition[test.df.g$condition=="Young"])/3*2 + 0.5  # Between each group
    gene_barplot <- ggplot(test.df.g, aes(fill = Cell_Type2, y = Expression_level, x = group)) + 
      geom_bar(stat = "identity",color = "black") + # Use 'position_dodge' to place bars side by side
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.4, position = position_dodge(.9)) +
      scale_x_discrete(labels = function(x) gsub("^Expected ||^Observed ", "", x)) + # Clean up x-axis labels (removing the 'PICs' part)
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            plot.margin = unit(c(1, 1, 1, 6), "lines")) + 
      scale_fill_manual(values = default_colors) +
      labs(fill = "") +
      xlab("PICs") +
      ylab(g) +
      geom_vline(xintercept = line_positions, linetype = "dashed", color = "grey",size=1)  # Add dashed lines
    ggsave(plot=gene_barplot,paste0(g,"_PICs_bar.pdf"), width = 16, height = 8)
    combined_plot.ObsvsExp.bar[[g]] <- gene_barplot
  } else {print(paste0("Abscent gene: ",g))}
}

library(patchwork)
combined_plot.ObsvsExp.bar.comb <- wrap_plots(combined_plot.ObsvsExp.bar, ncol = 2)
print(combined_plot.ObsvsExp.bar.comb)
ggsave("genes_combined_plot_ObsvsExp.pdf",width=30,height=35)






combined_data.exp.all.ls0 <- combined_data.exp.all.ls # backup

### generate DEGs of observed vs expected (AD) - major cell types ###

# function for Wilcoxon Rank Sum Test for DEGs
computeStats <- function(geneIndex) {
  # Extract the gene name
  geneName <- rownames(data_PICs_obs)[geneIndex]
  # Extract expression values for the current gene and ensure they are numeric
  expression_T <- as.numeric(data_PICs_obs[geneIndex, ])
  expression_B <- as.numeric(data_AD_exp[geneIndex, ])
  # Perform Wilcoxon Rank Sum test
  test_result <- tryCatch({
    wilcox.test(expression_T, expression_B, exact = FALSE)
  }, error = function(e) {
    # Return NA values in case of an error (e.g., all values are identical)
    list(p.value = NA)
  })
  p_value <- test_result$p.value
  # Calculate log2FC
  avg_expr_T <- mean(expression_T)
  avg_expr_B <- mean(expression_B)
  fold_change <- (avg_expr_T + 1e-6) / (avg_expr_B + 1e-6)  # Avoid division by zero
  log2FC <- log2(fold_change)
  # Return a named vector or list including the gene name
  return(list(gene = geneName, log2FC = log2FC, p_value = p_value))
}

# calculate the total count per cell of PICs
PICs_UMI <- colSums(PICs.combined[["RNA"]]$counts)

# preprocess the sparse matrix: convert to dense dataframe
T_RNA_df <- as.data.frame(T_cell.combined[["RNA"]]$counts)
B_RNA_df <- as.data.frame(B_cell.combined[["RNA"]]$counts)
PICs_RNA_df <- as.data.frame(PICs.combined[["RNA"]]$counts)

rownames_TB <- intersect(rownames(T_RNA_df),rownames(B_RNA_df))
rownames_all <- intersect(rownames_TB,rownames(PICs_RNA_df))
T_RNA_df <- T_RNA_df[rownames_all,]
B_RNA_df <- B_RNA_df[rownames_all,]
PICs_RNA_df <- PICs_RNA_df[rownames_all,]

library(parallel)
# number of cores for parallel computing
no_cores <- 20
# Initiate cluster
cl <- makeCluster(no_cores)

randSize <- 5000
ObsvsExp.DEGs.ls <- list()
for (t in major_T){
  t <- names(translate_T[translate_T==t])
  print(t)
  t_ct <- names(T_cell.combined$cell_type[T_cell.combined$cell_type==t]) # select cell IDs of a T cell type
  for (b in major_B){
    b <- names(translate_B[translate_B==b])
    print(b)
    AD_label <- paste0(t," and ",b)
    print(paste0("generate artificial doublets (AD): ",AD_label))
    b_ct <- names(B_cell.combined$cell_type[B_cell.combined$cell_type==b]) # select cell IDs of a B cell type
    # random sampling T and B cells
    t_rand <- sample(t_ct,randSize,replace = T)
    b_rand <- sample(b_ct,randSize,replace = T)
    # get gene expression matrix
    t_data <- T_RNA_df[,t_rand]
    b_data <- B_RNA_df[,b_rand]
    # get total counts per cell for selected cells
    T_counts_per_cell <- T_UMI[t_rand]
    B_counts_per_cell <- B_UMI[b_rand]
    
    # get expected logCPM
    data_AD_exp <- log1p((t_data + b_data)/(T_counts_per_cell + B_counts_per_cell)*1e6) # logCPM of AD
    colnames(data_AD_exp) <- paste0(AD_label,"_",seq(1,ncol(data_AD_exp)))
    
    # Calculate the number of non-zero entries for each gene
    non_zero_counts <- rowSums(data_AD_exp > 0)
    # Create a logical vector to identify genes with more than 3 non-zero entries; and min.pct = 0.01
    genes_with_non_zero <- non_zero_counts > 3 & non_zero_counts > ncol(data_AD_exp)*0.01
    # Subset the original data frame to include only those genes
    data_AD_exp <- data_AD_exp[genes_with_non_zero, ]
    
    # get observed logCPM
    PICs_ct <- names(PICs.combined$predicted_PICs[PICs.combined$predicted_PICs==AD_label])
    if (length(PICs_ct)>=3){ # only consider cell number >= 3
      PICs_data <- PICs_RNA_df[,PICs_ct]
      PICs_counts_per_cell <- PICs_UMI[PICs_ct]
      data_PICs_obs <- log1p(PICs_data/PICs_counts_per_cell*1e6)
      
      # Calculate the number of non-zero entries for each gene
      non_zero_counts <- rowSums(data_PICs_obs > 0)
      # Create a logical vector to identify genes with more than 3 non-zero entries; and min.pct = 0.01
      genes_with_non_zero <- non_zero_counts > 3 & non_zero_counts > ncol(data_PICs_obs)*0.01
      # Subset the original data frame to include only those genes
      data_PICs_obs <- data_PICs_obs[genes_with_non_zero, ]
      
      gene_filtered <- intersect(rownames(data_AD_exp),rownames(data_PICs_obs))
      data_AD_exp <- data_AD_exp[gene_filtered,]
      data_PICs_obs <- data_PICs_obs[gene_filtered,]
      # Perform Wilcoxon Rank Sum Test for Each Gene
      # Calculate Log2 Fold Change for Each Gene
      clusterExport(cl, varlist = c("data_PICs_obs", "data_AD_exp"))
      
      # Execute in parallel
      results <- parLapply(cl, 1:nrow(data_PICs_obs), computeStats)
      
      # Convert results list to a data frame
      results_df <- do.call(rbind, results)
      results_df <- data.frame(gene=as.character(results_df[,1]),
                               log2FC=as.numeric(results_df[,2]),
                               p_value=as.numeric(results_df[,3]))
      # Adjust p-values
      results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")
      
      # add the DEG result to a list
      ObsvsExp.DEGs.ls[[AD_label]] <- results_df
    }
  }
}
# stop cluster
stopCluster(cl)

ObsvsExp.DEGs.ls0 <- ObsvsExp.DEGs.ls # backup

# adjust the column name of each element in the list
for (o in names(ObsvsExp.DEGs.ls)){
  ObsvsExp.DEGs <- ObsvsExp.DEGs.ls[[o]]
  colnames(ObsvsExp.DEGs) <- c("gene",paste0(names(ObsvsExp.DEGs)[-1],"_",o)) # change column names
  ObsvsExp.DEGs.ls[[o]] <- ObsvsExp.DEGs
}

ObsvsExp.DEGs.df <- purrr::reduce(ObsvsExp.DEGs.ls, full_join, by = "gene")


# grep() returns the indices of column names that match the pattern
# The pattern "log2FC" is looked for in the names of the dataframe columns
log2FC_columns <- grep("log2FC", names(ObsvsExp.DEGs.df), value = TRUE)
# Selecting columns that contain "log2FC" in their names
ObsvsExp.DEGs.df.log <- ObsvsExp.DEGs.df[, c("gene",log2FC_columns)]

# Calculate the standard deviation for each gene across all samples, ignoring the first column
std_dev_log2FC_per_gene <- apply(ObsvsExp.DEGs.df.log[,-1], 1, sd, na.rm = TRUE)
names(std_dev_log2FC_per_gene) <- ObsvsExp.DEGs.df.log[,1]
std_dev_log2FC_per_gene <- sort(std_dev_log2FC_per_gene,decreasing=T)
top10_dev <- std_dev_log2FC_per_gene[1:10] # select 10 most deviant genes

absSum_per_gene <- rowSums(abs(ObsvsExp.DEGs.df.log[-1]),na.rm = T)
names(absSum_per_gene) <- ObsvsExp.DEGs.df.log[,1]
absSum_per_gene <- sort(absSum_per_gene,decreasing=T)
top10_absSum <- absSum_per_gene[1:10] # select 10 most log2FC genes




# volcano plot for each PICs vs AD
dir.create("ObsvsExp_Volcano")
major_cell4comp <- names(ObsvsExp.DEGs.ls)
ObsvsExp.vol.plot.ls <- list()
for (p in major_cell4comp){
  ObsvsExp_one_PIC <- ObsvsExp.DEGs.ls[[p]]
  # #markers4plot <- read.csv(paste0("DEGs_PICs_age/",d))
  # # sort genes
  # ObsvsExp_one_PIC$rank <- abs(ObsvsExp_one_PIC$log2FC_T0_0_B0*(-log10(ObsvsExp_one_PIC$adjusted_p_value_T0_0_B0)))
  # ObsvsExp_one_PIC <- ObsvsExp_one_PIC %>% arrange(desc(rank))
  # topgenes <- head(ObsvsExp_one_PIC$gene,n=20)
  topgenes_u <- ObsvsExp_one_PIC %>% filter(adjusted_p_value <0.05 & log2FC >0) %>% arrange(desc(log2FC)) %>% head(n=20)
  topgenes_u <- topgenes_u$gene
  topgenes_d <- ObsvsExp_one_PIC %>% filter(adjusted_p_value <0.05 & log2FC <0) %>% arrange(log2FC) %>% head(n=20)
  topgenes_d <- topgenes_d$gene
  topgenes <- c(topgenes_u,topgenes_d)
  # Create a label column, label only selected genes
  #ObsvsExp_one_PIC$label <- ifelse(ObsvsExp_one_PIC$gene %in% topgenes, ObsvsExp_one_PIC$gene, "")
  vol.plot <- EnhancedVolcano(ObsvsExp_one_PIC,
                              lab = ObsvsExp_one_PIC$gene,
                              selectLab=topgenes,
                              x = colnames(ObsvsExp_one_PIC)[2],
                              y = colnames(ObsvsExp_one_PIC)[4],
                              #xlim =c(-1.5,2),
                              ylab = bquote(~-Log[10] ~ italic(adj.P)),
                              legendLabels = c("NS", expression(Log[2] ~ FC), "adj.p-value", expression(adj.p - value ~ and
                                                                                                        ~ log[2] ~ FC)),
                              FCcutoff = 0.25,
                              pCutoff = 0.05,
                              labSize = 5, 
                              drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                              #max.overlaps = 20,
                              title = p,
                              subtitle = bquote(italic("Observed vs Expected")))
  ggsave(plot=vol.plot,paste0("ObsvsExp_Volcano/",p,"_ObsvsExp_Volcano.pdf"), width=9,height=8)
  ObsvsExp.vol.plot.ls[[p]] <- vol.plot
}
library(patchwork)
combined_plot.ObsvsExp <- wrap_plots(ObsvsExp.vol.plot.ls, ncol = 5)
print(combined_plot.ObsvsExp)
ggsave("ObsvsExp_Volcano/Volcano_PICs_ObsvsExp.pdf",width=30,height=35)

ObsvsExp.vol.plot.ls1 <- ObsvsExp.vol.plot.ls[!names(ObsvsExp.vol.plot.ls) %in% 
                                                c("TFH Cells and GC B Cells","Treg and GC B Cells","CD8+ EM T Cells and GC B Cells",
                                                  "CX4CR1_hi CD4+ EM T Cells and GC B Cells","CD8+ Naive T Cells and ABCs","CD8+ Naive T Cells and GC B Cells")]
combined_plot.ObsvsExp1 <- wrap_plots(ObsvsExp.vol.plot.ls1, ncol = 5)
ggsave(plot=combined_plot.ObsvsExp1,"ObsvsExp_Volcano/Volcano_PICs_ObsvsExp1.pdf",width=35,height=30)

# save the DEG list
dir.create("ObsvsExp_DEGs")
for (d in names(ObsvsExp.DEGs.ls)){
  write.csv(ObsvsExp.DEGs.ls[[d]],paste0("ObsvsExp_DEGs/",d,"_ObservedvsExpected.csv"))
}


### gene expression for young and aged separately ###
# young version of combined_data.exp.ls
# number of cores for parallel computing
no_cores <- 20
# Initiate cluster
cl <- makeCluster(no_cores)

set.seed(1234)
randSize <- 5000
ObsvsExp.DEGs.ls_Young <- list()
ObsvsExp.DEGs.ls_Aged <- list()
for (t in major_T){
  t <- names(translate_T[translate_T==t])
  print(t)
  # select cell IDs of a T cell type
  t_ct_Young <- T_cell.combined@meta.data %>% filter(cell_type==t & orig.ident=="Young_T_Cells") %>% row.names() # young
  t_ct_Aged <- T_cell.combined@meta.data %>% filter(cell_type==t & orig.ident=="Aged_T_Cells") %>% row.names() # aged
  for (b in major_B){
    b <- names(translate_B[translate_B==b])
    print(b)
    # select cell IDs of a B cell type
    b_ct_Young <- B_cell.combined@meta.data %>% filter(cell_type==b & orig.ident=="Young_B_Cells") %>% row.names() # young
    b_ct_Aged <- B_cell.combined@meta.data %>% filter(cell_type==b & orig.ident=="Aged_B_Cells") %>% row.names() # aged
    # if non of t_ct_Young, t_ct_Aged, b_ct_Young, and b_ct_Aged is zero
    if (length(t_ct_Young)>0 & length(t_ct_Aged)>0 & length(b_ct_Young)>0 & length(b_ct_Aged)>0){
      AD_label <- paste0(t," and ",b)
      print(paste0("generate artificial doublets (AD): ",AD_label))
      # random sampling T and B cells
      # young
      t_rand_Y <- sample(t_ct_Young,randSize,replace = T)
      b_rand_Y <- sample(b_ct_Young,randSize,replace = T)
      # aged
      t_rand_A <- sample(t_ct_Aged,randSize,replace = T)
      b_rand_A <- sample(b_ct_Aged,randSize,replace = T)
      # get gene expression matrix
      # young
      t_data_Y <- T_RNA_df[,t_rand_Y]
      b_data_Y <- B_RNA_df[,b_rand_Y]
      # aged
      t_data_A <- T_RNA_df[,t_rand_A]
      b_data_A <- B_RNA_df[,b_rand_A]
      
      # get total counts per cell for selected cells
      # young
      T_counts_per_cell_Y <- T_UMI[t_rand_Y]
      B_counts_per_cell_Y <- B_UMI[b_rand_Y]
      # aged
      T_counts_per_cell_A <- T_UMI[t_rand_A]
      B_counts_per_cell_A <- B_UMI[b_rand_A]
      
      # get expected logCPM
      # young
      data_AD_exp_Y <- log1p((t_data_Y + b_data_Y)/(T_counts_per_cell_Y + B_counts_per_cell_Y)*1e6) # logCPM of AD
      colnames(data_AD_exp_Y) <- paste0(AD_label,"_",seq(1,ncol(data_AD_exp_Y)))
      # aged
      data_AD_exp_A <- log1p((t_data_A + b_data_A)/(T_counts_per_cell_A + B_counts_per_cell_A)*1e6) # logCPM of AD
      colnames(data_AD_exp_A) <- paste0(AD_label,"_",seq(1,ncol(data_AD_exp_A)))
      
      # Calculate the number of non-zero entries for each gene
      # young
      non_zero_counts_Y <- rowSums(data_AD_exp_Y > 0)
      # aged
      non_zero_counts_A <- rowSums(data_AD_exp_A > 0)
      
      # Create a logical vector to identify genes with more than 3 non-zero entries; and min.pct = 0.01
      # young
      genes_with_non_zero_Y <- non_zero_counts_Y > 3 & non_zero_counts_Y > ncol(data_AD_exp_Y)*0.01
      # aged
      genes_with_non_zero_A <- non_zero_counts_A > 3 & non_zero_counts_A > ncol(data_AD_exp_A)*0.01
      
      # Subset the original data frame to include only those genes
      # young
      data_AD_exp_Y <- data_AD_exp_Y[genes_with_non_zero_Y, ]
      # aged
      data_AD_exp_A <- data_AD_exp_A[genes_with_non_zero_A, ]
      
      # get observed logCPM
      # young
      PICs_ct_Young <- PICs.combined@meta.data %>% filter(predicted_PICs==AD_label & orig.ident=="Young_PICS") %>% row.names()
      # aged
      PICs_ct_Aged <- PICs.combined@meta.data %>% filter(predicted_PICs==AD_label & orig.ident=="Aged_PICS") %>% row.names()
      
      # Young: DEGs in Obs vs Exp
      if (length(PICs_ct_Young)>=3){ # only consider cell number >= 3
        PICs_data <- PICs_RNA_df[,PICs_ct_Young]
        PICs_counts_per_cell <- PICs_UMI[PICs_ct_Young]
        data_PICs_obs <- log1p(PICs_data/PICs_counts_per_cell*1e6)
        
        # Calculate the number of non-zero entries for each gene
        non_zero_counts <- rowSums(data_PICs_obs > 0)
        # Create a logical vector to identify genes with more than 3 non-zero entries; and min.pct = 0.01
        genes_with_non_zero <- non_zero_counts > 3 & non_zero_counts > ncol(data_PICs_obs)*0.01
        # Subset the original data frame to include only those genes
        data_PICs_obs <- data_PICs_obs[genes_with_non_zero, ]
        
        gene_filtered <- intersect(rownames(data_AD_exp_Y),rownames(data_PICs_obs))
        data_AD_exp <- data_AD_exp_Y[gene_filtered,] # assign data_AD_exp_Y to data_AD_exp for parallel computeStats function
        data_PICs_obs <- data_PICs_obs[gene_filtered,]
        # Perform Wilcoxon Rank Sum Test for Each Gene
        # Calculate Log2 Fold Change for Each Gene
        print("Calculate Log2 Fold Change for Each Gene in Young")
        clusterExport(cl, varlist = c("data_PICs_obs", "data_AD_exp"))
        
        # Execute in parallel
        print("Execute in parallel")
        results <- parLapply(cl, 1:nrow(data_PICs_obs), computeStats)
        
        # Convert results list to a data frame
        results_df <- do.call(rbind, results)
        results_df <- data.frame(gene=as.character(results_df[,1]),
                                 log2FC=as.numeric(results_df[,2]),
                                 p_value=as.numeric(results_df[,3]))
        # Adjust p-values
        results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")
        
        # add the DEG result to a list
        ObsvsExp.DEGs.ls_Young[[AD_label]] <- results_df
      }
      
      # Aged: DEGs in Obs vs Exp
      if (length(PICs_ct_Aged)>=3){ # only consider cell number >= 3
        PICs_data <- PICs_RNA_df[,PICs_ct_Aged]
        PICs_counts_per_cell <- PICs_UMI[PICs_ct_Aged]
        data_PICs_obs <- log1p(PICs_data/PICs_counts_per_cell*1e6)
        
        # Calculate the number of non-zero entries for each gene
        non_zero_counts <- rowSums(data_PICs_obs > 0)
        # Create a logical vector to identify genes with more than 3 non-zero entries; and min.pct = 0.01
        genes_with_non_zero <- non_zero_counts > 3 & non_zero_counts > ncol(data_PICs_obs)*0.01
        # Subset the original data frame to include only those genes
        data_PICs_obs <- data_PICs_obs[genes_with_non_zero, ]
        
        gene_filtered <- intersect(rownames(data_AD_exp_A),rownames(data_PICs_obs))
        data_AD_exp <- data_AD_exp_A[gene_filtered,] # assign data_AD_exp_A to data_AD_exp for parallel computeStats function
        data_PICs_obs <- data_PICs_obs[gene_filtered,]
        # Perform Wilcoxon Rank Sum Test for Each Gene
        # Calculate Log2 Fold Change for Each Gene
        print("Calculate Log2 Fold Change for Each Gene in Aged")
        clusterExport(cl, varlist = c("data_PICs_obs", "data_AD_exp"))
        
        # Execute in parallel
        print("Execute in parallel")
        results <- parLapply(cl, 1:nrow(data_PICs_obs), computeStats)
        
        # Convert results list to a data frame
        results_df <- do.call(rbind, results)
        results_df <- data.frame(gene=as.character(results_df[,1]),
                                 log2FC=as.numeric(results_df[,2]),
                                 p_value=as.numeric(results_df[,3]))
        # Adjust p-values
        results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")
        
        # add the DEG result to a list
        ObsvsExp.DEGs.ls_Aged[[AD_label]] <- results_df
      }
    }
  }
}
# stop cluster
stopCluster(cl)

# save the DEG list
# young
dir.create("ObsvsExp_DEGs_Young")
for (d in names(ObsvsExp.DEGs.ls_Young)){
  write.csv(ObsvsExp.DEGs.ls_Young[[d]],paste0("ObsvsExp_DEGs_Young/",d,"_ObservedvsExpected_Young.csv"))
}
# aged
dir.create("ObsvsExp_DEGs_Aged")
for (d in names(ObsvsExp.DEGs.ls_Aged)){
  write.csv(ObsvsExp.DEGs.ls_Aged[[d]],paste0("ObsvsExp_DEGs_Aged/",d,"_ObservedvsExpected_Aged.csv"))
}

### overall comparison ###
# number of cores for parallel computing
no_cores <- 20
# Initiate cluster
cl <- makeCluster(no_cores)

set.seed(1234)
randSize <- 5000
# select cell IDs of young T cells
t_ct_Young.overall <- T_cell.combined@meta.data %>% filter(orig.ident=="Young_T_Cells") %>% row.names() # young
# select cell IDs of aged T cells
t_ct_Aged.overall <- T_cell.combined@meta.data %>% filter(orig.ident=="Aged_T_Cells") %>% row.names() # aged
# select cell IDs of young B cells
b_ct_Young.overall <- B_cell.combined@meta.data %>% filter(orig.ident=="Young_B_Cells") %>% row.names() # young
# select cell IDs of aged B cells
b_ct_Aged.overall <- B_cell.combined@meta.data %>% filter(orig.ident=="Aged_B_Cells") %>% row.names() # aged

AD_label <- "Overall_T_and_B"
print(paste0("generate artificial doublets (AD): ",AD_label))
# random sampling T and B cells
# young
t_rand_Y.overall <- sample(t_ct_Young.overall,randSize,replace = T)
b_rand_Y.overall <- sample(b_ct_Young.overall,randSize,replace = T)
# aged
t_rand_A.overall <- sample(t_ct_Aged.overall,randSize,replace = T)
b_rand_A.overall <- sample(b_ct_Aged.overall,randSize,replace = T)
# get gene expression matrix
# young
t_data_Y.overall <- T_RNA_df[,t_rand_Y.overall]
b_data_Y.overall <- B_RNA_df[,b_rand_Y.overall]
# aged
t_data_A.overall <- T_RNA_df[,t_rand_A.overall]
b_data_A.overall <- B_RNA_df[,b_rand_A.overall]

# get total counts per cell for selected cells
# young
T_counts_per_cell_Y.overall <- T_UMI[t_rand_Y.overall]
B_counts_per_cell_Y.overall <- B_UMI[b_rand_Y.overall]
# aged
T_counts_per_cell_A.overall <- T_UMI[t_rand_A.overall]
B_counts_per_cell_A.overall <- B_UMI[b_rand_A.overall]

# get expected logCPM
# young
data_AD_exp_Y.overall <- log1p((t_data_Y.overall + b_data_Y.overall)/(T_counts_per_cell_Y.overall + B_counts_per_cell_Y.overall)*1e6) # logCPM of AD
colnames(data_AD_exp_Y.overall) <- paste0(AD_label,"_",seq(1,ncol(data_AD_exp_Y.overall)))
# aged
data_AD_exp_A.overall <- log1p((t_data_A.overall + b_data_A.overall)/(T_counts_per_cell_A.overall + B_counts_per_cell_A.overall)*1e6) # logCPM of AD
colnames(data_AD_exp_A.overall) <- paste0(AD_label,"_",seq(1,ncol(data_AD_exp_A.overall)))

# Calculate the number of non-zero entries for each gene
# young
non_zero_counts_Y.overall <- rowSums(data_AD_exp_Y.overall > 0)
# aged
non_zero_counts_A.overall <- rowSums(data_AD_exp_A.overall > 0)

# Create a logical vector to identify genes with more than 3 non-zero entries; and min.pct = 0.01
# young
genes_with_non_zero_Y.overall <- non_zero_counts_Y.overall > 3 & non_zero_counts_Y.overall > ncol(data_AD_exp_Y.overall)*0.01
# aged
genes_with_non_zero_A.overall <- non_zero_counts_A.overall > 3 & non_zero_counts_A.overall > ncol(data_AD_exp_A.overall)*0.01

# Subset the original data frame to include only those genes
# young
data_AD_exp_Y.overall <- data_AD_exp_Y.overall[genes_with_non_zero_Y.overall, ]
# aged
data_AD_exp_A.overall <- data_AD_exp_A.overall[genes_with_non_zero_A.overall, ]

# get observed logCPM
# young
PICs_ct_Young.overall <- PICs.combined@meta.data %>% filter(orig.ident=="Young_PICS") %>% row.names()
# aged
PICs_ct_Aged.overall <- PICs.combined@meta.data %>% filter(orig.ident=="Aged_PICS") %>% row.names()


# Young: overall DEGs in Obs vs Exp
PICs_data.overall <- PICs_RNA_df[,PICs_ct_Young.overall]
PICs_counts_per_cell.overall <- PICs_UMI[PICs_ct_Young.overall]
data_PICs_obs.overall <- log1p(PICs_data.overall/PICs_counts_per_cell.overall*1e6)

# Calculate the number of non-zero entries for each gene
non_zero_counts.overall <- rowSums(data_PICs_obs.overall > 0)
# Create a logical vector to identify genes with more than 3 non-zero entries; and min.pct = 0.01
genes_with_non_zero.overall <- non_zero_counts.overall > 3 & non_zero_counts.overall > ncol(data_PICs_obs.overall)*0.01
# Subset the original data frame to include only those genes
data_PICs_obs.overall <- data_PICs_obs.overall[genes_with_non_zero.overall, ]

gene_filtered.overall <- intersect(rownames(data_AD_exp_Y.overall),rownames(data_PICs_obs.overall))
data_AD_exp <- data_AD_exp_Y.overall[gene_filtered.overall,] # assign data_AD_exp_Y to data_AD_exp for parallel computeStats function
data_PICs_obs <- data_PICs_obs.overall[gene_filtered.overall,]
# Perform Wilcoxon Rank Sum Test for Each Gene
# Calculate Log2 Fold Change for Each Gene
print("Calculate Log2 Fold Change for Each Gene in Young")
clusterExport(cl, varlist = c("data_PICs_obs", "data_AD_exp"))

# Execute in parallel
print("Execute in parallel")
results <- parLapply(cl, 1:nrow(data_PICs_obs), computeStats)

# Convert results list to a data frame
results_df <- do.call(rbind, results)
results_df <- data.frame(gene=as.character(results_df[,1]),
                         log2FC=as.numeric(results_df[,2]),
                         p_value=as.numeric(results_df[,3]))
# Adjust p-values
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")

# add the DEG result to a list
ObsvsExp.DEGs.Young.overall.df <- results_df
write.csv(ObsvsExp.DEGs.Young.overall.df,"ObsvsExp.DEGs.Young.overall.csv")



# Aged: overall DEGs in Obs vs Exp
PICs_data.overall <- PICs_RNA_df[,PICs_ct_Aged.overall]
PICs_counts_per_cell.overall <- PICs_UMI[PICs_ct_Aged.overall]
data_PICs_obs.overall <- log1p(PICs_data.overall/PICs_counts_per_cell.overall*1e6)

# Calculate the number of non-zero entries for each gene
non_zero_counts.overall <- rowSums(data_PICs_obs.overall > 0)
# Create a logical vector to identify genes with more than 3 non-zero entries; and min.pct = 0.01
genes_with_non_zero.overall <- non_zero_counts.overall > 3 & non_zero_counts.overall > ncol(data_PICs_obs.overall)*0.01
# Subset the original data frame to include only those genes
data_PICs_obs.overall <- data_PICs_obs.overall[genes_with_non_zero.overall, ]

gene_filtered.overall <- intersect(rownames(data_AD_exp_A.overall),rownames(data_PICs_obs.overall))
data_AD_exp <- data_AD_exp_A.overall[gene_filtered.overall,] # assign data_AD_exp_A to data_AD_exp for parallel computeStats function
data_PICs_obs <- data_PICs_obs.overall[gene_filtered.overall,]
# Perform Wilcoxon Rank Sum Test for Each Gene
# Calculate Log2 Fold Change for Each Gene
print("Calculate Log2 Fold Change for Each Gene in Aged")
clusterExport(cl, varlist = c("data_PICs_obs", "data_AD_exp"))

# Execute in parallel
print("Execute in parallel")
results <- parLapply(cl, 1:nrow(data_PICs_obs), computeStats)

# Convert results list to a data frame
results_df <- do.call(rbind, results)
results_df <- data.frame(gene=as.character(results_df[,1]),
                         log2FC=as.numeric(results_df[,2]),
                         p_value=as.numeric(results_df[,3]))
# Adjust p-values
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")

# add the DEG result to a list
ObsvsExp.DEGs.Aged.overall.df <- results_df
write.csv(ObsvsExp.DEGs.Aged.overall.df,"ObsvsExp.DEGs.Aged.overall.csv")





# UMIs: Create a scatter plot with smaller dots all cells
pdf("PICs_UMI.pdf", width = 10, height = 4)
plot(x = seq_along(PICs_UMI), y = log2(PICs_UMI), 
     main = paste0(length(PICs_UMI)," cells"), 
     xlab = "", ylab = "UMIs(log)", 
     pch = 19,  # solid circles
     cex = 0.15) # smaller dots
dev.off()

pdf("T_UMI.pdf", width = 10, height = 4)
plot(x = seq_along(T_UMI), y = T_UMI, 
     main = paste0(length(T_UMI)," cells"), 
     xlab = "", ylab = "UMIs(log)", 
     pch = 19,  # solid circles
     cex = 0.15) # smaller dots
dev.off()

pdf("B_UMI.pdf", width = 10, height = 4)
plot(x = seq_along(B_UMI), y = B_UMI, 
     main = paste0(length(B_UMI)," cells"), 
     xlab = "", ylab = "UMIs(log)", 
     pch = 19,  # solid circles
     cex = 0.15) # smaller dots
dev.off()