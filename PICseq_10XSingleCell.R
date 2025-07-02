# require "glmGamPoi", version="1.14.0" for Seurat v5 Sctransform v2 under R 4.2.3

## double check the R library path
# .libPaths("/opt/home/buckcenter.org/fwu/R/x86_64-pc-linux-gnu-library/4.2")
# .libPaths()

### Setup the Seurat objects ###
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggsci)
library(sctransform)
library(future)
library(dplyr)

# # change the current plan to access parallelization
# plan("multisession", workers = 4)
# options(future.globals.maxSize= +Inf) # increase maxSize of future function

#dir.create("revision")
setwd("~/PICseq/revision")
# import 10X results
# T and B single cells
sample_list <- c("Young_T_Cells","Young_B_Cells","Aged_T_Cells","Aged_B_Cells")
Sample.list.cmb <- list()
for (s in sample_list){
  sample <- Read10X(data.dir = paste0("/data/array2/fwu/saad/",s))
  sample <- CreateSeuratObject(sample,project=s)
  Sample.list.cmb[[s]] <- sample
}

# create seuratobject for QC
Sample.list.qc.cmb <- Sample.list.cmb

# function for Pre-process Seurat object: QC, PCA and UMAP
SCT_m<-function(l){
  l <- PercentageFeatureSet(l,pattern = "^MT-", col.name = "percent.mt")
  l <- subset(l, subset= nFeature_RNA>250 & percent.mt < 10)
  # recognize BCR and TCR genes
  pattern <- "Ig[hkl]v|Ig[kl]j|Ig[kl]c|Igh[adegm]|^Tr[abdg][vjc]"
  mat <- l[["RNA"]]$counts
  gene2keep <- grep(pattern, rownames(mat), value=T, invert = T)
  l <- subset(l, features = gene2keep) # remove BCR and TCR genes
  l <- SCTransform(l, vst.flavor = "v2", method = "glmGamPoi",verbose = FALSE)
  l <- RunPCA(l, npcs = 20, verbose = FALSE)
  l <- RunUMAP(l,dims = 1:20, verbose = FALSE)
  l <- FindNeighbors(l, dims = 1:20, verbose = FALSE)
  l <- FindClusters(l, verbose = FALSE)
}

# SCTransform for each dataset independently
library(parallel)
Sample.list.qc.cmb<-mclapply(Sample.list.qc.cmb,SCT_m, mc.cores=4)

# Remove doublet
library(DoubletFinder)
for (i in names(Sample.list.qc.cmb)){
  # pK Identification (no ground-truth)
  sweep.res.list_sample <- paramSweep(Sample.list.qc.cmb[[i]], PCs = 1:30, sct = T)
  sweep.stats_sample <- summarizeSweep(sweep.res.list_sample, GT = FALSE)
  bcmvn_sample <- find.pK(sweep.stats_sample)
  pK<-as.numeric(as.character(bcmvn_sample$pK))[bcmvn_sample$BCmetric==max(bcmvn_sample$BCmetric)]
  ## Homotypic Doublet Proportion Estimate
  annotations <- Sample.list.qc.cmb[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(Sample.list.qc.cmb[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies
  Sample.list.qc.cmb[[i]] <- doubletFinder(Sample.list.qc.cmb[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  re.pAnn<-names(Sample.list.qc.cmb[[i]]@meta.data)[(length(Sample.list.qc.cmb[[i]]@meta.data)-1)]
  Sample.list.qc.cmb[[i]] <- doubletFinder(Sample.list.qc.cmb[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = re.pAnn, sct = T)
}

save.image(file = "PICseq_Seurat_revision.RData")
# filter samples based on results of QC & Doublet
for (i in names(Sample.list.qc.cmb)){
  sample.meta.data<-Sample.list.qc.cmb[[i]]@meta.data
  singlet<-row.names(sample.meta.data)[sample.meta.data[length(sample.meta.data)]=="Singlet"]
  # Get the indices of the desired columns
  column_indices <- which(colnames(Sample.list.cmb[[i]]) %in% singlet)
  # Select the columns from the sparse matrix
  Sample.list.cmb[[i]]<-Sample.list.cmb[[i]][,column_indices]
}
# create seuratobject for integration
for (i in names(Sample.list.cmb)){
  Sample.list.cmb[[i]] <- CreateSeuratObject(Sample.list.cmb[[i]], project=i)
  Sample.list.cmb[[i]] <- SCTransform(Sample.list.cmb[[i]], vst.flavor = "v2", method = "glmGamPoi",verbose = F)
}
# create SCT seuratobjects for integration
for (i in names(Sample.list.cmb)){
  # recognize BCR and TCR genes
  pattern <- "Ig[hkl]v|Ig[kl]j|Ig[kl]c|Igh[adegm]|^Tr[abdg][vjc]"
  mat <- Sample.list.cmb[[i]][["RNA"]]$counts
  gene2keep <- grep(pattern, rownames(mat), value=T, invert = T)
  Sample.list.cmb[[i]] <- subset(Sample.list.cmb[[i]], features = gene2keep) # remove BCR and TCR genes
  Sample.list.cmb[[i]] <- SCTransform(Sample.list.cmb[[i]], vst.flavor = "v2", method = "glmGamPoi",verbose = F)
}


# Analysis T and B cell separately

T_cell.list.cmb <- Sample.list.cmb[c("Young_T_Cells","Aged_T_Cells")]
B_cell.list.cmb <- Sample.list.cmb[c("Young_B_Cells","Aged_B_Cells")]

### --- Integrate the T cells --- ###
T.features.cmb <- SelectIntegrationFeatures(object.list = T_cell.list.cmb, nfeatures = 3000)
T_cell.list.cmb <- PrepSCTIntegration(object.list = T_cell.list.cmb, anchor.features = T.features.cmb)

# find the anchor
T.Anchors <- FindIntegrationAnchors(object.list = T_cell.list.cmb, normalization.method = "SCT",
                                  anchor.features = T.features.cmb)
T_cell.combined <- IntegrateData(anchorset = T.Anchors, normalization.method = "SCT")

# Dimensional reduction
T_cell.combined <- RunPCA(T_cell.combined, verbose = FALSE)
# ElbowPlot
ElbowPlot(T_cell.combined,ndims = 50)
T_cell.combined <- RunUMAP(T_cell.combined, reduction = "pca", dims = 1:20)
# Cluster the cells
T_cell.combined <- FindNeighbors(T_cell.combined, dims = 1:20)
for (res in c(0.1, 0.2, 0.8)){
  T_cell.combined <- FindClusters(T_cell.combined, resolution = res)
}
colnames(T_cell.combined@meta.data)
apply(T_cell.combined@meta.data[,grep("snn_res",colnames(T_cell.combined@meta.data))],2,table)
Idents(T_cell.combined) <- "integrated_snn_res.0.2"

# Prepare object to run differential expression on SCT assay with multiple models
T_cell.combined <- PrepSCTFindMarkers(T_cell.combined)
# find markers for every cluster compared to all remaining cells
T_cell.markers.all <- FindAllMarkers(T_cell.combined, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25)

write.csv(T_cell.markers.all,"T_cell.markers.all_revision.csv")

Idents(T_cell.combined) <- "orig.ident"
T.overall.DEG.singlet.age <- FindMarkers(T_cell.combined, assay = "SCT", ident.1 ="Aged_T_Cells", ident.2 ="Young_T_Cells",
                                         min.cells.group = 3,min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
write.csv(T.overall.DEG.singlet.age, "Overall_T_Cells_AgedvsYoung.csv")

# find subclust for T cell
Idents(T_cell.combined) <- "integrated_snn_res.0.2"
T_cell.combined <- FindSubCluster(T_cell.combined,cluster=0,graph.name = "integrated_snn",resolution=0.1)
Idents(T_cell.combined) <- "sub.cluster"
T_cell.combined <- FindSubCluster(T_cell.combined,cluster=1,graph.name = "integrated_snn",resolution=0.1)
table(T_cell.combined$sub.cluster)
Idents(T_cell.combined) <- "sub.cluster"
T_cell.combined <- FindSubCluster(T_cell.combined,cluster=7,graph.name = "integrated_snn",resolution=0.05)
table(T_cell.combined$sub.cluster)
Idents(T_cell.combined) <- "sub.cluster"

# annotate T cell types
T_cell.combined$celltype <- case_when(T_cell.combined$sub.cluster ==
                                        "0_0" ~ "0 - CD8+ Ex T Cells",
                                      T_cell.combined$sub.cluster ==
                                        "0_1" ~ "1 - CD4+ Ex T Cells",
                                      T_cell.combined$sub.cluster ==
                                        "1" ~ "2 - CD8+ EM T Cells",
                                      T_cell.combined$sub.cluster ==
                                        "2" ~ "3 - CD4+ Naive T Cells",
                                      T_cell.combined$sub.cluster ==
                                        "3" ~ "4 - CD4+ EM T Cells",
                                      T_cell.combined$sub.cluster ==
                                        "4" ~ "5 - TFH Cells",
                                      T_cell.combined$sub.cluster ==
                                        "5" ~ "6 - CD8+ Naive T Cells",
                                      T_cell.combined$sub.cluster ==
                                        "6" ~ "7 - Tregs",
                                      T_cell.combined$sub.cluster ==
                                        "7" ~ "8 - CD8+ CM T Cells",
                                      T_cell.combined$sub.cluster ==
                                        "8" ~ "9 - Myeloid-like T Cells",
                                      T_cell.combined$sub.cluster ==
                                        "9" ~ "10 - Ki67+ T Cells",
                                      T_cell.combined$sub.cluster ==
                                        "10" ~ "11 - NKT Cells",
                                      T_cell.combined$sub.cluster ==
                                        "11" ~ "12 - CD4+ IFITM+ CISH+ T Cells",
                                      .default = as.character(T_cell.combined$sub.cluster))

T_cell.combined$celltype <- factor(T_cell.combined$celltype, levels = c("0 - CD8+ Ex T Cells","1 - CD4+ Ex T Cells","2 - CD8+ EM T Cells","3 - CD4+ Naive T Cells",
                                                                        "4 - CD4+ EM T Cells","5 - TFH Cells","6 - CD8+ Naive T Cells",
                                                                        "7 - Tregs","8 - CD8+ CM T Cells","9 - Myeloid-like T Cells",
                                                                        "10 - Ki67+ T Cells","11 - NKT Cells",
                                                                        "12 - CD4+ IFITM+ CISH+ T Cells"))
#celltype_backup <- T_cell.combined$celltype

# further find subclust for T cell cluster - 0
Idents(T_cell.combined) <- "celltype"
T_cell.combined <- FindSubCluster(T_cell.combined,cluster="0 - CD8+ Ex/Mem T Cells",graph.name = "integrated_snn",resolution=0.075)
table(T_cell.combined$sub.cluster)
Idents(T_cell.combined) <- "sub.cluster"
T_cell.markers.all.anno <- FindAllMarkers(T_cell.combined, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(T_cell.markers.all.anno,"T_cell.markers.all.annotation_revision.csv")

# findmarkers for all clusters
Idents(T_cell.combined) <- "celltype"
T_cell.markers.all.anno <- FindAllMarkers(T_cell.combined, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(T_cell.markers.all.anno,"T_cell.markers.all.annotation_revision.csv")




T_cell.combined$age_cell <- paste0(T_cell.combined$celltype,"_",T_cell.combined$orig.ident)
Idents(T_cell.combined) <- "age_cell"
# DEGs for each cell types in T cells
dir.create("Singlet_T_DEGs")
for (c in unique(T_cell.combined$celltype)){
  try({
    DEG.singlet.age <- FindMarkers(T_cell.combined, assay = "SCT", ident.1 =paste0(c,"_Aged_T_Cells"), ident.2 =paste0(c,"_Young_T_Cells"),
                           min.cells.group = 3,min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
    write.csv(DEG.singlet.age, paste0("Singlet_T_DEGs/",c,"_AgedvsYoung_revision.csv"))
  })
}

# top_gene_per_cluster.T
top_gene_per_cluster.T <- T_cell.markers.all %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(T_cell.combined, features = top10$gene,assay = "RNA")
ggsave("T_top10_markers.pdf",width = 12, height = 14)


Idents(T_cell.combined) <- "celltype"
DimPlot(T_cell.combined, reduction = "umap",label = T)
ggsave("umap_T_cell_clusters_anno.pdf", width=10, height=9)

# Cell type ratio change stack barplot
sub.prop.all<-data.frame()
for (l in unique(T_cell.combined$orig.ident)){
  sub.treat<-T_cell.combined@meta.data[T_cell.combined$orig.ident==l,]
  sub.prop<-data.frame(table(sub.treat$celltype)/sum(table(sub.treat$celltype)))
  sub.prop$sample<-l
  sub.prop.all<-rbind(sub.prop.all,sub.prop)
}

# cell type ratio change
ggplot(sub.prop.all, aes(x = Var1, y = Freq, fill = sample)) +
  geom_bar(stat = "identity",position = "dodge") +
  ggtitle("Cell type frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(x="Cell types", y="Proportion", title="Cell type frequency") +
  guides(fill=guide_legend(title='Groups')) +
  scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave("Cell_type_frequency_T.pdf",width=10,height=6)

library(dplyr)
# All sample groups: cell type ratio change stack barplot
## set the levels in order we want
sub.prop.all$sample<-factor(sub.prop.all$sample, 
                            levels=c("Young_T_Cells","Aged_T_Cells"))
## prep cumulative sums for line links
cell.stack.wide1<-sub.prop.all %>%
  pivot_wider(names_from=sample, values_from=Freq) %>% 
  arrange(by=desc(Var1)) %>% 
  replace(is.na(.), 0) %>% # remove NA
  mutate(y=cumsum(`Young_T_Cells`),
         yend=cumsum(Aged_T_Cells))
## set the levels in order we want
sub.prop.all$labs<-round(sub.prop.all$Freq,3) # prepare cell proportion label
sub.prop.all$labs[sub.prop.all$labs<0.05]<-""
ggplot(sub.prop.all, aes(fill=Var1, y=Freq, x=sample)) + 
  geom_bar(stat="identity",width=0.5, col="black") +
  geom_segment(data = cell.stack.wide1, aes(x=1.25, xend=1.75, y=y, yend=yend))+
  geom_text(aes(label= labs), size = 3, position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(x="Sample", y="Proportion", title="T cells") +
  guides(fill=guide_legend(title='Cell types'))
ggsave("T_cells_stacked_barplot.pdf",width=6,height=5)

# check T cell marker genes
# DefaultAssay(T_cell.combined) <- "SCT"
# T_cell.combined$sub.cluster <- factor(T_cell.combined$sub.cluster,
#                                       levels=c("0_0","0_1","1_0","1_1","1_2","2","3","4","5","6","7_0","7_1","8","9","10"))
# Idents(T_cell.combined) <- "sub.cluster"
# manually selected marker genes for T cells
T_markers <- unique(c(t_cell_markers_mouse <- c(
  "Cd8a", "Gzmk", "Cd38", "Pdcd1", "Ifng", "Tox", "Tigit", "Cd44", "Sell", "Fas", "Ccr7", "Ly6c1", "Dapl1", "Igfbp4", "Ccl3", "Ccl4", "Ccl5",
  "Cd4", "Gzmk", "Cd38", "Pdcd1", "Ifng", "Tox", "Tigit", "Cd44", "Sell", "Fas", "Ccr7", "Ly6c1", "Dapl1", "Igfbp4", "Ccl3", "Ccl4", "Ccl5",
  "Cd8b1", "Cd44", "Gzmk", "Cd38", "Ifng", "Tox", "Tigit", "Fas", "Ccl5", "Ccr7", "Dapl1", "Izumo1r", "Igfbp4", "Cx3cr1", "Tbx21",
  "Cd4", "Cd44", "Ccr7", "Sell", "Igfbp4",
  "Cd4", "Cd44", "Ifng", "Sell", "Pdcd1", "Maf",
  "Cd4", "Il21", "Tox2", "Pdcd1", "Maf", "Bcl6",
  "Cd8a", "Ccr7", "Sell", "Cd44",
  "Cd4", "Foxp3", "Areg", "Il2ra", "Ikzf2", "Ctla4",
  "Cd8a", "Cd44", "Sell",
  "Apoe", "C1qb", "H2-Aa", "Il1b",
  "Mki67", "Top2a", "Pclaf",
  "Klrb1c", "Zbtb16",
  "Ifitm1", "Ifitm2", "Ifitm3", "Cish", "Cd4"
)
))
Idents(T_cell.combined) <- "CT_code"
VlnPlot(T_cell.combined, features = T_markers,
        #split.by="orig.ident",split.plot = T,
        alpha = 0.2, assay= "SCT")
ggsave("T_cells_marker_VlnPlot_0.2_subs_revision2.pdf",width=20,height=26)

FeaturePlot(T_cell.combined, features = T_markers,
            #split.by = "orig.ident",
            ncol=4)
ggsave("T_cells_marker_FeaturePlot_0.2_subs2.pdf",width=20,height=40,limitsize = FALSE)

VlnPlot(T_cell.combined, features = c("Eomes", "Ccl3", "Ccl4", "Cish"),
        ncol=2,
        alpha = 0.2, assay= "SCT",
)
ggsave("T_cells_marker_VlnPlot_0.2_sub3_0311.pdf",width=10,height=8)
FeaturePlot(T_cell.combined, features = c("Eomes", "Ccl3", "Ccl4", "Cish"),
            ncol=2)
ggsave("T_cells_marker_FeaturePlot_0.2_sub3_0311.pdf",width=10,height=9)


# log normalization and scale data
DefaultAssay(T_cell.combined) <- "RNA"
T_cell.combined <- NormalizeData(T_cell.combined)
T_cell.combined <- ScaleData(T_cell.combined)


### --- Integrate the B cells --- ###
B.features.cmb <- SelectIntegrationFeatures(object.list = B_cell.list.cmb, nfeatures = 3000)
B_cell.list.cmb <- PrepSCTIntegration(object.list = B_cell.list.cmb, anchor.features = B.features.cmb)

# find the anchor
B.Anchors <- FindIntegrationAnchors(object.list = B_cell.list.cmb, normalization.method = "SCT",
                                    anchor.features = B.features.cmb)
B_cell.combined <- IntegrateData(anchorset = B.Anchors, normalization.method = "SCT")

# Dimensional reduction
B_cell.combined <- RunPCA(B_cell.combined, verbose = FALSE)
# ElbowPlot
ElbowPlot(B_cell.combined,ndims = 50)
B_cell.combined <- RunUMAP(B_cell.combined, reduction = "pca", dims = 1:20)
# Cluster the cells
B_cell.combined <- FindNeighbors(B_cell.combined, dims = 1:20)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8)){
  B_cell.combined <- FindClusters(B_cell.combined, resolution = res)
}
colnames(B_cell.combined@meta.data)
apply(B_cell.combined@meta.data[,grep("snn_res",colnames(B_cell.combined@meta.data))],2,table)
Idents(B_cell.combined) <- "integrated_snn_res.0.2"

# Prepare object to run differential expression on SCT assay with multiple models
B_cell.combined <- PrepSCTFindMarkers(B_cell.combined)
# find markers for every cluster compared to all remaining cells
B_cell.markers.all <- FindAllMarkers(B_cell.combined, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25)

write.csv(B_cell.markers.all,"B_cell.markers.all.csv")

B_cell.markers.C3 <- FindMarkers(B_cell.combined, assay = "SCT", ident.1 ="3_0", ident.2 ="3_1")
write.csv(B_cell.markers.C3, "B_cell.markers.C3.csv")

# find subclust for B cell
B_cell.combined <- FindSubCluster(B_cell.combined,cluster=3,graph.name = "integrated_snn",resolution=0.1)
Idents(B_cell.combined) <- "sub.cluster"
DimPlot(B_cell.combined,split.by = "orig.ident",label = T) # split by age
ggsave("umap_split_age_B_cell_clusters_sub3_revision.pdf", width=14, height=7)

# remove cluster 1 as it represent unhealthy cells
# B_cell.combined0 <- B_cell.combined
# Idents(B_cell.combined)
B_cell.combined <- subset(x = B_cell.combined, idents = c("1"), invert = TRUE)

# annotate B cell types
B_cell.combined$celltype <- case_when(B_cell.combined$sub.cluster ==
                                        "0" ~ "0 - FOBs",
                                      B_cell.combined$sub.cluster ==
                                        "2" ~ "1 - MZBs",
                                      B_cell.combined$sub.cluster ==
                                        "3_0" ~ "2 - ABCs",
                                      B_cell.combined$sub.cluster ==
                                        "3_1" ~ "3 - B1 Cells (1)",
                                      B_cell.combined$sub.cluster ==
                                        "4" ~ "4 - B1 Cells (2)",
                                      B_cell.combined$sub.cluster ==
                                        "5" ~ "5 - SRMhi B cells",
                                      B_cell.combined$sub.cluster ==
                                        "6" ~ "6 - GC B Cells",
                                      B_cell.combined$sub.cluster ==
                                        "7" ~ "0 - FOBs",
                                      B_cell.combined$sub.cluster ==
                                        "8" ~ "7 - EF PBs")

B_cell.combined$celltype <- factor(B_cell.combined$celltype, levels = c("0 - FOBs", "1 - MZBs",
                                                                        "2 - ABCs","3 - B1 Cells (1)","4 - B1 Cells (2)",
                                                                        "5 - SRMhi B cells","6 - GC B Cells","7 - EF PBs"))

DimPlot(B_cell.combined, reduction = "umap", group.by = "celltype",label=T)
ggsave("UMAP_B_revision.pdf",width=10, height=7)

save.image(file = "PICseq_Seurat_revision.RData")


# split column to code cell types
library(tidyr)
T_cell.combined@meta.data <- T_cell.combined@meta.data %>% 
  separate(celltype, into=c("CT_code", "CT_name"), sep=" - ", convert = T, remove = F)
B_cell.combined@meta.data <- B_cell.combined@meta.data %>% 
  separate(celltype, into=c("CT_code", "CT_name"), sep=" - ", convert = T, remove = F)

T_cell.combined$CT_code <- factor(T_cell.combined$CT_code, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12"))

B_cell.combined$CT_code <- factor(B_cell.combined$CT_code, levels = c("0","1","2","3","4","5","6","7"))


# findmarkers for all clusters
Idents(B_cell.combined) <- "celltype"
B_cell.combined <- PrepSCTFindMarkers(B_cell.combined)
B_cell.markers.all.anno <- FindAllMarkers(B_cell.combined, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(B_cell.markers.all.anno,"B_cell.markers.all.annotation_revision.csv")

B_markers <- unique(c("Fcer2a", "Fchsd2", "Icosl", "Ms4a4c", "Zfp318",
                      "Asb2", "Myof", "Cd1d1", "S1pr3", "Pik3r4", "Dusp16",
                      "Zeb2", "Adgre1", "Itgax", "Nt5e", "Itgam", "Tbx21",
                      "Ahnak", "Ctla4", "Zbtb32",
                      "Ahnak", "Ahnak2", "Ctla4", "Zbtb32", "Il10", "Mcoln2", "Tox",
                      "Shmt1", "Srm", "Pusl1", "Pus7",
                      "Gcsam", "S1pr2", "Mef2b", "Jchain", "Mki67", "Mybl1", "Aicda", "Slc41a2", "Bcl6", "Zbtb20",
                      "Mki67", "Zbtb20"))

VlnPlot(B_cell.combined, features = B_markers, alpha = 0.2, assay= "SCT")
ggsave("B_cells_marker_VlnPlot_0.2_sub2.pdf",width=18,height=40)

FeaturePlot(B_cell.combined, features = B_markers, ncol=4)
ggsave("B_cells_marker_FeaturePlot_0.2_sub2.pdf",width=15,height=40)

# log normalize and scale data
DefaultAssay(B_cell.combined) <- "RNA"
B_cell.combined <- NormalizeData(B_cell.combined)
B_cell.combined <- ScaleData(B_cell.combined)

save.image(file = "PICseq_Seurat_revision.RData")


# steps for re-annotate the clusters
Idents(T_cell.combined) <- "CT_code"
VlnPlot(T_cell.combined, pt.size=0,features = c("Cd8a", "Gzmk", "Cd38", "Pdcd1", "Ifng", "Tox", "Tigit", "Cd44", "Sell", "Fas", "Ccr7","Ly6c1", "Dapl1", "Igfbp4", "Ccl3", "Ccl4", "Ccl5",
                                                "Cd8b1", "Cd44", "Gzmk", "Cd38", "Ifng", "Tox", "Tigit", "Fas", "Ccl5", "Ccr7", "Dapl1", "Izumo1r", "Igfbp4", "Cx3cr1", "Tbx21",
                                                "Cd4", "Cd44", "Ccr7", "Sell", "Igfbp4",
                                                "Cd4", "Cd44", "Ifng", "Sell", "Pdcd1", "Maf",
                                                "Cd4", "Il21","Tox2", "Pdcd1", "Maf", "Bcl6",
                                                "Cd8a", "Ccr7", "Sell", "Cd44",
                                                "Cd4", "Foxp3", "Areg", "Il2ra", "Ikzf2", "Ctla4",
                                                "Cd8a", "Cd44", "Sell",
                                                "Apoe", "C1qb", "H2-Aa", "Il1b",
                                                "Mki67", "Top2a", "Pclaf",
                                                "Klrb1c", "Zbtb16",
                                                "Ifitm1", "Ifitm2", "Ifitm3", "Cish", "Cd4"))
ggsave("T_marker_genes.pdf",width=20,height=20)
FeaturePlot(T_cell.combined, features = c("Cd8a", "Gzmk", "Cd38", "Pdcd1", "Ifng", "Tox", "Tigit", "Cd44", "Sell", "Fas", "Ccr7","Ly6c1", "Dapl1", "Igfbp4", "Ccl3", "Ccl4", "Ccl5",
                                          "Cd8b1", "Cd44", "Gzmk", "Cd38", "Ifng", "Tox", "Tigit", "Fas", "Ccl5", "Ccr7", "Dapl1", "Izumo1r", "Igfbp4", "Cx3cr1", "Tbx21",
                                          "Cd4", "Cd44", "Ccr7", "Sell", "Igfbp4",
                                          "Cd4", "Cd44", "Ifng", "Sell", "Pdcd1", "Maf",
                                          "Cd4", "Il21","Tox2", "Pdcd1", "Maf", "Bcl6",
                                          "Cd8a", "Ccr7", "Sell", "Cd44",
                                          "Cd4", "Foxp3", "Areg", "Il2ra", "Ikzf2", "Ctla4",
                                          "Cd8a", "Cd44", "Sell",
                                          "Apoe", "C1qb", "H2-Aa", "Il1b",
                                          "Mki67", "Top2a", "Pclaf",
                                          "Klrb1c", "Zbtb16",
                                          "Ifitm1", "Ifitm2", "Ifitm3", "Cish", "Cd4"))



Idents(B_cell.combined) <- "CT_code"
VlnPlot(B_cell.combined, pt.size=0,features = c("Fcer2a", "Fchsd2", "Icosl", "Ms4a4c", "Zfp318",
                                      "Asb2", "Myof", "Cd1d1", "S1pr3", "Pik3r4", "Dusp16",
                                      "Zeb2", "Adgre1", "Itgax", "Nt5e", "Itgam", "Tbx21",
                                      "Ahnak", "Ctla4", "Zbtb32",
                                      "Ahnak", "Ahnak2", "Ctla4", "Zbtb32", "Il10", "Mcoln2", "Tox",
                                      "Shmt1", "Srm", "Pusl1", "Pus7",
                                      "Gcsam", "S1pr2", "Mef2b", "Jchain", "Mki67", "Mybl1", "Aicda", "Slc41a2","Bcl6", "Zbtb20",
                                      "Mki67", "Zbtb20"))
ggsave("B_marker_genes.pdf",width=20,height=20)
FeaturePlot(B_cell.combined, features = c("Fcer2a", "Fchsd2", "Icosl", "Ms4a4c", "Zfp318",
                                          "Asb2", "Myof", "Cd1d1", "S1pr3", "Pik3r4", "Dusp16",
                                          "Zeb2", "Adgre1", "Itgax", "Nt5e", "Itgam", "Tbx21",
                                          "Ahnak", "Ctla4", "Zbtb32",
                                          "Ahnak", "Ahnak2", "Ctla4", "Zbtb32", "Il10", "Mcoln2", "Tox",
                                          "Shmt1", "Srm", "Pusl1", "Pus7",
                                          "Gcsam", "S1pr2", "Mef2b", "Jchain", "Mki67", "Mybl1", "Aicda", "Slc41a2","Bcl6", "Zbtb20",
                                          "Mki67", "Zbtb20"))


B1_compare <- FindMarkers(B_cell.combined,ident.1="4",ident.2 = "5")
write.csv(B1_compare,"B1_compare_4vs5.csv")


## go to PICseq_script.R


















B_cell.combined$age_cell <- paste0(B_cell.combined$celltype,"_",B_cell.combined$orig.ident)
Idents(B_cell.combined) <- "age_cell"
# DEGs for each cell types in B cells
dir.create("Singlet_B_DEGs")
for (c in unique(B_cell.combined$celltype)){
  try({
    DEG.singlet.age <- FindMarkers(B_cell.combined, assay = "SCT", ident.1 =paste0(c,"_Aged_B_Cells"), ident.2 =paste0(c,"_Young_B_Cells"),
                                   min.cells.group = 3,min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
    write.csv(DEG.singlet.age, paste0("Singlet_B_DEGs/",c,"_AgedvsYoung_revision.csv"))
  })
}

Idents(B_cell.combined) <- "orig.ident"
B.overall.DEG.singlet.age <- FindMarkers(B_cell.combined, assay = "SCT", ident.1 ="Aged_B_Cells", ident.2 ="Young_B_Cells",
                               min.cells.group = 3,min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")
write.csv(B.overall.DEG.singlet.age, "Overall_B_Cells_AgedvsYoung.csv")
# 
# # volcano plot
library(EnhancedVolcano)
dir.create("Singlet_B_Volcano")
B_cells_DEG <- list.files("/opt/home/buckcenter.org/fwu/PICseq/Singlet_B")
for(ct in B_cells_DEG){
  B_CT_DEG <- read.csv(paste0("/opt/home/buckcenter.org/fwu/PICseq/Singlet_B/",ct))
  B_CT_DEG <- tibble::column_to_rownames(B_CT_DEG, "X")
  ct_name <- gsub("_AgedvsYoung.csv","",ct)
  vol.plot <- EnhancedVolcano(B_CT_DEG,
                              lab = row.names(B_CT_DEG),
                              #selectLab=topgenes,
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              #xlim =c(-1.5,2),
                              ylab = bquote(~-Log[10] ~ italic(adj.P)),
                              legendLabels = c("NS", expression(Log[2] ~ FC), "adj.p-value", expression(adj.p - value ~ and
                                                                                                        ~ log[2] ~ FC)),
                              FCcutoff = 0.5,
                              #pCutoff = 0.05,
                              labSize = 5,
                              drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                              #max.overlaps = 20,
                              title = ct_name,
                              subtitle = bquote(italic("Aged vs Young")))
  ggsave(plot=vol.plot,paste0("Singlet_B_Volcano/",ct,"_AgedvsYoung_Volcano.pdf"), width=9,height=8)
}
# 
# 
top_gene_per_cluster.B <- B_cell.markers.all %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(B_cell.combined, features = top10$gene,assay = "RNA")
ggsave("B_top10_markers.pdf",width = 12, height = 14)
# 
# 
# # arrange the order of sample
# B_cell.combined$orig.ident <- factor(B_cell.combined$orig.ident,
#                                      levels=c("Young_B_Cells","Aged_B_Cells"))
# # save the UMAP
# DimPlot(B_cell.combined, reduction = "umap",label = T)
# ggsave("umap_B_cell_clusters.pdf", width=8, height=7)

Idents(B_cell.combined) <- "celltype"
DimPlot(B_cell.combined, reduction = "umap",label = T)
ggsave("umap_B_cell_clusters_anno.pdf", width=10, height=9)

# DimPlot(B_cell.combined,split.by = "orig.ident",label = T) # split by age
# ggsave("umap_split_age_B_cell_clusters_0.2.pdf", width=14, height=7)
# # Idents(B_cell.combined) <- "integrated_snn_res.0.8" # res=0.8
# # DimPlot(B_cell.combined,split.by = "orig.ident",label = T) # split by age
# # ggsave("umap_split_age_B_cell_clusters_0.8.pdf", width=14, height=7)
# 
# # find subclust for B cell
# B_cell.combined <- FindSubCluster(B_cell.combined,cluster=3,graph.name = "integrated_snn",resolution=0.1)
# Idents(B_cell.combined) <- "sub.cluster"
# DimPlot(B_cell.combined,split.by = "orig.ident",label = T) # split by age
# ggsave("umap_split_age_B_cell_clusters_sub3.pdf", width=14, height=7)

# update all makers after subclustering
B_cell.markers.all <- FindAllMarkers(B_cell.combined, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(B_cell.markers.all,"B_cell.markers.all.csv")

# export the Seurat objects
saveRDS(T_cell.combined,"T_cell.combined.rds")
saveRDS(B_cell.combined,"B_cell.combined.rds")

# Cell type ratio change stack barplot
sub.prop.all<-data.frame()
for (l in unique(B_cell.combined$orig.ident)){
  sub.treat<-B_cell.combined@meta.data[B_cell.combined$orig.ident==l,]
  sub.prop<-data.frame(table(sub.treat$celltype)/sum(table(sub.treat$celltype)))
  sub.prop$sample<-l
  sub.prop.all<-rbind(sub.prop.all,sub.prop)
}

# cell type ratio change
ggplot(sub.prop.all, aes(x = Var1, y = Freq, fill = sample)) +
  geom_bar(stat = "identity",position = "dodge") +
  ggtitle("Cell type frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(x="Cell types", y="Proportion", title="Cell type frequency") +
  guides(fill=guide_legend(title='Groups'))
ggsave("Cell_type_frequency_B.pdf",width=10,height=6)

library(dplyr)
# All sample groups: cell type ratio change stack barplot
## set the levels in order we want
sub.prop.all$sample<-factor(sub.prop.all$sample, 
                            levels=c("Young_B_Cells","Aged_B_Cells"))
## prep cumulative sums for line links
cell.stack.wide1<-sub.prop.all %>%
  pivot_wider(names_from=sample, values_from=Freq) %>% 
  arrange(by=desc(Var1)) %>% 
  replace(is.na(.), 0) %>% # remove NA
  mutate(y=cumsum(`Young_B_Cells`),
         yend=cumsum(Aged_B_Cells))
## set the levels in order we want
sub.prop.all$labs<-round(sub.prop.all$Freq,3) # prepare cell proportion label
sub.prop.all$labs[sub.prop.all$labs<0.05]<-""
ggplot(sub.prop.all, aes(fill=Var1, y=Freq, x=sample)) + 
  geom_bar(stat="identity",width=0.5, col="black") +
  geom_segment(data = cell.stack.wide1, aes(x=1.25, xend=1.75, y=y, yend=yend))+
  geom_text(aes(label= labs), size = 3, position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(x="Sample", y="Proportion", title="B cells") +
  guides(fill=guide_legend(title='Cell types'))
ggsave("B_cells_stacked_barplot.pdf",width=6,height=5)

# check B cell marker genes 1
DefaultAssay(B_cell.combined) <- "SCT"
VlnPlot(B_cell.combined, features = c("Fcer2a", "Cr2", "Icosl", "Fchsd2", "Ms4a4c", "Zfp318", "Asb2", 
                                      "Myof", "Cd1d1", "S1pr3", "Pik3r4", "Iglv1", "Iglv2", "Iglv3", 
                                      "Fam129c", "Vpreb3", "Cd24a", "Zeb2", "Adgre1", "Fcer1g", "Nr4a2", 
                                      "Hes1", "Id3", "Cdkn2a", "Itgax", "Nt5e", "Itgam", "Tbx21", "Ahnak", 
                                      "Ahnak2", "Ctla4", "Zbtb32", "Il10", "Mcoln2", "Apoe", "Tacstd2", 
                                      "Igkv1-135", "Igkv10-96", "Igkv1-110", "Shmt1", "Srm", "Pusl1", 
                                      "Pus7", "Gcsam", "S1pr2", "Mef2b", "Jchain", "Mki67", "Mybl1", 
                                      "Aicda", "Slc41a2", "Bcl6", "Zbtb20", "Ighg2b", "Ighg2c", "Igha", 
                                      "Ighm", "Ighd", "Basp1", "Ccr6", "Pdcd1lg2", "Cd80", "Plxnb2", 
                                      "Cd40", "Igkv1-117", "Lyz2", "Ifitm2", "Ccr2", "Ltb4r1", "Igkv1-99", "Spn", "Cd5"),
        #split.by="orig.ident",split.plot = T,
        alpha = 0.2, assay= "SCT",
        #cols=c("#4DBBD5FF","#E64B35FF")
        )
ggsave("B_cells_marker_VlnPlot.pdf",width=20,height=40)

Idents(B_cell.combined) <- "CT_code"
VlnPlot(B_cell.combined, features = c("Cd19"),alpha = 0.2, assay= "SCT")
ggsave("B_cells_Cd19_VlnPlot.pdf",width=6,height=4)

VlnPlot(B_cell.combined, features = c("Fcer2a"),alpha = 0.2, assay= "SCT")
ggsave("B_cells_Fcer2a_VlnPlot.pdf",width=6,height=4)

VlnPlot(B_cell.combined, features = c("Cr2"),alpha = 0.2, assay= "SCT")
ggsave("B_cells_Cr2_VlnPlot.pdf",width=6,height=4)

VlnPlot(B_cell.combined, features = c("Zeb2"),alpha = 0.2, assay= "SCT")
ggsave("B_cells_Zeb2_VlnPlot.pdf",width=6,height=4)

VlnPlot(B_cell.combined, features = c("Itgax"),alpha = 0.2, assay= "SCT")
ggsave("B_cells_Itgax_VlnPlot.pdf",width=6,height=4)

VlnPlot(B_cell.combined, features = c("Tbx21"),alpha = 0.2, assay= "SCT")
ggsave("B_cells_Tbx21_VlnPlot.pdf",width=6,height=4)



# --- Volcano plot --- #
ggplotThemeBlock = function(slantXLabels = FALSE, theme="classic", display_legend=TRUE) {
  
  if(theme == "bw") {
    themeBlock = theme_bw()
  } else {
    themeBlock = theme_classic()
  }
  
  themeBlock = themeBlock + theme(axis.title.x=element_text(size=22),
                                  axis.text.x=element_text(size=20),
                                  axis.title.y=element_text(size=22),
                                  axis.text.y=element_text(size=20),
                                  strip.text=element_text(size=22),
                                  plot.title=element_text(size=24, hjust=0.5, vjust=0.5),
                                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  
  if(slantXLabels) {
    themeBlock = themeBlock + theme(axis.text.x=element_text(size=20, angle=45, hjust=1, vjust=1))
  } else {
    themeBlock = themeBlock + theme(axis.text.x=element_text(size=20, hjust=0.5))
  }
  
  # Either remove legend completely, or format appropriately
  if(display_legend) {
    themeBlock = themeBlock + theme(legend.title=element_text(size=22),
                                    legend.text=element_text(size=20),
                                    legend.text.align=0)  
    
  } else {
    themeBlock = themeBlock + theme(legend.position="none")
  }
  
  return(themeBlock)
}

library(ggrepel)
library(RColorBrewer)
volcanoLabelsOverall = read.delim("/opt/home/buckcenter.org/fwu/PICseq/volcano_plot_genes_global_categories.txt", header=T, as.is=T, check.names=F)

markersAgedvsYoung <- read.csv("Overall_B_Cells_AgedvsYoung.csv",as.is=T, check.names=F)
colnames(markersAgedvsYoung)[1]<-"symbol"
minP = min(markersAgedvsYoung %>% dplyr::select(p_val) %>% filter(p_val > 0))
plotDatGlobal = markersAgedvsYoung %>% mutate(p_val = ifelse(p_val == 0, minP, p_val))
plotDatGlobal = left_join(plotDatGlobal, volcanoLabelsOverall, by="symbol")
plotDatGlobal = plotDatGlobal %>% mutate(highlight=ifelse(is.na(category), "small", "large"))
plotDatGlobal$highlight = factor(plotDatGlobal$highlight, levels=c("small", "large"))
plotDatHighlight = plotDatGlobal %>% filter(!is.na(category))

categories = sort(unique(plotDatHighlight$category))
#categories = c(setdiff(categories, "NicheNet"), "NicheNet")

plotDatHighlight$category = factor(plotDatHighlight$category, levels=categories)

paletteColors = brewer.pal(7, "Accent")[-4]
paletteColors = c(paletteColors[c(6)], "orange1", "cadetblue", paletteColors[1])
names(paletteColors) = categories

# logFcBreaks = c(-2,-log2(3), -1, -log2(1.5), 0, log2(1.5), 1)
# logFcLabels = c(-4, -2, -1.5, "same", 1.5, 2)

upInAged = "Up Aged"
upInYoung = "Up Young"

pl = ggplot(data=plotDatGlobal, aes(x=avg_log2FC, y=-log10(p_val)))
pl = pl + geom_point(color="lightgrey", size=2)
pl = pl + geom_point(data=plotDatHighlight, aes(color=category), size=4)
pl = pl + geom_text_repel(data=plotDatHighlight, aes(label=symbol, color=category), fontface="bold", force=70, max.overlaps=200, size=6, segment.linetype=2)
pl = pl + scale_y_continuous(limits=c(0, 450), breaks=c(0,50,100,150,200,250,300, 350))
#pl = pl + scale_x_continuous(limits=c(-2.5, 1.3), breaks=logFcBreaks, labels=logFcLabels)
pl = pl + geom_vline(xintercept=0, linetype="dotted")
pl = pl + annotate(geom="text", label=upInYoung, x=-0.5, y=420, size=7)
pl = pl + annotate(geom="text", label=upInAged, parse=F, x=0.5, y=420, size=7)
pl = pl + scale_color_manual(name="Category", values=paletteColors)
pl = pl + ggplotThemeBlock()
pl = pl + xlab(expression("log2 Fold Change (Aged / Young)"))
pl = pl + ylab("-log10 P-Value")
pl
ggsave("overall_B_Aged_vs_Young_differential_genes_by_category.pdf", pl, h=9, w=18, useDingbats=FALSE)


# Genes to highlight for each cluster-specific Aged vs Young analysis volcano plot
##volcanoLabelsPerCluster
# DE files per cluster
##AgedVsYoungDE
library(stringr)
volcanoLabelsPerCluster = read.delim("/opt/home/buckcenter.org/fwu/PICseq/volcano_plot_genes_per_cluster.txt", header=T, as.is=T, check.names=F)
AgedVsYoungFiles = dir("/opt/home/buckcenter.org/fwu/PICseq/B_cell_singlet_DEGs", full.names=T)

#tcrColors = brewer.pal(n=3, name="Set2")[2:1]

AgedVsYoungDE = list()
for(currFile in AgedVsYoungFiles) {
  currCluster = str_replace(currFile, ".*B_cell_singlet_DEGs/([0-9]{1,2}).*", "\\1")
  print(currCluster)
  currDE = read.csv(currFile, header=T, as.is=T, check.names=F)
  AgedVsYoungDE[[currCluster]] = currDE
}

for(currCluster in names(AgedVsYoungDE)) {
  
  currMarkers = AgedVsYoungDE[[currCluster]]
  names(currMarkers)[1]<-"symbol"
  currHighlights = volcanoLabelsPerCluster %>% filter(cluster == currCluster) %>% dplyr::select(symbol, direction)
  currentMarkerHighlights = inner_join(currMarkers, currHighlights)
  
  minP = min(currMarkers %>% dplyr::select(p_val) %>% filter(p_val > 0))
  currPlotDat = currMarkers %>% mutate(p_val = ifelse(p_val == 0, minP, p_val))
  currPlotDat = left_join(currPlotDat, volcanoLabelsOverall, by="symbol")
  
  currentMarkerHighlights$direction = factor(currentMarkerHighlights$direction, levels=c("up_Aged", "down_Aged"))
  
  #paletteColorsDuo = c("up_Aged"=tcrColors[2], "down_Aged"=tcrColors[1])
  paletteColorsDuo = c("up_Aged"="#e34a33", "down_Aged"="#636363")
  
  #logFcBreaks = c(-2,-log2(3), -1, -log2(1.5), 0, log2(1.5), 1)
  #logFcLabels = c(-4, -3, -2, -1.5, "same", 1.5, 2)
  
  #  upInKO = "up in muMT ->"
  #  upInWT = "<- up in WT"
  
  pl = ggplot(data=currPlotDat, aes(x=avg_log2FC, y=-log10(p_val)))
  pl = pl + geom_point(color="lightgrey", size=2)
  pl = pl + geom_point(data=currentMarkerHighlights, aes(color=direction), size=6)
  pl = pl + geom_text_repel(data=currentMarkerHighlights, aes(label=symbol, color=direction), fontface="bold", force=70, max.overlaps=200, size=8, segment.linetype=2)
  #  pl = pl + scale_y_continuous(limits=c(0, 400))
  # pl = pl + scale_x_continuous(breaks=logFcBreaks, labels=logFcLabels)
  pl = pl + geom_vline(xintercept=0, linetype="dotted")
  pl = pl + scale_color_manual(name="", values=paletteColorsDuo, labels=c(expression("Up in Aged"), "Up in Young"))
  pl = pl + ggplotThemeBlock()
  pl = pl + xlab(expression("log2 Fold Change (Aged / Young)"))
  pl = pl + ylab("-log10 P-Value")  
  pl
  
  volcanoPlotFilePDF = str_c("volcano_plot_Aged_vs_Young_differential_genes-cluster_", currCluster  ,".pdf")
  ggsave(volcanoPlotFilePDF, pl, h=7, w=10)
}




# # check T cell marker genes 1
# DefaultAssay(T_cell.combined) <- "SCT"
# VlnPlot(T_cell.combined, features = c("Cd3e", "Cd4", "Il21", "Tox2", "Pdcd1", "Maf", "Bcl6", "Foxp3", 
#                                       "Areg", "Il2ra", "Il10", "Ctla4", "Cd8a", "Gzmk", "Eomes", "Cd38", 
#                                       "Tox", "Tigit", "Ifng", "Cd44", "Cx3cr1", "Cd38", "Ccl3", "Ccl4", 
#                                       "Ccl5", "Ccr7", "Sell", "Ccr9", "Lars2", "Malat1", "Apoe", 
#                                       "C1qb", "H2-aa", "Il1b", "Mki67", "Top2a", "Cish", "Ifitm1", 
#                                       "Ifitm3", "Ifitm2", "Ifitm7", "Ifitm1", "Ifitm3", "Ifitm2", "Ifitm7", 
#                                       "Ifitm1", "Ifitm3", "Ifitm2", "Ifitm7", "Klrb1c", "Zbtb16"),
#         #split.by="orig.ident",split.plot = T,
#         alpha = 0.2, assay= "SCT",
#         #cols=c("#4DBBD5FF","#E64B35FF")
# )
# ggsave("T_cells_marker_VlnPlot.pdf",width=20,height=35)
Idents(T_cell.combined) <- "CT_code"

VlnPlot(T_cell.combined, features = c("Cd3e"),alpha = 0.2, assay= "SCT")
ggsave("T_cells_Cd3e_VlnPlot.pdf",width=6,height=4)

VlnPlot(T_cell.combined, features = c("Cd4"),alpha = 0.2, assay= "SCT")
ggsave("T_cells_Cd4_VlnPlot.pdf",width=6,height=4)

VlnPlot(T_cell.combined, features = c("Cd8a"),alpha = 0.2, assay= "SCT")
ggsave("T_cells_Cd8a_VlnPlot.pdf",width=6,height=4)

VlnPlot(T_cell.combined, features = c("Pdcd1"),alpha = 0.2, assay= "SCT")
ggsave("T_cells_Pdcd1_VlnPlot.pdf",width=6,height=4)

VlnPlot(T_cell.combined, features = c("Cd44"),alpha = 0.2, assay= "SCT")
ggsave("T_cells_Cd44_VlnPlot.pdf",width=6,height=4)

VlnPlot(T_cell.combined, features = c("Sell"),alpha = 0.2, assay= "SCT")
ggsave("T_cells_Sell_VlnPlot.pdf",width=6,height=4)

VlnPlot(T_cell.combined, features = c("Gzmk"),alpha = 0.2, assay= "SCT")
ggsave("T_cells_Gzmk_VlnPlot.pdf",width=6,height=4)


# volcano plot for T cells
volcanoLabelsOverall = read.delim("/opt/home/buckcenter.org/fwu/PICseq/volcano_plot_genes_global_categories_T.txt", header=T, as.is=T, check.names=F)

markersAgedvsYoung <- read.csv("Overall_T_Cells_AgedvsYoung.csv",as.is=T, check.names=F)
colnames(markersAgedvsYoung)[1]<-"symbol"
minP = min(markersAgedvsYoung %>% dplyr::select(p_val) %>% filter(p_val > 0))
plotDatGlobal = markersAgedvsYoung %>% mutate(p_val = ifelse(p_val == 0, minP, p_val))
plotDatGlobal = left_join(plotDatGlobal, volcanoLabelsOverall, by="symbol")
plotDatGlobal = plotDatGlobal %>% mutate(highlight=ifelse(is.na(category), "small", "large"))
plotDatGlobal$highlight = factor(plotDatGlobal$highlight, levels=c("small", "large"))
plotDatHighlight = plotDatGlobal %>% filter(!is.na(category))

categories = sort(unique(plotDatHighlight$category))
#categories = c(setdiff(categories, "NicheNet"), "NicheNet")

plotDatHighlight$category = factor(plotDatHighlight$category, levels=categories)

paletteColors = brewer.pal(7, "Accent")
paletteColors = c("orange1","hotpink2","skyblue3",paletteColors[1])
names(paletteColors) = categories

# logFcBreaks = c(-2,-log2(3), -1, -log2(1.5), 0, log2(1.5), 1)
# logFcLabels = c(-4, -2, -1.5, "same", 1.5, 2)

upInAged = "Up Aged"
upInYoung = "Up Young"

pl = ggplot(data=plotDatGlobal, aes(x=avg_log2FC, y=-log10(p_val)))
pl = pl + geom_point(color="lightgrey", size=2)
pl = pl + geom_point(data=plotDatHighlight, aes(color=category), size=4)
pl = pl + geom_text_repel(data=plotDatHighlight, aes(label=symbol, color=category), fontface="bold", force=70, max.overlaps=200, size=6, segment.linetype=2)
pl = pl + scale_y_continuous(limits=c(0, 450), breaks=c(0,50,100,150,200,250,300, 350))
#pl = pl + scale_x_continuous(limits=c(-2.5, 1.3), breaks=logFcBreaks, labels=logFcLabels)
pl = pl + geom_vline(xintercept=0, linetype="dotted")
pl = pl + annotate(geom="text", label=upInYoung, x=-0.6, y=450, size=7)
pl = pl + annotate(geom="text", label=upInAged, parse=F, x=0.6, y=450, size=7)
pl = pl + scale_color_manual(name="Category", values=paletteColors)
pl = pl + ggplotThemeBlock()
pl = pl + xlab(expression("log2 Fold Change (Aged / Young)"))
pl = pl + ylab("-log10 P-Value")
pl
ggsave("overall_T_Aged_vs_Young_differential_genes_by_category_T_cell.pdf", pl, h=9, w=18, useDingbats=FALSE)


# Genes to highlight for each cluster-specific Aged vs Young analysis volcano plot
##volcanoLabelsPerCluster
# DE files per cluster
##AgedVsYoungDE
library(stringr)
volcanoLabelsPerCluster = read.delim("/opt/home/buckcenter.org/fwu/PICseq/volcano_plot_genes_per_cluster_T.txt", header=T, as.is=T, check.names=F)
AgedVsYoungFiles = dir("/opt/home/buckcenter.org/fwu/PICseq/T_cell_singlet_DEGs", full.names=T)

#tcrColors = brewer.pal(n=3, name="Set2")[2:1]

AgedVsYoungDE = list()
for(currFile in AgedVsYoungFiles) {
  currCluster = str_replace(currFile, ".*T_cell_singlet_DEGs/([0-9]{1,2}).*", "\\1")
  print(currCluster)
  currDE = read.csv(currFile, header=T, as.is=T, check.names=F)
  AgedVsYoungDE[[currCluster]] = currDE
}

for(currCluster in names(AgedVsYoungDE)) {
  
  currMarkers = AgedVsYoungDE[[currCluster]]
  names(currMarkers)[1]<-"symbol"
  currHighlights = volcanoLabelsPerCluster %>% filter(cluster == currCluster) %>% dplyr::select(symbol, direction)
  currentMarkerHighlights = inner_join(currMarkers, currHighlights)
  
  minP = min(currMarkers %>% dplyr::select(p_val) %>% filter(p_val > 0))
  currPlotDat = currMarkers %>% mutate(p_val = ifelse(p_val == 0, minP, p_val))
  currPlotDat = left_join(currPlotDat, volcanoLabelsOverall, by="symbol")
  
  currentMarkerHighlights$direction = factor(currentMarkerHighlights$direction, levels=c("up_Aged", "down_Aged"))
  
  #paletteColorsDuo = c("up_Aged"=tcrColors[2], "down_Aged"=tcrColors[1])
  paletteColorsDuo = c("up_Aged"="#e34a33", "down_Aged"="#636363")
  
  #logFcBreaks = c(-2,-log2(3), -1, -log2(1.5), 0, log2(1.5), 1)
  #logFcLabels = c(-4, -3, -2, -1.5, "same", 1.5, 2)
  
  #  upInKO = "up in muMT ->"
  #  upInWT = "<- up in WT"
  
  pl = ggplot(data=currPlotDat, aes(x=avg_log2FC, y=-log10(p_val)))
  pl = pl + geom_point(color="lightgrey", size=2)
  pl = pl + geom_point(data=currentMarkerHighlights, aes(color=direction), size=6)
  pl = pl + geom_text_repel(data=currentMarkerHighlights, aes(label=symbol, color=direction), fontface="bold", force=70, max.overlaps=200, size=8, segment.linetype=2)
  #  pl = pl + scale_y_continuous(limits=c(0, 400))
  # pl = pl + scale_x_continuous(breaks=logFcBreaks, labels=logFcLabels)
  pl = pl + geom_vline(xintercept=0, linetype="dotted")
  pl = pl + scale_color_manual(name="", values=paletteColorsDuo, labels=c(expression("Up in Aged"), "Up in Young"))
  pl = pl + ggplotThemeBlock()
  pl = pl + xlab(expression("log2 Fold Change (Aged / Young)"))
  pl = pl + ylab("-log10 P-Value")  
  pl
  
  volcanoPlotFilePDF = str_c("volcano_plot_Aged_vs_Young_differential_genes-cluster_", currCluster  ,".pdf")
  ggsave(volcanoPlotFilePDF, pl, h=7, w=10)
}




















# # FeaturePlot(B_cell.combined, features = c("Cd19", "Ms4a1", "Cr2", "Fcer2a", "Itgax", "Tbx21", "Zeb2", "Cxcr4", "Fas", "Cd22", "Cd1d1","Ly77"),
# #             split.by = "orig.ident",ncol=4)
# FeaturePlot(B_cell.combined, features = c("Cd19", "Fcer2a", "Cr2", "Tbx21", "Itgax"),
#             ncol=2)
# ggsave("B_cells_marker_FeaturePlot_0.2_sub3.pdf",width=10,height=15)
# 
# 
# VlnPlot(B_cell.combined, features = c("Igha", "Ighg2b", "Ighg2c"),
#         ncol=2,
#         alpha = 0.2, assay= "SCT",
# )
# ggsave("B_cells_marker_VlnPlot_0.2_sub3_0311.pdf",width=10,height=8)
# FeaturePlot(B_cell.combined, features = c("Igha", "Ighg2b", "Ighg2c"),
#             ncol=2)
# ggsave("B_cells_marker_FeaturePlot_0.2_sub3_0311.pdf",width=10,height=9)
# 
# # check B cell marker genes 2
# B_cell.combined$sub.cluster <- factor(B_cell.combined$sub.cluster,
#                                       levels=c("0", "1", "2", "3_0", "3_1", "4", "5", "6", "7", "8", "9","10", "11", "12", "13", "14"))
# Idents(B_cell.combined) <- "sub.cluster"
# DefaultAssay(B_cell.combined) <- "SCT"



# ### update meatadata ###
# # cell type annotations
# T_cell.combined@meta.data$cell_type <- case_when(T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "0_0" ~ "TFH Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "0_1" ~ "Tregs",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "1_0" ~ "CD8+ Ex/Mem T Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "1_1" ~ "CD4+ Ex/Mem T Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "1_2" ~ "CCL3+ CCL4+ Ex/Mem T Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "2" ~ "CD4+ Naive T Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "3" ~ "CD8+ Naive T Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "4" ~ "CD8+ CM T Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "5" ~ "Lars2hi Malat1hi T cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "6" ~ "Trbv3+ Naive T Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "7_0" ~ "Myeloid-like T Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "7_1" ~ "Ki67+ T Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "8" ~ "CD4+ IFITM+ CISH+ T Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "9" ~ "NKT Cells",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "10" ~ "Trbv4+ T Cells")
# 
# B_cell.combined$sub.cluster.char <- as.character(B_cell.combined$sub.cluster)
# B_cell.combined@meta.data$cell_type <- case_when(B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "0" ~ "FOBs",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "1" ~ "MZBs",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "2" ~ "Iglv1-3+ Immature B Cells",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "3_0" ~ "Zeb2lo B Cells",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "3_1" ~ "ABCs",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "4" ~ "B1 Cells",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "5" ~ "Igkv1-135+ B cells",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "6" ~ "Igkv10-96+ B Cells",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "7" ~ "Igkv1-110+ B Cells",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "8" ~ "Srmhi B Cells",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "9" ~ "GC B Cells",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "10" ~ "MBCs",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "11" ~ "Igkv1-117+ B Cells",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "12" ~ "PBs",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "13" ~ "Myeloid-Like B Cells",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "14" ~ "Igkv1-99+ B Cells")
# 
# 
# # Subcluster cell type annotations
# T_cell.combined@meta.data$clusters <- case_when(T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "0_0" ~ "0",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "0_1" ~ "1",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "1_0" ~ "2",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "1_1" ~ "3",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "1_2" ~ "4",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "2" ~ "5",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "3" ~ "6",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "4" ~ "7",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "5" ~ "8",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "6" ~ "9",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "7_0" ~ "10",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "7_1" ~ "11",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "8" ~ "12",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "9" ~ "13",
#                                                  T_cell.combined@meta.data$sub.cluster.char ==
#                                                    "10" ~ "14")
# 
# B_cell.combined@meta.data$clusters <- case_when(B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "0" ~ "0",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "1" ~ "1",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "2" ~ "2",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "3_0" ~ "3",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "3_1" ~ "4",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "4" ~ "5",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "5" ~ "6",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "6" ~ "7",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "7" ~ "8",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "8" ~ "9",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "9" ~ "10",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "10" ~ "11",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "11" ~ "12",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "12" ~ "13",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "13" ~ "14",
#                                                  B_cell.combined@meta.data$sub.cluster.char ==
#                                                    "14" ~ "15")
# 
# 
# # finalized the cluster name
# T_cell.combined$cluster_name <- paste0(T_cell.combined$clusters,"_",T_cell.combined$cell_type)
# B_cell.combined$cluster_name <- paste0(B_cell.combined$clusters,"_",B_cell.combined$cell_type)
# 
# # UMAP for publication
# # T cells
# T_cell.combined$cluster_name <- factor(T_cell.combined$cluster_name,
#                                       levels=c("0_TFH Cells","1_Tregs","2_CD8+ Ex/Mem T Cells","3_CD4+ Ex/Mem T Cells",
#                                                "4_CCL3+ CCL4+ Ex/Mem T Cells","5_CD4+ Naive T Cells","6_CD8+ Naive T Cells",
#                                                "7_CD8+ CM T Cells","8_Lars2hi Malat1hi T cells","9_Trbv3+ Naive T Cells",
#                                                "10_Myeloid-like T Cells","11_Ki67+ T Cells","12_CD4+ IFITM+ CISH+ T Cells",
#                                                "13_NKT Cells","14_Trbv4+ T Cells"))
# Idents(T_cell.combined) <- "cluster_name"
# DimPlot(T_cell.combined, reduction = "umap",label = T)
# ggsave("umap_T_cell_clusters_anno_fin.pdf", width=10, height=9)
# 
# # B cells
# B_cell.combined$cluster_name <- factor(B_cell.combined$cluster_name,
#                                        levels=c("0_FOBs","1_MZBs","2_Iglv1-3+ Immature B Cells",
#                                                 "3_Zeb2lo B Cells","4_ABCs","5_B1 Cells",
#                                                 "6_Igkv1-135+ B cells","7_Igkv10-96+ B Cells",
#                                                 "8_Igkv1-110+ B Cells","9_Srmhi B Cells","10_GC B Cells",
#                                                 "11_MBCs","12_Igkv1-117+ B Cells","13_PBs","14_Myeloid-Like B Cells",
#                                                 "15_Igkv1-99+ B Cells"))
# Idents(B_cell.combined) <- "cluster_name"
# DimPlot(B_cell.combined, reduction = "umap",label = T)
# ggsave("umap_B_cell_clusters_anno_fin.pdf", width=10, height=9)
# 
# # T cells numberic label
# T_cell.combined$clusters <- factor(T_cell.combined$clusters,
#                                        levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"))
# Idents(T_cell.combined) <- "clusters"
# DimPlot(T_cell.combined, reduction = "umap",label = T)
# ggsave("umap_T_cell_clusters_anno_fin_num.pdf", width=7, height=6)
# 
# # B cells numberic label
# B_cell.combined$clusters <- factor(B_cell.combined$clusters,
#                                    levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))
# Idents(B_cell.combined) <- "clusters"
# DimPlot(B_cell.combined, reduction = "umap",label = T)
# ggsave("umap_B_cell_clusters_anno_fin_num.pdf", width=7, height=6)







# ### --- PIC-seq analysis --- ###
# # load MetaCell package
# library("metacell")
# # set database for metacell analysis
# scdb_init("/opt/home/buckcenter.org/fwu/PICseq/PICdb", force_reinit=T)
# # import 10x data
# # singlet cell (T or B cells)
# mcell_import_multi_scmat_10x("sc", "dataset_table_SC.txt", "/data/array2/fwu/saad", force = FALSE)
# # doublet cell
# mcell_import_multi_scmat_10x("db", "dataset_table_DB.txt", "/data/array2/fwu/saad", force = FALSE)
# 
# save.image("/opt/home/buckcenter.org/fwu/PICseq/PICseq_Joint.RData")
# 
# # Before starting to analyze the data, we link the package to a figure directory:
# if(!dir.exists("PIC_figs")) dir.create("PIC_figs/")
# scfigs_init("PIC_figs/")
# 
# ### process singlet: ###
# # drop cells failed Seurat QC
# T.QC <- gsub("_1|_2","",rownames(T_cell.combined@meta.data))
# T.QC <- paste0(T_cell.combined$orig.ident,"_",T.QC)
# B.QC <- gsub("_1|_2","",rownames(B_cell.combined@meta.data))
# B.QC <- paste0(B_cell.combined$orig.ident,"_",B.QC)
# 
# mcell_mat_ignore_cells(new_mat_id = "sc", mat_id = "sc",ig_cells = c(T.QC,B.QC), reverse=T)
# 
# # Exploring and filtering the UMI matrix
# mcell_plot_umis_per_cell("sc")
# 
# sc_mat = scdb_mat("sc")
# print(dim(sc_mat@mat))
# 
# nms = c(rownames(sc_mat@mat), rownames(sc_mat@ignore_gmat))
# ig_genes = c(grep("^Igj", nms, v=T), 
#              grep("^Igh",nms,v=T),
#              grep("^Igk", nms, v=T), 
#              grep("^Igl", nms, v=T))
# bad_genes = unique(c(grep("ERCC", nms, v=T), grep("^mt-", nms, v=T), "Gm42418", "Gm26917", "Gm29216", ig_genes))
# 
# mcell_mat_ignore_genes(new_mat_id="sc", mat_id="sc", bad_genes, reverse=F)
# 
# # selecting feature genes
# mcell_add_gene_stat(gstat_id="sc", mat_id="sc", force=T)
# mcell_gset_filter_varmean(gset_id="sc_feats", gstat_id="sc", T_vm=0.08, force_new=T)
# mcell_gset_filter_cov(gset_id = "sc_feats", gstat_id="sc", T_tot=100, T_top3=2)
# mcell_plot_gstats(gstat_id="sc", gset_id="sc_feats")
# 
# # define cell cycle modules
# anchor_genes = c("Top2a", "Mki67", "Pcna", "Mcm4", "Cdk1", "Ube2c")
# 
# mcell_mat_rpt_cor_anchors(mat_id="sc", gene_anchors = anchor_genes, cor_thresh = 0.1, gene_anti = c(),
#                           tab_fn = "lateral_gmods.txt", sz_cor_thresh = 0.2)
# 
# gcor_mat = read.table("lateral_gmods.txt", header=T)
# head(gcor_mat)
# print(dim(gcor_mat))
# 
# foc_genes = apply(gcor_mat[, anchor_genes], 1, which.max)
# gset = gset_new_gset(sets = foc_genes, desc = "Cell cycle correlated genes")
# print(gset)
# 
# scdb_add_gset("sc_lateral", gset)
# 
# # Fine clustering and filtering of the gene set
# mcell_mat_ignore_genes(new_mat_id = "sc_lateral", mat_id = "sc", ig_genes = names(foc_genes), reverse = T)
# mcell_gset_split_by_dsmat(gset_id = "sc_lateral", mat_id = "sc_lateral", K = 5)
# gset = scdb_gset("sc_lateral")
# print(gset)
# mcell_plot_gset_cor_mats(gset_id = "sc_lateral", scmat_id = "sc_lateral")
# lateral_clusts = unique(gset@gene_set[ anchor_genes])
# mcell_gset_remove_clusts(gset_id = "sc_lateral", filt_clusts = lateral_clusts, new_id = "sc_lateral_filtered", reverse=T)
# 
# lateral_gset = scdb_gset("sc_lateral_filtered")
# print(lateral_gset)
# 
# ##################
# # Generating metacells while blacklisting lateral gene set
# # removing the lateral genes from the list of feature genes
# marker_gset.sc = scdb_gset("sc_feats")   
# marker_gset.sc = gset_new_restrict_gset(gset = marker_gset.sc, filt_gset = lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
# scdb_add_gset("sc_feats_filtered", marker_gset)
# 
# # generate the graph using the filtered list of feature genes:
# mcell_add_cgraph_from_mat_bknn(mat_id="sc", 
#                                gset_id = "sc_feats_filtered", 
#                                graph_id="sc_graph_filtered",
#                                K=100,
#                                dsamp=T)
# 
# # Resampling and generating the co-clustering graph
# mcell_coclust_from_graph_resamp(
#   coc_id="sc_coc500_filtered", 
#   graph_id="sc_graph_filtered",
#   min_mc_size=20, 
#   p_resamp=0.75, n_resamp=500)
# 
# mcell_mc_from_coclust_balanced(
#   coc_id="sc_coc500_filtered", 
#   mat_id= "sc",
#   mc_id= "sc_filtered_mc", 
#   K=30, min_mc_size=30, alpha=2)
# 
# # Removing outlier cells
# mcell_plot_outlier_heatmap(mc_id="sc_filtered_mc", mat_id = "sc", T_lfc=3)
# mcell_mc_split_filt(new_mc_id="sc_filtered_mc_f", 
#                     mc_id="sc_filtered_mc", 
#                     mat_id="sc",
#                     T_lfc=3, plot_mats=F)
# 
# # Selecting markers and coloring metacells
# mcell_gset_from_mc_markers(gset_id="sc_filtered_markers", mc_id="sc_filtered_mc_f", blacklist_gset_id="sc_lateral_filtered")
# mcell_gset_from_mc_markers(gset_id="sc_filtered_markers_lat", mc_id="sc_filtered_mc_f", filt_gset_id="sc_lateral_filtered")
# 
# mcell_mc_plot_marks(mc_id="sc_filtered_mc_f", gset_id="sc_filtered_markers", mat_id="sc", lateral_gset_id = "sc_filtered_markers_lat")
# # set default color
# mc_colorize_default("sc_filtered_mc_f")
# sc_filtered_mc_f <- scdb_mc("sc_filtered_mc_f")
# 
# # check some gene expression in metacells/clusters
# lfp = log2(sc_filtered_mc_f@mc_fp)
# tail(sort(lfp["Cr2",]))
# write.csv(lfp,"all_lfp.csv")
# plt = function(gene1, gene2, lfp, colors) 
# {
#   plot(lfp[gene1, ], lfp[gene2, ], pch=21, cex=3, bg=colors, xlab=gene1, ylab=gene2)
#   text(lfp[gene1, ], lfp[gene2, ], colnames(lfp))
#   
# }
# plt(gene1 = 'Cd19', gene2 = 'Cd3d', lfp = lfp, colors = sc_filtered_mc_f@colors)
# 
# # Projecting metacells and cells in 2D
# mcell_mc2d_force_knn(mc2d_id = "sc_filtered_mc_f", mc_id = "sc_filtered_mc_f", graph_id = "sc_graph_filtered")
# mcell_mc2d_plot(mc2d_id = "sc_filtered_mc_f")
# 
# # Visualizing the metacell confusion matrix
# mc_hc <- mcell_mc_hclust_confu(mc_id="sc_filtered_mc_f", graph_id="sc_graph_filtered")
# 
# mc_sup <- mcell_mc_hierarchy(mc_id="sc_filtered_mc_f",mc_hc=mc_hc, T_gap=0.04)
# 
# mcell_mc_plot_hierarchy(mc_id="sc_filtered_mc_f", 
#                         graph_id="sc_graph_filtered", 
#                         mc_order=mc_hc$order, 
#                         sup_mc = mc_sup, 
#                         width=2800, 
#                         height=4000, 
#                         min_nmc=2)

# ### process doublet: ###
# # Exploring and filtering the UMI matrix
# mcell_plot_umis_per_cell("db")
# 
# db_mat = scdb_mat("db")
# print(dim(db_mat@mat))
# 
# nms = c(rownames(db_mat@mat), rownames(db_mat@ignore_gmat))
# ig_genes = c(grep("^Igj", nms, v=T), 
#              grep("^Igh",nms,v=T),
#              grep("^Igk", nms, v=T), 
#              grep("^Igl", nms, v=T))
# bad_genes = unique(c(grep("ERCC", nms, v=T), grep("^mt-", nms, v=T), "Gm42418", "Gm26917", "Gm29216", ig_genes))
# 
# mcell_mat_ignore_genes(new_mat_id="db", mat_id="db", bad_genes, reverse=F)
# 
# # selecting feature genes
# mcell_add_gene_stat(gstat_id="db", mat_id="db", force=T)
# mcell_gset_filter_varmean(gset_id="db_feats", gstat_id="db", T_vm=0.08, force_new=T)
# mcell_gset_filter_cov(gset_id = "db_feats", gstat_id="db", T_tot=100, T_top3=2)
# mcell_plot_gstats(gstat_id="db", gset_id="db_feats")
# 
# # define cell cycle modules
# anchor_genes = c("Top2a", "Mki67", "Pcna", "Mcm4", "Cdk1", "Ube2c")
# 
# mcell_mat_rpt_cor_anchors(mat_id="db", gene_anchors = anchor_genes, cor_thresh = 0.1, gene_anti = c(),
#                           tab_fn = "lateral_gmods_db.txt", sz_cor_thresh = 0.2)
# 
# gcor_mat = read.table("lateral_gmods_db.txt", header=T)
# head(gcor_mat)
# print(dim(gcor_mat))
# 
# foc_genes = apply(gcor_mat[, anchor_genes], 1, which.max)
# gset = gset_new_gset(sets = foc_genes, desc = "Cell cycle correlated genes")
# print(gset)
# 
# scdb_add_gset("db_lateral", gset)
# 
# # Fine clustering and filtering of the gene set
# mcell_mat_ignore_genes(new_mat_id = "db_lateral", mat_id = "db", ig_genes = names(foc_genes), reverse = T)
# mcell_gset_split_by_dsmat(gset_id = "db_lateral", mat_id = "db_lateral", K = 5)
# gset = scdb_gset("db_lateral")
# print(gset)
# mcell_plot_gset_cor_mats(gset_id = "db_lateral", scmat_id = "db_lateral")
# lateral_clusts = unique(gset@gene_set[ anchor_genes])
# mcell_gset_remove_clusts(gset_id = "db_lateral", filt_clusts = lateral_clusts, new_id = "db_lateral_filtered", reverse=T)
# 
# lateral_gset = scdb_gset("db_lateral_filtered")
# print(lateral_gset)

# ##################
# # Generating metacells while blacklisting lateral gene set
# # removing the lateral genes from the list of feature genes
# marker_gset.db = scdb_gset("db_feats")  
# marker_gset.db = gset_new_restrict_gset(gset = marker_gset.db, filt_gset = lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
# scdb_add_gset("db_feats_filtered", marker_gset.db)
# 
# # generate the graph using the filtered list of feature genes:
# mcell_add_cgraph_from_mat_bknn(mat_id="db", 
#                                gset_id = "db_feats_filtered", 
#                                graph_id="db_graph_filtered",
#                                K=100,
#                                dsamp=T)
# 
# # Resampling and generating the co-clustering graph
# mcell_coclust_from_graph_resamp(
#   coc_id="db_coc500_filtered", 
#   graph_id="db_graph_filtered",
#   min_mc_size=20, 
#   p_resamp=0.75, n_resamp=500)
# 
# mcell_mc_from_coclust_balanced(
#   coc_id="db_coc500_filtered", 
#   mat_id= "db",
#   mc_id= "db_filtered_mc", 
#   K=30, min_mc_size=30, alpha=2)
# 
# # Removing outlier cells
# mcell_plot_outlier_heatmap(mc_id="db_filtered_mc", mat_id = "db", T_lfc=3)
# mcell_mc_split_filt(new_mc_id="db_filtered_mc_f", 
#                     mc_id="db_filtered_mc", 
#                     mat_id="db",
#                     T_lfc=3, plot_mats=F)
# 
# # Selecting markers and coloring metacells
# mcell_gset_from_mc_markers(gset_id="db_filtered_markers", mc_id="db_filtered_mc_f", blacklist_gset_id="db_lateral_filtered")
# mcell_gset_from_mc_markers(gset_id="db_filtered_markers_lat", mc_id="db_filtered_mc_f", filt_gset_id="db_lateral_filtered")
# 
# mcell_mc_plot_marks(mc_id="db_filtered_mc_f", gset_id="db_filtered_markers", mat_id="db", lateral_gset_id = "db_filtered_markers_lat")
# # set default color
# mc_colorize_default("db_filtered_mc_f")
# db_filtered_mc_f <- scdb_mc("db_filtered_mc_f")
# 
# # check some gene expression in metacells/clusters
# lfp_db = log2(db_filtered_mc_f@mc_fp)
# tail(sort(lfp_db["Cr2",]))
# tail(sort(lfp_db["Cd3d",]))
# write.csv(lfp_db,"all_lfp_db.csv")
# plt = function(gene1, gene2, lfp, colors) 
# {
#   plot(lfp[gene1, ], lfp[gene2, ], pch=21, cex=3, bg=colors, xlab=gene1, ylab=gene2)
#   text(lfp[gene1, ], lfp[gene2, ], colnames(lfp))
#   
# }
# plt(gene1 = 'Cd19', gene2 = 'Cd3d', lfp = lfp_db, colors = db_filtered_mc_f@colors)
# 
# # Projecting metacells and cells in 2D
# mcell_mc2d_force_knn(mc2d_id = "db_filtered_mc_f", mc_id = "db_filtered_mc_f", graph_id = "db_graph_filtered")
# mcell_mc2d_plot(mc2d_id = "db_filtered_mc_f")
# 
# # Visualizing the metacell confusion matrix
# mc_hc_db <- mcell_mc_hclust_confu(mc_id="db_filtered_mc_f", graph_id="db_graph_filtered")
# 
# mc_sup_db <- mcell_mc_hierarchy(mc_id="db_filtered_mc_f",mc_hc=mc_hc_db, T_gap=0.04)
# 
# mcell_mc_plot_hierarchy(mc_id="db_filtered_mc_f", 
#                         graph_id="db_graph_filtered", 
#                         mc_order=mc_hc_db$order, 
#                         sup_mc = mc_sup_db, 
#                         width=2800, 
#                         height=4000, 
#                         min_nmc=2)
# 
# 
# # test: convert to metacell object
# # for SCT
# T_cell.SCE <- as.SingleCellExperiment(T_cell.combined, assay= "SCT")
# T_cell.mat <- scm_import_sce_to_mat(T_cell.SCE, counts_slot = "logcounts")
# 
# B_cell.SCE <- as.SingleCellExperiment(B_cell.combined, assay= "SCT")
# B_cell.mat <- scm_import_sce_to_mat(B_cell.SCE, counts_slot = "logcounts")
# 
# TB_cells.mat <- scm_merge_mats(T_cell.mat, B_cell.mat)
# 
# # SCT normalized
# TB_cells.mat@mat[1:10,1:3]


######## --- PICs --- ########
# import 10X results
# PICs
sample_list_PICs <- c("Young_PICS","Aged_PICS")
Sample.list.cmb.PICS <- list()
for (s in sample_list_PICs){
  sample <- Read10X(data.dir = paste0("/data/array2/fwu/saad/",s))
  Sample.list.cmb.PICS[[s]] <- sample
}
# create seuratobject for integration
for (i in names(Sample.list.cmb.PICS)){
  Sample.list.cmb.PICS[[i]] <- CreateSeuratObject(Sample.list.cmb.PICS[[i]], project=i, min.cells = 3)
  Sample.list.cmb.PICS[[i]] <- PercentageFeatureSet(Sample.list.cmb.PICS[[i]], pattern = "^MT-", col.name = "percent.mt")
  Sample.list.cmb.PICS[[i]] <- subset(Sample.list.cmb.PICS[[i]], subset= nFeature_RNA>250 & percent.mt < 10)
  # recognize BCR and TCR genes
  pattern <- "Ig[hkl]v|Ig[kl]j|Ig[kl]c|Igh[adegm]|^Tr[abdg][vjc]"
  mat <- Sample.list.cmb.PICS[[i]][["RNA"]]$counts
  gene2keep <- grep(pattern, rownames(mat), value=T, invert = T)
  Sample.list.cmb.PICS[[i]] <- subset(Sample.list.cmb.PICS[[i]], features = gene2keep) # remove BCR and TCR genes
  Sample.list.cmb.PICS[[i]] <- SCTransform(Sample.list.cmb.PICS[[i]], vst.flavor = "v2", vars.to.regress = "percent.mt",method = "glmGamPoi",verbose = F)
}
### --- Integrate the PIC cells --- ###
PICs.features.cmb <- SelectIntegrationFeatures(object.list = Sample.list.cmb.PICS, nfeatures = 3000)
Sample.list.cmb.PICS <- PrepSCTIntegration(object.list = Sample.list.cmb.PICS, anchor.features = PICs.features.cmb)

# find the anchor
PICs.Anchors <- FindIntegrationAnchors(object.list = Sample.list.cmb.PICS, normalization.method = "SCT",
                                    anchor.features = PICs.features.cmb)
PICs.combined <- IntegrateData(anchorset = PICs.Anchors, normalization.method = "SCT")

# Dimensional reduction
PICs.combined <- RunPCA(PICs.combined, verbose = FALSE)
# ElbowPlot
ElbowPlot(PICs.combined,ndims = 50)
PICs.combined <- RunUMAP(PICs.combined, reduction = "pca", dims = 1:30)
# Cluster the cells
PICs.combined <- FindNeighbors(PICs.combined, dims = 1:30)

for (res in c(0.1, 0.2, 0.8)){
  PICs.combined <- FindClusters(PICs.combined, resolution = res)
}
colnames(PICs.combined@meta.data)
apply(PICs.combined@meta.data[,grep("snn_res",colnames(PICs.combined@meta.data))],2,table)
Idents(PICs.combined) <- "integrated_snn_res.0.2"

# Prepare object to run differential expression on SCT assay with multiple models
PICs.combined <- PrepSCTFindMarkers(PICs.combined)
# find markers for every cluster compared to all remaining cells
PICs.markers.all <- FindAllMarkers(PICs.combined, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25)

write.csv(PICs.markers.all,"PICs.markers.all.csv")
# save the UMAP
DimPlot(PICs.combined, reduction = "umap",label = T)
ggsave("umap_PICs_clusters.pdf", width=8, height=7)
DimPlot(PICs.combined,split.by = "orig.ident",label = T) # split by age
ggsave("umap_split_age_PICs_clusters_0.2.pdf", width=14, height=7)


# show the heatmap of markers
# T cells
T_top_gene_per_cluster <- unique(c(top_gene_per_cluster.T$gene, # top genes
                                 T_markers))
DoHeatmap(T_cell.combined, features = T_top_gene_per_cluster,assay = "RNA")
ggsave("T_markers.pdf",width = 12, height = 14)
# B cells
B_top_gene_per_cluster <- unique(c(top_gene_per_cluster.B$gene, # top genes
                                   B_markers))
DoHeatmap(B_cell.combined, features = B_top_gene_per_cluster,assay = "RNA")
ggsave("B_markers.pdf",width = 12, height = 14)


# find marker of T and B
T_SCT_df <- as.data.frame(T_cell.combined[["SCT"]]$data)
B_SCT_df <- as.data.frame(B_cell.combined[["SCT"]]$data)
rownames_TB <- intersect(rownames(T_SCT_df),rownames(B_SCT_df))
T_SCT_df <- T_SCT_df[rownames_TB,]
B_SCT_df <- B_SCT_df[rownames_TB,]

# Perform Wilcoxon Rank Sum Test for Each Gene
# Calculate Log2 Fold Change for Each Gene
library(parallel)
# Detect the number of CPU cores
no_cores <- 10  # Leave one core free for system processes
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, varlist = c("T_SCT_df", "B_SCT_df"))

computeStats <- function(geneIndex) {
  # Extract the gene name
  geneName <- rownames(T_SCT_df)[geneIndex]
  
  # Extract expression values for the current gene and ensure they are numeric
  expression_T <- as.numeric(T_SCT_df[geneIndex, ])
  expression_B <- as.numeric(B_SCT_df[geneIndex, ])
  
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

# Execute in parallel
results <- parLapply(cl, 1:nrow(T_SCT_df), computeStats)

stopCluster(cl)

# Convert results list to a data frame
results_df <- do.call(rbind, results)
results_df <- data.frame(gene=as.character(results_df[,1]),
                         log2FC=as.numeric(results_df[,2]),
                         p_value=as.numeric(results_df[,3]))

# Adjust p-values
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "BH")

T_DEGs <- results_df %>% filter(adjusted_p_value<0.05 & (log2FC>1 | log2FC<(-1))) %>% arrange(desc(log2FC))


# # double-check some genes
# expression_T <- as.numeric(T_SCT_df[14142, ])
# expression_B <- as.numeric(B_SCT_df[14142, ])
# test_result <- wilcox.test(expression_T, expression_B, exact = FALSE)
# p_value <- test_result$p.value
# 
# avg_expr_T <- mean(expression_T)
# avg_expr_B <- mean(expression_B)
# fold_change <- (avg_expr_T + 1e-6) / (avg_expr_B + 1e-6)  # Avoid division by zero
# log2FC <- log2(fold_change)
