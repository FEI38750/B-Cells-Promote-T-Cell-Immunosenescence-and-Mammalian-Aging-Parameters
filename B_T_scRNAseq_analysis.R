### Setup the Seurat objects ###
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(future)

setwd("~/Winer_Saad")
# import 10X results
  # WT
WT <- Read10X("/opt/home/buckcenter.org/fwu/Winer_Saad/WT/filtered_feature_bc_matrix")
  # KO
KO <- Read10X("/opt/home/buckcenter.org/fwu/Winer_Saad/KO/filtered_feature_bc_matrix")

WT<-CreateSeuratObject(WT,project="WT")
KO<-CreateSeuratObject(KO,project="μMT")
# Visualize QC metrics as a violin plot
VlnPlot(KO, features = c("nFeature_RNA"), y.max = 500)
# create seuratobject list
Sample.list <- list("WT"=WT, "μMT"=KO)
# create seuratobject list for QC
Sample.list.qc <- Sample.list

# function for Pre-process Seurat object: QC, PCA and UMAP
SCT_m<-function(l){
  l <- PercentageFeatureSet(l,pattern = "^MT-", col.name = "percent.mt")
  l <- subset(l, subset= nFeature_RNA>250 & percent.mt < 10)
  l <- SCTransform(l, vst.flavor = "v2", method = "glmGamPoi",verbose = FALSE)
  l <- RunPCA(l, npcs = 20, verbose = FALSE)
  l <- RunUMAP(l,dims = 1:20, verbose = FALSE)
  l <- FindNeighbors(l, dims = 1:20, verbose = FALSE)
  l <- FindClusters(l, verbose = FALSE)
}
# SCTransform for each dataset independently
library(parallel)
Sample.list.qc<-mclapply(Sample.list.qc,SCT_m, mc.cores=8)
# Remove doublet
library(DoubletFinder)
for (i in names(Sample.list.qc)){
  # pK Identification (no ground-truth)
  sweep.res.list_sample <- paramSweep_v3(Sample.list.qc[[i]], PCs = 1:20, sct = T)
  sweep.stats_sample <- summarizeSweep(sweep.res.list_sample, GT = FALSE)
  bcmvn_sample <- find.pK(sweep.stats_sample)
  pK<-as.numeric(as.character(bcmvn_sample$pK))[bcmvn_sample$BCmetric==max(bcmvn_sample$BCmetric)]
  ## Homotypic Doublet Proportion Estimate
  annotations <- Sample.list.qc[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(Sample.list.qc[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies
  Sample.list.qc[[i]] <- doubletFinder_v3(Sample.list.qc[[i]], PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  re.pAnn<-names(Sample.list.qc[[i]]@meta.data)[(length(Sample.list.qc[[i]]@meta.data)-1)]
  Sample.list.qc[[i]] <- doubletFinder_v3(Sample.list.qc[[i]], PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = re.pAnn, sct = T)
}

# filter samples based on results of QC & Doublet
for (i in names(Sample.list.qc)){
  sample.meta.data<-Sample.list.qc[[i]]@meta.data
  singlet<-row.names(sample.meta.data)[sample.meta.data[length(sample.meta.data)]=="Singlet"]
  Sample.list[[i]]<-Sample.list[[i]][,singlet]
}
# create SCT seuratobjects for integration
for (i in names(Sample.list)){
  Sample.list[[i]] <- SCTransform(Sample.list[[i]], vst.flavor = "v2", method = "glmGamPoi",verbose = F)
}
# Perform integration
features <- SelectIntegrationFeatures(object.list = Sample.list, nfeatures = 2000)
plan("multisession", workers = 8)
options(future.globals.maxSize= +Inf) # increase maxSize of future function
Sample.list <- PrepSCTIntegration(object.list = Sample.list, anchor.features = features)
Sample.list <- lapply(X = Sample.list, FUN = RunPCA, features = features, npcs = 20) # run PCA ahead
Anchors <- FindIntegrationAnchors(object.list = Sample.list, normalization.method = "SCT",
                                  anchor.features = features,
                                  dims = 1:20, reduction = "rpca", k.anchor = 4)
Sample.combined <- IntegrateData(anchorset = Anchors, normalization.method = "SCT",dims = 1:20)

# Dimensional reduction
Sample.combined <- RunPCA(Sample.combined, npcs = 20, verbose = FALSE)
# ElbowPlot(Sample.combined,ndims = 50)
Sample.combined <- RunUMAP(Sample.combined, reduction = "pca", dims = 1:20)
# Cluster the cells
Sample.combined <- FindNeighbors(Sample.combined, dims = 1:20)
plan("multisession", workers = 8)
Sample.combined <- FindClusters(Sample.combined, resolution = 0.5)
DimPlot(Sample.combined, reduction = "umap",label = T,repel=F)
ggsave("UMAP.pdf",width=6, height=5)
DimPlot(Sample.combined, reduction = "umap", group.by = "orig.ident")
ggsave("UMAP_conditions.pdf",width=6, height=5)


# roughly check cell proportion
Idents(Sample.combined)<-"orig.ident"
DefaultAssay(Sample.combined) <- "integrated"
Sample.combined.WT <- subset(Sample.combined, idents = "WT")
Sample.combined.KO <- subset(Sample.combined, idents = "μMT")

table(Sample.combined.WT$seurat_clusters)/8376
table(Sample.combined.KO$seurat_clusters)/7747

Idents(Sample.combined)<-"seurat_clusters"
# Prepare object to run differential expression on SCT assay with multiple models
Sample.combined <- PrepSCTFindMarkers(Sample.combined)
# find markers for every cluster compared to all remaining cells
plan("multisession", workers = 16)
markers.all <- FindAllMarkers(Sample.combined, assay = "SCT", min.pct = 0, logfc.threshold = 0.25)
write.csv(markers.all,"marker_substantial_differential_expression_per_cluster.csv")

# heatmap for marker in each group
markers.all %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Sample.combined, features = top10$gene, size=3) + 
  theme(text = element_text(size = 8.5))
ggsave("Heatmap_markers.pdf",width=12, height=12)

#DefaultAssay(Sample.combined) <- "SCT"
names(Sample.combined@meta.data)[1]<-"condition" # change name of a metadata

# ------ Cell type ratio change stack barplot ------ #
library(RColorBrewer)
sub.prop.all<-data.frame()
for (l in unique(Sample.combined$condition)){
  sub.treat<-Sample.combined@meta.data[Sample.combined$condition==l,]
  sub.prop<-data.frame(table(sub.treat$seurat_clusters)/sum(table(sub.treat$seurat_clusters)))
  sub.prop$sample<-l
  sub.prop.all<-rbind(sub.prop.all,sub.prop)
}
# cell type ratio change: all

sub.prop.all$sample <- factor(sub.prop.all$sample, levels = c("WT","μMT")) # change plot legend order 
ggplot(sub.prop.all, aes(x = Var1, y = Freq, fill = sample)) +
  geom_bar(stat = "identity",position = "dodge") +
  ggtitle("Cell type frequency") +
  scale_fill_brewer(palette="Set2",direction=-1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(x="Cell types", y="Fraction", title="Fraction of each cell type in each cluster") +
  guides(fill=guide_legend(title='Groups'))
ggsave("Cell_type_frequency.pdf",width=8,height=6, device=cairo_pdf)

library(tidyr)
# All sample groups: cell type ratio change stack barplot
## set the levels in order we want
sub.prop.all$sample<-factor(sub.prop.all$sample, 
                            levels=c("1G_unstimulated","uG_unstimulated","1G_stimulated","uG_stimulated"))
ggplot(sub.prop.all, aes(fill=Var1, y=Freq, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  labs(x="Sample", y="Proportion", title="PBMC") +
  guides(fill=guide_legend(title='Cell types'))
ggsave("stacked_barplot_all.pdf",width=8,height=6)

# comparison
sub.prop1<-sub.prop.all[sub.prop.all$sample %in% c("1G_unstimulated","uG_unstimulated"),]
## prep cumulative sums for line links
cell.stack.wide1<-sub.prop1 %>%
  pivot_wider(names_from=sample, values_from=Freq) %>% 
  arrange(by=desc(Var1)) %>% 
  replace(is.na(.), 0) %>% # remove NA
  mutate(y=cumsum(`1G_unstimulated`),
         yend=cumsum(uG_unstimulated))
## set the levels in order we want
sub.prop1$sample<-factor(sub.prop1$sample, 
                         levels=c("1G_unstimulated","uG_unstimulated"))
sub.prop1$labs<-round(sub.prop1$Freq,3) # prepare cell proportion label
sub.prop1$labs[sub.prop1$labs<0.05]<-""
ggplot(sub.prop1, aes(fill=Var1, y=Freq, x=sample)) + 
  geom_bar(stat="identity",width=0.5, col="black") +
  geom_segment(data = cell.stack.wide1, aes(x=1.25, xend=1.75, y=y, yend=yend))+
  geom_text(aes(label= labs), size = 3, position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(x="Sample", y="Proportion", title="PBMC") +
  guides(fill=guide_legend(title='Cell types'))
ggsave("stacked_barplot_1.pdf",width=7,height=6)

# violin plot to check maker genes provided by Saad
Idents(Sample.combined)
DefaultAssay(Sample.combined)<-"SCT"
VlnPlot(Sample.combined, features=c("Cd3e","Cd4","Cd8a","Cd8b1","Cd44","Pdcd1","Sell",
                                    "Ccr7","Dapl1","Lef1","Ly6c2","Gzmm","Ly6a","Tox",
                                    "Sostdc1","Tnfsf8","Klrk1","Klrb1c","Klra7","Igfbp4",
                                    "Ly6c1","Dusp10","Ccl4","Gzmk","Eomes","Tnfrsf4",
                                    "Foxp3","Areg","Lars2","Apoe","H2-Aa","Cxcr6","Cd226",
                                    "Isg15","Ifit3","Usp18","Rorc","Tmem176b","Gata3",
                                    "Trdc","Trgv2","Mki67","S100a9","Cd79a","Cd19"),
        pt.size = 0)
ggsave("Marker_check.pdf",width=30,height=30)

VlnPlot(Sample.combined, features=c("Malat1", "Mycbp", "Trbv3"),
        pt.size = 0)
ggsave("Marker_check2.pdf",width=16,height=7)


markers2plot<-c("Cd3e","Cd4","Cd8a","Cd8b1","Sell","Ccr7","Cd44",
  "Pdcd1","Ly6a","Tox","Eomes","Gzmk","Klrb1c",
  "Lars2","Malat1","Trbv3")
Idents(Sample.combined)<-"seurat_clusters"
markers_plot<-lapply(1:16,function(i)
  VlnPlot(Sample.combined,features=markers2plot[i],pt.size = 0) + NoLegend())
lapply(1:16, function(i)
  ggsave(filename=paste0("MTvsWT",i,".pdf"), plot=markers_plot[[i]],
         width = 6, height = 5,device=cairo_pdf))


# annotate clusters
Sample.combined$cell_types<-Sample.combined$seurat_clusters
Idents(Sample.combined)<-"cell_types"
cell_types <- c("0 - CD8+ Naive T Cells","1 - CD4+ EM T Cells","2 - CD8+ CM T Cells",
                "3 - CD4+ Naive T Cells","4 - Regulatory T Cells","5 - NKT Cells (1)",
                "6 - CD8+ EM T Cells","7 - Lars2hi T Cells","8 - Malat1hi T Cells",
                "9 - Trbv3+ Naive T Cells","10 - Myeloid-like CD3+ Cells","11 - NKT Cells (2)",
                "12 - IFN-Signature T Cells","13 - Mixed CD3+ Cells","14 - NKT Cells (3)",
                "15 - Ki67+ T Cells")
names(cell_types) <- levels(Sample.combined)
Sample.combined <- RenameIdents(Sample.combined, cell_types)
Sample.combined$cell_types<-Idents(Sample.combined)
DimPlot(Sample.combined, reduction = "umap", label = F)
ggsave("UMAP_annotation.pdf",width=6, height=5)

# DEGs
## overall WT vs μMT
Sample.combined$condition<-gsub("KO","μMT",Sample.combined$condition) #change to proper name
Idents(Sample.combined)<-"condition"
plan("multiprocess", workers = 16)
Overall_uMTvsWT_all <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="μMT", ident.2 ="WT",
                          min.pct=0, logfc.threshold = 0,
                          test.use = "MAST")
Overall_uMTvsWT_all<-Overall_uMTvsWT_all[order(Overall_uMTvsWT_all$avg_log2FC, decreasing = T),] # sort by log2FC
write.csv(Overall_uMTvsWT_all, "Overall_DEG_uMTvsWT_all.csv")

Overall_uMTvsWT <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="μMT", ident.2 ="WT",
                                   min.pct=0.05, logfc.threshold = 0.25,
                                   test.use = "MAST") # 347 DEGs
Overall_uMTvsWT<-Overall_uMTvsWT[order(Overall_uMTvsWT$avg_log2FC, decreasing = T),] # sort by log2FC
write.csv(Overall_uMTvsWT, "Overall_DEG_uMTvsWT.csv")

Overall_uMTvsWT_10K <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="μMT", ident.2 ="WT",
                                     min.pct=0.01, logfc.threshold = 0,
                                     test.use = "MAST") 
Overall_uMTvsWT_10K<-Overall_uMTvsWT_10K[order(Overall_uMTvsWT_10K$avg_log2FC, decreasing = T),] # sort by log2FC
write.csv(Overall_uMTvsWT_10K, "Overall_DEG_uMTvsWT_10K.csv")




# ------ find DEG for all cell types MT vs WT ------ #
Sample.combined$celltype.condition <- paste0(Sample.combined$cell_types,"_",Sample.combined$condition)
Idents(Sample.combined)<-"celltype.condition"
plan("multiprocess", workers = 16)
for (d in unique(Sample.combined$cell_types)){
  DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_μMT"), ident.2 =paste0(d,"_WT"),
                        min.pct=0.01, logfc.threshold = 0,
                        test.use = "MAST")
  write.csv(DEG.uG, paste0(d,"_μMTvsWT.csv"))
}

# ------ GSEA ------ #
library(clusterProfiler)
library(enrichplot)
library(ggnewscale)
library(DO.db)
library(msigdbr)
library(stringr)
library(dplyr)

## prerank genelist for all cell types
### H and C2.CP data base
msig.H <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
msig.C2.CP <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP") %>% 
  dplyr::select(gs_name, gene_symbol)
msig.HCP<-rbind(msig.H,msig.C2.CP)
setwd("~/Winer_Saad/Cell_types")
for (d in c("Global","0 - CD8+ Naive T Cells","1 - CD4+ EM T Cells","2 - CD8+ CM T Cells",
            "3 - CD4+ Naive T Cells","4 - Regulatory T Cells","5 - NKT Cells (1)",
            "6 - CD8+ EM T Cells","7 - Lars2hi T Cells")){
  DEG.uG.rnk <- read.csv(paste0(d,"_μMTvsWT.csv"))
  DEG.uG.rnk<-DEG.uG.rnk[order(DEG.uG.rnk$avg_log2FC,decreasing = T),]
  DEG.uG.rnk <- DEG.uG.rnk[DEG.uG.rnk$p_val <= 0.05,]
  DEG.uG.rnk4gsea<-DEG.uG.rnk[['avg_log2FC']]
  names(DEG.uG.rnk4gsea)<-DEG.uG.rnk[["X"]]
  genelist_rnk<-DEG.uG.rnk4gsea
  try(gsearun <- GSEA(genelist_rnk, TERM2GENE = msig.HCP, pvalueCutoff = 1))
  write.csv(gsearun@result,paste0(d,"_GSEA_μMTvsWT.csv"))
}

#### draw dot plot for GSEA results
cell4dot<-rev(c("Global","0 - CD8+ Naive T Cells","1 - CD4+ EM T Cells","2 - CD8+ CM T Cells",
            "3 - CD4+ Naive T Cells","4 - Regulatory T Cells","5 - NKT Cells (1)",
            "6 - CD8+ EM T Cells","7 - Lars2hi T Cells"))
ls.dot<-data.frame()
for (d in cell4dot){
  pathway <- read.csv(paste0(d,"_GSEA_μMTvsWT.csv"))
  if (nrow(pathway) != 0){
    pathway<-pathway[order(pathway$NES,decreasing = T),]
    pathway<-pathway[pathway$pvalue <= 0.05,]
    pathway$Cell_type<-d
    ls.dot<-rbind(ls.dot,pathway)
  }
}
pw<-ls.dot$ID %>% unique()
ls.dot$Cell_type<-factor(ls.dot$Cell_type, levels = unique(ls.dot$Cell_type)) # order the x axis
ls.dot<-ls.dot[order(ls.dot$NES,decreasing=F),]
ls.dot$ID<-factor(ls.dot$ID, levels = unique(ls.dot$ID))
ls.dot$direction <- with(ls.dot, ifelse(NES<0,"WT","μMT"))

ls.dot %>% filter(ID %in% pw) %>% 
  ggplot(aes(x=ID, y = Cell_type, color = direction, size = -log10(p.adjust))) + 
  geom_point() +
  scale_y_discrete(labels=function(y) str_wrap(y, width=70)) +
  ylab('Within-cluster μMT vs WT comparison') +
  xlab('GSEA pathway') +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle=45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=10))+
  theme(axis.line  = element_blank()) +
  theme(axis.ticks = element_blank()) +
  scale_color_brewer(palette="Set2",direction=-1,name = 'Pathway upregulated in')
  
ggsave("GSEA_μMTvsWT.pdf", width = 14, height = 6,limitsize = FALSE,device=cairo_pdf)


# iAge
library(scales)
library(ggpubr)

# read GE_iAge gene and coefficients
GE_iAge <- read.csv("GE_iAge_Exclude_age_mice.csv")
# diveide into to groups: 1G and uG
WT4iAge<-Sample.combined@assays[["SCT"]]@data[intersect(GE_iAge$Genes,row.names(Sample.combined@assays[["SCT"]]@data)),Sample.combined$condition=="WT"]
MT4iAge<-Sample.combined@assays[["SCT"]]@data[intersect(GE_iAge$Genes,row.names(Sample.combined@assays[["SCT"]]@data)),Sample.combined$condition=="μMT"]

# calculate iAge index for each cells
sumGene<-0
for (g in intersect(GE_iAge$Genes,row.names(Sample.combined@assays[["SCT"]]@data))){
  resGene<-rescale(WT4iAge[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
  sumGene<-sumGene+resGene
}
iAge_index_WT<-sumGene+10
summary(iAge_index_WT)

sumGene<-0
for (g in intersect(GE_iAge$Genes,row.names(Sample.combined@assays[["SCT"]]@data))){
  resGene<-rescale(MT4iAge[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
  sumGene<-sumGene+resGene
}
iAge_index_MT<-sumGene+10
summary(iAge_index_MT)

iAge.data<-rbind(data.frame(iAge_index=iAge_index_WT,group=rep("WT",length(iAge_index_WT))),
                 data.frame(iAge_index=iAge_index_MT,group=rep("μMT",length(iAge_index_MT))))


ggboxplot(iAge.data, x = "group", y = "iAge_index",
          color = "group", palette = "jco",outlier.shape = NA) +
  stat_compare_means (method = "wilcox.test", label = "p.signif") +
  coord_cartesian(ylim = c(0,18)) +
  scale_color_brewer(palette="Set2",direction=-1)
ggsave("GE_iAge_Mice_Saad.pdf", width = 6, height = 5,device=cairo_pdf)
write.csv(iAge.data,"GE_iAge_Mice_Saad.csv")
t.test.result<-t.test(iAge_index_MT,iAge_index_WT)
hist(iAge_index_WT); hist(iAge_index_MT)

# --- iAge for cell types --- #
# calculate iAge index for each cells
iAge.data.cell<- data.frame(matrix(ncol = 3, nrow = 0))
colnames(iAge.data.cell) <- c("iAge_index","seurat_clusters","group")
for (c in unique(Sample.combined$seurat_clusters)){
  cell4iage<-dplyr::filter(Sample.combined@meta.data,seurat_clusters==c & condition == "WT")
  cell4iage<-row.names(cell4iage)
  WT4iAge_cell<-WT4iAge[,cell4iage]
  sumGene<-0
  for (g in intersect(GE_iAge$Genes,row.names(Sample.combined@assays[["SCT"]]@data))){
    resGene<-rescale(WT4iAge_cell[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
    sumGene<-sumGene+resGene
  }
  iAge_index_WT<-sumGene+10
  iAge.data.cell<-rbind(iAge.data.cell,
                        data.frame(iAge_index=iAge_index_WT,seurat_clusters=rep(c,length(iAge_index_WT)),group=rep("WT",length(iAge_index_WT))))
}

iAge.data.cell2<- data.frame(matrix(ncol = 3, nrow = 0))
colnames(iAge.data.cell2) <- c("iAge_index","cell_type","group")
for (c in unique(Sample.combined$seurat_clusters)){
  cell4iage<-dplyr::filter(Sample.combined@meta.data,seurat_clusters==c & condition == "μMT")
  cell4iage<-row.names(cell4iage)
  MT4iAge_cell<-MT4iAge[,cell4iage]
  sumGene<-0
  for (g in intersect(GE_iAge$Genes,row.names(Sample.combined@assays[["SCT"]]@data))){
    resGene<-rescale(MT4iAge_cell[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
    sumGene<-sumGene+resGene
  }
  iAge_index_MT<-sumGene+10
  iAge.data.cell2<-rbind(iAge.data.cell2,
                         data.frame(iAge_index=iAge_index_MT,seurat_clusters=rep(c,length(iAge_index_MT)),group=rep("μMT",length(iAge_index_MT))))
}
iAge.data.cell<-rbind(iAge.data.cell,iAge.data.cell2)
write.csv(iAge.data.cell,"GE_iAge_mice_Saas_celltypes.csv")

iAge.data.cell$seurat_clusters <- factor(iAge.data.cell$seurat_clusters, levels = c("0","1","2","3",
                                                                                    "4","5","6","7",
                                                                                    "8","9","10","11",
                                                                                    "12","13","14","15")) # change plot legend order 
iAge.data.cell$group<-factor(iAge.data.cell$group, levels=c("WT","μMT"))
ggboxplot(iAge.data.cell, x = "seurat_clusters", y = "iAge_index",
          color = "group", palette = "jco",outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_brewer(palette="Set2",direction=-1)
ggsave("GEiage_mice_Saad_celltypes.pdf", width = 8, height = 6, device=cairo_pdf)


# calculate FC and P-value for each comparison in each cell type
iAge.data.cell<- data.frame(matrix(ncol = 3, nrow = 0))
t.test.result.df<- data.frame(matrix(ncol = 3, nrow = 0))
colnames(iAge.data.cell) <- c("iAge_index","seurat_clusters","group")
for (c in unique(Sample.combined$seurat_clusters)){
  cell4iage<-dplyr::filter(Sample.combined@meta.data,seurat_clusters==c & condition == "WT")
  cell4iage<-row.names(cell4iage)
  WT4iAge_cell<-WT4iAge[,cell4iage]
  sumGene<-0
  for (g in intersect(GE_iAge$Genes,row.names(Sample.combined@assays[["SCT"]]@data))){
    resGene<-rescale(WT4iAge_cell[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
    sumGene<-sumGene+resGene
  }
  iAge_index_WT<-sumGene+10
  
  cell4iage<-dplyr::filter(Sample.combined@meta.data,seurat_clusters==c & condition == "μMT")
  cell4iage<-row.names(cell4iage)
  MT4iAge_cell<-MT4iAge[,cell4iage]
  sumGene<-0
  for (g in intersect(GE_iAge$Genes,row.names(Sample.combined@assays[["SCT"]]@data))){
    resGene<-rescale(MT4iAge_cell[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
    sumGene<-sumGene+resGene
  }
  iAge_index_MT<-sumGene+10
  # non-parametric t test
  t.test.result <- wilcox.test(iAge_index_MT,iAge_index_WT)
  Log2FC<-log2(median(iAge_index_MT)/median(iAge_index_WT))
  t.test.result.df<-rbind(t.test.result.df,data.frame(cell_type=c,Log2FC=Log2FC,p.value=t.test.result$p.value))
}
write.csv(t.test.result.df,"GE_iAge_MTvsWT_Saad_mice_celltypes.csv")


## SenMayo
geneset_features <- list(SenMayo_mice)
### scoring
SC_module <- AddModuleScore(
  object = Sample.combined,
  features = geneset_features,
  name = 'SenMayo')
ggplot(SC_module@meta.data, aes(x=SenMayo1, fill=condition)) +
  geom_density(alpha=0.4) + theme_bw() + xlim(-0.15,0.25)
ggsave("SenMayo.pdf", width = 6, height = 5,device=cairo_pdf)
ggplot(SC_module@meta.data) + 
  ggridges::geom_density_ridges(aes(x=SenMayo1, y= cell_types, fill=condition), alpha=0.5) + theme_bw()
ggsave("SenMayo_celltypes.pdf", width = 6, height = 6,device=cairo_pdf)
t.test(filter(SC_module@meta.data,condition=="WT")$SenMayo1,
       filter(SC_module@meta.data,condition=="μMT")$SenMayo1)
# make a dataframe - dat for violin and boxplot
c1<-c(rep("WT",length(filter(SC_module@meta.data,condition=="WT")$SenMayo1)),
      rep("μMT",length(filter(SC_module@meta.data,condition=="μMT")$SenMayo1)))
c2<-c(filter(SC_module@meta.data,condition=="WT")$SenMayo1,filter(SC_module@meta.data,condition=="μMT")$SenMayo1)
dat<-data.frame(Group=c1,Score=c2)
# violin plot
ggplot(dat,aes(Group,Score)) + 
  geom_violin(aes(fill = Group)) +
  theme_bw() + scale_fill_brewer(palette = 'Set2',direction = -1) +
  labs(title="SenMayo") +
  geom_boxplot(width=0.1,outlier.shape = NA) +
  ggpubr::stat_compare_means(method = "t.test", 
                             comparisons=list(c("WT","μMT")), 
                             label = 'p.signif')
ggsave("MTvsWT_SenMayo_violin.pdf", width = 5, height = 7,device=cairo_pdf)
# boxplot
ggboxplot(dat, x="Group",y="Score",
          color = "Group", palette = "jco",outlier.shape = NA) + coord_cartesian(ylim = c(-0.12,0.15)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",label.y=0.11,
                     comparisons=list(c("1G","uG"))) +
  scale_color_brewer(palette="Set2",direction=-1)
ggsave("MTvsWT_SenMayo_box.pdf", width = 6, height = 5,device=cairo_pdf)


# --- VDJ analysis --- #
# load the TCR sequencing data
tcrDir = "~/Winer_Saad/winer_mouse_tcr/data/2022_winer_mouse_tcr/KO/TCR/"
clonotypeTcrFile = str_c(tcrDir, "clonotypes.csv")
cellTcrFile = str_c(tcrDir, "filtered_contig_annotations.csv")

clonotypes = read.csv(clonotypeTcrFile, check.names=FALSE, as.is=T)
cellTcrs = read.csv(cellTcrFile, check.names=FALSE, as.is=T)

toUniqueString = function(v) {return(str_c(sort(unique(v)), collapse=","))}
perCell = cellTcrs %>% group_by(barcode) %>% summarise_each(list(toUniqueString))

cellTcrDat = perCell[,c("barcode", "raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id")]

tcrDat = merge(clonotypes[,c("clonotype_id", "frequency", "proportion")], cellTcrDat, by.x="clonotype_id", by.y="raw_clonotype_id")
tcrDat$cell_id = str_c(tcrDat$barcode,"_2")

length(intersect(tcrDat$cell_id, rownames(Sample.combined@meta.data)))

tcrDatKO = tcrDat


tcrDir = "~/Winer_Saad/winer_mouse_tcr/data/2022_winer_mouse_tcr/WT/TCR/"
clonotypeTcrFile = str_c(tcrDir, "clonotypes.csv")
cellTcrFile = str_c(tcrDir, "filtered_contig_annotations.csv")

clonotypes = read.csv(clonotypeTcrFile, check.names=FALSE, as.is=T)
cellTcrs = read.csv(cellTcrFile, check.names=FALSE, as.is=T)

toUniqueString = function(v) {return(str_c(sort(unique(v)), collapse=","))}
perCell = cellTcrs %>% group_by(barcode) %>% summarise_each(list(toUniqueString))

cellTcrDat = perCell[,c("barcode", "raw_clonotype_id", "raw_consensus_id", "exact_subclonotype_id")]

tcrDat = merge(clonotypes[,c("clonotype_id", "frequency", "proportion")], cellTcrDat, by.x="clonotype_id", by.y="raw_clonotype_id")
tcrDat$cell_id = str_c(tcrDat$barcode,"_1")

length(intersect(tcrDat$cell_id, rownames(Sample.combined@meta.data)))

tcrDatWT = tcrDat

tcrDatAll = rbind(tcrDatKO, tcrDatWT)

Sample.combined$cell_id <- rownames(Sample.combined@meta.data)
meta = Sample.combined@meta.data
length(intersect(meta$cell_id, tcrDatAll$cell_id))
length(setdiff(meta$cell_id, tcrDatAll$cell_id))
length(setdiff(tcrDatAll$cell_id, meta$cell_id))

metaTcr = merge(meta, tcrDatAll, by="cell_id", all.x=TRUE, sort=FALSE)
rownames(metaTcr) = metaTcr$cell_id

setdiff(rownames(meta), rownames(metaTcr))

Sample.combined@meta.data = metaTcr

meta = Sample.combined@meta.data
coordUmap = as.data.frame(Sample.combined[["umap"]]@cell.embeddings)

meta<-meta[order(row.names(meta)), ]
coordUmap<-coordUmap[order(row.names(coordUmap)), ]
table(rownames(meta) == rownames(coordUmap))

dat = cbind.data.frame(meta, coordUmap, stringsAsFactors=F)
dat = dat %>% filter(!is.na(clonotype_id))

cloneDat = dat %>% 
  #            filter(frequency > 1) %>% 
  group_by(seurat_clusters, condition, clonotype_id) %>% 
  summarize(clone_count=n(),
            clone_x=mean(UMAP_1, na.rm=T),
            clone_y=mean(UMAP_2, na.rm=T))

cloneDatUmap = cloneDat %>%
  filter(clone_count > 2) %>%
  arrange(clone_count)


tcrColors = brewer.pal(n=3, name="Set2")[2:1]
legendLabels = c("WT", expression(mu~"MT"))

pl = ggplot(dat, aes(x=UMAP_1, y=UMAP_2))
pl = pl + geom_point(color="darkgrey")
pl = pl + geom_point(data=cloneDatUmap, aes(x=clone_x, clone_y, size=clone_count, fill=condition), pch=21)
pl = pl + scale_size_area(trans=scales::identity_trans(), max_size=10, n.breaks=10)
pl = pl + ggplotThemeBlock()
pl = pl + scale_fill_manual(values=tcrColors, labels=legendLabels)
#pl = pl + scale_fill_brewer(palette="Set2")
pl = pl + labs(fill="", size="Clonality")
pl

ggsave("umap_with_tcr_clone_size_overlay.pdf", pl, h=6, w=10)

plByCondition = pl + facet_wrap(~ condition, ncol=2)
plByCondition

ggsave("umap_with_tcr_clone_size_overlay-by_condition.pdf", plByCondition, h=6, w=20, device=cairo_pdf)

plByCluster = pl + facet_wrap(~ seurat_clusters, ncol=2)

ggsave("umap_with_tcr_clone_size_overlay-by_cluster.pdf", plByCluster, h=40, w=18, limitsize=F,device=cairo_pdf)

cloneDatCum = cloneDat %>% arrange(desc(clone_count)) %>% group_by(condition) %>% mutate(cumulative_count=cumsum(clone_count),
                                                                                         cumulative_index=row_number()
)

pl = ggplot(data=cloneDatCum, aes(x=cumulative_index, y=cumulative_count, color=condition))
pl = pl + geom_point()
pl = pl + geom_line(aes(group=condition))
pl = pl + ggplotThemeBlock()
pl = pl + scale_color_manual(name="", values=tcrColors, labels=legendLabels)
pl = pl + xlab("Unique TCR sequence count")
pl = pl + ylab("Unique cell count")
#pl = pl + scale_fill_brewer(palette="Set2")
#pl = pl + scale_x_log10()
#pl = pl + scale_y_log10()
pl

cloneDatTab = cloneDat %>% group_by(condition, clone_count) %>% summarize(clonality_count=n())

# Correct log-log plot for different # of WT/KO cells
cloneDatTotals = cloneDatTab %>% group_by(condition) %>% summarize(total_cells = sum(clone_count*clonality_count))

koRatio = cloneDatTotals %>% filter(condition == "WT") %>% dplyr::select(total_cells) / cloneDatTotals %>% filter(condition == "μMT") %>% dplyr::select(total_cells)
koRatio = as.numeric(koRatio)

cloneDatLogLog = cloneDatTab %>% mutate(clonality_count_adjusted = round(ifelse(condition == "μMT", clonality_count * koRatio, clonality_count)))
cloneDatLogLog$condition = factor(cloneDatLogLog$condition, levels=c("μMT", "WT"))


pl = ggplot(data=cloneDatLogLog, aes(x=clone_count, y=clonality_count_adjusted, fill=condition))
pl = pl + geom_point(size=6, pch=21, alpha=0.9)
pl = pl + scale_x_log10(breaks=c(1, 2, 3, 4, 5, 6, 7, 10, 15, 25, 40, 55, 110),
                        labels=c("1\n(Unique Clones)", as.character(c(2, 3, 4, 5, 6, 7, 10, 15, 25, 40, 55, 110))))
pl = pl + scale_y_log10(breaks=c(1,3,5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000))
pl = pl + geom_vline(xintercept=1.5, linetype="dashed")
#pl = pl + annotate(geom="text", x=3, y=5000, label="Evidence of\nClonal Expansion ->", size=10)
pl = pl + ggplotThemeBlock()
pl = pl + scale_fill_brewer(palette="Set2", labels=c(expression(mu~"MT"), "WT"))
pl = pl + labs(fill="")
pl = pl + xlab("Clonality\n(# cells with specific TCR sequence)")
pl = pl + ylab("Number of distinct TCR sequences")
pl

ggsave("clonality_size_frequency_log_log_dot_plot.pdf", pl, h=8, w=12, useDingbats=F)

# Test differences between the 2 distributions of TCR clone counts using 2 different tests (Wilcoxon Ranksum, and Kolmogorov Smirnov test for differences in distributions)
ks.test(cloneDatCum$clone_count[cloneDatCum$condition == "WT"], cloneDatCum$clone_count[cloneDatCum$condition == "μMT"])

wilcox.test(cloneDatCum$clone_count[cloneDatCum$condition == "WT"], cloneDatCum$clone_count[cloneDatCum$condition == "μMT"])

# Trajectory analysis
# ------ Monocle3 ------ #
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

DefaultAssay(Sample.combined)<-"integrated"
cds <- as.cell_data_set(Sample.combined)
cds <- cluster_cells(cds)

cds <- learn_graph(cds)

cds <- order_cells(cds)
plot_cells(cds)
ggsave("Trej_all.pdf", width =7, height = 6)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE)
ggsave("Trej_all_pseudotime.pdf", width =7, height = 6)

# trajectory for WT
Idents(Sample.combined)<-"condition"
cds <- as.cell_data_set(subset(Sample.combined, idents = "WT"))
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds)
ggsave("Trej_WT.pdf", width =7, height = 6)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE)
ggsave("Trej_WT_pseudotime.pdf", width =7, height = 6)

# trajectory for MT
Idents(Sample.combined)<-"condition"
cds <- as.cell_data_set(subset(Sample.combined, idents = "μMT"))
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds)
ggsave("Trej_MT.pdf", width =7, height = 6)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE)
ggsave("Trej_MT_pseudotime.pdf", width =7, height = 6)

# --- Volcano plot --- #
library(ggrepel)
volcanoLabelsOverall = read.delim("~/Winer_Saad/volcano_plot_genes_global_categories.txt", header=T, as.is=T, check.names=F)

markersKOvsWT <- read.csv("~/Winer_Saad/Cell_types/Global_μMTvsWT.csv",as.is=T, check.names=F)
colnames(markersKOvsWT)[1]<-"symbol"
minP = min(markersKOvsWT %>% dplyr::select(p_val) %>% filter(p_val > 0))
plotDatGlobal = markersKOvsWT %>% mutate(p_val = ifelse(p_val == 0, minP, p_val))
plotDatGlobal = left_join(plotDatGlobal, volcanoLabelsOverall, by="symbol")
plotDatGlobal = plotDatGlobal %>% mutate(highlight=ifelse(is.na(category), "small", "large"))
plotDatGlobal$highlight = factor(plotDatGlobal$highlight, levels=c("small", "large"))
plotDatHighlight = plotDatGlobal %>% filter(!is.na(category))

categories = sort(unique(plotDatHighlight$category))
categories = c(setdiff(categories, "NicheNet"), "NicheNet")

plotDatHighlight$category = factor(plotDatHighlight$category, levels=categories)

paletteColors = brewer.pal(7, "Accent")[-4]
paletteColors = c(paletteColors[c(4:6)], "orange1", "cadetblue", paletteColors[1])
names(paletteColors) = categories

# logFcBreaks = c(-2,-log2(3), -1, -log2(1.5), 0, log2(1.5), 1)
# logFcLabels = c(-4, -2, -1.5, "same", 1.5, 2)

upInKO = "Up~mu~MT"
upInWT = "Up WT"

pl = ggplot(data=plotDatGlobal, aes(x=avg_log2FC, y=-log10(p_val)))
pl = pl + geom_point(color="lightgrey", size=2)
pl = pl + geom_point(data=plotDatHighlight, aes(color=category), size=6)
pl = pl + geom_text_repel(data=plotDatHighlight, aes(label=symbol, color=category), fontface="bold", force=70, max.overlaps=200, size=8, segment.linetype=2)
pl = pl + scale_y_continuous(limits=c(0, 450), breaks=c(0,50,100,150,200,250,300, 350))
#pl = pl + scale_x_continuous(limits=c(-2.5, 1.3), breaks=logFcBreaks, labels=logFcLabels)
pl = pl + geom_vline(xintercept=0, linetype="dotted")
pl = pl + annotate(geom="text", label=upInWT, x=-0.5, y=420, size=7)
pl = pl + annotate(geom="text", label=upInKO, parse=TRUE, x=0.5, y=420, size=7)
pl = pl + scale_color_manual(name="Category", values=paletteColors)
pl = pl + ggplotThemeBlock()
pl = pl + xlab(expression("log2 Fold Change ("~mu~"MT / WT)"))
pl = pl + ylab("-log10 P-Value")
pl
ggsave("overall_muMT_vs_WT_differential_genes_by_category.pdf", pl, h=9, w=18, useDingbats=FALSE)

plStub = ggplot(data=plotDatGlobal, aes(x=avg_log2FC, y=-log10(p_val)))
plStub = plStub + geom_point(color="lightgrey", size=2)
plStub = plStub + scale_y_continuous(limits=c(0, 450), breaks=c(0,50,100,150,200,250,300, 350))
#plStub = plStub + scale_x_continuous(limits=c(-2.5, 1.3), breaks=logFcBreaks, labels=logFcLabels)
plStub = plStub + geom_vline(xintercept=0, linetype="dotted")
plStub = plStub + annotate(geom="text", label=upInWT, x=-0.5, y=420, size=9)
plStub = plStub + annotate(geom="text", label=upInKO, parse=TRUE, x=0.5, y=420, size=9)
plStub = plStub + ggplotThemeBlock()
plStub = plStub + xlab(expression("log2 Fold Change ("~mu~"MT / WT)"))
plStub = plStub + ylab("-log10 P-Value")
#plStub

uniqueCategories = unique(plotDatGlobal$category)
uniqueCategories = uniqueCategories[!is.na(uniqueCategories)]

for(currCategory in uniqueCategories) {
  
  plotDatHighlightCategory = plotDatHighlight %>% filter(category==currCategory)
  categoryColor = paletteColors[currCategory]
  
  plCat = plStub + geom_point(data=plotDatHighlightCategory, aes(color=category), size=6)
  plCat = plCat + geom_text_repel(data=plotDatHighlightCategory, aes(label=symbol, color=category), fontface="bold", force=70, max.overlaps=200, size=8, segment.linetype=2)
  plCat = plCat + scale_color_manual(name="Category", values=categoryColor)  
  
  catFormatted = str_replace_all(currCategory, "/", "_")
  
  outfileV = str_c("overall_muMT_vs_WT_differential_genes_by_category-", catFormatted ,".pdf")
  ggsave(outfileV, plCat, h=9, w=18, useDingbats=FALSE)
}


# Genes to highlight for each cluster-specific KO vs WT analysis volcano plot
##volcanoLabelsPerCluster
# DE files per cluster
##koVsWtDE
volcanoLabelsPerCluster = read.delim("~/Winer_Saad/volcano_plot_genes_per_cluster.txt", header=T, as.is=T, check.names=F)
koVsWtFiles = dir("~/Winer_Saad/Volcano", full.names=T)

koVsWtDE = list()
for(currFile in koVsWtFiles) {
  currCluster = str_replace(currFile, ".*Volcano/([0-9]).*", "\\1")
  print(currCluster)
  currDE = read.csv(currFile, header=T, as.is=T, check.names=F)
  koVsWtDE[[currCluster]] = currDE
}

for(currCluster in names(koVsWtDE)) {
  
  currMarkers = koVsWtDE[[currCluster]]
  names(currMarkers)[1]<-"symbol"
  currHighlights = volcanoLabelsPerCluster %>% filter(cluster == currCluster) %>% dplyr::select(symbol, direction)
  currentMarkerHighlights = inner_join(currMarkers, currHighlights)
  
  minP = min(currMarkers %>% dplyr::select(p_val) %>% filter(p_val > 0))
  currPlotDat = currMarkers %>% mutate(p_val = ifelse(p_val == 0, minP, p_val))
  currPlotDat = left_join(currPlotDat, volcanoLabelsOverall, by="symbol")
  
  currentMarkerHighlights$direction = factor(currentMarkerHighlights$direction, levels=c("up_KO", "down_KO"))
  
  paletteColorsDuo = c("up_KO"=tcrColors[2], "down_KO"=tcrColors[1])
  
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
  pl = pl + scale_color_manual(name="", values=paletteColorsDuo, labels=c(expression("Up in "~mu~"MT"), "Up in WT"))
  pl = pl + ggplotThemeBlock()
  pl = pl + xlab(expression("log2 Fold Change ("~mu~"MT / WT)"))
  pl = pl + ylab("-log10 P-Value")  
  pl
  
  volcanoPlotFilePDF = str_c("volcano_plot_muMT_vs_WT_differential_genes-cluster_", currCluster  ,".pdf")
  ggsave(volcanoPlotFilePDF, pl, h=7, w=10)
}



