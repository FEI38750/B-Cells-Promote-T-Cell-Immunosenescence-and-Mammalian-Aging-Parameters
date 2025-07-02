# Generate dataframe of young and aged artificial doublets (AD) and true PICs for DEG analysis (not for training; after training and annotation steps)
# To integrate with Seurat pipeline
# run this sciript after PICseq_10XSingleCell.R

# Step 1: extract the raw counts from singlet A and B and true PICs
rawCounts_T <- T_cell.combined[["RNA"]]$counts
rawCounts_B <- B_cell.combined[["RNA"]]$counts
rawCounts_PICs <- PICs.combined[["RNA"]]$counts

# Step 2: select cell IDs of singlet A and B and true PICs by conditions
# select cell IDs of young T cells
t_ct_Young.overall <- T_cell.combined@meta.data %>% filter(orig.ident=="Young_T_Cells") %>% row.names() # young
# select cell IDs of aged T cells
t_ct_Aged.overall <- T_cell.combined@meta.data %>% filter(orig.ident=="Aged_T_Cells") %>% row.names() # aged
# select cell IDs of young B cells
b_ct_Young.overall <- B_cell.combined@meta.data %>% filter(orig.ident=="Young_B_Cells") %>% row.names() # young
# select cell IDs of aged B cells
b_ct_Aged.overall <- B_cell.combined@meta.data %>% filter(orig.ident=="Aged_B_Cells") %>% row.names() # aged
# select cell IDs of young PICs cells
PICs_ct_Young.overall <- PICs.combined@meta.data %>% filter(orig.ident=="Young_PICS") %>% row.names() # young
# select cell IDs of aged PICs cells
PICs_ct_Aged.overall <- PICs.combined@meta.data %>% filter(orig.ident=="Aged_PICS") %>% row.names() # aged

# Make sure the same order of feature names in singlet A and B for AD generation
comfeatures <- intersect(rownames(rawCounts_T),rownames(rawCounts_B))
rawCounts_T <- rawCounts_T[comfeatures,]
rawCounts_B <- rawCounts_B[comfeatures,]

# Step 3: subset raw count matrix for conditions
# subset for young
rawCounts_T.young <- rawCounts_T[,t_ct_Young.overall]
rawCounts_B.young <- rawCounts_B[,b_ct_Young.overall]
rawCounts_PICs.young <- rawCounts_PICs[,PICs_ct_Young.overall]
# subset for aged
rawCounts_T.aged <- rawCounts_T[,t_ct_Aged.overall]
rawCounts_B.aged <- rawCounts_B[,b_ct_Aged.overall]
rawCounts_PICs.aged <- rawCounts_PICs[,PICs_ct_Aged.overall]

# Step 4: Generate AD matrix
# function to generate AD matrix
AD_generation <- function(rawCounts_T,rawCounts_B,major_T,major_B,T_cell.combined_major,B_cell.combined_major){
  randSize <- 5000
  MajorAD.df.ls <- list()
  for (t in major_T){
    print(t)
    #T_cell.combined_major <- T_cell.combined$major
    t_ct <- names(T_cell.combined_major[T_cell.combined_major==t]) # select cell IDs of a T cell type
    for (b in major_B){
      print(b)
      AD_label <- paste0("T",t,"_","B",b)
      print(paste0("generate artificial doublets (AD): ",AD_label))
      #B_cell.combined_major <- B_cell.combined$major
      b_ct <- names(B_cell.combined_major[B_cell.combined_major==b]) # select cell IDs of a B cell type
      # random sampling T and B cells
      t_rand <- sample(t_ct,randSize,replace = T)
      b_rand <- sample(b_ct,randSize,replace = T)
      # get gene expression matrix
      t_data <- rawCounts_T[,t_rand]
      b_data <- rawCounts_B[,b_rand]
      # raw counts of artificial doublets (AD)
      data_AD <- t_data + b_data
      colnames(data_AD) <- paste0(AD_label,"_",seq(1,ncol(data_AD)))
      
      MajorAD.df.ls[[AD_label]] <- data_AD
    }
  }
  MajorAD.df <- purrr::reduce(MajorAD.df.ls, cbind)
  return(MajorAD.df)
}

# generate AD count matrix for singlet A and B of major cell types in each condtions
# Young
T_cell.combined_major <- T_cell.combined$major[t_ct_Young.overall] # select cell id with cell type code
B_cell.combined_major <- B_cell.combined$major[b_ct_Young.overall]
MajorAD.young <- AD_generation(rawCounts_T=rawCounts_T.young,
                               rawCounts_B=rawCounts_B.young,
                               major_T=major_T, major_B=major_B,
                               T_cell.combined_major=T_cell.combined_major,
                               B_cell.combined_major=B_cell.combined_major)

# Aged
T_cell.combined_major <- T_cell.combined$major[t_ct_Aged.overall]
B_cell.combined_major <- B_cell.combined$major[b_ct_Aged.overall]
MajorAD.aged <- AD_generation(rawCounts_T=rawCounts_T.aged,
                               rawCounts_B=rawCounts_B.aged,
                               major_T=major_T, major_B=major_B,
                               T_cell.combined_major=T_cell.combined_major,
                               B_cell.combined_major=B_cell.combined_major)

# generate true PICs count matrix of major cell types in each conditions
# subset true PICs by conditions
# subset true PICs for Young and Aged in Major cell types
# select cell IDs in major cell types
PICs_Major_CellID <- PICs_anno$Cell_ID[PICs_anno$pred != "unknown"]
PICs_Major_CellID.young <- intersect(colnames(rawCounts_PICs.young),PICs_Major_CellID)
PICs_Major_CellID.aged <- intersect(colnames(rawCounts_PICs.aged),PICs_Major_CellID)
# subset true PICs by conditions
MajorPICs.young <- rawCounts_PICs.young[,PICs_Major_CellID.young]
MajorPICs.aged <- rawCounts_PICs.aged[,PICs_Major_CellID.aged]

# Step 5: Create new Seurat objects by including AD and true PICs
# sample_list.major.AD <- c("MajorAD.young","MajorAD.aged","MajorPICs.young","MajorPICs.aged")
# Sample.list.major.AD <- list()
# for (s in sample_list.major.AD){
#   Sample.list.major.AD[[s]] <- get(s)
# }
# # create seuratobject for integration
# for (i in names(Sample.list.major.AD)){
#   Sample.list.major.AD[[i]] <- CreateSeuratObject(Sample.list.major.AD[[i]], project=i,
#                                                   min.cells = 3, min.features = 200)
#   Sample.list.major.AD[[i]] <- SCTransform(Sample.list.major.AD[[i]], vst.flavor = "v2", method = "glmGamPoi",verbose = F)
# }


# convert to annotated cell types
# function for translating cell types code
Translate_name <- function(toTranslate_name,major_T,major_B,translate_T,translate_B){
  for (t in major_T){
    for (b in major_B){
      T_name <- names(translate_T[translate_T==t])
      B_name <- names(translate_B[translate_B==b])
      toTranslate_name <- gsub(paste0("T",t,"_","B",b), paste0(T_name," and ",B_name), toTranslate_name)
    }
  }
  return(toTranslate_name)
}

# Create a MajorAD young Seurat object
MajorAD_CellID <- colnames(MajorAD.young)
MajorAD_CellID <- gsub("_\\d+$", "", MajorAD_CellID)
MajorAD_CellNames <- Translate_name(toTranslate_name=MajorAD_CellID,
                                    major_T=major_T,major_B=major_B,
                                    translate_T=translate_T,translate_B=translate_B)

MajorAD.young.meta <- data.frame(MajorAD_CellNames)
rownames(MajorAD.young.meta) <- colnames(MajorAD.young)

MajorAD.young.obj <- CreateSeuratObject(counts = MajorAD.young, project = "MajorAD_Young", min.cells = 3, min.features = 200,
                              meta.data = MajorAD.young.meta)
MajorAD.young.obj$orig.ident <- "MajorAD_Young"

# Create a MajorAD aged Seurat object
MajorAD_CellID <- colnames(MajorAD.aged)
MajorAD_CellID <- gsub("_\\d+$", "", MajorAD_CellID)
MajorAD_CellNames <- Translate_name(toTranslate_name=MajorAD_CellID,
                                    major_T=major_T,major_B=major_B,
                                    translate_T=translate_T,translate_B=translate_B)

MajorAD.aged.meta <- data.frame(MajorAD_CellNames)
rownames(MajorAD.aged.meta) <- colnames(MajorAD.aged)

MajorAD.aged.obj <- CreateSeuratObject(counts = MajorAD.aged, project = "MajorAD_Aged", min.cells = 3, min.features = 200,
                                        meta.data = MajorAD.aged.meta)
MajorAD.aged.obj$orig.ident <- "MajorAD_Aged"

# Create a true PICs young Seurat object
MajorPICs_CellID <- colnames(MajorPICs.young)
MajorPICs.young.meta <- data.frame(MajorPICs_CellNames=PICs.combined$predicted_PICs[MajorPICs_CellID])
MajorPICs.young.obj <- CreateSeuratObject(counts = MajorPICs.young, project = "MajorPICs_Young", min.cells = 3, min.features = 200,
                                        meta.data = MajorPICs.young.meta)

# Create a true PICs aged Seurat object
MajorPICs_CellID <- colnames(MajorPICs.aged)
MajorPICs.aged.meta <- data.frame(MajorPICs_CellNames=PICs.combined$predicted_PICs[MajorPICs_CellID])
MajorPICs.aged.obj <- CreateSeuratObject(counts = MajorPICs.aged, project = "MajorPICs_Aged", min.cells = 3, min.features = 200,
                                          meta.data = MajorPICs.aged.meta)




Sample.list.ADPICs <- list(MajorAD_young=MajorAD.young.obj, MajorAD_aged=MajorAD.aged.obj, 
                           MajorPICs_young=MajorPICs.young.obj, MajorPICs_aged=MajorPICs.aged.obj)

for (i in names(Sample.list.ADPICs)){
  Sample.list.ADPICs[[i]] <- SCTransform(Sample.list.ADPICs[[i]], vst.flavor = "v2", method = "glmGamPoi",verbose = F)
}
# Perform integration
features <- SelectIntegrationFeatures(object.list = Sample.list.ADPICs, nfeatures = 3000)
Sample.list.ADPICs <- PrepSCTIntegration(object.list = Sample.list.ADPICs, anchor.features = features)
Sample.list.ADPICs <- lapply(X = Sample.list.ADPICs, FUN = RunPCA, features = features)

# find the anchor by using reference
Anchors.ADPICs <- FindIntegrationAnchors(object.list = Sample.list.ADPICs, normalization.method = "SCT",
                                      anchor.features = features, reduction = "rpca")
Sample.ADPICs <- IntegrateData(anchorset = Anchors.ADPICs, normalization.method = "SCT")

save.image(file = "PICseq_Seurat_revision.RData")

library(stringr)
# # Update the dataframe column
# MajorAD_CellNames <- str_replace_all(
#   Sample.ADPICs$MajorAD_CellNames,
#   c("CD8\\+ EM T Cells" = "CD8+ Ex/Mem T Cells",
#     "CX4CR1_hi CD4\\+ EM T Cells" = "CD4+ Ex/Mem T Cells",
#     "Nlrp6\\+ CD4\\+ EM T Cells" = "CD4+ IFITM+ CISH+ T Cells",
#     "Treg"="Tregs"))
# 
# Sample.ADPICs$MajorAD_CellNames <- MajorAD_CellNames

cell_types <- ifelse(is.na(Sample.ADPICs$MajorAD_CellNames), Sample.ADPICs$MajorPICs_CellNames, Sample.ADPICs$MajorAD_CellNames)
#cell_types <- gsub("^CX4","CX3",cell_types)
Sample.ADPICs$cell_types <- cell_types

Idents(Sample.ADPICs) <- "cell_types"
# Dimensional reduction
Sample.ADPICs <- RunPCA(Sample.ADPICs, npcs = 30 ,verbose = FALSE)
ElbowPlot(Sample.ADPICs,ndims = 30)
Sample.ADPICs <- RunUMAP(Sample.ADPICs, reduction = "pca", dims = 1:30)
# Cluster the cells
Sample.ADPICs <- FindNeighbors(Sample.ADPICs, dims = 1:30)
Sample.ADPICs <- FindClusters(Sample.ADPICs, resolution = 0.3)
Idents(Sample.ADPICs) <- "cell_types"
DimPlot(Sample.ADPICs, reduction = "umap",split.by = "orig.ident",label =F,ncol=2)
ggsave("AD_PICs_combined.pdf",width=16,height=12)

# Prepare object to run differential expression on SCT assay with multiple models
Sample.ADPICs <- PrepSCTFindMarkers(Sample.ADPICs)

AD_PICs <- ifelse(Sample.ADPICs$orig.ident %in% c("MajorAD_Young","MajorAD_Aged"), "AD", "PICs")
Sample.ADPICs$AD_PICs <- AD_PICs

# find markers for every cluster compared to all remaining cells
Idents(Sample.ADPICs) <- "orig.ident"
PICs_vs_Expect_Young <- FindMarkers(Sample.ADPICs, assay = "SCT", ident.1 = "MajorPICs_Young", ident.2 = "MajorAD_Young",
                           min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="AD_PICs")
write.csv(PICs_vs_Expect_Young,"Overall_Obs_vs_Expect_Young.csv")

PICs_vs_Expect_Aged <- FindMarkers(Sample.ADPICs, assay = "SCT", ident.1 = "MajorPICs_Aged", ident.2 = "MajorAD_Aged",
                                    min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="AD_PICs")
write.csv(PICs_vs_Expect_Aged,"Overall_Obs_vs_Expect_Aged.csv")

save.image("/data/array2/fwu/PICs/PICseq_Seurat_revision.RData")


Sample.ADPICs$age <- ifelse(Sample.ADPICs$orig.ident %in% c("MajorAD_Young","MajorPICs_Young"), "Young", "Aged")
Sample.ADPICs$CT_groups <- paste0(Sample.ADPICs$orig.ident,"_",Sample.ADPICs$cell_types)

# DEGs for each PICs vs Expected for young and aged
Idents(Sample.ADPICs) <- "CT_groups"
# young
dir.create("DEGs_PICs_vs_Exp_young")
for (c in unique(Sample.ADPICs$cell_types)){
  print(c)
  try({
  PICs_vs_Expect <- FindMarkers(Sample.ADPICs, assay = "SCT", ident.1 = paste0("MajorPICs_Young","_",c), ident.2 = paste0("MajorAD_Young","_",c),
                                     min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="AD_PICs")
  c <- gsub("/","_or_",c); print(c)
  write.csv(PICs_vs_Expect, paste0("DEGs_PICs_vs_Exp_young/",c,"_Obs_vs_Expect_Young.csv"))
  })
}

# aged
dir.create("DEGs_PICs_vs_Exp_aged")
for (c in unique(Sample.ADPICs$cell_types)){
  print(c)
  try({
    PICs_vs_Expect <- FindMarkers(Sample.ADPICs, assay = "SCT", ident.1 = paste0("MajorPICs_Aged","_",c), ident.2 = paste0("MajorAD_Aged","_",c),
                                  min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST", latent.vars="AD_PICs")
    c <- gsub("/","_or_",c); print(c)
    write.csv(PICs_vs_Expect, paste0("DEGs_PICs_vs_Exp_aged/",c,"_Obs_vs_Expect_Aged.csv"))
  })
}

save.image("/data/array2/fwu/PICs/PICseq_Seurat_revision.RData")

### Observed vs Expected gene expression ###
t_data.exp.all <- T_cell.combined[["RNA"]]$counts
b_data.exp.all <- B_cell.combined[["RNA"]]$counts

SeuratObj_metadata <- Sample.ADPICs@meta.data
SeuratObj_data <- Sample.ADPICs[["SCT"]]$data

# function for gene expression levels
obs_vs_exp_gene <- function(condition, Cell1data, Cell2data, SeuratObj_metadata,SeuratObj_data,randSize,groupAD,groupPICs,CT,AD_label=AD_label,t=t,b=b){
  # random sampling T and B cells
  t_rand <- sample(Cell1data,randSize,replace = T)
  b_rand <- sample(Cell2data,randSize,replace = T)
  # get gene expression matrix
  t_data <- t_data.exp.all[,t_rand]
  b_data <- b_data.exp.all[,b_rand]
  
  # proportion of each gene in the doublet
  T_prop <- rowMeans(t_data)/(rowMeans(t_data) + rowMeans(b_data))
  B_prop <- 1-T_prop
  
  # get expected cell ID
  if (AD_label=="T and B cells"){
    MajorAD_IDs <- SeuratObj_metadata %>% filter(orig.ident==groupAD) %>% rownames()
  } else {
    MajorAD_IDs <- SeuratObj_metadata %>% filter(CT_groups==paste0(groupAD,"_",CT)) %>% rownames()
  }
  # get expected gene expression
  data_AD_exp <- SeuratObj_data[,MajorAD_IDs]
  # get mean expression for each of gene
  expected_expression <- rowMeans(data_AD_exp)
  
    # get observed PICs gene expression
  if (AD_label=="T and B cells"){
    MajorPICs_IDs <- SeuratObj_metadata %>% filter(orig.ident==groupPICs) %>% rownames()
  } else {
    MajorPICs_IDs <- SeuratObj_metadata %>% filter(CT_groups==paste0(groupPICs,"_",CT)) %>% rownames()
  }
    data_PICs_obs <- SeuratObj_data[,MajorPICs_IDs]
    # get mean expression for each of gene
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
    PICs_ct_Young <- Sample.ADPICs@meta.data %>% filter(CT_groups==paste0("MajorPICs_Young_",AD_label)) %>% row.names() # young
    PICs_ct_Aged <- Sample.ADPICs@meta.data %>% filter(CT_groups==paste0("MajorPICs_Aged_",AD_label)) %>% row.names() # aged
    PICs_ct_overall <- Sample.ADPICs@meta.data %>% filter(AD_PICs == "PICs" & cell_types==AD_label) %>% row.names() # overall
    
    # Young
    if ((length(t_ct_Young) >= 3 | length(b_ct_Young) >= 3) & length(PICs_ct_Young) >= 3){
      print(paste0("generate artificial doublets (AD): ",AD_label," Young"))
      obs_vs_exp_gene_young <- obs_vs_exp_gene(condition="Young",Cell1data=t_ct_Young, Cell2data=b_ct_Young,randSize=randSize, 
                                               SeuratObj_metadata=SeuratObj_metadata,SeuratObj_data=SeuratObj_data,
                                               groupAD="MajorAD_Young", groupPICs="MajorPICs_Young", CT=AD_label,AD_label=AD_label)
    }
    # Aged
    if ((length(t_ct_Aged) >= 3 | length(b_ct_Aged) >= 3) & length(PICs_ct_Aged) >= 3){
      print(paste0("generate artificial doublets (AD): ",AD_label," Aged"))
      obs_vs_exp_gene_aged <- obs_vs_exp_gene(condition="Aged",Cell1data=t_ct_Aged, Cell2data=b_ct_Aged,randSize=randSize,
                                              SeuratObj_metadata=SeuratObj_metadata,SeuratObj_data=SeuratObj_data,
                                              groupAD="MajorAD_Aged", groupPICs="MajorPICs_Aged", CT=AD_label,AD_label=AD_label)
    }
    # # overall
    # if ((length(t_ct_overall) >= 3 | length(b_ct_overall) >= 3) & length(PICs_ct_overall) >= 3){
    #   print(paste0("generate artificial doublets (AD): ",AD_label," overall"))
    #   obs_vs_exp_gene_overall <- obs_vs_exp_gene(condition="Overall",Cell1data=t_ct_overall, Cell2data=b_ct_overall,randSize=randSize,
    #                                              SeuratObj_metadata=SeuratObj_metadata,SeuratObj_data=SeuratObj_data,
    #                                              groupAD="MajorAD_Young", groupPICs="MajorPICs_Young", CT=AD_label)
    #}
    # add to a list
    combined_data.exp.all.ls[[AD_label]] <- rbind(obs_vs_exp_gene_young,obs_vs_exp_gene_aged)
  }
}


combined_data.exp.major.ls <- combined_data.exp.all.ls[names(combined_data.exp.all.ls) %in% unique(Sample.ADPICs$cell_types)]
combined_data.exp.major.df <- purrr::reduce(combined_data.exp.major.ls, bind_rows)

# Plotting with error bars added only for "PIC"
test.df <- combined_data.exp.major.df %>% filter(condition != "Overall") %>%
  distinct()

test.df <- test.df %>%
  mutate(
    condition = factor(condition, levels = c("Young", "Aged")),
    PICs = factor(PICs, levels = unique(test.df$PICs)),
    group = interaction(PICs, condition) # Create interaction for grouping
  )

# Get the default color set
default_colors <- scales::hue_pal()(3)
# Change the third color to grey
default_colors[3] <- "darkgrey"

gene2look <- c("Ccl3","Tnf","Hacd4","Stat5a","Tgfb1","Hmgb1","Tlr2")

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



# volcano plot
library(EnhancedVolcano)
vol.plot.young <- EnhancedVolcano(PICs_vs_Expect_Young,
                lab = row.names(PICs_vs_Expect_Young),
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
                title = "Observed vs Expect Young",
                subtitle = bquote(italic("Observed vs Expected")))
ggsave(plot=vol.plot.young,paste0("ObsvsExp_Volcano_Young_overall.pdf"), width=9,height=8)

vol.plot.aged <- EnhancedVolcano(PICs_vs_Expect_Aged,
                            lab = row.names(PICs_vs_Expect_Aged),
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
                            title = "Observed vs Expect Aged",
                            subtitle = bquote(italic("Observed vs Expected")))
ggsave(plot=vol.plot.aged,paste0("ObsvsExp_Volcano_Aged_overall.pdf"), width=9,height=8)


# overall Observed vs Expected gene expression
randSize <- 5000
t_all_Young <- T_cell.combined@meta.data %>% filter(orig.ident=="Young_T_Cells") %>% row.names() # young
t_all_Aged <- T_cell.combined@meta.data %>% filter(orig.ident=="Aged_T_Cells") %>% row.names() # aged

b_all_Young <- B_cell.combined@meta.data %>% filter(orig.ident=="Young_B_Cells") %>% row.names() # young
b_all_Aged <- B_cell.combined@meta.data %>% filter(orig.ident=="Aged_B_Cells") %>% row.names() # aged

# create name of doublet
AD_label_all <- "T and B cells"

# select cell IDs of observed PICs
PICs_all_Young <- Sample.ADPICs@meta.data %>% filter(orig.ident=="MajorPICs_Young") %>% row.names() # young
PICs_all_Aged <- Sample.ADPICs@meta.data %>% filter(orig.ident=="MajorPICs_Aged") %>% row.names() # aged

obs_vs_exp_gene_all_young <- obs_vs_exp_gene(condition="Young",Cell1data=t_all_Young, Cell2data=b_all_Young,randSize=randSize, 
                                             SeuratObj_metadata=SeuratObj_metadata,SeuratObj_data=SeuratObj_data,
                                             groupAD="MajorAD_Young", groupPICs="MajorPICs_Young", CT=AD_label_all,AD_label=AD_label_all,
                                             t="T cell",b="B cell")
obs_vs_exp_gene_all_aged <- obs_vs_exp_gene(condition="Aged",Cell1data=t_all_Aged, Cell2data=b_all_Aged,randSize=randSize, 
                                            SeuratObj_metadata=SeuratObj_metadata,SeuratObj_data=SeuratObj_data,
                                            groupAD="MajorAD_Aged", groupPICs="MajorPICs_Aged", CT=AD_label_all,AD_label=AD_label_all,
                                            t="T cell",b="B cell")
combined_data.exp.overall <- rbind(obs_vs_exp_gene_all_young,obs_vs_exp_gene_all_aged)

# Plotting with error bars added only for "PIC"
combined_data.exp.overall$group <- paste0(combined_data.exp.overall$PICs,"_",combined_data.exp.overall$condition)
test.df <- combined_data.exp.overall
test.df$group <- factor(test.df$group,
                        levels = c("Expected T and B cells_Young","Observed T and B cells_Young",
                                   "Expected T and B cells_Aged","Observed T and B cells_Aged"))

# Get the default color set
default_colors <- scales::hue_pal()(3)
# Change the third color to grey
default_colors[3] <- "darkgrey"

gene2look <- c("Ccl3","Tnf","Hacd4","Stat5a","Tgfb1","Hmgb1","Tlr2")

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
    ggsave(plot=gene_barplot,paste0(g,"_PICs_overall_bar.pdf"), width = 8, height = 8)
    combined_plot.ObsvsExp.bar[[g]] <- gene_barplot
  } else {print(paste0("Abscent gene: ",g))}
}



















randSize <- 5000
MajorAD.df.ls <- list()
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
    t_data <- rawCounts_T[,t_rand]
    b_data <- rawCounts_B[,b_rand]
    # raw counts of artificial doublets (AD)
    data_AD <- t_data + b_data
    colnames(data_AD) <- paste0(AD_label,"_",seq(1,ncol(data_AD)))
    
    MajorAD.df.ls[[AD_label]] <- data_AD
  }
}
MajorAD.df <- purrr::reduce(MajorAD.df.ls, cbind)

MajorAD_CellID <- colnames(MajorAD.df)
MajorAD_CellID <- gsub("_\\d+$", "", MajorAD_CellID)