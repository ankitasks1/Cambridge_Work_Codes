library(Seurat)

setwd("/mnt/beegfs6/home3/reid/av638/rnaseq/fengzhu_lab/")

# Load raw expression data
human_raw_reads <- readRDS("/mnt/beegfs6/home3/reid/av638/rnaseq/fengzhu_lab/human/human-gastrula-shiny-main/raw_reads.rds")


# Load UMAP coordinates and cell type annotations
human_umap <- readRDS("/mnt/beegfs6/home3/reid/av638/rnaseq/fengzhu_lab/human/human-gastrula-shiny-main/umap.rds")
rownames(human_umap) <- human_umap$cell_name
#Give cell names to rows of reads
rownames(human_raw_reads) <- human_umap$cell_name

#transpose the data
human_raw_reads <- t(human_raw_reads)
dim(human_raw_reads)


# Create a Seurat object
human_seurat_obj <- CreateSeuratObject(human_raw_reads)
# Normalize and filter the data
# human_seurat_obj <- NormalizeData(human_seurat_obj)
human_seurat_obj <- SCTransform(human_seurat_obj, verbose = TRUE, variable.features.n = 3000)

table(rowSums(human_seurat_obj@assays$RNA) > 0)

human_seurat_obj@assays
human_seurat_obj <- FindVariableFeatures(human_seurat_obj, selection.method = "vst", nfeatures = 3000)
human_seurat_obj <- ScaleData(human_seurat_obj)
human_seurat_obj <- RunPCA(human_seurat_obj)
# human_seurat_obj <- FindNeighbors(human_seurat_obj, dims = 1:20)
# human_seurat_obj <- FindClusters(human_seurat_obj)
# Skip RunUMAP  becuase umap is already available
# human_seurat_obj <- RunUMAP(human_seurat_obj, reduction = "pca", dims = 1:30)


# Add UMAP coordinates to the Seurat object
nrow(human_umap) == ncol(human_seurat_obj@assays$RNA@counts) # should return TRUE
human_UMAP_coordinates <- human_umap[,c(2,3)]

# conver to matrix
human_UMAP_coordinates_mat <- as(human_UMAP_coordinates, "matrix")
colnames(UMAP_coordinates_mat) <- c("UMAP_1", "UMAP_2")
# Create DimReducObject and add to object
human_seurat_obj[['umap']] <- CreateDimReducObject(embeddings = human_UMAP_coordinates_mat, key = "UMAP_", global = T, assay = "RNA")

#Add cluster id
human_seurat_obj[['cell_type']]<- human_umap$cluster_id

#Add subcluster information
human_seurat_obj[['sub_cluster']]<- human_umap$sub_cluster

# Visualize the data
DimPlot(human_seurat_obj, reduction = "umap", group.by = "orig.ident")
DimPlot(human_seurat_obj, reduction = "umap", group.by = "cell_type")
DimPlot(human_seurat_obj, reduction = "umap", group.by = "sub_cluster")
# FeaturePlot(human_seurat_obj, features = c("gene1", "gene2"), pt.size = 0.1)

FeaturePlot(human_seurat_obj, reduction = "umap", feature = c("nFeature_SCT"))

# change default assay
DefaultAssay(human_seurat_obj) <- "SCT"

# Make a vector of gene names for the markers we know about
# Let's pretend we don't know about NKT cell markers here
features <- c("FRZB","ANXA1","POSTN","RGS16","DCN","IGFBP7","IGF2","AC132217.1")

# Plot the expression values of each markers on the UMAP
FeaturePlot(human_seurat_obj, reduction = "umap", features = features, pt.size = 0.1, label = TRUE)

# Draw violin plots of the distribution of expression values
# for each marker in each cluster 
VlnPlot(human_seurat_obj, features = features)

# Differentially expressed genes between samples
human_ExE_Mesoderm_Epiblast_markers <- FindMarkers(human_seurat_obj, 
                                    assay = "SCT", 
                                    ident.1 = "ExE Mesoderm", 
                                    ident.2 = c("Epiblast", "Erythroblasts"), 
                                    group.by="cell_type", 
                                    min.pct = 0.5)

head(human_ExE_Mesoderm_Epiblast_markers)

diffgenes <- c("H19")
FeaturePlot(human_seurat_obj, reduction = "umap", features = diffgenes, pt.size = 0.5, label = FALSE)

setgenes <- c("RAB20")
FeaturePlot(human_seurat_obj, reduction = "umap", features = setgenes, pt.size = 0.5, label = FALSE)


# DLL1, DOC2B, DKK1 and WNT8C

#Chicken Pre-primitive chick embryo lee et al.
library(data.table)
fpkm_lee_chick_embryo <- fread("./chicken/fpkm_lee_chick_embryo.txt")
dim(fpkm_lee_chick_embryo)
head(fpkm_lee_chick_embryo)
fpkm_lee_chick_embryo <- data.frame(fpkm_lee_chick_embryo)
fpkm_lee_chick_embryo["AO_MZ"] <- fpkm_lee_chick_embryo$HC1_WAO / fpkm_lee_chick_embryo$HC2_WMZ
fpkm_lee_chick_embryo["AO_AP"] <- fpkm_lee_chick_embryo$HC1_WAO / fpkm_lee_chick_embryo$HC3_WAP
fpkm_lee_chick_embryo["AAO_AMZ"] <- fpkm_lee_chick_embryo$HC4_AAO / fpkm_lee_chick_embryo$HC5_AMZ
fpkm_lee_chick_embryo["AAO_AAP"] <- fpkm_lee_chick_embryo$HC4_AAO / fpkm_lee_chick_embryo$HC6_AAP
fpkm_lee_chick_embryo["PAO_PMZ"] <- fpkm_lee_chick_embryo$HC7_PAO / fpkm_lee_chick_embryo$HC8_PMZ
fpkm_lee_chick_embryo["PAO_PAP"] <- fpkm_lee_chick_embryo$HC7_PAO / fpkm_lee_chick_embryo$HC9_PAP

#filter FC> 3 , dplyr way checked

library(dplyr)
fpkm_lee_chick_embryo_filtFC <- fpkm_lee_chick_embryo %>%
  filter(AO_MZ > 3 & AO_AP > 3 & AAO_AMZ > 3 & AAO_AAP > 3 & PAO_PMZ > 3 & PAO_PAP > 3)

fpkm_lee_chick_embryo_filtFC_FPKMfilt <- fpkm_lee_chick_embryo_filtFC %>%
  filter(HC1_WAO >= 10 & HC4_AAO >= 10 & HC7_PAO >= 10 & HC12_GWM >= 10)

#Monkey
macaca_sc_data <- Read10X(data.dir = "/mnt/home3/reid/av638/rnaseq/fengzhu_lab/monkey/GSE193007_MFE56636-processed/MFE56636MTX/", gene.column=1)
macaca_seurat_obj <- CreateSeuratObject(counts = macaca_sc_data, project = "macaca_sc")

macaca_seurat_obj
head(macaca_seurat_obj@meta.data)
macaca_seurat_obj@assays
DefaultAssay(macaca_seurat_obj) <- 'RNA'
macaca_seurat_obj@active.assay

# # Skip these steps as I need to add metadata, otheriwse these steps are required for processing step-by-step data
# # Identify percentage of mitochondrial reads in each cell
# macaca_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(macaca_seurat_obj, pattern = "^M")
# head(macaca_seurat_obj@meta.data)

#Explore plots
VlnPlot(macaca_seurat_obj, features = c("nCount_RNA"), log=TRUE)
VlnPlot(macaca_seurat_obj, features = c("nFeature_RNA"))
VlnPlot(macaca_seurat_obj, features = c("percent.mt"))

#How many gene detected per sample
# table(rowSums(macaca_seurat_obj@assays$RNA) > 0)

# filter the cells 

# macaca_seurat_obj <- subset(macaca_seurat_obj, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
# # Subsample the filtered cells to reduce dataset size
# cells= 800
# macaca_seurat_obj <- macaca_seurat_obj[, sample(colnames(macaca_seurat_obj), size = cells, replace=F)]

# macaca_seurat_objplot1 <- FeatureScatter(macaca_seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
# macaca_seurat_objplot2 <- FeatureScatter(macaca_seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# macaca_seurat_objplot1 + macaca_seurat_objplot2

#Filter cells (dead cells , apototic)
# macaca_cells_to_filter <- rownames(subset(macaca_seurat_obj, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)@meta.data)
# macaca_seurat_obj$keep <- rownames(macaca_seurat_obj@meta.data) %in% macaca_cells_to_filter


#How many cells are we keeping
# table(macaca_seurat_obj$keep)

# #Let's check which cells we are filtering
# #Raw counts per sample
# VlnPlot(macaca_seurat_obj, features = c("nCount_RNA"), log=TRUE, split.by="keep")
# 
# # Total features per sample
# VlnPlot(macaca_seurat_obj, features = c("nFeature_RNA"), split.by="keep")
# 
# #Percent mitochondrial reads per sample
# VlnPlot(macaca_seurat_obj, features = c("percent.mt"), split.by="keep")
# 

# If cutoffs are satisfactory for the filtering, we can remove those cells which didn't pass. Note that features here are genes and samples are cells.
# macaca_seurat_obj <- subset(macaca_seurat_obj, subset = keep == TRUE)
# macaca_seurat_obj
# Normalize and filter the data
macaca_seurat_obj <- SCTransform(macaca_seurat_obj, verbose = TRUE, variable.features.n = 3000)

# table(rowSums(macaca_seurat_obj@assays$RNA) > 0)

macaca_seurat_obj@assays
macaca_seurat_obj <- FindVariableFeatures(macaca_seurat_obj, selection.method = "vst", nfeatures = 3000)
macaca_seurat_obj <- ScaleData(macaca_seurat_obj)
macaca_seurat_obj <- RunPCA(macaca_seurat_obj)
# macaca_seurat_obj <- RunUMAP(macaca_seurat_obj, reduction = "pca", dims = 1:30)

macaca_meta.data <- read.csv('./monkey/GSE193007_MFE56636-processed/MFE56636-meta.csv', row.names = 1) #row names should contain cell barcode
# seurat <- add.meta.data(seurat, meta.data)
macaca_seurat_obj <- AddMetaData(
  object = macaca_seurat_obj,
  metadata = macaca_meta.data)

head(x = macaca_seurat_obj[[]])

# # Add UMAP coordinates to the Seurat object
# nrow(macaca_meta.data) == ncol(macaca_seurat_obj@assays$RNA@counts) # should return TRUE
# macaca_UMAP_coordinates <- macaca_meta.data[,c(9,10)]
# 
# # conver to matrix
# macaca_UMAP_coordinates_mat <- as(macaca_UMAP_coordinates, "matrix")
# colnames(macaca_UMAP_coordinates_mat) <- c("UMAP_1", "UMAP_2")

# or directly add
# Create DimReducObject and add to object
macaca_seurat_obj[['umap']] <- CreateDimReducObject(embeddings =  as(macaca_meta.data[,c(9,10)], "matrix"), key = "UMAP_", global = T, assay = "RNA")

# make sure we are still using the integrated assay
# DefaultAssay(call_int) <- 'integrated'
# 
# # build the nearest-neighbour graph and do clustering on it
# call_int <- FindNeighbors(call_int, dims = 1:10)
# 
# call_int <- FindClusters(call_int, resolution = 0.3, algorithm=2)
# 
# c8_markers <- FindMarkers(call_int, ident.1 = 7)

# Visualize the data
DimPlot(macaca_seurat_obj, reduction = "umap", group.by = "orig.ident")
DimPlot(macaca_seurat_obj, reduction = "umap", group.by = "cell_type")
DimPlot(macaca_seurat_obj, reduction = "umap", group.by = "stage")
# FeaturePlot(macaca_seurat_obj, features = c("gene1", "gene2"), pt.size = 0.1)

FeaturePlot(macaca_seurat_obj, reduction = "umap", feature = c("nFeature_SCT"))

# change default assay
DefaultAssay(macaca_seurat_obj) <- "SCT"

# Make a vector of gene names for the markers we know about
# Let's pretend we don't know about NKT cell markers here
features <- c("FRZB","ANXA1","POSTN","RGS16","DCN","IGFBP7")

# Plot the expression values of each markers on the UMAP
FeaturePlot(macaca_seurat_obj, reduction = "umap", features = features, pt.size = 0.1, label = TRUE)

# Draw violin plots of the distribution of expression values
# for each marker in each cluster 
VlnPlot(macaca_seurat_obj, features = features)

# Differentially expressed genes between samples
macaca_ExE_Mesoderm_Epiblast_markers <- FindMarkers(macaca_seurat_obj, 
                                                   assay = "SCT", 
                                                   ident.1 = "ExE.Meso", 
                                                   ident.2 = c("EPI", "Ery2","Ery1"), 
                                                   group.by="cell_type", 
                                                   min.pct = 0.5)

head(macaca_ExE_Mesoderm_Epiblast_markers)

diffgenes <- c("SOX2")
FeaturePlot(macaca_seurat_obj, reduction = "umap", features = diffgenes, pt.size = 0.5, label = FALSE)

setgenes <- c("POU5F1")
FeaturePlot(macaca_seurat_obj, reduction = "umap", features = setgenes, pt.size = 0.5, label = FALSE)

VlnPlot(macaca_seurat_obj, features = diffgenes)
VlnPlot(macaca_seurat_obj, features = diffgenes, group.by= "cell_type")


# Developmental programs in chicken early embryos by whole transcriptome analysis
mRNAcounts_gse86592_list <- list()
path_dir_gse86592 <- list.files("/mnt/home3/reid/av638/rnaseq/fengzhu_lab/chicken/Hwang_GSE86592_RAW")
for (files in path_dir_gse86592){
  print(files)
  mRNAcounts_gse86592 <- read.table(paste0("/mnt/home3/reid/av638/rnaseq/fengzhu_lab/chicken/Hwang_GSE86592_RAW/",files))
  mRNAcounts_gse86592_filt <- mRNAcounts_gse86592 %>% filter(grepl('ENSGAL', V1))
  mRNAcounts_gse86592_filt <- mRNAcounts_gse86592_filt[order(mRNAcounts_gse86592_filt$V1),]
  mRNAcounts_gse86592_list[[gsub("_mRNAcount.txt","", files)]] <- mRNAcounts_gse86592_filt
}
