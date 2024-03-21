library(Seurat)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(plyr)
library(ggpubr)
library(splitstackshape)
library(gprofiler2)
library(UpSetR)
library(tidyr)
library(org.Hs.eg.db)
library(stringr)
library(ggrepel)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringi)
library(dplyr)
library(gprofiler2)

1. One to one ortholog (bi-directional)
2. Use only ensembl chicken,human,monkey
3, Intersect using ensembl ID 
setwd("/mnt/beegfs6/home3/reid/av638/rnaseq/fengzhu_lab/")

#read annotations
#use ensembl biomart for homology between chicken and human
chicken_gGal6_human_mart_export <- fread("/mnt/home3/reid/av638/rnaseq/fengzhu_lab/chicken/chicken_gGal6_human_mart_export.txt", header = T)
chicken_gGal6_human_mart_export <- data.frame(chicken_gGal6_human_mart_export)
chicken_gGal6_human_mart_export_sub <- chicken_gGal6_human_mart_export[,c(7,8,1,4)] %>% distinct()
colnames(chicken_gGal6_human_mart_export_sub) <- colnames(human_chicken_hcop_gene_ortholog)
dim(chicken_gGal6_human_mart_export_sub) #29175     4
length(unique(chicken_gGal6_human_mart_export_sub$chicken_ensembl_gene)) #24356 unique ensembl gene ids means there will be some human genes which have multiple chicken genes 

# Remove genes which have no chicken ortholog (ensembl id) in humans
chicken_gGal6_human_mart_export_sub_nohuman_orth <- chicken_gGal6_human_mart_export_sub[which(chicken_gGal6_human_mart_export_sub$human_ensembl_gene != ""),]

# collapse other columns by keeping chicken_ensembl_gene as base so that wherever chicken_ensembl_gene is repeated the other rows will be combined as comma
chicken_gGal6_human_mart_export_sub_one_to <- ddply(chicken_gGal6_human_mart_export_sub_nohuman_orth, .(chicken_ensembl_gene), summarize, human_ensembl_gene = toString(human_ensembl_gene), human_symbol = toString(human_symbol), chicken_symbol = toString(chicken_symbol))
chicken_gGal6_human_mart_export_sub_one_to$human_ensembl_gene <- gsub(" ", "", chicken_gGal6_human_mart_export_sub_one_to$human_ensembl_gene)
chicken_gGal6_human_mart_export_sub_one_to$human_symbol <- gsub(" ", "", chicken_gGal6_human_mart_export_sub_one_to$human_symbol)
chicken_gGal6_human_mart_export_sub_one_to$chicken_symbol <- gsub(" ", "", chicken_gGal6_human_mart_export_sub_one_to$chicken_symbol)

# Count number of human orthologs for chicken ensemble gene
chicken_gGal6_human_mart_export_sub_one_to_one_filt <- chicken_gGal6_human_mart_export_sub_one_to
chicken_gGal6_human_mart_export_sub_one_to_one_filt["human_ensembl_gene_num"] <- unlist(lapply(strsplit(chicken_gGal6_human_mart_export_sub_one_to$human_ensembl_gene, ","), function(x) length(x)))

# Remove human ensembl gene duplicates (eg. ENSG125, ENSG127) of single chicken ensembl id (eg. ENSGAL123)
chicken_gGal6_human_mart_export_sub_one_to_one_filt <- chicken_gGal6_human_mart_export_sub_one_to_one_filt[which(chicken_gGal6_human_mart_export_sub_one_to_one_filt$human_ensembl_gene_num <= 1),]


#I am using the chicken_gGal6_human_mart_export_sub_one_to_one_filt  to remove multiple chicken orthologs  (eg. ENSGAL129, ENSGAL121) to single human ensembl id (eg. ENSG124)
human_gGal6_chicken_mart_export_sub_humanone_to <- ddply(chicken_gGal6_human_mart_export_sub_one_to_one_filt, .(human_ensembl_gene), summarize, chicken_ensembl_gene = toString(chicken_ensembl_gene), human_symbol = toString(human_symbol), chicken_symbol = toString(chicken_symbol), human_ensembl_gene_num = toString(human_ensembl_gene_num))
human_gGal6_chicken_mart_export_sub_humanone_to$chicken_ensembl_gene <- gsub(" ", "", human_gGal6_chicken_mart_export_sub_humanone_to$chicken_ensembl_gene)
human_gGal6_chicken_mart_export_sub_humanone_to$human_symbol <- gsub(" ", "", human_gGal6_chicken_mart_export_sub_humanone_to$human_symbol)
human_gGal6_chicken_mart_export_sub_humanone_to$chicken_symbol <- gsub(" ", "", human_gGal6_chicken_mart_export_sub_humanone_to$chicken_symbol)
human_gGal6_chicken_mart_export_sub_humanone_to_one_filt <- human_gGal6_chicken_mart_export_sub_humanone_to
human_gGal6_chicken_mart_export_sub_humanone_to_one_filt["chicken_ensembl_gene_num"] <- unlist(lapply(strsplit(human_gGal6_chicken_mart_export_sub_humanone_to_one_filt$chicken_ensembl_gene, ","), function(x) length(x)))
human_gGal6_chicken_mart_export_sub_humanone_to_one_filt <- human_gGal6_chicken_mart_export_sub_humanone_to_one_filt[which(human_gGal6_chicken_mart_export_sub_humanone_to_one_filt$chicken_ensembl_gene_num <= 1),]


# ===================================== Human ===================================#
# Dataset 1
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
features <- c("PKP1")

# Plot the expression values of each markers on the UMAP
FeaturePlot(human_seurat_obj, reduction = "umap", features = features, pt.size = 0.1, label = TRUE)

# Draw violin plots of the distribution of expression values
# for each marker in each cluster 
VlnPlot(human_seurat_obj, features = features)

# Differentially expressed genes between samples
# Extra-embryonic erythroblast, hemogenic endothelial, exe mesoderm.
human_differentially_exp_markers <- FindMarkers(human_seurat_obj, 
                                    assay = "SCT", 
                                    ident.1 = as.character(unique(human_seurat_obj@meta.data$cell_type))[(as.character(unique(human_seurat_obj@meta.data$cell_type)) %in% c("ExE Mesoderm", "Hemogenic Endothelium Progenitors")) == TRUE], 
                                    ident.2 = as.character(unique(human_seurat_obj@meta.data$cell_type))[(as.character(unique(human_seurat_obj@meta.data$cell_type)) %in% c("ExE Mesoderm", "Hemogenic Endothelium Progenitors")) == FALSE], 
                                    group.by="cell_type", 
                                    min.pct = 0.5)
human_Extraembryonic_markers <- human_differentially_exp_markers %>% dplyr::filter(p_val_adj < 0.00001)
human_Extraembryonic_markers <- human_differentially_exp_markers %>% dplyr::filter(avg_log2FC > 0)
#get non-extraembryonic markers
human_Extraembryonic_markers <- human_Extraembryonic_markers[order(-human_Extraembryonic_markers$avg_log2FC),]

head(human_Extraembryonic_markers)
dim(human_Extraembryonic_markers)
human_Extraembryonic_markers["human_symbol"] <- rownames(human_Extraembryonic_markers)
top_expressed_gene_human_Extraembryonic <- rownames(head(human_Extraembryonic_markers[order(-human_Extraembryonic_markers$avg_log2FC),],3))
setgenes <- top_expressed_gene_human_Extraembryonic
FeaturePlot(human_seurat_obj, reduction = "umap", features = setgenes, pt.size = 0.5, label = FALSE)

human_Extraembryonic_markers_chicken_orths <- merge(human_Extraembryonic_markers, human_gGal6_chicken_mart_export_sub_humanone_to_one_filt, by.x="human_symbol", by.y="human_symbol")

#============================ Chicken Pre-primitive chick embryo lee et al. ==========================#

#I downloaded raw data and reprocessed it using nextflow
# Dataset 2
load("/mnt/home3/reid/av638/rnaseq/fengzhu_lab/chicken/PRJEB32476/fastqs/fastq/outfolder_ens/star_salmon/deseq2_qc/deseq2.dds.RData")

dds_prjeb32476_ensembl <- dds
rm(dds)
prjeb32476_ensembl_countdata <- counts(dds_prjeb32476_ensembl, normalized=FALSE)

prjeb32476_ensembl_coldata <- fread("/mnt/home3/reid/av638/rnaseq/fengzhu_lab/chicken/PRJEB32476/fastqs/fastq/outfolder_ens/coldata.txt")
prjeb32476_ensembl_coldata <- data.frame(prjeb32476_ensembl_coldata %>% dplyr::select(c("Extract.Name", "ENA_EXPERIMENT", "organism.part", "sampling.site")))

prjeb32476_ensembl_coldata$sampling.site[prjeb32476_ensembl_coldata$sampling.site == "not.applicable"] <- "total"
prjeb32476_ensembl_coldata$organism.part[prjeb32476_ensembl_coldata$organism.part == "Kollers.sickle"] <- "Ksickle"
prjeb32476_ensembl_coldata$Extract.Name[prjeb32476_ensembl_coldata$Extract.Name == "Kollers.Sickle"] <- "Ksickle"
prjeb32476_ensembl_coldata$organism.part[prjeb32476_ensembl_coldata$organism.part == "hypoblast"] <- "other"
prjeb32476_ensembl_coldata$organism.part[prjeb32476_ensembl_coldata$organism.part == "Ksickle"] <- "other"
prjeb32476_ensembl_coldata["extraembryonic"] <- prjeb32476_ensembl_coldata$organism.part

#add extraembryonic specific column

prjeb32476_ensembl_coldata <- prjeb32476_ensembl_coldata %>%
  mutate(extraembryonic = case_when(
    organism.part == "area.opaca" ~ "yes",
    organism.part == "area.pellucida" ~ "no",
    organism.part == "marginal.zone.of.embryo" ~ "no",
    organism.part == "other" ~ "no",
    organism.part == "germ.wall" ~ "germ",
  ))


#remove duplicates
prjeb32476_ensembl_coldata <- prjeb32476_ensembl_coldata %>% distinct()
rownames(prjeb32476_ensembl_coldata) <- colData(dds_prjeb32476_ensembl)$sample
head(prjeb32476_ensembl_coldata)
prjeb32476_ensembl_coldata$Extract.Name <- factor(prjeb32476_ensembl_coldata$Extract.Name)
prjeb32476_ensembl_coldata$organism.part  <- factor(prjeb32476_ensembl_coldata$organism.part)
prjeb32476_ensembl_coldata$ENA_EXPERIMENT <- factor(prjeb32476_ensembl_coldata$ENA_EXPERIMENT)
prjeb32476_ensembl_coldata$sampling.site  <- factor(prjeb32476_ensembl_coldata$sampling.site)
prjeb32476_ensembl_coldata$extraembryonic  <- factor(prjeb32476_ensembl_coldata$extraembryonic)
all(rownames(prjeb32476_ensembl_coldata) == colnames(prjeb32476_ensembl_countdata)) #should print TRUE
prjeb32476_ensembl_ddsO <- DESeqDataSetFromMatrix(countData =prjeb32476_ensembl_countdata, colData = prjeb32476_ensembl_coldata, design = ~ extraembryonic)
prjeb32476_ensembl_keep <- rowSums(counts(prjeb32476_ensembl_ddsO)) > 10
prjeb32476_ensembl_ddsO <- prjeb32476_ensembl_ddsO[prjeb32476_ensembl_keep,]
colData(prjeb32476_ensembl_ddsO)

prjeb32476_ensembl_ddsO <- DESeq(prjeb32476_ensembl_ddsO)
prjeb32476_ensembl_vsdO <- DESeq2::vst(prjeb32476_ensembl_ddsO, blind = FALSE)
DESeq2::plotPCA(prjeb32476_ensembl_vsdO, intgroup = "organism.part")

prjeb32476_ensembl_res_extraembryonic <- results(prjeb32476_ensembl_ddsO, contrast=c("extraembryonic","yes","no"))
prjeb32476_ensembl_res_extraembryonic <- data.frame(prjeb32476_ensembl_res_extraembryonic)
prjeb32476_ensembl_res_extraembryonic["Chicken.ens.name"] <- rownames(prjeb32476_ensembl_res_extraembryonic) 
prjeb32476_ensembl_res_extraembryonic$threshold <- as.logical(prjeb32476_ensembl_res_extraembryonic$padj < 0.05)
prjeb32476_ensembl_res_extraembryonic_0.05 <- prjeb32476_ensembl_res_extraembryonic[which(prjeb32476_ensembl_res_extraembryonic$threshold == TRUE),]
prjeb32476_ensembl_res_extraembryonic_0.05_up_down <- data.frame(prjeb32476_ensembl_res_extraembryonic_0.05) %>% dplyr::filter(log2FoldChange > 0 | log2FoldChange < 0)
prjeb32476_ensembl_res_extraembryonic_0.05_up <- data.frame(prjeb32476_ensembl_res_extraembryonic_0.05) %>% dplyr::filter(log2FoldChange > 0)
prjeb32476_ensembl_res_extraembryonic_0.05_down <- data.frame(prjeb32476_ensembl_res_extraembryonic_0.05) %>% dplyr::filter(log2FoldChange <  0)



# assign orthologs one to one
prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_up_down_gene <- merge(prjeb32476_ensembl_res_extraembryonic_0.05_up_down, human_gGal6_chicken_mart_export_sub_humanone_to_one_filt, by.x= "Chicken.ens.name", by.y="chicken_ensembl_gene")
prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_up_gene <- merge(prjeb32476_ensembl_res_extraembryonic_0.05_up, human_gGal6_chicken_mart_export_sub_humanone_to_one_filt, by.x= "Chicken.ens.name", by.y="chicken_ensembl_gene")
prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_down_gene <- merge(prjeb32476_ensembl_res_extraembryonic_0.05_down, human_gGal6_chicken_mart_export_sub_humanone_to_one_filt, by.x= "Chicken.ens.name", by.y="chicken_ensembl_gene")

write.xlsx(prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_up_down_gene,"/mnt/beegfs6/home3/reid/av638/rnaseq/fengzhu_lab/prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_up_down_gene.xlsx")


#go analys one to one
gprofiler_prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_up_down <- gost(prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_up_down_gene$chicken_symbol, organism="ggallus")
gprofiler_prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_up <- gost(prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_up_gene$chicken_symbol, organism="ggallus")
gprofiler_prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_down <- gost(prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_down_gene$chicken_symbol, organism="ggallus")

# Developmental programs in chicken early embryos by whole transcriptome analysis
#Also perform analysis using the count given by Lab
#use ensmble biomart chicken provided by Lab itself, Galgal 5 version
Galgal5_Biomart <- fread("/mnt/home3/reid/av638/rnaseq/fengzhu_lab/chicken/GSE100798_GSE86592_RAW/Galgal5_Biomart.txt", header = T)
Galgal5_Biomart_df <- data.frame(Galgal5_Biomart)

#use ensmble biomart for homology between chicken and human
# chicken_gGal5_human_mart_export <- fread("/mnt/home3/reid/av638/rnaseq/fengzhu_lab/chicken/chicken_gGal5_human_mart_export.txt", header = T)
# chicken_gGal5_human_mart_export <- data.frame(chicken_gGal5_human_mart_export)
# chicken_gGal5_human_mart_export_sub <- chicken_gGal5_human_mart_export[,c(1,4,8)] %>% distinct()
# colnames(chicken_gGal5_human_mart_export_sub) <- c("Ensembl.Gene.ID", "Human.gene.name", "Chicken.gene.name")

# Dataset 3
#Read raw counts given by the lab
sorted_mRNAcounts_Galgal5 <- fread("/mnt/home3/reid/av638/rnaseq/fengzhu_lab/chicken/GSE100798_GSE86592_RAW/Raw_count_Galgal5.txt", header = T)
sorted_mRNAcounts_Galgal5_df <- data.frame(sorted_mRNAcounts_Galgal5)
rownames(sorted_mRNAcounts_Galgal5_df) <- sorted_mRNAcounts_Galgal5_df$Ensembl_ID
sorted_mRNAcounts_Galgal5_df <- sorted_mRNAcounts_Galgal5_df[,-1]
#Get only required columns
mRNAcounts_Galgal5_bulked_countData = as.matrix(sorted_mRNAcounts_Galgal5_df[,c(1:21)])
Galgal5_bulked_coldata <- data.frame(Sample = colnames(mRNAcounts_Galgal5_bulked_countData), 
                                     condition = unlist(lapply(strsplit(colnames(mRNAcounts_Galgal5_bulked_countData), "_"), function(x) x[1])),
                                     replicate = unlist(lapply(strsplit(colnames(mRNAcounts_Galgal5_bulked_countData), "_"), function(x) x[2])))
rownames(Galgal5_bulked_coldata) <- colnames(mRNAcounts_Galgal5_bulked_countData)
Galgal5_bulked_coldata$condition <- factor(Galgal5_bulked_coldata$condition)
Galgal5_bulked_coldata$replicate <- factor(Galgal5_bulked_coldata$replicate)
all(rownames(Galgal5_bulked_coldata) == colnames(mRNAcounts_Galgal5_bulked_countData)) #should print TRUE
Galgal5_bulked_dds <- DESeqDataSetFromMatrix(countData =mRNAcounts_Galgal5_bulked_countData, colData = Galgal5_bulked_coldata, design = ~ condition)
Galgal5_bulked_dds
Galgal5_bulked_keep <- rowSums(counts(Galgal5_bulked_dds)) >= 10
Galgal5_bulked_dds <- Galgal5_bulked_dds[Galgal5_bulked_keep,]
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

Galgal5_bulked_dds <- DESeq(Galgal5_bulked_dds)
Galgal5_bulked_vsd <- DESeq2::vst(Galgal5_bulked_dds)
plotPCA(Galgal5_bulked_vsd)

#DEG
#See size factors estimated by DESeq function: sizeFactors(dds)
res_Galgal5_bulked_EGKVIII_EGKVI <- results(Galgal5_bulked_dds, contrast=c("condition","EGK.VIII","EGK.VI"))
res_Galgal5_bulked_EGKVIII_EGKVI <- data.frame(res_Galgal5_bulked_EGKVIII_EGKVI)
res_Galgal5_bulked_EGKVIII_EGKVI["Ensembl.Gene.ID"] <- rownames(res_Galgal5_bulked_EGKVIII_EGKVI)
res_Galgal5_bulked_EGKVIII_EGKVI$threshold <- as.logical(res_Galgal5_bulked_EGKVIII_EGKVI$padj < 0.05)
res_Galgal5_bulked_EGKVIII_EGKVI_biomart <- merge(res_Galgal5_bulked_EGKVIII_EGKVI, Galgal5_Biomart_df, by="Ensembl.Gene.ID")
deseq2_results_res_Galgal5_bulked_EGKVIII_EGKVI_0.05 <- res_Galgal5_bulked_EGKVIII_EGKVI_biomart[which(res_Galgal5_bulked_EGKVIII_EGKVI_biomart$threshold == TRUE),]
# deseq2_results_res_Galgal5_bulked_EGKVIII_EGKVI_0.05_up_down <- data.frame(deseq2_results_res_Galgal5_bulked_EGKVIII_EGKVI_0.05) %>% dplyr::filter(log2FoldChange > 0 | log2FoldChange < 0)
# deseq2_results_res_Galgal5_bulked_EGKVIII_EGKVI_0.05_up <- data.frame(deseq2_results_res_Galgal5_bulked_EGKVIII_EGKVI_0.05) %>% dplyr::filter(log2FoldChange > 0)
# deseq2_results_res_Galgal5_bulked_EGKVIII_EGKVI_0.05_down <- data.frame(deseq2_results_res_Galgal5_bulked_EGKVIII_EGKVI_0.05) %>% dplyr::filter(log2FoldChange < 0)
# deseq2_results_res_Galgal5_bulked_EGKVIII_EGKVI_0.05_up_down_gene <- merge(deseq2_results_res_Galgal5_bulked_EGKVIII_EGKVI_0.05_up_down, chicken_gGal5_human_mart_export_sub, by.x= "Ensembl.Gene.ID", by.y="Ensembl.Gene.ID")

res_Galgal5_bulked_EGKX_EGKVIII <- results(Galgal5_bulked_dds, contrast=c("condition","EGK.X","EGK.VIII"))
res_Galgal5_bulked_EGKX_EGKVIII <- data.frame(res_Galgal5_bulked_EGKX_EGKVIII)
res_Galgal5_bulked_EGKX_EGKVIII["Ensembl.Gene.ID"] <- rownames(res_Galgal5_bulked_EGKX_EGKVIII)
res_Galgal5_bulked_EGKX_EGKVIII$threshold <- as.logical(res_Galgal5_bulked_EGKX_EGKVIII$padj < 0.05)
res_Galgal5_bulked_EGKX_EGKVIII_biomart <- merge(res_Galgal5_bulked_EGKX_EGKVIII, Galgal5_Biomart_df, by="Ensembl.Gene.ID")
deseq2_results_res_Galgal5_bulked_EGKX_EGKVIII_0.05 <- res_Galgal5_bulked_EGKX_EGKVIII_biomart[which(res_Galgal5_bulked_EGKX_EGKVIII_biomart$threshold == TRUE),]
deseq2_results_res_Galgal5_bulked_EGKX_EGKVIII_0.05_up_down <- data.frame(deseq2_results_res_Galgal5_bulked_EGKX_EGKVIII_0.05) %>% dplyr::filter(log2FoldChange > 0 | log2FoldChange < 0)
deseq2_results_res_Galgal5_bulked_EGKX_EGKVIII_0.05_up <- data.frame(deseq2_results_res_Galgal5_bulked_EGKX_EGKVIII_0.05) %>% dplyr::filter(log2FoldChange > 0)
deseq2_results_res_Galgal5_bulked_EGKX_EGKVIII_0.05_down <- data.frame(deseq2_results_res_Galgal5_bulked_EGKX_EGKVIII_0.05) %>% dplyr::filter(log2FoldChange < 0)

DE_expressed_gene_chicken_stage_df_one_to_one_human_orth <- merge(deseq2_results_res_Galgal5_bulked_EGKX_EGKVIII_0.05_up_down, human_gGal6_chicken_mart_export_sub_humanone_to_one_filt, by.x="Ensembl.Gene.ID", by.y="chicken_ensembl_gene")
write.xlsx(DE_expressed_gene_chicken_stage_df_one_to_one_human_orth,"/mnt/beegfs6/home3/reid/av638/rnaseq/fengzhu_lab/DE_expressed_gene_chicken_stage_df_one_to_one_human_orth.xlsx")


#Get only required columns single cell embryonic RNAseq data
# Dataset 4
mRNAcounts_Galgal5_sce_countData = as.matrix(sorted_mRNAcounts_Galgal5_df[,-c(1:21)])
Galgal5_sce_coldata <- data.frame(Sample = colnames(mRNAcounts_Galgal5_sce_countData), 
                                  condition = unlist(lapply(strsplit(colnames(mRNAcounts_Galgal5_sce_countData), "_"), function(x) x[1])),
                                  replicate = unlist(lapply(strsplit(colnames(mRNAcounts_Galgal5_sce_countData), "_"), function(x) x[2])))
rownames(Galgal5_sce_coldata) <- colnames(mRNAcounts_Galgal5_sce_countData)
Galgal5_sce_coldata$condition <- factor(Galgal5_sce_coldata$condition)
Galgal5_sce_coldata$replicate <- factor(Galgal5_sce_coldata$replicate)
all(rownames(Galgal5_sce_coldata) == colnames(mRNAcounts_Galgal5_sce_countData)) #should print TRUE
Galgal5_sce_dds <- DESeqDataSetFromMatrix(countData =mRNAcounts_Galgal5_sce_countData, colData = Galgal5_sce_coldata, design = ~ condition)
Galgal5_sce_dds
Galgal5_sce_keep <- rowSums(counts(Galgal5_sce_dds)) >= 10
Galgal5_sce_dds <- Galgal5_sce_dds[Galgal5_sce_keep,]

Galgal5_sce_dds <- DESeq(Galgal5_sce_dds)
Galgal5_sce_vsd <- DESeq2::vst(Galgal5_sce_dds)
plotPCA(Galgal5_sce_vsd)

#DEG
#See size factors estimated by DESeq function: sizeFactors(dds)
res_Galgal5_sce_EGK.X_Zygote <- data.frame(res_Galgal5_sce_EGK.X_Zygote)
res_Galgal5_sce_EGK.X_Zygote["Ensembl.Gene.ID"] <- rownames(res_Galgal5_sce_EGK.X_Zygote)
res_Galgal5_sce_EGK.X_Zygote$threshold <- as.logical(res_Galgal5_sce_EGK.X_Zygote$padj < 0.05)
res_Galgal5_sce_EGK.X_Zygote_biomart <- merge(res_Galgal5_sce_EGK.X_Zygote, Galgal5_Biomart_df, by="Ensembl.Gene.ID")
deseq2_results_res_Galgal5_sce_EGK.X_Zygote_0.05 <- res_Galgal5_sce_EGK.X_Zygote_biomart[which(res_Galgal5_sce_EGK.X_Zygote_biomart$threshold == TRUE),]
deseq2_results_res_Galgal5_sce_EGK.X_Zygote_0.05_up_down <- data.frame(deseq2_results_res_Galgal5_sce_EGK.X_Zygote_0.05) %>% dplyr::filter(log2FoldChange > 0.5 | log2FoldChange < -0.5)
deseq2_results_res_Galgal5_sce_EGK.X_Zygote_0.05_up <- data.frame(deseq2_results_res_Galgal5_sce_EGK.X_Zygote_0.05) %>% dplyr::filter(log2FoldChange > 0.5)
deseq2_results_res_Galgal5_sce_EGK.X_Zygote_0.05_down <- data.frame(deseq2_results_res_Galgal5_sce_EGK.X_Zygote_0.05) %>% dplyr::filter(log2FoldChange < -0.5)
deseq2_results_res_Galgal5_sce_EGK.X_Zygote_0.05_up_down_gene <- merge(deseq2_results_res_Galgal5_sce_EGK.X_Zygote_0.05_up_down, chicken_gGal5_human_mart_export_sub, by.x= "Ensembl.Gene.ID", by.y="Ensembl.Gene.ID")


#======================================= Monkey ===================================#
# Dataset 5
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
DimPlot(macaca_seurat_obj, reduction = "umap", group.by = "cell_cluster")
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
macaca_differentially_exp_markers <- FindMarkers(macaca_seurat_obj, 
                                                   assay = "SCT", 
                                                   ident.1 = as.character(unique(macaca_seurat_obj@meta.data$cell_type))[(as.character(unique(macaca_seurat_obj@meta.data$cell_type)) %in% c("EC","BP","Mac","Exe.Meso","ys.Meso1","ys.Meso2","ys.Endo1","ys.Endo2")) == TRUE], 
                                                   ident.2 = as.character(unique(macaca_seurat_obj@meta.data$cell_type))[(as.character(unique(macaca_seurat_obj@meta.data$cell_type)) %in% c("EC","BP","Mac","Exe.Meso","ys.Meso1","ys.Meso2","ys.Endo1","ys.Endo2")) == FALSE], 
                                                   group.by="cell_type", 
                                                   min.pct = 0.5)
macaca_Extraembryonic_markers <- macaca_differentially_exp_markers %>% dplyr::filter( p_val_adj< 0.00001)
macaca_Extraembryonic_markers <- macaca_differentially_exp_markers %>% dplyr::filter(avg_log2FC > 0)
head(macaca_Extraembryonic_markers)
dim(macaca_Extraembryonic_markers)
top_expressed_gene_macaca_Extraembryonic <- rownames(head(macaca_Extraembryonic_markers[order(-macaca_Extraembryonic_markers$avg_log2FC),],3))
macsetgenes <- top_expressed_gene_macaca_Extraembryonic
FeaturePlot(macaca_seurat_obj, reduction = "umap", features = macsetgenes, pt.size = 0.5, label = FALSE)


monkey_mmuI_human_chicken_ensembl_mart <- fread("/mnt/beegfs6/home3/reid/av638/rnaseq/fengzhu_lab/monkey/monkey_mmuI_human_chicken_ensembl_mart.txt")
colnames(monkey_mmuI_human_chicken_ensembl_mart) <- c("monkey_ensembl_gene", "monkey_symbol", "chicken_ensembl_gene", "chicken_symbol", "human_ensembl_gene", "human_symbol")
monkey_mmuI_human_chicken_ensembl_mart <- monkey_mmuI_human_chicken_ensembl_mart[which(monkey_mmuI_human_chicken_ensembl_mart$chicken_ensembl_gene != ""),]
monkey_mmuI_human_chicken_ensembl_mart <- monkey_mmuI_human_chicken_ensembl_mart[which(monkey_mmuI_human_chicken_ensembl_mart$human_ensembl_gene != ""),]

monkey_mmuI_human_chicken_ensembl_mart_filt <- monkey_mmuI_human_chicken_ensembl_mart
monkey_mmuI_human_chicken_ensembl_mart_filt <- ddply(monkey_mmuI_human_chicken_ensembl_mart, .(monkey_ensembl_gene), summarize, monkey_symbol = toString(monkey_symbol), chicken_ensembl_gene= toString(chicken_ensembl_gene), chicken_symbol= toString(chicken_symbol), human_ensembl_gene= toString(human_ensembl_gene), human_symbol=toString(human_symbol))
monkey_mmuI_human_chicken_ensembl_mart_filt["chicken_monkey_ensembl_gene_num"] <- unlist(lapply(strsplit(monkey_mmuI_human_chicken_ensembl_mart_filt$chicken_ensembl_gene, "_"), function(x) length(x)))
monkey_mmuI_human_chicken_ensembl_mart_filt["human_monkey_ensembl_gene_num"] <- unlist(lapply(strsplit(monkey_mmuI_human_chicken_ensembl_mart_filt$human_ensembl_gene, "_"), function(x) length(x)))

monkey_mmuI_human_chicken_ensembl_mart_filt <- monkey_mmuI_human_chicken_ensembl_mart_filt[which(monkey_mmuI_human_chicken_ensembl_mart_filt$chicken_monkey_ensembl_gene_num <=1),]
monkey_mmuI_human_chicken_ensembl_mart_filt <- monkey_mmuI_human_chicken_ensembl_mart_filt[which(monkey_mmuI_human_chicken_ensembl_mart_filt$human_monkey_ensembl_gene_num <=1),]


monkey_mmuI_human_chicken_ensembl_mart_filt2 <- ddply(monkey_mmuI_human_chicken_ensembl_mart_filt, .(chicken_ensembl_gene), summarize,  monkey_ensembl_gene= toString(monkey_ensembl_gene), monkey_symbol = toString(monkey_symbol), chicken_symbol= toString(chicken_symbol), human_ensembl_gene= toString(human_ensembl_gene), human_symbol=toString(human_symbol))

monkey_mmuI_human_chicken_ensembl_mart_filt2["monkey_chicken_ensembl_gene_num"] <- unlist(lapply(strsplit(monkey_mmuI_human_chicken_ensembl_mart_filt2$chicken_ensembl_gene, "_"), function(x) length(x)))
monkey_mmuI_human_chicken_ensembl_mart_filt2["monkey_human_ensembl_gene_num"] <- unlist(lapply(strsplit(monkey_mmuI_human_chicken_ensembl_mart_filt2$human_ensembl_gene, "_"), function(x) length(x)))

monkey_mmuI_human_chicken_ensembl_mart_filt2 <- monkey_mmuI_human_chicken_ensembl_mart_filt2[which(monkey_mmuI_human_chicken_ensembl_mart_filt2$monkey_chicken_ensembl_gene_num <=1),]
monkey_mmuI_human_chicken_ensembl_mart_filt2 <- monkey_mmuI_human_chicken_ensembl_mart_filt2[which(monkey_mmuI_human_chicken_ensembl_mart_filt2$monkey_human_ensembl_gene_num <=1),]


monkey_mmuI_human_chicken_ensembl_mart_filt3 <- ddply(monkey_mmuI_human_chicken_ensembl_mart_filt2, .(human_ensembl_gene), summarize,  monkey_ensembl_gene= toString(monkey_ensembl_gene), monkey_symbol = toString(monkey_symbol), chicken_symbol= toString(chicken_symbol), chicken_ensembl_gene= toString(chicken_ensembl_gene), human_symbol=toString(human_symbol))

monkey_mmuI_human_chicken_ensembl_mart_filt3["monkey_chicken_ensembl_gene_num"] <- unlist(lapply(strsplit(monkey_mmuI_human_chicken_ensembl_mart_filt3$chicken_ensembl_gene, "_"), function(x) length(x)))
monkey_mmuI_human_chicken_ensembl_mart_filt3["monkey_human_ensembl_gene_num"] <- unlist(lapply(strsplit(monkey_mmuI_human_chicken_ensembl_mart_filt3$human_ensembl_gene, "_"), function(x) length(x)))

monkey_mmuI_human_chicken_ensembl_mart_filt3 <- monkey_mmuI_human_chicken_ensembl_mart_filt3[which(monkey_mmuI_human_chicken_ensembl_mart_filt3$monkey_chicken_ensembl_gene_num <=1),]
monkey_mmuI_human_chicken_ensembl_mart_filt3 <- monkey_mmuI_human_chicken_ensembl_mart_filt3[which(monkey_mmuI_human_chicken_ensembl_mart_filt3$monkey_human_ensembl_gene_num <=1),]


#Human and Monkey Ortholog
macaca_Extraembryonic_markers["monkey_symbol"] <- rownames(macaca_Extraembryonic_markers)
DE_expressed_gene_macaca_Extraembryonic_df_one_to_one_human_orth <- merge(macaca_Extraembryonic_markers, monkey_mmuI_human_chicken_ensembl_mart_filt3, by="monkey_symbol")


#Integration of Dataset 1 to 4
# 
# #ensembl + ncbi
# common_extraembryonic_genes_of_interest <- intersect(intersect(intersect(rownames(human_Extraembryonic_markers),
#                               DE_expressed_gene_macaca_Extraembryonic_df_human_orth$human_symbol),
#                     DE_expressed_gene_chicken_Extraembryonic_df_human_orth$human_symbol),
#           c(prjeb32476_ensembl_res_extraembryonic_0.05_up_gene$human_symbol,
#             prjeb32476_res_extraembryonic_0.05_up_gene$human_symbol))
# gprofiler_common_extraembryonic_genes_of_interest <- gost(common_extraembryonic_genes_of_interest, organism="hsapiens")

#ensembl no monkey
#Upset plot if you want to get data out
common_extraembryonic_goi_listInput_one_to_one <- list(human_extraembryonic = unique(human_Extraembryonic_markers_chicken_orths$human_symbol), 
                                                       chicken_stage = DE_expressed_gene_chicken_stage_df_one_to_one_human_orth$human_symbol,
                                                       chicken_extraembryonic = prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_up_gene$human_symbol)
UpSetR::upset(fromList(common_extraembryonic_goi_listInput_one_to_one), order.by = "freq")
#common_extraembryonic_goi_listInput_one_to_one.png
source("https://raw.githubusercontent.com/ankitasks1/utilities/main/upset_out.R")
upset_extraembryonic_goi_listoutput_one_to_one <- Upsetout(common_extraembryonic_goi_listInput_one_to_one)

gprofiler_human_extraembryonic_chicken_stage_chicken_extraembryonic_genes_of_interest_one_to_one <- gost(rownames(upset_extraembryonic_goi_listoutput_one_to_one$human_extraembryonic_chicken_stage_chicken_extraembryonic), organism="hsapiens")

gprofiler_human_extraembryonic_chicken_stage_0_genes_of_interest_one_to_one <- gost(rownames(upset_extraembryonic_goi_listoutput_one_to_one$human_extraembryonic_chicken_stage_0), organism="hsapiens")

gprofiler_0_chicken_stage_chicken_extraembryonic_genes_of_interest_one_to_one <- gost(rownames(upset_extraembryonic_goi_listoutput_one_to_one$`0_chicken_stage_chicken_extraembryonic`), organism="hsapiens")

#--------------- summarize outputs ------------------#
head(human_gGal6_chicken_mart_export_sub_humanone_to_one_filt)

#Data1
length(unique(human_Extraembryonic_markers$human_symbol))
length(unique(human_Extraembryonic_markers_chicken_orths$human_symbol))

#Data2
length(unique(prjeb32476_ensembl_res_extraembryonic_0.05_up_down$Chicken.ens.name))
length(unique(prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_up_down_gene$Chicken.ens.name))

length(unique(prjeb32476_ensembl_res_extraembryonic_0.05_up$Chicken.ens.name))
length(unique(prjeb32476_ensembl_res_extraembryonic_0.05_one_to_one_up_gene$Chicken.ens.name))

#Data3
length(unique(deseq2_results_res_Galgal5_bulked_EGKX_EGKVIII_0.05_up_down$Ensembl.Gene.ID))
length(unique(DE_expressed_gene_chicken_stage_df_one_to_one_human_orth$Ensembl.Gene.ID))

upset_extraembryonic_goi_listoutput_one_to_one_sel <- upset_extraembryonic_goi_listoutput_one_to_one[names(upset_extraembryonic_goi_listoutput_one_to_one) %in% c("0_chicken_stage_chicken_extraembryonic", "human_extraembryonic_chicken_stage_0", "human_extraembryonic_chicken_stage_chicken_extraembryonic")]

upset_extraembryonic_goi_listoutput_one_to_one_sel_df <- do.call(rbind.data.frame, upset_extraembryonic_goi_listoutput_one_to_one_sel)
upset_extraembryonic_goi_listoutput_one_to_one_sel_df["human_symbol"] <- rownames(upset_extraembryonic_goi_listoutput_one_to_one_sel_df)
upset_extraembryonic_goi_listoutput_one_to_one_sel_df["human_symbol"] <- unlist(lapply(strsplit(upset_extraembryonic_goi_listoutput_one_to_one_sel_df$human_symbol, "\\."), function(x) x[2]))
rownames(upset_extraembryonic_goi_listoutput_one_to_one_sel_df) <- upset_extraembryonic_goi_listoutput_one_to_one_sel_df$human_symbol
upset_extraembryonic_goi_one_to_one_sel_de_stage <- merge(upset_extraembryonic_goi_listoutput_one_to_one_sel_df, DE_expressed_gene_chicken_stage_df_one_to_one_human_orth, by= "human_symbol")
upset_extraembryonic_goi_one_to_one_sel_de_stage <- upset_extraembryonic_goi_one_to_one_sel_de_stage[,c(21,6,5,7:12,20,1,22:23)]
# -> CsCe = Chicken stage and Chicken extraembryonic
# -> HeCs = Chicken stage and Human extarembryonic
# -> HeCsCe =  Chicken stage and Human extarembryonic and Chicken extraembryonic

upset_extraembryonic_goi_one_to_one_sel_de_stage$combos[which(upset_extraembryonic_goi_one_to_one_sel_de_stage$combos == "0_chicken_stage_chicken_extraembryonic")] <-  "CsCe"
upset_extraembryonic_goi_one_to_one_sel_de_stage$combos[which(upset_extraembryonic_goi_one_to_one_sel_de_stage$combos == "human_extraembryonic_chicken_stage_0")] <-  "HeCs"
upset_extraembryonic_goi_one_to_one_sel_de_stage$combos[which(upset_extraembryonic_goi_one_to_one_sel_de_stage$combos == "human_extraembryonic_chicken_stage_chicken_extraembryonic")] <- "HeCsCe"
upset_extraembryonic_goi_one_to_one_sel_de_stage <- upset_extraembryonic_goi_one_to_one_sel_de_stage[,-c(12:13)]
write.xlsx(upset_extraembryonic_goi_one_to_one_sel_de_stage,"/mnt/beegfs6/home3/reid/av638/rnaseq/fengzhu_lab/upset_extraembryonic_goi_one_to_one_sel_de_stage.xlsx")



upset_extraembryonic_goi_listoutput_one_to_one_sel_HeCe <- upset_extraembryonic_goi_listoutput_one_to_one[names(upset_extraembryonic_goi_listoutput_one_to_one) %in% c("human_extraembryonic_0_chicken_extraembryonic")]
rownames(upset_extraembryonic_goi_listoutput_one_to_one_sel_HeCe$human_extraembryonic_0_chicken_extraembryonic)
