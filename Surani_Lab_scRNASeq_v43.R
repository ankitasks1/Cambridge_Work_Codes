library(Seurat)
library(sctransform)
library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)
library(scran)
library(tidyverse)
library(BiocParallel)
library(patchwork)
bpp <- MulticoreParam(7)
setwd("/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/seurat/")

# Read data
C8.data <- Read10X(data.dir = "/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/cellranger_outs_v43/run_count_SITTC8/outs/filtered_feature_bc_matrix/")
D8.data <- Read10X(data.dir = "/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/cellranger_outs_v43/run_count_SITTD8/outs/filtered_feature_bc_matrix/")
E1.data <- Read10X(data.dir = "/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/cellranger_outs_v43/run_count_SITTE1/outs/filtered_feature_bc_matrix/")
F1.data <- Read10X(data.dir = "/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/cellranger_outs_v43/run_count_SITTF1/outs/filtered_feature_bc_matrix/")
G1.data <- Read10X(data.dir = "/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/cellranger_outs_v43/run_count_SITTG1/outs/filtered_feature_bc_matrix/")
G9.data <- Read10X(data.dir = "/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/cellranger_outs_v43/run_count_SITTG9/outs/filtered_feature_bc_matrix/")
H10.data <- Read10X(data.dir = "/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/cellranger_outs_v43/run_count_SITTH10/outs/filtered_feature_bc_matrix/")

# Make Seurat object
C8_obj <- CreateSeuratObject(counts = C8.data, project = "C8")
D8_obj <- CreateSeuratObject(counts = D8.data, project = "D8")
E1_obj <- CreateSeuratObject(counts = E1.data, project = "E1")
F1_obj <- CreateSeuratObject(counts = F1.data, project = "F1")
G1_obj <- CreateSeuratObject(counts = G1.data, project = "G1")
G9_obj <- CreateSeuratObject(counts = G9.data, project = "G9")
H10_obj <- CreateSeuratObject(counts = H10.data, project = "H10")

SLscp <- merge(C8_obj, y = c(D8_obj, E1_obj, F1_obj, G1_obj, G9_obj, H10_obj), 
              add.cell.ids = c("C8", "D8", "E1", "F1", "G1", "G9", "H10"), project = "SLscp")


# Examine merged Seurat object
SLscp

# @meta.data slot
head(SLscp@meta.data)


# @assays slot
SLscp@assays

# Default assay
DefaultAssay(SLscp) <- 'RNA'
SLscp@active.assay

# Quality Control (reads mapping to mitochondria)
# Percent of reads mapping to the mitochondrial genome

SLscp[["percent.mt"]] <- PercentageFeatureSet(SLscp, pattern = "^MT-")

head(SLscp@meta.data)

# Plot raw counts per sample
VlnPlot(SLscp, features = c("nCount_RNA"), log=TRUE)

# Total features per sample
VlnPlot(SLscp, features = c("nFeature_RNA"))

# Percent mitochondrial reads per sample
VlnPlot(SLscp, features = c("percent.mt"))

# How many genes were detected in the whole dataset
table(rowSums(SLscp@assays$RNA) > 0)

# Filtering cells
plot1 <- FeatureScatter(SLscp, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SLscp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filtering out outlying cells
cells_to_filter <- rownames(subset(SLscp, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)@meta.data)
SLscp$keep <- rownames(SLscp@meta.data) %in% cells_to_filter

# How many cells are we keeping?
table(SLscp$keep)

# Let's check which cells we are filtering
# Raw counts per sample
VlnPlot(SLscp, features = c("nCount_RNA"), log=TRUE, split.by="keep")

# Total features per sample
VlnPlot(SLscp, features = c("nFeature_RNA"), split.by="keep")

# Percent mitochondrial reads per sample
VlnPlot(SLscp, features = c("percent.mt"), split.by="keep")

# Remove low quality cells from Seurat object
# If cut-off are satisfactory, remove cells
SLscp <- subset(SLscp, subset = keep == TRUE)
SLscp

# Normalisation
SLscp <- SCTransform(SLscp, vars.to.regress = "percent.mt", verbose = TRUE, variable.features.n = 3000)

# The normalised data are now under $SCT in the assay slot
SLscp@assays

## Dimension reduction

# Variable features
top10 <- head(VariableFeatures(SLscp), 10)
plot1 <- VariableFeaturePlot(SLscp)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Linear dimension reduction (PCA)
SLscp <- RunPCA(SLscp, features = VariableFeatures(object = SLscp))

# examine pca plots
SLscp@reductions

# Plot the first two PCA components (classic PCA plot)
DimPlot(SLscp, reduction = "pca", group.by="orig.ident")


# Picking the number of components to use for UMAP projection
# We need to pick a number of PCs which capture most of the variance in the data. These will be used to generate a UMAP plot.
ElbowPlot(SLscp)

# Non-linear dimension reduction - UMAP
# A large proportion of the variation is captured by 16 components and there is a drop off thereafter
SLscp <- RunUMAP(SLscp, dims = 1:16)
head(SLscp@meta.data)
DimPlot(SLscp, group.by="orig.ident")

FeaturePlot(SLscp, features = c('CD79A'))

# Differentially expressed genes between samples
# SLscp_markers <- FindMarkers(SLscp, 
#                              assay = "SCT", 
#                              ident.1 = "ExE.Meso", 
#                              ident.2 = c("EPI", "Ery2","Ery1"),
#                              group.by="cell_type", 
#                              min.pct = 0.5)

# QC and exploratory analysis
samplesheet <- read.csv("/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/SLX-23228.HYC7MDSX5.s_4.contents.csv")
samplesheet["Sample"] <- samplesheet$Barcode
# parallelization for using multicore
bp.params <- MulticoreParam(workers = 7)

# load counts 
C8.counts <- read10xCounts("/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/cellranger_outs_v43/run_count_SITTC8/outs/filtered_feature_bc_matrix/", col.names=TRUE, BPPARAM = bp.params)
dim(counts(C8.counts))

counts(C8.counts)[1:10, 1:10]

# Details about features
rowData(C8.counts)

# Droplet annotation
colData(C8.counts)
colnames(counts(C8.counts))[1:6]

# Properties of scRNA-seq data
# Number of genes detected per cell
C8.genesPerCell <- colSums(counts(C8.counts) > 0)
plot(density(C8.genesPerCell), main="", xlab="Genes per cell")


# Total UMI for a gene versus the number of times detected
plot(rowSums(counts(C8.counts)) / rowSums(counts(C8.counts) > 0),
     rowMeans(counts(C8.counts) > 0),
     log = "x",
     xlab="Mean UMIs per cell",
     ylab="proportion of cells expressing the gene"
)

# Distribution of counts for a gene across cells
# top 20 genes
C8.rel_expression <- t( t(counts(C8.counts)) / colSums(counts(C8.counts))) * 100
rownames(C8.rel_expression) <- rowData(C8.counts)$Symbol
C8.most_expressed <- sort(rowSums( C8.rel_expression ), decreasing = T)[20:1]
C8.plot_data <- as.matrix(t(C8.rel_expression[names(C8.most_expressed),]))

boxplot(C8.plot_data, cex=0.1, las=1, xlab="% total count per cell", horizontal=TRUE)



# Load multiple samples
samples <- samplesheet$Barcode
list_of_files <- str_c("/mnt/home3/reid/av638/scrnaseq/surani_lab/Sabine/SLX-23228/cellranger_outs/", 
                       paste0("run_count_", samples), 
                       "/outs/filtered_feature_bc_matrix")
names(list_of_files) <- samples
list_of_files

SLscp.counts <- read10xCounts(list_of_files, col.names=TRUE, BPPARAM = bp.params)
SLscp.counts
colData(SLscp.counts)
SLscp.counts$Barcode <- rownames(colData(SLscp.counts))
colData(SLscp.counts) <- merge(colData(SLscp.counts), samplesheet, by="Sample", sort=FALSE)
rownames(colData(SLscp.counts)) <- SLscp.counts$Barcode
# Although the data has 36601 genes but there might be many of these will not have been detected in any droplet.

detected_genes <- rowSums(counts(SLscp.counts)) > 0
table(detected_genes)

# We can remove these before proceeding in order to reduce the size of the single cell experiment object.
SLscp.counts <- SLscp.counts[detected_genes,]

# Annotate genes
AnnotationHub <- AnnotationHub()
ens.hs.107 <- query(AnnotationHub, c("Homo sapiens", "EnsDb", 107))[[1]] 

genes <- rowData(SLscp.counts)$ID
gene_annot <- AnnotationDbi::select(ens.hs.107, 
                                    keys = genes,
                                    keytype = "GENEID",
                                    columns = c("GENEID", "SEQNAME")) %>%
  set_names(c("ID", "Chromosome"))
rowData(SLscp.counts) <- merge(rowData(SLscp.counts), gene_annot, by = "ID", sort=FALSE)
rownames(rowData(SLscp.counts)) <- rowData(SLscp.counts)$ID

rowData(SLscp.counts)

# Add per cell QC metrics
is.mito <- which(rowData(SLscp.counts)$Chromosome=="MT")

SLscp.counts <- addPerCellQC(SLscp.counts, subsets=list(Mito=is.mito), BPPARAM = bp.params)
colData(SLscp.counts)
# QC metric distribution

plotColData(SLscp.counts, x="Sample", y="sum") + 
  scale_y_log10() + 
  ggtitle("Total count")

plotColData(SLscp.counts, x="Sample", y="detected") + 
  scale_y_log10() + 
  ggtitle("Detected features")

plotColData(SLscp.counts, x="Sample", y="subsets_Mito_percent") + 
  ggtitle("Mito percent")

# Why normalize ? 
table(SLscp.counts$Sample)

colData(SLscp.counts) %>% 
  as.data.frame() %>% 
  arrange(subsets_Mito_percent) %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(colour = subsets_Mito_percent > 10)) + 
  facet_wrap(vars(Sample.name))



oneSamTab <- colData(SLscp.counts) %>% 
  as.data.frame() %>% 
  dplyr::select(Sample,Barcode.x, sum) %>% # Barcode.x is cell Barcode.y the samplesheet has this name so I kept it as such but can be change to avoid confusion with real cell Barcode
  mutate(cell_num = 1:n())

p_before_nom <- ggplot(data=oneSamTab, aes(x=cell_num, y=sum)) +
  geom_bar(stat = 'identity') +
  labs( x= 'Cell Index',
        y='Cell UMI counts',
        title = "Data: Before Normalization" ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size=20, color = 'red')
  )

p_before_nom

# Cluster cells

set.seed(100)
clust <- quickCluster(SLscp.counts, BPPARAM=bpp)
table(clust)


# Compute size factors
SLscp.counts <- computePooledFactors(SLscp.counts,
                            clusters = clust,
                            min.mean = 0.1,
                            BPPARAM = bpp)
deconv.sf <- sizeFactors(SLscp.counts)
summary(deconv.sf)


lib.sf <- librarySizeFactors(SLscp.counts)
data.frame(LibrarySizeFactors = lib.sf, 
           DeconvolutionSizeFactors = deconv.sf,
           SampleGroup = sce$SampleGroup) %>%
  ggplot(aes(x=LibrarySizeFactors, y=DeconvolutionSizeFactors)) +
  geom_point(aes(col=SampleGroup)) +
  geom_abline(slope = 1, intercept = 0)


#Apply size factors
SLscp.counts <- logNormCounts(SLscp.counts)
assayNames(SLscp.counts)


