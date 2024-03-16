### Load packages
library(devtools)
install_github("ankitasks1/Seqwalker", force=TRUE)
library(Seqwalker)

library(DESeq2)
library(DiffBind)
library(edgeR)
library(ggplot2)
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
library(stringi)
library(data.table)
library(openxlsx)
library(ggfortify)
library(ggvenn)
library(gprofiler2)
library(pheatmap)
library(ComplexHeatmap)
library(ggbreak)
library(Rsamtools)
library(Rsubread)
library(RColorBrewer)


setwd("/mnt/beegfs6/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm")


###########################################
####  common_files among all datasets   ###
###########################################
integration_features_list <- list()
# include upstream of a gene https://support.bioconductor.org/p/78652/
integration_features_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_gencode_out.sorted.chr.txt")
integration_features_list[["gene_ext"]] <- makeGRangesFromDataFrame(integration_features_list[["gene"]], keep.extra.columns=TRUE)
# extend the upstream by 2kb (to capture promoter along with gene body)
integration_features_list$gene_ext <- data.frame(resize(integration_features_list$gene_ext, width = width(integration_features_list$gene_ext) + 2000, fix = "end"))
integration_features_list$gene_ext <- integration_features_list$gene_ext[which(integration_features_list$gene_ext$start > 0),]
write.table(integration_features_list$gene_ext[,c(1,2,3,5:8)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/","gene_gencode_human_upstream_2kb.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)
system("sort -k1,1 -k2,2n gene_gencode_human_upstream_2kb.txt > gene_gencode_human_upstream_2kb.sorted.txt")
# bins5kb
system("bedtools makewindows -g hg38.chrom.sizes -w 5000 | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"\"bin_\"NR\"\\t\"$1\"%\"$2\"%\"$3\"\\t\"\".\"}' > hg38_5kb.txt")
# bins10kb
system("bedtools makewindows -g hg38.chrom.sizes -w 10000 | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"\"bin_\"NR\"\\t\"$1\"%\"$2\"%\"$3\"\\t\"\".\"}' > hg38_10kb.txt")
# bins80kb
system("bedtools makewindows -g hg38.chrom.sizes -w 80000 | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"\"bin_\"NR\"\\t\"$1\"%\"$2\"%\"$3\"\\t\"\".\"}' > hg38_80kb.txt")

# promoter
integration_features_list[["promoter"]] <- data.frame(GenomicRanges::promoters(makeGRangesFromDataFrame(integration_features_list[["gene"]], keep.extra.columns=TRUE), upstream=2000, downstream=200))
# filter chrM
integration_features_list$promoter <- integration_features_list$promoter[which(integration_features_list$promoter$seqnames != "chrM"),]
write.table(integration_features_list$promoter[,c(1,2,3,7,8,5)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/","gene_gencodev41_promoter.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

# tss of genes
integration_features_list[["transcript"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "transcript_v41_human_transcript_grch38.sorted.txt")
integration_features_list[["tss"]] <- data.frame(GenomicRanges::resize(makeGRangesFromDataFrame(integration_features_list[["transcript"]], keep.extra.columns=TRUE), width=2, fix='start'))
write.table(integration_features_list$tss[,c(1,2,3,7,9,5)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/","transcript_gencodev41_tss.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

integration_features_list[["gene_tss"]] <- data.frame(GenomicRanges::resize(makeGRangesFromDataFrame(integration_features_list[["gene"]], keep.extra.columns=TRUE), width=2, fix='start'))
write.table(integration_features_list$gene_tss[,c(1,2,3)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/","gene_gencodev41_genetss.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

write.table(integration_features_list$gene_ext[,c(1,2,3,7,8)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/","gene_gencodev41_gene.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

# rnaseqkd_gene_groups
integration_path <- "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/"


# region of interest RM
# genebody: gene_gencode_human_gencode_out.sorted.chr.txt

# promoter (+2000 , 0): gene_gencodev41_promoter_ideal.txt
# promoter ideal (since 200 bp are interfering with gene body I assumed promoter to be 2000 bp upstream only )
integration_features_list[["promoter_i"]] <- data.frame(GenomicRanges::promoters(makeGRangesFromDataFrame(integration_features_list[["gene"]], keep.extra.columns=TRUE), upstream=2000, downstream=0))
# filter chrM
integration_features_list$promoter_i <- integration_features_list$promoter_i[which(integration_features_list$promoter_i$seqnames != "chrM"),]
write.table(integration_features_list$promoter_i[,c(1,2,3,7,8,5)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/","gene_gencodev41_promoter_ideal.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

system("sort -k1,1 -k2,2n hg38.chrom.sizes > hg38.chrom.sizes.bed")

################################################
###       RNA-Seq data: Knockdown            ###
################################################
rnaseqkd_path <- "/mnt/home3/reid/av638/rnaseq/iva_lab_oct23/"

rnaseqkd_list <- list()

rnaseqkd_list[["dds"]] <- load_rdata(rnaseqkd_path, "outfolder/star_salmon/deseq2_qc/deseq2.dds.RData")
rnaseqkd_list[["gene"]] <- read_genefile(rnaseqkd_path, "gene_gencode_human__gencode_out_chr.txt")
rnaseqkd_list[["counts"]] <- dds_counts(rnaseqkd_list$dds)
rnaseqkd_list[["samples_selected"]] <- sample_filter(rnaseqkd_list$counts, samplefilter = TRUE, c(7:9))
rnaseqkd_list[["processed_counts"]] <- reassign_genenames(rnaseqkd_list$samples_selected, rnaseqkd_list$gene, c(2:10))
rnaseqkd_list[["coldata"]] <- make_coldata_dds(rnaseqkd_list$dds, samplefilter = TRUE, c(7:9))

# check
all(rownames(rnaseqkd_list$coldata) == colnames(rnaseqkd_list$processed_counts)) #should print TRUE

# de analysis
rnaseqkd_list[["deseq2"]] <- deseq2_de(rnaseqkd_list$processed_counts, rnaseqkd_list$coldata, 0.05, 2)

# with(rnaseqkd_list$deseq2$ma_plots$shC1_c1509$data, plot(mean, lfc, pch=19, cex = 0.3, main="Ma plot", xlim=c(0,15), col="grey", ylim=c(-8,8)))
# #Upregulated
# with(subset(rnaseqkd_list$deseq2$ma_plots$shC1_c1509$data, padj<.05 & lfc > 1), points(mean, lfc, pch=19, cex = 0.4, xlim=c(0,15), ylim=c(-8,8),lwd = 0.4, col="black",bg="grey"))
# #Downregulated
# with(subset(rnaseqkd_list$deseq2$ma_plots$shC1_c1509$data, padj<.05 & lfc < -1), points(mean, lfc, pch=19, cex = 0.4, xlim=c(0,15), ylim=c(-8,8), lwd = 0.4,col="black", bg="grey"))

# Overlap between cramp1 and suz12 target genes
# Venn diagram
list_rnaseqkd_de_shC1_c1509_shS_c1509 <- list(shC1_c1509 = unique(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$de$feature_id), shS_c1509=unique(rnaseqkd_list$deseq2$de_analysis$shS_c1509$de$feature_id))
plot_venn(list_rnaseqkd_de_shC1_c1509_shS_c1509, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_rnaseqkd_up_shC1_c1509_shS_c1509 <- list(shC1_c1509 = unique(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$up$feature_id), shS_c1509=unique(rnaseqkd_list$deseq2$de_analysis$shS_c1509$up$feature_id))
plot_venn(list_rnaseqkd_up_shC1_c1509_shS_c1509, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_rnaseqkd_down_shC1_c1509_shS_c1509 <- list(shC1_c1509 = unique(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$down$feature_id), shS_c1509=unique(rnaseqkd_list$deseq2$de_analysis$shS_c1509$down$feature_id))
plot_venn(list_rnaseqkd_down_shC1_c1509_shS_c1509, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

# extract unique genes to each category
# take the rownames or ensID_Gene instead of direct gene symbol, I took rownames which is same as ensID_Gene
set_rnaseqde_uniq_shC1_up <- data.frame(ensID_Gene=setdiff(unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shC1_c1509), unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shS_c1509)), label="only_shC1_up")
set_rnaseqde_uniq_shS_up <- data.frame(ensID_Gene=setdiff(unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shS_c1509), unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shC1_c1509)), label="only_shS_up")
set_rnaseqde_both_shC1_shS_up <- data.frame(ensID_Gene=intersect(unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shC1_c1509), unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shS_c1509)), label="both_shC1_shS_up")
set_rnaseqde_none_up <- data.frame(ensID_Gene=setdiff(setdiff(rnaseqkd_list$gene$ens_gene, unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shC1_c1509)), unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shS_c1509)), label="none_up")
set_rnaseqde_merged_shC1_shS_up <- rbind.data.frame(set_rnaseqde_uniq_shC1_up, set_rnaseqde_uniq_shS_up, set_rnaseqde_both_shC1_shS_up, set_rnaseqde_none_up)

list_set_rnaseqde_merged_shC1_shS_up <- list(shC1_up=set_rnaseqde_uniq_shC1_up$ensID_Gene, shS_up=set_rnaseqde_uniq_shS_up$ensID_Gene, shC1_shS_up=set_rnaseqde_both_shC1_shS_up$ensID_Gene, none_up=set_rnaseqde_none_up$ensID_Gene)
ggvenn::ggvenn(
  list_set_rnaseqde_merged_shC1_shS_up, 
  fill_color = c("#0073C2FF", "#EFC000FF","#999999", "#FF9999", "#56B4E9"),
  stroke_size = 0.5, set_name_size = 4
)

set_rnaseqde_merged_shC1_shS_up_pos <- merge(set_rnaseqde_merged_shC1_shS_up, integration_features_list$gene, by.x="ensID_Gene", by.y="ens_gene")
write.table(set_rnaseqde_merged_shC1_shS_up_pos[,c(3:5,1,2,6)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/","set_rnaseqde_merged_shC1_shS_up_pos.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

# export files
write.xlsx(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$de, file = "rnaseqkd_list_deseq2_de_analysis_shC1_c1509_de.xlsx", colNames = TRUE, rowNames = TRUE)
write.xlsx(rnaseqkd_list$deseq2$de_analysis$shS_c1509$de, file = "rnaseqkd_list_deseq2_de_analysis_shS_c1509_de.xlsx", colNames = TRUE, rowNames = TRUE)

# go term analysis
# de
# list_rnaseqkd_de_shC1_c1509_shS_c1509
rnaseqkd_list[["goterm"]][["de"]] <- go_term_analysis(list_rnaseqkd_de_shC1_c1509_shS_c1509, "hsapiens")

# list_rnaseqkd_up_shC1_c1509_shS_c1509
rnaseqkd_list[["goterm"]][["up"]] <- go_term_analysis(list_rnaseqkd_up_shC1_c1509_shS_c1509, "hsapiens")

# list_rnaseqkd_down_shC1_c1509_shS_c1509
rnaseqkd_list[["goterm"]][["down"]] <- go_term_analysis(list_rnaseqkd_down_shC1_c1509_shS_c1509, "hsapiens")

# unique category
# remove none_up
list_set_rnaseqde_merged_shC1_shS_up_filt <- list_set_rnaseqde_merged_shC1_shS_up
list_set_rnaseqde_merged_shC1_shS_up_filt$none_up <- NULL
rnaseqkd_list[["goterm"]][["unique"]] <- go_term_analysis(list_set_rnaseqde_merged_shC1_shS_up_filt, "hsapiens")

# ma plot shrinkage
rnaseqkd_list[["deseq2"]][["maplots_shrinkage"]][["shC1_c1509"]] <- maplot_general_shrinkage(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$all, c(1,3,6,7), 0.05, 1)

################################################
###       ATAC-Seq data: Knockdown  PE       ###
################################################

# deseq2 with consensus peaks (direct from nextflow)
atacseqkd_deseq2_path <- "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/"

atacseqkd_deseq2_list <- list()

atacseqkd_deseq2_list[["dds"]] <- load_rdata(atacseqkd_deseq2_path, "boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/deseq2/consensus_peaks.mLb.clN.dds.RData")
atacseqkd_deseq2_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_upstream_2kb.sorted.txt")
atacseqkd_deseq2_list[["counts"]] <- dds_counts(atacseqkd_deseq2_list$dds)
atacseqkd_deseq2_list[["samples_selected"]] <- sample_filter(atacseqkd_deseq2_list$counts, samplefilter = FALSE)
atacseqkd_deseq2_list[["processed_counts"]] <- atacseqkd_deseq2_list[["samples_selected"]][,c(1:6)]
atacseqkd_deseq2_list[["coldata"]] <- make_coldata_dds(atacseqkd_deseq2_list$dds, samplefilter = FALSE)
atacseqkd_deseq2_list[["consensus"]] <- data.frame(fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.bed"))
# check
all(rownames(atacseqkd_deseq2_list$coldata) == colnames(atacseqkd_deseq2_list$processed_counts)) #should print TRUE

# de analysis
atacseqkd_deseq2_list[["dar_analysis"]] <- deseq2_de(atacseqkd_deseq2_list$processed_counts, atacseqkd_deseq2_list$coldata, 0.05, 2)

# get positions of dars
atacseqkd_deseq2_list[["annotation"]] <- consensus_gene_anno("atacseqkd", "deseq2", atacseqkd_deseq2_list$dar_analysis$de_analysis, atacseqkd_deseq2_list$consensus, atacseqkd_deseq2_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

atacseqkd_deseq2_list[["chipseeker"]][["annotation"]] <- chipseeker_region_anno("atacseqkd", "deseq2", atacseqkd_deseq2_list$dar_analysis$de_analysis, atacseqkd_deseq2_list$consensus[,c(1:4)], "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2) 

atacseqkd_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl[["maplot"]] <- maplot(atacseqkd_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$all, c(1,2,6,7), 0.05, 2)
atacseqkd_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl[["maplot"]] <- maplot(atacseqkd_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$all, c(1,2,6,7), 0.05, 2)

write.table(atacseqkd_deseq2_list$annotation$shCRAMP1_shControl$all$de[,c(8:10,3,7,1,11:12)] %>% distinct(), "atacseqkd_deseq2_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_deseq2_list$annotation$shSUZ12_shControl$all$de[,c(8:10,3,7,1,11:12)] %>% distinct(), "atacseqkd_deseq2_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_deseq2_de_shCRAMP1_shSUZ12_shControl <- list(shCRAMP1 = unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$de$feature_id), shSUZ12=unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$de$feature_id))
plot_venn(list_atacseqkd_deseq2_de_shCRAMP1_shSUZ12_shControl, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_atacseqkd_deseq2_up_shCRAMP1_shSUZ12_shControl <- list(shCRAMP1 = unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$up$feature_id), shSUZ12=unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$up$feature_id))
plot_venn(list_atacseqkd_deseq2_up_shCRAMP1_shSUZ12_shControl, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotPCA(atacseqkd_deseq2_list$dar_analysis$vstO, intgroup="sample", returnData=FALSE)


atacseqkd_deseq2_list[["aggregate_pergene"]][["shCRAMP1_shControl"]][["nearest"]][["de"]] <- aggregate_feature(atacseqkd_deseq2_list$annotation$shCRAMP1_shControl$nearest$de, c(3), "ens_gene", mean)
atacseqkd_deseq2_list[["aggregate_pergene"]][["shSUZ12_shControl"]][["nearest"]][["de"]] <- aggregate_feature(atacseqkd_deseq2_list$annotation$shSUZ12_shControl$nearest$de, c(3), "ens_gene", mean)

colnames(atacseqkd_deseq2_list$aggregate_pergene$shCRAMP1_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(atacseqkd_deseq2_list$aggregate_pergene$shSUZ12_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")

# edgeR with consensus peaks (direct from nextflow)
atacseqkd_edger_list <- list()

atacseqkd_edger_list[["dds"]] <- load_rdata(atacseqkd_deseq2_path, "boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/deseq2/consensus_peaks.mLb.clN.dds.RData")
atacseqkd_edger_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_upstream_2kb.sorted.txt")
atacseqkd_edger_list[["counts"]] <- dds_counts(atacseqkd_edger_list$dds)
atacseqkd_edger_list[["samples_selected"]] <- sample_filter(atacseqkd_edger_list$counts, samplefilter = FALSE)
atacseqkd_edger_list[["processed_counts"]] <- atacseqkd_edger_list[["samples_selected"]][,c(6,4,5,1,3,2)] # Specific requirement here
atacseqkd_edger_list[["coldata"]] <- make_coldata_dds(atacseqkd_edger_list$dds, samplefilter = FALSE)
colnames(atacseqkd_edger_list$coldata)[colnames(atacseqkd_edger_list$coldata) == "condition"] <- "Condition" # Specific requirement here
atacseqkd_edger_list$coldata$Condition <- factor(atacseqkd_edger_list$coldata$Condition)
atacseqkd_edger_list$coldata <- atacseqkd_edger_list$coldata[c(6,4,5,1,3,2),c(1:3)] # Specific requirement here

atacseqkd_edger_list[["dge_obj"]] <- edger_create(atacseqkd_edger_list$processed_counts, atacseqkd_edger_list$coldata)
atacseqkd_edger_list[["dar_analysis"]] <- edger_de(atacseqkd_edger_list$dge_obj$lmfit, atacseqkd_edger_list$dge_obj$design, 0.05, 2)
atacseqkd_edger_list[["consensus"]] <- data.frame(fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.bed"))
atacseqkd_edger_list[["annotation"]] <- consensus_gene_anno("atacseqkd", "edger", atacseqkd_edger_list$dar_analysis$dars, atacseqkd_edger_list$consensus, atacseqkd_edger_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)
# atacseqkd_edger_list[["chipseeker"]][["plots"]] <- run_chipseeker("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/", "atacseqkd", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "*.broadPeak")
atacseqkd_edger_list[["chipseeker"]][["annotation"]] <- chipseeker_region_anno("atacseqkd", "edger", atacseqkd_edger_list$dar_analysis$dars, atacseqkd_edger_list$consensus[,c(1:4)], "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2) 
atacseqkd_edger_list$dar_analysis$dars$coef_shCRAMP1[["maplot"]] <- maplot(atacseqkd_edger_list$dar_analysis$dars$coef_shCRAMP1$all, c(2,1,5,6), 0.05, 2)
atacseqkd_edger_list$dar_analysis$dars$coef_shSUZ12[["maplot"]] <- maplot(atacseqkd_edger_list$dar_analysis$dars$coef_shSUZ12$all, c(2,1,5,6), 0.05, 2)

# get positions of dars
write.table(atacseqkd_edger_list$annotation$coef_shCRAMP1$all$de[,c(7:9,2,6,1,10:11)] %>% distinct(), "atacseqkd_edger_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_edger_list$annotation$coef_shSUZ12$all$de[,c(7:9,2,6,1,10:11)] %>% distinct(), "atacseqkd_edger_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_edger_de_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_edger_list$dar_analysis$dars$coef_shCRAMP1$de$feature_id), shSUZ12=unique(atacseqkd_edger_list$dar_analysis$dars$coef_shSUZ12$de$feature_id))
plot_venn(list_atacseqkd_edger_de_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_atacseqkd_edger_up_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_edger_list$dar_analysis$dars$coef_shCRAMP1$up$feature_id), shSUZ12=unique(atacseqkd_edger_list$dar_analysis$dars$coef_shSUZ12$up$feature_id))
plot_venn(list_atacseqkd_edger_up_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotMDS(atacseqkd_edger_list$dge_obj$mds,top = 100000, pch = 16, cex = 1, dim.plot = c(1,2), col=c(rep("black",2), rep("red",2), rep("blue",2)))
plotMDS(atacseqkd_edger_list$dge_obj$mds)

atacseqkd_edger_list[["aggregate_pergene"]][["coef_shCRAMP1"]][["nearest"]][["de"]] <- aggregate_feature(atacseqkd_edger_list$annotation$coef_shCRAMP1$nearest$de, c(3), "ens_gene", mean)
atacseqkd_edger_list[["aggregate_pergene"]][["coef_shSUZ12"]][["nearest"]][["de"]] <- aggregate_feature(atacseqkd_edger_list$annotation$coef_shSUZ12$nearest$de, c(3), "ens_gene", mean)

colnames(atacseqkd_edger_list$aggregate_pergene$coef_shCRAMP1$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(atacseqkd_edger_list$aggregate_pergene$coef_shSUZ12$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")

# edgeR with consensus peaks (direct from nextflow)
atacseqkd_limma_list <- list()

atacseqkd_limma_list[["dds"]] <- load_rdata(atacseqkd_deseq2_path, "boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/deseq2/consensus_peaks.mLb.clN.dds.RData")
atacseqkd_limma_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_upstream_2kb.sorted.txt")
atacseqkd_limma_list[["counts"]] <- dds_counts(atacseqkd_limma_list$dds)
atacseqkd_limma_list[["samples_selected"]] <- sample_filter(atacseqkd_limma_list$counts, samplefilter = FALSE)
atacseqkd_limma_list[["processed_counts"]] <- atacseqkd_limma_list[["samples_selected"]][,c(6,4,5,1,3,2)] # Specific requirement here
atacseqkd_limma_list[["coldata"]] <- make_coldata_dds(atacseqkd_limma_list$dds, samplefilter = FALSE)
colnames(atacseqkd_limma_list$coldata)[colnames(atacseqkd_limma_list$coldata) == "condition"] <- "Condition" # Specific requirement here
atacseqkd_limma_list$coldata$Condition <- factor(atacseqkd_limma_list$coldata$Condition)
atacseqkd_limma_list$coldata <- atacseqkd_limma_list$coldata[c(6,4,5,1,3,2),c(1:3)] # Specific requirement here


atacseqkd_limma_list[["dge_obj"]] <- limma_create(atacseqkd_limma_list$processed_counts, atacseqkd_limma_list$coldata)

# create contrasts
atacseqkd_limma_list[["contrasts"]] <- limma::makeContrasts(paste0(colnames(atacseqkd_limma_list$dge_obj$design)[c(2,1)], collapse = "-"),
                                                            paste0(colnames(atacseqkd_limma_list$dge_obj$design)[c(3,1)], collapse = "-"), levels = atacseqkd_limma_list$dge_obj$design)

atacseqkd_limma_list[["dar_analysis"]] <- limma_de(atacseqkd_limma_list$dge_obj$fit, atacseqkd_limma_list$contrasts, 0.05, 2)

atacseqkd_limma_list[["consensus"]] <- data.frame(fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.bed"))
# scaffolds with non-chr annotation wil be filtered
atacseqkd_limma_list[["annotation"]] <- consensus_gene_anno("atacseqkd", "limma", atacseqkd_limma_list$dar_analysis$dars,atacseqkd_limma_list$consensus, atacseqkd_limma_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)
atacseqkd_limma_list[["chipseeker"]][["annotation"]] <- chipseeker_region_anno("atacseqkd", "limma", atacseqkd_limma_list$dar_analysis$dars, atacseqkd_limma_list$consensus[,c(1:4)], "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2) 

# Create MA plot
atacseqkd_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl[["maplot"]] <- maplot_general(atacseqkd_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl$all, c(2,1,5,7), 0.05, 1)
atacseqkd_limma_list$dar_analysis$dars$coef_shSUZ12_shControl[["maplot"]] <- maplot_general(atacseqkd_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$all, c(2,1,5,7), 0.05, 1)

# get positions of dars
write.table(atacseqkd_limma_list$annotation$coef_shCRAMP1$all$de[,c(8:10,2,6,1,11:12)] %>% distinct(), "atacseqkd_limma_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_limma_list$annotation$coef_shSUZ12$all$de[,c(8:10,2,6,1,11:12)] %>% distinct(), "atacseqkd_limma_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_limma_de_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl$de$feature_id), shSUZ12=unique(atacseqkd_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$de$feature_id))
plot_venn(list_atacseqkd_limma_de_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_atacseqkd_limma_up_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl$up$feature_id), shSUZ12=unique(atacseqkd_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$up$feature_id))
plot_venn(list_atacseqkd_limma_up_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotMDS(atacseqkd_limma_list$dge_obj$dgelist,top = 100000, pch = 16, cex = 1, dim.plot = c(1,2), col=c(rep("black",2), rep("red",2), rep("blue",2)))
plotMDS(atacseqkd_limma_list$dge_obj$dgelist)

atacseqkd_limma_list[["aggregate_pergene"]][["coef_shCRAMP1_shControl"]][["nearest"]][["de"]] <- aggregate_feature(atacseqkd_limma_list$annotation$coef_shCRAMP1_shControl$nearest$de, c(3), "ens_gene", mean)
atacseqkd_limma_list[["aggregate_pergene"]][["coef_shSUZ12_shControl"]][["nearest"]][["de"]] <- aggregate_feature(atacseqkd_limma_list$annotation$coef_shSUZ12_shControl$nearest$de, c(3), "ens_gene", mean)

colnames(atacseqkd_limma_list$aggregate_pergene$coef_shCRAMP1_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(atacseqkd_limma_list$aggregate_pergene$coef_shSUZ12_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")

#--------------------
atacseqkd_venn_consensus_list_CRAMP1 <- list(
  e_CRAMP1_up = unique(do.call(paste, c(atacseqkd_edger_list$annotation$coef_shCRAMP1$all$up[,c(7:9)], sep="%"))),
  e_CRAMP1_down = unique(do.call(paste, c(atacseqkd_edger_list$annotation$coef_shCRAMP1$all$down[,c(7:9)], sep="%"))),
  l_CRAMP1_up = unique(do.call(paste, c(atacseqkd_limma_list$annotation$coef_shCRAMP1_shControl$all$up[,c(8:10)], sep="%"))),
  l_CRAMP1_down = unique(do.call(paste, c(atacseqkd_limma_list$annotation$coef_shCRAMP1_shControl$all$down[,c(8:10)], sep="%"))),
  d_CRAMP1_up = unique(do.call(paste, c(atacseqkd_deseq2_list$annotation$shCRAMP1_shControl$all$up[,c(8:10)], sep="%"))),
  d_CRAMP1_down = unique(do.call(paste, c(atacseqkd_deseq2_list$annotation$shCRAMP1_shControl$all$down[,c(8:10)], sep="%")))
)
upset(fromList(atacseqkd_venn_consensus_list_CRAMP1), order.by = "freq", nsets = 6)
#--------------------
# diffbind
diffbind_atacseqkd_samplesheet <- read.csv("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/diffbind_atac_seq_samplesheet.csv")
atacseqkd_diffbind_list <- list()
atacseqkd_diffbind_list[["prep"]] <- diffbind_prep(diffbind_atacseqkd_samplesheet)
plot(atacseqkd_diffbind_list$prep$dba_obj)
plot(atacseqkd_diffbind_list$prep$ovp_rate,type ='b',ylab='# peaks', xlab='Overlap at least this many peaksets')



atacseqkd_diffbind_list[["counts"]] <- diffbind_count(atacseqkd_diffbind_list$prep$dba_obj, 100)
atacseqkd_diffbind_list[["norm"]] <-  diffbind_norm(atacseqkd_diffbind_list$counts$dba_obj)

# check for a match before and after dars analysis
atacseqkd_diffbind_list[["counts"]][["reads"]]  <- dba.count(atacseqkd_diffbind_list$counts$dba_obj, peaks=NULL, score=DBA_SCORE_READS)
atacseqkd_diffbind_list[["counts"]][["counts"]] <- dba.peakset(atacseqkd_diffbind_list$counts$reads, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
atacseqkd_diffbind_counts <- data.frame(atacseqkd_diffbind_list$counts$counts)
rownames(atacseqkd_diffbind_counts) <- paste0(atacseqkd_diffbind_counts$CHR, "%",atacseqkd_diffbind_counts$START,"%", atacseqkd_diffbind_counts$END)
atacseqkd_diffbind_counts <- atacseqkd_diffbind_counts[,-c(1:3)]
atacseqkd_diffbind_list[["regions"]] <- data.frame(atacseqkd_diffbind_list$counts$counts[,c(1:3)], interval=unique(do.call(paste, c(atacseqkd_diffbind_list$counts$counts[,1:3], sep="%"))))

atacseqkd_diffbind_list[["dar_analysis"]][["edgeR"]] <- diffbind_de(atacseqkd_diffbind_list$norm$norm, 0.05, 2, "edgeR")
atacseqkd_diffbind_list[["dar_analysis"]][["deseq2"]] <- diffbind_de(atacseqkd_diffbind_list$norm$norm, 0.05, 2, "deseq2")

atacseqkd_diffbind_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_upstream_2kb.sorted.txt")

atacseqkd_diffbind_list[["annotation"]][["deseq2"]] <- diffbind_gene_anno("atacseqkd", "diffbind",atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts, atacseqkd_diffbind_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)
atacseqkd_diffbind_list[["annotation"]][["edgeR"]] <- diffbind_gene_anno("atacseqkd", "diffbind",atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts, atacseqkd_diffbind_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

atacseqkd_diffbind_list[["chipseeker"]][["annotation"]][["deseq2"]] <- chipseeker_region_anno("atacseqkd", "diffbind",atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts, atacseqkd_diffbind_list$regions, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2) 
atacseqkd_diffbind_list[["chipseeker"]][["annotation"]][["edgeR"]] <- chipseeker_region_anno("atacseqkd", "diffbind",atacseqkd_diffbind_list$dar_analysis$edger$contrasts, atacseqkd_diffbind_list$regions, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2) 

# since we need only deseq2 results for downstream analysis, I perform mean analysis with deseq2
atacseqkd_diffbind_list[["aggregate_pergene"]][["deseq2"]][["contrast_shCRAMP1_shControl"]][["nearest"]][["de"]] <- aggregate_feature(atacseqkd_diffbind_list$annotation$deseq2$contrast_shCRAMP1_shControl$nearest$de, c(9), "ens_gene", mean)
atacseqkd_diffbind_list[["aggregate_pergene"]][["deseq2"]][["contrast_shSUZ12_shControl"]][["nearest"]][["de"]] <- aggregate_feature(atacseqkd_diffbind_list$annotation$deseq2$contrast_shSUZ12_shControl$nearest$de, c(9), "ens_gene", mean)

colnames(atacseqkd_diffbind_list$aggregate_pergene$deseq2$contrast_shCRAMP1_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(atacseqkd_diffbind_list$aggregate_pergene$deseq2$contrast_shSUZ12_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")

# Create MA plot
atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shCRAMP1_shControl[["maplot"]] <- maplot(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shCRAMP1_shControl$all, c(6,9,11,12), 0.05, 2)
atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shSUZ12_shControl[["maplot"]] <- maplot(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shSUZ12_shControl$all, c(6,9,11,12), 0.05, 2)
atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shControl[["maplot"]] <- maplot(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shControl$all, c(6,9,11,12), 0.05, 2)
atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl[["maplot"]] <- maplot(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl$all, c(6,9,11,12), 0.05, 2)

# e= edger and d= deseq2
write.table(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shCRAMP1_shControl$de[,c(1:3,9,11:12)] %>% distinct(), "atacseqkd_diffbind_e_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shSUZ12_shControl$de[,c(1:3,9,11:12)] %>% distinct(), "atacseqkd_diffbind_e_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shControl$de[,c(1:3,9,11:12)] %>% distinct(), "atacseqkd_diffbind_d_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl$de[,c(1:3,9,11:12)] %>% distinct(), "atacseqkd_diffbind_d_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# extract dars in specified format
write.table(cbind(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shControl$de[,c(1:3,12)] %>% distinct(), dars="dars", strand="."), "atacseqkd_diffbind_d_shCRAMP1_shControl_de_re.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(cbind(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl$de[,c(1:3,12)] %>% distinct(), dars="dars", strand="."), "atacseqkd_diffbind_d_shSUZ12_shControl_de_re.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# up
write.table(cbind(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shControl$up[,c(1:3,12)] %>% distinct(), dars="dars", strand="."), "atacseqkd_diffbind_d_shCRAMP1_shControl_up_re.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(cbind(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl$up[,c(1:3,12)] %>% distinct(), dars="dars", strand="."), "atacseqkd_diffbind_d_shSUZ12_shControl_up_re.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# down
write.table(cbind(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shControl$down[,c(1:3,12)] %>% distinct(), dars="dars", strand="."), "atacseqkd_diffbind_d_shCRAMP1_shControl_down_re.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(cbind(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl$down[,c(1:3,12)] %>% distinct(), dars="dars", strand="."), "atacseqkd_diffbind_d_shSUZ12_shControl_down_re.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_diffbind_e_de_shCRAMP1_vs_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shCRAMP1_shContro$de$feature_id), shSUZ12=unique(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shSUZ12_shControl$de$feature_id))
plot_venn(list_atacseqkd_diffbind_e_de_shCRAMP1_vs_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))
list_atacseqkd_diffbind_d_de_shCRAMP1_vs_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shContro$de$feature_id), shSUZ12=unique(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl$de$feature_id))
plot_venn(list_atacseqkd_diffbind_d_de_shCRAMP1_vs_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

#PCA
dba.plotPCA(atacseqkd_diffbind_list$norm$norm, method=DBA_ALL_METHODS, attributes=DBA_FACTOR, label=DBA_ID)



# atacseqkd_diffbind_list[["binding_mat"]] <- data.frame(atacseqkd_diffbind_list$counts$dba_obj$binding[,c(1:3)])
# atacseqkd_diffbind_list$binding_mat$CHR <- paste0("chr",atacseqkd_diffbind_list$binding_mat$CHR)

length(intersect(
do.call(paste, c(atacseqkd_diffbind_list$counts$counts[,c(1:3)], sep="%")),
do.call(paste, c(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shCRAMP1_shControl$all[,c(1:3)], sep="%"))
))

# length(intersect(
#   do.call(paste, c(atacseqkd_diffbind_list$counts$counts[,c(1:3)], sep="%")),
#   do.call(paste, c(atacseqkd_diffbind_list$binding_mat, sep="%"))
# ))

# diffbind_deseq2 
atacseqkd_diffbind_deseq2_list <- list()

atacseqkd_diffbind_deseq2_list[["counts"]] <- atacseqkd_diffbind_counts
atacseqkd_diffbind_deseq2_list[["coldata"]] <- atacseqkd_diffbind_list$counts$info
rownames(atacseqkd_diffbind_deseq2_list$coldata) <- atacseqkd_diffbind_deseq2_list$coldata$ID
atacseqkd_diffbind_deseq2_list$coldata$condition <- factor(atacseqkd_diffbind_deseq2_list$coldata$Condition)
atacseqkd_diffbind_deseq2_list$coldata <- atacseqkd_diffbind_deseq2_list$coldata[,c(1,8,3)]
all(rownames(atacseqkd_diffbind_deseq2_list$coldata) == colnames(atacseqkd_diffbind_deseq2_list$counts)) #should print TRUE
# de analysis
atacseqkd_diffbind_deseq2_list[["dar_analysis"]] <- deseq2_de(atacseqkd_diffbind_deseq2_list$counts, atacseqkd_diffbind_deseq2_list$coldata, 0.05, 2)

atacseqkd_diffbind_deseq2_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_upstream_2kb.sorted.txt")
# get positions of dars
atacseqkd_diffbind_deseq2_list[["annotation"]] <- diffbind_and_packages_gene_anno("atacseqkd", "diffbind_deseq2", atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis, atacseqkd_diffbind_deseq2_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl[["maplot"]] <- maplot(atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$all, c(1,2,6,7), 0.05, 2)
atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl[["maplot"]] <- maplot(atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$all, c(1,2,6,7), 0.05, 2)

write.table(atacseqkd_diffbind_deseq2_list$annotation$shCRAMP1_shControl$all$de[,c(1:3,5,9)] %>% distinct(), "atacseqkd_diffbind_deseq2_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_deseq2_list$annotation$shSUZ12_shControl$all$de[,c(1:3,5,9)] %>% distinct(), "atacseqkd_diffbind_deseq2_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_diffbind_deseq2_de_shCRAMP1_shSUZ12_shControl <- list(shCRAMP1 = unique(atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$de$feature_id), shSUZ12=unique(atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$de$feature_id))
plot_venn(list_atacseqkd_diffbind_deseq2_de_shCRAMP1_shSUZ12_shControl, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotPCA(atacseqkd_diffbind_deseq2_list$dar_analysis$vstO, intgroup="ID", returnData=FALSE)


# DE deseq2 manually from diffbind mean the logFC per peaks to specify genes
atacseqkd_diffbind_deseq2_list[["aggregate_pergene"]][["deseq2"]][["shCRAMP1_shControl"]][["nearest"]][["de"]] <- aggregate_feature(atacseqkd_diffbind_deseq2_list$annotation$shCRAMP1_shControl$nearest$de, c(5), "ens_gene", mean)
colnames(atacseqkd_diffbind_deseq2_list$aggregate_pergene$deseq2$shCRAMP1_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")

atacseqkd_diffbind_deseq2_list[["aggregate_pergene"]][["deseq2"]][["shSUZ12_shControl"]][["nearest"]][["de"]] <- aggregate_feature(atacseqkd_diffbind_deseq2_list$annotation$shSUZ12_shControl$nearest$de, c(5), "ens_gene", mean)
colnames(atacseqkd_diffbind_deseq2_list$aggregate_pergene$deseq2$shSUZ12_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")


# edgeR analysis
atacseqkd_diffbind_edger_list <- list()

atacseqkd_diffbind_edger_list[["counts"]] <- atacseqkd_diffbind_counts
atacseqkd_diffbind_edger_list[["coldata"]] <- atacseqkd_diffbind_list$counts$info
atacseqkd_diffbind_edger_list$coldata$Condition <- factor(atacseqkd_diffbind_edger_list$coldata$Condition)
atacseqkd_diffbind_edger_list[["dge_obj"]] <- edger_create(atacseqkd_diffbind_edger_list$counts, atacseqkd_diffbind_edger_list$coldata)

atacseqkd_diffbind_edger_list[["dar_analysis"]] <- edger_de(atacseqkd_diffbind_edger_list$dge_obj$lmfit, atacseqkd_diffbind_edger_list$dge_obj$design, 0.05, 2)

atacseqkd_diffbind_edger_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_upstream_2kb.sorted.txt")

atacseqkd_diffbind_edger_list[["annotation"]] <- diffbind_and_packages_gene_anno("atacseqkd", "diffbind_edger", atacseqkd_diffbind_edger_list$dar_analysis$dars, atacseqkd_diffbind_edger_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shCRAMP1[["maplot"]] <- maplot(atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shCRAMP1$all, c(2,1,5,6), 0.05, 2)
atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shSUZ12[["maplot"]] <- maplot(atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shSUZ12$all, c(2,1,5,6), 0.05, 2)

# get positions of dars
write.table(atacseqkd_diffbind_edger_list$annotation$coef_shCRAMP1$all$de[,c(1:4,8)] %>% distinct(), "atacseqkd_diffbind_edger_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_edger_list$annotation$coef_shSUZ12$all$de[,c(1:4,8)] %>% distinct(), "atacseqkd_diffbind_edger_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_diffbind_edger_de_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shCRAMP1$de$feature_id), shSUZ12=unique(atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shSUZ12$de$feature_id))
plot_venn(list_atacseqkd_diffbind_edger_de_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotMDS(atacseqkd_diffbind_edger_list$dge_obj$mds,top = 100000, pch = 16, cex = 1, dim.plot = c(1,2), col=c(rep("black",2), rep("red",2), rep("blue",2)))
plotMDS(atacseqkd_diffbind_edger_list$dge_obj$mds)

# limma analysis
atacseqkd_diffbind_limma_list <- list()

atacseqkd_diffbind_limma_list[["counts"]] <- atacseqkd_diffbind_counts
atacseqkd_diffbind_limma_list[["coldata"]] <- atacseqkd_diffbind_list$counts$info
atacseqkd_diffbind_limma_list$coldata$Condition <- factor(atacseqkd_diffbind_limma_list$coldata$Condition)
atacseqkd_diffbind_limma_list[["dge_obj"]] <- limma_create(atacseqkd_diffbind_limma_list$counts, atacseqkd_diffbind_limma_list$coldata)

# create contrasts
atacseqkd_diffbind_limma_list[["contrasts"]] <- limma::makeContrasts(paste0(colnames(atacseqkd_diffbind_limma_list$dge_obj$design)[c(2,1)], collapse = "-"),
                                                                     paste0(colnames(atacseqkd_diffbind_limma_list$dge_obj$design)[c(3,1)], collapse = "-"), levels = atacseqkd_diffbind_limma_list$dge_obj$design)

atacseqkd_diffbind_limma_list[["dar_analysis"]] <- limma_de(atacseqkd_diffbind_limma_list$dge_obj$fit, atacseqkd_diffbind_limma_list$contrasts, 0.05, 2)

atacseqkd_diffbind_limma_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_upstream_2kb.sorted.txt")

atacseqkd_diffbind_limma_list[["annotation"]] <- diffbind_and_packages_gene_anno("atacseqkd", "diffbind_limma", atacseqkd_diffbind_limma_list$dar_analysis$dars, atacseqkd_diffbind_limma_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

# Create MA plot
atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl[["maplot"]] <- maplot(atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl$all, c(2,1,5,7), 0.05, 2)
atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shSUZ12_shControl[["maplot"]] <- maplot(atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$all, c(2,1,5,7), 0.05, 2)

# get positions of dars
write.table(atacseqkd_diffbind_limma_list$annotation$coef_shCRAMP1$all$de[,c(1:4,8)] %>% distinct(), "atacseqkd_diffbind_limma_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_limma_list$annotation$coef_shSUZ12$all$de[,c(1:4,8)] %>% distinct(), "atacseqkd_diffbind_limma_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_diffbind_limma_de_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl$de$feature_id), shSUZ12=unique(atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$de$feature_id))
plot_venn(list_atacseqkd_diffbind_limma_de_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotMDS(atacseqkd_diffbind_limma_list$dge_obj$dgelist,top = 100000, pch = 16, cex = 1, dim.plot = c(1,2), col=c(rep("black",2), rep("red",2), rep("blue",2)))
plotMDS(atacseqkd_diffbind_limma_list$dge_obj$dgelist)

#--------------------
atacseqkd_diffbind_venn_packages_list_CRAMP1 <- list(
  e_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shCRAMP1$all$up[,c(1:3)], sep="%"))),
  e_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shCRAMP1$all$down[,c(1:3)], sep="%"))),
  l_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shCRAMP1_shControl$all$up[,c(1:3)], sep="%"))),
  l_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shCRAMP1_shControl$all$down[,c(1:3)], sep="%"))),
  d_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shCRAMP1_shControl$all$up[,c(1:3)], sep="%"))),
  d_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shCRAMP1_shControl$all$down[,c(1:3)], sep="%")))
)
upset(fromList(atacseqkd_diffbind_venn_packages_list_CRAMP1), order.by = "freq", nsets = 6)
#--------------------

# Since diffbind regions are shared I can compare softwares up and down intervals
### for diffbind edgeR contrast control - sample so take down i.e they are upregulated in CRAMP1
### for limma edgeR contrast sample - control so take down i.e they are upregulated in CRAMP1 (atacseqkd_limma_list$contrasts)
# Note: I labelled down as up and up as down for diffbind outputs
#--------------------
atacseqkd_diffbind_with_packages_list_CRAMP1 <- list(
  e_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shCRAMP1$all$up[,c(1:3)], sep="%"))),
  e_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shCRAMP1$all$down[,c(1:3)], sep="%"))),
  l_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shCRAMP1_shControl$all$up[,c(1:3)], sep="%"))),
  l_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shCRAMP1_shControl$all$down[,c(1:3)], sep="%"))),
  d_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shCRAMP1_shControl$all$up[,c(1:3)], sep="%"))),
  d_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shCRAMP1_shControl$all$down[,c(1:3)], sep="%"))),
  db_e_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$edgeR$contrast_shCRAMP1_shControl$all$up[,c(1:3)], sep="%"))),
  db_e_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$edgeR$contrast_shCRAMP1_shControl$all$down[,c(1:3)], sep="%"))),
  db_d_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$deseq2$contrast_shCRAMP1_shControl$all$up[,c(1:3)], sep="%"))),
  db_d_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$deseq2$contrast_shCRAMP1_shControl$all$down[,c(1:3)], sep="%")))
)
upset(fromList(atacseqkd_diffbind_with_packages_list_CRAMP1), order.by = "freq", nsets = 10)
#--------------------
t(do.call(cbind.data.frame,lapply(atacseqkd_diffbind_with_packages_list_CRAMP1, function(x) length(x))))


# similar analysis for SUZ12
#--------------------
atacseqkd_venn_consensus_list_SUZ12 <- list(
  e_SUZ12_up = unique(do.call(paste, c(atacseqkd_edger_list$annotation$coef_shSUZ12$all$up[,c(7:9)], sep="%"))),
  e_SUZ12_down = unique(do.call(paste, c(atacseqkd_edger_list$annotation$coef_shSUZ12$all$down[,c(7:9)], sep="%"))),
  l_SUZ12_up = unique(do.call(paste, c(atacseqkd_limma_list$annotation$coef_shSUZ12_shControl$all$up[,c(8:10)], sep="%"))),
  l_SUZ12_down = unique(do.call(paste, c(atacseqkd_limma_list$annotation$coef_shSUZ12_shControl$all$down[,c(8:10)], sep="%"))),
  d_SUZ12_up = unique(do.call(paste, c(atacseqkd_deseq2_list$annotation$shSUZ12_shControl$all$up[,c(8:10)], sep="%"))),
  d_SUZ12_down = unique(do.call(paste, c(atacseqkd_deseq2_list$annotation$shSUZ12_shControl$all$down[,c(8:10)], sep="%")))
)
upset(fromList(atacseqkd_venn_consensus_list_SUZ12), order.by = "freq", nsets = 6)
#--------------------

#--------------------
atacseqkd_diffbind_venn_packages_list_SUZ12 <- list(
  e_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shSUZ12$all$up[,c(1:3)], sep="%"))),
  e_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shSUZ12$all$down[,c(1:3)], sep="%"))),
  l_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shSUZ12_shControl$all$up[,c(1:3)], sep="%"))),
  l_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shSUZ12_shControl$all$down[,c(1:3)], sep="%"))),
  d_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shSUZ12_shControl$all$up[,c(1:3)], sep="%"))),
  d_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shSUZ12_shControl$all$down[,c(1:3)], sep="%")))
)
upset(fromList(atacseqkd_diffbind_venn_packages_list_SUZ12), order.by = "freq", nsets = 6)
#--------------------

#--------------------
atacseqkd_diffbind_with_packages_list_SUZ12 <- list(
  e_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shSUZ12$all$up[,c(1:3)], sep="%"))),
  e_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shSUZ12$all$down[,c(1:3)], sep="%"))),
  l_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shSUZ12_shControl$all$up[,c(1:3)], sep="%"))),
  l_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shSUZ12_shControl$all$down[,c(1:3)], sep="%"))),
  d_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shSUZ12_shControl$all$up[,c(1:3)], sep="%"))),
  d_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shSUZ12_shControl$all$down[,c(1:3)], sep="%"))),
  db_e_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$edgeR$contrast_shSUZ12_shControl$all$up[,c(1:3)], sep="%"))),
  db_e_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$edgeR$contrast_shSUZ12_shControl$all$down[,c(1:3)], sep="%"))),
  db_d_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$deseq2$contrast_shSUZ12_shControl$all$up[,c(1:3)], sep="%"))),
  db_d_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$deseq2$contrast_shSUZ12_shControl$all$down[,c(1:3)], sep="%")))
)
upset(fromList(atacseqkd_diffbind_with_packages_list_SUZ12), order.by = "freq", nsets = 10)
#--------------------
t(do.call(cbind.data.frame,lapply(atacseqkd_diffbind_with_packages_list_SUZ12, function(x) length(x))))

#-------------------------------------------


#--------------------
atacseqkd_diffbind_with_packages_list_CRAMP1_SUZ12 <- list(
  e_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shSUZ12$all$up[,c(1:3)], sep="%"))),
  e_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shSUZ12$all$down[,c(1:3)], sep="%"))),
  l_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shSUZ12_shControl$all$up[,c(1:3)], sep="%"))),
  l_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shSUZ12_shControl$all$down[,c(1:3)], sep="%"))),
  d_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shSUZ12_shControl$all$up[,c(1:3)], sep="%"))),
  d_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shSUZ12_shControl$all$down[,c(1:3)], sep="%"))),
  db_e_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$edgeR$contrast_shSUZ12_shControl$all$up[,c(1:3)], sep="%"))),
  db_e_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$edgeR$contrast_shSUZ12_shControl$all$down[,c(1:3)], sep="%"))),
  db_d_SUZ12_down = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$deseq2$contrast_shSUZ12_shControl$all$up[,c(1:3)], sep="%"))),
  db_d_SUZ12_up = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$deseq2$contrast_shSUZ12_shControl$all$down[,c(1:3)], sep="%"))),
  e_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shCRAMP1$all$up[,c(1:3)], sep="%"))),
  e_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_edger_list$annotation$coef_shCRAMP1$all$down[,c(1:3)], sep="%"))),
  l_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shCRAMP1_shControl$all$up[,c(1:3)], sep="%"))),
  l_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_limma_list$annotation$coef_shCRAMP1_shControl$all$down[,c(1:3)], sep="%"))),
  d_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shCRAMP1_shControl$all$up[,c(1:3)], sep="%"))),
  d_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_deseq2_list$annotation$shCRAMP1_shControl$all$down[,c(1:3)], sep="%"))),
  db_e_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$edgeR$contrast_shCRAMP1_shControl$all$up[,c(1:3)], sep="%"))),
  db_e_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$edgeR$contrast_shCRAMP1_shControl$all$down[,c(1:3)], sep="%"))),
  db_d_CRAMP1_down = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$deseq2$contrast_shCRAMP1_shControl$all$up[,c(1:3)], sep="%"))),
  db_d_CRAMP1_up = unique(do.call(paste, c(atacseqkd_diffbind_list$annotation$deseq2$contrast_shCRAMP1_shControl$all$down[,c(1:3)], sep="%")))
)

upset(fromList(atacseqkd_diffbind_with_packages_list_CRAMP1_SUZ12), order.by = "freq", nsets = 20)
#--------------------
#-------------------------------------------
atacseqkd_venn_packages_list_CRAMP1_upsetobj <- Upsetout(atacseqkd_venn_packages_list_CRAMP1)
names(atacseqkd_venn_packages_list_CRAMP1_upsetobj) <-
   gsub("_0","",gsub("0_","", names(atacseqkd_venn_packages_list_CRAMP1_upsetobj)))

atacseqkd_venn_packages_list_CRAMP1_df <- NULL
for(i in c("e_CRAMP1_up", "l_CRAMP1_up", "e_CRAMP1_down", "l_CRAMP1_down", "d_CRAMP1_up")){
  print(i)
  subset <- intersect(atacseqkd_venn_packages_list_CRAMP1[[i]], rownames(atacseqkd_venn_packages_list_CRAMP1_upsetobj[[i]]))
  df <- data.frame(interval=subset, category=i)
  atacseqkd_venn_packages_list_CRAMP1_df <- rbind.data.frame(atacseqkd_venn_packages_list_CRAMP1_df, df)
}

# quantify over various features
atacseqkd_quantify_list <- list()
atacseqkd_quantify_list[["featurecounts"]][["featureslist"]]=list(gene_ext=c(".chr_1.5kb.txt", "gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), bins5kb=c("hg38_5kb.txt", "hg38_5kb.txt"), bins10kb=c("hg38_10kb.txt", "hg38_10kb.txt"), promoter=c("_gencodev41_promoter.txt", "gene_gencodev41_promoter.txt"))
atacseqkd_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_groups"]] <-  c("_rkd_gene_groups.txt", "selected_rkd_gene_groups.txt")
atacseqkd_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_uniq_diff"]] <-  c("_shC1_shS_up_pos.txt", "set_rnaseqde_merged_shC1_shS_up_pos.txt")
atacseqkd_quantify_list[["featurecounts"]][["pairs"]] <- list(shCRAMP1_shControl_1 = c("shCRAMP1_REP1.mLb.clN.sorted.bam","shControl_REP1.mLb.clN.sorted.bam"), shCRAMP1_shControl_2 = c("shCRAMP1_REP2.mLb.clN.sorted.bam", "shControl_REP2.mLb.clN.sorted.bam") , shSUZ12_shControl_1 = c("shSUZ12_REP1.mLb.clN.sorted.bam", "shControl_REP1.mLb.clN.sorted.bam"), shSUZ12_shControl_2 = c("shSUZ12_REP2.mLb.clN.sorted.bam", "shControl_REP2.mLb.clN.sorted.bam"))
atacseqkd_quantify_list[["featurecounts"]][["featuresmatrix"]] <- quantify_featurecounts_multifeature("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "atacseqkd", ".mLb.clN.sorted.bam", atacseqkd_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "shControl",pe=TRUE, diff="multi", atacseqkd_quantify_list$featurecounts$pairs)
atacseqkd_quantify_list[["featurecounts"]][["featuresmatrix_plots"]] <- meta_plot(atacseqkd_quantify_list$featurecounts$featuresmatrix, c(7:10), rep=1, splitside=3)

colnames(atacseqkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$counts) <- gsub(".mLb.clN.sorted.bam", "", colnames(atacseqkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$counts))

DESeq2::plotPCA(edgeR::cpm(atacseqkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$counts), label=TRUE)

ggplot(atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")

atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st[["merge_rep"]] <- sapply(strsplit(as.character(atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st$col), "_"), function(x) paste0(x[1:2], collapse = "_"))
atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st[["merge_rep_category"]] <- paste0(sapply(strsplit(as.character(atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st$col), "_"), function(x) paste0(x[1:2], collapse = "_")),"_", 
                                                                                                                 as.character(atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st$id))

# remove none_up category
atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub <- dplyr::filter(atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, !grepl('none_up', row))

ggplot(atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub, aes(x=merge_rep, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")

ggplot(atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub, aes(x=merge_rep, y=value, col =id)) + 
  geom_violin(aes(fill = id), color="black", trim = FALSE, position = position_dodge(0.9), size=0.25, width = 1) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9))  + theme_classic() + 
  geom_hline(yintercept = 0, linetype="dotted")+ labs(title="ATAC-Seq KD", x ="Samples", y = "Sample-Control (normalized values)") 

pairwise.wilcox.test(atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub$value, atacseqkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub$merge_rep_category,
                     p.adjust.method = "BH")


# histone marks
# quantify bam files in a given histone mark peak coordinates
atacseqkd_quantify_list[["featurecounts"]][["histonemarks"]] <- quantify_featurecounts("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/", "atacseqkd","\\.bam$", "\\_peaks_id.bed$","histone_marks", c(12,4,1:3,7), pairedend=TRUE, "hg38", "_",merge_sites_files=FALSE)

# quantify total reads in a bam file, used Rsamtools
atacseqkd_quantify_list[["bams"]] <- quantify_bams("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/", "atacseqkd","\\.bam$")
atacseqkd_labelpath = "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/histone_marks_label.txt"
atacseqkd_quantify_list[["featurecounts"]][["histonemarks"]][["summed"]] <- summedreads_per_feature(atacseqkd_quantify_list$featurecounts$histonemarks$histone_marks_countmatrix, atacseqkd_quantify_list$bams$total_reads, atacseqkd_labelpath)

# plot heatmap 
atacseqkd_quantify_list$featurecounts$histonemarks$summed$all_features_ratio$sample <- gsub(".mLb.clN.sorted.bam","",atacseqkd_quantify_list$featurecounts$histonemarks$summed$all_features_ratio$sample)
atacseqkd_quantify_list$featurecounts$histonemarks$summed[["susbet_all_features_ratio"]] <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$all_features_ratio[,c(2,3,7)]
atacseqkd_quantify_list$featurecounts$histonemarks$summed[["table_afr"]] <- data.frame(tidyr::pivot_wider(atacseqkd_quantify_list$featurecounts$histonemarks$summed$susbet_all_features_ratio, names_from = sample, values_from = Ratio))
rownames(atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr) <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr$label
atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr[,-1]
write.table(atacseqkd_quantify_list$featurecounts$histonemarks$summed$all_features, "atacseqkd_histonemarks_all_features.txt", sep="\t", quote = F, append = F, row.names = F, col.names = T)

# subtraction of respective control
atacseqkd_quantify_list$featurecounts$histonemarks$summed[["diff_table_afr"]] <- NULL
for (i in c("shCRAMP1_REP1", "shCRAMP1_REP2","shSUZ12_REP1","shSUZ12_REP2")){
  if(i == "shCRAMP1_REP1"){
    atacseqkd_quantify_list$featurecounts$histonemarks$summed[["diff_table_afr"]][[i]] <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr[grepl(i, colnames(atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr))] - atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr[["shControl_REP1"]]
  }else if(i == "shCRAMP1_REP2"){
    atacseqkd_quantify_list$featurecounts$histonemarks$summed[["diff_table_afr"]][[i]] <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr[grepl(i, colnames(atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr))] - atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr[["shControl_REP2"]]
  }else if(i == "shSUZ12_REP1"){
    atacseqkd_quantify_list$featurecounts$histonemarks$summed[["diff_table_afr"]][[i]] <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr[grepl(i, colnames(atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr))] - atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr[["shControl_REP1"]]
  }else if(i == "shSUZ12_REP2"){
    atacseqkd_quantify_list$featurecounts$histonemarks$summed[["diff_table_afr"]][[i]] <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr[grepl(i, colnames(atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr))] - atacseqkd_quantify_list$featurecounts$histonemarks$summed$table_afr[["shControl_REP2"]]
  }
}
atacseqkd_quantify_list$featurecounts$histonemarks$summed[["merged_diff_table_afr"]]<- do.call(cbind.data.frame, atacseqkd_quantify_list$featurecounts$histonemarks$summed$diff_table_afr)
atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr <- t(atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr)
atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr[order(rownames(atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr)),]
breaksListz <- seq(0, 0.01, by = 0.001)

pheatmap::pheatmap(as.matrix(atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr), na_col = "grey",breaks = breaksListz,
                   color = colorRampPalette(c("white", "orange"))(length(breaksListz)),
                   clustering_distance_cols = "euclidean", cluster_rows = T, cluster_cols = T, clustering_method = "ward.D")

# mean the histone marks
atacseqkd_quantify_list$featurecounts$histonemarks$summed[["merged_diff_table_afr_st"]] <- data.frame(stack(as.matrix(do.call(cbind.data.frame, atacseqkd_quantify_list$featurecounts$histonemarks$summed$diff_table_afr))))
atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr_st[["id"]] <- paste0(sapply(strsplit(as.character(atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr_st$row), "_"),function(x) x[1]), "%", as.character(atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr_st$col))
atacseqkd_quantify_list$featurecounts$histonemarks$summed[["aggregated"]] <- aggregate_feature(atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr_st, c(3), "id", mean)
atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated[["histonemarks"]] <- sapply(strsplit(as.character(atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated$Group.1), "%"),function(x) x[1])
atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated[["sample"]] <- sapply(strsplit(as.character(atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated$Group.1), "%"),function(x) x[2])
atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated[,-1]
atacseqkd_quantify_list$featurecounts$histonemarks$summed[["histone_mark_df"]] <- data.frame(tidyr::pivot_wider(atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated, names_from = sample, values_from = x))
rownames(atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df) <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df$histonemarks
atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df[,-1]
atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df["shCRAMP1"] <- base::rowMeans(atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df[,1:2])
atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df["shSUZ12"] <- base::rowMeans(atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df[,3:4])
atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df["shCRAMP1_median"] <- matrixStats::rowMedians(as.matrix(atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df[,1:2]))
atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df["shSUZ12_median"] <- matrixStats::rowMedians(as.matrix(atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df[,3:4]))
pheatmap::pheatmap(as.matrix(t(atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df[,c(5:6)])), na_col = "grey",breaks = breaksListz,
                   color = colorRampPalette(c("white", "darkred"))(length(breaksListz)),
                   clustering_distance_cols = "euclidean", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D")


# All samples
atacseqkd_bams_path <- "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/"
atacseqkd_correlation_heatmap_divlog2 <- fread(paste0(atacseqkd_bams_path, "atacseqkd_correlation_heatmap_divlog2.txt"), header = TRUE)
colnames(atacseqkd_correlation_heatmap_divlog2) <- gsub("'","",colnames(atacseqkd_correlation_heatmap_divlog2))
atacseqkd_correlation_heatmap_divlog2$V1 <- gsub("'","",atacseqkd_correlation_heatmap_divlog2$V1)
atacseqkd_correlation_heatmap_divlog2 <- data.frame(atacseqkd_correlation_heatmap_divlog2)
rownames(atacseqkd_correlation_heatmap_divlog2) <- atacseqkd_correlation_heatmap_divlog2$V1
atacseqkd_correlation_heatmap_divlog2 <- atacseqkd_correlation_heatmap_divlog2[,-1]
breaksList_cor = seq(-1, 1, by = 0.01)
pheatmap(as.matrix(atacseqkd_correlation_heatmap_divlog2), display_numbers = TRUE,  number_color = "black", breaks = breaksList_cor, legend_breaks = c(-1, 0, 1), fontsize_number = 9, border_color = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList_cor)), clustering_distance_cols = "euclidean", clustering_method = "ward.D")


# Annotation peaks
atacseqkd_quantify_list[["chipseeker"]][["plots"]] <- run_chipseeker("/mnt/beegfs6/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "atacseqkd", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "*_re.bed")

atacseqkd_quantify_list[["chipseeker"]][["plots"]][["pie_annodf"]] <- replot_pie_chipseeker(atacseqkd_quantify_list$chipseeker$plots$gr_anno)

atacseqkd_quantify_list$chipseeker$plots$pie_annodf$df_st$sample <- gsub("atacseqkd_diffbind_d_|_re.bed","",atacseqkd_quantify_list$chipseeker$plots$pie_annodf$df_st$sample)

ggplot(atacseqkd_quantify_list$chipseeker$plots$pie_annodf$df_st, aes(x="Feature", y=Frequency, fill=Feature))+ facet_wrap( ~ sample, ncol=2, nrow=3) +
  geom_bar(width = 1, stat = "identity", color="white") + coord_polar("y", start=0) + theme_void() + scale_fill_manual(values=c("#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251", "#B565A7", "#009B77", "#DD4124", "#D65076")) 

# RM told me to further simplify the pie chart to promoter, gene body and intergenic regions
category_chipseeker_df <- data.frame(id= c("Promoter" , "Gene_body", "Gene_body", "Gene_body", "Gene_body", "Gene_body", "Gene_body", "Downstream","Intergenic"), Feature=c(unique(as.character(atacseqkd_quantify_list$chipseeker$plots$pie_annodf$df_st$Feature))))

atacseqkd_quantify_list[["chipseeker"]][["plots"]][["id_pie_annodf"]] <- merge(atacseqkd_quantify_list$chipseeker$plots$pie_annodf$df_st, category_chipseeker_df, by="Feature")

atacseqkd_quantify_list[["chipseeker"]][["plots"]][["id_pie_annodf"]]["id_to_agg"] <- do.call(paste, c(atacseqkd_quantify_list[["chipseeker"]][["plots"]][["id_pie_annodf"]][,c(3:4)], sep="%"))

atacseqkd_quantify_list[["chipseeker"]][["plots"]][["aggregated"]] <- aggregate_feature(atacseqkd_quantify_list$chipseeker$plots$id_pie_annodf, c(2), "id_to_agg", sum)

atacseqkd_quantify_list$chipseeker$plots$aggregated$aggregated <- cSplit(atacseqkd_quantify_list$chipseeker$plots$aggregated$aggregated, "Group.1", "%")
colnames(atacseqkd_quantify_list$chipseeker$plots$aggregated$aggregated) <- c("Frequency", "sample","Feature")

ggplot(atacseqkd_quantify_list$chipseeker$plots$aggregated$aggregated, aes(x="Feature", y=Frequency, fill=Feature))+ facet_wrap( ~ sample, ncol=2, nrow=3) +
  geom_bar(width = 1, stat = "identity", color="white") + coord_polar("y", start=0) + theme_void() + scale_fill_manual(values=c("#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251", "#B565A7", "#009B77", "#DD4124", "#D65076")) 

atacseqkd_quantify_list[["chipseeker"]][["plots"]][["assigned"]][["shCRAMP1_shControl"]] <- as.data.frame(atacseqkd_quantify_list$chipseeker$plots$gr_anno$atacseqkd_diffbind_d_shCRAMP1_shControl_de_re.bed@anno)
atacseqkd_quantify_list[["chipseeker"]][["plots"]][["assigned"]][["shSUZ12_shControl"]] <- as.data.frame(atacseqkd_quantify_list$chipseeker$plots$gr_anno$atacseqkd_diffbind_d_shSUZ12_shControl_de_re.bed@anno)

write.table(atacseqkd_quantify_list$chipseeker$plots$assigned$shCRAMP1_shControl, "atacseqkd_shCRAMP1_shControl_chipseeker_all_features.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_quantify_list$chipseeker$plots$assigned$shSUZ12_shControl, "atacseqkd_shSUZ12_shControl_chipseeker_all_features.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)

write.table(atacseqkd_quantify_list$chipseeker$plots$assigned$shCRAMP1_shControl[which(atacseqkd_quantify_list$chipseeker$plots$assigned$shCRAMP1_shControl$annotation == "Distal Intergenic"),], "atacseqkd_shCRAMP1_shControl_chipseeker_intergenic.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_quantify_list$chipseeker$plots$assigned$shSUZ12_shControl[which(atacseqkd_quantify_list$chipseeker$plots$assigned$shSUZ12_shControl$annotation == "Distal Intergenic"),], "atacseqkd_shSUZ12_shControl_chipseeker_intergenic.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)

write.table(atacseqkd_quantify_list$chipseeker$plots$assigned$shCRAMP1_shControl[which(atacseqkd_quantify_list$chipseeker$plots$assigned$shCRAMP1_shControl$annotation != "Distal Intergenic"),], "atacseqkd_shCRAMP1_shControl_chipseeker_Nonintergenic.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_quantify_list$chipseeker$plots$assigned$shSUZ12_shControl[which(atacseqkd_quantify_list$chipseeker$plots$assigned$shSUZ12_shControl$annotation != "Distal Intergenic"),], "atacseqkd_shSUZ12_shControl_chipseeker_Nonintergenic.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)

################################################
###         CUT&RUN data: Knockout SE        ###
################################################

# deseq2 with consensus peaks (direct from nextflow)
cutnrunko_deseq2_path <- "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/"

cutnrunko_deseq2_list <- list()

cutnrunko_deseq2_list[["dds"]] <- load_rdata(cutnrunko_deseq2_path, "outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/deseq2/H3K27me3.consensus_peaks.dds.RData")
cutnrunko_deseq2_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_upstream_2kb.sorted.txt")
cutnrunko_deseq2_list[["counts"]] <- dds_counts(cutnrunko_deseq2_list$dds)
cutnrunko_deseq2_list[["samples_selected"]] <- sample_filter(cutnrunko_deseq2_list$counts, samplefilter = FALSE)
cutnrunko_deseq2_list[["processed_counts"]] <- cutnrunko_deseq2_list[["samples_selected"]][,c(1:4)]
cutnrunko_deseq2_list[["coldata"]] <- data.frame(colData(cutnrunko_deseq2_list$dds))
cutnrunko_deseq2_list$coldata <- data.frame(sample=cutnrunko_deseq2_list$coldata$sample, condition=factor(gsub('[0-9]+', '', unlist(lapply(strsplit(cutnrunko_deseq2_list$coldata$sample, "_"), function(x) x[2])))), replicate=factor(c(3,1,1,2)), sizeFactors=cutnrunko_deseq2_list$coldata$sizeFactor)
rownames(cutnrunko_deseq2_list$coldata) <- cutnrunko_deseq2_list$coldata$sample
cutnrunko_deseq2_list[["consensus"]] <- data.frame(fread("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/H3K27me3.consensus_peaks.bed"))
# check
all(rownames(cutnrunko_deseq2_list$coldata) == colnames(cutnrunko_deseq2_list$processed_counts)) #should print TRUE

# de analysis
cutnrunko_deseq2_list[["dar_analysis"]] <- deseq2_de(cutnrunko_deseq2_list$processed_counts, cutnrunko_deseq2_list$coldata, 0.05, 2)

# get positions of dars
cutnrunko_deseq2_list[["annotation"]] <- consensus_gene_anno("cutnrunko", "deseq2", cutnrunko_deseq2_list$dar_analysis$de_analysis, cutnrunko_deseq2_list$consensus, cutnrunko_deseq2_list$gene, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/", "txt", 0.05, 2)

plotPCA(cutnrunko_deseq2_list$dar_analysis$vstO, intgroup="sample", returnData=FALSE)

cutnrunko_deseq2_list[["aggregate_pergene"]][["KO_WT"]][["nearest"]][["de"]] <- aggregate_feature(cutnrunko_deseq2_list$annotation$KO_WT$nearest$de, c(3), "ens_gene", mean)
cutnrunko_deseq2_list[["aggregate_pergene"]][["KO_WT"]][["nearest"]][["up"]] <- aggregate_feature(cutnrunko_deseq2_list$annotation$KO_WT$nearest$up, c(3), "ens_gene", mean)
cutnrunko_deseq2_list[["aggregate_pergene"]][["KO_WT"]][["nearest"]][["down"]] <- aggregate_feature(cutnrunko_deseq2_list$annotation$KO_WT$nearest$down, c(3), "ens_gene", mean)

colnames(cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest$up$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest$down$aggregated) <- c("ensID_Gene", "mean_log2FC")

# quantify over various features
# now since I have information for groups to take
selected_rkd_gene_groups <- list.files(integration_path, pattern="cluster*", full.names = T)
integration_features_list[["selected_rkd_gene_groups"]] <- rearrange_groups(selected_rkd_gene_groups)
write.table(integration_features_list$selected_rkd_gene_groups[,c(2:4,8,9,5)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/","selected_rkd_gene_groups.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

cutnrunko_quantify_list <- list()
cutnrunko_quantify_list[["featurecounts"]][["featureslist"]]=list(gene_ext=c(".chr_1.5kb.txt", "gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), bins5kb=c("hg38_5kb.txt", "hg38_5kb.txt"), bins10kb=c("hg38_10kb.txt", "hg38_10kb.txt"), promoter=c("_gencodev41_promoter.txt", "gene_gencodev41_promoter.txt"))
cutnrunko_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_groups"]] <-  c("_rkd_gene_groups.txt", "selected_rkd_gene_groups.txt")
cutnrunko_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_uniq_diff"]] <-  c("_shC1_shS_up_pos.txt", "set_rnaseqde_merged_shC1_shS_up_pos.txt")
cutnrunko_quantify_list[["featurecounts"]][["pairs"]] <- list(H3K27me3_KO1_IgG=c("H3K27me3_KO1.mLb.clN.sorted.bam", "IgG_KO1.mLb.clN.sorted.bam"), H3K27me3_KO2_IgG=c("H3K27me3_KO2.mLb.clN.sorted.bam", "IgG_KO2.mLb.clN.sorted.bam"), H3K27me3_KO3_IgG=c("H3K27me3_KO3.mLb.clN.sorted.bam", "IgG_KO3.mLb.clN.sorted.bam"), H3K27me3_WT_IgG=c("H3K27me3_WT.mLb.clN.sorted.bam", "IgG_WT.mLb.clN.sorted.bam"))
cutnrunko_quantify_list[["featurecounts"]][["featuresmatrix"]] <- quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "cutnrunko", ".mLb.clN.sorted.bam", cutnrunko_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=FALSE, diff="multi", cutnrunko_quantify_list$featurecounts$pairs)
cutnrunko_quantify_list[["featurecounts"]][["featuresmatrix_plots"]] <- meta_plot(cutnrunko_quantify_list$featurecounts$featuresmatrix, c(9:12), rep=1, splitside=3)

# quantify overs atacseq dars
cutnrunko_quantify_list[["featurecounts"]][["atacseqdars"]] <- list(shCRAMP1_shControl=c("_diffbind_d_shCRAMP1_shControl_de_re.bed", "atacseqkd_diffbind_d_shCRAMP1_shControl_de_re.bed"), shSUZ12_shControl=c("_diffbind_d_shSUZ12_shControl_de_re.bed", "atacseqkd_diffbind_d_shSUZ12_shControl_de_re.bed"))
cutnrunko_quantify_list[["featurecounts"]][["pairs"]] <- list(H3K27me3_KO1_IgG=c("H3K27me3_KO1.mLb.clN.sorted.bam", "IgG_KO1.mLb.clN.sorted.bam"), H3K27me3_KO2_IgG=c("H3K27me3_KO2.mLb.clN.sorted.bam", "IgG_KO2.mLb.clN.sorted.bam"), H3K27me3_KO3_IgG=c("H3K27me3_KO3.mLb.clN.sorted.bam", "IgG_KO3.mLb.clN.sorted.bam"), H3K27me3_WT_IgG=c("H3K27me3_WT.mLb.clN.sorted.bam", "IgG_WT.mLb.clN.sorted.bam"))
cutnrunko_quantify_list[["featurecounts"]][["darsmatrix"]] <- quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "cutnrunko", ".mLb.clN.sorted.bam", cutnrunko_quantify_list$featurecounts$atacseqdars, c(5,4,1:3,6), "hg38", "IgG",pe=FALSE, diff="multi", cutnrunko_quantify_list$featurecounts$pairs)

ggplot(cutnrunko_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")
ggplot(cutnrunko_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, aes(x=value, color=col, linetype=id)) + geom_density() + scale_linetype_manual(values = c(rep("solid",1),rep("longdash",1),rep("dashed",1),rep("dotted",1))) + geom_vline(aes(xintercept=0), color="grey", linetype="dashed", size=0.2) + scale_color_manual(values = c("blue", "blue", "blue", "orange"))+theme_classic()
# 
# ggplot(cutnrunko_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, aes(x=col, y=value, color=id))  + geom_violin(trim=FALSE, position = position_dodge(0.4)) + theme_classic() + geom_hline(yintercept = 0, linetype="dotted")

ggplot(cutnrunko_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, aes(x=col, y=value, col =id)) + 
  geom_violin(aes(color = id), trim = FALSE, position = position_dodge(0.9), size=0.25, width = 1.5) +
  geom_boxplot(aes(color = id), width = 0.05, position = position_dodge(0.9))  + theme_classic() + geom_hline(yintercept = 0, linetype="dotted")

# remove none_up category
cutnrunko_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub <- dplyr::filter(cutnrunko_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, !grepl('none_up', row))

cutnrunko_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub$col <- factor(cutnrunko_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub$col, levels = c("H3K27me3_WT_IgG","H3K27me3_KO1_IgG","H3K27me3_KO2_IgG","H3K27me3_KO3_IgG"))

ggplot(cutnrunko_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")

ggplot(cutnrunko_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub, aes(x=col, y=value, col =id)) + 
  geom_violin(aes(fill = id), color="black", trim = FALSE, position = position_dodge(0.9), size=0.25, width = 1) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9))  + theme_classic() + 
  geom_hline(yintercept = 0, linetype="dotted")+ labs(title="CUT&RUN KO", x ="Samples", y = "H3K27me3 level")

# intersect with DE genes from RNA-Seq categories
overlap_cutnrunko_rnaseqkd_de_genes_mat <- matrix(0,4,4)
overlap_cutnrunko_rnaseqkd_de_genes_mat_exp <- matrix(0,4,4)
total_genes_expressed_in_rnaseqkd <- rnaseqkd_list$deseq2$de_analysis$c1509_shC1$all$ensID_Gene # take any categories after filtering the genes will be same, i.e expressed genes in c1509_shC1 before any filtering will same as c1509_shS
for (i in names(list_set_rnaseqde_merged_shC1_shS_up)){
  # print(list_set_rnaseqde_merged_shC1_shS_up[[i]])
  index_i <- match(i, names(list_set_rnaseqde_merged_shC1_shS_up))
  # print(index_i)
  # print(head(get(i)),2)
  # print(dim(get(i))[1])
  overlap_cutnrunko_rnaseqkd_de_genes_mat[index_i,1] <- length(list_set_rnaseqde_merged_shC1_shS_up[[i]])
  total_genes_expressed_in_rnaseqkd_categorised <- intersect(list_set_rnaseqde_merged_shC1_shS_up[[i]], total_genes_expressed_in_rnaseqkd)
  overlap_cutnrunko_rnaseqkd_de_genes_mat_exp[index_i,1] <- length(total_genes_expressed_in_rnaseqkd_categorised)
  rnaseqkd_de_genes_info <- list_set_rnaseqde_merged_shC1_shS_up[[i]]
  for (j in names(cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest)){
    index_j <- match(j, names(cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest))
    # print(index_j)
    print(paste0(index_i," vs ", index_j))
    print(paste0(i," vs ", j))
    print("performing for  all genes")
    cutnrunko_dbr_nearest_genes_info <- cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest[[j]][["aggregated"]]
    overlap_cutnrunko_rnaseqkd_de_genes <- intersect(rnaseqkd_de_genes_info, cutnrunko_dbr_nearest_genes_info$ensID_Gene)
    print(length(overlap_cutnrunko_rnaseqkd_de_genes))
    overlap_cutnrunko_rnaseqkd_de_genes_mat[index_i,index_j+1] <- length(overlap_cutnrunko_rnaseqkd_de_genes)
    print("performing for only expressed genes")
    # take only 
    cutnrunko_dbr_nearest_exp_genes_info <- intersect(total_genes_expressed_in_rnaseqkd, cutnrunko_dbr_nearest_genes_info$ensID_Gene)
    overlap_cutnrunko_rnaseqkd_exp_de_genes <- intersect(rnaseqkd_de_genes_info, cutnrunko_dbr_nearest_exp_genes_info)
    print(length(overlap_cutnrunko_rnaseqkd_exp_de_genes))
    overlap_cutnrunko_rnaseqkd_de_genes_mat_exp[index_i,index_j+1] <- length(overlap_cutnrunko_rnaseqkd_exp_de_genes)
  }
}
rownames(overlap_cutnrunko_rnaseqkd_de_genes_mat) <- names(list_set_rnaseqde_merged_shC1_shS_up)
colnames(overlap_cutnrunko_rnaseqkd_de_genes_mat) <- c("all","de", "up", "down")

rownames(overlap_cutnrunko_rnaseqkd_de_genes_mat_exp) <- names(list_set_rnaseqde_merged_shC1_shS_up)
colnames(overlap_cutnrunko_rnaseqkd_de_genes_mat_exp) <- c("all","de", "up", "down")

# de 
#shC1_up
1-phyper(13,277,19393-277,1404)
# shS_up
1-phyper(99,1054,19393-1054,1404)
# shC1_shS_up
1-phyper(79,443,19393-443,1404)

# up 
#shC1_up
1-phyper(3,277,19393-277,311)
# shS_up
1-phyper(13,1054,19393-1054,311)
# shC1_shS_up
1-phyper(7,443,19393-443,311)

# down 
#shC1_up
1-phyper(10,277,19393-277,1102)
# shS_up
1-phyper(88,1054,19393-1054,1102)
# shC1_shS_up
1-phyper(72,443,19393-443,1102)
# 

# for expresed genes
# de 
#shC1_up
1-phyper(13,277,19393-277,557)
# shS_up
1-phyper(99,1054,19393-1054,557)
# shC1_shS_up
1-phyper(79,443,19393-443,557)

# up 
#shC1_up
1-phyper(3,277,19393-277,75)
# shS_up
1-phyper(13,1054,19393-1054,75)
# shC1_shS_up
1-phyper(7,443,19393-443,75)

# down 
#shC1_up
1-phyper(10,277,19393-277,488)
# shS_up
1-phyper(88,1054,19393-1054,488)
# shC1_shS_up
1-phyper(72,443,19393-443,488)
#
# total =61852
# total white balls=277
# no. of white balls drawn = 13
# no. of black balls=61852-277 = 61575
# no. of balls drawn=1404
# 1-phyper(13,277,61575,1404)
# 1-phyper(79,443,61409,1404)
# 1-phyper(99,1054,60798,1404)
# 1-phyper(1213,60078,1774,1404)

################################################
###         CUT&RUN data: Knockdown SE       ###
################################################

# deseq2 with consensus peaks (direct from nextflow)
cutnrunkd_deseq2_path <- "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/"

cutnrunkd_deseq2_list <- list()

cutnrunkd_deseq2_list[["dds"]] <- load_rdata(cutnrunkd_deseq2_path, "outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/deseq2/H3K27me3.consensus_peaks.dds.RData")
cutnrunkd_deseq2_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_upstream_2kb.sorted.txt")
cutnrunkd_deseq2_list[["counts"]] <- dds_counts(cutnrunkd_deseq2_list$dds)
cutnrunkd_deseq2_list[["samples_selected"]] <- sample_filter(cutnrunkd_deseq2_list$counts, samplefilter = FALSE)
cutnrunkd_deseq2_list[["processed_counts"]] <- cutnrunkd_deseq2_list[["samples_selected"]][,c(1:3)]
cutnrunkd_deseq2_list[["consensus"]] <- data.frame(fread("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/H3K27me3.consensus_peaks.bed"))
cutnrunkd_deseq2_list[["normalized_counts"]] <- data.frame(edgeR::cpm(cutnrunkd_deseq2_list$processed_counts))
cutnrunkd_deseq2_list$normalized_counts["ensID_Gene"] <- rownames(cutnrunkd_deseq2_list$normalized_counts)

# get positions of peaks
cutnrunkd_deseq2_list[["annotation"]] <- consensus_gene_count_anno("cutnrunko", "deseq2", cutnrunkd_deseq2_list$normalized_counts, cutnrunkd_deseq2_list$consensus, cutnrunkd_deseq2_list$gene, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/", "txt")


# pca 
cutnrunkd_deseq2_list[["vst"]] <- DESeq2::vst(cutnrunkd_deseq2_list$dds)
plotPCA(cutnrunkd_deseq2_list$vst, intgroup="sample", returnData=FALSE)

cutnrunkd_deseq2_list[["aggregate_pergene"]][["nearest"]] <- aggregate_feature(cutnrunkd_deseq2_list$annotation$nearest, c(2:4), "ens_gene", sum)

rownames(cutnrunkd_deseq2_list$aggregate_pergene$nearest$aggregated) <- cutnrunkd_deseq2_list$aggregate_pergene$nearest$aggregated$Group.1
cutnrunkd_deseq2_list$aggregate_pergene$nearest$aggregated <- cutnrunkd_deseq2_list$aggregate_pergene$nearest$aggregated[,-1]
cutnrunkd_deseq2_list[["aggregate_pergene"]][["nearest"]][["log_sum_counts"]] <- log(data.frame(cutnrunkd_deseq2_list$aggregate_pergene$nearest$aggregated) + 1,2)


# subtract sample to shControl
for (i in colnames(cutnrunkd_deseq2_list$aggregate_pergene$nearest$log_sum_counts)){
  if (i %notlike% "shControl"){
    print(i)
    cutnrunkd_deseq2_list$aggregate_pergene$nearest$log_sum_counts[paste0(gsub("_","",gsub("H3K27me3","",i)),"_","shControl")] <- cutnrunkd_deseq2_list$aggregate_pergene$nearest$log_sum_counts[i] - cutnrunkd_deseq2_list$aggregate_pergene$nearest$log_sum_counts["shControl_H3K27me3"]
  }
}
cutnrunkd_deseq2_list$aggregate_pergene$nearest$log_sum_counts["ensID_Gene"] <- rownames(cutnrunkd_deseq2_list$aggregate_pergene$nearest$log_sum_counts)

cutnrunkd_quantify_list <- list()
cutnrunkd_quantify_list[["featurecounts"]][["featureslist"]]=list(gene_ext=c(".chr_1.5kb.txt", "gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), bins5kb=c("hg38_5kb.txt", "hg38_5kb.txt"), bins10kb=c("hg38_10kb.txt", "hg38_10kb.txt"), promoter=c("_gencodev41_promoter.txt", "gene_gencodev41_promoter.txt"))
cutnrunkd_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_groups"]] <-  c("_rkd_gene_groups.txt", "selected_rkd_gene_groups.txt")
cutnrunkd_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_uniq_diff"]] <-  c("_shC1_shS_up_pos.txt", "set_rnaseqde_merged_shC1_shS_up_pos.txt")
cutnrunkd_quantify_list[["featurecounts"]][["pairs"]] <- list(
  shControl_IgG = c("shControl_H3K27me3.mLb.clN.sorted.bam", "shControl_IgG.mLb.clN.sorted.bam"),
  shCRAMP1_IgG = c("shCRAMP1_H3K27me3.mLb.clN.sorted.bam","shCRAMP1_IgG.mLb.clN.sorted.bam"),
  shSUZ12_IgG = c("shSUZ12_H3K27me3.mLb.clN.sorted.bam", "shSUZ12_IgG.mLb.clN.sorted.bam"))
cutnrunkd_quantify_list[["featurecounts"]][["featuresmatrix"]] <- quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "cutnrunkd", ".mLb.clN.sorted.bam", cutnrunkd_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "shControl",pe=FALSE, diff="multi", cutnrunkd_quantify_list$featurecounts$pairs)
cutnrunkd_quantify_list[["featurecounts"]][["featuresmatrix_plots"]] <- meta_plot(cutnrunkd_quantify_list$featurecounts$featuresmatrix, c(7:9), rep=1, splitside=3)

ggplot(cutnrunkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")

# remove none_up category
cutnrunkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub <- dplyr::filter(cutnrunkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, !grepl('none_up', row))

ggplot(cutnrunkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")

ggplot(cutnrunkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub, aes(x=col, y=value, col =id)) + 
  geom_violin(aes(fill = id), color="black", trim = FALSE, position = position_dodge(0.9), size=0.25, width = 1) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9))  + theme_classic() + 
  geom_hline(yintercept = 0, linetype="dotted")+ labs(title="CUT&RUN KD", x ="Samples", y = "H3K27me3 level") +facet_wrap(~id)

cutnrunkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff[["comparisons"]] <- list( c("shControl_IgG", "shCRAMP1_IgG"), c("shControl_IgG", "shSUZ12_IgG"), c("shSUZ12_IgG", "shCRAMP1_IgG") )

# ggviolin(cutnrunkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st_sub, x = "col",
#          y = "value",
#          combine = TRUE, 
#          color = "id", palette = "jco",
#          add = "boxplot", width = 0.6) + stat_compare_means(comparisons = cutnrunkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$comparisons)

# Add global the p-value
################################################
###         CUT&TAG data: Wildtype  SE       ###
################################################
cutntagwt_quantify_list <- list()
cutntagwt_quantify_list[["featurecounts"]] <- quantify_featurecounts("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary", "cutntagwt","\\.bam$", "\\_peaks_id.bed$","histone_marks", c(12,4,1:3,7), pairedend=FALSE, "hg38", "_",merge_sites_files=FALSE)

cutntagwt_quantify_list[["bams"]] <- quantify_bams("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary", "cutntagwt","\\.bam$")
cutntagwt_labelpath = "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/histone_marks_label.txt"
cutntagwt_quantify_list[["summed"]] <- summedreads_per_feature(cutntagwt_quantify_list$featurecounts$histone_marks_countmatrix, cutntagwt_quantify_list$bams$total_reads, cutntagwt_labelpath)

cutntagwt_quantify_list$summed$all_features_ratio$sample <- gsub(".bam","",gsub(".mLb.clN.sorted.bam", "", cutntagwt_quantify_list$summed$all_features_ratio$sample))
cutntagwt_quantify_list$summed$all_features_ratio$sample <- factor(cutntagwt_quantify_list$summed$all_features_ratio$sample, levels = c("IgG","ENCFF508LLH","H1_4","H3K27me3","panH1"))

for (i in as.character(unique(cutntagwt_quantify_list$summed$all_features_ratio$sample))){
  if (i != "IgG"){
    print(i)
    cutntagwt_mysubset <- cutntagwt_quantify_list$summed$all_features_ratio[cutntagwt_quantify_list$summed$all_features_ratio$sample %in% c(i, "IgG"), ]
    cutntagwt_mysubset_plot <- ggplot(cutntagwt_mysubset, aes(x = label, y = Ratio)) +
      geom_col(aes(color = sample, fill = sample), position = position_dodge(0.7), width = 0.6) + scale_color_manual(values = c("black","black"))+
      scale_fill_manual(values = c("#d3d3d3", "#fde396"))+theme_classic()+
      theme(axis.text.x = element_text(angle = 90), text = element_text(size=25), plot.title = element_text(hjust = 0.5))+
      labs(title = i, y = "ratio = Reads_Inside/Reads_outside ", x = "Histone Marks")
    cutntagwt_quantify_list[["summed"]][["plot"]][[i]] <- cutntagwt_mysubset_plot
  }
}

# Combine the plots using patchwork
gridExtra::grid.arrange(cutntagwt_quantify_list$summed$plot$ENCFF508LLH, cutntagwt_quantify_list$summed$plot$H3K27me3, cutntagwt_quantify_list$summed$plot$panH1, cutntagwt_quantify_list$summed$plot$H1_4, ncol =1)

#save as  cutntagwt_combined_plot.pdf at pdf 40X23 portrait

# assign genes
cutntagwt_quantify_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_upstream_2kb.sorted.txt")

# plot heatmap 
cutntagwt_quantify_list$summed[["susbet_all_features_ratio"]] <- cutntagwt_quantify_list$summed$all_features_ratio[,c(2,3,7)]
cutntagwt_quantify_list$summed[["table_afr"]] <- data.frame(tidyr::pivot_wider(cutntagwt_quantify_list$summed$susbet_all_features_ratio, names_from = sample, values_from = Ratio))
rownames(cutntagwt_quantify_list$summed$table_afr) <- cutntagwt_quantify_list$summed$table_afr$label
cutntagwt_quantify_list$summed$table_afr <- cutntagwt_quantify_list$summed$table_afr[,-1]
cutntagwt_quantify_list$summed$table_afr <- cutntagwt_quantify_list$summed$table_afr - cutntagwt_quantify_list$summed$table_afr$IgG
cutntagwt_quantify_list$summed$table_afr <- t(cutntagwt_quantify_list$summed$table_afr)
breaksList <- seq(-0.4, 0.4, by = 0.01)
pheatmap::pheatmap(as.matrix(cutntagwt_quantify_list$summed$table_afr), na_col = "black",breaks = breaksList,
                   color = colorRampPalette(c("blue","white", "orange"))(length(breaksList)),
                   clustering_distance_cols = "euclidean", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D")

write.table(cutntagwt_quantify_list$summed$all_features, "cutntagwt_histonemarks_all_features.txt", sep="\t", quote = F, append = F, row.names = F, col.names = T)

# make a notlike function
`%notlike%` <- Negate('%like%')

# quantify over various features
cutntagwt_quantify_list[["featurecounts"]][["featureslist"]]=list(gene_ext=c(".chr_1.5kb.txt", "gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), bins5kb=c("hg38_5kb.txt", "hg38_5kb.txt"), bins10kb=c("hg38_10kb.txt", "hg38_10kb.txt"), promoter=c("_gencodev41_promoter.txt", "gene_gencodev41_promoter.txt"))
cutntagwt_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_groups"]] <-  c("_rkd_gene_groups.txt", "selected_rkd_gene_groups.txt")
cutntagwt_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_uniq_diff"]] <-  c("_shC1_shS_up_pos.txt", "set_rnaseqde_merged_shC1_shS_up_pos.txt")
cutntagwt_quantify_list[["featurecounts"]][["featuresmatrix"]] <- quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "cutntagwt", ".mLb.clN.sorted.bam", cutntagwt_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=FALSE, diff="single")
cutntagwt_quantify_list[["featurecounts"]][["featuresmatrix_plots"]] <- meta_plot(cutntagwt_quantify_list$featurecounts$featuresmatrix, c(5:7), rep=1, splitside=3)


# for check if two df are exactly same and to check if quantify_featurecounts_multifeature is created properly and doing the same thing as individual steps look at Check Notes

# quantify overs atacseq dars
cutntagwt_quantify_list[["featurecounts"]][["atacseqdars"]] <- list(shCRAMP1_shControl=c("_diffbind_d_shCRAMP1_shControl_de_re.bed", "atacseqkd_diffbind_d_shCRAMP1_shControl_de_re.bed"), shSUZ12_shControl=c("_diffbind_d_shSUZ12_shControl_de_re.bed", "atacseqkd_diffbind_d_shSUZ12_shControl_de_re.bed"))
cutntagwt_quantify_list[["featurecounts"]][["darsmatrix"]] <- quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "cutntagwt", ".mLb.clN.sorted.bam", cutntagwt_quantify_list$featurecounts$atacseqdars, c(5,4,1:3,6), "hg38", "IgG", pe=FALSE, diff="single")


ggplot(cutntagwt_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")

ggplot(cutntagwt_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, aes(x=col, y=value, col =id)) + 
  geom_violin(aes(color = id), trim = FALSE, position = position_dodge(0.9), size=0.25, width = 10) +
  geom_boxplot(aes(color = id), width = 0.05, position = position_dodge(0.9))  + theme_classic() + geom_hline(yintercept = 0, linetype="dotted")

##############################################################
###         CUT&TAG data: WT and KD and FLAG (PE)        ###
##############################################################
# ls *.fastq.gz -1 |  grep R1 | sort -k1,1 | awk -F'_' '{print $2"\t"$1"_"$2"_"$3"_"$4"_""R1""_"$6"\t"$1"_"$2"_"$3"_"$4"_""R2""_"$6}' > fastq_info.txt
# add header vim fastq_info.txt
# samplesheet prep for nextflow
# Rachael told me control for the FLAG labelled samples is "FLAG_H14_old".
# The one labelled as FLAG_old is actually the FLAG_H1.4_old sample. So I changed it, But we are not sure still
cutntagwkd_quantify_list <- list()
cutntagwkd_quantify_list[["nextflow"]][["info"]] <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_dec23/info.txt", header = T))
cutntagwkd_quantify_list[["nextflow"]][["sheet"]] <- data.frame(ID=rep(c(cutntagwkd_quantify_list$nextflow$info$ID),2),
                                                                replicate=rep(1,length(cutntagwkd_quantify_list$nextflow$info$ID) * 2),
                                                                group=gsub("-","_",rep(c(cutntagwkd_quantify_list$nextflow$info$Sample),2)),
                                                                control=rep(c("shControl_IgG","shControl_IgG","","shCRAMP1_IgG","shCRAMP1_IgG","","shSUZ12_IgG","shSUZ12_IgG","","WT_IgG","WT_IgG","WT_IgG","WT_IgG","WT_IgG","","FLAG_H14_old","FLAG_H14_old","FLAG_H14_old","FLAG_H14_old","FLAG_H14_old","",""),2)
)

cutntagwkd_quantify_list[["nextflow"]][["fastqs"]] <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_dec23/fastq_info.txt", header = T))
cutntagwkd_quantify_list[["nextflow"]][["samplesheet"]] <- merge(cutntagwkd_quantify_list$nextflow$sheet, cutntagwkd_quantify_list$nextflow$fastqs, by="ID")
cutntagwkd_quantify_list$nextflow$samplesheet <- cutntagwkd_quantify_list$nextflow$samplesheet %>% dplyr::distinct()
cutntagwkd_quantify_list$nextflow$samplesheet <- cutntagwkd_quantify_list$nextflow$samplesheet[order(cutntagwkd_quantify_list$nextflow$samplesheet$group),][,c(3,2,5,6,4)]
# write.table(cutntagwkd_quantify_list$nextflow$samplesheet, "/mnt/home3/reid/av638/cutntag/iva_lab_dec23/samplesheet.csv", sep=",", col.names = T, row.names = F, quote = F, append = F)

# histone marks
cutntagwkd_quantify_list[["featurecounts"]] <- quantify_featurecounts("/mnt/home3/reid/av638/cutntag/iva_lab_dec23/outfolder/02_alignment/bowtie2/target/markdup", "cutntagwkd","\\.bam$", "\\_peaks_id.bed$","histone_marks", c(12,4,1:3,7), pairedend=TRUE, "hg38", "_",merge_sites_files=FALSE)

cutntagwkd_quantify_list[["bams"]] <- quantify_bams("/mnt/home3/reid/av638/cutntag/iva_lab_dec23/outfolder/02_alignment/bowtie2/target/markdup", "cutntagwkd","\\.bam$")
cutntagwkd_labelpath = "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/histone_marks_label.txt"
cutntagwkd_quantify_list[["summed"]] <- summedreads_per_feature(cutntagwkd_quantify_list$featurecounts$histone_marks_countmatrix, cutntagwkd_quantify_list$bams$total_reads, cutntagwkd_labelpath)

# plot heatmap 
cutntagwkd_quantify_list$summed$all_features_ratio$sample <- gsub("_R1.target.markdup.sorted.bam","",cutntagwkd_quantify_list$summed$all_features_ratio$sample)
cutntagwkd_quantify_list$summed[["susbet_all_features_ratio"]] <- cutntagwkd_quantify_list$summed$all_features_ratio[,c(2,3,7)]
cutntagwkd_quantify_list$summed[["table_afr"]] <- data.frame(tidyr::pivot_wider(cutntagwkd_quantify_list$summed$susbet_all_features_ratio, names_from = sample, values_from = Ratio))
rownames(cutntagwkd_quantify_list$summed$table_afr) <- cutntagwkd_quantify_list$summed$table_afr$label
cutntagwkd_quantify_list$summed$table_afr <- cutntagwkd_quantify_list$summed$table_afr[,-1]
write.table(cutntagwkd_quantify_list$summed$all_features, "cutntagwkd_histonemarks_all_features.txt", sep="\t", quote = F, append = F, row.names = F, col.names = T)


cutntagwkd_quantify_list$summed[["diff_table_afr"]] <- NULL
for (i in c("FLAG", "WT", "shCRAMP1", "shSUZ12", "shControl")){
  if(i == "FLAG"){
    cutntagwkd_quantify_list$summed[["diff_table_afr"]][["FLAG"]] <- cutntagwkd_quantify_list$summed$table_afr[grepl(i, colnames(cutntagwkd_quantify_list$summed$table_afr))] - cutntagwkd_quantify_list$summed$table_afr[["FLAG_old"]]
  }else if(i == "WT"){
    cutntagwkd_quantify_list$summed[["diff_table_afr"]][["WT"]] <- cutntagwkd_quantify_list$summed$table_afr[grepl(i, colnames(cutntagwkd_quantify_list$summed$table_afr))] - cutntagwkd_quantify_list$summed$table_afr[["WT_IgG"]]
  }else if(i == "shCRAMP1"){
    cutntagwkd_quantify_list$summed[["diff_table_afr"]][["shCRAMP1"]] <- cutntagwkd_quantify_list$summed$table_afr[grepl(i, colnames(cutntagwkd_quantify_list$summed$table_afr))] - cutntagwkd_quantify_list$summed$table_afr[["shCRAMP1_IgG"]]
  }else if(i == "shSUZ12"){
    cutntagwkd_quantify_list$summed[["diff_table_afr"]][["shSUZ12"]] <- cutntagwkd_quantify_list$summed$table_afr[grepl(i, colnames(cutntagwkd_quantify_list$summed$table_afr))] - cutntagwkd_quantify_list$summed$table_afr[["shSUZ12_IgG"]]
  }else if(i == "shControl"){
    cutntagwkd_quantify_list$summed[["diff_table_afr"]][["shControl"]] <- cutntagwkd_quantify_list$summed$table_afr[grepl(i, colnames(cutntagwkd_quantify_list$summed$table_afr))] - cutntagwkd_quantify_list$summed$table_afr[["shControl_IgG"]]
  }
}
cutntagwkd_quantify_list$summed[["merged_diff_table_afr"]]<- do.call(cbind.data.frame, cutntagwkd_quantify_list$summed$diff_table_afr)
cutntagwkd_quantify_list$summed$merged_diff_table_afr <- t(cutntagwkd_quantify_list$summed$merged_diff_table_afr)
rownames(cutntagwkd_quantify_list$summed$merged_diff_table_afr) <- sapply(strsplit(rownames(cutntagwkd_quantify_list$summed$merged_diff_table_afr), "\\."), function(x) x[2])
cutntagwkd_quantify_list$summed$merged_diff_table_afr <- cutntagwkd_quantify_list$summed$merged_diff_table_afr[order(rownames(cutntagwkd_quantify_list$summed$merged_diff_table_afr)),]
breaksListp <- seq(0, 0.2, by = 0.01)
pheatmap::pheatmap(as.matrix(cutntagwkd_quantify_list$summed$merged_diff_table_afr), na_col = "grey",breaks = breaksListp,
                   color = colorRampPalette(c("white", "orange"))(length(breaksListp)),
                   clustering_distance_cols = "euclidean", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D")

pheatmap::pheatmap(as.matrix(cutntagwkd_quantify_list$summed$merged_diff_table_afr[c(6,7:17,19:22),]), na_col = "grey",breaks = breaksListp,
                   color = colorRampPalette(c("white", "orange"))(length(breaksListp)),
                   clustering_distance_cols = "euclidean", cluster_rows = T, cluster_cols = T, clustering_method = "ward.D")

corrplot::corrplot(cor(cutntagwkd_quantify_list$summed$merged_diff_table_afr))
cutntagwkd_quantify_list$summed[["pca"]] <- prcomp(cutntagwkd_quantify_list$summed$merged_diff_table_afr, center = TRUE, scale. = TRUE) 
biplot(cutntagwkd_quantify_list$summed$pca)
factoextra::fviz_pca_ind(cutntagwkd_quantify_list$summed$pca, geom = c("point", "text"))

# quantify over various features
cutntagwkd_quantify_list[["featurecounts"]][["featureslist"]]=list(gene_ext=c(".chr_1.5kb.txt", "gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), bins5kb=c("hg38_5kb.txt", "hg38_5kb.txt"), bins10kb=c("hg38_10kb.txt", "hg38_10kb.txt"), promoter=c("_gencodev41_promoter.txt", "gene_gencodev41_promoter.txt"))
cutntagwkd_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_groups"]] <-  c("_rkd_gene_groups.txt", "selected_rkd_gene_groups.txt")
cutntagwkd_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_uniq_diff"]] <-  c("_shC1_shS_up_pos.txt", "set_rnaseqde_merged_shC1_shS_up_pos.txt")
cutntagwkd_quantify_list[["featurecounts"]][["pairsinfo"]] <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_dec23/informationpairs.txt"))
cutntagwkd_quantify_list$featurecounts$pairsinfo <- cutntagwkd_quantify_list$featurecounts$pairsinfo %>% distinct()
cutntagwkd_quantify_list$featurecounts$pairsinfo["id"] <- paste0(cutntagwkd_quantify_list$featurecounts$pairsinfo$group, "_",cutntagwkd_quantify_list$featurecounts$pairsinfo$control)
cutntagwkd_quantify_list$featurecounts$pairsinfo$group <- paste0(cutntagwkd_quantify_list$featurecounts$pairsinfo$group, "_R1.target.markdup.sorted.bam")
cutntagwkd_quantify_list$featurecounts$pairsinfo$control <- paste0(cutntagwkd_quantify_list$featurecounts$pairsinfo$control, "_R1.target.markdup.sorted.bam")
cutntagwkd_quantify_list[["featurecounts"]][["pairs"]] <- lapply(split(cutntagwkd_quantify_list$featurecounts$pairsinfo[, c("group", "control")], cutntagwkd_quantify_list$featurecounts$pairsinfo$id), function(x) as.vector(unlist(x)))
cutntagwkd_quantify_list[["featurecounts"]][["featuresmatrix"]] <- 
  quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutntag/iva_lab_dec23/outfolder/02_alignment/bowtie2/target/markdup/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "cutntagwkd", "_R1.target.markdup.sorted.bam", 
                                         cutntagwkd_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=TRUE, diff="multi", 
                                         cutntagwkd_quantify_list$featurecounts$pairs)

cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt[["cpm"]] <- edgeR::cpm(cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$counts)
cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb[["pca_counts"]] <- make_pca(cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$cpm, "_R1.target.markdup.sorted.bam", "pca")

cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb[["pca_diff"]] <- make_pca(cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts[,c(23:38)], "_R1.target.markdup.sorted.bam", "pca")

cutntagwkd_quantify_list[["featurecounts"]][["featuresmatrix_plots"]] <- meta_plot(cutntagwkd_quantify_list$featurecounts$featuresmatrix, c(23:38), rep=1, splitside=3)
ggplot(cutntagwkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")+ coord_flip()
# save as cutntagwkd_quantify_list_featuresmatrix_plots.pdf

cutntagwkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_sub_st <- dplyr::filter(cutntagwkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, !grepl('FLAG|WT', col))
ggplot(cutntagwkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_sub_st, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")+ coord_flip()

cutntagwkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_subH1_4_st <- dplyr::filter(cutntagwkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, !grepl('FLAG|WT|H1_', col))
ggplot(cutntagwkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_subH1_4_st, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")+ coord_flip()

cutntagwkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_subH1_st <- dplyr::filter(cutntagwkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, !grepl('FLAG|WT|H14', col))
ggplot(cutntagwkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_subH1_st, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")+ coord_flip()

##############################################################
###          CUT&TAG data: WT RM mapped (PE)               ###
##############################################################
cutntagrmwt_quantify_list <- list()

cutntagrmwt_path <- "/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams/"

# histone marks
cutntagrmwt_quantify_list[["featurecounts"]] <- quantify_featurecounts("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams", "cutntagrmwt","\\.bam$", "\\_peaks_id.bed$","histone_marks", c(12,4,1:3,7), pairedend=TRUE, "hg38", "_",merge_sites_files=FALSE)

cutntagrmwt_quantify_list[["bams"]] <- quantify_bams("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams", "cutntagrmwt","\\.bam$")
cutntagrmwt_labelpath = "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/histone_marks_label.txt"
cutntagrmwt_quantify_list[["summed"]] <- summedreads_per_feature(cutntagrmwt_quantify_list$featurecounts$histone_marks_countmatrix, cutntagrmwt_quantify_list$bams$total_reads, cutntagrmwt_labelpath)

# plot heatmap 
cutntagrmwt_quantify_list$summed$all_features_ratio$sample <- gsub("ITit257-RM-|_sorted.bam","",cutntagrmwt_quantify_list$summed$all_features_ratio$sample)
cutntagrmwt_quantify_list$summed[["susbet_all_features_ratio"]] <- cutntagrmwt_quantify_list$summed$all_features_ratio[,c(2,3,7)]
cutntagrmwt_quantify_list$summed[["table_afr"]] <- data.frame(tidyr::pivot_wider(cutntagrmwt_quantify_list$summed$susbet_all_features_ratio, names_from = sample, values_from = Ratio))
rownames(cutntagrmwt_quantify_list$summed$table_afr) <- cutntagrmwt_quantify_list$summed$table_afr$label
cutntagrmwt_quantify_list$summed$table_afr <- cutntagrmwt_quantify_list$summed$table_afr[,-1]
write.table(cutntagrmwt_quantify_list$summed$all_features, "cutntagrmwt_histonemarks_all_features.txt", sep="\t", quote = F, append = F, row.names = F, col.names = T)

# subtraction of respective control
cutntagrmwt_quantify_list$summed[["diff_table_afr"]] <- NULL
for (i in c("F", "V5", "F|V5")){
  if(i == "F"){
    cutntagrmwt_quantify_list$summed[["diff_table_afr"]][["FLAG"]] <- cutntagrmwt_quantify_list$summed$table_afr[grepl(i, colnames(cutntagrmwt_quantify_list$summed$table_afr))] - cutntagrmwt_quantify_list$summed$table_afr[["F_S6"]]
  }else if(i == "V5"){
    cutntagrmwt_quantify_list$summed[["diff_table_afr"]][["V5"]] <- cutntagrmwt_quantify_list$summed$table_afr[grepl(i, colnames(cutntagrmwt_quantify_list$summed$table_afr))] - cutntagrmwt_quantify_list$summed$table_afr[["V5_S13"]]
  }else if(i == "F|V5"){
    #taking not equal strategy (!grepl for either F or V5)
    cutntagrmwt_quantify_list$summed[["diff_table_afr"]][["WT"]] <- cutntagrmwt_quantify_list$summed$table_afr[!grepl(i, colnames(cutntagrmwt_quantify_list$summed$table_afr))] - cutntagrmwt_quantify_list$summed$table_afr[["IgG_S22"]]
  }
}
cutntagrmwt_quantify_list$summed[["merged_diff_table_afr"]]<- do.call(cbind.data.frame, cutntagrmwt_quantify_list$summed$diff_table_afr)
cutntagrmwt_quantify_list$summed$merged_diff_table_afr <- t(cutntagrmwt_quantify_list$summed$merged_diff_table_afr)
rownames(cutntagrmwt_quantify_list$summed$merged_diff_table_afr) <- sapply(strsplit(rownames(cutntagrmwt_quantify_list$summed$merged_diff_table_afr), "_"), function(x) x[1])
cutntagrmwt_quantify_list$summed$merged_diff_table_afr <- cutntagrmwt_quantify_list$summed$merged_diff_table_afr[order(rownames(cutntagrmwt_quantify_list$summed$merged_diff_table_afr)),]
breaksListp <- seq(0, 0.2, by = 0.01)
pheatmap::pheatmap(as.matrix(cutntagrmwt_quantify_list$summed$merged_diff_table_afr[c(1,7,19,2:6,8:18,20:22),]), na_col = "grey",breaks = breaksListp,
                   color = colorRampPalette(c("white", "orange"))(length(breaksListp)),
                   clustering_distance_cols = "euclidean", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D")

pheatmap::pheatmap(as.matrix(cutntagrmwt_quantify_list$summed$merged_diff_table_afr), na_col = "grey",breaks = breaksListp,
                   color = colorRampPalette(c("white", "orange"))(length(breaksListp)),
                   clustering_distance_cols = "euclidean", cluster_rows = T, cluster_cols = T, clustering_method = "ward.D")

# mean the histone marks
cutntagrmwt_quantify_list$summed[["merged_diff_table_afr_st"]] <- data.frame(stack(as.matrix(do.call(cbind.data.frame, cutntagrmwt_quantify_list$summed$diff_table_afr))))
cutntagrmwt_quantify_list$summed$merged_diff_table_afr_st[["id"]] <- paste0(sapply(strsplit(as.character(cutntagrmwt_quantify_list$summed$merged_diff_table_afr_st$row), "_"),function(x) x[1]), "%", as.character(cutntagrmwt_quantify_list$summed$merged_diff_table_afr_st$col))
cutntagrmwt_quantify_list$summed[["aggregated"]] <- aggregate_feature(cutntagrmwt_quantify_list$summed$merged_diff_table_afr_st, c(3), "id", mean)
cutntagrmwt_quantify_list$summed$aggregated$aggregated[["histonemarks"]] <- sapply(strsplit(as.character(cutntagrmwt_quantify_list$summed$aggregated$aggregated$Group.1), "%"),function(x) x[1])
cutntagrmwt_quantify_list$summed$aggregated$aggregated[["sample"]] <- sapply(strsplit(as.character(cutntagrmwt_quantify_list$summed$aggregated$aggregated$Group.1), "%"),function(x) x[2])
cutntagrmwt_quantify_list$summed$aggregated$aggregated <- cutntagrmwt_quantify_list$summed$aggregated$aggregated[,-1]
cutntagrmwt_quantify_list$summed[["histone_mark_df"]] <- data.frame(tidyr::pivot_wider(cutntagrmwt_quantify_list$summed$aggregated$aggregated, names_from = sample, values_from = x))
rownames(cutntagrmwt_quantify_list$summed$histone_mark_df) <- cutntagrmwt_quantify_list$summed$histone_mark_df$histonemarks
cutntagrmwt_quantify_list$summed$histone_mark_df <- cutntagrmwt_quantify_list$summed$histone_mark_df[,-1]
# remove H1.5new (RM suggested me)
pheatmap::pheatmap(as.matrix(t(cutntagrmwt_quantify_list$summed$histone_mark_df)[c(1,7,19,2:6,8:15,17:18,20:22),c(1,3,10:11,12,2,7:9,4,5,6)]), na_col = "grey",breaks = breaksListp,
                   color = colorRampPalette(c("white", "darkred"))(length(breaksListp)),
                   clustering_distance_cols = "euclidean", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D")
# remove H3K27me3 (RM/Iva suggested as it is too dark) and H1Xs H1.5new
breaksListp1 <- seq(0, 0.05, by = 0.01)
pheatmap::pheatmap(as.matrix(t(cutntagrmwt_quantify_list$summed$histone_mark_df)[c(1,7,19,2:5,8:11,13:15,17),c(1,3,10:11,12,2,7:9,4,5,6)]), na_col = "grey",breaks = breaksListp1,border_color = "grey",
                   color = colorRampPalette(c("white", "darkred"))(length(breaksListp1)),
                   clustering_distance_cols = "euclidean", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D")

# for correlation remove controls as it interfere with stdev calculation
corrplot::corrplot(cor(cutntagrmwt_quantify_list$summed$merged_diff_table_afr[-c(1,7,19),]))
corrplot::corrplot(cor(t(cutntagrmwt_quantify_list$summed$merged_diff_table_afr[-c(1,7,19),])), method = 'square', addCoef.col = 'black', col= colorRampPalette(c("darkred", "yellow", "navy"))(20))
cutntagrmwt_quantify_list$summed[["pca"]] <- prcomp(cutntagrmwt_quantify_list$summed$merged_diff_table_afr, center = TRUE, scale. = TRUE) 
biplot(cutntagrmwt_quantify_list$summed$pca)
factoextra::fviz_pca_ind(cutntagrmwt_quantify_list$summed$pca, geom = c("point", "text"))

# genome wide and other features assessment
# quantify over various features
cutntagrmwt_quantify_list[["featurecounts"]][["featureslist"]]=list(gene_ext=c(".chr_1.5kb.txt", "gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), bins5kb=c("hg38_5kb.txt", "hg38_5kb.txt"), bins10kb=c("hg38_10kb.txt", "hg38_10kb.txt"), bins80kb=c("hg38_80kb.txt", "hg38_80kb.txt"), promoter=c("_gencodev41_promoter.txt", "gene_gencodev41_promoter.txt"))
cutntagrmwt_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_groups"]] <-  c("_rkd_gene_groups.txt", "selected_rkd_gene_groups.txt")
cutntagrmwt_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_uniq_diff"]] <-  c("_shC1_shS_up_pos.txt", "set_rnaseqde_merged_shC1_shS_up_pos.txt")
cutntagrmwt_quantify_list[["featurecounts"]][["pairsinfo"]] <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams/informationpairs.txt"))
cutntagrmwt_quantify_list$featurecounts$pairsinfo <- cutntagrmwt_quantify_list$featurecounts$pairsinfo %>% distinct()
cutntagrmwt_quantify_list$featurecounts$pairsinfo["id"] <- paste0(cutntagrmwt_quantify_list$featurecounts$pairsinfo$group, "_",cutntagrmwt_quantify_list$featurecounts$pairsinfo$control)

# replacing all - by . is required
cutntagrmwt_quantify_list$featurecounts$pairsinfo <- mutate_all(cutntagrmwt_quantify_list$featurecounts$pairsinfo, ~gsub("-", ".", .))

cutntagrmwt_quantify_list[["featurecounts"]][["pairs"]] <- lapply(split(cutntagrmwt_quantify_list$featurecounts$pairsinfo[, c("groupbam", "controlbam")], cutntagrmwt_quantify_list$featurecounts$pairsinfo$id), function(x) as.vector(unlist(x)))
cutntagrmwt_quantify_list[["featurecounts"]][["featuresmatrix"]] <- 
  quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "cutntagrmwt", "_sorted.bam", 
                                         cutntagrmwt_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=TRUE, diff="multi", 
                                         cutntagrmwt_quantify_list$featurecounts$pairs)

# plot and add the regression line
ggplot(cutntagrmwt_quantify_list$featurecounts$featuresmatrix$bins80kb$log_normalized_counts, aes(x=H14_S17_IgG_S22, y=FH14_S3_F_S6)) + 
  geom_point(shape=18, color="#2c7bb0", size=1)+
  geom_smooth(method=lm, color="darkred") + theme_classic()

cutntagrmwt_quantify_list[["featurecounts"]][["featuresmatrixdiv"]] <- 
  quantify_featurecounts_multifeature_div("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "cutntagrmwt", "_sorted.bam", 
                                         cutntagrmwt_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=TRUE, diff="multi", 
                                         cutntagrmwt_quantify_list$featurecounts$pairs)


# plot and add the regression line
ggplot(cutntagrmwt_quantify_list$featurecounts$featuresmatrixdiv$bins80kb$log_normalized_counts, aes(x=H14_S17_IgG_S22, y=FH14_S3_F_S6)) + 
  geom_point(shape=18, color="#2c7bb0", size=1)+
  geom_smooth(method=lm, color="darkred") + theme_classic()+
  stat_cor(method = "pearson", label.x = -3, label.y = 8, r.digits = 3, aes(label = ..r.label..))

cor(cutntagrmwt_quantify_list$featurecounts$featuresmatrixdiv$bins80kb$log_normalized_counts$H14_S17_IgG_S22, cutntagrmwt_quantify_list$featurecounts$featuresmatrixdiv$bins80kb$log_normalized_counts$FH14_S3_F_S6, method = 'pearson')

# deeptools was for used correlation analysis
# import output of plotCorrealtion
# Heatmap/Correlation plot
# ----------- divlog2 method  -----------------#
# All samples
cutntagrmwt_correlation_heatmap_divlog2 <- fread(paste0(cutntagrmwt_path, "cutntagrmwt_correlation_heatmap_divlog2.txt"), header = TRUE)
colnames(cutntagrmwt_correlation_heatmap_divlog2) <- gsub("'","",colnames(cutntagrmwt_correlation_heatmap_divlog2))
cutntagrmwt_correlation_heatmap_divlog2$V1 <- gsub("'","",cutntagrmwt_correlation_heatmap_divlog2$V1)
cutntagrmwt_correlation_heatmap_divlog2 <- data.frame(cutntagrmwt_correlation_heatmap_divlog2)
rownames(cutntagrmwt_correlation_heatmap_divlog2) <- cutntagrmwt_correlation_heatmap_divlog2$V1
cutntagrmwt_correlation_heatmap_divlog2 <- cutntagrmwt_correlation_heatmap_divlog2[,-1]
breaksList_cor = seq(-1, 1, by = 0.01)
pheatmap(as.matrix(cutntagrmwt_correlation_heatmap_divlog2), display_numbers = TRUE,  number_color = "black", breaks = breaksList_cor, legend_breaks = c(-1, 0, 1), fontsize_number = 9, border_color = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList_cor)), clustering_distance_cols = "euclidean", clustering_method = "ward.D")

# All ENCODE and GEO mix
all_ENC_cutntagrmwt_correlation_heatmap_divlog2 <- fread(paste0(cutntagrmwt_path, "all_ENC_cutntagrmwt_correlation_heatmap_divlog2.txt"), header = TRUE)
colnames(all_ENC_cutntagrmwt_correlation_heatmap_divlog2) <- gsub("'","",colnames(all_ENC_cutntagrmwt_correlation_heatmap_divlog2))
all_ENC_cutntagrmwt_correlation_heatmap_divlog2$V1 <- gsub("'","",all_ENC_cutntagrmwt_correlation_heatmap_divlog2$V1)
all_ENC_cutntagrmwt_correlation_heatmap_divlog2 <- data.frame(all_ENC_cutntagrmwt_correlation_heatmap_divlog2)
rownames(all_ENC_cutntagrmwt_correlation_heatmap_divlog2) <- all_ENC_cutntagrmwt_correlation_heatmap_divlog2$V1
all_ENC_cutntagrmwt_correlation_heatmap_divlog2 <- all_ENC_cutntagrmwt_correlation_heatmap_divlog2[,-1]
breaksList_cor = seq(-1, 1, by = 0.01)
pheatmap(as.matrix(all_ENC_cutntagrmwt_correlation_heatmap_divlog2), display_numbers = TRUE,  number_color = "black", breaks = breaksList_cor, legend_breaks = c(-1, 0, 1), fontsize_number = 6, border_color = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList_cor)), clustering_distance_cols = "euclidean", clustering_method = "ward.D")

# H1s endo
H1endo_cutntagrmwt_correlation_heatmap_divlog2 <- fread(paste0(cutntagrmwt_path, "H1endo_cutntagrmwt_correlation_heatmap_divlog2.txt"), header = TRUE)
colnames(H1endo_cutntagrmwt_correlation_heatmap_divlog2) <- gsub("'","",colnames(H1endo_cutntagrmwt_correlation_heatmap_divlog2))
H1endo_cutntagrmwt_correlation_heatmap_divlog2$V1 <- gsub("'","",H1endo_cutntagrmwt_correlation_heatmap_divlog2$V1)
H1endo_cutntagrmwt_correlation_heatmap_divlog2 <- data.frame(H1endo_cutntagrmwt_correlation_heatmap_divlog2)
rownames(H1endo_cutntagrmwt_correlation_heatmap_divlog2) <- H1endo_cutntagrmwt_correlation_heatmap_divlog2$V1
H1endo_cutntagrmwt_correlation_heatmap_divlog2 <- H1endo_cutntagrmwt_correlation_heatmap_divlog2[,-1]
breaksList_cor = seq(-1, 1, by = 0.01)
pheatmap(as.matrix(H1endo_cutntagrmwt_correlation_heatmap_divlog2), display_numbers = TRUE,  number_color = "black", breaks = breaksList_cor, legend_breaks = c(-1, 0, 1), fontsize_number = 9, border_color = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList_cor)), clustering_distance_cols = "euclidean", clustering_method = "ward.D")

# H1s endo ENCODE and GEO mix
H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2 <- fread(paste0(cutntagrmwt_path, "H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2.txt"), header = TRUE)
colnames(H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2) <- gsub("'","",colnames(H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2))
H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2$V1 <- gsub("'","",H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2$V1)
H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2 <- data.frame(H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2)
rownames(H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2) <- H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2$V1
H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2 <- H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2[,-1]
breaksList_cor = seq(-1, 1, by = 0.01)
pheatmap(as.matrix(H1endo_ENC_cutntagrmwt_correlation_heatmap_divlog2), display_numbers = TRUE,  number_color = "black", breaks = breaksList_cor, legend_breaks = c(-1, 0, 1), fontsize_number = 9, border_color = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList_cor)), clustering_distance_cols = "euclidean", clustering_method = "ward.D")

# H1s endo spearman
H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp <- fread(paste0(cutntagrmwt_path, "H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp.txt"), header = TRUE)
colnames(H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp) <- gsub("'","",colnames(H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp))
H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp$V1 <- gsub("'","",H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp$V1)
H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp <- data.frame(H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp)
rownames(H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp) <- H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp$V1
H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp <- H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp[,-1]
breaksList_sp = seq(0, 1, by = 0.01)
pheatmap(as.matrix(H1endo_cutntagrmwt_correlation_heatmap_divlog2_sp), display_numbers = TRUE,  number_color = "black", breaks = breaksList_sp, legend_breaks = c(0, 1), fontsize_number = 9, border_color = NA, color = colorRampPalette((brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList_sp)), clustering_distance_cols = "euclidean", clustering_method = "ward.D")

# ----------- subtraction method -----------------#
# All samples
cutntagrmwt_correlation_heatmap_subtract <- fread(paste0(cutntagrmwt_path, "cutntagrmwt_correlation_heatmap_sub.txt"), header = TRUE)
colnames(cutntagrmwt_correlation_heatmap_subtract) <- gsub("'","",colnames(cutntagrmwt_correlation_heatmap_subtract))
cutntagrmwt_correlation_heatmap_subtract$V1 <- gsub("'","",cutntagrmwt_correlation_heatmap_subtract$V1)
cutntagrmwt_correlation_heatmap_subtract <- data.frame(cutntagrmwt_correlation_heatmap_subtract)
rownames(cutntagrmwt_correlation_heatmap_subtract) <- cutntagrmwt_correlation_heatmap_subtract$V1
cutntagrmwt_correlation_heatmap_subtract <- cutntagrmwt_correlation_heatmap_subtract[,-1]
breaksList_cor = seq(-1, 1, by = 0.01)
pheatmap(as.matrix(cutntagrmwt_correlation_heatmap_subtract), display_numbers = TRUE,  number_color = "black", breaks = breaksList_cor, legend_breaks = c(-1, 0, 1), fontsize_number = 9, border_color = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList_cor)), clustering_distance_cols = "euclidean", clustering_method = "ward.D")

# H1s samples
H1endo_cutntagrmwt_correlation_heatmap_subtract <- fread(paste0(cutntagrmwt_path, "H1endo_cutntagrmwt_correlation_heatmap_sub.txt"), header = TRUE)
colnames(H1endo_cutntagrmwt_correlation_heatmap_subtract) <- gsub("'","",colnames(H1endo_cutntagrmwt_correlation_heatmap_subtract))
H1endo_cutntagrmwt_correlation_heatmap_subtract$V1 <- gsub("'","",H1endo_cutntagrmwt_correlation_heatmap_subtract$V1)
H1endo_cutntagrmwt_correlation_heatmap_subtract <- data.frame(H1endo_cutntagrmwt_correlation_heatmap_subtract)
rownames(H1endo_cutntagrmwt_correlation_heatmap_subtract) <- H1endo_cutntagrmwt_correlation_heatmap_subtract$V1
H1endo_cutntagrmwt_correlation_heatmap_subtract <- H1endo_cutntagrmwt_correlation_heatmap_subtract[,-1]
breaksList_cor = seq(-1, 1, by = 0.01)
pheatmap(as.matrix(H1endo_cutntagrmwt_correlation_heatmap_subtract), display_numbers = TRUE,  number_color = "black", breaks = breaksList_cor, legend_breaks = c(-1, 0, 1), fontsize_number = 9, border_color = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList_cor)), clustering_distance_cols = "euclidean", clustering_method = "ward.D")

# Scatter Plot
cutntagrmwt_correlation_scatter_divlog2 <- data.frame(fread(paste0(cutntagrmwt_path, "cutntagrmwt_correlation_scores_per_bin_divlog2.tab"), header = TRUE))
colnames(cutntagrmwt_correlation_scatter_divlog2) <- sapply(strsplit(colnames(cutntagrmwt_correlation_scatter_divlog2), "\\."), function(x) x[2])
rownames(cutntagrmwt_correlation_scatter_divlog2) <- do.call(paste, c(cutntagrmwt_correlation_scatter_divlog2[,c(1:3)], sep="%"))
cutntagrmwt_correlation_scatter_divlog2 <- cutntagrmwt_correlation_scatter_divlog2[,c(4:22)]
cutntagrmwt_quantify_list[["deeptools"]][["scatterplot"]] <- scatterplot(cutntagrmwt_correlation_scatter_divlog2)

# Annotation peaks
cutntagrmwt_quantify_list[["chipseeker"]][["plots"]] <- run_chipseeker("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams/endoH1peaks/", "cutntagrmwt", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "*.broadPeak")

cutntagrmwt_quantify_list[["chipseeker"]][["plots"]][["pie_annodf"]] <- replot_pie_chipseeker(cutntagrmwt_quantify_list$chipseeker$plots$gr_anno)
cutntagrmwt_quantify_list$chipseeker$plots$pie_annodf$df_st$sample <- gsub("endo|_peaks.broadPeak","",cutntagrmwt_quantify_list$chipseeker$plots$pie_annodf$df_st$sample)

ggplot(cutntagrmwt_quantify_list$chipseeker$plots$pie_annodf$df_st, aes(x="Feature", y=Frequency, fill=Feature))+ facet_wrap( ~ sample, ncol=2, nrow=3) +
  geom_bar(width = 1, stat = "identity", color="white") + coord_polar("y", start=0) + theme_void() + scale_fill_manual(values=c("#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251", "#B565A7", "#009B77", "#DD4124", "#D65076")) 

# RM told me to further simplify the pie chart to promoter, gene body and intergenic regions
category_chipseeker_df <- data.frame(id= c("Promoter" , "Gene_body", "Gene_body", "Gene_body", "Gene_body", "Gene_body", "Gene_body", "Intergenic", "Downstream"), Feature=c(unique(as.character(cutntagrmwt_quantify_list$chipseeker$plots$pie_annodf$df_st$Feature))))

cutntagrmwt_quantify_list[["chipseeker"]][["plots"]][["id_pie_annodf"]] <- merge(cutntagrmwt_quantify_list$chipseeker$plots$pie_annodf$df_st, category_chipseeker_df, by="Feature")

cutntagrmwt_quantify_list[["chipseeker"]][["plots"]][["id_pie_annodf"]]["id_to_agg"] <- do.call(paste, c(cutntagrmwt_quantify_list[["chipseeker"]][["plots"]][["id_pie_annodf"]][,c(3:4)], sep="%"))

cutntagrmwt_quantify_list[["chipseeker"]][["plots"]][["aggregated"]] <- aggregate_feature(cutntagrmwt_quantify_list$chipseeker$plots$id_pie_annodf, c(2), "id_to_agg", sum)

cutntagrmwt_quantify_list$chipseeker$plots$aggregated$aggregated <- cSplit(cutntagrmwt_quantify_list$chipseeker$plots$aggregated$aggregated, "Group.1", "%")
colnames(cutntagrmwt_quantify_list$chipseeker$plots$aggregated$aggregated) <- c("Frequency", "sample","Feature")

ggplot(cutntagrmwt_quantify_list$chipseeker$plots$aggregated$aggregated, aes(x="Feature", y=Frequency, fill=Feature))+ facet_wrap( ~ sample, ncol=2, nrow=3) +
  geom_bar(width = 1, stat = "identity", color="white") + coord_polar("y", start=0) + theme_void() + scale_fill_manual(values=c("#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251", "#B565A7", "#009B77", "#DD4124", "#D65076")) 

#################################
######## Public data ############
#################################
public_list <- list()
public_list[["peakfilelist"]] <- list.files("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "AHF_peaks_id.bed")
public_list[["gene"]] <- read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration_rm/", "gene_gencode_human_upstream_2kb.sorted.txt")
public_list[["peaks_anno"]] <- annotate_peaks_to_genes(public_list$peakfilelist, public_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration","chipseq", "encode", "txt")
write.table((public_list$peaks_anno$H3K27me3_ENCFF801AHF_peaks_id.bed$nearest$all[,c(13:15,19,18)] %>% dplyr::distinct()), "H3K27me3_ENCFF801AHF_peaks_id_neareast_anno.txt", sep="\t", quote = F, append = F, col.names = F, row.names = F)

# Integration of omics data
data_integration_list <- list()

# for gene
data_integration_list[["ensID_gene"]][["list"]] <- 
  list(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$de %>% select(r_shC1=2, ensID_Gene=7),
       rnaseqkd_list$deseq2$de_analysis$shS_c1509$de %>% select(r_shS=2, ensID_Gene=7),
       atacseqkd_deseq2_list$aggregate_pergene$shCRAMP1_shControl$nearest$de$aggregated %>% select(a_shCRAMP1=2, ensID_Gene=1),
       atacseqkd_deseq2_list$aggregate_pergene$shSUZ12_shControl$nearest$de$aggregated %>% select(a_shSUZ12=2, ensID_Gene=1),
       cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest$de$aggregated %>% select(cko_KO=2, ensID_Gene=1),
       cutnrunkd_deseq2_list$aggregate_pergene$nearest$log_sum_counts[,c(4:6)],
       cutntagwt_quantify_list$featurecounts$featuresmatrix$gene_ext$log_normalized_counts %>% select(5:7, ensID_Gene=8),
       cutntagwkd_quantify_list$featurecounts$featuresmatrix$gene_ext$log_normalized_counts %>% dplyr::select(23:38, ensID_Gene=39))

data_integration_list[["ensID_gene"]][["merged"]] <- integrate_featurecounts(data_integration_list$ensID_gene$list, "ensID_Gene")

# Some predictable
col_fun = circlize::colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
col_fun2 = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data_integration_list$ensID_gene$merged[["somepredictable"]] <- data_integration_list$ensID_gene$merged$df[complete.cases(data_integration_list$ensID_gene$merged$df[, c("r_shC1")]), ]
data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat1"]] <- ComplexHeatmap::Heatmap(as.matrix(data_integration_list$ensID_gene$merged$somepredictable[,c(1:5)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, col = col_fun)
data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat2"]] <- ComplexHeatmap::Heatmap(as.matrix(data_integration_list$ensID_gene$merged$somepredictable[,c(6:7)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, col = col_fun2)
data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat3"]] <- ComplexHeatmap::Heatmap(as.matrix(data_integration_list$ensID_gene$merged$somepredictable[,c(8:10)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, col = col_fun)
data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat4"]] <- ComplexHeatmap::Heatmap(as.matrix(data_integration_list$ensID_gene$merged$somepredictable[,c(11:26)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, col = col_fun)

data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat1"]]+ data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat2"]] + data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat3"]] + data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat4"]]

#complex_heatmap_data_integration_list_ensID_gene_merged_df_somepredictable_ch.pdf

data_integration_list$ensID_gene$merged[["pca"]][["somepredictable"]]  <- make_pca(data_integration_list$ensID_gene$merged$df, "none", "pca")

# for bins5kb from featurecounts
data_integration_list[["bins5kb"]][["list"]] <- 
  list(atacseqkd_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts %>% dplyr::select(7:10, bins5kb=11),
       cutnrunko_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts %>% dplyr::select(9:12, bins5kb=13),
       cutnrunkd_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts %>% dplyr::select(7:9, bins5kb=10),
       cutntagwt_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts %>% dplyr::select(5:7, bins5kb=8),
       cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts %>% dplyr::select(23:38, bins5kb=39))

data_integration_list[["bins5kb"]][["merged"]] <- integrate_featurecounts(data_integration_list$bins5kb$list, "bins5kb")

# PCA on bins of subset of samples using bins5kb diff counts
data_integration_list$bins5kb$merged[["bin_prcomp"]] <- prcomp(data_integration_list$bins5kb$merged$df[,c(1:4)], center = TRUE, scale. = TRUE)

data_integration_list$bins5kb$merged[["bin_pca"]] <- factoextra::fviz_pca_ind(data_integration_list$bins5kb$merged$bin_prcomp, geom = c("point"))


factoextra::fviz_pca_ind(data_integration_list$bins5kb$merged$bin_prcomp, 
             geom = c("point"), 
             col.ind = "contrib",  # Color by loadings
             gradient.cols = c("blue", "green"),  # Color gradient
             pointsize = 0.1,  # Reduce point size based on loadings
             title = "PCA Plot with Loadings Color and Reduced Point Size")

# for bins from bedtools coverage
intergration_omics_list <- list()
bin_cov_atacseq_path <- list.files(path="/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library", pattern = "*.txt_coverage.pe.bed", full.names = T)
bin_cov_cutnrun_KO_path <- list.files(path="/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary", pattern = "*.txt_coverage_se.bed", full.names = T)
bin_cov_cutnrun_KD_path <- list.files(path="/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary", pattern = "*.txt_coverage_se.bed", full.names = T)
bin_cov_cutntag_path <- list.files(path="/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary", pattern = "*.txt_coverage_se.bed", full.names = T)
bin_path_combined <- c(bin_cov_atacseq_path, bin_cov_cutnrun_KO_path, bin_cov_cutnrun_KD_path, bin_cov_cutntag_path)

intergration_omics_list[["bin5kb"]] <- integrate_bedtools(bin_path_combined, "5kb", c(1:3), seq(4,175,7))
intergration_omics_list[["bin10kb"]] <- integrate_bedtools(bin_path_combined, "10kb", c(1:3), seq(4,175,7))

#<--------------------------># 
# all cut&tag histone marks merge cutntagwt and cutntagwkd
intergration_omics_list[["histonemarks"]][["cutntag"]] <- rbind(cutntagwt_quantify_list$summed$table_afr, cutntagwkd_quantify_list$summed$merged_diff_table_afr)

pheatmap::pheatmap(as.matrix(intergration_omics_list$histonemarks$cutntag[c(1:5,11,12:22,24:27),]), na_col = "grey",breaks = breaksListp,
                   color = colorRampPalette(c("white", "orange"))(length(breaksListp)),
                   clustering_distance_cols = "euclidean", cluster_rows = T, cluster_cols = T, clustering_method = "ward.D")

intergration_omics_list[["featurecounts"]][["featuresmatrix_plots"]][["selected_rkd_uniq_diff"]][["cutntagwt_cutntagwkd"]] <- 
  rbind.data.frame(cutntagwt_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st, 
                   cutntagwkd_quantify_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$df_st)

ggplot(intergration_omics_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$cutntagwt_cutntagwkd, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")+ coord_flip()

intergration_omics_list[["featurecounts"]][["featuresmatrix_plots"]][["selected_rkd_uniq_diff"]][["cutntagwt_cutntagwkd_filt"]] <- dplyr::filter(intergration_omics_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$cutntagwt_cutntagwkd, !grepl('FLAG|WT', col))

ggplot(intergration_omics_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$cutntagwt_cutntagwkd_filt, aes(x=col, y=value, color=id))  + geom_boxplot(width=0.7) + theme_classic()+ geom_hline(yintercept = 0, linetype="dotted")+ coord_flip()
intergration_omics_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$cutntagwt_cutntagwkd_filt <- intergration_omics_list$featurecounts$featuresmatrix_plots$selected_rkd_uniq_diff$cutntagwt_cutntagwkd_filt[,c(3,2)]

# save as cutntagwkd_quantify_list_featuresmatrix_plots.pdf

# integrate older and new cutntag data

# combine
# cutntagrmwt_quantify_list$summed$merged_diff_table_afr
# cutntagwkd_quantify_list$summed$merged_diff_table_afr
intergration_omics_list$histonemarks[["cutntag_rmwt_wkd"]][["merged_diff_table_afr"]] <- rbind.data.frame(cutntagrmwt_quantify_list$summed$merged_diff_table_afr, cutntagwkd_quantify_list$summed$merged_diff_table_afr)

pheatmap::pheatmap(as.matrix(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$merged_diff_table_afr[c(1:19, 28:29, 43:44),]), na_col = "grey",breaks = breaksListp,
                   color = colorRampPalette(c("white", "darkred"))(length(breaksListp)),
                   clustering_distance_cols = "euclidean", cluster_rows = T, cluster_cols = T, clustering_method = "ward.D")

# mean the histone marks
intergration_omics_list$histonemarks$cutntag_rmwt_wkd[["merged_diff_table_afr_st"]] <- data.frame(stack(as.matrix(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$merged_diff_table_afr)))
intergration_omics_list$histonemarks$cutntag_rmwt_wkd$merged_diff_table_afr_st[["id"]] <- paste0(sapply(strsplit(as.character(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$merged_diff_table_afr_st$col), "_"),function(x) x[1]), "%", as.character(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$merged_diff_table_afr_st$row))
intergration_omics_list$histonemarks$cutntag_rmwt_wkd[["aggregated"]] <- aggregate_feature(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$merged_diff_table_afr_st, c(3), "id", mean)
intergration_omics_list$histonemarks$cutntag_rmwt_wkd$aggregated$aggregated[["histonemarks"]] <- sapply(strsplit(as.character(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$aggregated$aggregated$Group.1), "%"),function(x) x[1])
intergration_omics_list$histonemarks$cutntag_rmwt_wkd$aggregated$aggregated[["sample"]] <- sapply(strsplit(as.character(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$aggregated$aggregated$Group.1), "%"),function(x) x[2])
intergration_omics_list$histonemarks$cutntag_rmwt_wkd$aggregated$aggregated <- intergration_omics_list$histonemarks$cutntag_rmwt_wkd$aggregated$aggregated[,-1]
intergration_omics_list$histonemarks$cutntag_rmwt_wkd[["histone_mark_df"]] <- data.frame(tidyr::pivot_wider(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$aggregated$aggregated, names_from = sample, values_from = x))
rownames(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$histone_mark_df) <- intergration_omics_list$histonemarks$cutntag_rmwt_wkd$histone_mark_df$histonemarks
intergration_omics_list$histonemarks$cutntag_rmwt_wkd$histone_mark_df <- intergration_omics_list$histonemarks$cutntag_rmwt_wkd$histone_mark_df[,-1]

breaksListp1 <- seq(0, 0.1, by = 0.01)
pheatmap::pheatmap(as.matrix(t(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$histone_mark_df[,c(6:7, 8:13, 23:28, 33:41)])), na_col = "grey",breaks = breaksListp1,
                   color = colorRampPalette(c("white", "darkred"))(length(breaksListp1)),
                   clustering_distance_cols = "euclidean", cluster_rows = T, cluster_cols = T, clustering_method = "ward.D")



# H1s endo ENCODE and GEO mix cutntag new and old
H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2 <- fread(paste0(cutntagrmwt_path, "H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2.txt"), header = TRUE)
colnames(H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2) <- gsub("'","",colnames(H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2))
H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2$V1 <- gsub("'","",H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2$V1)
H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2 <- data.frame(H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2)
rownames(H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2) <- H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2$V1
H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2 <- H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2[,-1]
breaksList_cor = seq(-1, 1, by = 0.01)
pheatmap(as.matrix(H1endo_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2), display_numbers = TRUE,  number_color = "black", breaks = breaksList_cor, legend_breaks = c(-1, 0, 1), fontsize_number = 9, border_color = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList_cor)), clustering_distance_cols = "euclidean", clustering_method = "ward.D")

# All ENCODE and GEO mix cutntag new and old
all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2 <- fread(paste0(cutntagrmwt_path, "all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2.txt"), header = TRUE)
colnames(all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2) <- gsub("'","",colnames(all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2))
all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2$V1 <- gsub("'","",all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2$V1)
all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2 <- data.frame(all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2)
rownames(all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2) <- all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2$V1
all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2 <- all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2[,-1]
breaksList_cor = seq(-1, 1, by = 0.01)
pheatmap(as.matrix(all_ENC_cutntagrmwt_wkd_correlation_heatmap_divlog2), display_numbers = TRUE,  number_color = "black", breaks = breaksList_cor, legend_breaks = c(-1, 0, 1), fontsize_number = 7, border_color = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList_cor)), clustering_distance_cols = "euclidean", clustering_method = "ward.D")


# remove H1.5new (RM suggested me)
pheatmap::pheatmap(as.matrix(t(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$histone_mark_df)[c(1,7,19,2:6,8:15,17:18,20:22),c(1,3,10:11,12,2,7:9,4,5,6)]), na_col = "grey",breaks = breaksListp,
                   color = colorRampPalette(c("white", "darkred"))(length(breaksListp)),
                   clustering_distance_cols = "euclidean", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D")
# remove H3K27me3 (RM/Iva suggested as it is too dark) and H1Xs H1.5new
breaksListp1 <- seq(0, 0.05, by = 0.01)
pheatmap::pheatmap(as.matrix(t(intergration_omics_list$histonemarks$cutntag_rmwt_wkd$histone_mark_df)[c(1,7,19,2:5,8:11,13:15,17),c(1,3,10:11,12,2,7:9,4,5,6)]), na_col = "grey",breaks = breaksListp1,border_color = "grey",
                   color = colorRampPalette(c("white", "darkred"))(length(breaksListp1)),
                   clustering_distance_cols = "euclidean", cluster_rows = F, cluster_cols = F, clustering_method = "ward.D")


# integrate atacseq-dars and histone marks data 
histonemarks_list1 <- list.files("/mnt/beegfs6/home3/reid/av638/atacseq/iva_lab_gencode/integration/pairwise/", pattern = "_peaks_id.bed", full.names = TRUE)
atacseqkd_dars_list1 <- list.files("/mnt/beegfs6/home3/reid/av638/atacseq/iva_lab_gencode/integration/pairwise/", pattern = "atacseqkd*", full.names = TRUE)

enrichment_counter(histonemarks_list1, atacseqkd_dars_list1)

################################################ END OF ANALYSIS ###########################################################

#  Check Notes:
# > if (identical(cutntagwt_quantify_list$featurecounts$gene_ext$gene_ext_countmatrix$gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt$counts[,c(2:5)], cutntagwt_quantify_list$featurecounts$featuresmatrix$gene_ext$gene_ext_countmatrix$gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt$counts)) {
#   +     print("Data frames are exactly the same.")
#   + } else {
#     +     print("Data frames are different.")
#     + }
# [1] "Data frames are exactly the same."
# > 
#   > if (identical(cutntagwt_quantify_list$featurecounts$promoter$promoter_countmatrix$gene_gencodev41_promoter.txt$counts[,c(2:5)], cutntagwt_quantify_list$featurecounts$featuresmatrix$promoter$promoter_countmatrix$gene_gencodev41_promoter.txt$counts)) {
#     +     print("Data frames are exactly the same.")
#     + } else {
#       +     print("Data frames are different.")
#       + }
# [1] "Data frames are exactly the same."
# > 
#   > if (identical(cutntagwt_quantify_list$featurecounts$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$counts[,c(2:5)], cutntagwt_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$counts)) {
#     +     print("Data frames are exactly the same.")
#     + } else {
#       +     print("Data frames are different.")
#       + }
# [1] "Data frames are exactly the same."
# > 
#   > if (identical(cutntagwt_quantify_list$featurecounts$bins10kb$bins10kb_countmatrix$hg38_10kb.txt$counts[,c(2:5)], cutntagwt_quantify_list$featurecounts$featuresmatrix$bins10kb$bins10kb_countmatrix$hg38_10kb.txt$counts)) {
#     +     print("Data frames are exactly the same.")
#     + } else {
#       +     print("Data frames are different.")
#       + }
# [1] "Data frames are exactly the same."