### Load packages

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

setwd("/mnt/beegfs6/home3/reid/av638/atacseq/iva_lab_gencode/integration")

################################################
###          create functions.               ###
################################################

# load rdata
fn_load_rdata <- function(path, rdata) {
  load(paste0(path, rdata))
  return(dds)
}

# import genome file 
fn_read_genefile <- function(path, file){
  out <- read.table(paste0(path, file), header = F, stringsAsFactors = F)
  colnames(out) <- c("chr","start", "end", "strand", "type", "ensID", "gene")
  out["ens_gene"] <- paste0(out$ensID, "%", out$gene)
  return(out)
}

# deseq2 based analysis
fn_dds_counts <- function(dds){
  counts <- DESeq2::counts(dds, normalized=FALSE)
  return(counts)
}

fn_sample_filter <- function(countmatrix, samplefilter, columnsremove){
  if (samplefilter == TRUE){
    countmatrix <- countmatrix[,-columnsremove] # columns is vector eg, c(7:9)
    countmatrix <- data.frame(countmatrix)
    countmatrix["ensID"] <- rownames(countmatrix)
    return(countmatrix)
  }else {
    countmatrix <- data.frame(countmatrix)
    countmatrix["ensID"] <- rownames(countmatrix)
    return(countmatrix)
  }
}

fn_reassign_genenames <- function(countmatrix, genedf, columnskeep){
  genecounts <- merge(countmatrix, genedf ,by.x = "ensID", by.y ="ensID", all.y =T)
  rownames(genecounts) <- genecounts$ens_gene
  genecounts <- genecounts[,columnskeep]
  return(genecounts)
}

fn_make_coldata_dds <- function(dds, samplefilter,rowsremove){
  if(samplefilter == TRUE){
    coldata <- data.frame(colData(dds))
    colnames(coldata) <- c("sample", "condition", "replicate", "sizeFactor")
    coldata <- coldata[-rowsremove,] # columns is vector eg, c(7:9)
    coldata$condition <- factor(coldata$condition)
    coldata$replicate <- factor(coldata$replicate)
    return(coldata)
  }else{
    coldata <- data.frame(colData(dds))
    colnames(coldata) <- c("sample", "condition", "replicate", "sizeFactor")
    coldata$condition <- factor(coldata$condition)
    coldata$replicate <- factor(coldata$replicate)
    return(coldata)
  }
}

fn_deseq2 <- function(counts, coldata, fdr, fc){
  storelist <- list()
  ddsO <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = coldata, design = ~ condition)
  keep <- rowSums(counts(ddsO)) > 10
  ddsO <- ddsO[keep,]
  # de analysis
  storelist[["ddsO"]] <- DESeq2::DESeq(ddsO)
  # vst
  vstO <- DESeq2::vst(ddsO)
  storelist[["vstO"]] <- vstO
  # pca
  message("Plotting PCA...")
  storelist[["pca"]] <- DESeq2::plotPCA(storelist$vstO, intgroup="condition", returnData=TRUE)
  contrast_list <- list()
  ma_list <- list()
  message("Performing Deseq2 analysis...")
  for (cont1 in unique(coldata$condition)){
    for (cont2 in unique(coldata$condition)){
      if (cont1!=cont2){
        print(paste0(cont1,"_",cont2))
        assayseqde <- DESeq2::results(storelist$ddsO, contrast=c("condition", cont1, cont2))
        assayseqde = assayseqde[order(rownames(assayseqde)),]
        assayseqde["ensID_Gene"] <- rownames(assayseqde)
        assayseqdedf <- data.frame(assayseqde)
        assayseqdedf["ensID_Gene"] <- rownames(assayseqdedf)
        message("Applying fdr and logfc cutoff...")
        # padj filter
        assayseqde$threshold <- as.logical(assayseqde$padj < fdr)
        assayseqde0.05 <- data.frame(assayseqde[which(assayseqde$threshold == TRUE),])
        # logfc filter
        assayseqde0.05_de <- assayseqde0.05 %>% dplyr::filter((log2FoldChange > log(fc,2)) | (log2FoldChange < -log(fc,2)))
        assayseqde0.05_up <- assayseqde0.05 %>% dplyr::filter((log2FoldChange > log(fc,2)))
        assayseqde0.05_down <- assayseqde0.05 %>% dplyr::filter((log2FoldChange < -log(fc,2)))
        contrast_list[[paste0(cont1,"_",cont2)]] <- list(all=assayseqdedf, de=assayseqde0.05_de, up=assayseqde0.05_up, down=assayseqde0.05_down)
        message("Plotting MA plot...")
        ma_list[[paste0(cont1,"_",cont2)]] <- ggmaplot(assayseqdedf, fdr = fdr, fc = fc, size = 0.3, palette = c("#D22B2B", "#1465AC", "darkgray"),
                                                       genenames = unlist(lapply(strsplit(assayseqdedf$ensID_Gene, "%"), function(x) x[2])), legend = "top", top = 20, font.label = c("bold", 5),
                                                       font.legend = "bold", font.main = "bold", ggtheme = ggplot2::theme_minimal())
      }
    }
  }
  storelist[["de_analysis"]] <- contrast_list
  storelist[["ma_plots"]] <- ma_list
  return(storelist)
}


fn_plot_venn <- function(list, color){
  ggvenn_list_venn <- ggvenn(list, fill_color = color, stroke_size = 0.5, set_name_size = 4)
  return(ggvenn_list_venn)
}
# edger based analysis
fn_edger_create <- function(data, coldata){
  storelist <- list()
  message("creating group... ")
  group <- factor(coldata$Condition, levels = unique(coldata$Condition))
  
  message("creating dgelist object... ")
  storelist[["dgelist"]] <- edgeR::DGEList(data, group = group)
  keep <- filterByExpr(storelist$dgelist, storelist$design)
  storelist$dgelist <- storelist$dgelist[keep,,keep.lib.sizes=FALSE]
  
  message("performing normalisation... ")
  storelist$dgelist <- calcNormFactors(storelist$dgelist)
  
  message("creating design... ")
  storelist[["design"]] <- model.matrix(~group)
  colnames(storelist$design) <- levels(storelist$dgelist$samples$group)
  
  message("estimating dispersion... ")
  storelist$dgelist <- estimateGLMCommonDisp(storelist$dgelist, storelist$design)
  storelist$dgelist <- estimateGLMTrendedDisp(storelist$dgelist, storelist$design) 
  storelist$dgelist <- estimateGLMTagwiseDisp(storelist$dgelist, storelist$design)
  
  message("performing MDS analysis")
  storelist[["mds"]] <- limma::plotMDS(storelist$dgelist, col=c("red", "blue", "blue", "blue"), pch=6,cex = 3)
  storelist[["bcv"]] <- plotBCV(storelist$dgelist)
  
  message("performing model fitting using design... ")
  storelist[["lmfit"]] <- glmFit(storelist$dgelist, storelist$design)
  return(storelist)
}

fn_edger_de <- function(lmfit, design, fdr, fc){  
  storelist <- list()
  delist <- list()
  message("obtaining de results...")
  for (i in c(2,3)){
    lrt <- glmLRT(lmfit, coef=i)
    delist[["class"]] <- topTags(lrt, n = Inf, sort.by = "PValue", p.value = 1, adjust.method="BH")
    delist[["all"]] <- delist$class$table
    comparison <- unique(delist$class$comparison)
    delist$all["coordinates"] <- rownames(delist$all)
    delist$all <- data.frame(splitstackshape::cSplit(delist$all, "coordinates", "%"))
    delist$all <- delist$all[,c((dim(delist$all)[2] -2): (dim(delist$all)[2]), 1:(dim(delist$all)[2]-3))]
    delist[["de"]] <- delist$all %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2) | logFC < -log(fc,2)))
    delist[["up"]] <- delist$all %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2)))
    delist[["down"]] <- delist$all %>% dplyr::filter(FDR < fdr & (logFC < -log(fc,2)))
    storelist[["dars"]][[paste0("coef_", comparison)]] <- delist
  }
  return(storelist)
}

# diffbind based analysis
fn_diffbind_prep <- function(samplesheet){
  diffblist <- list()
  # Creating diffbind object 
  obj <- dba(sampleSheet=samplesheet,  scoreCol=5, minOverlap=1)
  diffblist[["dba_obj"]] <- obj
  # Plot of number of peaks that appear in atleast one, two , three so on ... until total samples together.
  ovp_rate <- dba.overlap(obj,mode=DBA_OLAP_RATE)
  diffblist[["ovp_rate"]] <- ovp_rate
  return(diffblist)
}

fn_diffbind_count <- function(obj, summit){
  diffblist <- list()
  # Calculate a binding matrix with scores based on read counts,  summits=100 for ATAC-Seq # https://www.biostars.org/p/9493721/
  diffblist[["dba_obj"]] <- dba.count(obj, minOverlap=1, summits=100, score = DBA_SCORE_READS)
  diffblist[["info"]]  <- dba.show(diffblist$dba_obj)
  return(diffblist)
}

fn_diffbind_de <- function(obj, fdr, fc){
  diffblist <- list()
  message("Performing normalisation...")
  diffblist[["norm"]] <- dba.normalize(obj, method=DBA_ALL_METHODS)
  message("Performing Diffbind analysis...")
  message("Performing Contrasts...")
  diffblist$norm <- dba.contrast(diffblist$norm, minMembers = 2,categories=DBA_CONDITION)
  contrasts <- dba.show(diffblist$norm, bContrasts=TRUE)
  diffblist$norm <- dba.analyze(diffblist$norm, method=DBA_ALL_METHODS)
  contrastlist <- list()
  for (i in as.numeric(rownames(contrasts))){
    print(paste0("performing analysis for contrast",i, " i.e. ", paste0(as.character(contrasts[i,][,c(4,2)]), collapse = "_")))
    storelist <- list()
    message("Fetching DE info...")
    # Extracting results, contrast 1 means , th means threshold for FDR, which if 1 give all sites 
    storelist[["all"]] <- data.frame(dba.report(diffblist$norm, method=DBA_EDGER, contrast = i, th=1))
    storelist[["de"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2) | Fold < -log(fc,2)))
    storelist[["up"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2)))
    storelist[["down"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold < -log(fc,2)))
    message("Performing PCA plot...")
    storelist[["pca"]] <- dba.plotPCA(diffblist$norm, contrast=i, method=DBA_EDGER, attributes=DBA_FACTOR, label=DBA_ID)
    message("Performing venn diagran...")
    storelist[["venn"]] <- dba.plotVenn(diffblist$norm, contrast=i, method=DBA_ALL_METHODS)
    message("Performing MA-plot...")
    storelist[["ma_plot"]] <- dba.plotMA(diffblist$norm, method=DBA_EDGER)
    contrastlist[[paste0("contrast_",paste0(as.character(contrasts[i,][,c(4,2)]), collapse = "_"))]] <- storelist
  }
  diffblist[["contrasts"]] <- contrastlist
  return(diffblist)
}


# limma based analysis
fn_limma_create <- function(data, coldata){
  storelist <- list()
  message(" creating design ... ")
  storelist[["design"]] <- model.matrix(~0 + coldata$Condition)
  colnames(storelist$design) <- levels(coldata$Condition)
  message(" creating dgelist object... ")
  storelist[["dgelist"]] <- edgeR::DGEList(data)
  keep <- filterByExpr(storelist$dgelist, storelist$design)
  storelist$dgelist <- storelist$dgelist[keep,,keep.lib.sizes=FALSE]
  storelist$limma <- calcNormFactors(storelist$dgelist)
  message(" performing voom transformation... ")
  storelist[["voom"]] <- limma::voom(storelist$limma, storelist$design)
  message(" performing model fitting using design... ")
  storelist[["fit"]] <- lmFit(storelist$voom, storelist$design)
  return(storelist)
}

fn_limma_de <- function(lmfit, contrasts, fdr, fc){
  message("Performing second fit...")
  storelist <- list()
  storelist[["fit2"]] <- contrasts.fit(lmfit, contrasts)
  storelist$fit2 <- eBayes(storelist$fit2)
  storelist[["decide"]] <- decideTests(storelist$fit2)
  storelist[["venn"]] <- vennDiagram(storelist$decide)
  
  message("obtaining results...")
  delist <- list()
  for (cont in colnames(contrasts)){
    delist[["all"]] <- topTable(storelist$fit2, number=Inf,coef=cont, adjust="BH")
    delist$all["coordinates"] <- rownames(delist$all)
    delist$all <- data.frame(splitstackshape::cSplit(delist$all, "coordinates", "%"))
    delist$all <- delist$all[,c((dim(delist$all)[2] -2): (dim(delist$all)[2]), 1:(dim(delist$all)[2]-3))]
    delist[["de"]] <- delist$all %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2) | logFC < -log(fc,2)))
    delist[["up"]] <- delist$all %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2)))
    delist[["down"]] <- delist$all %>% dplyr::filter(adj.P.Val < fdr & (logFC < -log(fc,2)))
    storelist[["dars"]][[paste0("coef_", gsub("-","_",cont))]] <- delist
  }
  return(storelist)
}


# annotation to genes and nearest distance filter for all packages
fn_nearest_gene_anno <- function(assaytype, software, contrasts, genefile, path, outformat, fdr, fc){
  contrastlist <- list()
  for (i in names(contrasts)){
    print(i)
    all_df <- data.frame(contrasts[[i]][["all"]])
    message("writing all data...")
    write.table(all_df, paste0(path,"/",assaytype, "_",software,"_",i,"_","all.",outformat), sep="\t", quote = F, append = F, row.names = F, col.names = F)
    write.table(genefile, paste0(path,"/","gene_gencode_v41_out.",outformat), sep="\t", quote = F, append = F, row.names = F, col.names = F)
    message("sorting...")
    system(paste0("sort -k1,1 -k2,2n ", paste0(path,"/",assaytype, "_",software,"_",i,"_","all.",outformat),  " | grep chr > ", paste0(path,"/",assaytype, "_",software,"_",i,"_","all_sorted.",outformat)))
    message("annotation to genes ...")
    all_anno <- data.frame(fread(cmd=paste0("bedtools closest -a ", paste0(path,"/",assaytype, "_",software,"_",i,"_","all_sorted.",outformat), " -b ", paste0(path,"/","gene_gencode_v41_out.",outformat), " -d")))
    colnames(all_anno) <- c(colnames(all_df), "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_type", "ensID", "gene", "ens_gene","Distance")
    storelist <- list()
    storelist[["all"]][["all"]] <- all_anno
    message("getting peaks nearest distance to genes...")
    nearest_distance = 5000
    all_anno_nearest <- all_anno %>% dplyr::filter(Distance < nearest_distance )
    storelist[["nearest"]][["all"]] <- all_anno_nearest
    if (software == "diffbind"){
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2) | Fold < -log(fc,2))) 
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (Fold < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2) | Fold < -log(fc,2))) 
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(FDR < fdr & (Fold < -log(fc,2)))
      contrastlist[[i]] <- storelist
    }else if (software == "limma") {
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2) | logFC < -log(fc,2))) 
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(adj.P.Val < fdr & (logFC < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2) | logFC < -log(fc,2))) 
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(adj.P.Val < fdr & (logFC < -log(fc,2)))
      contrastlist[[i]] <- storelist
    }else if (software == "edger") {
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2) | logFC < -log(fc,2))) 
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2) | logFC < -log(fc,2))) 
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC < -log(fc,2)))
      contrastlist[[i]] <- storelist
    }else if (software == "deseq2") {
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2) | log2FoldChange < -log(fc,2))) 
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(padj < fdr & (log2FoldChange < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2) | log2FoldChange < -log(fc,2))) 
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(padj < fdr & (log2FoldChange < -log(fc,2)))
      contrastlist[[i]] <- storelist
    }
  }
  return(contrastlist)
}

fn_consensus_gene_anno <- function(assaytype, software, contrasts, consensus_bed, genefile, path, outformat, fdr, fc){
  message("writing consensus peaks data...")
  write.table(consensus_bed, paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), sep="\t", quote = F, append = F, row.names = F, col.names = F)
  
  message("writing gene peaks data...")
  write.table(genefile, paste0(path,"/","gene_gencode_v41_out.",outformat), sep="\t", quote = F, append = F, row.names = F, col.names = F)
  
  message("sorting consensus and genefile...")
  system(paste0("sort -k1,1 -k2,2n ", paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), " | grep chr > ", paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), ".chr.bed"))
  system(paste0("sort -k1,1 -k2,2n ", paste0(path,"/","gene_gencode_v41_out.",outformat), " | grep chr > ", paste0(path,"/","gene_gencode_v41_out.",outformat), ".sorted.chr.txt"))
  
  message("annotating consensus peaks file to genes...")
  system(paste0("bedtools closest -a ", paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), ".chr.bed", " -b ", paste0(path,"/","gene_gencode_v41_out.",outformat), ".sorted.chr.txt", " -d > ", paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), ".chr.anno.bed"))

  consensus_to_gene <- data.frame(fread(paste0(paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), ".chr.anno.bed")))
  colnames(consensus_to_gene) <- c("peak_chr","peak_start", "peak_end","peak_interval","peak_score","peak_strand" ,"gene_chr", "gene_start", "gene_end", "gene_strand", "gene_type", "ensID", "gene", "ens_gene","Distance")
  
  contrastlist <- list()
  for (i in names(contrasts)){
    print(i)
    all_df <- data.frame(contrasts[[i]][["all"]])
    
    message("annotation to consensus gene ...")
    all_anno <- merge(all_df, consensus_to_gene, by.x="ensID_Gene", by.y="peak_interval")
    storelist <- list()
    storelist[["all"]][["all"]] <- all_anno
    
    message("getting peaks nearest distance to genes...")
    nearest_distance = 5000
    all_anno_nearest <- all_anno %>% dplyr::filter(Distance < nearest_distance )
    storelist[["nearest"]][["all"]] <- all_anno_nearest
    if (software == "deseq2"){
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2) | log2FoldChange < -log(fc,2))) 
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(padj < fdr & (log2FoldChange < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2) | log2FoldChange < -log(fc,2))) 
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(padj < fdr & (log2FoldChange < -log(fc,2)))
      contrastlist[[i]] <- storelist
    }
  }
  return(contrastlist)
}

fn_consensus_gene_count_anno <- function(assaytype, software, normcountmatrix, consensus_bed, genefile, path, outformat){
  message("writing consensus peaks data...")
  write.table(consensus_bed, paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), sep="\t", quote = F, append = F, row.names = F, col.names = F)
  
  message("writing gene peaks data...")
  write.table(genefile, paste0(path,"/","gene_gencode_v41_out.",outformat), sep="\t", quote = F, append = F, row.names = F, col.names = F)
  
  message("sorting consensus and genefile...")
  system(paste0("sort -k1,1 -k2,2n ", paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), " | grep chr > ", paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), ".chr.bed"))
  system(paste0("sort -k1,1 -k2,2n ", paste0(path,"/","gene_gencode_v41_out.",outformat), " | grep chr > ", paste0(path,"/","gene_gencode_v41_out.",outformat), ".sorted.chr.txt"))
  
  message("annotating consensus peaks file to genes...")
  system(paste0("bedtools closest -a ", paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), ".chr.bed", " -b ", paste0(path,"/","gene_gencode_v41_out.",outformat), ".sorted.chr.txt", " -d > ", paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), ".chr.anno.bed"))
  
  consensus_to_gene <- data.frame(fread(paste0(paste0(path,"/",assaytype, "_","consensus_peaks.",outformat), ".chr.anno.bed")))
  colnames(consensus_to_gene) <- c("peak_chr","peak_start", "peak_end","peak_interval","peak_score","peak_strand" ,"gene_chr", "gene_start", "gene_end", "gene_strand", "gene_type", "ensID", "gene", "ens_gene","Distance")
  
  storelist <- list()
  all_df <- data.frame(normcountmatrix)
  message("annotation to consensus gene ...")
  all_anno <- merge(all_df, consensus_to_gene, by.x="ensID_Gene", by.y="peak_interval")
  all_anno <- all_anno %>% dplyr::distinct()
  storelist[["all"]] <- all_anno
  message("getting peaks nearest distance to genes...")
  nearest_distance = 5000
  all_anno_nearest <- all_anno %>% dplyr::filter(Distance < nearest_distance )
  storelist[["nearest"]] <- all_anno_nearest
  return(storelist)
}

# export bed files for dars, dbrs
fn_export_bed <- function(path, annotation, type, de_pattern, outformat){
  for (names in names(annotation)){
    for (de_pattern in names(annotation[[names]][[type]])){
      write.table(annotation[[names]][[type]][[de_pattern]], paste0(path, "/",
                                                                    gsub("\\$", "_", deparse(substitute(annotation))),"_",names, "_",type, "_",de_pattern, ".",outformat), 
                  sep="\t", append = F, quote=F, col.names = T, row.names = F)
    }
  }
}

# run chipseeker
"peakformat=*.broadPeak"
"txdb=TxDb.Hsapiens.UCSC.hg38.knownGene"
"org_db=org.Hs.eg.db"
"assaytype=atacseqkd"
fn_run_chipseeker <- function(peaks_path, assaytype, txdb, org_db, peakformat){
  storelist <- list()
  peakfiles <- list.files(path=peaks_path,pattern = peakformat)
  for (peaksfile in peakfiles){
    print(peaksfile)
    message("reading peaklist...")
    peaks_temp <- read.table(paste0(peaks_path,peaksfile), header = F)
    print(class(peaks_temp))
    message("making GRanges object...")
    peaks_tempgr <- makeGRangesFromDataFrame(peaks_temp,keep.extra.columns=TRUE, seqnames.field=c("V1"),start.field=c("V2"),end.field=c("V3"))
    storelist[["gr_obj"]][[peaksfile]] <- peaks_tempgr
    message("annotating peaks...")
    print(txdb)
    print(org_db)
    print(head(peaks_tempgr))
    peaks_tempgr_anno <- ChIPseeker::annotatePeak(peaks_tempgr, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb=org_db)
    storelist[["gr_anno"]][[peaksfile]] <- peaks_tempgr_anno
    message("plotting pie, bar, disToTSS...")
    storelist[["plotAnnoPie"]][[peaksfile]] <- plotAnnoPie(peaks_tempgr_anno)
  }
  storelist[["plotAnnoBar"]] <- plotAnnoBar(storelist$gr_anno)
  storelist[["plotDistToTSS"]] <- plotDistToTSS(storelist$gr_anno)
  message("fetaching promoter...")
  promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
  message("generating tagmatrix...")
  storelist[["tagMatrixList"]] <- lapply(storelist[["gr_obj"]], getTagMatrix, windows=promoter)
  message("plotting average plot...")
  storelist[["plotAvgProf"]] <- plotAvgProf(storelist[["tagMatrixList"]], xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
  storelist[["plotPeakProf2"]] <- plotPeakProf2(storelist[["gr_obj"]], upstream = 3000, downstream = 3000, conf = 0.95,by = "gene", type = "start_site", TxDb = txdb, facet = "row", nbin = 800)
  return(storelist)
}

# quantify per feature in a given interval using bedtools coverage
# peaks_path = 
# bamfile = "\\.name\\.bam$"
# bedfile = "_chr.bedpe"
# sites_file = "atac_shControl_merged.txt"

fn_quantify_bedtools <- function(peaks_path, assaytype, bamfile, bedfile ,sites_file){
  bam_files <- list.files(peaks_path, pattern = bamfile, full.names = TRUE)
  for (bam_file in bam_files) {
    cat(paste("Processing BAM file:", bam_file, "\n"))
    output_file <- paste0(sub("\\.name\\.bam$", "", bam_file), "_ac_coverage.pe.bed")
    bedtools_command <- paste("bedtools coverage -a", sites_file, "-b", paste0(bam_file, bedfile), "| sort -k1,1 -k2,2n >", output_file)
    # system(bedtools_command, intern = TRUE)
    print(bedtools_command)
  }
}

# quantify per feature in a given interval using featurecounts
fn_quantify_featurecounts <- function(peaks_path, assaytype, bamfiles, sites_files, sites_type, sites_column_to_rearrange, pairedend=FALSE, refgenome, delim, merge_sites_files=FALSE){
  storelist <- list()
  bam_files <- list.files(peaks_path, pattern = bamfiles, full.names = TRUE)
  sites_files <- list.files(peaks_path, pattern = sites_files, full.names = TRUE)
  storelist[["bams"]] <- bam_files
  for (i in sites_files){
    print(basename(i))
    message("reading the peak file and reaaranging to saf format...")
    sites <- data.frame(fread(i, header = F))[,sites_column_to_rearrange]
    colnames(sites) <- c("id1","id2", "chr", "start", "end", "strand")
    sites["geneid"] <- paste0(sites$id1, delim ,sites$id2)
    saf_file <- sites[,c(7,3,4,5,6)]
    storelist[["saf"]][[basename(i)]] <- saf_file
    message("running featurecounts...")
    storelist[[paste0(sites_type, "_","countmatrix")]][[basename(i)]] <- Rsubread::featureCounts(bam_files, annot.inbuilt = refgenome, annot.ext = saf_file, isGTFAnnotationFile = FALSE, countMultiMappingReads = FALSE,
                            isPairedEnd = pairedend, nthreads = 12) # annot.ext will overide annot.inbuilt if both provided but to be on safe side I supply refgenome
  }
  return(storelist)
}

# quantify for multiple feature using featurecounts
fn_quantify_featurecounts_multifeature <- function(assay_bam_path, feature_path, assay, bam_extension, feature_list, columnstorearrange, ref, control,pe=FALSE, diff=FALSE){
  storelist <- list()
  for (i in names(feature_list)){
    message(paste0("copying ",feature_list[[i]][2], " to ", assay_bam_path, "..."))
    system(paste0("cp ", feature_path, feature_list[[i]][2], " ", assay_bam_path))
    storelist[[i]] <- fn_quantify_featurecounts(assay_bam_path, assay, paste0("\\",bam_extension,"$"), paste0("\\",feature_list[[i]][1],"$"), i, columnstorearrange, pairedend=pe, ref, "%",merge_sites_files=FALSE)
    storelist[[i]][["log_normalized_counts"]] <- log(data.frame(edgeR::cpm(storelist[[i]][[paste0(i, "_countmatrix")]][[feature_list[[i]][2]]][["counts"]])) + 1,2)
    if (diff==TRUE){
      for (j in colnames(storelist[[i]][["log_normalized_counts"]])){
        if (j %notlike% control){
          print(paste0("feature: ",i, ", sample: ",j))
          message(paste0("subtracting ", control, " from ", j,  "..."))
          storelist[[i]][["log_normalized_counts"]][paste0(gsub(bam_extension,"",j),"_",control)] <- storelist[[i]][["log_normalized_counts"]][j] - storelist[[i]][["log_normalized_counts"]][paste0(control,bam_extension)]
        }
      }
    }else{print("Not to calculate differences")}
    storelist[[i]][["log_normalized_counts"]]["id"] <- rownames(storelist[[i]][["log_normalized_counts"]])
  }
  return(storelist)
}


# quantify reads in a given bam files
fn_quantify_bams <- function(peaks_path, assaytype, bamfiles){
  storelist <- list()
  bam_files <- list.files(peaks_path, pattern = bamfiles, full.names = TRUE)
  for (bam_file_path in bam_files) {
    print(basename(bam_file_path))
    bam_file <- Rsamtools::BamFile(bam_file_path)
    message("counting reads...")
    readcounts_stats <- Rsamtools::countBam(bam_file)
    message("storing reads information...")
    storelist[["reads"]][[basename(bam_file_path)]] <- readcounts_stats$records
  }
  message("getting total reads information...")
  reads_df <- do.call(rbind.data.frame, storelist$reads)
  rownames(reads_df) <- names(storelist$reads)
  reads_df["sample"] <- rownames(reads_df)
  colnames(reads_df) <- c("reads", "sample")
  storelist[["total_reads"]] <- reads_df
  return(storelist)
}

# further process the featurecounts quantified matrix
fn_summedreads_per_feature <- function(fearturecountmatrix, bamcountsreads, labels){
  storelist <- list()
  featurelist <- list()
  for (features in names(fearturecountmatrix)){
    message("summing up ", features, "...")
    feature_colsum_df <- data.frame(colSums(fearturecountmatrix[[features]][["counts"]]))
    colnames(feature_colsum_df) <- c("reads")
    feature_colsum_df["sample"] <- rownames(feature_colsum_df)
    feature_colsum_df[["feature"]] <- features
    featurelist[[features]] <- feature_colsum_df
  }
  storelist_df <- do.call(rbind.data.frame, featurelist)
  feature_bam_count_df <- merge(bamcountsreads, storelist_df, by="sample")
  colnames(feature_bam_count_df) <- c("sample", "TotalReads", "ReadsInside", "feature")
  feature_bam_count_df["ReadsOutside"] <- feature_bam_count_df$TotalReads - feature_bam_count_df$ReadsInside
  feature_bam_count_df["Ratio"] <- feature_bam_count_df$ReadsInside / feature_bam_count_df$ReadsOutside
  storelist[["per_features"]] <- featurelist
  rownames(feature_bam_count_df) <- paste0(feature_bam_count_df$sample, "%",feature_bam_count_df$feature)
  storelist[["all_features"]] <- feature_bam_count_df
  labels <- data.frame(fread(labels, header = F))
  colnames(labels) <- c("col", "label")
  feature_bam_count_df_labeled <- merge(labels, feature_bam_count_df, by.x="col", by.y="feature")
  storelist[["all_features_ratio"]] <- feature_bam_count_df_labeled
  return(storelist)
}


# aggregate value per feature
fn_aggregate_feature <- function(annomatrix, columnstoaggregate, aggreagtebycolumn, operation){
  storelist <- list()
  storelist[["aggregated"]] <- stats::aggregate(annomatrix[,columnstoaggregate], by=list(annomatrix[[aggreagtebycolumn]]), operation)
  return(storelist)
}
# integrate multiple omics data on featurecounts
fn_integrate_featurecounts <- function(inputlist, sharedcolumn){
  storelist <- list()
  message("initializing the first dataframe...")
  df <- inputlist[[1]]
  message("merging datasets...")
  for (i in 2:length(inputlist)) {
    print(i)
    message(paste0("merging data... ", i))
    df <- merge(df, inputlist[[i]], by = sharedcolumn, all.x = TRUE, all.y = TRUE) 
  }
  rownames(df) <- df[[sharedcolumn]]
  df <- df[,-1]
  storelist[["df"]] <- df
  return(storelist)
}

# integrate multiple omics data on bedtools coverage output
fn_integrate_bedtools <- function(omics_data_path_lists, binsize, columns_pos, columns_samples){
  storelist <- list()
  omics_data_path_bin <- omics_data_path_lists[grepl(binsize, omics_data_path_lists)]
  bin_list <- list()
  for (files in omics_data_path_bin){
    filename <- sub("\\..*", "", basename(files))
    print(filename)
    bin_list[[filename]] <- data.frame(fread(files, header = F))
  }
  
  bin_df <- do.call(cbind.data.frame, bin_list)
  bin_df <- bin_df[,c(columns_pos, columns_samples)]
  storelist[["df"]] <- bin_df
  rownames(bin_df) <- as.character(apply(bin_df[,columns_pos], 1, function(df) paste0(df, collapse = "%")))
  bin_df <- bin_df[,-columns_pos]
  colnames(bin_df) <- gsub(".V4","", colnames(bin_df))
  
  # remove empty bins
  bin_df_filt <- bin_df[rowSums(bin_df) > 0,]
  bin_df_CPM <- edgeR::cpm(bin_df_filt)
  storelist[["cpm"]] <- bin_df_CPM
  bin_df_CPM_log <- data.frame(log(bin_df_CPM + 1, 2))
  storelist[["cpm_log"]] <- bin_df_CPM_log
  return(storelist)
} 

# make a notlike function
`%notlike%` <- Negate('%like%')

#### common_files among all datasets
# bins5kb
system("bedtools makewindows -g hg38.chrom.sizes -w 5000 | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"\"bin_\"NR\"\\t\"$1\"%\"$2\"%\"$3\"\\t\"\".\"}' > hg38_5kb.txt")
# bins10kb
system("bedtools makewindows -g hg38.chrom.sizes -w 10000 | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"\"bin_\"NR\"\\t\"$1\"%\"$2\"%\"$3\"\\t\"\".\"}' > hg38_10kb.txt")

# gene_extended
integration_features_list <- list()
integration_features_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_gencode_out.sorted.chr.txt")
integration_features_list[["gene_ext"]] <- makeGRangesFromDataFrame(integration_features_list[["gene"]], keep.extra.columns=TRUE)
# since it was complicated to extend only the one side of gene coordinates I extended both sides
integration_features_list$gene_ext <- data.frame(GenomicRanges::shift(integration_features_list$gene_ext, 1500))
write.table(integration_features_list$gene_ext[,c(1,2,3,7,8,5)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/","gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

# promoter
integration_features_list[["promoter"]] <- data.frame(GenomicRanges::promoters(makeGRangesFromDataFrame(integration_features_list[["gene"]], keep.extra.columns=TRUE), upstream=2000, downstream=200))
# filter chrM
integration_features_list$promoter <- integration_features_list$promoter[which(integration_features_list$promoter$seqnames != "chrM"),]
write.table(integration_features_list$promoter[,c(1,2,3,7,8,5)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/","gene_gencodev41_promoter.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)


################################################
###       RNA-Seq data: Knockdown            ###
################################################
rnaseqkd_path <- "/mnt/home3/reid/av638/rnaseq/iva_lab_oct23/"

rnaseqkd_list <- list()

rnaseqkd_list[["dds"]] <- fn_load_rdata(rnaseqkd_path, "outfolder/star_salmon/deseq2_qc/deseq2.dds.RData")
rnaseqkd_list[["gene"]] <- fn_read_genefile(rnaseqkd_path, "gene_gencode_human__gencode_out_chr.txt")
rnaseqkd_list[["counts"]] <- fn_dds_counts(rnaseqkd_list$dds)
rnaseqkd_list[["samples_selected"]] <- fn_sample_filter(rnaseqkd_list$counts, samplefilter = TRUE, c(7:9))
rnaseqkd_list[["processed_counts"]] <- fn_reassign_genenames(rnaseqkd_list$samples_selected, rnaseqkd_list$gene, c(2:10))
rnaseqkd_list[["coldata"]] <- fn_make_coldata_dds(rnaseqkd_list$dds, samplefilter = TRUE, c(7:9))

# check
all(rownames(rnaseqkd_list$coldata) == colnames(rnaseqkd_list$processed_counts)) #should print TRUE

# de analysis
rnaseqkd_list[["deseq2"]] <- fn_deseq2(rnaseqkd_list$processed_counts, rnaseqkd_list$coldata, 0.05, 2)

# Overlap between cramp1 and suz12 target genes
# Venn diagram
list_rnaseqkd_de_shC1_c1509_shS_c1509 <- list(shC1_c1509 = unique(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$de$ensID_Gene), shS_c1509=unique(rnaseqkd_list$deseq2$de_analysis$shS_c1509$de$ensID_Gene))
fn_plot_venn(list_rnaseqkd_de_shC1_c1509_shS_c1509, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_rnaseqkd_up_shC1_c1509_shS_c1509 <- list(shC1_c1509 = unique(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$up$ensID_Gene), shS_c1509=unique(rnaseqkd_list$deseq2$de_analysis$shS_c1509$up$ensID_Gene))
fn_plot_venn(list_rnaseqkd_up_shC1_c1509_shS_c1509, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

# take the rownames or ensID_Gene instead of direct gene symbol, I took rownames which is same as ensID_Gene
set_rnaseqde_uniq_shC1_up <- data.frame(Gene=setdiff(unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shC1_c1509), unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shS_c1509)), label="only_shC1_up")
set_rnaseqde_uniq_shS_up <- data.frame(Gene=setdiff(unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shS_c1509), unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shC1_c1509)), label="only_shS_up")
set_rnaseqde_both_shC1_shS_up <- data.frame(Gene=intersect(unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shC1_c1509), unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shS_c1509)), label="both_shC1_shS_up")
set_rnaseqde_none_up <- data.frame(Gene=setdiff(setdiff(rnaseqkd_list$gene$ens_gene, unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shC1_c1509)), unique(list_rnaseqkd_up_shC1_c1509_shS_c1509$shS_c1509)), label="none_up")
set_rnaseqde_merged_shC1_shS_up <- rbind.data.frame(set_rnaseqde_uniq_shC1_up, set_rnaseqde_uniq_shS_up, set_rnaseqde_both_shC1_shS_up, set_rnaseqde_none_up)


################################################
###       ATAC-Seq data: Knockdown           ###
################################################

# deseq2 with consensus peaks (direct from nextflow)
atacseqkd_deseq2_path <- "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/"

atacseqkd_deseq2_list <- list()

atacseqkd_deseq2_list[["dds"]] <- fn_load_rdata(atacseqkd_deseq2_path, "boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/deseq2/consensus_peaks.mLb.clN.dds.RData")
atacseqkd_deseq2_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/", "gene_gencode_human_gencode_out.sorted.chr.txt")
atacseqkd_deseq2_list[["counts"]] <- fn_dds_counts(atacseqkd_deseq2_list$dds)
atacseqkd_deseq2_list[["samples_selected"]] <- fn_sample_filter(atacseqkd_deseq2_list$counts, samplefilter = FALSE)
atacseqkd_deseq2_list[["processed_counts"]] <- atacseqkd_deseq2_list[["samples_selected"]][,c(1:6)]
atacseqkd_deseq2_list[["coldata"]] <- fn_make_coldata_dds(atacseqkd_deseq2_list$dds, samplefilter = FALSE)
atacseqkd_deseq2_list[["consensus"]] <- data.frame(fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.bed"))
# check
all(rownames(atacseqkd_deseq2_list$coldata) == colnames(atacseqkd_deseq2_list$processed_counts)) #should print TRUE

# de analysis
atacseqkd_deseq2_list[["dar_analysis"]] <- fn_deseq2(atacseqkd_deseq2_list$processed_counts, atacseqkd_deseq2_list$coldata, 0.05, 2)

# get positions of dars
atacseqkd_deseq2_list[["annotation"]] <- fn_consensus_gene_anno("atacseqkd", "deseq2", atacseqkd_deseq2_list$dar_analysis$de_analysis, atacseqkd_deseq2_list$consensus, atacseqkd_deseq2_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_deseq2_de_shCRAMP1_shSUZ12_shControl <- list(shCRAMP1 = unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$de$ensID_Gene), shSUZ12=unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$de$ensID_Gene))
fn_plot_venn(list_atacseqkd_deseq2_de_shCRAMP1_shSUZ12_shControl, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_atacseqkd_deseq2_up_shCRAMP1_shSUZ12_shControl <- list(shCRAMP1 = unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$up$ensID_Gene), shSUZ12=unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$up$ensID_Gene))
fn_plot_venn(list_atacseqkd_deseq2_up_shCRAMP1_shSUZ12_shControl, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotPCA(atacseqkd_deseq2_list$dar_analysis$vstO, intgroup="sample", returnData=FALSE)


atacseqkd_deseq2_list[["aggregate_pergene"]][["shCRAMP1_shControl"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_deseq2_list$annotation$shCRAMP1_shControl$nearest$de, c(3), "ens_gene", mean)
atacseqkd_deseq2_list[["aggregate_pergene"]][["shSUZ12_shControl"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_deseq2_list$annotation$shSUZ12_shControl$nearest$de, c(3), "ens_gene", mean)

colnames(atacseqkd_deseq2_list$aggregate_pergene$shCRAMP1_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(atacseqkd_deseq2_list$aggregate_pergene$shSUZ12_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")

# diffbind
diffbind_atacseqkd_samplesheet <- read.csv("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/diffbind_atac_seq_samplesheet.csv")
atacseqkd_diffbind_list <- list()
atacseqkd_diffbind_list[["prep"]] <- fn_diffbind_prep(diffbind_atacseqkd_samplesheet)
plot(atacseqkd_diffbind_list$prep$dba_obj)
plot(atacseqkd_diffbind_list$prep$ovp_rate,type ='b',ylab='# peaks', xlab='Overlap at least this many peaksets')

atacseqkd_diffbind_list[["counts"]] <- fn_diffbind_count(atacseqkd_diffbind_list$prep$dba_obj, 100)

atacseqkd_diffbind_list[["dar_analysis"]] <- fn_diffbind_de(atacseqkd_diffbind_list$counts$dba_obj, 0.05, 2)

atacseqkd_diffbind_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/", "gene_gencode_human_gencode_out.sorted.chr.txt")

atacseqkd_diffbind_list[["annotation"]] <- fn_nearest_gene_anno("atacseqkd", "diffbind",atacseqkd_diffbind_list$dar_analysis$contrasts, atacseqkd_diffbind_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

# limma
atacseqkd_limma_list <- list()
atacseqkd_counts <- data.frame(atacseqkd_diffbind_list$counts$dba_obj$binding)
rownames(atacseqkd_counts) <- paste0("chr",atacseqkd_counts$CHR, "%",atacseqkd_counts$START,"%", atacseqkd_counts$END)
atacseqkd_limma_list[["counts"]] <- atacseqkd_counts[,c(4:9)]
atacseqkd_limma_list[["coldata"]] <- atacseqkd_diffbind_list$counts$info
atacseqkd_limma_list$coldata$Condition <- factor(atacseqkd_limma_list$coldata$Condition)
atacseqkd_limma_list[["dge_obj"]] <- fn_limma_create(atacseqkd_limma_list$counts, atacseqkd_limma_list$coldata)

# create contrasts
atacseqkd_limma_list[["contrasts"]] <- limma::makeContrasts(paste0(colnames(atacseqkd_limma_list$dge_obj$design)[c(2,1)], collapse = "_vs_"),
                                          paste0(colnames(atacseqkd_limma_list$dge_obj$design)[c(3,1)], collapse = "-"), levels = atacseqkd_limma_list$dge_obj$design)

atacseqkd_limma_list[["dar_analysis"]] <- fn_limma_de(atacseqkd_limma_list$dge_obj$fit, atacseqkd_limma_list$contrasts, 0.05, 2)

atacseqkd_limma_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/", "gene_gencode_human_gencode_out.sorted.chr.txt")

atacseqkd_limma_list[["annotation"]] <- fn_nearest_gene_anno("atacseqkd", "limma", atacseqkd_limma_list$dar_analysis$dars, atacseqkd_limma_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)


# edgeR analysis
atacseqkd_edger_list <- list()
atacseqkd_counts <- data.frame(atacseqkd_diffbind_list$counts$dba_obj$binding)
rownames(atacseqkd_counts) <- paste0("chr",atacseqkd_counts$CHR, "%",atacseqkd_counts$START,"%", atacseqkd_counts$END)
atacseqkd_edger_list[["counts"]] <- atacseqkd_counts[,c(4:9)]
atacseqkd_edger_list[["coldata"]] <- atacseqkd_diffbind_list$counts$info
atacseqkd_edger_list$coldata$Condition <- factor(atacseqkd_edger_list$coldata$Condition)
atacseqkd_edger_list[["dge_obj"]] <- fn_edger_create(atacseqkd_edger_list$counts, atacseqkd_edger_list$coldata)

atacseqkd_edger_list[["dar_analysis"]] <- fn_edger_de(atacseqkd_edger_list$dge_obj$lmfit, atacseqkd_edger_list$dge_obj$design, 0.05, 2)

atacseqkd_edger_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/", "gene_gencode_human_gencode_out.sorted.chr.txt")

atacseqkd_edger_list[["annotation"]] <- fn_nearest_gene_anno("atacseqkd", "edger", atacseqkd_edger_list$dar_analysis$dars, atacseqkd_edger_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

atacseqkd_edger_list[["chipseeker"]] <- fn_run_chipseeker("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/", "atacseqkd", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "*.broadPeak")
  
#  export bed files
fn_export_bed("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", atacseqkd_edger_list$annotation, "all", "all", "txt")


################################################
###         CUT&RUN data: Knockout           ###
################################################

# deseq2 with consensus peaks (direct from nextflow)
cutnrunko_deseq2_path <- "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/"

cutnrunko_deseq2_list <- list()

cutnrunko_deseq2_list[["dds"]] <- fn_load_rdata(cutnrunko_deseq2_path, "outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/deseq2/H3K27me3.consensus_peaks.dds.RData")
cutnrunko_deseq2_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/", "gene_gencode_human_gencode_out.sorted.chr.txt")
cutnrunko_deseq2_list[["counts"]] <- fn_dds_counts(cutnrunko_deseq2_list$dds)
cutnrunko_deseq2_list[["samples_selected"]] <- fn_sample_filter(cutnrunko_deseq2_list$counts, samplefilter = FALSE)
cutnrunko_deseq2_list[["processed_counts"]] <- cutnrunko_deseq2_list[["samples_selected"]][,c(1:4)]
cutnrunko_deseq2_list[["coldata"]] <- data.frame(colData(cutnrunko_deseq2_list$dds))
cutnrunko_deseq2_list$coldata <- data.frame(sample=cutnrunko_deseq2_list$coldata$sample, condition=factor(gsub('[0-9]+', '', unlist(lapply(strsplit(cutnrunko_deseq2_list$coldata$sample, "_"), function(x) x[2])))), replicate=factor(c(3,1,1,2)), sizeFactors=cutnrunko_deseq2_list$coldata$sizeFactor)
rownames(cutnrunko_deseq2_list$coldata) <- cutnrunko_deseq2_list$coldata$sample
cutnrunko_deseq2_list[["consensus"]] <- data.frame(fread("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/H3K27me3.consensus_peaks.bed"))
# check
all(rownames(cutnrunko_deseq2_list$coldata) == colnames(cutnrunko_deseq2_list$processed_counts)) #should print TRUE

# de analysis
cutnrunko_deseq2_list[["dar_analysis"]] <- fn_deseq2(cutnrunko_deseq2_list$processed_counts, cutnrunko_deseq2_list$coldata, 0.05, 2)

# get positions of dars
cutnrunko_deseq2_list[["annotation"]] <- fn_consensus_gene_anno("cutnrunko", "deseq2", cutnrunko_deseq2_list$dar_analysis$de_analysis, cutnrunko_deseq2_list$consensus, cutnrunko_deseq2_list$gene, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/", "txt", 0.05, 2)

plotPCA(cutnrunko_deseq2_list$dar_analysis$vstO, intgroup="sample", returnData=FALSE)

cutnrunko_deseq2_list[["aggregate_pergene"]][["KO_WT"]][["nearest"]][["de"]] <- fn_aggregate_feature(cutnrunko_deseq2_list$annotation$KO_WT$nearest$de, c(3), "ens_gene", mean)

colnames(cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")


cutnrunko_quantify_list <- list()

# quantify over various features
cutnrunko_quantify_list[["featurecounts"]][["featureslist"]]=list(gene_ext=c(".chr_1.5kb.txt", "gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), bins5kb=c("hg38_5kb.txt", "hg38_5kb.txt"), bins10kb=c("hg38_10kb.txt", "hg38_10kb.txt"), promoter=c("_gencodev41_promoter.txt", "gene_gencodev41_promoter.txt"))

cutnrunko_quantify_list[["featurecounts"]][["featuresmatrix"]] <- fn_quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "cutnrunko", ".mLb.clN.sorted.bam", cutnrunko_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=FALSE, diff=FALSE)

################################################
###         CUT&RUN data: Knockdown          ###
################################################

# deseq2 with consensus peaks (direct from nextflow)
cutnrunkd_deseq2_path <- "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/"

cutnrunkd_deseq2_list <- list()

cutnrunkd_deseq2_list[["dds"]] <- fn_load_rdata(cutnrunkd_deseq2_path, "outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/deseq2/H3K27me3.consensus_peaks.dds.RData")
cutnrunkd_deseq2_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/", "gene_gencode_human_gencode_out.sorted.chr.txt")
cutnrunkd_deseq2_list[["counts"]] <- fn_dds_counts(cutnrunkd_deseq2_list$dds)
cutnrunkd_deseq2_list[["samples_selected"]] <- fn_sample_filter(cutnrunkd_deseq2_list$counts, samplefilter = FALSE)
cutnrunkd_deseq2_list[["processed_counts"]] <- cutnrunkd_deseq2_list[["samples_selected"]][,c(1:3)]
cutnrunkd_deseq2_list[["consensus"]] <- data.frame(fread("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/H3K27me3.consensus_peaks.bed"))
cutnrunkd_deseq2_list[["normalized_counts"]] <- data.frame(edgeR::cpm(cutnrunkd_deseq2_list$processed_counts))
cutnrunkd_deseq2_list$normalized_counts["ensID_Gene"] <- rownames(cutnrunkd_deseq2_list$normalized_counts)

# get positions of peaks
cutnrunkd_deseq2_list[["annotation"]] <- fn_consensus_gene_count_anno("cutnrunko", "deseq2", cutnrunkd_deseq2_list$normalized_counts, cutnrunkd_deseq2_list$consensus, cutnrunkd_deseq2_list$gene, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/", "txt")


# pca 
cutnrunkd_deseq2_list[["vst"]] <- DESeq2::vst(cutnrunkd_deseq2_list$dds)
plotPCA(cutnrunkd_deseq2_list$vst, intgroup="sample", returnData=FALSE)

cutnrunkd_deseq2_list[["aggregate_pergene"]][["nearest"]] <- fn_aggregate_feature(cutnrunkd_deseq2_list$annotation$nearest, c(2:4), "ens_gene", sum)

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

################################################
###         CUT&TAG data: Wildtype  SE       ###
################################################
cutntagwt_quantify_list <- list()
cutntagwt_quantify_list[["featurecounts"]] <- fn_quantify_featurecounts("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary", "cutntagwt","\\.bam$", "\\_peaks_id.bed$","histone_marks", c(12,4,1:3,7), pairedend=FALSE, "hg38", "_",merge_sites_files=FALSE)

cutntagwt_quantify_list[["bams"]] <- fn_quantify_bams("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary", "cutntagwt","\\.bam$")
cutntagwt_labelpath = "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/histone_marks_label.txt"
cutntagwt_quantify_list[["summed"]] <- fn_summedreads_per_feature(cutntagwt_quantify_list$featurecounts$histone_marks_countmatrix, cutntagwt_quantify_list$bams$total_reads, cutntagwt_labelpath)

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
cutntagwt_quantify_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/", "gene_gencode_human_gencode_out.sorted.chr.txt")

# quantify over various features
cutntagwt_quantify_list[["featurecounts"]][["featureslist"]]=list(gene_ext=c(".chr_1.5kb.txt", "gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), bins5kb=c("hg38_5kb.txt", "hg38_5kb.txt"), bins10kb=c("hg38_10kb.txt", "hg38_10kb.txt"), promoter=c("_gencodev41_promoter.txt", "gene_gencodev41_promoter.txt"))

cutntagwt_quantify_list[["featurecounts"]][["featuresmatrix"]] <- fn_quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "cutntagwt", ".mLb.clN.sorted.bam", cutntagwt_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=FALSE, diff=TRUE)

# for check if two df are exactly same and to check if fn_quantify_featurecounts_multifeature is created properly and doing the same thing as individual steps look at Check Notes

################################################
###         CUT&TAG data: wt and ko          ###
################################################
# ls *.fastq.gz -1 |  grep R1 | sort -k1,1 | awk -F'_' '{print $2"\t"$1"_"$2"_"$3"_"$4"_""R1""_"$6"\t"$1"_"$2"_"$3"_"$4"_""R2""_"$6}' > fastq_info.txt
# add header vim fastq_info.txt
# samplesheet prep for nextflow
# Rachael told me control for the FLAG labelled samples is "FLAG_H14_old".
# The one labelled as FLAG_old is actually the FLAG_H1.4_old sample. So I changed it
cutntagwko_quantify_list <- list()
cutntagwko_quantify_list[["nextflow"]][["info"]] <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_dec23/info.txt", header = T))
cutntagwko_quantify_list[["nextflow"]][["sheet"]] <- data.frame(ID=rep(c(cutntagwko_quantify_list$nextflow$info$ID),2),
                replicate=rep(1,length(cutntagwko_quantify_list$nextflow$info$ID) * 2),
                group=gsub("-","_",rep(c(cutntagwko_quantify_list$nextflow$info$Sample),2)),
                control=rep(c("shControl_IgG","shControl_IgG","","shCRAMP1_IgG","shCRAMP1_IgG","","shSUZ12_IgG","shSUZ12_IgG","","WT_IgG","WT_IgG","WT_IgG","WT_IgG","WT_IgG","","FLAG_H14_old","FLAG_H14_old","FLAG_H14_old","FLAG_H14_old","FLAG_H14_old","",""),2)
)

cutntagwko_quantify_list[["nextflow"]][["fastqs"]] <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_dec23/fastq_info.txt", header = T))
cutntagwko_quantify_list[["nextflow"]][["samplesheet"]] <- merge(cutntagwko_quantify_list$nextflow$sheet, cutntagwko_quantify_list$nextflow$fastqs, by="ID")
cutntagwko_quantify_list$nextflow$samplesheet <- cutntagwko_quantify_list$nextflow$samplesheet %>% dplyr::distinct()
cutntagwko_quantify_list$nextflow$samplesheet <- cutntagwko_quantify_list$nextflow$samplesheet[order(cutntagwko_quantify_list$nextflow$samplesheet$group),][,c(3,2,5,6,4)]
write.table(cutntagwko_quantify_list$nextflow$samplesheet, "/mnt/home3/reid/av638/cutntag/iva_lab_dec23/samplesheet.csv", sep=",", col.names = T, row.names = F, quote = F, append = F)

# Integration of omics data
# for gene
data_integration_list <- list()
data_integration_list[["ensID_gene"]][["list"]] <- 
  list(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$de %>% select(r_shC1=2, ensID_Gene=7),
       rnaseqkd_list$deseq2$de_analysis$shS_c1509$de %>% select(r_shS=2, ensID_Gene=7),
       atacseqkd_deseq2_list$aggregate_pergene$shCRAMP1_shControl$nearest$de$aggregated %>% select(a_shCRAMP1=2, ensID_Gene=1),
       atacseqkd_deseq2_list$aggregate_pergene$shSUZ12_shControl$nearest$de$aggregated %>% select(a_shSUZ12=2, ensID_Gene=1),
       cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest$de$aggregated %>% select(cko_KO=2, ensID_Gene=1),
       cutnrunkd_deseq2_list$aggregate_pergene$nearest$log_sum_counts[,c(4:6)],
       cutntagwt_quantify_list$featurecounts$gene_ext$log_normalized_counts[,c(7:10)])

data_integration_list[["ensID_gene"]][["merged"]] <- fn_integrate_featurecounts(data_integration_list$ensID_gene$list, "ensID_Gene")

# Some predictable
col_fun = circlize::colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
col_fun2 = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
data_integration_list$ensID_gene$merged[["somepredictable"]] <- data_integration_list$ensID_gene$merged$df[complete.cases(data_integration_list$ensID_gene$merged$df[, c("r_shC1")]), ]
data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat1"]] <- ComplexHeatmap::Heatmap(as.matrix(data_integration_list$ensID_gene$merged$somepredictable[,c(1:5)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, col = col_fun)
data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat2"]] <- ComplexHeatmap::Heatmap(as.matrix(data_integration_list$ensID_gene$merged$somepredictable[,c(6:7)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, col = col_fun2)
data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat3"]] <- ComplexHeatmap::Heatmap(as.matrix(data_integration_list$ensID_gene$merged$somepredictable[,c(8:10)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, col = col_fun)

data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat1"]]+ data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat2"]] + data_integration_list$ensID_gene$merged[["Heatmap"]][["somepredictable"]][["mat3"]]

#complex_heatmap_data_integration_list$ensID_gene$merged$df_somepredictable_ch.pdf
# for bins
intergration_omics_list <- list()
bin_cov_atacseq_path <- list.files(path="/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library", pattern = "*.txt_coverage.pe.bed", full.names = T)
bin_cov_cutnrun_KO_path <- list.files(path="/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary", pattern = "*.txt_coverage_se.bed", full.names = T)
bin_cov_cutnrun_KD_path <- list.files(path="/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary", pattern = "*.txt_coverage_se.bed", full.names = T)
bin_cov_cutntag_path <- list.files(path="/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary", pattern = "*.txt_coverage_se.bed", full.names = T)
bin_path_combined <- c(bin_cov_atacseq_path, bin_cov_cutnrun_KO_path, bin_cov_cutnrun_KD_path, bin_cov_cutntag_path)

intergration_omics_list[["bin5kb"]] <- fn_integrate_bedtools(bin_path_combined, "5kb", c(1:3), seq(4,175,7))
intergration_omics_list[["bin10kb"]] <- fn_integrate_bedtools(bin_path_combined, "10kb", c(1:3), seq(4,175,7))






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