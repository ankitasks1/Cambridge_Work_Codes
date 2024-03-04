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
library(RColorBrewer)


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
        assayseqde["feature_id"] <- rownames(assayseqde)
        assayseqdedf <- data.frame(assayseqde)
        assayseqdedf["feature_id"] <- rownames(assayseqdedf)
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
                                                       genenames = unlist(lapply(strsplit(assayseqdedf$feature_id, "%"), function(x) x[2])), legend = "top", top = 20, font.label = c("bold", 5),
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
    delist$all["feature_id"] <- rownames(delist$all)
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

fn_diffbind_norm <- function(obj){
  diffblist <- list()
  message("Performing normalisation...")
  diffblist[["norm"]] <- dba.normalize(obj, method=DBA_ALL_METHODS)
  message("Performing Diffbind analysis...")
  message("Performing Contrasts...")
  diffblist$norm <- dba.contrast(diffblist$norm, minMembers = 2,categories=DBA_CONDITION)
  diffblist$norm <- dba.analyze(diffblist$norm, method=DBA_ALL_METHODS)
  return(diffblist)
}

fn_diffbind_de <- function(normobj, fdr, fc, method){
  contrastlist <- list()
  diffblist <- list()
  contrasts <- dba.show(normobj, bContrasts=TRUE)
  for (i in as.numeric(rownames(contrasts))){
    print(paste0("performing analysis for contrast",i, " i.e. ", paste0(as.character(contrasts[i,][,c(4,2)]), collapse = "_")))
    storelist <- list()
    message("Fetching DE info...")
    # Extracting results, contrast 1 means , th means threshold for FDR, which if 1 give all sites 
    if(method=="edgeR"){
      storelist <- list()
      storelist[["all"]] <- data.frame(dba.report(normobj, method=DBA_EDGER, contrast = i, th=1))
      storelist$all["feature_id"] <- unique(do.call(paste, c(storelist$all[,c(1:3)], sep="%")))
      storelist[["de"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2) | Fold < -log(fc,2)))
      storelist[["up"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2)))
      storelist[["down"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold < -log(fc,2)))
      message("Performing PCA plot...")
      storelist[["pca"]] <- dba.plotPCA(normobj, contrast=i, method=DBA_EDGER, attributes=DBA_FACTOR, label=DBA_ID)
      message("Performing venn diagran...")
      storelist[["venn"]] <- dba.plotVenn(normobj, contrast=i, method=DBA_ALL_METHODS)
      message("Performing MA-plot...")
      storelist[["ma_plot"]] <- dba.plotMA(normobj, method=DBA_EDGER)
    }else if(method=="deseq2"){
      storelist <- list()
      storelist[["all"]] <- data.frame(dba.report(normobj, method=DBA_DESEQ2, contrast = i, th=1))
      storelist$all["feature_id"] <- unique(do.call(paste, c(storelist$all[,c(1:3)], sep="%")))
      storelist[["de"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2) | Fold < -log(fc,2)))
      storelist[["up"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2)))
      storelist[["down"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold < -log(fc,2)))
      message("Performing PCA plot...")
      storelist[["pca"]] <- dba.plotPCA(normobj, contrast=i, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)
      message("Performing venn diagran...")
      storelist[["venn"]] <- dba.plotVenn(normobj, contrast=i, method=DBA_ALL_METHODS)
      message("Performing MA-plot...")
      storelist[["ma_plot"]] <- dba.plotMA(normobj, method=DBA_DESEQ2)
    }
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
    delist$all["feature_id"] <- rownames(delist$all)
    delist[["de"]] <- delist$all %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2) | logFC < -log(fc,2)))
    delist[["up"]] <- delist$all %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2)))
    delist[["down"]] <- delist$all %>% dplyr::filter(adj.P.Val < fdr & (logFC < -log(fc,2)))
    storelist[["dars"]][[paste0("coef_", gsub("-","_",cont))]] <- delist
  }
  return(storelist)
}

# annotation to genes and nearest distance filter for all packages
fn_diffbind_gene_anno <- function(assaytype, software, contrasts, genefile, path, outformat, fdr, fc){
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
    nearest_distance = 0
    all_anno_nearest <- all_anno %>% dplyr::filter(Distance == nearest_distance)
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
    }
  }
  return(contrastlist)
}

fn_diffbind_and_packages_gene_anno <- function(assaytype, software, contrasts, genefile, path, outformat, fdr, fc){
  contrastlist <- list()
  for (i in names(contrasts)){
    print(i)
    all_df <- data.frame(contrasts[[i]][["all"]])
    message("adding coordinates...")
    all_df["coordinates"] <- rownames(all_df)
    all_df <- data.frame(splitstackshape::cSplit(all_df, "coordinates", "%"))
    all_df <- all_df[,c((dim(all_df)[2] -2): (dim(all_df)[2]), 1:(dim(all_df)[2]-3))]
    
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
    nearest_distance = 0
    all_anno_nearest <- all_anno %>% dplyr::filter(Distance == nearest_distance)
    storelist[["nearest"]][["all"]] <- all_anno_nearest
    if (software == "diffbind_limma") {
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2) | logFC < -log(fc,2))) 
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(adj.P.Val < fdr & (logFC < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2) | logFC < -log(fc,2))) 
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(adj.P.Val < fdr & (logFC < -log(fc,2)))
      contrastlist[[i]] <- storelist
    }else if (software == "diffbind_edger") {
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2) | logFC < -log(fc,2))) 
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2) | logFC < -log(fc,2))) 
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC < -log(fc,2)))
      contrastlist[[i]] <- storelist
    }else if (software == "diffbind_deseq2") {
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
    all_anno <- merge(all_df, consensus_to_gene, by.x="feature_id", by.y="peak_interval")
    storelist <- list()
    storelist[["all"]][["all"]] <- all_anno
    
    message("getting peaks nearest distance to genes...")
    nearest_distance = 0
    all_anno_nearest <- all_anno %>% dplyr::filter(Distance == nearest_distance)
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
    }else if (software == "edger") {
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2) | logFC < -log(fc,2))) 
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2) | logFC < -log(fc,2))) 
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC < -log(fc,2)))
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
  nearest_distance = 0
  all_anno_nearest <- all_anno %>% dplyr::filter(Distance == nearest_distance)
  storelist[["nearest"]] <- all_anno_nearest
  return(storelist)
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
    storelist[["anno_df"]][[peaksfile]] <- data.frame(peaks_tempgr_anno)
    message("plotting pie, bar, disToTSS...")
    storelist[["plotAnnoPie"]][[peaksfile]] <- ChIPseeker::plotAnnoPie(peaks_tempgr_anno)
  }
  storelist[["plotAnnoBar"]] <- ChIPseeker::plotAnnoBar(storelist$gr_anno)
  storelist[["plotDistToTSS"]] <- ChIPseeker::plotDistToTSS(storelist$gr_anno)
  message("fetaching promoter...")
  promoter <- ChIPseeker::getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
  message("generating tagmatrix...")
  storelist[["tagMatrixList"]] <- lapply(storelist[["gr_obj"]], ChIPseeker::getTagMatrix, windows=promoter)
  message("plotting average plot...")
  storelist[["plotAvgProf"]] <- ChIPseeker::plotAvgProf(storelist[["tagMatrixList"]], xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
  storelist[["plotPeakProf2"]] <- ChIPseeker::plotPeakProf2(storelist[["gr_obj"]], upstream = 3000, downstream = 3000, conf = 0.95,by = "gene", type = "start_site", TxDb = txdb, facet = "row", nbin = 800)
  return(storelist)
}

fn_chipseeker_region_anno <- function(assaytype, software, contrasts, region_bed, path, txdb, org_db, regionformat, fdr, fc){
  message("writing region bed data...", regionformat)
  write.table(region_bed, paste0(path, "/", assaytype, "_chipseeker_", software, "_", regionformat), sep="\t", quote = F, append = F, row.names = F, col.names = F)
  regionlist <- list()
  regionfiles <- list.files(path=path, pattern = paste0("_chipseeker_", software, "_", regionformat))
  for (regionfile in regionfiles){
    message("Region file in analysis: ",regionfile)
    message("reading region file...")
    region_temp <- read.table(paste0(path, regionfile), header = F)
    message("making GRanges object...")
    region_tempgr <- makeGRangesFromDataFrame(region_temp,keep.extra.columns=TRUE, seqnames.field=c("V1"),start.field=c("V2"),end.field=c("V3"))
    regionlist[["gr_obj"]][[regionfile]] <- region_tempgr
    message("annotating peaks...")
    region_tempgr_anno <- ChIPseeker::annotatePeak(region_tempgr, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb=org_db)
    regionlist[["gr_anno"]][[regionfile]] <- region_tempgr_anno
    regionlist[["anno_df"]][[regionfile]] <- data.frame(region_tempgr_anno)
  }
  region_to_gene <- regionlist[["anno_df"]][[regionfile]]
  contrastlist <- list() # create empty list to store ALL contrast data
  contrastlist[["region"]] <- regionlist
  for (i in names(contrasts)){
    print(i)
    all_df <- data.frame(contrasts[[i]][["all"]])
    message("annotation to region bed file ...")
    all_anno <- merge(all_df, region_to_gene, by.x="feature_id", by.y="V4")
    storelist <- list() # create empty list to store EACH contrast data
    storelist[["all"]][["all"]] <- all_anno

    message("getting peaks nearest distance to genes...")
    upstreamFromTSS = 2000
    downstreamFromTSS = 100
    all_anno_nearest <- all_anno %>% dplyr::filter(distanceToTSS < downstreamFromTSS & distanceToTSS > -upstreamFromTSS)
    storelist[["nearest"]][["all"]] <- all_anno_nearest
    if (software == "deseq2"){
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2) | log2FoldChange < -log(fc,2)))
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(padj < fdr & (log2FoldChange < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2) | log2FoldChange < -log(fc,2)))
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(padj < fdr & (log2FoldChange > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(padj < fdr & (log2FoldChange < -log(fc,2)))
      contrastlist[["contrasts"]][[i]] <- storelist
    }else if (software == "edger") {
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2) | logFC < -log(fc,2)))
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (logFC < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2) | logFC < -log(fc,2)))
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(FDR < fdr & (logFC < -log(fc,2)))
      contrastlist[["contrasts"]][[i]] <- storelist
    }else if (software == "limma") {
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2) | logFC < -log(fc,2)))
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(adj.P.Val < fdr & (logFC < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2) | logFC < -log(fc,2)))
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(adj.P.Val < fdr & (logFC > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(adj.P.Val < fdr & (logFC < -log(fc,2)))
      contrastlist[["contrasts"]][[i]] <- storelist
    }else if (software == "diffbind"){
      message(paste0("package used: ", software))
      storelist[["nearest"]][["de"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2) | Fold < -log(fc,2)))
      storelist[["nearest"]][["up"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2)))
      storelist[["nearest"]][["down"]] <- all_anno_nearest %>% dplyr::filter(FDR < fdr & (Fold < -log(fc,2)))
      storelist[["all"]][["de"]] <- all_anno %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2) | Fold < -log(fc,2)))
      storelist[["all"]][["up"]] <- all_anno %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2)))
      storelist[["all"]][["down"]] <- all_anno %>% dplyr::filter(FDR < fdr & (Fold < -log(fc,2)))
      contrastlist[["contrasts"]][[i]] <- storelist
    }
  }
  return(contrastlist)
}
# names=c("chr", "start", "end", "feature_id")
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

# quantify for multiple feature using featurecounts subtract
fn_quantify_featurecounts_multifeature <- function(assay_bam_path, feature_path, assay, bam_extension, feature_list, columnstorearrange, ref, control,pe=FALSE, diff="single", pairs){
  storelist <- list()
  for (i in names(feature_list)){
    message(paste0("copying ",feature_list[[i]][2], " to ", assay_bam_path, "..."))
    system(paste0("cp ", feature_path, feature_list[[i]][2], " ", assay_bam_path))
    message("performing featurecounts...")
    storelist[[i]] <- fn_quantify_featurecounts(assay_bam_path, assay, paste0("\\",bam_extension,"$"), paste0("\\",feature_list[[i]][1],"$"), i, columnstorearrange, pairedend=pe, ref, "%",merge_sites_files=FALSE)
    storelist[[i]][["log_normalized_counts"]] <- log(data.frame(edgeR::cpm(storelist[[i]][[paste0(i, "_countmatrix")]][[feature_list[[i]][2]]][["counts"]])) + 1,2)
    message("executing nomralization with respective controls...")
    if (diff=="single"){
      for (j in colnames(storelist[[i]][["log_normalized_counts"]])){
        if (j %notlike% control){
          print(paste0("feature: ",i, ", sample: ",j))
          message(paste0("subtracting ", control, " from ", j,  "..."))
          storelist[[i]][["log_normalized_counts"]][paste0(gsub(bam_extension,"",j),"_",control)] <- storelist[[i]][["log_normalized_counts"]][j] - storelist[[i]][["log_normalized_counts"]][paste0(control,bam_extension)]
        }
      }
    }else if (diff=="multi"){
      for (k in names(pairs)){
          print(paste0("feature: ",i, ", sample: ", pairs[[k]][1], " control: ", pairs[[k]][2]))
          message(paste0("subtracting ", pairs[[k]][2], " from ", pairs[[k]][1],  "..."))
          storelist[[i]][["log_normalized_counts"]][k] <- storelist[[i]][["log_normalized_counts"]][pairs[[k]][1]] - storelist[[i]][["log_normalized_counts"]][pairs[[k]][2]]

        }
      }
    storelist[[i]][["log_normalized_counts"]]["id"] <- rownames(storelist[[i]][["log_normalized_counts"]])
  }
  return(storelist)
}
# quantify for multiple feature using featurecounts_div
fn_quantify_featurecounts_multifeature_div <- function(assay_bam_path, feature_path, assay, bam_extension, feature_list, columnstorearrange, ref, control,pe=FALSE, diff="single", pairs){
  storelist <- list()
  for (i in names(feature_list)){
    message(paste0("copying ",feature_list[[i]][2], " to ", assay_bam_path, "..."))
    system(paste0("cp ", feature_path, feature_list[[i]][2], " ", assay_bam_path))
    message("performing featurecounts...")
    storelist[[i]] <- fn_quantify_featurecounts(assay_bam_path, assay, paste0("\\",bam_extension,"$"), paste0("\\",feature_list[[i]][1],"$"), i, columnstorearrange, pairedend=pe, ref, "%",merge_sites_files=FALSE)
    storelist[[i]][["normalized_counts"]] <- data.frame(edgeR::cpm(storelist[[i]][[paste0(i, "_countmatrix")]][[feature_list[[i]][2]]][["counts"]])) 
    storelist[[i]][["log_normalized_counts"]] <- log(data.frame(edgeR::cpm(storelist[[i]][[paste0(i, "_countmatrix")]][[feature_list[[i]][2]]][["counts"]])) + 1,2)
    message("executing nomralization with respective controls...")
    if (diff=="single"){
      for (j in colnames(storelist[[i]][["normalized_counts"]])){
        if (j %notlike% control){
          print(paste0("feature: ",i, ", sample: ",j))
          message(paste0("dividing ", control, " from ", j,  "..."))
          storelist[[i]][["log_normalized_counts"]][paste0(gsub(bam_extension,"",j),"_",control)] <- log(((storelist[[i]][["normalized_counts"]][j] + 1) / (storelist[[i]][["normalized_counts"]][paste0(control,bam_extension)] + 1)),2)
        }
      }
    }else if (diff=="multi"){
      for (k in names(pairs)){
        print(paste0("feature: ",i, ", sample: ", pairs[[k]][1], " control: ", pairs[[k]][2]))
        message(paste0("dividing ", pairs[[k]][2], " from ", pairs[[k]][1],  "..."))
        storelist[[i]][["log_normalized_counts"]][k] <- log(((storelist[[i]][["normalized_counts"]][pairs[[k]][1]] + 1) / (storelist[[i]][["normalized_counts"]][pairs[[k]][2]] + 1)), 2)
      }
    }
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
fn_aggregate_feature <- function(annomatrix, columnstoaggregate, aggregatebycolumn, operation){
  storelist <- list()
  message("aggregating...")
  storelist[["aggregated"]] <- stats::aggregate(annomatrix[,columnstoaggregate], by=list(annomatrix[[aggregatebycolumn]]), operation)
  return(storelist)
}

fn_make_pca <- function(counts, file_extension, attribute){
  storelist <- list()
  message("filtering the counts for zero rowsums...")
  storelist[["norm_counts_filtered"]] <- counts[rowSums(counts) > 0,]
  colnames(storelist$norm_counts_filtered) <- gsub(file_extension,"", colnames(storelist$norm_counts_filtered))
  message("generating pca")
  storelist[["prcomp"]] <- prcomp(t(storelist$norm_counts_filtered), center = TRUE, scale. = TRUE)
  storelist[[attribute]] <- factoextra::fviz_pca_ind(storelist$prcomp, geom = c("point", "text"), repel = TRUE)
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

source("/mnt/beegfs6/home3/reid/av638/atacseq/iva_lab_gencode/integration/upset_out.R")

# general MA plot using ggplot2
fn_maplot <- function(df, columnbyorder, fdr, fc){
  plotlist <- list()
  message("rearranging and plotting ma...")
  df <- df[,columnbyorder]
  colnames(df) <- c("baseMean","log2FoldChange","padj","feature")
  plotlist[["maplot"]] <- ggmaplot(df, fdr = fdr, fc = fc, size = 0.3, palette = c("#D22B2B", "#1465AC", "darkgray"),
           genenames = df$feature, legend = "top", top = 20,
           font.label = c("bold", 5), font.legend = "bold", font.main = "bold", ggtheme = ggplot2::theme_minimal())
}

fn_maplot_general <- function(df, columnbyorder, fdr, logfc){
  plotlist <- list()
  message("rearranging and plotting ma...")
  df <- df[,columnbyorder]
  colnames(df) <- c("baseMean","log2FoldChange","padj","feature")
  df$color <- ifelse(df$log2FoldChange > logfc & df$padj < fdr, "#B31B21", 
                     ifelse(df$log2FoldChange < -logfc & df$padj < fdr, "#1465AC", "grey"))
  plotlist[["maplot"]] <- ggplot(df, aes(x = baseMean, y = log2FoldChange, color = color)) + geom_point(alpha = 0.8, size = 0.4) +
    ggtitle("MA Plot") + xlab("Average Expression (AveExpr)") + ylab("Log-Fold Change (logFC)") +
    scale_color_identity() + theme_minimal()+
    annotate("text", x = max(df$baseMean), y = max(df$log2FoldChange), 
             label = paste("Up:", sum(df$log2FoldChange > logfc & df$padj < fdr)),
             hjust = 1, vjust = 1, color = "black") +
    annotate("text", x = max(df$baseMean), y = min(df$log2FoldChange), 
             label = paste("Down:", sum(df$log2FoldChange < -logfc & df$padj < fdr)),
             hjust = 1, vjust = 1, color = "black")
}

fn_maplot_general_shrinkage <- function(df, columnbyorder, fdr, logfc){
  plotlist <- list()
  message("rearranging and plotting ma...")
  df <- df[,columnbyorder]
  colnames(df) <- c("baseMean","lfcSE","padj","feature")
  df$color <- ifelse(df$lfcSE > logfc & df$padj < fdr, "#B31B21", 
                     ifelse(df$lfcSE < -logfc & df$padj < fdr, "#1465AC", "grey"))
  plotlist[["maplot"]] <- ggplot(df, aes(x = baseMean, y = lfcSE, color = color)) + geom_point(alpha = 0.8, size = 0.4) +
    ggtitle("MA Plot") + xlab("Average Expression (AveExpr)") + ylab("Log-Fold Change Shrinkage (logFCSE)") +
    scale_color_identity() + theme_minimal()+
    annotate("text", x = max(df$baseMean), y = max(df$lfcSE), 
             label = paste("Up:", sum(df$lfcSE > logfc & df$padj < fdr)),
             hjust = 1, vjust = 1, color = "black") +
    annotate("text", x = max(df$baseMean), y = min(df$lfcSE), 
             label = paste("Down:", sum(df$lfcSE < -logfc & df$padj < fdr)),
             hjust = 1, vjust = 1, color = "black")
}


fn_rearrange_groups <- function(grouplist){
  storelist <- list()
  for (i in grouplist){
    genes <- data.frame(fread(i, header = F))
    colnames(genes) <- "gene"
    message("Analysing gene groups...", basename(i))
    gene_group <- merge(genes, integration_features_list$gene, by="gene")
    gene_group["group"] <- basename(i)
    gene_group <- gene_group[which(gene_group$gene != "Metazoa_SRP"),]
    storelist[[basename(i)]] <- gene_group
  }
  storelist_df <- do.call(rbind.data.frame, storelist)
  return(storelist_df)
}


fn_meta_plot <- function(featurematrix, columnstoplot, multi=FALSE, value="value",color="col", linetype="col", rep, splitside=3, colorbyid=FALSE){
  storelist <- list()
  for (i in names(featurematrix)){
    message("rearranging matrix for ",i)
    df <- featurematrix[[i]][["log_normalized_counts"]][,columnstoplot]
    df_st <- data.frame(stack(as.matrix(df)))
    message("plotting meta plots for ",i)
    if (mean(sapply(strsplit(as.character(df_st$row), "%"), function(x) length(x))) <= 2){
      message("no name splitting required for ", i)
      storelist[[i]][["df"]] <- df
      storelist[[i]][["df_st"]] <- df_st
      storelist[[i]][["meta_plot"]] <- ggplot(df_st, aes_string(x=value, color=color, linetype=color)) + geom_density() +
        scale_linetype_manual(values = c(rep("solid",rep),rep("longdash",rep), rep("dashed",rep),rep("dotted",rep))) +
        geom_vline(aes(xintercept=0), color="grey", linetype="dashed", size=0.2) + theme_classic()
    }else if (mean(sapply(strsplit(as.character(df_st$row), "%"), function(x) length(x)))  > 2){
      message("splitting name required for ", i, ", splitting...")
      df_st["id"] <- sapply(strsplit(as.character(df_st$row), "%"), function(x) x[splitside])
      storelist[[i]][["df"]] <- df
      storelist[[i]][["df_st"]] <- df_st
      if (!grepl("bins", i, ignore.case = TRUE)){
        storelist[[i]][["meta_plot"]] <- ggplot(df_st, aes_string(x=value, color="id", linetype=color)) + geom_density() +
          scale_linetype_manual(values = c(rep("solid",rep),rep("longdash",rep),rep("dashed",rep),rep("dotted",rep))) +
          geom_vline(aes(xintercept=0), color="grey", linetype="dashed", size=0.2) + theme_classic()
      }else{
        storelist[[i]][["meta_plot"]] <- ggplot(df_st, aes_string(x=value, color=color, linetype=color)) + geom_density() +
          scale_linetype_manual(values = c(rep("solid",rep),rep("longdash",rep),rep("dashed",rep),rep("dotted",rep))) +
          geom_vline(aes(xintercept=0), color="grey", linetype="dashed", size=0.2) + theme_classic()
      }
    }
  }
  return(storelist)
}

fn_annotate_peaks_to_genes <- function(peakfileslist, genefile, path, assaytype, software, outformat){
  for (i in peakfileslist){
    print(i)
    all_df <- data.frame(fread(paste0(path, "/", i)))
    message("writing all data...")
    write.table(all_df, paste0(path,"/",assaytype, "_",software,"_",i,"_","all.",outformat), sep="\t", quote = F, append = F, row.names = F, col.names = F)
    write.table(genefile, paste0(path,"/","gene_gencode_v41_out.",outformat), sep="\t", quote = F, append = F, row.names = F, col.names = F)
    message("sorting...")
    system(paste0("sort -k1,1 -k2,2n ", paste0(path,"/",assaytype, "_",software,"_",i,"_","all.",outformat),  " | grep chr > ", paste0(path,"/",assaytype, "_",software,"_",i,"_","all_sorted.",outformat)))
    message("annotation to genes ...")
    all_anno <- data.frame(fread(cmd=paste0("bedtools closest -a ", paste0(path,"/",assaytype, "_",software,"_",i,"_","all_sorted.",outformat), " -b ", paste0(path,"/","gene_gencode_v41_out.",outformat), " -d")))
    colnames(all_anno) <- c(colnames(all_df), "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_type", "ensID", "gene", "ens_gene","Distance")
    storelist <- list()
    storelist[[i]][["all"]][["all"]] <- all_anno
    message("getting peaks nearest distance to genes...")
    nearest_distance = 0
    all_anno_nearest <- all_anno %>% dplyr::filter(Distance == nearest_distance)
    storelist[[i]][["nearest"]][["all"]] <- all_anno_nearest
    return(storelist)
  }
}

fn_go_term_analysis <- function(gene_list, organism){
  storelist <- list()
  for (names in names(gene_list)) {
    print(names)
    message("getting gene names ...")
    geneset <- sapply(strsplit(gene_list[[names]], "%"), function(x) x[2])
    message("total genes belonging to ", names, " are...")
    print(length(geneset))
    message("performing gost ...")
    gost <- gost(query=geneset, organism=organism, evcodes = TRUE)
    storelist[[names]][["gost"]] <- gost
    if (length(gost) > 0){
      gost <- gost$result[order(gost$result$p_value),]
      gost["term_name_collapse"] <- paste0(gost$term_id,"_",gost$source,"_" ,gsub(" ", ".", gost$term_name))
      storelist[[names]][["gost_re"]] <- gost
      message("extracting GO:BP ...")
      gost_gobp <- gost[which(gost$source == "GO:BP"),]
      gost_gobp_bar <- gost_gobp[,c(11,3)]
      gost_gobp_bar_top <- head(gost_gobp_bar,10)
      gost_gobp_bar_top$p_value <- -log10(gost_gobp_bar_top$p_value)
      storelist[[names]][["top"]] <- gost_gobp_bar_top
      message("plotting barplot ...")
      plotBar <- ggbarplot(gost_gobp_bar_top, x = "term_name", y = "p_value", color = "#5d93c4ff", fill = "#5d93c4ff" , sort.by.groups = FALSE,x.text.angle = 90, ylab = "-log10(p.value)", xlab = "Biological Process", legend.title = gsub("gost.","",names), lab.size = 9, sort.val = "asc", rotate = TRUE,  position = position_dodge(),ggtheme = theme_bw())
      storelist[[names]][["barplot"]] <- plotBar
    }
  }
  return(storelist)
}

fn_scatterplot <- function(df){
  storelist <- list()
  for (i in colnames(df)){
    for (j in colnames(df)){
      message("sample in analysis: ",i," ",j)
      message("plotting scatterplot...")
      storelist[[paste0(i,"_",j)]][["scatter"]] <- ggplot(df, aes_string(x=i, y=j)) + geom_point(shape=18, color="#2c7bb0", size=1)+
        geom_smooth(method=lm, color="red") + theme_classic()+
        stat_cor(method = "pearson", label.x = -4, label.y = 5, r.digits = 3, aes(label = ..r.label..))
      message("correlation analysis...")
      storelist[[paste0(i,"_",j)]][["cor"]] <- cor(df[[i]], df[[j]], method = 'pearson')
    }
  }
  return(storelist)
}


fn_replot_pie_chipseeker <- function(gr_anno){
  storelist <- list()
  pie_data_gr_anno_all <- c()
  for (i in names(gr_anno)){
    print(i)
    pie_data_gr_anno <- gr_anno[[i]]@annoStat
    pie_data_gr_anno["sample"] <- i
    pie_data_gr_anno_all <- rbind.data.frame(pie_data_gr_anno_all, pie_data_gr_anno)
  }
  storelist[["df_st"]] <- pie_data_gr_anno_all 
  return(storelist)
}


###########################################
####  common_files among all datasets   ###
###########################################
integration_features_list <- list()
# include upstream of a gene https://support.bioconductor.org/p/78652/
integration_features_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_gencode_out.sorted.chr.txt")
integration_features_list[["gene_ext"]] <- makeGRangesFromDataFrame(integration_features_list[["gene"]], keep.extra.columns=TRUE)
# extend the upstream by 2kb (to capture promoter)
integration_features_list$gene_ext <- data.frame(resize(integration_features_list$gene_ext, width = width(integration_features_list$gene_ext) + 2000, fix = "end"))
integration_features_list$gene_ext <- integration_features_list$gene_ext[which(integration_features_list$gene_ext$start > 0),]
write.table(integration_features_list$gene_ext[,c(1,2,3,5:8)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/","gene_gencode_human_upstream_2kb.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)
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
write.table(integration_features_list$promoter[,c(1,2,3,7,8,5)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/","gene_gencodev41_promoter.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

# tss of genes
integration_features_list[["transcript"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "transcript_v41_human_transcript_grch38.sorted.txt")
integration_features_list[["tss"]] <- data.frame(GenomicRanges::resize(makeGRangesFromDataFrame(integration_features_list[["transcript"]], keep.extra.columns=TRUE), width=2, fix='start'))
write.table(integration_features_list$tss[,c(1,2,3,7,9,5)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/","transcript_gencodev41_tss.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

integration_features_list[["gene_tss"]] <- data.frame(GenomicRanges::resize(makeGRangesFromDataFrame(integration_features_list[["gene"]], keep.extra.columns=TRUE), width=2, fix='start'))
write.table(integration_features_list$gene_tss[,c(1,2,3)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/","gene_gencodev41_genetss.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

write.table(integration_features_list$gene_ext[,c(1,2,3,7,8)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/","gene_gencodev41_gene.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

# rnaseqkd_gene_groups
integration_path <- "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration"


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

# with(rnaseqkd_list$deseq2$ma_plots$shC1_c1509$data, plot(mean, lfc, pch=19, cex = 0.3, main="Ma plot", xlim=c(0,15), col="grey", ylim=c(-8,8)))
# #Upregulated
# with(subset(rnaseqkd_list$deseq2$ma_plots$shC1_c1509$data, padj<.05 & lfc > 1), points(mean, lfc, pch=19, cex = 0.4, xlim=c(0,15), ylim=c(-8,8),lwd = 0.4, col="black",bg="grey"))
# #Downregulated
# with(subset(rnaseqkd_list$deseq2$ma_plots$shC1_c1509$data, padj<.05 & lfc < -1), points(mean, lfc, pch=19, cex = 0.4, xlim=c(0,15), ylim=c(-8,8), lwd = 0.4,col="black", bg="grey"))

# Overlap between cramp1 and suz12 target genes
# Venn diagram
list_rnaseqkd_de_shC1_c1509_shS_c1509 <- list(shC1_c1509 = unique(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$de$ensID_Gene), shS_c1509=unique(rnaseqkd_list$deseq2$de_analysis$shS_c1509$de$ensID_Gene))
fn_plot_venn(list_rnaseqkd_de_shC1_c1509_shS_c1509, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_rnaseqkd_up_shC1_c1509_shS_c1509 <- list(shC1_c1509 = unique(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$up$ensID_Gene), shS_c1509=unique(rnaseqkd_list$deseq2$de_analysis$shS_c1509$up$ensID_Gene))
fn_plot_venn(list_rnaseqkd_up_shC1_c1509_shS_c1509, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_rnaseqkd_down_shC1_c1509_shS_c1509 <- list(shC1_c1509 = unique(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$down$ensID_Gene), shS_c1509=unique(rnaseqkd_list$deseq2$de_analysis$shS_c1509$down$ensID_Gene))
fn_plot_venn(list_rnaseqkd_down_shC1_c1509_shS_c1509, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

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
write.table(set_rnaseqde_merged_shC1_shS_up_pos[,c(3:5,1,2,6)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/","set_rnaseqde_merged_shC1_shS_up_pos.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

# export files
write.xlsx(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$de, file = "rnaseqkd_list_deseq2_de_analysis_shC1_c1509_de.xlsx", colNames = TRUE, rowNames = TRUE)
write.xlsx(rnaseqkd_list$deseq2$de_analysis$shS_c1509$de, file = "rnaseqkd_list_deseq2_de_analysis_shS_c1509_de.xlsx", colNames = TRUE, rowNames = TRUE)

# go term analysis
# de
# list_rnaseqkd_de_shC1_c1509_shS_c1509
rnaseqkd_list[["goterm"]][["de"]] <- fn_go_term_analysis(list_rnaseqkd_de_shC1_c1509_shS_c1509, "hsapiens")

# list_rnaseqkd_up_shC1_c1509_shS_c1509
rnaseqkd_list[["goterm"]][["up"]] <- fn_go_term_analysis(list_rnaseqkd_up_shC1_c1509_shS_c1509, "hsapiens")

# list_rnaseqkd_down_shC1_c1509_shS_c1509
rnaseqkd_list[["goterm"]][["down"]] <- fn_go_term_analysis(list_rnaseqkd_down_shC1_c1509_shS_c1509, "hsapiens")

# unique category
# remove none_up
list_set_rnaseqde_merged_shC1_shS_up_filt <- list_set_rnaseqde_merged_shC1_shS_up
list_set_rnaseqde_merged_shC1_shS_up_filt$none_up <- NULL
rnaseqkd_list[["goterm"]][["unique"]] <- fn_go_term_analysis(list_set_rnaseqde_merged_shC1_shS_up_filt, "hsapiens")

# ma plot shrinkage
rnaseqkd_list[["deseq2"]][["maplots_shrinkage"]][["shC1_c1509"]] <- fn_maplot_general_shrinkage(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$all, c(1,3,6,7), 0.05, 1)

################################################
###       ATAC-Seq data: Knockdown  PE       ###
################################################

# deseq2 with consensus peaks (direct from nextflow)
atacseqkd_deseq2_path <- "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/"

atacseqkd_deseq2_list <- list()

atacseqkd_deseq2_list[["dds"]] <- fn_load_rdata(atacseqkd_deseq2_path, "boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/deseq2/consensus_peaks.mLb.clN.dds.RData")
atacseqkd_deseq2_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_upstream_2kb.sorted.txt")
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

atacseqkd_deseq2_list[["chipseeker"]][["annotation"]] <- fn_chipseeker_region_anno("atacseqkd", "deseq2", atacseqkd_deseq2_list$dar_analysis$de_analysis, atacseqkd_deseq2_list$consensus[,c(1:4)], "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2) 

atacseqkd_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl[["maplot"]] <- fn_maplot(atacseqkd_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$all, c(1,2,6,7), 0.05, 2)
atacseqkd_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl[["maplot"]] <- fn_maplot(atacseqkd_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$all, c(1,2,6,7), 0.05, 2)

write.table(atacseqkd_deseq2_list$annotation$shCRAMP1_shControl$all$de[,c(8:10,3,7,1,11:12)] %>% distinct(), "atacseqkd_deseq2_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_deseq2_list$annotation$shSUZ12_shControl$all$de[,c(8:10,3,7,1,11:12)] %>% distinct(), "atacseqkd_deseq2_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_deseq2_de_shCRAMP1_shSUZ12_shControl <- list(shCRAMP1 = unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$de$feature_id), shSUZ12=unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$de$feature_id))
fn_plot_venn(list_atacseqkd_deseq2_de_shCRAMP1_shSUZ12_shControl, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_atacseqkd_deseq2_up_shCRAMP1_shSUZ12_shControl <- list(shCRAMP1 = unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$up$feature_id), shSUZ12=unique(atacseqkd_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$up$feature_id))
fn_plot_venn(list_atacseqkd_deseq2_up_shCRAMP1_shSUZ12_shControl, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotPCA(atacseqkd_deseq2_list$dar_analysis$vstO, intgroup="sample", returnData=FALSE)


atacseqkd_deseq2_list[["aggregate_pergene"]][["shCRAMP1_shControl"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_deseq2_list$annotation$shCRAMP1_shControl$nearest$de, c(3), "ens_gene", mean)
atacseqkd_deseq2_list[["aggregate_pergene"]][["shSUZ12_shControl"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_deseq2_list$annotation$shSUZ12_shControl$nearest$de, c(3), "ens_gene", mean)

colnames(atacseqkd_deseq2_list$aggregate_pergene$shCRAMP1_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(atacseqkd_deseq2_list$aggregate_pergene$shSUZ12_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")

# edgeR with consensus peaks (direct from nextflow)
atacseqkd_edger_list <- list()

atacseqkd_edger_list[["dds"]] <- fn_load_rdata(atacseqkd_deseq2_path, "boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/deseq2/consensus_peaks.mLb.clN.dds.RData")
atacseqkd_edger_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_upstream_2kb.sorted.txt")
atacseqkd_edger_list[["counts"]] <- fn_dds_counts(atacseqkd_edger_list$dds)
atacseqkd_edger_list[["samples_selected"]] <- fn_sample_filter(atacseqkd_edger_list$counts, samplefilter = FALSE)
atacseqkd_edger_list[["processed_counts"]] <- atacseqkd_edger_list[["samples_selected"]][,c(6,4,5,1,3,2)] # Specific requirement here
atacseqkd_edger_list[["coldata"]] <- fn_make_coldata_dds(atacseqkd_edger_list$dds, samplefilter = FALSE)
colnames(atacseqkd_edger_list$coldata)[colnames(atacseqkd_edger_list$coldata) == "condition"] <- "Condition" # Specific requirement here
atacseqkd_edger_list$coldata$Condition <- factor(atacseqkd_edger_list$coldata$Condition)
atacseqkd_edger_list$coldata <- atacseqkd_edger_list$coldata[c(6,4,5,1,3,2),c(1:3)] # Specific requirement here

atacseqkd_edger_list[["dge_obj"]] <- fn_edger_create(atacseqkd_edger_list$processed_counts, atacseqkd_edger_list$coldata)
atacseqkd_edger_list[["dar_analysis"]] <- fn_edger_de(atacseqkd_edger_list$dge_obj$lmfit, atacseqkd_edger_list$dge_obj$design, 0.05, 2)
atacseqkd_edger_list[["consensus"]] <- data.frame(fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.bed"))
atacseqkd_edger_list[["annotation"]] <- fn_consensus_gene_anno("atacseqkd", "edger", atacseqkd_edger_list$dar_analysis$dars, atacseqkd_edger_list$consensus, atacseqkd_edger_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)
# atacseqkd_edger_list[["chipseeker"]][["plots"]] <- fn_run_chipseeker("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/", "atacseqkd", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "*.broadPeak")
atacseqkd_edger_list[["chipseeker"]][["annotation"]] <- fn_chipseeker_region_anno("atacseqkd", "edger", atacseqkd_edger_list$dar_analysis$dars, atacseqkd_edger_list$consensus[,c(1:4)], "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2) 
atacseqkd_edger_list$dar_analysis$dars$coef_shCRAMP1[["maplot"]] <- fn_maplot(atacseqkd_edger_list$dar_analysis$dars$coef_shCRAMP1$all, c(2,1,5,6), 0.05, 2)
atacseqkd_edger_list$dar_analysis$dars$coef_shSUZ12[["maplot"]] <- fn_maplot(atacseqkd_edger_list$dar_analysis$dars$coef_shSUZ12$all, c(2,1,5,6), 0.05, 2)

# get positions of dars
write.table(atacseqkd_edger_list$annotation$coef_shCRAMP1$all$de[,c(7:9,2,6,1,10:11)] %>% distinct(), "atacseqkd_edger_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_edger_list$annotation$coef_shSUZ12$all$de[,c(7:9,2,6,1,10:11)] %>% distinct(), "atacseqkd_edger_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_edger_de_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_edger_list$dar_analysis$dars$coef_shCRAMP1$de$feature_id), shSUZ12=unique(atacseqkd_edger_list$dar_analysis$dars$coef_shSUZ12$de$feature_id))
fn_plot_venn(list_atacseqkd_edger_de_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_atacseqkd_edger_up_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_edger_list$dar_analysis$dars$coef_shCRAMP1$up$feature_id), shSUZ12=unique(atacseqkd_edger_list$dar_analysis$dars$coef_shSUZ12$up$feature_id))
fn_plot_venn(list_atacseqkd_edger_up_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotMDS(atacseqkd_edger_list$dge_obj$mds,top = 100000, pch = 16, cex = 1, dim.plot = c(1,2), col=c(rep("black",2), rep("red",2), rep("blue",2)))
plotMDS(atacseqkd_edger_list$dge_obj$mds)

atacseqkd_edger_list[["aggregate_pergene"]][["coef_shCRAMP1"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_edger_list$annotation$coef_shCRAMP1$nearest$de, c(3), "ens_gene", mean)
atacseqkd_edger_list[["aggregate_pergene"]][["coef_shSUZ12"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_edger_list$annotation$coef_shSUZ12$nearest$de, c(3), "ens_gene", mean)

colnames(atacseqkd_edger_list$aggregate_pergene$coef_shCRAMP1$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(atacseqkd_edger_list$aggregate_pergene$coef_shSUZ12$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")

# edgeR with consensus peaks (direct from nextflow)
atacseqkd_limma_list <- list()

atacseqkd_limma_list[["dds"]] <- fn_load_rdata(atacseqkd_deseq2_path, "boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/deseq2/consensus_peaks.mLb.clN.dds.RData")
atacseqkd_limma_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_upstream_2kb.sorted.txt")
atacseqkd_limma_list[["counts"]] <- fn_dds_counts(atacseqkd_limma_list$dds)
atacseqkd_limma_list[["samples_selected"]] <- fn_sample_filter(atacseqkd_limma_list$counts, samplefilter = FALSE)
atacseqkd_limma_list[["processed_counts"]] <- atacseqkd_limma_list[["samples_selected"]][,c(6,4,5,1,3,2)] # Specific requirement here
atacseqkd_limma_list[["coldata"]] <- fn_make_coldata_dds(atacseqkd_limma_list$dds, samplefilter = FALSE)
colnames(atacseqkd_limma_list$coldata)[colnames(atacseqkd_limma_list$coldata) == "condition"] <- "Condition" # Specific requirement here
atacseqkd_limma_list$coldata$Condition <- factor(atacseqkd_limma_list$coldata$Condition)
atacseqkd_limma_list$coldata <- atacseqkd_limma_list$coldata[c(6,4,5,1,3,2),c(1:3)] # Specific requirement here


atacseqkd_limma_list[["dge_obj"]] <- fn_limma_create(atacseqkd_limma_list$processed_counts, atacseqkd_limma_list$coldata)

# create contrasts
atacseqkd_limma_list[["contrasts"]] <- limma::makeContrasts(paste0(colnames(atacseqkd_limma_list$dge_obj$design)[c(2,1)], collapse = "-"),
                                                            paste0(colnames(atacseqkd_limma_list$dge_obj$design)[c(3,1)], collapse = "-"), levels = atacseqkd_limma_list$dge_obj$design)

atacseqkd_limma_list[["dar_analysis"]] <- fn_limma_de(atacseqkd_limma_list$dge_obj$fit, atacseqkd_limma_list$contrasts, 0.05, 2)

atacseqkd_limma_list[["consensus"]] <- data.frame(fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.bed"))
# scaffolds with non-chr annotation wil be filtered
atacseqkd_limma_list[["annotation"]] <- fn_consensus_gene_anno("atacseqkd", "limma", atacseqkd_limma_list$dar_analysis$dars,atacseqkd_limma_list$consensus, atacseqkd_limma_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)
atacseqkd_limma_list[["chipseeker"]][["annotation"]] <- fn_chipseeker_region_anno("atacseqkd", "limma", atacseqkd_limma_list$dar_analysis$dars, atacseqkd_limma_list$consensus[,c(1:4)], "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2) 

# Create MA plot
atacseqkd_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl[["maplot"]] <- fn_maplot_general(atacseqkd_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl$all, c(2,1,5,7), 0.05, 1)
atacseqkd_limma_list$dar_analysis$dars$coef_shSUZ12_shControl[["maplot"]] <- fn_maplot_general(atacseqkd_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$all, c(2,1,5,7), 0.05, 1)

# get positions of dars
write.table(atacseqkd_limma_list$annotation$coef_shCRAMP1$all$de[,c(8:10,2,6,1,11:12)] %>% distinct(), "atacseqkd_limma_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_limma_list$annotation$coef_shSUZ12$all$de[,c(8:10,2,6,1,11:12)] %>% distinct(), "atacseqkd_limma_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_limma_de_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl$de$feature_id), shSUZ12=unique(atacseqkd_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$de$feature_id))
fn_plot_venn(list_atacseqkd_limma_de_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

list_atacseqkd_limma_up_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl$up$feature_id), shSUZ12=unique(atacseqkd_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$up$feature_id))
fn_plot_venn(list_atacseqkd_limma_up_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotMDS(atacseqkd_limma_list$dge_obj$dgelist,top = 100000, pch = 16, cex = 1, dim.plot = c(1,2), col=c(rep("black",2), rep("red",2), rep("blue",2)))
plotMDS(atacseqkd_limma_list$dge_obj$dgelist)

atacseqkd_limma_list[["aggregate_pergene"]][["coef_shCRAMP1_shControl"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_limma_list$annotation$coef_shCRAMP1_shControl$nearest$de, c(3), "ens_gene", mean)
atacseqkd_limma_list[["aggregate_pergene"]][["coef_shSUZ12_shControl"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_limma_list$annotation$coef_shSUZ12_shControl$nearest$de, c(3), "ens_gene", mean)

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
atacseqkd_diffbind_list[["prep"]] <- fn_diffbind_prep(diffbind_atacseqkd_samplesheet)
plot(atacseqkd_diffbind_list$prep$dba_obj)
plot(atacseqkd_diffbind_list$prep$ovp_rate,type ='b',ylab='# peaks', xlab='Overlap at least this many peaksets')

atacseqkd_diffbind_list[["counts"]] <- fn_diffbind_count(atacseqkd_diffbind_list$prep$dba_obj, 100)
atacseqkd_diffbind_list[["regions"]] <- data.frame(atacseqkd_diffbind_list$counts$counts[,c(1:3)], interval=unique(do.call(paste, c(atacseqkd_diffbind_list$counts$counts[,1:3], sep="%"))))
atacseqkd_diffbind_list[["norm"]] <-  fn_diffbind_norm(atacseqkd_diffbind_list$counts$dba_obj)

atacseqkd_diffbind_list[["dar_analysis"]][["edgeR"]] <- fn_diffbind_de(atacseqkd_diffbind_list$norm$norm, 0.05, 2, "edgeR")
atacseqkd_diffbind_list[["dar_analysis"]][["deseq2"]] <- fn_diffbind_de(atacseqkd_diffbind_list$norm$norm, 0.05, 2, "deseq2")

atacseqkd_diffbind_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_upstream_2kb.sorted.txt")

atacseqkd_diffbind_list[["annotation"]][["deseq2"]] <- fn_diffbind_gene_anno("atacseqkd", "diffbind",atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts, atacseqkd_diffbind_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)
atacseqkd_diffbind_list[["annotation"]][["edgeR"]] <- fn_diffbind_gene_anno("atacseqkd", "diffbind",atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts, atacseqkd_diffbind_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

atacseqkd_diffbind_list[["chipseeker"]][["annotation"]][["deseq2"]] <- fn_chipseeker_region_anno("atacseqkd", "diffbind",atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts, atacseqkd_diffbind_list$regions, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2) 
atacseqkd_diffbind_list[["chipseeker"]][["annotation"]][["edgeR"]] <- fn_chipseeker_region_anno("atacseqkd", "diffbind",atacseqkd_diffbind_list$dar_analysis$edger$contrasts, atacseqkd_diffbind_list$regions, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2) 

# since we need only deseq2 results for downstream analysis, I perform mean analysis with deseq2
atacseqkd_diffbind_list[["aggregate_pergene"]][["deseq2"]][["contrast_shCRAMP1_shControl"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_diffbind_list$annotation$deseq2$contrast_shCRAMP1_shControl$nearest$de, c(9), "ens_gene", mean)
atacseqkd_diffbind_list[["aggregate_pergene"]][["deseq2"]][["contrast_shSUZ12_shControl"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_diffbind_list$annotation$deseq2$contrast_shSUZ12_shControl$nearest$de, c(9), "ens_gene", mean)

colnames(atacseqkd_diffbind_list$aggregate_pergene$deseq2$contrast_shCRAMP1_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(atacseqkd_diffbind_list$aggregate_pergene$deseq2$contrast_shSUZ12_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")

# Create MA plot
atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shCRAMP1_shControl[["maplot"]] <- fn_maplot(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shCRAMP1_shControl$all, c(6,9,11,12), 0.05, 2)
atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shSUZ12_shControl[["maplot"]] <- fn_maplot(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shSUZ12_shControl$all, c(6,9,11,12), 0.05, 2)
atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shControl[["maplot"]] <- fn_maplot(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shControl$all, c(6,9,11,12), 0.05, 2)
atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl[["maplot"]] <- fn_maplot(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl$all, c(6,9,11,12), 0.05, 2)

# e= edger and d= deseq2
write.table(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shCRAMP1_shControl$de[,c(1:3,9,11:12)] %>% distinct(), "atacseqkd_diffbind_e_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shSUZ12_shControl$de[,c(1:3,9,11:12)] %>% distinct(), "atacseqkd_diffbind_e_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shControl$de[,c(1:3,9,11:12)] %>% distinct(), "atacseqkd_diffbind_d_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl$de[,c(1:3,9,11:12)] %>% distinct(), "atacseqkd_diffbind_d_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# extract dars in specified format
write.table(cbind(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shControl$de[,c(1:3,12)] %>% distinct(), dars="dars", strand="."), "atacseqkd_diffbind_d_shCRAMP1_shControl_de_re.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(cbind(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl$de[,c(1:3,12)] %>% distinct(), dars="dars", strand="."), "atacseqkd_diffbind_d_shSUZ12_shControl_de_re.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_diffbind_e_de_shCRAMP1_vs_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shCRAMP1_shContro$de$feature_id), shSUZ12=unique(atacseqkd_diffbind_list$dar_analysis$edgeR$contrasts$contrast_shSUZ12_shControl$de$feature_id))
fn_plot_venn(list_atacseqkd_diffbind_e_de_shCRAMP1_vs_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))
list_atacseqkd_diffbind_d_de_shCRAMP1_vs_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shCRAMP1_shContro$de$feature_id), shSUZ12=unique(atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts$contrast_shSUZ12_shControl$de$feature_id))
fn_plot_venn(list_atacseqkd_diffbind_d_de_shCRAMP1_vs_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

#PCA
dba.plotPCA(atacseqkd_diffbind_list$norm$norm, method=DBA_ALL_METHODS, attributes=DBA_FACTOR, label=DBA_ID)

# check for a match before and after dars analysis
atacseqkd_diffbind_list[["counts"]][["reads"]]  <- dba.count(atacseqkd_diffbind_list$counts$dba_obj, peaks=NULL, score=DBA_SCORE_READS)
atacseqkd_diffbind_list[["counts"]][["counts"]] <- dba.peakset(atacseqkd_diffbind_list$counts$reads, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
atacseqkd_diffbind_counts <- data.frame(atacseqkd_diffbind_list$counts$counts)
rownames(atacseqkd_diffbind_counts) <- paste0(atacseqkd_diffbind_counts$CHR, "%",atacseqkd_diffbind_counts$START,"%", atacseqkd_diffbind_counts$END)
atacseqkd_diffbind_counts <- atacseqkd_diffbind_counts[,-c(1:3)]

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
atacseqkd_diffbind_deseq2_list[["dar_analysis"]] <- fn_deseq2(atacseqkd_diffbind_deseq2_list$counts, atacseqkd_diffbind_deseq2_list$coldata, 0.05, 2)

atacseqkd_diffbind_deseq2_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_upstream_2kb.sorted.txt")
# get positions of dars
atacseqkd_diffbind_deseq2_list[["annotation"]] <- fn_diffbind_and_packages_gene_anno("atacseqkd", "diffbind_deseq2", atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis, atacseqkd_diffbind_deseq2_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl[["maplot"]] <- fn_maplot(atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$all, c(1,2,6,7), 0.05, 2)
atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl[["maplot"]] <- fn_maplot(atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$all, c(1,2,6,7), 0.05, 2)

write.table(atacseqkd_diffbind_deseq2_list$annotation$shCRAMP1_shControl$all$de[,c(1:3,5,9)] %>% distinct(), "atacseqkd_diffbind_deseq2_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_deseq2_list$annotation$shSUZ12_shControl$all$de[,c(1:3,5,9)] %>% distinct(), "atacseqkd_diffbind_deseq2_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_diffbind_deseq2_de_shCRAMP1_shSUZ12_shControl <- list(shCRAMP1 = unique(atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shCRAMP1_shControl$de$feature_id), shSUZ12=unique(atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis$shSUZ12_shControl$de$feature_id))
fn_plot_venn(list_atacseqkd_diffbind_deseq2_de_shCRAMP1_shSUZ12_shControl, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotPCA(atacseqkd_diffbind_deseq2_list$dar_analysis$vstO, intgroup="ID", returnData=FALSE)


# DE deseq2 manually from diffbind mean the logFC per peaks to specify genes
atacseqkd_diffbind_deseq2_list[["aggregate_pergene"]][["deseq2"]][["shCRAMP1_shControl"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_diffbind_deseq2_list$annotation$shCRAMP1_shControl$nearest$de, c(5), "ens_gene", mean)
colnames(atacseqkd_diffbind_deseq2_list$aggregate_pergene$deseq2$shCRAMP1_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")

atacseqkd_diffbind_deseq2_list[["aggregate_pergene"]][["deseq2"]][["shSUZ12_shControl"]][["nearest"]][["de"]] <- fn_aggregate_feature(atacseqkd_diffbind_deseq2_list$annotation$shSUZ12_shControl$nearest$de, c(5), "ens_gene", mean)
colnames(atacseqkd_diffbind_deseq2_list$aggregate_pergene$deseq2$shSUZ12_shControl$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")


# edgeR analysis
atacseqkd_diffbind_edger_list <- list()

atacseqkd_diffbind_edger_list[["counts"]] <- atacseqkd_diffbind_counts
atacseqkd_diffbind_edger_list[["coldata"]] <- atacseqkd_diffbind_list$counts$info
atacseqkd_diffbind_edger_list$coldata$Condition <- factor(atacseqkd_diffbind_edger_list$coldata$Condition)
atacseqkd_diffbind_edger_list[["dge_obj"]] <- fn_edger_create(atacseqkd_diffbind_edger_list$counts, atacseqkd_diffbind_edger_list$coldata)

atacseqkd_diffbind_edger_list[["dar_analysis"]] <- fn_edger_de(atacseqkd_diffbind_edger_list$dge_obj$lmfit, atacseqkd_diffbind_edger_list$dge_obj$design, 0.05, 2)

atacseqkd_diffbind_edger_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_upstream_2kb.sorted.txt")

atacseqkd_diffbind_edger_list[["annotation"]] <- fn_diffbind_and_packages_gene_anno("atacseqkd", "diffbind_edger", atacseqkd_diffbind_edger_list$dar_analysis$dars, atacseqkd_diffbind_edger_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shCRAMP1[["maplot"]] <- fn_maplot(atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shCRAMP1$all, c(2,1,5,6), 0.05, 2)
atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shSUZ12[["maplot"]] <- fn_maplot(atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shSUZ12$all, c(2,1,5,6), 0.05, 2)

# get positions of dars
write.table(atacseqkd_diffbind_edger_list$annotation$coef_shCRAMP1$all$de[,c(1:4,8)] %>% distinct(), "atacseqkd_diffbind_edger_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_edger_list$annotation$coef_shSUZ12$all$de[,c(1:4,8)] %>% distinct(), "atacseqkd_diffbind_edger_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_diffbind_edger_de_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shCRAMP1$de$feature_id), shSUZ12=unique(atacseqkd_diffbind_edger_list$dar_analysis$dars$coef_shSUZ12$de$feature_id))
fn_plot_venn(list_atacseqkd_diffbind_edger_de_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

plotMDS(atacseqkd_diffbind_edger_list$dge_obj$mds,top = 100000, pch = 16, cex = 1, dim.plot = c(1,2), col=c(rep("black",2), rep("red",2), rep("blue",2)))
plotMDS(atacseqkd_diffbind_edger_list$dge_obj$mds)

# limma analysis
atacseqkd_diffbind_limma_list <- list()

atacseqkd_diffbind_limma_list[["counts"]] <- atacseqkd_diffbind_counts
atacseqkd_diffbind_limma_list[["coldata"]] <- atacseqkd_diffbind_list$counts$info
atacseqkd_diffbind_limma_list$coldata$Condition <- factor(atacseqkd_diffbind_limma_list$coldata$Condition)
atacseqkd_diffbind_limma_list[["dge_obj"]] <- fn_limma_create(atacseqkd_diffbind_limma_list$counts, atacseqkd_diffbind_limma_list$coldata)

# create contrasts
atacseqkd_diffbind_limma_list[["contrasts"]] <- limma::makeContrasts(paste0(colnames(atacseqkd_diffbind_limma_list$dge_obj$design)[c(2,1)], collapse = "-"),
                                                                     paste0(colnames(atacseqkd_diffbind_limma_list$dge_obj$design)[c(3,1)], collapse = "-"), levels = atacseqkd_diffbind_limma_list$dge_obj$design)

atacseqkd_diffbind_limma_list[["dar_analysis"]] <- fn_limma_de(atacseqkd_diffbind_limma_list$dge_obj$fit, atacseqkd_diffbind_limma_list$contrasts, 0.05, 2)

atacseqkd_diffbind_limma_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_upstream_2kb.sorted.txt")

atacseqkd_diffbind_limma_list[["annotation"]] <- fn_diffbind_and_packages_gene_anno("atacseqkd", "diffbind_limma", atacseqkd_diffbind_limma_list$dar_analysis$dars, atacseqkd_diffbind_limma_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)

# Create MA plot
atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl[["maplot"]] <- fn_maplot(atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl$all, c(2,1,5,7), 0.05, 2)
atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shSUZ12_shControl[["maplot"]] <- fn_maplot(atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$all, c(2,1,5,7), 0.05, 2)

# get positions of dars
write.table(atacseqkd_diffbind_limma_list$annotation$coef_shCRAMP1$all$de[,c(1:4,8)] %>% distinct(), "atacseqkd_diffbind_limma_shCRAMP1_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseqkd_diffbind_limma_list$annotation$coef_shSUZ12$all$de[,c(1:4,8)] %>% distinct(), "atacseqkd_diffbind_limma_shSUZ12_shControl_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# Overlap between cramp1 and suz12 target intervals
# Venn diagram
list_atacseqkd_diffbind_limma_de_shCRAMP1_coef_shSUZ12 <- list(shCRAMP1 = unique(atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shCRAMP1_shControl$de$feature_id), shSUZ12=unique(atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$de$feature_id))
fn_plot_venn(list_atacseqkd_diffbind_limma_de_shCRAMP1_coef_shSUZ12, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))

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
atacseqkd_quantify_list[["featurecounts"]][["featuresmatrix"]] <- fn_quantify_featurecounts_multifeature("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "atacseqkd", ".mLb.clN.sorted.bam", atacseqkd_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "shControl",pe=TRUE, diff="multi", atacseqkd_quantify_list$featurecounts$pairs)
atacseqkd_quantify_list[["featurecounts"]][["featuresmatrix_plots"]] <- fn_meta_plot(atacseqkd_quantify_list$featurecounts$featuresmatrix, c(7:10), rep=1, splitside=3)

colnames(atacseqkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$counts) <- gsub(".mLb.clN.sorted.bam", "", colnames(atacseqkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$counts))

plotPCA(edgeR::cpm(atacseqkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$counts), label=TRUE)

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
atacseqkd_quantify_list[["featurecounts"]][["histonemarks"]] <- fn_quantify_featurecounts("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/", "atacseqkd","\\.bam$", "\\_peaks_id.bed$","histone_marks", c(12,4,1:3,7), pairedend=TRUE, "hg38", "_",merge_sites_files=FALSE)

# quantify total reads in a bam file, used Rsamtools
atacseqkd_quantify_list[["bams"]] <- fn_quantify_bams("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/", "atacseqkd","\\.bam$")
atacseqkd_labelpath = "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/histone_marks_label.txt"
atacseqkd_quantify_list[["featurecounts"]][["histonemarks"]][["summed"]] <- fn_summedreads_per_feature(atacseqkd_quantify_list$featurecounts$histonemarks$histone_marks_countmatrix, atacseqkd_quantify_list$bams$total_reads, atacseqkd_labelpath)

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
atacseqkd_quantify_list$featurecounts$histonemarks$summed[["aggregated"]] <- fn_aggregate_feature(atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr_st, c(3), "id", mean)
atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated[["histonemarks"]] <- sapply(strsplit(as.character(atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated$Group.1), "%"),function(x) x[1])
atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated[["sample"]] <- sapply(strsplit(as.character(atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated$Group.1), "%"),function(x) x[2])
atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated[,-1]
atacseqkd_quantify_list$featurecounts$histonemarks$summed[["histone_mark_df"]] <- data.frame(tidyr::pivot_wider(atacseqkd_quantify_list$featurecounts$histonemarks$summed$aggregated$aggregated, names_from = sample, values_from = x))
rownames(atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df) <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df$histonemarks
atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df <- atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df[,-1]
atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df["shCRAMP1"] <- rowMeans(atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df[,1:2])
atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df["shSUZ12"] <- rowMeans(atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df[,3:4])
atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df["shCRAMP1_median"] <- rowMedians(atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df[,1:2], cols=T)
atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df["shSUZ12_median"] <- rowMedians(atacseqkd_quantify_list$featurecounts$histonemarks$summed$histone_mark_df[,3:4])
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

################################################
###         CUT&RUN data: Knockout SE        ###
################################################

# deseq2 with consensus peaks (direct from nextflow)
cutnrunko_deseq2_path <- "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/"

cutnrunko_deseq2_list <- list()

cutnrunko_deseq2_list[["dds"]] <- fn_load_rdata(cutnrunko_deseq2_path, "outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/deseq2/H3K27me3.consensus_peaks.dds.RData")
cutnrunko_deseq2_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_upstream_2kb.sorted.txt")
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
cutnrunko_deseq2_list[["aggregate_pergene"]][["KO_WT"]][["nearest"]][["up"]] <- fn_aggregate_feature(cutnrunko_deseq2_list$annotation$KO_WT$nearest$up, c(3), "ens_gene", mean)
cutnrunko_deseq2_list[["aggregate_pergene"]][["KO_WT"]][["nearest"]][["down"]] <- fn_aggregate_feature(cutnrunko_deseq2_list$annotation$KO_WT$nearest$down, c(3), "ens_gene", mean)

colnames(cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest$de$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest$up$aggregated) <- c("ensID_Gene", "mean_log2FC")
colnames(cutnrunko_deseq2_list$aggregate_pergene$KO_WT$nearest$down$aggregated) <- c("ensID_Gene", "mean_log2FC")

# quantify over various features
# now since I have information for groups to take
selected_rkd_gene_groups <- list.files(integration_path, pattern="cluster*", full.names = T)
integration_features_list[["selected_rkd_gene_groups"]] <- fn_rearrange_groups(selected_rkd_gene_groups)
write.table(integration_features_list$selected_rkd_gene_groups[,c(2:4,8,9,5)], paste0("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/","selected_rkd_gene_groups.txt"), sep="\t", quote=F, append = F, row.names = F, col.names = F)

cutnrunko_quantify_list <- list()
cutnrunko_quantify_list[["featurecounts"]][["featureslist"]]=list(gene_ext=c(".chr_1.5kb.txt", "gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), bins5kb=c("hg38_5kb.txt", "hg38_5kb.txt"), bins10kb=c("hg38_10kb.txt", "hg38_10kb.txt"), promoter=c("_gencodev41_promoter.txt", "gene_gencodev41_promoter.txt"))
cutnrunko_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_groups"]] <-  c("_rkd_gene_groups.txt", "selected_rkd_gene_groups.txt")
cutnrunko_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_uniq_diff"]] <-  c("_shC1_shS_up_pos.txt", "set_rnaseqde_merged_shC1_shS_up_pos.txt")
cutnrunko_quantify_list[["featurecounts"]][["pairs"]] <- list(H3K27me3_KO1_IgG=c("H3K27me3_KO1.mLb.clN.sorted.bam", "IgG_KO1.mLb.clN.sorted.bam"), H3K27me3_KO2_IgG=c("H3K27me3_KO2.mLb.clN.sorted.bam", "IgG_KO2.mLb.clN.sorted.bam"), H3K27me3_KO3_IgG=c("H3K27me3_KO3.mLb.clN.sorted.bam", "IgG_KO3.mLb.clN.sorted.bam"), H3K27me3_WT_IgG=c("H3K27me3_WT.mLb.clN.sorted.bam", "IgG_WT.mLb.clN.sorted.bam"))
cutnrunko_quantify_list[["featurecounts"]][["featuresmatrix"]] <- fn_quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "cutnrunko", ".mLb.clN.sorted.bam", cutnrunko_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=FALSE, diff="multi", cutnrunko_quantify_list$featurecounts$pairs)
cutnrunko_quantify_list[["featurecounts"]][["featuresmatrix_plots"]] <- fn_meta_plot(cutnrunko_quantify_list$featurecounts$featuresmatrix, c(9:12), rep=1, splitside=3)

# quantify overs atacseq dars
cutnrunko_quantify_list[["featurecounts"]][["atacseqdars"]] <- list(shCRAMP1_shControl=c("_diffbind_d_shCRAMP1_shControl_de_re.bed", "atacseqkd_diffbind_d_shCRAMP1_shControl_de_re.bed"), shSUZ12_shControl=c("_diffbind_d_shSUZ12_shControl_de_re.bed", "atacseqkd_diffbind_d_shSUZ12_shControl_de_re.bed"))
cutnrunko_quantify_list[["featurecounts"]][["pairs"]] <- list(H3K27me3_KO1_IgG=c("H3K27me3_KO1.mLb.clN.sorted.bam", "IgG_KO1.mLb.clN.sorted.bam"), H3K27me3_KO2_IgG=c("H3K27me3_KO2.mLb.clN.sorted.bam", "IgG_KO2.mLb.clN.sorted.bam"), H3K27me3_KO3_IgG=c("H3K27me3_KO3.mLb.clN.sorted.bam", "IgG_KO3.mLb.clN.sorted.bam"), H3K27me3_WT_IgG=c("H3K27me3_WT.mLb.clN.sorted.bam", "IgG_WT.mLb.clN.sorted.bam"))
cutnrunko_quantify_list[["featurecounts"]][["darsmatrix"]] <- fn_quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "cutnrunko", ".mLb.clN.sorted.bam", cutnrunko_quantify_list$featurecounts$atacseqdars, c(5,4,1:3,6), "hg38", "IgG",pe=FALSE, diff="multi", cutnrunko_quantify_list$featurecounts$pairs)

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

cutnrunkd_deseq2_list[["dds"]] <- fn_load_rdata(cutnrunkd_deseq2_path, "outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/deseq2/H3K27me3.consensus_peaks.dds.RData")
cutnrunkd_deseq2_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_upstream_2kb.sorted.txt")
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

cutnrunkd_quantify_list <- list()
cutnrunkd_quantify_list[["featurecounts"]][["featureslist"]]=list(gene_ext=c(".chr_1.5kb.txt", "gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), bins5kb=c("hg38_5kb.txt", "hg38_5kb.txt"), bins10kb=c("hg38_10kb.txt", "hg38_10kb.txt"), promoter=c("_gencodev41_promoter.txt", "gene_gencodev41_promoter.txt"))
cutnrunkd_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_groups"]] <-  c("_rkd_gene_groups.txt", "selected_rkd_gene_groups.txt")
cutnrunkd_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_uniq_diff"]] <-  c("_shC1_shS_up_pos.txt", "set_rnaseqde_merged_shC1_shS_up_pos.txt")
cutnrunkd_quantify_list[["featurecounts"]][["pairs"]] <- list(
  shControl_IgG = c("shControl_H3K27me3.mLb.clN.sorted.bam", "shControl_IgG.mLb.clN.sorted.bam"),
  shCRAMP1_IgG = c("shCRAMP1_H3K27me3.mLb.clN.sorted.bam","shCRAMP1_IgG.mLb.clN.sorted.bam"),
  shSUZ12_IgG = c("shSUZ12_H3K27me3.mLb.clN.sorted.bam", "shSUZ12_IgG.mLb.clN.sorted.bam"))
cutnrunkd_quantify_list[["featurecounts"]][["featuresmatrix"]] <- fn_quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "cutnrunkd", ".mLb.clN.sorted.bam", cutnrunkd_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "shControl",pe=FALSE, diff="multi", cutnrunkd_quantify_list$featurecounts$pairs)
cutnrunkd_quantify_list[["featurecounts"]][["featuresmatrix_plots"]] <- fn_meta_plot(cutnrunkd_quantify_list$featurecounts$featuresmatrix, c(7:9), rep=1, splitside=3)

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
cutntagwt_quantify_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_upstream_2kb.sorted.txt")

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

# quantify over various features
cutntagwt_quantify_list[["featurecounts"]][["featureslist"]]=list(gene_ext=c(".chr_1.5kb.txt", "gene_gencode_human_gencode_out.sorted.chr_1.5kb.txt"), bins5kb=c("hg38_5kb.txt", "hg38_5kb.txt"), bins10kb=c("hg38_10kb.txt", "hg38_10kb.txt"), promoter=c("_gencodev41_promoter.txt", "gene_gencodev41_promoter.txt"))
cutntagwt_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_groups"]] <-  c("_rkd_gene_groups.txt", "selected_rkd_gene_groups.txt")
cutntagwt_quantify_list[["featurecounts"]][["featureslist"]][["selected_rkd_uniq_diff"]] <-  c("_shC1_shS_up_pos.txt", "set_rnaseqde_merged_shC1_shS_up_pos.txt")
cutntagwt_quantify_list[["featurecounts"]][["featuresmatrix"]] <- fn_quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "cutntagwt", ".mLb.clN.sorted.bam", cutntagwt_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=FALSE, diff="single")
cutntagwt_quantify_list[["featurecounts"]][["featuresmatrix_plots"]] <- fn_meta_plot(cutntagwt_quantify_list$featurecounts$featuresmatrix, c(5:7), rep=1, splitside=3)


# for check if two df are exactly same and to check if fn_quantify_featurecounts_multifeature is created properly and doing the same thing as individual steps look at Check Notes

# quantify overs atacseq dars
cutntagwt_quantify_list[["featurecounts"]][["atacseqdars"]] <- list(shCRAMP1_shControl=c("_diffbind_d_shCRAMP1_shControl_de_re.bed", "atacseqkd_diffbind_d_shCRAMP1_shControl_de_re.bed"), shSUZ12_shControl=c("_diffbind_d_shSUZ12_shControl_de_re.bed", "atacseqkd_diffbind_d_shSUZ12_shControl_de_re.bed"))
cutntagwt_quantify_list[["featurecounts"]][["darsmatrix"]] <- fn_quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "cutntagwt", ".mLb.clN.sorted.bam", cutntagwt_quantify_list$featurecounts$atacseqdars, c(5,4,1:3,6), "hg38", "IgG", pe=FALSE, diff="single")


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
cutntagwkd_quantify_list[["featurecounts"]] <- fn_quantify_featurecounts("/mnt/home3/reid/av638/cutntag/iva_lab_dec23/outfolder/02_alignment/bowtie2/target/markdup", "cutntagwkd","\\.bam$", "\\_peaks_id.bed$","histone_marks", c(12,4,1:3,7), pairedend=TRUE, "hg38", "_",merge_sites_files=FALSE)

cutntagwkd_quantify_list[["bams"]] <- fn_quantify_bams("/mnt/home3/reid/av638/cutntag/iva_lab_dec23/outfolder/02_alignment/bowtie2/target/markdup", "cutntagwkd","\\.bam$")
cutntagwkd_labelpath = "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/histone_marks_label.txt"
cutntagwkd_quantify_list[["summed"]] <- fn_summedreads_per_feature(cutntagwkd_quantify_list$featurecounts$histone_marks_countmatrix, cutntagwkd_quantify_list$bams$total_reads, cutntagwkd_labelpath)

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
  fn_quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutntag/iva_lab_dec23/outfolder/02_alignment/bowtie2/target/markdup/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "cutntagwkd", "_R1.target.markdup.sorted.bam", 
                                         cutntagwkd_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=TRUE, diff="multi", 
                                         cutntagwkd_quantify_list$featurecounts$pairs)

cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt[["cpm"]] <- edgeR::cpm(cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$counts)
cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb[["pca_counts"]] <- fn_make_pca(cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$cpm, "_R1.target.markdup.sorted.bam", "pca")

cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb[["pca_diff"]] <- fn_make_pca(cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts[,c(23:38)], "_R1.target.markdup.sorted.bam", "pca")

cutntagwkd_quantify_list[["featurecounts"]][["featuresmatrix_plots"]] <- fn_meta_plot(cutntagwkd_quantify_list$featurecounts$featuresmatrix, c(23:38), rep=1, splitside=3)
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
cutntagrmwt_quantify_list[["featurecounts"]] <- fn_quantify_featurecounts("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams", "cutntagrmwt","\\.bam$", "\\_peaks_id.bed$","histone_marks", c(12,4,1:3,7), pairedend=TRUE, "hg38", "_",merge_sites_files=FALSE)

cutntagrmwt_quantify_list[["bams"]] <- fn_quantify_bams("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams", "cutntagrmwt","\\.bam$")
cutntagrmwt_labelpath = "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/histone_marks_label.txt"
cutntagrmwt_quantify_list[["summed"]] <- fn_summedreads_per_feature(cutntagrmwt_quantify_list$featurecounts$histone_marks_countmatrix, cutntagrmwt_quantify_list$bams$total_reads, cutntagrmwt_labelpath)

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
cutntagrmwt_quantify_list$summed[["aggregated"]] <- fn_aggregate_feature(cutntagrmwt_quantify_list$summed$merged_diff_table_afr_st, c(3), "id", mean)
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
  fn_quantify_featurecounts_multifeature("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "cutntagrmwt", "_sorted.bam", 
                                         cutntagrmwt_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=TRUE, diff="multi", 
                                         cutntagrmwt_quantify_list$featurecounts$pairs)

# plot and add the regression line
ggplot(cutntagrmwt_quantify_list$featurecounts$featuresmatrix$bins80kb$log_normalized_counts, aes(x=H14_S17_IgG_S22, y=FH14_S3_F_S6)) + 
  geom_point(shape=18, color="#2c7bb0", size=1)+
  geom_smooth(method=lm, color="darkred") + theme_classic()

cutntagrmwt_quantify_list[["featurecounts"]][["featuresmatrixdiv"]] <- 
  fn_quantify_featurecounts_multifeature_div("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "cutntagrmwt", "_sorted.bam", 
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
cutntagrmwt_quantify_list[["deeptools"]][["scatterplot"]] <- fn_scatterplot(cutntagrmwt_correlation_scatter_divlog2)

# Annotation peaks
cutntagrmwt_quantify_list[["chipseeker"]][["plots"]] <- fn_run_chipseeker("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams/endoH1peaks/", "cutntagrmwt", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "*.broadPeak")

cutntagrmwt_quantify_list[["chipseeker"]][["plots"]][["pie_annodf"]] <- fn_replot_pie_chipseeker(cutntagrmwt_quantify_list$chipseeker$plots$gr_anno)
cutntagrmwt_quantify_list$chipseeker$plots$pie_annodf$df_st$sample <- gsub("endo|_peaks.broadPeak","",cutntagrmwt_quantify_list$chipseeker$plots$pie_annodf$df_st$sample)

ggplot(cutntagrmwt_quantify_list$chipseeker$plots$pie_annodf$df_st, aes(x="Feature", y=Frequency, fill=Feature))+ facet_wrap( ~ sample, ncol=2, nrow=3) +
  geom_bar(width = 1, stat = "identity", color="white") + coord_polar("y", start=0) + theme_void() + scale_fill_manual(values=c("#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251", "#B565A7", "#009B77", "#DD4124", "#D65076")) 

#################################
######## Public data ############
#################################
public_list <- list()
public_list[["peakfilelist"]] <- list.files("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration", "AHF_peaks_id.bed")
public_list[["gene"]] <- fn_read_genefile("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "gene_gencode_human_upstream_2kb.sorted.txt")
public_list[["peaks_anno"]] <- fn_annotate_peaks_to_genes(public_list$peakfilelist, public_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration","chipseq", "encode", "txt")
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
       cutntagwt_quantify_list$featurecounts$gene_ext$log_normalized_counts[,c(7:10)],
       cutntagwkd_quantify_list$featurecounts$featuresmatrix$gene_ext$log_normalized_counts %>% dplyr::select(23:38, ensID_Gene=39))

data_integration_list[["ensID_gene"]][["merged"]] <- fn_integrate_featurecounts(data_integration_list$ensID_gene$list, "ensID_Gene")

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

data_integration_list$ensID_gene$merged[["pca"]][["somepredictable"]]  <- fn_make_pca(data_integration_list$ensID_gene$merged$df, "none", "pca")

# for bins5kb from featurecounts
data_integration_list[["bins5kb"]][["list"]] <- 
  list(atacseqkd_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts %>% dplyr::select(7:10, bins5kb=11),
       cutnrunko_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts %>% dplyr::select(9:12, bins5kb=13),
       cutnrunkd_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts %>% dplyr::select(7:9, bins5kb=10),
       cutntagwt_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts %>% dplyr::select(5:7, bins5kb=8),
       cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb$log_normalized_counts %>% dplyr::select(23:38, bins5kb=39))

data_integration_list[["bins5kb"]][["merged"]] <- fn_integrate_featurecounts(data_integration_list$bins5kb$list, "bins5kb")

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

intergration_omics_list[["bin5kb"]] <- fn_integrate_bedtools(bin_path_combined, "5kb", c(1:3), seq(4,175,7))
intergration_omics_list[["bin10kb"]] <- fn_integrate_bedtools(bin_path_combined, "10kb", c(1:3), seq(4,175,7))


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