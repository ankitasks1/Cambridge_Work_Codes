setwd("/mnt/home3/reid/av638/atacseq/surani_lab/wg")
library(DiffBind)
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
# Prepare diffbind object
diffbind_prep <- function(samplesheet, minoverlap){
  diffblist <- list()
  # Creating diffbind object
  print(paste0('minOverlap:', minoverlap))
  obj <- dba(sampleSheet=samplesheet,  scoreCol=5, minOverlap=minoverlap)
  diffblist[["dba_obj"]] <- obj
  # Plot of number of peaks that appear in atleast one, two , three so on ... until total samples together.
  ovp_rate <- dba.overlap(obj,mode=DBA_OLAP_RATE)
  diffblist[["ovp_rate"]] <- ovp_rate
  return(diffblist)
}

# Create a count matrix using peak coordinates (set range from summit to 100)
# If peaks parameter is missing or a mask, generates a consensus peakset using minOverlap parameter (after applying the mask if present)
diffbind_count <- function(obj, minoverlap, summit){
  diffblist <- list()
  print(paste0('minOverlap:', minoverlap))
  print(paste0('summits: ', summit))
  # Calculate a binding matrix with scores based on read counts,  summits=100 for ATAC-Seq # https://rdrr.io/bioc/DiffBind/man/DiffBind3.html, mapQCth=15 default
  diffblist[["dba_obj"]] <- dba.count(obj, minOverlap=minoverlap, summits=summit, score = DBA_SCORE_READS) # DBA_SCORE_READS = raw read count for interval using only reads from ChIP
  diffblist[["info"]]  <- dba.show(diffblist$dba_obj)
  return(diffblist)
}


diffbind_norm <- function(obj, mimember){
  diffblist <- list()
  print("Performing normalisation...")
  diffblist[["norm"]] <- dba.normalize(obj, method=DBA_DESEQ2, library = DBA_LIBSIZE_FULL) # Normalize only with deseq2 method, # DBA_LIBSIZE_FULL default
  print("Performing Diffbind analysis...")
  print("Performing Contrasts...") # categories based
  print(paste0('mimember: ', mimember))
  diffblist$norm <- dba.contrast(diffblist$norm, minMembers = mimember, categories=DBA_CONDITION)
  diffblist$norm <- dba.analyze(diffblist$norm)
  return(diffblist)
}


diffbind_de <- function(normobj, fdr, fc, method){
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
    }else if(method=="deseq2"){
      storelist <- list()
      storelist[["all"]] <- data.frame(dba.report(normobj, method=DBA_DESEQ2, contrast = i, th=1))
      storelist$all["feature_id"] <- unique(do.call(paste, c(storelist$all[,c(1:3)], sep="%")))
      storelist[["de"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2) | Fold < -log(fc,2)))
      storelist[["up"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold > log(fc,2)))
      storelist[["down"]] <- storelist$all %>% dplyr::filter(FDR < fdr & (Fold < -log(fc,2)))
    }
    contrastlist[[paste0("contrast_",paste0(as.character(contrasts[i,][,c(4,2)]), collapse = "_"))]] <- storelist
  }
  diffblist[["contrasts"]] <- contrastlist
  return(diffblist)
}


diffbind_gene_anno <- function(assaytype, software, contrasts, genefile, path, outformat, fdr, fc){
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


run_chipseeker <- function(peaks_path, assaytype, txdb, org_db, peakformat){
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

read_genefile <- function(path, file){
  out <- read.table(paste0(path, file), header = F, stringsAsFactors = F)
  colnames(out) <- c("chr","start", "end", "strand", "type", "ensID", "gene")
  out["ens_gene"] <- paste0(out$ensID, "%", out$gene)
  return(out)
}

chipseeker_region_anno <- function(assaytype, software, contrasts, region_bed, path, txdb, org_db, regionformat, fdr, fc){
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

aggregate_feature <- function(annomatrix, columnstoaggregate, aggregatebycolumn, operation){
  storelist <- list()
  message("aggregating...")
  storelist[["aggregated"]] <- stats::aggregate(annomatrix[,columnstoaggregate], by=list(annomatrix[[aggregatebycolumn]]), operation)
  return(storelist)
}

maplot <- function(df, columnbyorder, fdr, fc){
  plotlist <- list()
  message("rearranging and plotting ma...")
  df <- df[,columnbyorder]
  colnames(df) <- c("baseMeanLog2","log2FoldChange","padj","feature")
  plotlist[["maplot"]] <- ggmaplot(df, fdr = fdr, fc = fc, size = 0.3, palette = c("#D22B2B", "#1465AC", "darkgray"),
                                   genenames = df$feature, legend = "top", top = 20,
                                   font.label = c("bold", 5), font.legend = "bold", font.main = "bold", ggtheme = ggplot2::theme_minimal())
}


quantify_featurecounts <- function(peaks_path, assaytype, bamfiles, sites_files, sites_type, sites_column_to_rearrange, pairedend=FALSE, refgenome, delim, merge_sites_files=FALSE){
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

quantify_mq_featurecounts <- function(peaks_path, assaytype, bamfiles, sites_files, sites_type, sites_column_to_rearrange, pairedend=FALSE, refgenome, delim, merge_sites_files=FALSE, mq){
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
                                                                                                 isPairedEnd = pairedend, nthreads = 12, minMQS = mq) # annot.ext will overide annot.inbuilt if both provided but to be on safe side I supply refgenome
  }
  return(storelist)
}

deseq2_de <- function(counts, coldata, fdr, fc){
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

deseq2_de_like_diffbind <- function(counts, coldata, fdr, fc, db_size_factor, fittype = "local", shrinkage="ashr", dars_analysis="no"){
  storelist <- list()
  
  ddsO <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = coldata, design = ~ condition)
  
  # filtering performed
  keep <- rowSums(counts(ddsO)) > 10
  ddsO <- ddsO[keep,]
  
  # de analysis
  message("Performing Deseq2 analysis...")
  
  # set size factor from diffbind
  sizeFactors(ddsO) <- db_size_factor
  
  message("estimate dispersion, fitType = ", fittype)
  ddsO <- estimateDispersions(ddsO, fitType = fittype)
  # DESeq2::nbinomWaldTest, with defaults
  ddsO <- nbinomWaldTest(ddsO)
  
  storelist[["ddsO"]]  <- ddsO
  
  # de analysis
  contrast_list <- list()
  ma_list <- list()
  if (dars_analysis == "yes"){
    for (cont1 in unique(coldata$condition)){
      for (cont2 in unique(coldata$condition)){
        if (cont1!=cont2){
          message("getting results for the contrast:")
          print(paste0(cont1,"_",cont2))
          
          db_res_temp <- DESeq2::results(storelist$ddsO, contrast=c("condition", cont1, cont2))
          
          message("performing shrinkage as diffbind does for logfoldchange, shrinkage type: ", shrinkage)
          db_res <- lfcShrink(ddsO, contrast=c("condition", cont1, cont2), res=db_res_temp, type=shrinkage) 
          db_res = db_res[order(rownames(db_res)),]
          db_res["feature_id"] <- rownames(db_res)
          db_resdf <- data.frame(db_res)
          db_resdf["feature_id"] <- rownames(db_resdf)
          message("Applying fdr and logfc cutoff...")
          # padj filter
          db_res$threshold <- as.logical(db_res$padj < fdr)
          db_res0.05 <- data.frame(db_res[which(db_res$threshold == TRUE),])
          # logfc filter
          db_res0.05_de <- db_res0.05 %>% dplyr::filter((log2FoldChange > log(fc,2)) | (log2FoldChange < -log(fc,2)))
          db_res0.05_up <- db_res0.05 %>% dplyr::filter((log2FoldChange > log(fc,2)))
          db_res0.05_down <- db_res0.05 %>% dplyr::filter((log2FoldChange < -log(fc,2)))
          contrast_list[[paste0(cont1,"_",cont2)]] <- list(all=db_resdf, de=db_res0.05_de, up=db_res0.05_up, down=db_res0.05_down)
          message("Plotting MA plot...")
          ma_list[[paste0(cont1,"_",cont2)]] <- ggmaplot(db_resdf, fdr = fdr, fc = fc, size = 0.3, palette = c("#D22B2B", "#1465AC", "darkgray"),
                                                         genenames = unlist(lapply(strsplit(db_resdf$feature_id, "%"), function(x) x[2])), legend = "top", top = 20, font.label = c("bold", 5),
                                                         font.legend = "bold", font.main = "bold", ggtheme = ggplot2::theme_minimal())
        }
      }
    }
  } else {
    contrast_list <- NA
    ma_list <- NA
  }
  
  storelist[["de_analysis"]] <- contrast_list
  storelist[["ma_plots"]] <- ma_list
  return(storelist)
}



consensus_gene_anno <- function(assaytype, software, contrasts, consensus_bed, genefile, path, outformat, fdr, fc){
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

# PCA
pcaplot <- function(counts, colorvec, dotvec){
  storelist <- list()
  countsfilt <- counts[rowSums(counts) >= 10,]
  tcountsfilt = t(countsfilt)
  tcountsfilt = data.frame(tcountsfilt)
  tcountsfilt["Color"] <-  rownames(tcountsfilt)
  dim(tcountsfilt)
  
  dfx <- tcountsfilt[,-c(dim(tcountsfilt)[2])]
  PC <- prcomp(dfx, scale. = T)
  PCi <- data.frame(PC$x,Color=tcountsfilt$Color)
  percentage <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
  percentage <- paste( colnames(PCi), "(", paste( as.character(percentage), "%", ")", sep="") )
  theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
  
  out <- ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
    theme + xlab(percentage[1]) + ylab(percentage[2])+
    geom_point(size=4,alpha=1,aes(shape=Color))+
    scale_color_manual(values = colorvec)+
    scale_shape_manual(values=dotvec) + theme_bw()
  storelist[["pca"]] <- out
  storelist[["countsfilt"]] <- countsfilt
  storelist[["transposed_countsfilt"]] <- tcountsfilt
  return(storelist)
}


# data is paired-end  as rectified from samtools view -H BAMfile

# diffbind
diffbind_atacseq_samplesheet <- read.csv("/mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg/samplesheet.csv")
# remove ME_1
diffbind_atacseq_samplesheet <- diffbind_atacseq_samplesheet[-c(14),]
atacseq_diffbind_list <- list()
# atacseq_diffbind_list[["prep"]][["dba_obj"]] <- dba(sampleSheet=diffbind_atacseq_samplesheet,  scoreCol=5, minOverlap=1)
# atacseq_diffbind_list[["prep"]][["ovp_rate"]] <- dba.overlap(atacseq_diffbind_list$prep$dba_obj, mode=DBA_OLAP_RATE)

atacseq_diffbind_list[["prep"]] <- diffbind_prep(diffbind_atacseq_samplesheet, 2) # #function(ss, minoverlap) 
plot(atacseq_diffbind_list$prep$dba_obj)
plot(atacseq_diffbind_list$prep$ovp_rate,type ='b',ylab='# peaks', xlab='Overlap at least this many peaksets')

# atacseq_diffbind_list[["counts"]][["dba_obj"]] <- dba.count(atacseq_diffbind_list$prep$dba_obj, minOverlap=1, score=DBA_SCORE_READS, summits=FALSE, mapQCth=0)
# atacseq_diffbind_list[["counts"]][["info"]]  <- dba.show(atacseq_diffbind_list$counts$dba_obj)

# mQCth = 15 is default
atacseq_diffbind_list[["counts"]] <- diffbind_count(atacseq_diffbind_list$prep$dba_obj, 2, 200) #function(obj, minoverlap, summit)

# Normalize
atacseq_diffbind_list[["norm"]] <- diffbind_norm(atacseq_diffbind_list$counts$dba_obj, 2) # function(obj, mimember)
# atacseq_diffbind_list[["norm"]] <- dba.normalize(atacseq_diffbind_list$counts$dba_obj)
# atacseq_diffbind_list$norm <- dba.contrast(atacseq_diffbind_list$norm, minMembers = 2,categories=DBA_CONDITION)
# atacseq_diffbind_list$norm <- dba.analyze(atacseq_diffbind_list$norm)

atacseq_diffbind_list$norm$norm$norm$DESeq2$norm.calc
atacseq_diffbind_list$norm$norm$norm$DESeq2$norm.method
atacseq_diffbind_list$norm$norm$norm$DESeq2$lib.calc
atacseq_diffbind_list$norm$norm$norm$DESeq2$lib.method
atacseq_diffbind_list$norm$norm$norm$DESeq2$lib.sizes
atacseq_diffbind_list$norm$norm$norm$DESeq2$norm.facs

# sizefactor diffbind
atacseq_diffbind_list$norm$norm$DESeq2$facs

# PCA
dba.plotPCA(atacseq_diffbind_list$norm$norm, attributes=DBA_CONDITION, label=DBA_REPLICATE, method= DBA_DESEQ2, score=DBA_SCORE_READS)

dba.plotPCA(atacseq_diffbind_list$norm$norm, attributes=DBA_CONDITION, label=DBA_REPLICATE, method= DBA_DESEQ2, score=DBA_SCORE_NORMALIZED)


# Differential analysis with counts
atacseq_diffbind_list[["dar_analysis"]][["deseq2"]] <- diffbind_de(atacseq_diffbind_list$norm$norm, 0.05, 2, "deseq2")


# since this data is also hg38 the file gene_gencode_human_upstream_2kb.sorted.txt (promoter (2kb) + gene coordinates) can be used
atacseq_diffbind_list[["gene"]] <- read_genefile("/mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg/", "gene_gencode_human_upstream_2kb.sorted.txt")

atacseq_diffbind_list[["annotation"]][["deseq2"]] <- diffbind_gene_anno("atacseq", "diffbind",atacseq_diffbind_list$dar_analysis$deseq2$contrasts, atacseq_diffbind_list$gene, "/mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg/", "txt", 0.05, 4)

# since we need only deseq2 results for downstream analysis, I perform mean analysis with deseq2
# comparison required
# contrast_hESC_PreME
# contrast_PreME_ME
# contrast_ME_DE
# contrast_PreME_hPGCLC_d2
# contrast_hPGCLC_d2_hPGCLC_d4
# contrast_hPGCLC_d4_hPGC_wk7_M
# contrast_hPGC_wk7_M_hPGC_wk7_F
# contrast_hPGC_wk7_M_Soma_wk7_M
# contrast_hPGCLC_d2_DE

atacseq_diffbind_list[["dar_analysis"]][["deseq2"]][["maplot"]][["contrast_PreME_hESC"]] <- maplot(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_hESC$all, c(6,9,11,12),0.05, 2)
atacseq_diffbind_list[["dar_analysis"]][["deseq2"]][["maplot"]][["contrast_PreME_ME"]] <- maplot(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_ME$all, c(6,9,11,12),0.05, 2)
atacseq_diffbind_list[["dar_analysis"]][["deseq2"]][["maplot"]][["contrast_ME_DE"]] <- maplot(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_ME_DE$all, c(6,9,11,12),0.05, 2)
atacseq_diffbind_list[["dar_analysis"]][["deseq2"]][["maplot"]][["contrast_PreME_hPGCLC_d2"]] <- maplot(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_hPGCLC_d2$all, c(6,9,11,12),0.05, 2)
atacseq_diffbind_list[["dar_analysis"]][["deseq2"]][["maplot"]][["contrast_hPGCLC_d4_hPGCLC_d2"]] <- maplot(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGCLC_d4_hPGCLC_d2$all, c(6,9,11,12),0.05, 2)
atacseq_diffbind_list[["dar_analysis"]][["deseq2"]][["maplot"]][["contrast_hPGC_wk7_M_hPGCLC_d4"]] <- maplot(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGC_wk7_M_hPGCLC_d4$all, c(6,9,11,12),0.05, 2)
atacseq_diffbind_list[["dar_analysis"]][["deseq2"]][["maplot"]][["contrast_hPGC_wk7_M_hPGC_wk7_F"]] <- maplot(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGC_wk7_M_hPGC_wk7_F$all, c(6,9,11,12),0.05, 2)
atacseq_diffbind_list[["dar_analysis"]][["deseq2"]][["maplot"]][["contrast_Soma_wk7_M_hPGC_wk7_M"]] <- maplot(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_Soma_wk7_M_hPGC_wk7_M$all, c(6,9,11,12),0.05, 2)
atacseq_diffbind_list[["dar_analysis"]][["deseq2"]][["maplot"]][["contrast_hPGCLC_d2_DE"]] <- maplot(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGCLC_d2_DE$all, c(6,9,11,12),0.05, 2)

dim(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_hESC$de)
dim(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_ME$de)
dim(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_ME_DE$de)
dim(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_hPGCLC_d2$de)
dim(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGCLC_d4_hPGCLC_d2$de)
dim(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGC_wk7_M_hPGCLC_d4$de)
dim(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGC_wk7_M_hPGC_wk7_F$de)
dim(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_Soma_wk7_M_hPGC_wk7_M$de)
dim(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGCLC_d2_DE$de)

write.table(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_hESC$de, "contrast_PreME_hESC_de.bed", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_ME$de, "contrast_PreME_ME_de.bed" , sep = "\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_ME_DE$de, "contrast_ME_DE_de.bed" , sep = "\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_hPGCLC_d2$de, "contrast_PreME_hPGCLC_d2_de.bed" , sep = "\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGCLC_d4_hPGCLC_d2$de, "contrast_hPGCLC_d4_hPGCLC_d2_de.bed" , sep = "\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGC_wk7_M_hPGCLC_d4$de, "contrast_hPGC_wk7_M_hPGCLC_d4_de.bed" , sep = "\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGC_wk7_M_hPGC_wk7_F$de, "contrast_hPGC_wk7_M_hPGC_wk7_F_de.bed" , sep = "\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_Soma_wk7_M_hPGC_wk7_M$de, "contrast_Soma_wk7_M_hPGC_wk7_M_de.bed" , sep = "\t", quote = F, append = F, row.names = F, col.names = F)
write.table(atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGCLC_d2_DE$de, "contrast_hPGCLC_d2_DE_de.bed" , sep = "\t", quote = F, append = F, row.names = F, col.names = F)

atacseq_diffbind_matrix <- matrix(0, nrow=length(names(atacseq_diffbind_list$dar_analysis$deseq2$contrasts)), ncol=2)

for (i in names(atacseq_diffbind_list$dar_analysis$deseq2$contrasts)){
  print(paste0(gsub("contrast_","",i),"  ",as.numeric(length(atacseq_diffbind_list$dar_analysis$deseq2$contrasts[[i]]$de$feature_id)[1])))
  atacseq_diffbind_matrix[which(i == names(atacseq_diffbind_list$dar_analysis$deseq2$contrasts)), 1] <- gsub("contrast_","",i)
  atacseq_diffbind_matrix[which(i == names(atacseq_diffbind_list$dar_analysis$deseq2$contrasts)), 2] <- as.numeric(length(atacseq_diffbind_list$dar_analysis$deseq2$contrasts[[i]]$de$feature_id)[1])
}

atacseq_diffbind_matrix_df <- data.frame(atacseq_diffbind_matrix)
colnames(atacseq_diffbind_matrix_df) <- c("sample", "dars")
atacseq_diffbind_matrix_df$dars <- as.numeric(atacseq_diffbind_matrix_df$dars)

# upset plot all
# no elelement in

atacseq_upset_list_diffbind <- list(PreME_hESC = atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_hESC$de$feature_id,
                                    PreME_ME = atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_ME$de$feature_id,
                                    ME_DE = atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_ME_DE$de$feature_id,
                                    PreME_hPGCLC_d2 = atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_PreME_hPGCLC_d2$de$feature_id,
                                    hPGCLC_d4_hPGCLC_d2 = atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGCLC_d4_hPGCLC_d2$de$feature_id,
                                    hPGC_wk7_M_hPGCLC_d4 = atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGC_wk7_M_hPGCLC_d4$de$feature_id,
                                    hPGC_wk7_M_hPGC_wk7_F = atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGC_wk7_M_hPGC_wk7_F$de$feature_id,
                                    Soma_wk7_M_hPGC_wk7_M = atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_Soma_wk7_M_hPGC_wk7_M$de$feature_id,
                                    hPGCLC_d2_DE = atacseq_diffbind_list$dar_analysis$deseq2$contrasts$contrast_hPGCLC_d2_DE$de$feature_id)

upset(fromList(atacseq_upset_list_diffbind), order.by = "freq")

# --------------------------  END OF ANALYSIS --------------------------------------- #

# atacseq_diffbind_list[["chipseeker"]][["annotation"]][["deseq2"]] <- chipseeker_region_anno("atacseq", "diffbind",atacseq_diffbind_list$dar_analysis$deseq2$contrasts, atacseq_diffbind_list$regions, "/mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2) 

# 
# # deseq2
# # prepare the consensus peaks file using bedtools
# # run the following in /mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg/ATACseq/macs2_narrow
# system("cat /mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg/ATACseq/macs2_narrow/*.narrowPeak > all_Peaks.bed")
# system("sort -k1,1 -k2,2n all_Peaks.bed > all_Peaks_sorted.bed")
# system("bedtools merge -i all_Peaks_sorted.bed > merged_Peaks.bed")
# merged_Peaks <- data.frame(fread("merged_Peaks.bed"))
# colnames(merged_Peaks) <- c("chr", "start", "end")
# merged_Peaks["peakid"] <- paste0("consensus_", rownames(merged_Peaks))
# merged_Peaks["featureid"] <- paste0(merged_Peaks$chr, "%", merged_Peaks$start, "%", merged_Peaks$end)
# merged_Peaks["peakid_featureid"] <- paste0(merged_Peaks$peakid, "_", merged_Peaks$featureid)
# merged_Peaks["score"] <- 0
# merged_Peaks["strand"] <- "."
# 
# write.table(merged_Peaks[,c(1:5,8)], "ATACseq/bam_clean/merged_Peaks_sorted.bed", sep="\t", append = F, quote = F, col.names = F, row.names = F)
# write.table(merged_Peaks[,c(1:3,6:8)], "consensus_Peaks_sorted.bed", sep="\t", append = F, quote = F, col.names = F, row.names = F)
# 
# # Extract consensus peaks obtained from diffbind (note: only the consensus peak coordinates) # peak coordinates will be as per diffbind
# head(dba.peakset(atacseq_diffbind_list$prep$dba_obj, bRetrieve=TRUE, DataType=DBA_DATA_FRAME))
# merged_Peaks_db <- data.frame(dba.peakset(atacseq_diffbind_list$prep$dba_obj, bRetrieve=TRUE, DataType=DBA_DATA_FRAME))[,c(1:3)]
# colnames(merged_Peaks_db) <- c("chr", "start", "end")
# merged_Peaks_db["peakid"] <- paste0("consensus_", rownames(merged_Peaks_db))
# merged_Peaks_db["featureid"] <- paste0(merged_Peaks_db$chr, "%", merged_Peaks_db$start, "%", merged_Peaks_db$end)
# merged_Peaks_db["peakid_featureid"] <- paste0(merged_Peaks_db$peakid, "_", merged_Peaks_db$featureid)
# merged_Peaks_db["score"] <- 0
# merged_Peaks_db["strand"] <- "."
# 
# write.table(merged_Peaks_db[,c(1:5,8)], "ATACseq/bam_clean/merged_Peaks_db_sorted.bed", sep="\t", append = F, quote = F, col.names = F, row.names = F)
# write.table(merged_Peaks_db[,c(1:3,6:8)], "consensus_Peaks_from_db_sorted.bed", sep="\t", append = F, quote = F, col.names = F, row.names = F)
# 
# 
# # empty list store deseq2 analysis
# atacseq_quantify_list <- list()
# 
# atacseq_quantify_list[["featurecounts"]][["consensus"]] <- quantify_mq_featurecounts("/mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg/ATACseq/bam_clean/", "atacseq","\\.bam$", "\\_Peaks_db_sorted.bed$","consensus", c(4,5,1:3,6), pairedend=TRUE, "hg38", "_",merge_sites_files=FALSE, 15)
# head(atacseq_quantify_list$featurecounts$consensus$consensus_countmatrix$merged_Peaks_db_sorted.bed$counts, 2)
# 
# colnames(atacseq_quantify_list$featurecounts$consensus$consensus_countmatrix$merged_Peaks_db_sorted.bed$counts) <- gsub("_clean.bam", "", gsub("ATAC_", "", colnames(atacseq_quantify_list$featurecounts$consensus$consensus_countmatrix$merged_Peaks_db_sorted.bed$counts)))
# 
# atacseq_quantify_list$featurecounts$consensus$consensus_countmatrix$merged_Peaks_db_sorted.bed$counts <- atacseq_quantify_list$featurecounts$consensus$consensus_countmatrix$merged_Peaks_db_sorted.bed$counts[,c(1:4, 10:13,5:9,14:20 )]
# 
# # empty list to store deseq2 output
# atacseq_deseq2_list <- list()
# 
# atacseq_deseq2_list[["counts"]] <- atacseq_quantify_list$featurecounts$consensus$consensus_countmatrix$merged_Peaks_db_sorted.bed$counts
# atacseq_deseq2_list[["coldata"]] <- data.frame(fread("coldata.txt", header = T))
# rownames(atacseq_deseq2_list$coldata) <- atacseq_deseq2_list$coldata$SampleID
# atacseq_deseq2_list$coldata$condition <- factor(atacseq_deseq2_list$coldata$Condition)
# atacseq_deseq2_list$coldata$replicate <- factor(atacseq_deseq2_list$coldata$Replicate)
# atacseq_deseq2_list$coldata <- atacseq_deseq2_list$coldata[,c(1,2,5,6)]
# all(rownames(atacseq_deseq2_list$coldata) == colnames(atacseq_deseq2_list$counts)) #should print TRUE
# 
# 
# # pre checks
# atacseq_deseq2_list[["prechecks"]][["atacseq_ddsO"]] <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(atacseq_deseq2_list$counts), colData = atacseq_deseq2_list$coldata, design = ~ condition)
# atacseq_deseq2_list[["prechecks"]][["atacseq_keep"]] <- rowSums(counts(atacseq_deseq2_list$prechecks$atacseq_ddsO)) > 10
# atacseq_deseq2_list[["prechecks"]][["atacseq_ddsO"]] <- atacseq_deseq2_list$prechecks$atacseq_ddsO[atacseq_deseq2_list$prechecks$atacseq_keep,]
# atacseq_deseq2_list[["prechecks"]][["atacseq_vstO"]] <- DESeq2::vst(atacseq_deseq2_list$prechecks$atacseq_ddsO) 
# atacseq_deseq2_list[["prechecks"]][["pcaData"]] <- DESeq2::plotPCA(atacseq_deseq2_list$prechecks$atacseq_vstO, intgroup=c("condition", "replicate"), returnData=TRUE)
# atacseq_deseq2_list[["prechecks"]][["percentVar"]] <- round(100 * attr(atacseq_deseq2_list$prechecks$pcaData, "percentVar"))
# ggplot(atacseq_deseq2_list$prechecks$pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",atacseq_deseq2_list$rechecks$percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",atacseq_deseq2_list$rechecks$percentVar[2],"% variance")) +
#   coord_fixed() + theme_bw()
# 
# 
# atacseq_deseq2_list[["prechecks"]][["norm"]] <- deseq2_de_like_diffbind(atacseq_deseq2_list$counts, atacseq_deseq2_list$coldata, 0.05, 2, atacseq_diffbind_list$norm$norm$DESeq2$facs, fittype = "local", shrinkage="ashr", dars_analysis="no")
# 
# 
# # raw pca
# atacseq_deseq2_list[["prechecks"]][["full_d_pca_res"]]  <- prcomp(t(counts(atacseq_deseq2_list$prechecks$norm$ddsO)), center = TRUE, scale. = TRUE)
# factoextra::fviz_pca_ind(atacseq_deseq2_list[["prechecks"]][["full_d_pca_res"]], geom = c("point", "text"), repel = TRUE, labelsize =4)
# 
# # normalized pca
# atacseq_deseq2_list[["prechecks"]][["full_d_pca_resn"]] <- prcomp(t(counts(atacseq_deseq2_list$prechecks$norm$ddsO, normalized =T)), center = TRUE, scale. = TRUE)
# factoextra::fviz_pca_ind(atacseq_deseq2_list[["prechecks"]][["full_d_pca_resn"]], geom = c("point", "text"), repel = TRUE, labelsize = 4)
# 
# atacseq_deseq2_list[["prechecks"]][["colorvector"]] <-  c("blue","blue","magenta","magenta","green","green","green","orange","orange","#0090B6","#0090B6","darkgrey","darkgrey","brown", "brown","brown","black","black", "red", "red")
# atacseq_deseq2_list[["prechecks"]][["dotvector"]] <-  c(1,2,1,2,1,2,3,1,2,1,2,1,2,1,2,3,1,2,1,2)
# atacseq_deseq2_list[["prechecks"]][["raw_pca"]] <- pcaplot(counts(atacseq_deseq2_list$prechecks$norm$ddsO), atacseq_deseq2_list$prechecks$colorvector, atacseq_deseq2_list$prechecks$dotvector)
# atacseq_deseq2_list[["prechecks"]][["norm_pca"]] <- pcaplot(counts(atacseq_deseq2_list$prechecks$norm$ddsO, normalized =T), atacseq_deseq2_list$prechecks$colorvector, atacseq_deseq2_list$prechecks$dotvector)
# 
# 
# # # ME_1 need to be removed ..
# # # rechecks
# # atacseq_deseq2_list$counts <- atacseq_deseq2_list$counts[,-c(14)]
# # atacseq_deseq2_list$coldata <- atacseq_deseq2_list$coldata[-c(14),]
# # atacseq_deseq2_list[["rechecks"]][["atacseq_ddsO"]] <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(atacseq_deseq2_list$counts), colData = atacseq_deseq2_list$coldata, design = ~ condition)
# # atacseq_deseq2_list[["rechecks"]][["atacseq_keep"]] <- rowSums(counts(atacseq_deseq2_list$rechecks$atacseq_ddsO)) > 10
# # atacseq_deseq2_list[["rechecks"]][["atacseq_ddsO"]] <- atacseq_deseq2_list$rechecks$atacseq_ddsO[atacseq_deseq2_list$rechecks$atacseq_keep,]
# # atacseq_deseq2_list[["rechecks"]][["atacseq_vstO"]] <- DESeq2::vst(atacseq_deseq2_list$rechecks$atacseq_ddsO) 
# # atacseq_deseq2_list[["rechecks"]][["pcaData"]] <- DESeq2::plotPCA(atacseq_deseq2_list$rechecks$atacseq_vstO, intgroup=c("condition", "replicate"), returnData=TRUE)
# # atacseq_deseq2_list[["rechecks"]][["percentVar"]] <- round(100 * attr(atacseq_deseq2_list$rechecks$pcaData, "percentVar"))
# # ggplot(atacseq_deseq2_list$rechecks$pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
# #   geom_point(size=3) +
# #   xlab(paste0("PC1: ",atacseq_deseq2_list$rechecks$percentVar[1],"% variance")) +
# #   ylab(paste0("PC2: ",atacseq_deseq2_list$rechecks$percentVar[2],"% variance")) + 
# #   coord_fixed() + theme_bw()
# # 
# # de analysis
# atacseq_deseq2_list[["dar_analysis"]] <- deseq2_de_like_diffbind(atacseq_deseq2_list$counts, atacseq_deseq2_list$coldata, 0.05, 2, atacseq_diffbind_list$norm$norm$DESeq2$facs, fittype = "local", shrinkage="ashr", dars_analysis="yes")
# 
# atacseq_deseq2_matrix <- matrix(0, nrow=length(names(atacseq_deseq2_list$dar_analysis$de_analysis)), ncol=2)
# 
# for (i in names(atacseq_deseq2_list$dar_analysis$de_analysis)){
#   print(paste0(gsub("contrast_","",i),"  ",as.numeric(length(atacseq_deseq2_list$dar_analysis$de_analysis[[i]]$de$feature_id)[1])))
#   atacseq_deseq2_matrix[which(i == names(atacseq_deseq2_list$dar_analysis$de_analysis)), 1] <- gsub("contrast_","",i)
#   atacseq_deseq2_matrix[which(i == names(atacseq_deseq2_list$dar_analysis$de_analysis)), 2] <- as.numeric(length(atacseq_deseq2_list$dar_analysis$de_analysis[[i]]$de$feature_id)[1])
# }
# 
# atacseq_deseq2_matrix_df <- data.frame(atacseq_deseq2_matrix)
# colnames(atacseq_deseq2_matrix_df) <- c("sample", "dars")
# atacseq_deseq2_matrix_df$dars <- as.numeric(atacseq_deseq2_matrix_df$dars)
# 
# 
# atacseq_deseq2_list[["gene"]] <- read_genefile("/mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg/", "gene_gencode_human_upstream_2kb.sorted.txt")
# 
# atacseq_deseq2_list[["consensus"]] <- data.frame(fread("/mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg/consensus_Peaks_sorted.bed"))
# 
# # get positions of dars
# atacseq_deseq2_list[["annotation"]] <- consensus_gene_anno("atacseq", "deseq2", atacseq_deseq2_list$dar_analysis$de_analysis, atacseq_deseq2_list$consensus, atacseq_deseq2_list$gene, "/mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg", "txt", 0.05, 2)
# 
# atacseq_deseq2_list[["chipseeker"]][["annotation"]] <- chipseeker_region_anno("atacseq", "deseq2", atacseq_deseq2_list$dar_analysis$de_analysis, atacseq_deseq2_list$consensus[,c(1:4)], "/mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2)
# 
# # 
# # write.table(atacseq_deseq2_list$annotation$hESC_PreME$all$de[,c(8:10,3,7,1,11:12)] %>% distinct(), "atacseq_deseq2_hESC_PreME_de.bed", sep="\t", quote = F, append = F, row.names = F, col.names = F)
# 
# # intersect dars from deseq2 and diffbind
# atacseq_deseq2_diffbind_int_matrix <- matrix(0, nrow=36, ncol=2)
# 
# counter = 0
# for (i in names(atacseq_diffbind_list$dar_analysis$deseq2$contrasts)){
#   for (j in names(atacseq_deseq2_list$dar_analysis$de_analysis)){
#     if (gsub("contrast_","",i) == j){
#       counter = counter + 1
#       print(paste0(counter," ", gsub("contrast_","",i)," ", j ," ", length(intersect(sapply(strsplit(atacseq_deseq2_list$dar_analysis$de_analysis[[j]]$de$feature_id, "_"), function(x) x[3]), atacseq_diffbind_list$dar_analysis$deseq2$contrasts[[i]]$de$feature_id))))
#       atacseq_deseq2_diffbind_int_matrix[counter, 1] <- j
#       atacseq_deseq2_diffbind_int_matrix[counter, 2] <- length(intersect(sapply(strsplit(atacseq_deseq2_list$dar_analysis$de_analysis[[j]]$de$feature_id, "_"), function(x) x[3]), atacseq_diffbind_list$dar_analysis$deseq2$contrasts[[i]]$de$feature_id))
#     }
#   }
# }
# 
# atacseq_deseq2_diffbind_int_matrix_df <- data.frame(atacseq_deseq2_diffbind_int_matrix)
# colnames(atacseq_deseq2_diffbind_int_matrix_df) <- c("sample", "dars")
# atacseq_deseq2_diffbind_int_matrix_df$dars <- as.numeric(atacseq_deseq2_diffbind_int_matrix_df$dars)
# 
# # merge deseq2 and diffbind
# atacseq_deseq2_diffbind_matrix_df <- merge(merge(atacseq_deseq2_matrix_df, atacseq_diffbind_matrix_df, by="sample"), atacseq_deseq2_diffbind_int_matrix_df, by="sample")
# 
# 
# 





