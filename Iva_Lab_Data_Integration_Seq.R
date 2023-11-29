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

fn_sample_filter <- function(countmatrix, columnsremove){
  countmatrix <- countmatrix[,-columnsremove] # columns is vector eg, c(7:9)
  countmatrix <- data.frame(countmatrix)
  countmatrix["ensID"] <- rownames(countmatrix)
  return(countmatrix)
}

fn_reassign_genenames <- function(countmatrix, genedf, columnskeep){
  genecounts <- merge(countmatrix, genedf ,by.x = "ensID", by.y ="ensID", all.y =T)
  rownames(genecounts) <- genecounts$ens_gene
  genecounts <- genecounts[,columnskeep]
  return(genecounts)
}

fn_make_coldata_dds <- function(dds, rowsremove){
  coldata <- data.frame(colData(dds))
  colnames(coldata) <- c("sample", "condition", "replicate", "sizeFactor")
  coldata <- coldata[-rowsremove,] # columns is vector eg, c(7:9)
  coldata$condition <- factor(coldata$condition)
  coldata$replicate <- factor(coldata$replicate)
  return(coldata)
}


fn_deseq2 <- function(counts, coldata, fdr, fc, contrasts){
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
        rnaseqde <- DESeq2::results(storelist$ddsO, contrast=c("condition", cont1, cont2))
        rnaseqde = rnaseqde[order(rownames(rnaseqde)),]
        rnaseqde["ensID_Gene"] <- rownames(rnaseqde)
        rnaseqdedf <- data.frame(rnaseqde)
        rnaseqdedf["ensID_Gene"] <- rownames(rnaseqdedf)
        message("Applying fdr and logfc cutoff...")
        # padj filter
        rnaseqde$threshold <- as.logical(rnaseqde$padj < fdr)
        rnaseqde0.05 <- data.frame(rnaseqde[which(rnaseqde$threshold == TRUE),])
        # logfc filter
        rnaseqde0.05_de <- rnaseqde0.05 %>% dplyr::filter((log2FoldChange > log(fc,2)) | (log2FoldChange < -log(fc,2)))
        rnaseqde0.05_up <- rnaseqde0.05 %>% dplyr::filter((log2FoldChange > log(fc,2)))
        rnaseqde0.05_down <- rnaseqde0.05 %>% dplyr::filter((log2FoldChange < -log(fc,2)))
        contrast_list[[paste0(cont1,"_",cont2)]] <- list(all=rnaseqdedf, de=rnaseqde0.05_de, up=rnaseqde0.05_up, down=rnaseqde0.05_down)
        message("Plotting MA plot...")
        ma_list[[paste0(cont1,"_",cont2)]] <- ggmaplot(rnaseqdedf, fdr = fdr, fc = fc, size = 0.3, palette = c("#D22B2B", "#1465AC", "darkgray"),
                 genenames = unlist(lapply(strsplit(rnaseqdedf$ensID_Gene, "%"), function(x) x[2])), legend = "top", top = 20, font.label = c("bold", 5),
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
    }
  }
  return(contrastlist)
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
    storelist[["plotAnnoBar"]] <- plotAnnoBar(storelist$gr_anno)
    storelist[["plotDistToTSS"]] <- plotDistToTSS(storelist$gr_anno)
    message("fetaching promoter...")
    promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
    message("generating tagmatrix...")
    storelist[["tagMatrixList"]] <- lapply(storelist[["gr_obj"]], getTagMatrix, windows=promoter)
    message("plotting average plot...")
    storelist[["plotAvgProf"]] <- plotAvgProf(storelist[["tagMatrixList"]], xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
    ## binning method
    storelist[["plotPeakProf2"]] <- plotPeakProf2(storelist[["gr_obj"]], upstream = 3000, downstream = 3000, conf = 0.95,by = "gene", type = "start_site", TxDb = txdb, facet = "row", nbin = 800)
  }
  return(storelist)
}


################################################
###       RNA-Seq data: Knockdown            ###
################################################
rnaseqkd_path <- "/mnt/home3/reid/av638/rnaseq/iva_lab_oct23/"

rnaseqkd_list <- list()

rnaseqkd_list[["dds"]] <- fn_load_rdata(rnaseqkd_path, "outfolder/star_salmon/deseq2_qc/deseq2.dds.RData")
rnaseqkd_list[["gene"]] <- fn_read_genefile(rnaseqkd_path, "gene_gencode_human__gencode_out_chr.txt")
rnaseqkd_list[["counts"]] <- fn_dds_counts(rnaseqkd_list$dds)
rnaseqkd_list[["samples_selected"]] <- fn_sample_filter(rnaseqkd_list$counts, c(7:9))
rnaseqkd_list[["processed_counts"]] <- fn_reassign_genenames(rnaseqkd_list$samples_selected, rnaseqkd_list$gene, c(2:10))
rnaseqkd_list[["coldata"]] <- fn_make_coldata_dds(rnaseqkd_list$dds, c(7:9))

# check
all(rownames(rnaseqkd_list$coldata) == colnames(rnaseqkd_list$processed_counts)) #should print TRUE

# de analysis
rnaseqkd_list[["deseq2"]] <- fn_deseq2(rnaseqkd_list$processed_counts, rnaseqkd_list$coldata, 0.05, 2, rnaseqkd_contrasts)

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

atacseqkd_edger_list[["chipseeker"]] <- fn_run_chipseeker("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/", "atacseqkd", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "*.broadPeak")
  
#  export bed files
fn_export_bed("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", atacseqkd_edger_list$annotation, "all", "all", "txt")



#####==================================#####

# ATAC-Seq data

load("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/deseq2/consensus_peaks.mLb.clN.dds.RData")

atacseqkd_dds <- dds
rm(dds)
atacseqkd_precountdata <- counts(atacseqkd_dds, normalized=FALSE)
head(atacseqkd_precountdata)
dim(atacseqkd_precountdata)
colSums(atacseqkd_precountdata)
atacseqkd_precountdata <- data.frame(atacseqkd_precountdata)
write.table(atacseqkd_precountdata, "atacseqkd_precountdata.txt", quote = F, append = F, sep="\t")

atacseqkd_vstNorm <- DESeq2::vst(atacseqkd_dds) 
atacseqkd_vstNorm_pca <- plotPCA(atacseqkd_vstNorm, intgroup="sample", returnData=TRUE)
plotPCA(atacseqkd_vstNorm, intgroup="Group1", returnData=FALSE)


# get countdata
atacseqkd_countdata <- atacseqkd_precountdata[,c(1:dim(atacseqkd_precountdata)[2])]


# create sample-specific coldata object
atacseqkd_dds_coldata <- data.frame(colData(atacseqkd_dds))
colnames(atacseqkd_dds_coldata) <- c("sample", "condition", "replicate", "sizeFactor")
atacseqkd_dds_coldata$condition <- factor(atacseqkd_dds_coldata$condition)
atacseqkd_dds_coldata$replicate <- factor(atacseqkd_dds_coldata$replicate)

# prepare count matrix 
atacseqkd_countdata = as.matrix(atacseqkd_countdata) 
all(rownames(atacseqkd_dds_coldata) == colnames(atacseqkd_countdata)) #should print TRUE
atacseqkd_ddsO <- DESeqDataSetFromMatrix(countData = atacseqkd_countdata, colData = atacseqkd_dds_coldata, design = ~ condition)
atacseqkd_keep <- rowSums(counts(atacseqkd_ddsO)) > 10
atacseqkd_ddsO <- atacseqkd_ddsO[atacseqkd_keep,]

# pca 
atacseqkd_vst0Norm <- DESeq2::vst(atacseqkd_ddsO)
atacseqkd_vst0Norm_pca <- plotPCA(atacseqkd_vst0Norm, intgroup="condition", returnData=TRUE)
plotPCA(atacseqkd_vst0Norm, intgroup="condition", returnData=FALSE)

# de analysis
atacseqkd_ddsO <- DESeq(atacseqkd_ddsO)


# shCRAMP1 vs shControl
atacseqkd_de_shCRAMP1_shControl <- results(atacseqkd_ddsO, contrast=c("condition", "shCRAMP1", "shControl"))
atacseqkd_de_shCRAMP1_shControl = atacseqkd_de_shCRAMP1_shControl[order(rownames(atacseqkd_de_shCRAMP1_shControl)),]
atacseqkd_de_shCRAMP1_shControl$threshold <- as.logical(atacseqkd_de_shCRAMP1_shControl$padj < 0.05)
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05 <- atacseqkd_de_shCRAMP1_shControl[which(atacseqkd_de_shCRAMP1_shControl$threshold == TRUE),]
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_df <- data.frame(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05)
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_df["interval"] <- rownames(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_df)
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_df <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_df[,c(8,1:7)]
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_df <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5) | (log2FoldChange < -1.5))
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_up <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5))
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_down <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange < -1.5))


# log2(1.5) = 2.83
ggmaplot(atacseqkd_de_shCRAMP1_shControl,
         fdr = 0.05, fc = 2.83, size = 0.3,
         palette = c("#D22B2B", "#1465AC", "darkgray"),
         genenames = rownames(atacseqkd_de_shCRAMP1_shControl),
         legend = "top", top = 20,
         font.label = c("bold", 5),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())

# shSUZ12 vs shControl
atacseqkd_de_shSUZ12_shControl <- results(atacseqkd_ddsO, contrast=c("condition", "shSUZ12", "shControl"))
atacseqkd_de_shSUZ12_shControl = atacseqkd_de_shSUZ12_shControl[order(rownames(atacseqkd_de_shSUZ12_shControl)),]
atacseqkd_de_shSUZ12_shControl$threshold <- as.logical(atacseqkd_de_shSUZ12_shControl$padj < 0.05)
deseq2_atacseqkd_de_shSUZ12_shControl_0.05 <- atacseqkd_de_shSUZ12_shControl[which(atacseqkd_de_shSUZ12_shControl$threshold == TRUE),]
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_df <- data.frame(deseq2_atacseqkd_de_shSUZ12_shControl_0.05)
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_df["interval"] <- rownames(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_df)
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_df <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_df[,c(8,1:7)]
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_df <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5) | (log2FoldChange < -1.5))
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_up <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5))
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_down <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange < -1.5))


ggmaplot(atacseqkd_de_shSUZ12_shControl,
         fdr = 0.05, fc = 2.83, size = 0.3,
         palette = c("#D22B2B", "#1465AC", "darkgray"),
         genenames = rownames(atacseqkd_de_shCRAMP1_shControl),
         legend = "top", top = 20,
         font.label = c("bold", 5),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())


# annotate 
atacseqkd_consensus_peaks.mLb.clN.annotatePeaks <- data.frame(fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.annotatePeaks.txt", header = T))
colnames(atacseqkd_consensus_peaks.mLb.clN.annotatePeaks) <- c("interval",colnames(atacseqkd_consensus_peaks.mLb.clN.annotatePeaks)[2:length(atacseqkd_consensus_peaks.mLb.clN.annotatePeaks)])


deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann <- merge(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_df, atacseqkd_consensus_peaks.mLb.clN.annotatePeaks, by="interval")
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann <- merge(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_df, atacseqkd_consensus_peaks.mLb.clN.annotatePeaks, by="interval")
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_up_ann <- merge(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_up, atacseqkd_consensus_peaks.mLb.clN.annotatePeaks, by="interval")
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_down_ann <- merge(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_down, atacseqkd_consensus_peaks.mLb.clN.annotatePeaks, by="interval")
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_up_ann <- merge(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_up, atacseqkd_consensus_peaks.mLb.clN.annotatePeaks, by="interval")
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_down_ann <- merge(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_down, atacseqkd_consensus_peaks.mLb.clN.annotatePeaks, by="interval")

write.xlsx(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann, "deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann.xlsx")
write.xlsx(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann, "deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann.xlsx")

deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann_nearest <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann 
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann_nearest <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann 
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_up_ann_nearest <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_up_ann 
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_down_ann_nearest <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_down_ann 
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_up_ann_nearest <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_up_ann 
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_down_ann_nearest <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_down_ann 

#Use manual annotation
sort -k1,1 -k2,2n consensus_peaks.mLb.clN.bed | grep chr > consensus_peaks.mLb.clN.sorted.chr.bed
sort -k1,1 -k2,2n /mnt/home3/reid/av638/rnaseq/iva_lab_oct23/gene_gencode_human__gencode_out.txt | grep chr > gene_gencode_human_gencode_out.sorted.chr.txt
bedtools closest -a consensus_peaks.mLb.clN.sorted.chr.bed -b gene_gencode_human_gencode_out.sorted.chr.txt -d > consensus_peaks.mLb.clN.sorted.chr.ann.bed

atacseqkd_consensus_peaks.gene <- data.frame(fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.sorted.chr.ann.bed", header = F))
atacseqkd_consensus_peaks.gene <- atacseqkd_consensus_peaks.gene[,c(4,12,13,14)]
colnames(atacseqkd_consensus_peaks.gene) <- c("interval", "ensID","Gene", "Distance")
atacseqkd_consensus_peaks.gene["ensID_Gene"] <- paste0(atacseqkd_consensus_peaks.gene$ensID, "%", atacseqkd_consensus_peaks.gene$Gene)

deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_anno <- merge(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_df, atacseqkd_consensus_peaks.gene, by="interval")
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_anno <- merge(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_df, atacseqkd_consensus_peaks.gene, by="interval")
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_up_anno <- merge(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_up, atacseqkd_consensus_peaks.gene, by="interval")
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_down_anno <- merge(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_down, atacseqkd_consensus_peaks.gene, by="interval")
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_up_anno <- merge(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_up, atacseqkd_consensus_peaks.gene, by="interval")
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_down_anno <- merge(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_down, atacseqkd_consensus_peaks.gene, by="interval")

# For extracting DARs coordinates HOMER outputs can serve thd purpose, only peak coordinates which are DARs need to be extracted
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_annotated <- merge(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_df, atacseqkd_consensus_peaks.mLb.clN.annotatePeaks, by="interval")
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_annotated <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_annotated[,c(9:11,1)]
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_annotated_ext <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_annotated
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_annotated_ext$Start <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_annotated_ext$Start - 1000
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_annotated_ext$End <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_annotated_ext$End + 1000

write.table(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_annotated_ext, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_annotated_ext.bed", sep = "\t", append = F, quote = F, row.names = F, col.names = F)

# Take only chr  atacseq peaks i.e. consensus_peaks.mLb.clN.sorted.chr.bed
atacseq_consensus_peaks.mLb.clN.bed <- data.frame(fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/consensus_peaks.mLb.clN.sorted.chr.bed", header = F))
colnames(atacseq_consensus_peaks.mLb.clN.bed) <- c("chr", "Start", "End", "interval", "q", "strand")
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_peaks <- merge(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_df, atacseq_consensus_peaks.mLb.clN.bed, by="interval")
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_peaks <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_peaks[,c(9:11,1)]
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_peaks_ext <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_peaks
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_peaks_ext$Start <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_peaks_ext$Start - 1000
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_peaks_ext$End <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_peaks_ext$End + 1000

write.table(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_peaks_ext, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_peaks_ext.bed", sep = "\t", append = F, quote = F, row.names = F, col.names = F)


deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_anno_nearest <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_anno %>% dplyr::filter(Distance < 5000 )
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_anno_nearest <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_anno %>% dplyr::filter(Distance < 5000 )
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_up_anno_nearest <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_up_anno %>% dplyr::filter(Distance < 5000 )
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_down_anno_nearest <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_down_anno %>% dplyr::filter(Distance < 5000 )
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_up_anno_nearest <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_up_anno %>% dplyr::filter(Distance < 5000 )
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_down_anno_nearest <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_down_anno %>% dplyr::filter(Distance < 5000 )



# prepare bed file (requested by Iva)

deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann.bed <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann[,c(9:11,1)]
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann.bed <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann[,c(9:11,1)]

write.table(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann.bed, "deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann.bed", sep = "\t", append = F, quote = F, row.names = F, col.names = F)
write.table(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann.bed, "deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann.bed", sep = "\t", append = F, quote = F, row.names = F, col.names = F)


deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann.UP.bed <- deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann[which(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann$log2FoldChange > 1.5),][,c(9:11,1)]
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann.UP.bed <- deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann[which(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann$log2FoldChange > 1.5),][,c(9:11,1)]

write.table(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann.UP.bed, "deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann.UP.bed", sep = "\t", append = F, quote = F, row.names = F, col.names = F)
write.table(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann.UP.bed, "deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann.UP.bed", sep = "\t", append = F, quote = F, row.names = F, col.names = F)




list_atacseq_venn <- list(shCRAMP1_shControl = unique(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_df$interval), shSUZ12_shControl=unique(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_df$interval))
library(ggvenn)
ggvenn(
  list_atacseq_venn, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

list_atacseq_venn.gene <- list(shCRAMP1_shControl = unique(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann_nearest$Gene.Name), shSUZ12_shControl=unique(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_ann_nearest$Gene.Name))
library(ggvenn)
ggvenn(
  list_atacseq_venn.gene, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)
# add differentially bound intervals annotation
#
# intersect shSUZ12 vs shCRAMP1
atacseqkd_intersected_shSUZ12_shCRAMP1 <- data.frame(interval = intersect(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_df$interval, deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_df$interval))
atacseqkd_intersected_shSUZ12_shCRAMP1_ann <- merge(atacseqkd_intersected_shSUZ12_shCRAMP1, atacseqkd_consensus_peaks.mLb.clN.annotatePeaks, by="interval")

# subtract shSUZ12 from shCRAMP1
atacseqkd_subtracted_shSUZ12_shCRAMP1 <- data.frame(interval = setdiff(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_df$interval, deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_df$interval))
atacseqkd_subtracted_shSUZ12_shCRAMP1_ann <- merge(atacseqkd_subtracted_shSUZ12_shCRAMP1, atacseqkd_consensus_peaks.mLb.clN.annotatePeaks, by="interval")

# subtract shCRAMP1 from shSUZ12
atacseqkd_subtracted_shCRAMP1_shSUZ12 <- data.frame(interval = setdiff(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_fc_df$interval, deseq2_atacseqkd_de_shSUZ12_shControl_0.05_fc_df$interval))
atacseqkd_subtracted_shCRAMP1_shSUZ12_ann <- merge(atacseqkd_subtracted_shCRAMP1_shSUZ12, atacseqkd_consensus_peaks.mLb.clN.annotatePeaks, by="interval")
# atacseqkd_subtracted_shCRAMP1_shSUZ12_ann_nearest <- atacseqkd_subtracted_shCRAMP1_shSUZ12_ann 


# # go term analysis
# gost_atacseqkd_intersected_shSUZ12_shCRAMP1_ann <- gost(query=atacseqkd_intersected_shSUZ12_shCRAMP1_ann$Gene.Name, organism="hsapiens")
# gost_atacseqkd_subtracted_shSUZ12_shCRAMP1_ann <- gost(query=atacseqkd_subtracted_shSUZ12_shCRAMP1_ann$Gene.Name, organism="hsapiens")
# gost_atacseqkd_subtracted_shCRAMP1_shSUZ12_ann <- gost(query=atacseqkd_subtracted_shCRAMP1_shSUZ12_ann$Gene.Name, organism="hsapiens")
# 
# ### Bar plotting with GO:BP terms 
# goterm_analysis_list = list(intersected_shSUZ12_shCRAMP1 = gost_atacseqkd_intersected_shSUZ12_shCRAMP1_ann,  atacseqkd_subtracted_shSUZ12_shCRAMP1 = gost_atacseqkd_subtracted_shSUZ12_shCRAMP1_ann,  atacseqkd_subtracted_shCRAMP1_shSUZ12 = gost_atacseqkd_subtracted_shCRAMP1_shSUZ12_ann)
# 
# goterm_barplot_list <- list()
# gost_gobp_bar_top_list <- list()
# for (names in names(goterm_analysis_list)){
#   print(names)
#   gost.temp4 <- goterm_analysis_list[[names]]
#   gost.temp4.res.sorted <- gost.temp4$result[order(gost.temp4$result$p_value),]
#   gost.temp4.res.sorted["term_name_collapse"] <- paste0(gost.temp4.res.sorted$term_id,"_",gost.temp4.res.sorted$source,"_" ,gsub(" ", ".", gost.temp4.res.sorted$term_name))
#   gost.temp4.res.sorted_gobp <- gost.temp4.res.sorted[which(gost.temp4.res.sorted$source == "GO:BP"),]
#   gost.temp4.res.sorted_gobp_bar <- gost.temp4.res.sorted_gobp[,c(11,3)]
#   gost.temp4.res.sorted_gobp_bar_top <- head(gost.temp4.res.sorted_gobp_bar,10)
#   gost.temp4.res.sorted_gobp_bar_top$p_value <- -log10(gost.temp4.res.sorted_gobp_bar_top$p_value)
#   gost_gobp_bar_top_list[[names]] <- gost.temp4.res.sorted_gobp_bar_top
#   plotBar_temp4 <- ggbarplot(gost.temp4.res.sorted_gobp_bar_top, x = "term_name", y = "p_value", color = "#5d93c4ff", fill = "#5d93c4ff" , sort.by.groups = FALSE,x.text.angle = 90, ylab = "-log10(p.value)", xlab = "Biological Process", legend.title = gsub("gost.","",names), lab.size = 9, sort.val = "asc", rotate = TRUE,  position = position_dodge(),ggtheme = theme_bw())
#   goterm_barplot_list[[names]] <- plotBar_temp4
# }
# 
# pdf(file = "gobp_barlot_plot_list_rnaseq.pdf", height = 10, width = 20)
# ggpubr::ggarrange(plotlist = goterm_barplot_list, ncol=2, nrow=3, common.legend = F, labels=names(goterm_barplot_list), vjust = 1,hjust=-0.5,font.label = list(size = 8, color = "black", face = "bold", family = NULL))
# dev.off()

# atac-peaks annotation

peaks_path = "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/"
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(gprofiler2)
library(ggpubr)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


mypeakfile <- list.files(path=peaks_path,pattern = "*.broadPeak")
mypeakfile_list_gr <- list()
mypeakfile_anno_list <- list()
for (peaksfile in mypeakfile){
  print(peaksfile)
  peaks_temp <- read.table(paste0(peaks_path,peaksfile), header = F)
  peaks_tempgr <- makeGRangesFromDataFrame(peaks_temp,keep.extra.columns=TRUE, seqnames.field=c("V1"),start.field=c("V2"),end.field=c("V3"))
  mypeakfile_list_gr[[peaksfile]] <- peaks_tempgr
  peaks_tempgr_anno <- annotatePeak(peaks_tempgr, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
  mypeakfile_anno_list[[peaksfile]] <- peaks_tempgr_anno
}

# count peaks
macs2_peak.mLb.clN.summary <- fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/qc/macs2_peak.mLb.clN.summary.txt")
macs2_peak.mLb.clN.summary <- data.frame(macs2_peak.mLb.clN.summary[,c(7,9)])
macs2_peak.mLb.clN.summary <- macs2_peak.mLb.clN.summary %>% distinct()
rownames(macs2_peak.mLb.clN.summary) <- macs2_peak.mLb.clN.summary$sample


ggplot(macs2_peak.mLb.clN.summary, aes(x=sample, y=num_peaks, fill=as.factor(sample))) + 
  geom_bar(stat = "identity") + coord_flip()+
  scale_fill_manual(values = c("grey","grey", "blue","blue", "green", "green") )


# pie chart peak distribution
plotAnnoPie(mypeakfile_anno_list$shControl_REP1.mLb.clN_peaks.broadPeak)
plotAnnoPie(mypeakfile_anno_list$shControl_REP2.mLb.clN_peaks.broadPeak)


plotAnnoBar(mypeakfile_anno_list)
plotDistToTSS(mypeakfile_anno_list)
upsetplot(mypeakfile_anno_list)

#average plots
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(mypeakfile_list_gr, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
## binning method 
plotPeakProf2(mypeakfile_list_gr, upstream = 3000, downstream = 3000, conf = 0.95,
              by = "gene", type = "start_site", TxDb = txdb,
              facet = "row", nbin = 800)
tagHeatmap(tagMatrixList)

# cutoffed
mypeak_cutoffedfile <- list.files(path=peaks_path,pattern = "*.broadPeak")
mypeak_cutoffedfile_list_gr <- list()
mypeak_cutoffedfile_anno_list <- list()
for (peak_cutoffedsfile in mypeak_cutoffedfile){
  print(peak_cutoffedsfile)
  peak_cutoffeds_temp <- read.table(paste0(peak_cutoffeds_path,peak_cutoffedsfile), header = F)
  peak_cutoffeds_tempgr <- makeGRangesFromDataFrame(peak_cutoffeds_temp,keep.extra.columns=TRUE, seqnames.field=c("V1"),start.field=c("V2"),end.field=c("V3"))
  mypeak_cutoffedfile_list_gr[[peak_cutoffedsfile]] <- peak_cutoffeds_tempgr
  peak_cutoffeds_tempgr_anno <- annotatepeak_cutoffed(peak_cutoffeds_tempgr, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
  mypeak_cutoffedfile_anno_list[[peak_cutoffedsfile]] <- peak_cutoffeds_tempgr_anno
}

# count peak_cutoffeds
atacseqkd_macs2_peak_cutoffed.mLb.clN.summary <- fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak_cutoffed/qc/macs2_peak_cutoffed.mLb.clN.summary.txt")
atacseqkd_macs2_peak_cutoffed.mLb.clN.summary <- data.frame(atacseqkd_macs2_peak_cutoffed.mLb.clN.summary[,c(7,9)])
atacseqkd_macs2_peak_cutoffed.mLb.clN.summary <- atacseqkd_macs2_peak_cutoffed.mLb.clN.summary %>% distinct()
rownames(atacseqkd_macs2_peak_cutoffed.mLb.clN.summary) <- atacseqkd_macs2_peak_cutoffed.mLb.clN.summary$sample


ggplot(atacseqkd_macs2_peak_cutoffed.mLb.clN.summary, aes(x=sample, y=num_peak_cutoffeds, fill=as.factor(sample))) + 
  geom_bar(stat = "identity") + coord_flip()+
  scale_fill_manual(values = c("grey","grey", "blue","blue", "green", "green") )


# pie chart peak_cutoffed distribution
plotAnnoPie(mypeak_cutoffedfile_anno_list$shControl_REP1.mLb.clN_peak_cutoffeds.broadpeak_cutoffed)
plotAnnoPie(mypeak_cutoffedfile_anno_list$shControl_REP2.mLb.clN_peak_cutoffeds.broadpeak_cutoffed)


plotAnnoBar(mypeak_cutoffedfile_anno_list)
plotDistToTSS(mypeak_cutoffedfile_anno_list)
upsetplot(mypeak_cutoffedfile_anno_list)

# heatmap: tornado plot: performed by deeptools
grep start refTSS_v4.1_human_coordinate.hg38.bed.txt -v | awk '{print $1"\t"$2"\t"$3}' > refTSS_v4.1_human_coordinate.hg38.bed_chr.txt
computeMatrix reference-point -S shControl.mRp.clN.bigWig  shCRAMP1.mRp.clN.bigWig  shSUZ12.mRp.clN.bigWig -R refTSS_v4.1_human_coordinate.hg38.bed_chr.txt -o AroundTSSmatrix --outFileNameMatrix AroundTSSmatrix.txt --outFileSortedRegions AroundTSSmatrix_nosort.txt -a 3000 -b 3000 -bs 25 -p 15 --missingDataAsZero


# motif calling
# meme
bedtools getfasta -fi /mnt/home3/reid/av638/atacseq/iva_lab_gencode/GRCh38.p13.genome.fa -bed shControl.mRp.clN_peaks.broadPeak -fo shControl.mRp.clN_peaks.fa
/opt/meme-5.5.0/bin/meme-chip -meme-nmotifs 5 shControl.mRp.clN_peaks.fa -o  shControl
bedtools getfasta -fi /mnt/home3/reid/av638/atacseq/iva_lab_gencode/GRCh38.p13.genome.fa -bed shCRAMP1.mRp.clN_peaks.broadPeak -fo shCRAMP1.mRp.clN_peaks.fa
/opt/meme-5.5.0/bin/meme-chip -meme-nmotifs 5 shCRAMP1.mRp.clN_peaks.fa -o shCRAMP1
bedtools getfasta -fi /mnt/home3/reid/av638/atacseq/iva_lab_gencode/GRCh38.p13.genome.fa -bed shSUZ12.mRp.clN_peaks.broadPeak -fo shSUZ12.mRp.clN_peaks.fa
/opt/meme-5.5.0/bin/meme-chip -meme-nmotifs 5 shSUZ12.mRp.clN_peaks.fa -o shSUZ12

# ENCODE
H3K27me3_ENCFF801AHF <- data.frame(fread("/mnt/home3/reid/av638/ENCODE/histone/H3K27me3_ENCFF801AHF_peaks_id.bed"))
H3K27me3_ENCFF323WOT <- data.frame(fread("/mnt/home3/reid/av638/ENCODE/histone/H3K27me3_ENCFF323WOT_peaks_id.bed"))
H3K9me3_ENCFF963GZJ <- data.frame(fread("/mnt/home3/reid/av638/ENCODE/histone/H3K9me3_ENCFF963GZJ_peaks_id.bed"))

combined_mypeakfile <- list.files(path=peaks_path,full.names = TRUE, pattern = "*.broadPeak")
combined_mypeakfile <- append(combined_mypeakfile, "/mnt/home3/reid/av638/ENCODE/histone/H3K27me3_ENCFF801AHF_peaks_id.bed")
combined_mypeakfile <- append(combined_mypeakfile, "/mnt/home3/reid/av638/ENCODE/histone/H3K27me3_ENCFF323WOT_peaks_id.bed")
combined_mypeakfile <-append(combined_mypeakfile, "/mnt/home3/reid/av638/ENCODE/histone/H3K9me3_ENCFF963GZJ_peaks_id.bed")


combined_mypeakfile_list_gr <- list()
combined_mypeakfile_anno_list <- list()
for (peaksfile in combined_mypeakfile){
  filename <- unlist(strsplit(peaksfile,"/"))
  filename <- filename[length(filename)]
  print(peaksfile)
  print(filename)
  peaks_temp <- read.table(peaksfile, header = F)
  # print(head(peaks_temp,1))
  peaks_tempgr <- makeGRangesFromDataFrame(peaks_temp,keep.extra.columns=TRUE, seqnames.field=c("V1"),start.field=c("V2"),end.field=c("V3"))
  combined_mypeakfile_list_gr[[filename]] <- peaks_tempgr
  peaks_tempgr_anno <- annotatePeak(peaks_tempgr, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
  combined_mypeakfile_anno_list[[filename]] <- peaks_tempgr_anno
}

# average plots
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
combined_tagMatrixList <- lapply(combined_mypeakfile_list_gr, getTagMatrix, windows=promoter)
plotAvgProf(combined_tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(combined_tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
## binning method 
plotPeakProf2(combined_mypeakfile_list_gr, upstream = 3000, downstream = 3000, conf = 0.95,
              by = "gene", type = "start_site", TxDb = txdb,
              facet = "row", nbin = 800)
tagHeatmap(combined_tagMatrixList)
plotAnnoBar(combined_mypeakfile_anno_list)

# Diffbind for ATAC-Seq
library(DiffBind)
diffbind_atac_seq_samplesheet <- read.csv("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/diffbind_atac_seq_samplesheet.csv")
diffbind_atac_seq_samplesheet
# Taking all the peaks with minOverlap=1
# Creating diffbind object 
atacseqkd_diffbind_obj <- dba(sampleSheet=diffbind_atac_seq_samplesheet,  scoreCol=5, minOverlap=1)
plot(atacseqkd_diffbind_obj)

# Plot of number of peaks that appear in atleast one, two , three so on ... until total samples together.
atacseqkd_diffbind_obj_ovp_rate <- dba.overlap(atacseqkd_diffbind_obj,mode=DBA_OLAP_RATE)
plot(atacseqkd_diffbind_obj_ovp_rate,type ='b',ylab='# peaks', xlab='Overlap at least this many peaksets')
#  Calculate a binding matrix with scores based on read counts,  summits=100 for ATC-Seq # https://www.biostars.org/p/9493721/
atacseqkd_diffbind_obj <- dba.count(atacseqkd_diffbind_obj, minOverlap=1, summits=100, score = DBA_SCORE_READS)
atacseqkd_diffbind_info <- dba.show(atacseqkd_diffbind_obj)
plot(atacseqkd_diffbind_obj)
atacseqkd_diffbind_norm <- dba.normalize(atacseqkd_diffbind_obj, method=DBA_ALL_METHODS)
atacseqkd_diffbind_norm_shCRAMP1_shControl <- dba.contrast(atacseqkd_diffbind_norm, contrast=c("Condition","shCRAMP1","shControl"))
atacseqkd_diffbind_norm_shSUZ12_shControl <- dba.contrast(atacseqkd_diffbind_norm, contrast=c("Condition","shSUZ12","shControl"))


# analyse using both DESeq2 and edgeR
atacseqkd_diffbind_norm_shCRAMP1_shControl <- dba.analyze(atacseqkd_diffbind_norm_shCRAMP1_shControl, method=DBA_ALL_METHODS)
atacseqkd_diffbind_norm_shSUZ12_shControl <- dba.analyze(atacseqkd_diffbind_norm_shSUZ12_shControl, method=DBA_ALL_METHODS)


#in-one step
atacseqkd_diffbind_norm_cont <- dba.contrast(atacseqkd_diffbind_norm, minMembers = 2,categories=DBA_CONDITION)
dba.show(atacseqkd_diffbind_norm_cont, bContrasts=TRUE)
atacseqkd_diffbind_norm_cont <- dba.analyze(atacseqkd_diffbind_norm_cont, method=DBA_ALL_METHODS)
dba.plotPCA(atacseqkd_diffbind_norm_cont, contrast=1, method=DBA_EDGER, attributes=DBA_FACTOR, label=DBA_ID)
atacseqkd_diffbind_norm_cont_shCRAMP1_shControl_edgeR <- dba.report(atacseqkd_diffbind_norm_cont, method=DBA_EDGER, contrast = 1, th=1)
atacseqkd_diffbind_norm_cont_shCRAMP1_shControl_edgeR_df <- data.frame(atacseqkd_diffbind_norm_cont_shCRAMP1_shControl_edgeR)
atacseqkd_diffbind_norm_cont_shCRAMP1_shControl_edgeR_df_de <- atacseqkd_diffbind_norm_cont_shCRAMP1_shControl_edgeR_df %>% dplyr::filter(FDR < 0.05 & (Fold > 1 | Fold < -1)) 

atacseqkd_diffbind_norm_cont_shSUZ12_shControl_edgeR <- dba.report(atacseqkd_diffbind_norm_cont, method=DBA_EDGER, contrast = 2, th=1)
atacseqkd_diffbind_norm_cont_shSUZ12_shControl_edgeR_df <- data.frame(atacseqkd_diffbind_norm_cont_shSUZ12_shControl_edgeR)
atacseqkd_diffbind_norm_cont_shSUZ12_shControl_edgeR_df_de <- atacseqkd_diffbind_norm_cont_shSUZ12_shControl_edgeR_df %>% dplyr::filter(FDR < 0.05 & (Fold > 1 | Fold < -1)) 


# USE THIS CONTRAST way as it is more short and same as individual
# shCRAMP1 summary of results for each tool we can use dba.show
dba.show(atacseqkd_diffbind_norm_shCRAMP1_shControl, bContrasts=T)	

dba.plotPCA(atacseqkd_diffbind_norm_shCRAMP1_shControl, contrast=1, method=DBA_EDGER, attributes=DBA_FACTOR, label=DBA_ID)

# Venn
dba.plotVenn(atacseqkd_diffbind_norm_shCRAMP1_shControl,contrast=1,method=DBA_ALL_METHODS)

# MA-plot
dba.plotMA(atacseqkd_diffbind_norm_shCRAMP1_shControl, method=DBA_DESEQ2)
dba.plotMA(atacseqkd_diffbind_norm_shCRAMP1_shControl, method=DBA_EDGER)

# Volcano plots
dba.plotVolcano(atacseqkd_diffbind_norm_shCRAMP1_shControl)

# Profile plots
dba.plotProfile(atacseqkd_diffbind_norm_shCRAMP1_shControl)

# shSUZ12 summary of results for each tool we can use dba.show
dba.show(atacseqkd_diffbind_norm_shSUZ12_shControl, bContrasts=T)	

dba.plotPCA(atacseqkd_diffbind_norm_shSUZ12_shControl, contrast=1, method=DBA_EDGER, attributes=DBA_FACTOR, label=DBA_ID)

# Venn
dba.plotVenn(atacseqkd_diffbind_norm_shSUZ12_shControl,contrast=1,method=DBA_ALL_METHODS)

# MA-plot
dba.plotMA(atacseqkd_diffbind_norm_shSUZ12_shControl, method=DBA_DESEQ2)
dba.plotMA(atacseqkd_diffbind_norm_shSUZ12_shControl, method=DBA_EDGER)

# Volcano plots
dba.plotVolcano(atacseqkd_diffbind_norm_shSUZ12_shControl, method=DBA_DESEQ2)
dba.plotVolcano(atacseqkd_diffbind_norm_shSUZ12_shControl, method=DBA_EDGER)

# Profile plots
dba.plotProfile(atacseqkd_diffbind_norm_shSUZ12_shControl)

# Extracting results, contrast 1 means , th means threshold for FDR, which if 1 give all sites 
atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR <- dba.report(atacseqkd_diffbind_norm_shCRAMP1_shControl, method=DBA_EDGER, contrast = 1, th=1)
atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df <- data.frame(atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR)
atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df["interval"] <- paste0("interval_",c(1:dim(atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df)[1]))
write.table(atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df, file="atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df.txt", sep="\t", append=F, quote=F, row.names=F, col.names = F)

# Create bed files for each keeping only significant peaks (p < 0.05)
atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df_de <- atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df %>% dplyr::filter(FDR < 0.05 & (Fold > 1 | Fold < -1)) 
atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df_up <- atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df %>% dplyr::filter(FDR < 0.05 & Fold > 1) 
atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df_down <- atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df %>% dplyr::filter(FDR < 0.05 & Fold < -1) 

write.table(atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df_up, file="atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df_up.txt", sep="\t", append=F, quote=F, row.names=F)
write.table(atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df_down, file="atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df_up.txt", sep="\t", append=F, quote=F, row.names=F)

# annotate
system("sort -k1,1 -k2,2n /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df.txt | grep chr > /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_sorted.txt")
atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_anno <- data.frame(fread("bedtools closest -a /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_sorted.txt -b /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/gene_gencode_human_gencode_out.sorted.chr.txt -d"))
colnames(atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_anno) <- c(colnames(atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_df), "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_type", "ensID", "gene", "Distance")
atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_anno_nearest <- atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_anno %>% dplyr::filter(Distance < 5000 )

atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_anno_nearest_de <- atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_anno_nearest %>% dplyr::filter(FDR < 0.05 & (Fold > 1 | Fold < -1)) 
atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_anno_nearest_up <- atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_anno_nearest %>% dplyr::filter(FDR < 0.05 & Fold > 1) 
atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_anno_nearest_down <- atacseqkd_diffbind_norm_shCRAMP1_shControl_edgeR_anno_nearest %>% dplyr::filter(FDR < 0.05 & Fold < -1) 

atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR <- dba.report(atacseqkd_diffbind_norm_shSUZ12_shControl, method=DBA_EDGER, contrast = 1, th=1)
atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df <- data.frame(atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR)
atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df["interval"] <- paste0("interval_",c(1:dim(atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df)[1]))
write.table(atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df, file="atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df.txt", sep="\t", append=F, quote=F, row.names=F, col.names = F)

# Create bed files for each keeping only significant peaks (p < 0.05)
atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df_de <- atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df %>% dplyr::filter(FDR < 0.05 & (Fold > 1 | Fold < -1)) 
atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df_up <- atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df %>% dplyr::filter(FDR < 0.05 & Fold > 1) 
atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df_down <- atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df %>% dplyr::filter(FDR < 0.05 & Fold < -1) 

write.table(atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df_up, file="atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df_up.txt", sep="\t", append=F, quote=F, row.names=F)
write.table(atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df_down, file="atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df_up.txt", sep="\t", append=F, quote=F, row.names=F)

# annotate
system("sort -k1,1 -k2,2n /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df.txt | grep chr > /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_sorted.txt")
atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_anno <- data.frame(fread("bedtools closest -a /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_sorted.txt -b /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/gene_gencode_human_gencode_out.sorted.chr.txt -d"))
colnames(atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_anno) <- c(colnames(atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_df), "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_type", "ensID", "gene", "Distance")
atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_anno_nearest <- atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_anno %>% dplyr::filter(Distance < 5000 )

atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_anno_nearest_de <- atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_anno_nearest %>% dplyr::filter(FDR < 0.05 & (Fold > 1 | Fold < -1)) 
atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_anno_nearest_up <- atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_anno_nearest %>% dplyr::filter(FDR < 0.05 & Fold > 1) 
atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_anno_nearest_down <- atacseqkd_diffbind_norm_shSUZ12_shControl_edgeR_anno_nearest %>% dplyr::filter(FDR < 0.05 & Fold < -1) 


#---------------------------- ATAC-Seq Analysis after annotation to genes --------------------------------#

# # get countdata for homer annotated (genic)
# atacseqkd_genic_precountdata_int <- atacseqkd_precountdata
# 
# atacseqkd_genic_precountdata_int["interval"] <- rownames(atacseqkd_genic_precountdata_int)
# 
# atacseqkd_genic_precountdata_int_ann <- merge(atacseqkd_genic_precountdata_int, atacseqkd_consensus_peaks.mLb.clN.annotatePeaks, by="interval")
# atacseqkd_genic_precountdata_int_ann_nearest <- atacseqkd_genic_precountdata_int_ann 
# 
# atacseqkd_genic_geneagg_countdata <- aggregate(atacseqkd_genic_precountdata_int_ann_nearest[,c(2:7)],by=list(atacseqkd_genic_precountdata_int_ann_nearest$Gene.Name), sum)
# rownames(atacseqkd_genic_geneagg_countdata) <- atacseqkd_genic_geneagg_countdata$Group.1 
# atacseqkd_genic_geneagg_countdata <- atacseqkd_genic_geneagg_countdata[,-1]
# 
# # create sample-specific coldata object
# atacseqkd_genic_dds_coldata <- atacseqkd_dds_coldata
# 
# # prepare count matrix 
# atacseqkd_genic_countdata = as.matrix(atacseqkd_genic_geneagg_countdata) 
# all(rownames(atacseqkd_genic_dds_coldata) == colnames(atacseqkd_genic_countdata)) #should print TRUE
# atacseqkd_genic_ddsO <- DESeqDataSetFromMatrix(countData = atacseqkd_genic_countdata, colData = atacseqkd_genic_dds_coldata, design = ~ condition)
# atacseqkd_genic_keep <- rowSums(counts(atacseqkd_genic_ddsO)) > 10
# atacseqkd_genic_ddsO <- atacseqkd_genic_ddsO[atacseqkd_genic_keep,]
# 
# # pca 
# atacseqkd_genic_vst0Norm <- DESeq2::vst(atacseqkd_genic_ddsO)
# atacseqkd_genic_vst0Norm_pca <- plotPCA(atacseqkd_genic_vst0Norm, intgroup="condition", returnData=TRUE)
# plotPCA(atacseqkd_genic_vst0Norm, intgroup="condition", returnData=FALSE)
# 
# # de analysis
# atacseqkd_genic_ddsO <- DESeq(atacseqkd_genic_ddsO)
# 
# 
# # shCRAMP1 vs shControl
# atacseqkd_genic_de_shCRAMP1_shControl <- results(atacseqkd_genic_ddsO, contrast=c("condition", "shCRAMP1", "shControl"))
# atacseqkd_genic_de_shCRAMP1_shControl = atacseqkd_genic_de_shCRAMP1_shControl[order(rownames(atacseqkd_genic_de_shCRAMP1_shControl)),]
# atacseqkd_genic_de_shCRAMP1_shControl$threshold <- as.logical(atacseqkd_genic_de_shCRAMP1_shControl$padj < 0.05)
# deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05 <- atacseqkd_genic_de_shCRAMP1_shControl[which(atacseqkd_genic_de_shCRAMP1_shControl$threshold == TRUE),]
# deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_df <- data.frame(deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05)
# deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_df["Gene"] <- rownames(deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_df)
# deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_df <- deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_df[,c(8,1:7)]
# deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_fc_df <- deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5) | (log2FoldChange < -1.5))
# deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_fc_up <- deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5))
# deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_fc_down <- deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange < -1.5))
# 
# 
# # log2(1.5) = 2.83
# ggmaplot(atacseqkd_genic_de_shCRAMP1_shControl,
#          fdr = 0.05, fc = 2.83, size = 0.3,
#          palette = c("#D22B2B", "#1465AC", "darkgray"),
#          genenames = rownames(atacseqkd_genic_de_shCRAMP1_shControl),
#          legend = "top", top = 20,
#          font.label = c("bold", 5),
#          font.legend = "bold",
#          font.main = "bold",
#          ggtheme = ggplot2::theme_minimal())
# 
# # shSUZ12 vs shControl
# atacseqkd_genic_de_shSUZ12_shControl <- results(atacseqkd_genic_ddsO, contrast=c("condition", "shSUZ12", "shControl"))
# atacseqkd_genic_de_shSUZ12_shControl = atacseqkd_genic_de_shSUZ12_shControl[order(rownames(atacseqkd_genic_de_shSUZ12_shControl)),]
# atacseqkd_genic_de_shSUZ12_shControl$threshold <- as.logical(atacseqkd_genic_de_shSUZ12_shControl$padj < 0.05)
# deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05 <- atacseqkd_genic_de_shSUZ12_shControl[which(atacseqkd_genic_de_shSUZ12_shControl$threshold == TRUE),]
# deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05_df <- data.frame(deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05)
# deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05_df["Gene"] <- rownames(deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05_df)
# deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05_df <- deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05_df[,c(8,1:7)]
# deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05_fc_df <- deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5) | (log2FoldChange < -1.5))
# deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05_fc_up <- deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5))
# deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05_fc_down <- deseq2_atacseqkd_genic_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange < -1.5))
# 
# 
# ggmaplot(atacseqkd_genic_de_shSUZ12_shControl,
#          fdr = 0.05, fc = 2.83, size = 0.3,
#          palette = c("#D22B2B", "#1465AC", "darkgray"),
#          genenames = rownames(atacseqkd_genic_de_shCRAMP1_shControl),
#          legend = "top", top = 20,
#          font.label = c("bold", 5),
#          font.legend = "bold",
#          font.main = "bold",
#          ggtheme = ggplot2::theme_minimal())
# 
# 
# 
# 
# 
# # get countdata for manually annotated gene (mg)
# atacseqkd_mg_precountdata_int <- atacseqkd_precountdata
# 
# atacseqkd_mg_precountdata_int["interval"] <- rownames(atacseqkd_mg_precountdata_int)
# 
# atacseqkd_mg_precountdata_int_ann <- merge(atacseqkd_mg_precountdata_int, atacseqkd_consensus_peaks.gene, by="interval")
# atacseqkd_mg_precountdata_int_ann_nearest <- atacseqkd_mg_precountdata_int_ann %>% dplyr::filter(Distance < 5000)
# 
# atacseqkd_mg_geneagg_countdata <- aggregate(atacseqkd_mg_precountdata_int_ann_nearest[,c(2:7)],by=list(atacseqkd_mg_precountdata_int_ann_nearest$Gene), sum)
# rownames(atacseqkd_mg_geneagg_countdata) <- atacseqkd_mg_geneagg_countdata$Group.1 
# atacseqkd_mg_geneagg_countdata <- atacseqkd_mg_geneagg_countdata[,-1]
# 
# # create sample-specific coldata object
# atacseqkd_mg_dds_coldata <- atacseqkd_dds_coldata
# 
# # prepare count matrix 
# atacseqkd_mg_countdata = as.matrix(atacseqkd_mg_geneagg_countdata) 
# all(rownames(atacseqkd_mg_dds_coldata) == colnames(atacseqkd_mg_countdata)) #should print TRUE
# atacseqkd_mg_ddsO <- DESeqDataSetFromMatrix(countData = atacseqkd_mg_countdata, colData = atacseqkd_mg_dds_coldata, design = ~ condition)
# atacseqkd_mg_keep <- rowSums(counts(atacseqkd_mg_ddsO)) > 10
# atacseqkd_mg_ddsO <- atacseqkd_mg_ddsO[atacseqkd_mg_keep,]
# 
# # pca 
# atacseqkd_mg_vst0Norm <- DESeq2::vst(atacseqkd_mg_ddsO)
# atacseqkd_mg_vst0Norm_pca <- plotPCA(atacseqkd_mg_vst0Norm, intgroup="condition", returnData=TRUE)
# plotPCA(atacseqkd_mg_vst0Norm, intgroup="condition", returnData=FALSE)
# 
# # de analysis
# atacseqkd_mg_ddsO <- DESeq(atacseqkd_mg_ddsO)
# 
# 
# # shCRAMP1 vs shControl
# atacseqkd_mg_de_shCRAMP1_shControl <- results(atacseqkd_mg_ddsO, contrast=c("condition", "shCRAMP1", "shControl"))
# atacseqkd_mg_de_shCRAMP1_shControl = atacseqkd_mg_de_shCRAMP1_shControl[order(rownames(atacseqkd_mg_de_shCRAMP1_shControl)),]
# atacseqkd_mg_de_shCRAMP1_shControl$threshold <- as.logical(atacseqkd_mg_de_shCRAMP1_shControl$padj < 0.05)
# deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05 <- atacseqkd_mg_de_shCRAMP1_shControl[which(atacseqkd_mg_de_shCRAMP1_shControl$threshold == TRUE),]
# deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_df <- data.frame(deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05)
# deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_df["Gene"] <- rownames(deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_df)
# deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_df <- deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_df[,c(8,1:7)]
# deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_fc_df <- deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5) | (log2FoldChange < -1.5))
# deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_fc_up <- deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5))
# deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_fc_down <- deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange < -1.5))
# 
# 
# # log2(1.5) = 2.83
# ggmaplot(atacseqkd_mg_de_shCRAMP1_shControl,
#          fdr = 0.05, fc = 2.83, size = 0.3,
#          palette = c("#D22B2B", "#1465AC", "darkgray"),
#          genenames = rownames(atacseqkd_mg_de_shCRAMP1_shControl),
#          legend = "top", top = 20,
#          font.label = c("bold", 5),
#          font.legend = "bold",
#          font.main = "bold",
#          ggtheme = ggplot2::theme_minimal())
# 
# # shSUZ12 vs shControl
# atacseqkd_mg_de_shSUZ12_shControl <- results(atacseqkd_mg_ddsO, contrast=c("condition", "shSUZ12", "shControl"))
# atacseqkd_mg_de_shSUZ12_shControl = atacseqkd_mg_de_shSUZ12_shControl[order(rownames(atacseqkd_mg_de_shSUZ12_shControl)),]
# atacseqkd_mg_de_shSUZ12_shControl$threshold <- as.logical(atacseqkd_mg_de_shSUZ12_shControl$padj < 0.05)
# deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05 <- atacseqkd_mg_de_shSUZ12_shControl[which(atacseqkd_mg_de_shSUZ12_shControl$threshold == TRUE),]
# deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05_df <- data.frame(deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05)
# deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05_df["Gene"] <- rownames(deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05_df)
# deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05_df <- deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05_df[,c(8,1:7)]
# deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05_fc_df <- deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5) | (log2FoldChange < -1.5))
# deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05_fc_up <- deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1.5))
# deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05_fc_down <- deseq2_atacseqkd_mg_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange < -1.5))
# 
# 
# ggmaplot(atacseqkd_mg_de_shSUZ12_shControl,
#          fdr = 0.05, fc = 2.83, size = 0.3,
#          palette = c("#D22B2B", "#1465AC", "darkgray"),
#          genenames = rownames(atacseqkd_mg_de_shCRAMP1_shControl),
#          legend = "top", top = 20,
#          font.label = c("bold", 5),
#          font.legend = "bold",
#          font.main = "bold",
#          ggtheme = ggplot2::theme_minimal())
# 
# 
# 


#####============================================================================================================================#####

# cutnrun KO
load("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/deseq2/H3K27me3.consensus_peaks.dds.RData")
cutnrunKOdds <- dds
rm(dds)
cutnrunKOprecountdata <- counts(cutnrunKOdds, normalized=FALSE)
head(cutnrunKOprecountdata)
# rearrange
cutnrunKOprecountdata <- cutnrunKOprecountdata[,c(3,2,4,1)]
dim(cutnrunKOprecountdata)
colSums(cutnrunKOprecountdata)
cutnrunKOprecountdata <- data.frame(cutnrunKOprecountdata)
write.table(cutnrunKOprecountdata, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/cutnrunKOprecountdata.txt", quote = F, append = F, sep="\t")


# NOISeq based analysis
# create sample-specific coldata object
library(NOISeq)
#NOISeq-real
cutnrunKOdds_factors <- data.frame(colData(cutnrunKOdds))

cutnrunKOdds_factors["condition"] <- gsub('[0-9]+', '', unlist(lapply(strsplit(cutnrunKOdds_factors$sample, "_"), function(x) x[2])))
cutnrunKOdds_factors["replicate"] <- c(3,1,1,2)


# rearrange
cutnrunKOdds_factors <- cutnrunKOdds_factors[c(3,2,4,1),c(1,3,4,2)]

cutnrunKOdds_factors$condition <- factor(cutnrunKOdds_factors$condition)
cutnrunKOdds_factors$replicate <- factor(cutnrunKOdds_factors$replicate)


# annotate 
cutnrunKO_consensus_peaks_annotatePeaks <- data.frame(fread("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/H3K27me3.consensus_peaks.annotatePeaks.txt", header = T))
colnames(cutnrunKO_consensus_peaks_annotatePeaks) <- c("interval",colnames(cutnrunKO_consensus_peaks_annotatePeaks)[2:length(cutnrunKO_consensus_peaks_annotatePeaks)])
rownames(cutnrunKO_consensus_peaks_annotatePeaks) <- cutnrunKO_consensus_peaks_annotatePeaks$interval
cutnrunKO_consensus_peaks_annotatePeaks["interval_length"] <- cutnrunKO_consensus_peaks_annotatePeaks$End - cutnrunKO_consensus_peaks_annotatePeaks$Start
cutnrunKOprecountdata["interval"] <- rownames(cutnrunKOprecountdata)
cutnrunKOprecountdata_ann <- merge(cutnrunKOprecountdata, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
rownames(cutnrunKOprecountdata_ann) <- cutnrunKOprecountdata_ann$interval

# filter 
cutnrunKOprecountdata_ann_filt <- cutnrunKOprecountdata_ann[which(rowSums(cutnrunKOprecountdata_ann[,c(2:5)]) > 10),]

cutnrunKOcountdata <- cutnrunKOprecountdata_ann_filt[,c(2,3,4,5)]


cutnrunKO_length <-  data.frame(length=cutnrunKOprecountdata_ann_filt$interval_length)
rownames(cutnrunKO_length) <- rownames(cutnrunKOprecountdata_ann_filt)

cutnrunKO_biotype <-  data.frame(length=cutnrunKOprecountdata_ann_filt$Gene.Type)
rownames(cutnrunKO_biotype) <- rownames(cutnrunKOprecountdata_ann_filt)

#make noiseq object
cutnrunKO_nsobj <- readData(data = (cutnrunKOcountdata),  length = cutnrunKO_length$length, biotype = cutnrunKO_biotype, factors = cutnrunKOdds_factors)
cutnrunKO_myPCA = dat(cutnrunKO_nsobj, type = "PCA")
explo.plot(cutnrunKO_myPCA, factor = "condition")

#get RPKM
cutnrunKO_myRPKM = rpkm(assayData(cutnrunKO_nsobj)$exprs, long = cutnrunKO_length$length, k = 0, lc = 1)
cutnrunKO_myRPKM <- data.frame(cutnrunKO_myRPKM)
cutnrunKO_myRPKM["interval"] <- rownames(cutnrunKO_myRPKM)
cutnrunKO_myRPKM_ann <- merge(cutnrunKO_myRPKM, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")



# #NOIseq-real: Two replicates in one of the experimental conditions are enough to run the algorithm
# #NOISeq BIO cannot be run becuase it two replicates in both condition
# cutnrunKO_myRPKM_noiseq <- noiseq(cutnrunKO_nsobj, k = 0.5, norm = "rpkm", factor = "condition", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "technical")
# 
# head(cutnrunKO_myRPKM_noiseq@results[[1]])
# dim(cutnrunKO_myRPKM_noiseq@results[[1]])
# 
# ### select the differentially expressed features
# cutnrunKO_myRPKM_noiseq.deg = degenes(cutnrunKO_myRPKM_noiseq, q = 0.8, M = NULL)
# cutnrunKO_myRPKM_noiseq.deg.up = degenes(cutnrunKO_myRPKM_noiseq, q = 0.8, M = "up")
# cutnrunKO_myRPKM_noiseq.deg.down = degenes(cutnrunKO_myRPKM_noiseq, q = 0.8, M = "down")
# 
# # Expression plot
# DE.plot(cutnrunKO_myRPKM_noiseq, q = 0.8, graphic = "expr", log.scale = TRUE)
# 
# # MD plot
# DE.plot(cutnrunKO_myRPKM_noiseq, q = 0.8, graphic = "MD")
# 
# 
# #annotate
# cutnrunKO_myRPKM_noiseq.deg.up["interval"] <- rownames(cutnrunKO_myRPKM_noiseq.deg.up)
# cutnrunKO_myRPKM_noiseq.deg.up_ann <- merge(cutnrunKO_myRPKM_noiseq.deg.up, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
# 
# cutnrunKO_myRPKM_noiseq.deg.down["interval"] <- rownames(cutnrunKO_myRPKM_noiseq.deg.down)
# cutnrunKO_myRPKM_noiseq.deg.down_ann <- merge(cutnrunKO_myRPKM_noiseq.deg.down, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
# 
# 
# # NOISeq-slim can be run if I treat each individual replicates as separate
# cutnrunKOdds_factors_slim <- data.frame(colData(cutnrunKOdds))
# 
# cutnrunKOdds_factors_slim["condition"] <- unlist(lapply(strsplit(cutnrunKOdds_factors_slim$sample, "_"), function(x) x[2]))
# 
# 
# # rearrange
# cutnrunKOdds_factors_slim <- cutnrunKOdds_factors_slim[c(3,2,4,1),c(1,3,2)]
# 
# cutnrunKOdds_factors_slim$condition <- factor(cutnrunKOdds_factors_slim$condition)
# 
# #make noiseq object
# cutnrunKO_slim_nsobj <- readData(data = (cutnrunKOcountdata),  length = cutnrunKO_length$length, biotype = cutnrunKO_biotype, factors = cutnrunKOdds_factors_slim)
# cutnrunKO_slim_myPCA = dat(cutnrunKO_slim_nsobj, type = "PCA")
# explo.plot(cutnrunKO_slim_myPCA, factor = "condition")
# 
# 
# cutnrunKO_slim_KO1_WT_myRPKM_noiseq <- noiseq(cutnrunKO_slim_nsobj, factor = "condition", conditions = c("KO1", "WT"), k = NULL, norm = "rpkm", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "no")
# cutnrunKO_slim_KO2_WT_myRPKM_noiseq <- noiseq(cutnrunKO_slim_nsobj, factor = "condition", conditions = c("KO2", "WT"), k = NULL, norm = "rpkm", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "no")
# cutnrunKO_slim_KO3_WT_myRPKM_noiseq <- noiseq(cutnrunKO_slim_nsobj, factor = "condition", conditions = c("KO3", "WT"), k = NULL, norm = "rpkm", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "no")
# 
# head(cutnrunKO_slim_KO1_WT_myRPKM_noiseq@results[[1]])
# dim(cutnrunKO_slim_KO1_WT_myRPKM_noiseq@results[[1]])
# 
# head(cutnrunKO_slim_KO2_WT_myRPKM_noiseq@results[[1]])
# dim(cutnrunKO_slim_KO2_WT_myRPKM_noiseq@results[[1]])
# 
# head(cutnrunKO_slim_KO3_WT_myRPKM_noiseq@results[[1]])
# dim(cutnrunKO_slim_KO3_WT_myRPKM_noiseq@results[[1]])
# 
# 
# ### select the differentially expressed features
# cutnrunKO_slim_KO1_WT_myRPKM_noiseq.deg = degenes(cutnrunKO_slim_KO1_WT_myRPKM_noiseq, q = 0.9, M = NULL)
# cutnrunKO_slim_KO1_WT_myRPKM_noiseq.deg.up = degenes(cutnrunKO_slim_KO1_WT_myRPKM_noiseq, q = 0.9, M = "up")
# cutnrunKO_slim_KO1_WT_myRPKM_noiseq.deg.down = degenes(cutnrunKO_slim_KO1_WT_myRPKM_noiseq, q = 0.9, M = "down")
# 
# 
# cutnrunKO_slim_KO2_WT_myRPKM_noiseq.deg = degenes(cutnrunKO_slim_KO2_WT_myRPKM_noiseq, q = 0.9, M = NULL)
# cutnrunKO_slim_KO2_WT_myRPKM_noiseq.deg.up = degenes(cutnrunKO_slim_KO2_WT_myRPKM_noiseq, q = 0.9, M = "up")
# cutnrunKO_slim_KO2_WT_myRPKM_noiseq.deg.down = degenes(cutnrunKO_slim_KO2_WT_myRPKM_noiseq, q = 0.9, M = "down")
# 
# 
# cutnrunKO_slim_KO3_WT_myRPKM_noiseq.deg = degenes(cutnrunKO_slim_KO3_WT_myRPKM_noiseq, q = 0.9, M = NULL)
# cutnrunKO_slim_KO3_WT_myRPKM_noiseq.deg.up = degenes(cutnrunKO_slim_KO3_WT_myRPKM_noiseq, q = 0.9, M = "up")
# cutnrunKO_slim_KO3_WT_myRPKM_noiseq.deg.down = degenes(cutnrunKO_slim_KO3_WT_myRPKM_noiseq, q = 0.9, M = "down")
# 
# 
# #annotate
# cutnrunKO_slim_KO1_WT_myRPKM_noiseq.deg.up["interval"] <- rownames(cutnrunKO_slim_KO1_WT_myRPKM_noiseq.deg.up)
# cutnrunKO_slim_KO1_WT_myRPKM_noiseq.deg.up_ann <- merge(cutnrunKO_slim_KO1_WT_myRPKM_noiseq.deg.up, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
# 
# cutnrunKO_slim_KO1_WT_myRPKM_noiseq.deg.down["interval"] <- rownames(cutnrunKO_slim_KO1_WT_myRPKM_noiseq.deg.down)
# cutnrunKO_slim_KO1_WT_myRPKM_noiseq.deg.down_ann <- merge(cutnrunKO_slim_KO1_WT_myRPKM_noiseq.deg.down, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
# 
# 
# #annotate
# cutnrunKO_slim_KO2_WT_myRPKM_noiseq.deg.up["interval"] <- rownames(cutnrunKO_slim_KO2_WT_myRPKM_noiseq.deg.up)
# cutnrunKO_slim_KO2_WT_myRPKM_noiseq.deg.up_ann <- merge(cutnrunKO_slim_KO2_WT_myRPKM_noiseq.deg.up, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
# 
# cutnrunKO_slim_KO2_WT_myRPKM_noiseq.deg.down["interval"] <- rownames(cutnrunKO_slim_KO2_WT_myRPKM_noiseq.deg.down)
# cutnrunKO_slim_KO2_WT_myRPKM_noiseq.deg.down_ann <- merge(cutnrunKO_slim_KO2_WT_myRPKM_noiseq.deg.down, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
# 
# 
# #annotate
# cutnrunKO_slim_KO3_WT_myRPKM_noiseq.deg.up["interval"] <- rownames(cutnrunKO_slim_KO3_WT_myRPKM_noiseq.deg.up)
# cutnrunKO_slim_KO3_WT_myRPKM_noiseq.deg.up_ann <- merge(cutnrunKO_slim_KO3_WT_myRPKM_noiseq.deg.up, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
# 
# cutnrunKO_slim_KO3_WT_myRPKM_noiseq.deg.down["interval"] <- rownames(cutnrunKO_slim_KO3_WT_myRPKM_noiseq.deg.down)
# cutnrunKO_slim_KO3_WT_myRPKM_noiseq.deg.down_ann <- merge(cutnrunKO_slim_KO3_WT_myRPKM_noiseq.deg.down, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
# 
# 
# 
# # MD plot
# DE.plot(cutnrunKO_slim_KO1_WT_myRPKM_noiseq, q = 0.9, graphic = "MD")
# # MD plot
# DE.plot(cutnrunKO_slim_KO2_WT_myRPKM_noiseq, q = 0.9, graphic = "MD")
# # MD plot
# DE.plot(cutnrunKO_slim_KO3_WT_myRPKM_noiseq, q = 0.9, graphic = "MD")


# DESEq2 based analysis
# Use DESeq2 and KO1=KO2=KO3=KO
# prepare count matrix
# create sample-specific coldata object
cutnrunKOdds_coldata <- data.frame(colData(cutnrunKOdds))
cutnrunKOdds_coldata["condition"] <- gsub('[0-9]+', '', unlist(lapply(strsplit(cutnrunKOdds_coldata$sample, "_"), function(x) x[2])))
# cutnrunKOdds_coldata["condition"] <- unlist(lapply(strsplit(cutnrunKOdds_coldata$sample, "_"), function(x) x[2]))

cutnrunKOdds_coldata["replicate"] <- c(3,1,1,2)

# rearrange
cutnrunKOdds_coldata <- cutnrunKOdds_coldata[c(3,2,4,1),c(1,3,4,2)]


cutnrunKOdds_coldata$condition <- factor(cutnrunKOdds_coldata$condition)
cutnrunKOdds_coldata$replicate <- factor(cutnrunKOdds_coldata$replicate)

all(rownames(cutnrunKOdds_coldata) == colnames(cutnrunKOcountdata)) #should print TRUE
cutnrunKOddsO <- DESeqDataSetFromMatrix(countData = cutnrunKOcountdata, colData = cutnrunKOdds_coldata, design = ~condition)
cutnrunKOkeep <- rowSums(counts(cutnrunKOddsO)) > 10
cutnrunKOddsO <- cutnrunKOddsO[cutnrunKOkeep,]

# pca
cutnrunKOvst0Norm <- DESeq2::vst(cutnrunKOddsO)
cutnrunKOvst0Norm_pca <- plotPCA(cutnrunKOvst0Norm, intgroup="condition", returnData=TRUE)
plotPCA(cutnrunKOvst0Norm, intgroup="condition", returnData=FALSE)


# de analysis
cutnrunKOddsO <- DESeq(cutnrunKOddsO)

# KO vs WT
cutnrunKOde_KO_WT <- results(cutnrunKOddsO, contrast=c("condition", "KO", "WT"))
cutnrunKOde_KO_WT = cutnrunKOde_KO_WT[order(rownames(cutnrunKOde_KO_WT)),]
cutnrunKOde_KO_WT$threshold <- as.logical(cutnrunKOde_KO_WT$padj < 0.05)
deseq2_cutnrunKOde_KO_WT_0.05 <- cutnrunKOde_KO_WT[which(cutnrunKOde_KO_WT$threshold == TRUE),]
deseq2_cutnrunKOde_KO_WT_0.05_df <- data.frame(deseq2_cutnrunKOde_KO_WT_0.05)
deseq2_cutnrunKOde_KO_WT_0.05_df["interval"] <- rownames(deseq2_cutnrunKOde_KO_WT_0.05_df)
deseq2_cutnrunKOde_KO_WT_0.05_df <- deseq2_cutnrunKOde_KO_WT_0.05_df[,c(8,1:7)]
deseq2_cutnrunKOde_KO_WT_0.05_fc_df <- deseq2_cutnrunKOde_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange > 1.5) | (log2FoldChange < -1.5))

cutnrunKOde_KO_WT_df <- data.frame(cutnrunKOde_KO_WT)
cutnrunKOde_KO_WT_df["interval"] <- rownames(cutnrunKOde_KO_WT_df) 

deseq2_cutnrunKOde_KO_WT_0.05_fc.up <- deseq2_cutnrunKOde_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange > 1.5) )
deseq2_cutnrunKOde_KO_WT_0.05_fc.down <- deseq2_cutnrunKOde_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange < -1.5) )




# annotate all cutrun KO peaks data
cutnrunKOde_KO_WT_re <- cutnrunKOde_KO_WT
cutnrunKOde_KO_WT_re["interval"] <- rownames(cutnrunKOde_KO_WT_re)
cutnrunKOde_KO_WT_re_df <- data.frame(cutnrunKOde_KO_WT_re)
cutnrunKOde_KO_WT_re_anno <- merge(cutnrunKOde_KO_WT_re_df, cutnrunKO_consensus_peaks.gene, by="interval")
cutnrunKOde_KO_WT_re_anno_nearest <- cutnrunKOde_KO_WT_re_anno %>% dplyr::filter(Distance < 5000 )


#annotate DE data
deseq2_cutnrunKOde_KO_WT_0.05_fc_df["interval"] <- rownames(deseq2_cutnrunKOde_KO_WT_0.05_fc_df)
deseq2_cutnrunKOde_KO_WT_0.05_fc_ann <- merge(deseq2_cutnrunKOde_KO_WT_0.05_fc_df, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")

deseq2_cutnrunKOde_KO_WT_0.05_fc.up["interval"] <- rownames(deseq2_cutnrunKOde_KO_WT_0.05_fc.up)
deseq2_cutnrunKOde_KO_WT_0.05_fc.up_ann <- merge(deseq2_cutnrunKOde_KO_WT_0.05_fc.up, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")

deseq2_cutnrunKOde_KO_WT_0.05_fc.down["interval"] <- rownames(deseq2_cutnrunKOde_KO_WT_0.05_fc.down)
deseq2_cutnrunKOde_KO_WT_0.05_fc.down_ann <- merge(deseq2_cutnrunKOde_KO_WT_0.05_fc.down, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")


# log2(1.5) = 2.83
ggmaplot(cutnrunKOde_KO_WT_df,
         fdr = 0.05, fc = 2.83, size = 0.3,
         palette = c("#D22B2B", "#1465AC", "darkgray"),
         genenames = cutnrunKOde_KO_WT_df$interval,
         legend = "top", top = 20,
         font.label = c("bold", 5),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())



# # edgeR based analysis
# 
# # create sample-specific coldata object
# cutnrunKO_meta <- data.frame(colData(cutnrunKOdds))
# cutnrunKO_meta["condition"] <- gsub('[0-9]+', '', unlist(lapply(strsplit(cutnrunKO_meta$sample, "_"), function(x) x[2])))
# cutnrunKO_meta <- cutnrunKO_meta[c(3,2,4,1),c(1,3,2)]
# 
# cutnrunKO_group <- factor(cutnrunKO_meta$condition, levels = c("WT", "KO"))
# cutnrunKO_eObj <- DGEList(counts = cutnrunKOcountdata, group = cutnrunKO_group)
# eObjkeep <- filterByExpr(cutnrunKO_eObj, group = cutnrunKO_group)
# cutnrunKO_eObj <- cutnrunKO_eObj[eObjkeep,,keep.lib.sizes=FALSE]
# 
# #normalize
# cutnrunKO_eObj <- calcNormFactors(cutnrunKO_eObj)
# 
# # Perform MDS analysis
# limma::plotMDS(cutnrunKO_eObj, col=c("red", "blue", "blue", "blue"), pch=6,cex = 3)
# 
# #design
# cutnrunKO_design <- model.matrix(~0+cutnrunKO_group, data=cutnrunKO_eObj$samples)
# colnames(cutnrunKO_design) <- levels(cutnrunKO_eObj$samples$group)
# 
# #estimate dispersion
# cutnrunKO_eObj <- estimateDisp(cutnrunKO_eObj, cutnrunKO_design)
# 
# 
# plotBCV(cutnrunKO_eObj)
# 
# cutnrunKO_eObj_et <- exactTest(cutnrunKO_eObj)
# 
# 
# # Extract significant genes
# # dim(rawdata)[[1]] contains the total number of genes (number of lines)
# cutnrunKO_eOBJ_top_intervals <- topTags(cutnrunKO_eObj_et, n=dim(cutnrunKOcountdata)[[1]])
# 
# #How many are significant between the cell lines?
# table(cutnrunKO_eOBJ_top_intervals$table$FDR < 0.05)
# 
# 
# cutnrunKO_eOBJ_sig_intervals <- cutnrunKO_eOBJ_top_intervals$table[which(cutnrunKO_eOBJ_top_intervals$table$FDR < 0.05),]
# cutnrunKO_eOBJ_sig_intervals_all <- cutnrunKO_eOBJ_top_intervals$table[which(cutnrunKO_eOBJ_top_intervals$table$FDR < 1),]
# 
# # Filter for upregulated
# cutnrunKO_eOBJ_sig_intervals.up <- cutnrunKO_eOBJ_sig_intervals[which(cutnrunKO_eOBJ_sig_intervals$logFC > 1.5),]
# 
# # Filter for downregulated genes
# cutnrunKO_eOBJ_sig_intervals.down <- cutnrunKO_eOBJ_sig_intervals[which(cutnrunKO_eOBJ_sig_intervals$logFC < -1.5),]
# 
# plot(cutnrunKO_eOBJ_sig_intervals_all$logCPM, cutnrunKO_eOBJ_sig_intervals_all$logFC, cex=0.2)
# 
# # Make a basic MA plot 
# with(cutnrunKO_eOBJ_sig_intervals_all, plot(logCPM, logFC, pch=21, cex = 0.3, main="MA plot", col="grey", xlim=c(0,12), ylim=c(-10,10)))
# #Coloroffsets
# #Upregulated
# with(subset(cutnrunKO_eOBJ_sig_intervals_all, FDR < 0.05 & logFC > 1.5), points(logCPM, logFC, pch=21, cex = 0.4, xlim=c(0,12), ylim=c(-10,10),lwd = 0.4, col="black",bg="red"))
# #Downregulated
# with(subset(cutnrunKO_eOBJ_sig_intervals_all, FDR < 0.05 & logFC < -1.5), points(logCPM, logFC, pch=21, cex = 0.4, xlim=c(0,12), ylim=c(-10,10), lwd = 0.4,col="black", bg="blue"))
# 
# 
# #annotate
# cutnrunKO_eOBJ_sig_intervals.up["interval"] <- rownames(cutnrunKO_eOBJ_sig_intervals.up)
# cutnrunKO_eOBJ_sig_intervals.up_ann <- merge(cutnrunKO_eOBJ_sig_intervals.up, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
# 
# cutnrunKO_eOBJ_sig_intervals.down["interval"] <- rownames(cutnrunKO_eOBJ_sig_intervals.down)
# cutnrunKO_eOBJ_sig_intervals.down_ann <- merge(cutnrunKO_eOBJ_sig_intervals.down, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
# 
# 
# #NOISeq annotate
# cutnrunKO_myRPKM_noiseq.deg.up_ann
# cutnrunKO_myRPKM_noiseq.deg.down_ann
# 
# #DESEq2 annotate
# deseq2_cutnrunKOde_KO_WT_0.05_fc.up_ann
# deseq2_cutnrunKOde_KO_WT_0.05_fc.down_ann
# 
# #edgeR annotate
# cutnrunKO_eOBJ_sig_intervals.up_ann
# cutnrunKO_eOBJ_sig_intervals.down_ann
# 
# # intervals
# list_cutnrunKO_venn.up <- list(noiseq_cutrunKO_up = (cutnrunKO_myRPKM_noiseq.deg.up_ann$interval),
#                             deseq2_cutrunKO_up = (deseq2_cutnrunKOde_KO_WT_0.05_fc.up_ann$interval),
#                             edgeR_cutrunKO_up = (cutnrunKO_eOBJ_sig_intervals.up_ann$interval)
#                             )
# list_cutnrunKO_venn.down <- list(
#                                noiseq_cutrunKO_down = (cutnrunKO_myRPKM_noiseq.deg.down_ann$interval),
#                                deseq2_cutrunKO_down = (deseq2_cutnrunKOde_KO_WT_0.05_fc.down_ann$interval),
#                                edgeR_cutrunKO_down = (cutnrunKO_eOBJ_sig_intervals.down_ann$interval)
# )
# 
library(ggvenn)
ggvenn(
  list_cutnrunKO_venn.up, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#FF6961"),
  stroke_size = 0.5, set_name_size = 4
)
ggvenn(
  list_cutnrunKO_venn.down, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#FF6961"),
  stroke_size = 0.5, set_name_size = 4
)

# nearest genes
cutnrunKO_myRPKM_noiseq.deg.up_ann_nearest <- cutnrunKO_myRPKM_noiseq.deg.up_ann 
cutnrunKO_myRPKM_noiseq.deg.down_ann_nearest <- cutnrunKO_myRPKM_noiseq.deg.down_ann 
deseq2_cutnrunKOde_KO_WT_0.05_fc.up_ann_nearest <- deseq2_cutnrunKOde_KO_WT_0.05_fc.up_ann 
deseq2_cutnrunKOde_KO_WT_0.05_fc.down_ann_nearest <- deseq2_cutnrunKOde_KO_WT_0.05_fc.down_ann 
cutnrunKO_eOBJ_sig_intervals.up_ann_nearest <- cutnrunKO_eOBJ_sig_intervals.up_ann 
cutnrunKO_eOBJ_sig_intervals.down_ann_nearest <- cutnrunKO_eOBJ_sig_intervals.down_ann 

deseq2_cutnrunKOde_KO_WT_0.05_fc_ann_nearest <- deseq2_cutnrunKOde_KO_WT_0.05_fc_ann 

list_cutnrunKO_venn.up.gene <- list(noiseq_cutrunKO_up = (cutnrunKO_myRPKM_noiseq.deg.up_ann_nearest$Gene.Name),
                                    deseq2_cutrunKO_up = (deseq2_cutnrunKOde_KO_WT_0.05_fc.up_ann_nearest$Gene.Name),
                                    edgeR_cutrunKO_up = (cutnrunKO_eOBJ_sig_intervals.up_ann_nearest$Gene.Name)
)
list_cutnrunKO_venn.down.gene <- list(
  noiseq_cutrunKO_down = (cutnrunKO_myRPKM_noiseq.deg.down_ann_nearest$Gene.Name),
  deseq2_cutrunKO_down = (deseq2_cutnrunKOde_KO_WT_0.05_fc.down_ann_nearest$Gene.Name),
  edgeR_cutrunKO_down = (cutnrunKO_eOBJ_sig_intervals.down_ann_nearest$Gene.Name)
)

library(ggvenn)
ggvenn(
  list_cutnrunKO_venn.up.gene, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#FF6961"),
  stroke_size = 0.5, set_name_size = 4
)
ggvenn(
  list_cutnrunKO_venn.down.gene, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#FF6961"),
  stroke_size = 0.5, set_name_size = 4
)

#Use manual annotation
sort -k1,1 -k2,2n H3K27me3.consensus_peaks.bed| grep chr > H3K27me3.consensus_peaks_chr.bed
sort -k1,1 -k2,2n /mnt/home3/reid/av638/rnaseq/iva_lab_oct23/gene_gencode_human__gencode_out.txt | grep chr > gene_gencode_human_gencode_out.sorted.chr.txt
bedtools closest -a H3K27me3.consensus_peaks_chr.bed -b gene_gencode_human_gencode_out.sorted.chr.txt -d > H3K27me3.consensus_peaks_chr.ann.bed

cutnrunKO_consensus_peaks.gene <- data.frame(fread("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/H3K27me3.consensus_peaks_chr.ann.bed", header = F))
cutnrunKO_consensus_peaks.gene <- cutnrunKO_consensus_peaks.gene[,c(4,12,13,14)]
colnames(cutnrunKO_consensus_peaks.gene) <- c("interval", "ensID","Gene", "Distance")
cutnrunKO_consensus_peaks.gene["ensID_Gene"] <- paste0(cutnrunKO_consensus_peaks.gene$ensID, "%", cutnrunKO_consensus_peaks.gene$Gene)

#annotate
deseq2_cutnrunKOde_KO_WT_0.05_fc_anno <- merge(deseq2_cutnrunKOde_KO_WT_0.05_fc_df, cutnrunKO_consensus_peaks.gene, by="interval")

deseq2_cutnrunKOde_KO_WT_0.05_fc.up_anno <- merge(deseq2_cutnrunKOde_KO_WT_0.05_fc.up, cutnrunKO_consensus_peaks.gene , by="interval")
deseq2_cutnrunKOde_KO_WT_0.05_fc.down_anno <- merge(deseq2_cutnrunKOde_KO_WT_0.05_fc.down, cutnrunKO_consensus_peaks.gene , by="interval")

deseq2_cutnrunKOde_KO_WT_0.05_fc_anno_nearest <- deseq2_cutnrunKOde_KO_WT_0.05_fc_anno %>% dplyr::filter(Distance < 5000 )

deseq2_cutnrunKOde_KO_WT_0.05_fc.up_anno_nearest <- deseq2_cutnrunKOde_KO_WT_0.05_fc.up_anno %>% dplyr::filter(Distance < 5000 )
deseq2_cutnrunKOde_KO_WT_0.05_fc.down_anno_nearest <- deseq2_cutnrunKOde_KO_WT_0.05_fc.down_anno %>% dplyr::filter(Distance < 5000 )



# cutnrun_KO-peaks annotation

peaks_path = "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/"
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(gprofiler2)
library(ggpubr)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


cutnrunKO_mypeakfile <- list.files(path=peaks_path,pattern = "*.broadPeak")
cutnrunKO_mypeakfile_list_gr <- list()
cutnrunKO_mypeakfile_anno_list <- list()
for (peaksfile in cutnrunKO_mypeakfile){
  print(peaksfile)
  peaks_temp <- read.table(paste0(peaks_path,peaksfile), header = F)
  peaks_tempgr <- makeGRangesFromDataFrame(peaks_temp,keep.extra.columns=TRUE, seqnames.field=c("V1"),start.field=c("V2"),end.field=c("V3"))
  cutnrunKO_mypeakfile_list_gr[[peaksfile]] <- peaks_tempgr
  peaks_tempgr_anno <- annotatePeak(peaks_tempgr, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
  cutnrunKO_mypeakfile_anno_list[[peaksfile]] <- peaks_tempgr_anno
}

# count peaks
cutnrunKO_macs2_peak.summary <- fread("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/qc/macs2_peak.summary.txt")
cutnrunKO_macs2_peak.summary <- data.frame(cutnrunKO_macs2_peak.summary[,c(7,9)])
cutnrunKO_macs2_peak.summary <- cutnrunKO_macs2_peak.summary %>% distinct()
rownames(cutnrunKO_macs2_peak.summary) <- cutnrunKO_macs2_peak.summary$sample


ggplot(cutnrunKO_macs2_peak.summary, aes(x=sample, y=num_peaks, fill=as.factor(sample))) + 
  geom_bar(stat = "identity") + coord_flip()+
  scale_fill_manual(values = c("orange","orange", "orange", "navy") ) + theme_bw()


# pie chart peak distribution
plotAnnoPie(cutnrunKO_mypeakfile_anno_list$H3K27me3_WT_peaks.broadPeak)
plotAnnoPie(cutnrunKO_mypeakfile_anno_list$H3K27me3_KO1_peaks.broadPeak)
plotAnnoPie(cutnrunKO_mypeakfile_anno_list$H3K27me3_KO2_peaks.broadPeak)
plotAnnoPie(cutnrunKO_mypeakfile_anno_list$H3K27me3_KO3_peaks.broadPeak)

plotAnnoBar(cutnrunKO_mypeakfile_anno_list)
plotDistToTSS(cutnrunKO_mypeakfile_anno_list)
upsetplot(cutnrunKO_mypeakfile_list_gr)


# Diffbind for ATAC-Seq
library(DiffBind)
diffbind_cutnrunKO_samplesheet <- read.csv("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/diffbind_cutnrunKO_samplesheet.csv")
diffbind_cutnrunKO_samplesheet
# Taking all the peaks with minOverlap=1
# Creating diffbind object 
cutnrunKO_diffbind_obj <- dba(sampleSheet=diffbind_cutnrunKO_samplesheet,  scoreCol=5, minOverlap=1)
plot(cutnrunKO_diffbind_obj)

# Plot of number of peaks that appear in atleast one, two , three so on ... until total samples together.
cutnrunKO_diffbind_obj_ovp_rate <- dba.overlap(cutnrunKO_diffbind_obj,mode=DBA_OLAP_RATE)
plot(cutnrunKO_diffbind_obj_ovp_rate,type ='b',ylab='# peaks', xlab='Overlap at least this many peaksets')
#  Calculate a binding matrix with scores based on read counts,  summits=100 for ATC-Seq # https://www.biostars.org/p/9493721/
cutnrunKO_diffbind_obj <- dba.count(cutnrunKO_diffbind_obj, minOverlap=1, summits=200, score = DBA_SCORE_READS)
cutnrunKO_diffbind_info <- dba.show(cutnrunKO_diffbind_obj)
plot(cutnrunKO_diffbind_obj)
cutnrunKO_diffbind_norm <- dba.normalize(cutnrunKO_diffbind_obj, method=DBA_ALL_METHODS)
cutnrunKO_diffbind_norm_KO_WT <- dba.contrast(cutnrunKO_diffbind_norm, contrast=c("Condition","KO","WT"))


# analyse using both DESeq2 and edgeR
cutnrunKO_diffbind_norm_KO_WT <- dba.analyze(cutnrunKO_diffbind_norm_KO_WT, method=DBA_ALL_METHODS)


#in-one step
cutnrunKO_diffbind_norm_cont <- dba.contrast(cutnrunKO_diffbind_norm, minMembers = 1,categories=DBA_CONDITION)
dba.show(cutnrunKO_diffbind_norm_cont, bContrasts=TRUE)
cutnrunKO_diffbind_norm_cont <- dba.analyze(cutnrunKO_diffbind_norm_cont, method=DBA_ALL_METHODS)
dba.plotPCA(cutnrunKO_diffbind_norm_cont, contrast=1, method=DBA_EDGER, attributes=DBA_FACTOR, label=DBA_ID)
cutnrunKO_diffbind_norm_cont_KO_WT_edgeR <- dba.report(cutnrunKO_diffbind_norm_cont, method=DBA_EDGER, contrast = 1, th=1)
cutnrunKO_diffbind_norm_cont_KO_WT_edgeR_df <- data.frame(cutnrunKO_diffbind_norm_cont_KO_WT_edgeR)
cutnrunKO_diffbind_norm_cont_KO_WT_edgeR_df_de <- cutnrunKO_diffbind_norm_cont_KO_WT_edgeR_df %>% dplyr::filter(FDR < 0.05 & (Fold > 1 | Fold < -1)) 



# USE THIS CONTRAST way as it is more short and same as individual
# KO summary of results for each tool we can use dba.show
dba.show(cutnrunKO_diffbind_norm_KO_WT, bContrasts=T)	

dba.plotPCA(cutnrunKO_diffbind_norm_KO_WT, contrast=1, method=DBA_EDGER, attributes=DBA_FACTOR, label=DBA_ID)

# Venn
dba.plotVenn(cutnrunKO_diffbind_norm_KO_WT,contrast=1,method=DBA_ALL_METHODS)

# MA-plot
dba.plotMA(cutnrunKO_diffbind_norm_KO_WT, method=DBA_DESEQ2)
dba.plotMA(cutnrunKO_diffbind_norm_KO_WT, method=DBA_EDGER)

# Volcano plots
dba.plotVolcano(cutnrunKO_diffbind_norm_KO_WT)

# Profile plots
dba.plotProfile(cutnrunKO_diffbind_norm_KO_WT)

#
# Extracting results, contrast 1 means , th means threshold for FDR, which if 1 give all sites 
cutnrunKO_diffbind_norm_KO_WT_edgeR <- dba.report(cutnrunKO_diffbind_norm_KO_WT, method=DBA_EDGER, contrast = 1, th=1)
cutnrunKO_diffbind_norm_KO_WT_edgeR_df <- data.frame(cutnrunKO_diffbind_norm_KO_WT_edgeR)
cutnrunKO_diffbind_norm_KO_WT_edgeR_df["interval"] <- paste0("interval_",c(1:dim(cutnrunKO_diffbind_norm_KO_WT_edgeR_df)[1]))
write.table(cutnrunKO_diffbind_norm_KO_WT_edgeR_df, file="cutnrunKO_diffbind_norm_KO_WT_edgeR_df.txt", sep="\t", append=F, quote=F, row.names=F, col.names = F)

# Create bed files for each keeping only significant peaks (p < 0.05)
cutnrunKO_diffbind_norm_KO_WT_edgeR_df_de <- cutnrunKO_diffbind_norm_KO_WT_edgeR_df %>% dplyr::filter(FDR < 0.05 & (Fold > 1 | Fold < -1)) 
cutnrunKO_diffbind_norm_KO_WT_edgeR_df_up <- cutnrunKO_diffbind_norm_KO_WT_edgeR_df %>% dplyr::filter(FDR < 0.05 & Fold > 1) 
cutnrunKO_diffbind_norm_KO_WT_edgeR_df_down <- cutnrunKO_diffbind_norm_KO_WT_edgeR_df %>% dplyr::filter(FDR < 0.05 & Fold < -1) 

write.table(cutnrunKO_diffbind_norm_KO_WT_edgeR_df_up, file="cutnrunKO_diffbind_norm_KO_WT_edgeR_df_up.txt", sep="\t", append=F, quote=F, row.names=F)
write.table(cutnrunKO_diffbind_norm_KO_WT_edgeR_df_down, file="cutnrunKO_diffbind_norm_KO_WT_edgeR_df_up.txt", sep="\t", append=F, quote=F, row.names=F)

# annotate
system("sort -k1,1 -k2,2n /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/cutnrunKO_diffbind_norm_KO_WT_edgeR_df.txt | grep chr > /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/cutnrunKO_diffbind_norm_KO_WT_edgeR_sorted.txt")
cutnrunKO_diffbind_norm_KO_WT_edgeR_anno <- data.frame(fread("bedtools closest -a /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/cutnrunKO_diffbind_norm_KO_WT_edgeR_sorted.txt -b /mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/macs2/broad_peak/consensus/gene_gencode_human_gencode_out.sorted.chr.txt -d"))
colnames(cutnrunKO_diffbind_norm_KO_WT_edgeR_anno) <- c(colnames(cutnrunKO_diffbind_norm_KO_WT_edgeR_df), "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_type", "ensID", "gene", "Distance")
cutnrunKO_diffbind_norm_KO_WT_edgeR_anno_nearest <- cutnrunKO_diffbind_norm_KO_WT_edgeR_anno %>% dplyr::filter(Distance < 5000 )

cutnrunKO_diffbind_norm_KO_WT_edgeR_anno_nearest_de <- cutnrunKO_diffbind_norm_KO_WT_edgeR_anno_nearest %>% dplyr::filter(FDR < 0.05 & (Fold > 1 | Fold < -1)) 
cutnrunKO_diffbind_norm_KO_WT_edgeR_anno_nearest_up <- cutnrunKO_diffbind_norm_KO_WT_edgeR_anno_nearest %>% dplyr::filter(FDR < 0.05 & Fold > 1) 
cutnrunKO_diffbind_norm_KO_WT_edgeR_anno_nearest_down <- cutnrunKO_diffbind_norm_KO_WT_edgeR_anno_nearest %>% dplyr::filter(FDR < 0.05 & Fold < -1) 





#---------------------------- CUT&RUN-KO-Seq Analysis after annotation to genes --------------------------------#

# # get countdata for homer annotated (genic)
# cutnrunKO_genic_precountdata_int <- cutnrunKOprecountdata
# 
# cutnrunKO_genic_precountdata_int["interval"] <- rownames(cutnrunKO_genic_precountdata_int)
# 
# cutnrunKO_genic_precountdata_int_ann <- merge(cutnrunKO_genic_precountdata_int, cutnrunKO_consensus_peaks_annotatePeaks, by="interval")
# cutnrunKO_genic_precountdata_int_ann_nearest <- cutnrunKO_genic_precountdata_int_ann 
# 
# cutnrunKO_genic_geneagg_countdata <- aggregate(cutnrunKO_genic_precountdata_int_ann_nearest[,c(2:5)],by=list(cutnrunKO_genic_precountdata_int_ann_nearest$Gene.Name), sum)
# rownames(cutnrunKO_genic_geneagg_countdata) <- cutnrunKO_genic_geneagg_countdata$Group.1 
# cutnrunKO_genic_geneagg_countdata <- cutnrunKO_genic_geneagg_countdata[,-1]
# 
# # create sample-specific coldata object
# cutnrunKO_genic_dds_coldata <- cutnrunKOdds_coldata
# 
# # prepare count matrix 
# cutnrunKO_genic_countdata = as.matrix(cutnrunKO_genic_geneagg_countdata) 
# all(rownames(cutnrunKO_genic_dds_coldata) == colnames(cutnrunKO_genic_countdata)) #should print TRUE
# cutnrunKO_genic_ddsO <- DESeqDataSetFromMatrix(countData = cutnrunKO_genic_countdata, colData = cutnrunKO_genic_dds_coldata, design = ~ condition)
# cutnrunKO_genic_keep <- rowSums(counts(cutnrunKO_genic_ddsO)) > 10
# cutnrunKO_genic_ddsO <- cutnrunKO_genic_ddsO[cutnrunKO_genic_keep,]
# 
# # pca 
# cutnrunKO_genic_vst0Norm <- DESeq2::vst(cutnrunKO_genic_ddsO)
# cutnrunKO_genic_vst0Norm_pca <- plotPCA(cutnrunKO_genic_vst0Norm, intgroup="condition", returnData=TRUE)
# plotPCA(cutnrunKO_genic_vst0Norm, intgroup="condition", returnData=FALSE)
# 
# # de analysis
# cutnrunKO_genic_ddsO <- DESeq(cutnrunKO_genic_ddsO)
# 
# 
# # KO vs WT
# cutnrunKO_genic_de_KO_WT <- results(cutnrunKO_genic_ddsO, contrast=c("condition", "KO", "WT"))
# cutnrunKO_genic_de_KO_WT = cutnrunKO_genic_de_KO_WT[order(rownames(cutnrunKO_genic_de_KO_WT)),]
# cutnrunKO_genic_de_KO_WT$threshold <- as.logical(cutnrunKO_genic_de_KO_WT$padj < 0.05)
# deseq2_cutnrunKO_genic_de_KO_WT_0.05 <- cutnrunKO_genic_de_KO_WT[which(cutnrunKO_genic_de_KO_WT$threshold == TRUE),]
# deseq2_cutnrunKO_genic_de_KO_WT_0.05_df <- data.frame(deseq2_cutnrunKO_genic_de_KO_WT_0.05)
# deseq2_cutnrunKO_genic_de_KO_WT_0.05_df["Gene"] <- rownames(deseq2_cutnrunKO_genic_de_KO_WT_0.05_df)
# deseq2_cutnrunKO_genic_de_KO_WT_0.05_df <- deseq2_cutnrunKO_genic_de_KO_WT_0.05_df[,c(8,1:7)]
# deseq2_cutnrunKO_genic_de_KO_WT_0.05_fc_df <- deseq2_cutnrunKO_genic_de_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange > 1.5) | (log2FoldChange < -1.5))
# deseq2_cutnrunKO_genic_de_KO_WT_0.05_fc_up <- deseq2_cutnrunKO_genic_de_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange > 1.5))
# deseq2_cutnrunKO_genic_de_KO_WT_0.05_fc_down <- deseq2_cutnrunKO_genic_de_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange < -1.5))
# 
# 
# # log2(1.5) = 2.83
# ggmaplot(cutnrunKO_genic_de_KO_WT,
#          fdr = 0.05, fc = 2.83, size = 0.3,
#          palette = c("#D22B2B", "#1465AC", "darkgray"),
#          genenames = rownames(cutnrunKO_genic_de_KO_WT),
#          legend = "top", top = 20,
#          font.label = c("bold", 5),
#          font.legend = "bold",
#          font.main = "bold",
#          ggtheme = ggplot2::theme_minimal())
# 
# 
# 
# 
# # get countdata for manually annotated gene (mg)
# cutnrunKO_mg_precountdata_int <- cutnrunKOprecountdata
# 
# cutnrunKO_mg_precountdata_int["interval"] <- rownames(cutnrunKO_mg_precountdata_int)
# 
# cutnrunKO_mg_precountdata_int_ann <- merge(cutnrunKO_mg_precountdata_int, cutnrunKO_consensus_peaks.gene, by="interval")
# cutnrunKO_mg_precountdata_int_ann_nearest <- cutnrunKO_mg_precountdata_int_ann %>% dplyr::filter(Distance < 5000)
# 
# cutnrunKO_mg_geneagg_countdata <- aggregate(cutnrunKO_mg_precountdata_int_ann_nearest[,c(2:5)],by=list(cutnrunKO_mg_precountdata_int_ann_nearest$Gene), sum)
# rownames(cutnrunKO_mg_geneagg_countdata) <- cutnrunKO_mg_geneagg_countdata$Group.1 
# cutnrunKO_mg_geneagg_countdata <- cutnrunKO_mg_geneagg_countdata[,-1]
# 
# # create sample-specific coldata object
# cutnrunKO_mg_dds_coldata <- cutnrunKOdds_coldata
# 
# # prepare count matrix 
# cutnrunKO_mg_countdata = as.matrix(cutnrunKO_mg_geneagg_countdata) 
# all(rownames(cutnrunKO_mg_dds_coldata) == colnames(cutnrunKO_mg_countdata)) #should print TRUE
# cutnrunKO_mg_ddsO <- DESeqDataSetFromMatrix(countData = cutnrunKO_mg_countdata, colData = cutnrunKO_mg_dds_coldata, design = ~ condition)
# cutnrunKO_mg_keep <- rowSums(counts(cutnrunKO_mg_ddsO)) > 10
# cutnrunKO_mg_ddsO <- cutnrunKO_mg_ddsO[cutnrunKO_mg_keep,]
# 
# # pca 
# cutnrunKO_mg_vst0Norm <- DESeq2::vst(cutnrunKO_mg_ddsO)
# cutnrunKO_mg_vst0Norm_pca <- plotPCA(cutnrunKO_mg_vst0Norm, intgroup="condition", returnData=TRUE)
# plotPCA(cutnrunKO_mg_vst0Norm, intgroup="condition", returnData=FALSE)
# 
# # de analysis
# cutnrunKO_mg_ddsO <- DESeq(cutnrunKO_mg_ddsO)
# 
# 
# # KO vs WT
# cutnrunKO_mg_de_KO_WT <- results(cutnrunKO_mg_ddsO, contrast=c("condition", "KO", "WT"))
# cutnrunKO_mg_de_KO_WT = cutnrunKO_mg_de_KO_WT[order(rownames(cutnrunKO_mg_de_KO_WT)),]
# cutnrunKO_mg_de_KO_WT$threshold <- as.logical(cutnrunKO_mg_de_KO_WT$padj < 0.05)
# deseq2_cutnrunKO_mg_de_KO_WT_0.05 <- cutnrunKO_mg_de_KO_WT[which(cutnrunKO_mg_de_KO_WT$threshold == TRUE),]
# deseq2_cutnrunKO_mg_de_KO_WT_0.05_df <- data.frame(deseq2_cutnrunKO_mg_de_KO_WT_0.05)
# deseq2_cutnrunKO_mg_de_KO_WT_0.05_df["Gene"] <- rownames(deseq2_cutnrunKO_mg_de_KO_WT_0.05_df)
# deseq2_cutnrunKO_mg_de_KO_WT_0.05_df <- deseq2_cutnrunKO_mg_de_KO_WT_0.05_df[,c(8,1:7)]
# deseq2_cutnrunKO_mg_de_KO_WT_0.05_fc_df <- deseq2_cutnrunKO_mg_de_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange > 1.5) | (log2FoldChange < -1.5))
# deseq2_cutnrunKO_mg_de_KO_WT_0.05_fc_up <- deseq2_cutnrunKO_mg_de_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange > 1.5))
# deseq2_cutnrunKO_mg_de_KO_WT_0.05_fc_down <- deseq2_cutnrunKO_mg_de_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange < -1.5))
# 
# 
# # log2(1.5) = 2.83
# ggmaplot(cutnrunKO_mg_de_KO_WT,
#          fdr = 0.05, fc = 2.83, size = 0.3,
#          palette = c("#D22B2B", "#1465AC", "darkgray"),
#          genenames = rownames(cutnrunKO_mg_de_KO_WT),
#          legend = "top", top = 20,
#          font.label = c("bold", 5),
#          font.legend = "bold",
#          font.main = "bold",
#          ggtheme = ggplot2::theme_minimal())
# 
# 


#####============================================================================================================================#####

# cutnrun KD
load("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/deseq2/H3K27me3.consensus_peaks.dds.RData")
cutnrunKDdds <- dds
rm(dds)
cutnrunKDprecountdata <- counts(cutnrunKDdds, normalized=FALSE)
head(cutnrunKDprecountdata)
# rearrange
cutnrunKDprecountdata <- cutnrunKDprecountdata[,c(3,2,1)]
dim(cutnrunKDprecountdata)
colSums(cutnrunKDprecountdata)
cutnrunKDprecountdata <- data.frame(cutnrunKDprecountdata)
write.table(cutnrunKDprecountdata, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/cutnrunKDprecountdata.txt", quote = F, append = F, sep="\t")

# NOISeq based analysis
# create sample-specific coldata object
library(NOISeq)
#NOISeq-real
cutnrunKDdds_factors <- data.frame(colData(cutnrunKDdds))

cutnrunKDdds_factors["condition"] <- unlist(lapply(strsplit(cutnrunKDdds_factors$sample, "_"), function(x) x[1]))
cutnrunKDdds_factors["replicate"] <- c(1,1,1)


# rearrange
cutnrunKDdds_factors <- cutnrunKDdds_factors[c(3,2,1),c(1,3,4,2)]

cutnrunKDdds_factors$condition <- factor(cutnrunKDdds_factors$condition)
cutnrunKDdds_factors$replicate <- factor(cutnrunKDdds_factors$replicate)


# annotate 
cutnrunKD_consensus_peaks_annotatePeaks <- data.frame(fread("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/H3K27me3.consensus_peaks.annotatePeaks.txt", header = T))
colnames(cutnrunKD_consensus_peaks_annotatePeaks) <- c("interval",colnames(cutnrunKD_consensus_peaks_annotatePeaks)[2:length(cutnrunKD_consensus_peaks_annotatePeaks)])
rownames(cutnrunKD_consensus_peaks_annotatePeaks) <- cutnrunKD_consensus_peaks_annotatePeaks$interval
cutnrunKD_consensus_peaks_annotatePeaks["interval_length"] <- cutnrunKD_consensus_peaks_annotatePeaks$End - cutnrunKD_consensus_peaks_annotatePeaks$Start
cutnrunKDprecountdata["interval"] <- rownames(cutnrunKDprecountdata)
cutnrunKDprecountdata_ann <- merge(cutnrunKDprecountdata, cutnrunKD_consensus_peaks_annotatePeaks, by="interval")
rownames(cutnrunKDprecountdata_ann) <- cutnrunKDprecountdata_ann$interval
cutnrunKDcountdata <- cutnrunKDprecountdata_ann[,c(2,3,4)]


cutnrunKD_length <-  data.frame(length=cutnrunKDprecountdata_ann$interval_length)
rownames(cutnrunKD_length) <- rownames(cutnrunKDprecountdata_ann)

cutnrunKD_biotype <-  data.frame(length=cutnrunKDprecountdata_ann$Gene.Type)
rownames(cutnrunKD_biotype) <- rownames(cutnrunKDprecountdata_ann)

#make noiseq object
cutnrunKD_nsobj <- readData(data = (cutnrunKDcountdata),  length = cutnrunKD_length$length, biotype = cutnrunKD_biotype, factors = cutnrunKDdds_factors)
cutnrunKD_myPCA = dat(cutnrunKD_nsobj, type = "PCA")
explo.plot(cutnrunKD_myPCA, factor = "condition")

cutnrunKD_myRPKM = rpkm(assayData(cutnrunKD_nsobj)$exprs, long = cutnrunKD_length$length, k = 0, lc = 1)
cutnrunKD_myRPKM <- data.frame(cutnrunKD_myRPKM)
cutnrunKD_myRPKM["interval"] <- rownames(cutnrunKD_myRPKM)
cutnrunKD_myRPKM_ann <- merge(cutnrunKD_myRPKM, cutnrunKD_consensus_peaks_annotatePeaks, by="interval")


cutnrunKD_myRPKM_ann["shCRAMP1_shControl"] <- log(cutnrunKD_myRPKM_ann$shCRAMP1_H3K27me3 / cutnrunKD_myRPKM_ann$shControl_H3K27me3, 2)
cutnrunKD_myRPKM_ann["shSUZ12_shControl"] <- log(cutnrunKD_myRPKM_ann$shSUZ12_H3K27me3 / cutnrunKD_myRPKM_ann$shControl_H3K27me3, 2)




#NOIseq-real: Two replicates in one of the experimental conditions are enough to run the algorithm, But I don't have two replicates either
#NOISeq BIO cannot be run becuase it two replicates in both condition
# NOISeq-slim can be run if I treat each individual replicates as separate
cutnrunKDdds_factors_slim <- data.frame(colData(cutnrunKDdds))

cutnrunKDdds_factors_slim["condition"] <- unlist(lapply(strsplit(cutnrunKDdds_factors_slim$sample, "_"), function(x) x[1]))


# rearrange
cutnrunKDdds_factors_slim <- cutnrunKDdds_factors_slim[c(3,2,1),c(1,3,2)]

cutnrunKDdds_factors_slim$condition <- factor(cutnrunKDdds_factors_slim$condition)

#make noiseq object
cutnrunKD_slim_nsobj <- readData(data = (cutnrunKDcountdata),  length = cutnrunKD_length$length, biotype = cutnrunKD_biotype, factors = cutnrunKDdds_factors_slim)
cutnrunKD_slim_myPCA = dat(cutnrunKD_slim_nsobj, type = "PCA")
explo.plot(cutnrunKD_slim_myPCA, factor = "condition")


cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq <- noiseq(cutnrunKD_slim_nsobj, factor = "condition", conditions = c("shCRAMP1", "shControl"), k = NULL, norm = "rpkm", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "no")
cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq <- noiseq(cutnrunKD_slim_nsobj, factor = "condition", conditions = c("shSUZ12", "shControl"), k = NULL, norm = "rpkm", pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "no")

head(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq@results[[1]])
dim(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq@results[[1]])

head(cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq@results[[1]])
dim(cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq@results[[1]])



### select the differentially expressed features
cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.deg = degenes(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq, q = 0.8, M = NULL)
cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.up = degenes(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq, q = 0.8, M = "up")
cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.down = degenes(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq, q = 0.8, M = "down")



### select the differentially expressed features
cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.deg = degenes(cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq, q = 0.8, M = NULL)
cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.up = degenes(cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq, q = 0.8, M = "up")
cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.down = degenes(cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq, q = 0.8, M = "down")


#annotate
cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.up["interval"] <- rownames(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.up)
cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.up_ann <- merge(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.up, cutnrunKD_consensus_peaks_annotatePeaks, by="interval")

cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.down["interval"] <- rownames(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.down)
cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.down_ann <- merge(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.down, cutnrunKD_consensus_peaks_annotatePeaks, by="interval")


#annotate
cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.up["interval"] <- rownames(cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.up)
cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.up_ann <- merge(cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.up, cutnrunKD_consensus_peaks_annotatePeaks, by="interval")

cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.down["interval"] <- rownames(cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.down)
cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.down_ann <- merge(cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.down, cutnrunKD_consensus_peaks_annotatePeaks, by="interval")

cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.up_ann_nearest <- cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.up_ann 
cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.down_ann_nearest <- cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.down_ann 

cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.up_ann_nearest <- cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.up_ann 
cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.down_ann_nearest <- cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq.down_ann 


# MD plot
DE.plot(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq, q = 0.8, graphic = "MD")
# MD plot
DE.plot(cutnrunKD_slim_shSUZ12_shControl_myRPKM_noiseq, q = 0.8, graphic = "MD")


#Use manual annotation
sort -k1,1 -k2,2n H3K27me3.consensus_peaks.bed | grep chr > H3K27me3.consensus_peaks_chr.bed
sort -k1,1 -k2,2n /mnt/home3/reid/av638/rnaseq/iva_lab_oct23/gene_gencode_human__gencode_out.txt | grep chr > gene_gencode_human_gencode_out.sorted.chr.txt
bedtools closest -a H3K27me3.consensus_peaks_chr.bed -b gene_gencode_human_gencode_out.sorted.chr.txt -d > H3K27me3.consensus_peaks_chr.ann.bed

cutnrunKD_consensus_peaks.gene <- data.frame(fread("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary/macs2/broadPeak/consensus/H3K27me3/H3K27me3.consensus_peaks_chr.ann.bed", header = F))
cutnrunKD_consensus_peaks.gene <- cutnrunKD_consensus_peaks.gene[,c(4,12,13,14)]
colnames(cutnrunKD_consensus_peaks.gene) <- c("interval","ensID", "Gene", "Distance")
cutnrunKD_consensus_peaks.gene["ensID_Gene"] <- paste0(cutnrunKD_consensus_peaks.gene$ensID, "%", cutnrunKD_consensus_peaks.gene$Gene)

#annotate

cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.up_anno <- merge(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.up, cutnrunKD_consensus_peaks.gene, by="interval")
cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.down_anno <- merge(cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.down, cutnrunKD_consensus_peaks.gene, by="interval")

cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.up_anno_nearest <- cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.up_anno %>% dplyr::filter(Distance < 5000)
cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.down_anno_nearest <- cutnrunKD_slim_shCRAMP1_shControl_myRPKM_noiseq.down_anno %>% dplyr::filter(Distance < 5000)

#calculate CPM


cutnrunKD_myCPM <- edgeR::cpm(cutnrunKDcountdata)
cutnrunKD_myCPM <- data.frame((cutnrunKD_myCPM))
cutnrunKD_myCPM["interval"] <- rownames(cutnrunKD_myCPM)
cutnrunKD_myCPM_ann <- merge(cutnrunKD_myCPM, cutnrunKD_consensus_peaks_annotatePeaks, by="interval")
cutnrunKD_myCPM_anno <- merge(cutnrunKD_myCPM, cutnrunKD_consensus_peaks.gene, by="interval")

cutnrunKD_myCPM_ann_nearest <- cutnrunKD_myCPM_ann 

cutnrunKD_myCPM_anno_nearest <- cutnrunKD_myCPM_anno %>% dplyr::filter(Distance < 5000)

#---------------------------- CUT&RUN-KD-Seq Analysis after annotation to genes --------------------------------#

# # get countdata for homer annotated (genic)
# cutnrunKD_genic_precountdata_int <- cutnrunKDprecountdata
# 
# cutnrunKD_genic_precountdata_int["interval"] <- rownames(cutnrunKD_genic_precountdata_int)
# 
# cutnrunKD_genic_precountdata_int_ann <- merge(cutnrunKD_genic_precountdata_int, cutnrunKD_consensus_peaks_annotatePeaks, by="interval")
# cutnrunKD_genic_precountdata_int_ann_nearest <- cutnrunKD_genic_precountdata_int_ann 
# 
# cutnrunKD_genic_geneagg_countdata <- aggregate(cutnrunKD_genic_precountdata_int_ann_nearest[,c(2:4)],by=list(cutnrunKD_genic_precountdata_int_ann_nearest$Gene.Name), sum)
# rownames(cutnrunKD_genic_geneagg_countdata) <- cutnrunKD_genic_geneagg_countdata$Group.1 
# 
# cutnrunKD_genic_countdata <- cutnrunKD_genic_geneagg_countdata[,-1]
# 
# cutnrunKD_genic_slim_myCPM <- edgeR::cpm(cutnrunKD_genic_countdata)
# cutnrunKD_genic_slim_myCPM <- data.frame((cutnrunKD_genic_slim_myCPM))
# cutnrunKD_genic_slim_myCPM["Gene"] <- rownames(cutnrunKD_genic_slim_myCPM)
# 
# #manual annotation
# cutnrunKD_mg_precountdata_int <- cutnrunKDprecountdata
# 
# cutnrunKD_mg_precountdata_int["interval"] <- rownames(cutnrunKD_mg_precountdata_int)
# 
# 
# cutnrunKD_mg_precountdata_int_ann <- merge(cutnrunKD_mg_precountdata_int, cutnrunKD_consensus_peaks.gene, by="interval")
# cutnrunKD_mg_precountdata_int_ann_nearest <- cutnrunKD_mg_precountdata_int_ann %>% dplyr::filter(Distance < 5000)
# 
# cutnrunKD_mg_geneagg_countdata <- aggregate(cutnrunKD_mg_precountdata_int_ann_nearest[,c(2:4)],by=list(cutnrunKD_mg_precountdata_int_ann_nearest$Gene), sum)
# rownames(cutnrunKD_mg_geneagg_countdata) <- cutnrunKD_mg_geneagg_countdata$Group.1 
# 
# cutnrunKD_mg_countdata <- cutnrunKD_mg_geneagg_countdata[,-1]
# 
# cutnrunKD_mg_slim_myCPM <- edgeR::cpm(cutnrunKD_mg_countdata)
# 
# cutnrunKD_mg_slim_myCPM <- data.frame((cutnrunKD_mg_slim_myCPM))
# cutnrunKD_mg_slim_myCPM["Gene"] <- rownames(cutnrunKD_mg_slim_myCPM)

####=========================================================================================####
cutntag_path <- "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/"

./quantify.sh


cutntag_bamtobed <- list.files(cutntag_path, pattern = ".bam_chr.bed")
cutntag_bamtobed_list <- list()
for (bamtobedfiles in cutntag_bamtobed){
  print(bamtobedfiles)
  bamtobedfiles_bed <- data.frame(fread(paste0(cutntag_path,bamtobedfiles)))
  cutntag_bamtobed_list[bamtobedfiles] <- dim(bamtobedfiles_bed)[1]
}
cutntag_bamtobed_list_df <- do.call(rbind.data.frame, cutntag_bamtobed_list)
rownames(cutntag_bamtobed_list_df) <- names(cutntag_bamtobed_list)
cutntag_bamtobed_list_df["sample"] <- gsub(".mLb.clN.sorted","",gsub(".bam_chr.bed","",rownames(cutntag_bamtobed_list_df)))
colnames(cutntag_bamtobed_list_df) <- c("reads", "sample")

histone_marks <- list.files(cutntag_path, pattern = "_peaks_id.bed")
cutntag_df_coverage_list <- list()
cutntag_df_ReadsInside_list <- list()
for (histone_file in histone_marks){
  print(histone_file)
  histone_file_tag <- gsub("_peaks_id.bed","", histone_file)
  cutntag_files <- list.files(cutntag_path, pattern = histone_file_tag)
  cutntag_files <- cutntag_files[-which(cutntag_files %in% c(histone_file))]
  # print(cutntag_files)
  cutntag_list <- list()
  for (bedfiles in cutntag_files){
    # print(bedfiles)
    bedfiles_name <- gsub(".bam","",gsub(".mLb.clN.sorted.bam","",gsub("_coverage.bed", "", bedfiles)))
    cutntag_list[[bedfiles_name]] <- data.frame(fread(paste0(cutntag_path,bedfiles)))
  }
  cutntag_df <- do.call(cbind.data.frame, cutntag_list)
  print(dim(cutntag_df))
  cutntag_df <- cutntag_df[,c(1:3, seq(13,80,16))]
  print(paste0(bedfiles_name, "_", histone_file_tag))
  print(head(cutntag_df,1))
  cutntag_column_name <- cutntag_df[,c(1,2,3)]
  colnames(cutntag_column_name) <- c("chr", "start", "end")
  rownames(cutntag_df) <- paste0(cutntag_column_name$chr, "%",cutntag_column_name$start, "%",cutntag_column_name$end)
  cutntag_df <- cutntag_df[,c(4:8)]
  cutntag_df_coverage_list[[paste0(bedfiles_name)]] <- cutntag_df
  cutntag_df_colsum <- data.frame(colSums(cutntag_df))
  cutntag_df_colsum["sample"] <- gsub(paste0("_",histone_file_tag, ".V13"),"",rownames(cutntag_df_colsum))
  cutntag_df_colsum["id"] <- rownames(cutntag_df_colsum)
  cutntag_bamtobed_colsum_df <- merge(cutntag_df_colsum, cutntag_bamtobed_list_df, by="sample")
  cutntag_bamtobed_colsum_df <- cutntag_bamtobed_colsum_df[,c(1,3,2,4)]
  colnames(cutntag_bamtobed_colsum_df) <- c("sample", "id", "ReadsInside", "TotalReads")
  cutntag_bamtobed_colsum_df["ReadsOutside"] <- cutntag_bamtobed_colsum_df$TotalReads - cutntag_bamtobed_colsum_df$ReadsInside
  cutntag_bamtobed_colsum_df["Ratio"] <- cutntag_bamtobed_colsum_df$ReadsInside / cutntag_bamtobed_colsum_df$ReadsOutside
  cutntag_df_ReadsInside_list[[paste0(histone_file_tag)]] <- cutntag_bamtobed_colsum_df
}

cutntag_df_ReadsInside_list_df <- do.call(cbind.data.frame, cutntag_df_ReadsInside_list)
cutntag_df_ReadsInside_sub_df <- cutntag_df_ReadsInside_list_df[,c(1, seq(6,114,6))]
rownames(cutntag_df_ReadsInside_sub_df) <- cutntag_df_ReadsInside_sub_df$H2AFZ_ENCFF213OTI.sample
cutntag_df_ReadsInside_sub_df <- cutntag_df_ReadsInside_sub_df[,-1]

cutntag_df_ReadsInside_sub_df_st <- data.frame(stack(as.matrix(cutntag_df_ReadsInside_sub_df)))
cutntag_df_ReadsInside_sub_df_st$col <- gsub(".Ratio","",cutntag_df_ReadsInside_sub_df_st$col)
cutntag_df_ReadsInside_sub_df_st_labels <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/histone_marks_label.txt", header = F))
cutntag_df_ReadsInside_sub_df_st_labels$V1 <- gsub("_peaks_id.bed", "", cutntag_df_ReadsInside_sub_df_st_labels$V1)
colnames(cutntag_df_ReadsInside_sub_df_st_labels) <- c("col", "label")
cutntag_df_ReadsInside_sub_df_st_re <- merge(cutntag_df_ReadsInside_sub_df_st_labels, cutntag_df_ReadsInside_sub_df_st, by="col")

cutntag_df_ReadsInside_sub_df_st_re$row <- factor(cutntag_df_ReadsInside_sub_df_st_re$row, levels = c("IgG","ENCFF508LLH","H1_4","H3K27me3","panH1"))

# subset by samples I want, viz. H1, panH1, H3K27me3
cutntag_df_ReadsInside_sub_df_st_H3K27me3 <- cutntag_df_ReadsInside_sub_df_st_re[cutntag_df_ReadsInside_sub_df_st_re$row %in% c("H3K27me3", "IgG"), ]
cutntag_df_ReadsInside_sub_df_st_H1_4 <- cutntag_df_ReadsInside_sub_df_st_re[cutntag_df_ReadsInside_sub_df_st_re$row %in% c("H1_4", "IgG"), ]
cutntag_df_ReadsInside_sub_df_st_panH1 <- cutntag_df_ReadsInside_sub_df_st_re[cutntag_df_ReadsInside_sub_df_st_re$row %in% c("panH1", "IgG"), ]

cutntag_plot1 <- ggplot(cutntag_df_ReadsInside_sub_df_st_H3K27me3, aes(x = label, y = value)) +
  geom_col(aes(color = row, fill = row), position = position_dodge(0.7), width = 0.6) +
  scale_color_manual(values = c("black","black"))+
  scale_fill_manual(values = c("#d3d3d3", "#fde396"))+theme_classic()+
  theme(axis.text.x = element_text(angle = 90), text = element_text(size=25), plot.title = element_text(hjust = 0.5))+
  labs(title = "H3K27me3", y = "ratio = Reads_Inside/Reads_outside ", x = "Histone Marks")


cutntag_plot2 <- ggplot(cutntag_df_ReadsInside_sub_df_st_H1_4, aes(x = label, y = value)) +
  geom_col(aes(color = row, fill = row), position = position_dodge(0.7), width = 0.6) +
  scale_color_manual(values = c("black","black"))+
  scale_fill_manual(values = c("#d3d3d3", "#ffdfba"))+theme_classic()+
  theme(axis.text.x = element_text(angle = 90), text = element_text(size=25), plot.title = element_text(hjust = 0.5))+
  labs(title = "H1_4", y = "ratio = Reads_Inside/Reads_outside ", x = "Histone Marks")

cutntag_plot3 <- ggplot(cutntag_df_ReadsInside_sub_df_st_panH1, aes(x = label, y = value)) +
  geom_col(aes(color = row, fill = row), position = position_dodge(0.7), width = 0.6) +
  scale_color_manual(values = c("black","black"))+
  scale_fill_manual(values = c("#d3d3d3", "#baffc9"))+theme_classic()+
  theme(axis.text.x = element_text(angle = 90), text = element_text(size=25), plot.title = element_text(hjust = 0.5))+
  labs(title = "panH1", y = "ratio = Reads_Inside/Reads_outside ", x = "Histone Marks")

# Combine the plots using patchwork
cutntag_combined_plot <- (cutntag_plot1 | cutntag_plot2 | cutntag_plot3)
cutntag_combined_plot <- cutntag_plot1 + cutntag_plot2 + cutntag_plot3
grid.arrange(cutntag_plot1, cutntag_plot2, cutntag_plot3, ncol = 1)
#save as  cutntag_combined_plot.pdf at pdf 30X20 portrait


# Over genes extended with promoter (1500bp)
#Annotation to regions
gene_v41_human_out <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/gene_v41_human_out.txt", header = F))
gene_v41_human_out_plus <- gene_v41_human_out[which(gene_v41_human_out$V4 == "+"),]
gene_v41_human_out_minus <- gene_v41_human_out[which(gene_v41_human_out$V4 == "-"),]

# extend TSS by 1500 bp
gene_v41_human_out_plus$V2 <- gene_v41_human_out_plus$V2 - 1500
gene_v41_human_out_minus$V3 <- gene_v41_human_out_minus$V3 + 1500

write.table(gene_v41_human_out_plus, "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/gene_v41_human_out_plus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)
write.table(gene_v41_human_out_minus, "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/gene_v41_human_out_minus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)

cat gene_v41_human_out_plus.txt gene_v41_human_out_minus.txt | sort -k1,1 -k2,2n > gene_v41_human_out_plus_minus.txt
grep "chrM" gene_v41_human_out_plus_minus.txt -v > gene_v41_human_out_plus_minus_freeM.txt

./quantify_over_genes.sh

cutntag_gene_cov_files <- list.files(cutntag_path, pattern = "_gene_coverage.bed")
cutntag_gene_coverage_list <- list()

for (gene_cov_file in cutntag_gene_cov_files){
  print(gene_cov_file)
  gene_cov_file_tag <- gsub(".bam_gene_coverage.bed","", gene_cov_file)
  print(gene_cov_file_tag)
  cutntag_gene_coverage_list[[gene_cov_file_tag]] <- data.frame(fread(paste0(cutntag_path,gene_cov_file)))
}

cutntag_gene_coverage_df <- do.call(cbind.data.frame, cutntag_gene_coverage_list)

cutntag_gene_precountdata <- cutntag_gene_coverage_df[,c(6,7, seq(8,dim(cutntag_gene_coverage_df)[2],11))]
colnames(cutntag_gene_precountdata) <- gsub(".mLb.clN.sorted.V8","",colnames(cutntag_gene_precountdata))
colnames(cutntag_gene_precountdata) <- gsub(".V8","",colnames(cutntag_gene_precountdata))
colnames(cutntag_gene_precountdata) <- c("ensID","Gene",colnames(cutntag_gene_precountdata)[3:7])
cutntag_gene_precountdata["ensID_Gene"] <- paste0(cutntag_gene_precountdata$ensID, "%", cutntag_gene_precountdata$Gene)
cutntag_gene_precountdata <- cutntag_gene_precountdata[,c(8,3:7)]
cutntag_gene_countdata <- aggregate(cutntag_gene_precountdata[,c(2:6)],by=list(cutntag_gene_precountdata$ensID_Gene), sum)

rownames(cutntag_gene_countdata) <- cutntag_gene_countdata$Group.1
cutntag_gene_countdata <- cutntag_gene_countdata[,-1]
#remove ENCODE H3K27me3
cutntag_gene_countdata <- cutntag_gene_countdata[,-1]

cutntag_gene_myCPM <- edgeR::cpm(cutntag_gene_countdata)
cutntag_gene_myCPM <- data.frame(cutntag_gene_myCPM)
cutntag_gene_myCPM["Gene"] <- rownames(cutntag_gene_myCPM)
cutntag_gene_myCPM_log <- log(cutntag_gene_myCPM[,c(1:4)] + 1,2)
cutntag_gene_myCPM_log["H1_4_IgG"] <- cutntag_gene_myCPM_log$H1_4 - cutntag_gene_myCPM_log$IgG
cutntag_gene_myCPM_log["H3K27me3_IgG"] <- cutntag_gene_myCPM_log$H3K27me3 - cutntag_gene_myCPM_log$IgG
cutntag_gene_myCPM_log["panH1_IgG"] <- cutntag_gene_myCPM_log$panH1 - cutntag_gene_myCPM_log$IgG
cutntag_gene_myCPM_log["Gene"] <- rownames(cutntag_gene_myCPM_log)
# cutntag_gene_myCPM_log_enriched <- cutntag_gene_myCPM_log[,c(5:7)]
# cutntag_gene_myCPM_log_overall_enriched <- cutntag_gene_myCPM_log_enriched[which(cutntag_gene_myCPM_log_enriched$H1_4_IgG > 0.3),]
# cutntag_gene_myCPM_log_overall_enriched <- cutntag_gene_myCPM_log_overall_enriched[which(cutntag_gene_myCPM_log_overall_enriched$H3K27me3_IgG > 0.3),]
# cutntag_gene_myCPM_log_overall_enriched <- cutntag_gene_myCPM_log_overall_enriched[which(cutntag_gene_myCPM_log_overall_enriched$panH1_IgG > 0.3),]
# dim(cutntag_gene_myCPM_log_overall_enriched)






# K27_peggy_Total <- 21061179
# H3K27me3_Total <- 3765724
# IgG_Total <- 768932
# H1_Total <- 3551411
# panH1_Total <- 2742886
# 
# K27_peggy_bound <- 5883369
# H3K27me3_bound <- 2371952
# IgG_bound <- 65557
# H1_bound <- 584840
# panH1_bound <- 269995
# 
# matrix_bound <- data.frame(ReadsInside = c(K27_peggy_bound, H3K27me3_bound, IgG_bound, H1_bound, panH1_bound),
#            Readsoutside = c(K27_peggy_Total - K27_peggy_bound, H3K27me3_Total - H3K27me3_bound, IgG_Total -IgG_bound, H1_Total - H1_bound, panH1_Total- panH1_bound))
# 
# rownames(matrix_bound) <- c("K27_peggy", "H3K27m3", "IgG", "H1", "panH1")
# test <- matrix_bound[c(4,3),]
# 
# chisq.test(test)


#pca on raw counts


cutntag_pca_res <- prcomp(t(cutntag_df), scale. = TRUE)

autoplot(cutntag_pca_res, label = TRUE, label.size = 3)

# # Integration of all omics data
# 
# # A). HOMER annotation, Peak DE first and then Gene, hpgd
# deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann_nearest     
# deseq2_rnaseqde_shC1_c1509_0.05_fc_df
# deseq2_cutnrunKOde_KO_WT_0.05_fc_ann_nearest
# cutnrunKD_myCPM_ann_nearest
# cutntag_gene_myCPM
# # consesus ann
# hpgd_genes_akd_rkd <- merge(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_ann_nearest, deseq2_rnaseqde_shC1_c1509_0.05_fc_df, by.x="Gene.Name", by.y="gene", all.x=TRUE, all.y=TRUE)
# hpgd_genes_akd_rkd <- hpgd_genes_akd_rkd[,c(1,4,28)]
# hpgd_genes_akd_rkd_cko <- merge(hpgd_genes_akd_rkd, deseq2_cutnrunKOde_KO_WT_0.05_fc_ann_nearest, by.x="Gene.Name", by.y="Gene.Name", all.x=TRUE, all.y=TRUE)
# hpgd_genes_akd_rkd_cko <- hpgd_genes_akd_rkd_cko[,c(1,2,3,6)]
# colnames(hpgd_genes_akd_rkd_cko) <- c("Gene", "akd_shCRAMP1_shControl", "rkd_shC1_c1509", "cko_KO_WT")
# 
# # In this strategy since different interval might get annotat to same genes yo will non unique gene name so rownames shoube adjusted as 
# rownames(hpgd_genes_akd_rkd_cko) <- paste0(rownames(hpgd_genes_akd_rkd_cko), "_",hpgd_genes_akd_rkd_cko$Gene)
# 
# #remove Y_RNA
# hpgd_genes_akd_rkd_cko <- hpgd_genes_akd_rkd_cko[which(hpgd_genes_akd_rkd_cko$Gene != "Y_RNA"),]
# 
# hpgd_genes_akd_rkd_cko <- hpgd_genes_akd_rkd_cko[,-1]
# 
# # remove data if no interpretation can be made about gene i.e. NA
# hpgd_genes_akd_rkd_cko <- na.omit(hpgd_genes_akd_rkd_cko)
# 
# pheatmap(as.matrix(hpgd_genes_akd_rkd_cko), fontsize_row = 1,cluster_rows = T,cluster_cols = F)
# 
# 
# #pheatmap_hpgd_genes_akd_rkd_cko.pdf
# 
# 
# #HOMER annotation, Gene annotation first and then DE, hgd
# deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_fc_df
# deseq2_rnaseqde_shC1_c1509_0.05_fc_df
# deseq2_cutnrunKO_genic_de_KO_WT_0.05_fc_df
# cutnrunKD_genic_slim_myCPM
# cutntag_gene_myCPM
# 
# hgd_genes_akd_rkd <- merge(deseq2_atacseqkd_genic_de_shCRAMP1_shControl_0.05_fc_df, deseq2_rnaseqde_shC1_c1509_0.05_fc_df, by.x="Gene", by.y="gene", all.x=TRUE, all.y=TRUE)
# hgd_genes_akd_rkd <- hgd_genes_akd_rkd[,c(1,3,10)]
# hgd_genes_akd_rkd_cko <- merge(hgd_genes_akd_rkd, deseq2_cutnrunKO_genic_de_KO_WT_0.05_fc_df, by.x="Gene", by.y="Gene", all.x=TRUE, all.y=TRUE)
# hgd_genes_akd_rkd_cko <- hgd_genes_akd_rkd_cko[,c(1,2,3,5)]
# colnames(hgd_genes_akd_rkd_cko) <- c("Gene", "akd_shCRAMP1_shControl", "rkd_shC1_c1509", "cko_KO_WT")
# 
# rownames(hgd_genes_akd_rkd_cko) <- hgd_genes_akd_rkd_cko$Gene
# 
# #remove Y_RNA
# hgd_genes_akd_rkd_cko <- hgd_genes_akd_rkd_cko[which(hgd_genes_akd_rkd_cko$Gene != "Y_RNA"),]
# 
# hgd_genes_akd_rkd_cko <- hgd_genes_akd_rkd_cko[,-1]
# 
# # remove data if no interpretation can be made about gene i.e. NA
# hgd_genes_akd_rkd_cko <- na.omit(hgd_genes_akd_rkd_cko)
# 
# pheatmap(as.matrix(hgd_genes_akd_rkd_cko), fontsize_row = 1,cluster_rows = T,cluster_cols = F)
# 
# 
# #pheatmap_hgd_genes_akd_rkd_cko.pdf
# 
# #MANUAL annotation, Peak DE first and then Gene, mpgd
# deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_anno_nearest
# deseq2_rnaseqde_shC1_c1509_0.05_fc_df
# deseq2_cutnrunKOde_KO_WT_0.05_fc_anno_nearest
# cutnrunKD_myCPM_anno_nearest
# cutntag_gene_myCPM
# 
# # consesus ann
# mpgd_genes_akd_rkd <- merge(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_anno_nearest, deseq2_rnaseqde_shC1_c1509_0.05_fc_df, by.x="Gene", by.y="gene", all.x=TRUE, all.y=TRUE)
# mpgd_genes_akd_rkd <- mpgd_genes_akd_rkd[,c(1,4,12)]
# mpgd_genes_akd_rkd_cko <- merge(mpgd_genes_akd_rkd, deseq2_cutnrunKOde_KO_WT_0.05_fc_anno_nearest, by.x="Gene", by.y="Gene", all.x=TRUE, all.y=TRUE)
# mpgd_genes_akd_rkd_cko <- mpgd_genes_akd_rkd_cko[,c(1,2,3,6)]
# colnames(mpgd_genes_akd_rkd_cko) <- c("Gene", "akd_shCRAMP1_shControl", "rkd_shC1_c1509", "cko_KO_WT")
# 
# # In this strategy since different interval might get annotat to same genes yo will non unique gene name so rownames shoube adjusted as 
# rownames(mpgd_genes_akd_rkd_cko) <- paste0(rownames(mpgd_genes_akd_rkd_cko), "_",mpgd_genes_akd_rkd_cko$Gene)
# 
# #remove Y_RNA
# mpgd_genes_akd_rkd_cko <- mpgd_genes_akd_rkd_cko[which(mpgd_genes_akd_rkd_cko$Gene != "Y_RNA"),]
# 
# mpgd_genes_akd_rkd_cko <- mpgd_genes_akd_rkd_cko[,-1]
# 
# # remove data if no interpretation can be made about gene i.e. NA
# mpgd_genes_akd_rkd_cko <- na.omit(mpgd_genes_akd_rkd_cko)
# 
# pheatmap(as.matrix(mpgd_genes_akd_rkd_cko), fontsize_row = 1,cluster_rows = T,cluster_cols = T)
# 
# 
# #pheatmap_mpgd_genes_akd_rkd_cko.pdf
# 
# 
# 
# #MANUAL annotation, Gene annotation first and then DE, mgd
# deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_fc_df
# deseq2_rnaseqde_shC1_c1509_0.05_fc_df
# deseq2_cutnrunKO_mg_de_KO_WT_0.05_fc_df
# cutnrunKD_mg_slim_myCPM
# cutntag_gene_myCPM
# 
# 
# mgd_genes_akd_rkd <- merge(deseq2_atacseqkd_mg_de_shCRAMP1_shControl_0.05_fc_df, deseq2_rnaseqde_shC1_c1509_0.05_fc_df, by.x="Gene", by.y="gene", all.x=TRUE, all.y=TRUE)
# mgd_genes_akd_rkd <- mgd_genes_akd_rkd[,c(1,3,10)]
# mgd_genes_akd_rkd_cko <- merge(mgd_genes_akd_rkd, deseq2_cutnrunKO_mg_de_KO_WT_0.05_fc_df, by.x="Gene", by.y="Gene", all.x=TRUE, all.y=TRUE)
# mgd_genes_akd_rkd_cko <- mgd_genes_akd_rkd_cko[,c(1,2,3,5)]
# colnames(mgd_genes_akd_rkd_cko) <- c("Gene", "akd_shCRAMP1_shControl", "rkd_shC1_c1509", "cko_KO_WT")
# 
# # In this strategy since different interval might get annotat to same genes yo will non unique gene name so rownames shoube adjusted as 
# rownames(mgd_genes_akd_rkd_cko) <- mgd_genes_akd_rkd_cko$Gene
# 
# #remove Y_RNA
# mgd_genes_akd_rkd_cko <- mgd_genes_akd_rkd_cko[which(mgd_genes_akd_rkd_cko$Gene != "Y_RNA"),]
# 
# mgd_genes_akd_rkd_cko <- mgd_genes_akd_rkd_cko[,-1]
# 
# # remove data if no interpretation can be made about gene i.e. NA
# mgd_genes_akd_rkd_cko <- na.omit(mgd_genes_akd_rkd_cko)
# 
# pheatmap(as.matrix(mgd_genes_akd_rkd_cko), fontsize_row = 1,cluster_rows = T,cluster_cols = F)
# rownames(na.omit(hgd_genes_akd_rkd_cko))
# 
# #pheatmap_mgd_genes_akd_rkd_cko.pdf
# 
# length(unique(unlist(lapply(strsplit(rownames(na.omit(hgd_genes_akd_rkd_cko)), "_"), function(x) x[2]))))
# length(unique(unlist(lapply(strsplit(rownames(na.omit(mgd_genes_akd_rkd_cko)), "_"), function(x) x[2]))))
# length(unique(unlist(lapply(strsplit(rownames(na.omit(hpgd_genes_akd_rkd_cko)), "_"), function(x) x[2]))))
# length(unique(unlist(lapply(strsplit(rownames(na.omit(mpgd_genes_akd_rkd_cko)), "_"), function(x) x[2]))))
# 
# 
# # Next merge this with ckd and ct
# #MANUAL annotation, Peak DE first and then Gene, mpgd
# deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_anno_nearest
# deseq2_rnaseqde_shC1_c1509_0.05_fc_df
# deseq2_cutnrunKOde_KO_WT_0.05_fc_anno_nearest
# cutnrunKD_myCPM_anno_nearest
# cutntag_gene_myCPM

#get single value per data types
#mean log2FC
deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_anno_nearest_agg <- aggregate(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_anno_nearest[,c(3)],by=list(deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_anno_nearest$ensID_Gene), mean)
deseq2_atacseqkd_de_shSUZ12_shControl_0.05_anno_nearest_agg <- aggregate(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_anno_nearest[,c(3)],by=list(deseq2_atacseqkd_de_shSUZ12_shControl_0.05_anno_nearest$ensID_Gene), mean)

#mean log2FC
deseq2_cutnrunKOde_KO_WT_0.05_fc_anno_nearest_agg <- aggregate(deseq2_cutnrunKOde_KO_WT_0.05_fc_anno_nearest[,c(3)],by=list(deseq2_cutnrunKOde_KO_WT_0.05_fc_anno_nearest$ensID_Gene), mean)
#sum CPM
cutnrunKD_myCPM_anno_nearest_agg <- aggregate(cutnrunKD_myCPM_anno_nearest[,c(2:4)],by=list(cutnrunKD_myCPM_anno_nearest$ensID_Gene), sum)
rownames(cutnrunKD_myCPM_anno_nearest_agg) <- cutnrunKD_myCPM_anno_nearest_agg$Group.1
cutnrunKD_myCPM_anno_nearest_agg_log <- log(cutnrunKD_myCPM_anno_nearest_agg[,c(2:4)] + 1,2)
cutnrunKD_myCPM_anno_nearest_agg_log["Group.1"] <- rownames(cutnrunKD_myCPM_anno_nearest_agg_log)
cutnrunKD_myCPM_anno_nearest_agg_log["H3K27me3_shCRAMP1_shControl"] <- cutnrunKD_myCPM_anno_nearest_agg_log$shCRAMP1_H3K27me3 - cutnrunKD_myCPM_anno_nearest_agg_log$shControl_H3K27me3
cutnrunKD_myCPM_anno_nearest_agg_log["H3K27me3_shSUZ12_shControl"] <- cutnrunKD_myCPM_anno_nearest_agg_log$shSUZ12_H3K27me3 - cutnrunKD_myCPM_anno_nearest_agg_log$shControl_H3K27me3

sel_mpgd_genes_rkdc_rkds <- merge(deseq2_rnaseqde_shC1_c1509_0.05_fc_df, deseq2_rnaseqde_shS_c1509_0.05_fc_df, by.x="ensID_Gene", by.y="ensID_Gene", all.x=TRUE, all.y=TRUE)
sel_mpgd_genes_rkdc_rkds <- sel_mpgd_genes_rkdc_rkds[,c(1,4,12)]

# consesus anno
sel_mpgd_genes_rkdc_rkds_akdc <- merge(sel_mpgd_genes_rkdc_rkds, deseq2_atacseqkd_de_shCRAMP1_shControl_0.05_anno_nearest_agg , by.x="ensID_Gene", by.y="Group.1", all.x=TRUE, all.y=TRUE)

sel_mpgd_genes_rkdc_rkds_akdc_akds <- merge(sel_mpgd_genes_rkdc_rkds_akdc, deseq2_atacseqkd_de_shSUZ12_shControl_0.05_anno_nearest_agg , by.x="ensID_Gene", by.y="Group.1", all.x=TRUE, all.y=TRUE)


sel_mpgd_genes_rkdc_rkds_akdc_akds_cko <- merge(sel_mpgd_genes_rkdc_rkds_akdc_akds, deseq2_cutnrunKOde_KO_WT_0.05_fc_anno_nearest_agg, by.x="ensID_Gene", by.y="Group.1", all.x=TRUE, all.y=TRUE)

colnames(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko) <- c("ensID_Gene","rkd_shC1_c1509","rkd_shS_c1509", "akd_shCRAMP1_shControl","akd_shSUZ12_shControl", "cko_KO_WT")

sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd <- merge(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko, cutnrunKD_myCPM_anno_nearest_agg_log, by.x="ensID_Gene", by.y="Group.1", all.x=TRUE, all.y=TRUE)

sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct <- merge(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd, cutntag_gene_myCPM_log, by.x="ensID_Gene", by.y="Gene", all.x=TRUE, all.y=TRUE)

rownames(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct) <- sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct$ensID_Gene

#remove Y_RNA
sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all <- sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct[which(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct$ensID_Gene != "Y_RNA"),]
#Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all[,c(2:6, 10:11,16:18)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE)


library(ComplexHeatmap)
# sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_mat1 <- Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all[,c(2:6)]), na_col = "white", row_names_gp = gpar(fontsize = 5), row_km = 3)
# sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_mat2 <- Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all[,c(10:11)]), na_col = "white", row_names_gp = gpar(fontsize = 5), row_km = 3)
# sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_mat3 <- Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all[,c(16:18)]), na_col = "white", row_names_gp = gpar(fontsize = 5), row_km = 3)
# 
# sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_mat1 + sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_mat2 + sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_mat3

#complex_heatmap_sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_ch.pdf

# sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative_mat <- Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative[,c(2:6, 10:11,16:18)]), na_col = "white", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE)
# Heatmap_sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative_mat.pdf


#All informative
sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative <- na.omit(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all)

Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative[,c(2:6, 10:11,16:18)]), na_col = "white", row_names_gp = gpar(fontsize = 5))

library(ComplexHeatmap)
sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative_mat1 <- Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative[,c(2:6)]), na_col = "white", row_names_gp = gpar(fontsize = 5))
sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative_mat2 <- Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative[,c(10:11)]), na_col = "white", row_names_gp = gpar(fontsize = 5))
sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative_mat3 <- Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative[,c(16:18)]), na_col = "white", row_names_gp = gpar(fontsize = 5))

sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative_mat1 + sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative_mat2 + sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative_mat3

#complex_heatmap_sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_informative_ch.pdf

col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
col_fun2 = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
# Some predictable
sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable <- sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all[complete.cases(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all[, c("rkd_shC1_c1509")]), ]
library(ComplexHeatmap)
sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable_mat1 <- Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable[,c(2:6)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, col = col_fun)
sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable_mat2 <- Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable[,c(10:11)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, col = col_fun2)
sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable_mat3 <- Heatmap(as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable[,c(16:18)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, col = col_fun)

sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable_mat1 + sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable_mat2 + sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable_mat3

#complex_heatmap_sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable_ch.pdf
#complex_heatmap_sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable_mat1.pdf
sel_tree_rkdc_rkds_akdc_akds_cko_sub <- as.matrix(sel_mpgd_genes_rkdc_rkds_akdc_akds_cko_ckd_ct_all_somepredictable[,c(2:6)])
breaksList3 = seq(-10, 10, by = 0.01)

# pheatmap_sel_tree_rkdc_rkds_akdc_akds_cko_sub <- pheatmap(sel_tree_rkdc_rkds_akdc_akds_cko_sub,
#          color = colorRampPalette(c("blue", "white", "red"))(length(breaksList3)),
#          breaks = breaksList3,
#          fontsize = 8,
#          cluster_cols = FALSE,
#          clustering_distance_cols = "euclidean", 
#          clustering_method = "ward.D", 
#          border_color = NA,
#          cutree_rows = 3)


# take genes of all 3 clusters manually
cluster1_genes <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/bigwig/cluster1_genes.txt", header = F))
cluster2_genes <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/bigwig/cluster2_genes.txt", header = F))
cluster3_genes <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/bigwig/cluster3_genes.txt", header = F))
colnames(cluster1_genes) <- "genes"
colnames(cluster2_genes) <- "genes"
colnames(cluster3_genes) <- "genes"

# awk '{print $1"\t"$2"\t"$3"\t"0"\t"$7"\t"$4}' gene_v41_human_out.txt > gene_v41_human_out_re.txt


gene_v41_human_out_re.bed <- data.frame(fread("/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/bigwig/gene_v41_human_out_re.txt", header = F))
colnames(gene_v41_human_out_re.bed) <- c("chr", "start", "end", "none","genes", "strand")
gene_v41_human_out_re_cluster1 <- merge(cluster1_genes, gene_v41_human_out_re.bed, by="genes")
gene_v41_human_out_re_cluster1 <- gene_v41_human_out_re_cluster1[which(gene_v41_human_out_re_cluster1$genes != "Metazoa_SRP"),]
write.table(gene_v41_human_out_re_cluster1[,c(2:4,1,5:6)], "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/bigwig/gene_v41_human_out_re_cluster1.txt", sep="\t", append = F, quote = F, col.names = F, row.names = F)
dim(gene_v41_human_out_re_cluster1)

gene_v41_human_out_re_cluster2 <- merge(cluster2_genes, gene_v41_human_out_re.bed, by="genes")
gene_v41_human_out_re_cluster2 <- gene_v41_human_out_re_cluster2[which(gene_v41_human_out_re_cluster2$genes != "Metazoa_SRP"),]
write.table(gene_v41_human_out_re_cluster2[,c(2:4,1,5:6)], "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/bigwig/gene_v41_human_out_re_cluster2.txt", sep="\t", append = F, quote = F, col.names = F, row.names = F)
dim(gene_v41_human_out_re_cluster2)

gene_v41_human_out_re_cluster3 <- merge(cluster3_genes, gene_v41_human_out_re.bed, by="genes")
gene_v41_human_out_re_cluster3 <- gene_v41_human_out_re_cluster3[which(gene_v41_human_out_re_cluster3$genes != "Metazoa_SRP"),]
write.table(gene_v41_human_out_re_cluster3[,c(2:4,1,5:6)], "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/bigwig/gene_v41_human_out_re_cluster3.txt", sep="\t", append = F, quote = F, col.names = F, row.names = F)
dim(gene_v41_human_out_re_cluster3)


# Notes on 30 Oct 2023
#cpm -> log2 cpm +1
# NA -> 0 -> make it white color instead of yellow
# Give uniform color to ckd and cutntag
# cutandtag convert to z-scores of CPM
# cutntag calculate enrichment over IgG
# add Suz12 kd rnaseq, atacseq, 
# add H1 rnaseq kd
# do the matrix for all genes without removing anything


#---------------- Bin-based analysis ------------#
# create genomic bins, 5kb, 10kb
bedtools makewindows -g hg38.chrom.sizes -w 5000 > hg38_5kb.txt
bedtools makewindows -g hg38.chrom.sizes -w 10000 > hg38_10kb.txt

# for PE, atac
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job_bb ./bamtobed.sh
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job_cb ./convert_bedpe_to_bed.sh
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job ./quantify_over_bins_pe.sh

# for all SE, cutnrun KO,KD and cutntag
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py ./bamtobed_se.sh
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job ./quantify_over_bins_se.sh
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job ./quantify_over_bins_se.sh
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job ./quantify_over_bins_se.sh

bin_cov_atacseq_path <- list.files(path="/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library", pattern = "*.txt_coverage.pe.bed", full.names = T)
bin_cov_cutnrun_KO_path <- list.files(path="/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary", pattern = "*.txt_coverage_se.bed", full.names = T)
bin_cov_cutnrun_KD_path <- list.files(path="/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary", pattern = "*.txt_coverage_se.bed", full.names = T)
bin_cov_cutntag_path <- list.files(path="/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary", pattern = "*.txt_coverage_se.bed", full.names = T)
bin_path_combined <- c(bin_cov_atacseq_path, bin_cov_cutnrun_KO_path, bin_cov_cutnrun_KD_path, bin_cov_cutntag_path)
bin_path_combined_5kb <- bin_path_combined[grepl("5kb", bin_path_combined)]
bin_path_combined_10kb <- bin_path_combined[grepl("10kb", bin_path_combined)]

bin_5kb_list <- list()
for (files in bin_path_combined_5kb){
  filename <- sub("\\..*", "", basename(files))
  print(filename)
  bin_5kb_list[[filename]] <- data.frame(fread(files, header = F))
}
bin_5kb_df <- do.call(cbind.data.frame, bin_5kb_list)
bin_5kb_df <- bin_5kb_df[,c(1:3, seq(4,175,7))]
rownames(bin_5kb_df) <- paste0(bin_5kb_df$shControl_REP1.V1, "%", bin_5kb_df$shControl_REP1.V2, "%", bin_5kb_df$shControl_REP1.V3)
bin_5kb_df <- bin_5kb_df[,-c(1:3)]
colnames(bin_5kb_df) <- gsub(".V4","", colnames(bin_5kb_df))
bin_5kb_df_filt <- bin_5kb_df[rowSums(bin_5kb_df) > 0,]
bin_5kb_df_CPM <- edgeR::cpm(bin_5kb_df_filt)
bin_5kb_df_CPM_log <- log(bin_5kb_df_CPM + 1,2)
bin_5kb_df_CPM_log <- data.frame(bin_5kb_df_CPM_log)
bin_5kb_df_CPM_log["akd_shCRAMP1_shControl_1"] <- bin_5kb_df_CPM_log$shCRAMP1_REP1 - bin_5kb_df_CPM_log$shControl_REP1
bin_5kb_df_CPM_log["akd_shCRAMP1_shControl_2"] <- bin_5kb_df_CPM_log$shCRAMP1_REP2 - bin_5kb_df_CPM_log$shControl_REP2
bin_5kb_df_CPM_log["akd_shSUZ12_shControl_1"] <- bin_5kb_df_CPM_log$shSUZ12_REP1 - bin_5kb_df_CPM_log$shControl_REP1
bin_5kb_df_CPM_log["akd_shSUZ12_shControl_2"] <- bin_5kb_df_CPM_log$shSUZ12_REP2 - bin_5kb_df_CPM_log$shControl_REP2
bin_5kb_df_CPM_log["cko_H3K27me3_KO1_IgG"] <- bin_5kb_df_CPM_log$H3K27me3_KO1 - bin_5kb_df_CPM_log$IgG_KO1
bin_5kb_df_CPM_log["cko_H3K27me3_KO2_IgG"] <- bin_5kb_df_CPM_log$H3K27me3_KO2 - bin_5kb_df_CPM_log$IgG_KO2
bin_5kb_df_CPM_log["cko_H3K27me3_KO3_IgG"] <- bin_5kb_df_CPM_log$H3K27me3_KO3 - bin_5kb_df_CPM_log$IgG_KO3
bin_5kb_df_CPM_log["cko_H3K27me3_WT_IgG"] <- bin_5kb_df_CPM_log$H3K27me3_WT - bin_5kb_df_CPM_log$IgG_WT
bin_5kb_df_CPM_log["ckd_shControl_H3K27me3_IgG"] <- bin_5kb_df_CPM_log$shControl_H3K27me3 - bin_5kb_df_CPM_log$shControl_IgG
bin_5kb_df_CPM_log["ckd_shCRAMP1_H3K27me3_IgG"] <- bin_5kb_df_CPM_log$shCRAMP1_H3K27me3 - bin_5kb_df_CPM_log$shCRAMP1_IgG
bin_5kb_df_CPM_log["ckd_shSUZ12_H3K27me3_IgG"] <- bin_5kb_df_CPM_log$shSUZ12_H3K27me3 - bin_5kb_df_CPM_log$shSUZ12_IgG
bin_5kb_df_CPM_log["ct_ENCFF508LLH_IgG"] <- bin_5kb_df_CPM_log$ENCFF508LLH - bin_5kb_df_CPM_log$IgG
bin_5kb_df_CPM_log["ct_H1_4_IgG"] <- bin_5kb_df_CPM_log$H1_4 - bin_5kb_df_CPM_log$IgG
bin_5kb_df_CPM_log["ct_H3K27me3_IgG"] <- bin_5kb_df_CPM_log$H3K27me3 - bin_5kb_df_CPM_log$IgG
bin_5kb_df_CPM_log["ct_panH1_IgG"] <- bin_5kb_df_CPM_log$panH1 - bin_5kb_df_CPM_log$IgG

bin_5kb_df_CPM_log_sub <- bin_5kb_df_CPM_log[,c(26:40)]

# 10 kb
bin_10kb_list <- list()
for (files in bin_path_combined_10kb){
  filename <- sub("\\..*", "", basename(files))
  print(filename)
  bin_10kb_list[[filename]] <- data.frame(fread(files, header = F))
}
bin_10kb_df <- do.call(cbind.data.frame, bin_10kb_list)
bin_10kb_df <- bin_10kb_df[,c(1:3, seq(4,175,7))]
rownames(bin_10kb_df) <- paste0(bin_10kb_df$shControl_REP1.V1, "%", bin_10kb_df$shControl_REP1.V2, "%", bin_10kb_df$shControl_REP1.V3)
bin_10kb_df <- bin_10kb_df[,-c(1:3)]
colnames(bin_10kb_df) <- gsub(".V4","", colnames(bin_10kb_df))
bin_10kb_df_filt <- bin_10kb_df[rowSums(bin_10kb_df) > 0,]
bin_10kb_df_CPM <- edgeR::cpm(bin_10kb_df_filt)
bin_10kb_df_CPM_log <- log(bin_10kb_df_CPM + 1,2)
bin_10kb_df_CPM_log <- data.frame(bin_10kb_df_CPM_log)
bin_10kb_df_CPM_log["akd_shCRAMP1_shControl_1"] <- bin_10kb_df_CPM_log$shCRAMP1_REP1 - bin_10kb_df_CPM_log$shControl_REP1
bin_10kb_df_CPM_log["akd_shCRAMP1_shControl_2"] <- bin_10kb_df_CPM_log$shCRAMP1_REP2 - bin_10kb_df_CPM_log$shControl_REP2
bin_10kb_df_CPM_log["akd_shSUZ12_shControl_1"] <- bin_10kb_df_CPM_log$shSUZ12_REP1 - bin_10kb_df_CPM_log$shControl_REP1
bin_10kb_df_CPM_log["akd_shSUZ12_shControl_2"] <- bin_10kb_df_CPM_log$shSUZ12_REP2 - bin_10kb_df_CPM_log$shControl_REP2
bin_10kb_df_CPM_log["cko_H3K27me3_KO1_IgG"] <- bin_10kb_df_CPM_log$H3K27me3_KO1 - bin_10kb_df_CPM_log$IgG_KO1
bin_10kb_df_CPM_log["cko_H3K27me3_KO2_IgG"] <- bin_10kb_df_CPM_log$H3K27me3_KO2 - bin_10kb_df_CPM_log$IgG_KO2
bin_10kb_df_CPM_log["cko_H3K27me3_KO3_IgG"] <- bin_10kb_df_CPM_log$H3K27me3_KO3 - bin_10kb_df_CPM_log$IgG_KO3
bin_10kb_df_CPM_log["cko_H3K27me3_WT_IgG"] <- bin_10kb_df_CPM_log$H3K27me3_WT - bin_10kb_df_CPM_log$IgG_WT
bin_10kb_df_CPM_log["ckd_shControl_H3K27me3_IgG"] <- bin_10kb_df_CPM_log$shControl_H3K27me3 - bin_10kb_df_CPM_log$shControl_IgG
bin_10kb_df_CPM_log["ckd_shCRAMP1_H3K27me3_IgG"] <- bin_10kb_df_CPM_log$shCRAMP1_H3K27me3 - bin_10kb_df_CPM_log$shCRAMP1_IgG
bin_10kb_df_CPM_log["ckd_shSUZ12_H3K27me3_IgG"] <- bin_10kb_df_CPM_log$shSUZ12_H3K27me3 - bin_10kb_df_CPM_log$shSUZ12_IgG
bin_10kb_df_CPM_log["ct_ENCFF508LLH_IgG"] <- bin_10kb_df_CPM_log$ENCFF508LLH - bin_10kb_df_CPM_log$IgG
bin_10kb_df_CPM_log["ct_H1_4_IgG"] <- bin_10kb_df_CPM_log$H1_4 - bin_10kb_df_CPM_log$IgG
bin_10kb_df_CPM_log["ct_H3K27me3_IgG"] <- bin_10kb_df_CPM_log$H3K27me3 - bin_10kb_df_CPM_log$IgG
bin_10kb_df_CPM_log["ct_panH1_IgG"] <- bin_10kb_df_CPM_log$panH1 - bin_10kb_df_CPM_log$IgG

bin_10kb_df_CPM_log_sub <- bin_10kb_df_CPM_log[,c(26:40)]


# DARs analysis using Bins
# bin_5kb_df_filt
atacseqkd_bin5kb_precountdata <- bin_5kb_df_filt[,c(1:6)]
head(atacseqkd_bin5kb_precountdata)
dim(atacseqkd_bin5kb_precountdata)
colSums(atacseqkd_bin5kb_precountdata)
atacseqkd_bin5kb_precountdata <- data.frame(atacseqkd_bin5kb_precountdata)
#write.table(atacseqkd_bin5kb_precountdata, "atacseqkd_bin5kb_precountdata.txt", quote = F, append = F, sep="\t")

# get countdata
atacseqkd_bin5kb_countdata <- atacseqkd_bin5kb_precountdata[,c(1:dim(atacseqkd_bin5kb_precountdata)[2])]

# create sample-specific coldata object
atacseqkd_bin5kb_dds_coldata <- data.frame(sample = colnames(atacseqkd_bin5kb_countdata), 
                                        condition = unlist(lapply(strsplit(colnames(atacseqkd_bin5kb_countdata), "_"), function(x) x[1])), 
                                        replicate = unlist(lapply(strsplit(colnames(atacseqkd_bin5kb_countdata), "_"), function(x) x[2])))

rownames(atacseqkd_bin5kb_dds_coldata) <- atacseqkd_bin5kb_dds_coldata$sample
atacseqkd_bin5kb_dds_coldata$condition <- factor(atacseqkd_bin5kb_dds_coldata$condition)
atacseqkd_bin5kb_dds_coldata$replicate <- factor(atacseqkd_bin5kb_dds_coldata$replicate)

# prepare count matrix 
atacseqkd_bin5kb_countdata = as.matrix(atacseqkd_bin5kb_countdata) 
all(rownames(atacseqkd_bin5kb_dds_coldata) == colnames(atacseqkd_bin5kb_countdata)) #should print TRUE
atacseqkd_bin5kb_ddsO <- DESeqDataSetFromMatrix(countData = atacseqkd_bin5kb_countdata, colData = atacseqkd_bin5kb_dds_coldata, design = ~ condition)
atacseqkd_bin5kb_keep <- rowSums(counts(atacseqkd_bin5kb_ddsO)) > 10
atacseqkd_bin5kb_ddsO <- atacseqkd_bin5kb_ddsO[atacseqkd_bin5kb_keep,]

# pca 
atacseqkd_bin5kb_vst0Norm <- DESeq2::vst(atacseqkd_bin5kb_ddsO)
atacseqkd_bin5kb_vst0Norm_pca <- plotPCA(atacseqkd_bin5kb_vst0Norm, intgroup="condition", returnData=TRUE)
plotPCA(atacseqkd_bin5kb_vst0Norm, intgroup="condition", returnData=FALSE)

# de analysis
atacseqkd_bin5kb_ddsO <- DESeq(atacseqkd_bin5kb_ddsO)


# shCRAMP1 vs shControl
atacseqkd_bin5kb_de_shCRAMP1_shControl <- results(atacseqkd_bin5kb_ddsO, contrast=c("condition", "shCRAMP1", "shControl"))
atacseqkd_bin5kb_de_shCRAMP1_shControl = atacseqkd_bin5kb_de_shCRAMP1_shControl[order(rownames(atacseqkd_bin5kb_de_shCRAMP1_shControl)),]
atacseqkd_bin5kb_de_shCRAMP1_shControl$threshold <- as.logical(atacseqkd_bin5kb_de_shCRAMP1_shControl$padj < 0.05)
deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05 <- atacseqkd_bin5kb_de_shCRAMP1_shControl[which(atacseqkd_bin5kb_de_shCRAMP1_shControl$threshold == TRUE),]
deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_df <- data.frame(deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05)
deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_df["interval"] <- rownames(deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_df)
deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_df <- deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_df[,c(8,1:7)]
deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_fc_df <- deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1) | (log2FoldChange < -1))
deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_fc_up <- deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1))
deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_fc_down <- deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_df %>% dplyr::filter((log2FoldChange < -1))


# log2(1.5) = 2.83
ggmaplot(atacseqkd_bin5kb_de_shCRAMP1_shControl,
         fdr = 0.05, fc = 2, size = 0.3,
         palette = c("#D22B2B", "#1465AC", "darkgray"),
         genenames = rownames(atacseqkd_bin5kb_de_shCRAMP1_shControl),
         legend = "top", top = 20,
         font.label = c("bold", 5),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())

# shSUZ12 vs shControl
atacseqkd_bin5kb_de_shSUZ12_shControl <- results(atacseqkd_bin5kb_ddsO, contrast=c("condition", "shSUZ12", "shControl"))
atacseqkd_bin5kb_de_shSUZ12_shControl = atacseqkd_bin5kb_de_shSUZ12_shControl[order(rownames(atacseqkd_bin5kb_de_shSUZ12_shControl)),]
atacseqkd_bin5kb_de_shSUZ12_shControl$threshold <- as.logical(atacseqkd_bin5kb_de_shSUZ12_shControl$padj < 0.05)
deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05 <- atacseqkd_bin5kb_de_shSUZ12_shControl[which(atacseqkd_bin5kb_de_shSUZ12_shControl$threshold == TRUE),]
deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_df <- data.frame(deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05)
deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_df["interval"] <- rownames(deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_df)
deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_df <- deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_df[,c(8,1:7)]
deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_fc_df <- deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1) | (log2FoldChange < -1))
deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_fc_up <- deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange > 1))
deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_fc_down <- deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_df %>% dplyr::filter((log2FoldChange < -1))


ggmaplot(atacseqkd_bin5kb_de_shSUZ12_shControl,
         fdr = 0.05, fc = 2, size = 0.3,
         palette = c("#D22B2B", "#1465AC", "darkgray"),
         genenames = rownames(atacseqkd_bin5kb_de_shCRAMP1_shControl),
         legend = "top", top = 20,
         font.label = c("bold", 5),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())


# Cut & Run KO


# DBRs analysis using Bins
# bin_5kb_df_filt
cutnrunKO_bin5kb_precountdata <- bin_5kb_df_filt[,c(7:14)]
head(cutnrunKO_bin5kb_precountdata)
dim(cutnrunKO_bin5kb_precountdata)
colSums(cutnrunKO_bin5kb_precountdata)
#write.table(cutnrunKO_bin5kb_precountdata, "cutnrunKO_bin5kb_precountdata.txt", quote = F, append = F, sep="\t")

# Since the bins are not IgG normalized it is important to correct for background 
cutnrunKO_bin5kb_precountdata["H3K27me3_KO1_IgG_KO1"] <- cutnrunKO_bin5kb_precountdata$H3K27me3_KO1 - cutnrunKO_bin5kb_precountdata$IgG_KO1
cutnrunKO_bin5kb_precountdata["H3K27me3_KO2_IgG_KO2"] <- cutnrunKO_bin5kb_precountdata$H3K27me3_KO2 - cutnrunKO_bin5kb_precountdata$IgG_KO2
cutnrunKO_bin5kb_precountdata["H3K27me3_KO3_IgG_KO3"] <- cutnrunKO_bin5kb_precountdata$H3K27me3_KO3 - cutnrunKO_bin5kb_precountdata$IgG_KO3
cutnrunKO_bin5kb_precountdata["H3K27me3_WT_IgG_WT"] <- cutnrunKO_bin5kb_precountdata$H3K27me3_WT - cutnrunKO_bin5kb_precountdata$IgG_WT

# Identify enriched consensus bins
# get regions that are enriched regions in atleast one of the H3K27me3 samples (KO1-3, WT)
# this will give bins like a consensus bins profile where the enrichemnt is observed in atleast one of the sample
# though they are low counts i.e. 1 and above, they are will be filtered on rowsums during DBR analysis
cutnrunKO_bin5kb_precountdata <- cutnrunKO_bin5kb_precountdata[which(cutnrunKO_bin5kb_precountdata$H3K27me3_KO1_IgG_KO1 > 0 
                                                                     | cutnrunKO_bin5kb_precountdata$H3K27me3_KO2_IgG_KO2 > 0 
                                                                     | cutnrunKO_bin5kb_precountdata$H3K27me3_KO3_IgG_KO3 > 0 
                                                                     | cutnrunKO_bin5kb_precountdata$H3K27me3_WT_IgG_WT > 0),]
# Now take only the H3K27me3 counts in bins
cutnrunKO_bin5kb_countdata <- cutnrunKO_bin5kb_precountdata[,c(1:4)]


# create sample-specific coldata object
cutnrunKO_bin5kb_dds_coldata <- data.frame(sample = colnames(cutnrunKO_bin5kb_countdata))

cutnrunKO_bin5kb_dds_coldata["condition"] <- gsub('[0-9]+', '', unlist(lapply(strsplit(cutnrunKO_bin5kb_dds_coldata$sample, "_"), function(x) x[2])))
cutnrunKO_bin5kb_dds_coldata["replicate"] <- c(1,2,3,1)

# rearrange
rownames(cutnrunKO_bin5kb_dds_coldata) <- cutnrunKO_bin5kb_dds_coldata$sample
cutnrunKO_bin5kb_dds_coldata$condition <- factor(cutnrunKO_bin5kb_dds_coldata$condition)
cutnrunKO_bin5kb_dds_coldata$replicate <- factor(cutnrunKO_bin5kb_dds_coldata$replicate)

# prepare count matrix 
cutnrunKO_bin5kb_countdata = as.matrix(cutnrunKO_bin5kb_countdata) 
all(rownames(cutnrunKO_bin5kb_dds_coldata) == colnames(cutnrunKO_bin5kb_countdata)) #should print TRUE
cutnrunKO_bin5kb_ddsO <- DESeqDataSetFromMatrix(countData = cutnrunKO_bin5kb_countdata, colData = cutnrunKO_bin5kb_dds_coldata, design = ~condition)
cutnrunKO_bin5kb_keep <- rowSums(counts(cutnrunKO_bin5kb_ddsO)) > 10
cutnrunKO_bin5kb_ddsO <- cutnrunKO_bin5kb_ddsO[cutnrunKO_bin5kb_keep,]

# pca
cutnrunKO_bin5kb_vst0Norm <- DESeq2::vst(cutnrunKO_bin5kb_ddsO)
cutnrunKO_bin5kb_vst0Norm_pca <- plotPCA(cutnrunKO_bin5kb_vst0Norm, intgroup="condition", returnData=TRUE)
plotPCA(cutnrunKO_bin5kb_vst0Norm, intgroup="condition", returnData=FALSE)


# de analysis
cutnrunKO_bin5kb_ddsO <- DESeq(cutnrunKO_bin5kb_ddsO)

# KO vs WT
cutnrunKO_bin5kb_de_KO_WT <- results(cutnrunKO_bin5kb_ddsO, contrast=c("condition", "KO", "WT"))
cutnrunKO_bin5kb_de_KO_WT = cutnrunKO_bin5kb_de_KO_WT[order(rownames(cutnrunKO_bin5kb_de_KO_WT)),]
cutnrunKO_bin5kb_de_KO_WT$threshold <- as.logical(cutnrunKO_bin5kb_de_KO_WT$padj < 0.05)
deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05 <- cutnrunKO_bin5kb_de_KO_WT[which(cutnrunKO_bin5kb_de_KO_WT$threshold == TRUE),]
deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_df <- data.frame(deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05)
deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_df["interval"] <- rownames(deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_df)
deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_df <- deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_df[,c(8,1:7)]
deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_fc_df <- deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange > 1) | (log2FoldChange < -1))

cutnrunKO_bin5kb_de_KO_WT_df <- data.frame(cutnrunKO_bin5kb_de_KO_WT)
cutnrunKO_bin5kb_de_KO_WT_df["interval"] <- rownames(cutnrunKO_bin5kb_de_KO_WT_df) 

deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_fc.up <- deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange > 1) )
deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_fc.down <- deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_df %>% dplyr::filter((log2FoldChange < -1) )

# log2(1) = 2
ggmaplot(cutnrunKO_bin5kb_de_KO_WT_df,
         fdr = 0.05, fc = 2, size = 0.3,
         palette = c("#D22B2B", "#1465AC", "darkgray"),
         genenames = cutnrunKO_bin5kb_de_KO_WT_df$interval,
         legend = "top", top = 20,
         font.label = c("bold", 5),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())

# Cut & Run KD and Cut & Tag data

# DBRs analysis using Bins
# bin_5kb_df_filt
cutnrunKD_bin5kb_precountdata <- bin_5kb_df_filt[,c(15:20)]
head(cutnrunKD_bin5kb_precountdata)
dim(cutnrunKD_bin5kb_precountdata)
colSums(cutnrunKD_bin5kb_precountdata)
#write.table(cutnrunKD_bin5kb_precountdata, "cutnrunKD_bin5kb_precountdata.txt", quote = F, append = F, sep="\t")

# Since the bins are not IgG normalized it is important to correct for background 
cutnrunKD_bin5kb_precountdata["shControl_H3K27me3_IgG"] <- cutnrunKD_bin5kb_precountdata$shControl_H3K27me3 - cutnrunKD_bin5kb_precountdata$shControl_IgG
cutnrunKD_bin5kb_precountdata["shCRAMP1_H3K27me3_IgG"] <- cutnrunKD_bin5kb_precountdata$shCRAMP1_H3K27me3 - cutnrunKD_bin5kb_precountdata$shCRAMP1_IgG
cutnrunKD_bin5kb_precountdata["shSUZ12_H3K27me3_IgG"] <- cutnrunKD_bin5kb_precountdata$shSUZ12_H3K27me3 - cutnrunKD_bin5kb_precountdata$shSUZ12_IgG

# Identify enriched consensus bins
# get regions that are enriched regions in atleast one of the cutandrun KD samples
# this will give bins like a consensus bins profile where the enrichemnt is observed in atleast one of the sample
# let Cut&Tag data as such
cutnrunKD_bin5kb_precountdata <- cutnrunKD_bin5kb_precountdata[which(cutnrunKD_bin5kb_precountdata$shControl_H3K27me3_IgG > 0 
                                                                     | cutnrunKD_bin5kb_precountdata$shCRAMP1_H3K27me3_IgG > 0 
                                                                     | cutnrunKD_bin5kb_precountdata$shSUZ12_H3K27me3_IgG > 0),]

# For Cut&RUn KD subtract sample-control just to get differential kind of profile
cutnrunKD_bin5kb_precountdata["shCRAMP1_H3K27me3_IgG_shControl"] <- cutnrunKD_bin5kb_precountdata$shCRAMP1_H3K27me3_IgG - cutnrunKD_bin5kb_precountdata$shControl_H3K27me3_IgG
cutnrunKD_bin5kb_precountdata["shSUZ12_H3K27me3_IgG_shControl"] <- cutnrunKD_bin5kb_precountdata$shSUZ12_H3K27me3_IgG - cutnrunKD_bin5kb_precountdata$shControl_H3K27me3_IgG

# Since there is no DBRs analysis here, I just take CPM log subtracted from respective IgG for Cut&Run KD and Cut&Tag data
cutnrunKD_bin5kb_countdata <- cutnrunKD_bin5kb_precountdata[,c(3,5)]

cutnrunKD_bin5kb_countdata_CPM <- edgeR::cpm(cutnrunKD_bin5kb_countdata)
cutnrunKD_bin5kb_countdata_CPM_log <- log(cutnrunKD_bin5kb_countdata_CPM + 1,2)
cutnrunKD_bin5kb_countdata_CPM_log <- data.frame(cutnrunKD_bin5kb_countdata_CPM_log)

# merge
sel_mpgd_genes_akdc_akds <- merge(deseq2_atacseqkd_bin5kb_de_shCRAMP1_shControl_0.05_df, deseq2_atacseqkd_bin5kb_de_shSUZ12_shControl_0.05_df, by.x="interval", by.y="interval", all.x=TRUE, all.y=TRUE)
sel_mpgd_genes_akdc_akds <- sel_mpgd_genes_akdc_akds[,c(1,3,10)]
sel_mpgd_genes_akdc_akds_cko <- merge(sel_mpgd_genes_akdc_akds, deseq2_cutnrunKO_bin5kb_de_KO_WT_0.05_fc_df, by.x="interval", by.y="interval", all.x=TRUE, all.y=TRUE)
sel_mpgd_genes_akdc_akds_cko <- sel_mpgd_genes_akdc_akds_cko[,c(1,2,3,5)]
colnames(sel_mpgd_genes_akdc_akds_cko) <- c("interval","akd_shCRAMP1_shControl","akd_shSUZ12_shControl", "cko_KO_WT")
rownames(sel_mpgd_genes_akdc_akds_cko) <- sel_mpgd_genes_akdc_akds_cko$interval
sel_mpgd_genes_akdc_akds_cko <- sel_mpgd_genes_akdc_akds_cko[,-1]
library(ComplexHeatmap)
sel_mpgd_genes_akdc_akds_cko_mat1 <- Heatmap(as.matrix(sel_mpgd_genes_akdc_akds_cko), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, cluster_rows = FALSE, col = col_fun)


#----------------- Over ATAC-Seq shCRAMP1 DARs peaks ---------------#
#Only DARs 
# create genomic atac_shCRAMP1_dars

# for PE, atac
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job_akdc ./quantify_over_atac_shCRAMP1_dars_pe.sh

# for all SE, cutnrun KO,KD and cutntag
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job_akdc ./quantify_over_atac_shCRAMP1_dars_se.sh
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job_akdc ./quantify_over_atac_shCRAMP1_dars_se.sh
/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job_akdc ./quantify_over_atac_shCRAMP1_dars_se.sh


atac_shCRAMP1_dars_cov_atacseq_path <- list.files(path="/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library", pattern = "*akdc_coverage.pe.bed", full.names = T)
atac_shCRAMP1_dars_cov_cutnrun_KO_path <- list.files(path="/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary", pattern = "*akdc_coverage_se.bed", full.names = T)
atac_shCRAMP1_dars_cov_cutnrun_KD_path <- list.files(path="/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/bowtie2/mergedLibrary", pattern = "*akdc_coverage_se.bed", full.names = T)
atac_shCRAMP1_dars_cov_cutntag_path <- list.files(path="/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary", pattern = "*akdc_coverage_se.bed", full.names = T)
atac_shCRAMP1_dars_path_combined <- c(atac_shCRAMP1_dars_cov_atacseq_path, atac_shCRAMP1_dars_cov_cutnrun_KO_path, atac_shCRAMP1_dars_cov_cutnrun_KD_path, atac_shCRAMP1_dars_cov_cutntag_path)

atac_shCRAMP1_dars_list <- list()
for (files in atac_shCRAMP1_dars_path_combined){
  filename <- sub("\\..*", "", basename(files))
  print(filename)
  atac_shCRAMP1_dars_list[[filename]] <- data.frame(fread(files, header = F))
}
atac_shCRAMP1_dars_df <- do.call(cbind.data.frame, atac_shCRAMP1_dars_list)
atac_shCRAMP1_dars_df <- atac_shCRAMP1_dars_df[,c(1:4, seq(5,200,8))]
rownames(atac_shCRAMP1_dars_df) <- paste0(atac_shCRAMP1_dars_df$shControl_REP1.V1, "%", atac_shCRAMP1_dars_df$shControl_REP1.V2, "%", atac_shCRAMP1_dars_df$shControl_REP1.V3, "%", atac_shCRAMP1_dars_df$shControl_REP1.V4)
atac_shCRAMP1_dars_df <- atac_shCRAMP1_dars_df[,-c(1:4)]
colnames(atac_shCRAMP1_dars_df) <- gsub(".V5","", colnames(atac_shCRAMP1_dars_df))
atac_shCRAMP1_dars_df_filt <- atac_shCRAMP1_dars_df[rowSums(atac_shCRAMP1_dars_df) > 0,]
atac_shCRAMP1_dars_df_CPM <- edgeR::cpm(atac_shCRAMP1_dars_df_filt)
atac_shCRAMP1_dars_df_CPM_log <- log(atac_shCRAMP1_dars_df_CPM + 1,2)
atac_shCRAMP1_dars_df_CPM_log <- data.frame(atac_shCRAMP1_dars_df_CPM_log)
atac_shCRAMP1_dars_df_CPM_log["akd_shCRAMP1_shControl_1"] <- atac_shCRAMP1_dars_df_CPM_log$shCRAMP1_REP1 - atac_shCRAMP1_dars_df_CPM_log$shControl_REP1
atac_shCRAMP1_dars_df_CPM_log["akd_shCRAMP1_shControl_2"] <- atac_shCRAMP1_dars_df_CPM_log$shCRAMP1_REP2 - atac_shCRAMP1_dars_df_CPM_log$shControl_REP2
atac_shCRAMP1_dars_df_CPM_log["akd_shSUZ12_shControl_1"] <- atac_shCRAMP1_dars_df_CPM_log$shSUZ12_REP1 - atac_shCRAMP1_dars_df_CPM_log$shControl_REP1
atac_shCRAMP1_dars_df_CPM_log["akd_shSUZ12_shControl_2"] <- atac_shCRAMP1_dars_df_CPM_log$shSUZ12_REP2 - atac_shCRAMP1_dars_df_CPM_log$shControl_REP2
atac_shCRAMP1_dars_df_CPM_log["cko_H3K27me3_KO1_IgG"] <- atac_shCRAMP1_dars_df_CPM_log$H3K27me3_KO1 - atac_shCRAMP1_dars_df_CPM_log$IgG_KO1
atac_shCRAMP1_dars_df_CPM_log["cko_H3K27me3_KO2_IgG"] <- atac_shCRAMP1_dars_df_CPM_log$H3K27me3_KO2 - atac_shCRAMP1_dars_df_CPM_log$IgG_KO2
atac_shCRAMP1_dars_df_CPM_log["cko_H3K27me3_KO3_IgG"] <- atac_shCRAMP1_dars_df_CPM_log$H3K27me3_KO3 - atac_shCRAMP1_dars_df_CPM_log$IgG_KO3
atac_shCRAMP1_dars_df_CPM_log["cko_H3K27me3_WT_IgG"] <- atac_shCRAMP1_dars_df_CPM_log$H3K27me3_WT - atac_shCRAMP1_dars_df_CPM_log$IgG_WT
atac_shCRAMP1_dars_df_CPM_log["ckd_shControl_H3K27me3_IgG"] <- atac_shCRAMP1_dars_df_CPM_log$shControl_H3K27me3 - atac_shCRAMP1_dars_df_CPM_log$shControl_IgG
atac_shCRAMP1_dars_df_CPM_log["ckd_shCRAMP1_H3K27me3_IgG"] <- atac_shCRAMP1_dars_df_CPM_log$shCRAMP1_H3K27me3 - atac_shCRAMP1_dars_df_CPM_log$shCRAMP1_IgG
atac_shCRAMP1_dars_df_CPM_log["ckd_shSUZ12_H3K27me3_IgG"] <- atac_shCRAMP1_dars_df_CPM_log$shSUZ12_H3K27me3 - atac_shCRAMP1_dars_df_CPM_log$shSUZ12_IgG
atac_shCRAMP1_dars_df_CPM_log["ct_ENCFF508LLH_IgG"] <- atac_shCRAMP1_dars_df_CPM_log$ENCFF508LLH - atac_shCRAMP1_dars_df_CPM_log$IgG
atac_shCRAMP1_dars_df_CPM_log["ct_H1_4_IgG"] <- atac_shCRAMP1_dars_df_CPM_log$H1_4 - atac_shCRAMP1_dars_df_CPM_log$IgG
atac_shCRAMP1_dars_df_CPM_log["ct_H3K27me3_IgG"] <- atac_shCRAMP1_dars_df_CPM_log$H3K27me3 - atac_shCRAMP1_dars_df_CPM_log$IgG
atac_shCRAMP1_dars_df_CPM_log["ct_panH1_IgG"] <- atac_shCRAMP1_dars_df_CPM_log$panH1 - atac_shCRAMP1_dars_df_CPM_log$IgG

atac_shCRAMP1_dars_df_CPM_log["cko_H3K27me3_KO1_IgG_WT"] <- atac_shCRAMP1_dars_df_CPM_log$cko_H3K27me3_KO1_IgG - atac_shCRAMP1_dars_df_CPM_log$cko_H3K27me3_WT_IgG
atac_shCRAMP1_dars_df_CPM_log["cko_H3K27me3_KO2_IgG_WT"] <- atac_shCRAMP1_dars_df_CPM_log$cko_H3K27me3_KO2_IgG - atac_shCRAMP1_dars_df_CPM_log$cko_H3K27me3_WT_IgG
atac_shCRAMP1_dars_df_CPM_log["cko_H3K27me3_KO3_IgG_WT"] <- atac_shCRAMP1_dars_df_CPM_log$cko_H3K27me3_KO3_IgG - atac_shCRAMP1_dars_df_CPM_log$cko_H3K27me3_WT_IgG

atac_shCRAMP1_dars_df_CPM_log["ckd_shCRAMP1_H3K27me3_IgG_shControl"] <- atac_shCRAMP1_dars_df_CPM_log$ckd_shCRAMP1_H3K27me3_IgG - atac_shCRAMP1_dars_df_CPM_log$ckd_shControl_H3K27me3_IgG
atac_shCRAMP1_dars_df_CPM_log["ckd_shSUZ12_H3K27me3_IgG_shControl"] <- atac_shCRAMP1_dars_df_CPM_log$ckd_shSUZ12_H3K27me3_IgG - atac_shCRAMP1_dars_df_CPM_log$ckd_shControl_H3K27me3_IgG

atac_shCRAMP1_dars_df_CPM_log_sub <- atac_shCRAMP1_dars_df_CPM_log[,c(26:29,41:45,37:40)]
dim(atac_shCRAMP1_dars_df_CPM_log_sub)


Heatmap(as.matrix(atac_shCRAMP1_dars_df_CPM_log_sub), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_columns = FALSE, use_raster = TRUE)

# atac_shCRAMP1_dars_df_CPM_log_sub_sorted <- atac_shCRAMP1_dars_df_CPM_log_sub[order(-atac_shCRAMP1_dars_df_CPM_log_sub$akd_shCRAMP1_shControl_1),]
atac_shCRAMP1_dars_df_CPM_log_sub_mat1 <- Heatmap(as.matrix(atac_shCRAMP1_dars_df_CPM_log_sub[,c(1:4)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_rows = TRUE, cluster_columns = FALSE, use_raster = TRUE, row_km = 3)
atac_shCRAMP1_dars_df_CPM_log_sub_mat2 <- Heatmap(as.matrix(atac_shCRAMP1_dars_df_CPM_log_sub[,c(5:9)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_rows = TRUE, cluster_columns = FALSE, use_raster = TRUE)
atac_shCRAMP1_dars_df_CPM_log_sub_mat3 <- Heatmap(as.matrix(atac_shCRAMP1_dars_df_CPM_log_sub[,c(10:13)]), na_col = "yellow", row_names_gp = gpar(fontsize = 5), cluster_rows = TRUE, cluster_columns = FALSE, use_raster = TRUE)

atac_shCRAMP1_dars_df_CPM_log_sub_mat1 + atac_shCRAMP1_dars_df_CPM_log_sub_mat2 + atac_shCRAMP1_dars_df_CPM_log_sub_mat3
 # notes H1, TRB1, PRC2 activity
library(M3C)
# do PCA
pca(atac_shCRAMP1_dars_df_CPM_log[,c(26:29,41:45,37:40)],colvec=c('gold'),printres=TRUE, text=colnames(atac_shCRAMP1_dars_df_CPM_log[,c(26:29,41:45,37:40)]))


# Over regions extended from gene coordinates  (1Kb , 2kb, 3Kb)
#Annotation to regions , Note only promoter need to be taken not gene body
# Over genes extended with promoter (1500bp)
#Annotation to regions
gene_v41_human_out <- data.frame(fread("/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/gene_v41_human_out.txt", header = F))
#1kb
gene_v41_human_TSS1kb <- gene_v41_human_out
# extend TSS code
gene_v41_human_TSS1kb_plus <- gene_v41_human_TSS1kb[which(gene_v41_human_TSS1kb$V4 == "+"),]

gene_v41_human_TSS1kb_plus$V2 <- gene_v41_human_TSS1kb_plus$V2 - 1000

gene_v41_human_TSS1kb_minus <- gene_v41_human_TSS1kb[which(gene_v41_human_TSS1kb$V4 == "-"),]

gene_v41_human_TSS1kb_minus$V3 <- gene_v41_human_TSS1kb_minus$V3 + 1000

# once written should write again untill run the above extend TSS code again, becuase re-subtraction should be avoided 
write.table(gene_v41_human_TSS1kb_plus, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/gene_v41_human_TSS1kb_plus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)
write.table(gene_v41_human_TSS1kb_minus, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/gene_v41_human_TSS1kb_minus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)

cat gene_v41_human_TSS1kb_plus.txt gene_v41_human_TSS1kb_minus.txt | sort -k1,1 -k2,2n > gene_v41_human_TSS1kb_plus_minus.txt
grep "chrM" gene_v41_human_TSS1kb_plus_minus.txt -v > gene_v41_human_TSS1kb_plus_minus_freeM.txt


#2kb
gene_v41_human_TSS2kb <- gene_v41_human_out
# extend TSS code
gene_v41_human_TSS2kb_plus <- gene_v41_human_TSS2kb[which(gene_v41_human_TSS2kb$V4 == "+"),]

gene_v41_human_TSS2kb_plus$V2 <- gene_v41_human_TSS2kb_plus$V2 - 2000

gene_v41_human_TSS2kb_minus <- gene_v41_human_TSS2kb[which(gene_v41_human_TSS2kb$V4 == "-"),]

gene_v41_human_TSS2kb_minus$V3 <- gene_v41_human_TSS2kb_minus$V3 + 2000

# once written should write again untill run the above extend TSS code again, becuase re-subtraction should be avoided 
write.table(gene_v41_human_TSS2kb_plus, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/gene_v41_human_TSS2kb_plus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)
write.table(gene_v41_human_TSS2kb_minus, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/gene_v41_human_TSS2kb_minus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)

cat gene_v41_human_TSS2kb_plus.txt gene_v41_human_TSS2kb_minus.txt | sort -k1,1 -k2,2n > gene_v41_human_TSS2kb_plus_minus.txt
grep "chrM" gene_v41_human_TSS2kb_plus_minus.txt -v > gene_v41_human_TSS2kb_plus_minus_freeM.txt

#3kb

gene_v41_human_TSS3kb <- gene_v41_human_out
# extend TSS code
gene_v41_human_TSS3kb_plus <- gene_v41_human_TSS3kb[which(gene_v41_human_TSS3kb$V4 == "+"),]

gene_v41_human_TSS3kb_plus$V2 <- gene_v41_human_TSS3kb_plus$V2 - 3000

gene_v41_human_TSS3kb_minus <- gene_v41_human_TSS3kb[which(gene_v41_human_TSS3kb$V4 == "-"),]

gene_v41_human_TSS3kb_minus$V3 <- gene_v41_human_TSS3kb_minus$V3 + 3000

# once written should write again untill run the above extend TSS code again, becuase re-subtraction should be avoided 
write.table(gene_v41_human_TSS3kb_plus, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/gene_v41_human_TSS3kb_plus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)
write.table(gene_v41_human_TSS3kb_minus, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/gene_v41_human_TSS3kb_minus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)

cat gene_v41_human_TSS3kb_plus.txt gene_v41_human_TSS3kb_minus.txt | sort -k1,1 -k2,2n > gene_v41_human_TSS3kb_plus_minus.txt
grep "chrM" gene_v41_human_TSS3kb_plus_minus.txt -v > gene_v41_human_TSS3kb_plus_minus_freeM.txt

--------

# rename as the names are too big
mv gene_v41_human_TSS3kb_plus_minus_freeM.txt gene_3kb_freeM.txt
mv gene_v41_human_TSS2kb_plus_minus_freeM.txt gene_2kb_freeM.txt
mv gene_v41_human_TSS1kb_plus_minus_freeM.txt gene_1kb_freeM.txt


/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job_gu ./quantify_over_genes_upstream_regions.sh

# for vars in *.bam
# do
# echo "${vars}"
# for genes_upstream in gene*kb*freeM.txt
# do
# echo "${genes_upstream}"
# bedtools coverage -a ${genes_upstream} -b ${vars}_chr.bed | sort -k1,1 -k2,2n > ${vars}_${genes_upstream}_coverage_se.bed
# done
# done


cutnrunKO_genesupstream_path <- "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/"
cutnrunKO_genesupstream_cov_files <- list.files(cutnrunKO_genesupstream_path, pattern = "_freeM.txt_coverage_se.bed")
cutnrunKO_genesupstream_coverage_list <- list()

for (gene_cov_file in cutnrunKO_genesupstream_cov_files){
  print(gene_cov_file)
  gene_cov_file_tag <- gsub("_freeM.txt_coverage_se.bed","",gsub(".mLb.clN.sorted.bam_gene","", gene_cov_file))
  print(gene_cov_file_tag)
  cutnrunKO_genesupstream_coverage_list[[gene_cov_file_tag]] <- data.frame(fread(paste0(cutnrunKO_genesupstream_path, gene_cov_file)))
}

cutnrunKO_genesupstream_coverage_df <- do.call(cbind.data.frame, cutnrunKO_genesupstream_coverage_list)

cutnrunKO_genesupstream_precountdata <- cutnrunKO_genesupstream_coverage_df[,c(6,7, seq(8,dim(cutnrunKO_genesupstream_coverage_df)[2],11))]
colnames(cutnrunKO_genesupstream_precountdata) <- c("ensID","Gene",colnames(cutnrunKO_genesupstream_precountdata)[3:dim(cutnrunKO_genesupstream_precountdata)[2]])
rownames(cutnrunKO_genesupstream_precountdata) <- paste0(cutnrunKO_genesupstream_precountdata$ensID, "%",cutnrunKO_genesupstream_precountdata$Gene)
cutnrunKO_genesupstream_precountdata <- cutnrunKO_genesupstream_precountdata[,-c(1,2)]
colnames(cutnrunKO_genesupstream_precountdata) <- gsub(".V8","",colnames(cutnrunKO_genesupstream_precountdata))


# use the default cpm function becuase it is a gene-count matrix so for ~60k gene we cannot expect the bias will occur fro any given sample and we do not need give separate lib size.
cutnrunKO_genesupstream_precountdata_myCPM <- data.frame(edgeR::cpm(cutnrunKO_genesupstream_precountdata))
cutnrunKO_genesupstream_precountdata_myCPM_log <- log(cutnrunKO_genesupstream_precountdata_myCPM + 1, 2)

# subtract sample - IgG

# Assuming your data is stored in a data frame named 'df'

library(dplyr)

# Columns to subtract
cutnrunKO_h3k27me3_columns <- colnames(cutnrunKO_genesupstream_precountdata_myCPM_log)[1:12]

cutnrunKO_igg_columns <- colnames(cutnrunKO_genesupstream_precountdata_myCPM_log)[13:24]


# Subtract H3K27me3 from IgG for each corresponding condition
for (i in seq_along(cutnrunKO_h3k27me3_columns)) {
  print(paste0(cutnrunKO_h3k27me3_columns[i], "  ",cutnrunKO_igg_columns[i]))
  cutnrunKO_genesupstream_precountdata_myCPM_log[paste0("diff_",cutnrunKO_h3k27me3_columns[i], "_",cutnrunKO_igg_columns[i])] <- cutnrunKO_genesupstream_precountdata_myCPM_log[cutnrunKO_h3k27me3_columns[i]] - cutnrunKO_genesupstream_precountdata_myCPM_log[cutnrunKO_igg_columns[i]]
}

# select only diff columns 
cutnrunKO_genesupstream_precountdata_myCPM_log_diff <- cutnrunKO_genesupstream_precountdata_myCPM_log %>% dplyr::select(contains('diff'))
cutnrunKO_genesupstream_precountdata_myCPM_log_diff_st <- data.frame(stack(as.matrix(cutnrunKO_genesupstream_precountdata_myCPM_log_diff)))
cutnrunKO_genesupstream_precountdata_myCPM_log_diff_st$col <- as.character(cutnrunKO_genesupstream_precountdata_myCPM_log_diff_st$col)
cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set <- merge(cutnrunKO_genesupstream_precountdata_myCPM_log_diff_st, set_rnaseqde_merged_shC1_shS_up, by.x="row", by.y="Gene", all.x=F, all.y=T)
cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set["id"] <- paste0(unlist(lapply(strsplit(cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set$col,"_"), function(x) paste0(x[2:4],collapse = "_"))))

cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_1kb <- cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set %>% dplyr::filter(grepl('1kb', col))
cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb <- cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set %>% dplyr::filter(grepl('2kb', col))
cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_3kb <- cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set %>% dplyr::filter(grepl('3kb', col))

ggplot(cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb, aes(x=value, color=label, linetype=id)) +
  geom_density()+
  scale_color_manual(values=c("#00D100", "black","orange","blue"))+
  scale_linetype_manual(values = c(rep("dotted",1), rep("longdash",1),rep("dashed",1),rep("solid",1))) + theme_classic() + xlim(c(-4,5)) + ylim(c(0,1)) +
  geom_vline(aes(xintercept=0),color="grey", linetype="dashed", size=0.2)


cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_WT <- cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb %>% dplyr::filter(grepl('WT', col))

plot_cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_WT <- ggplot(cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_WT, aes(x=value, color=label, linetype=id)) +
  geom_density()+
  scale_color_manual(values=c("#00D100", "black","orange","blue"))+
  scale_linetype_manual(values = c(rep("solid",1))) + theme_classic() + xlim(c(-4,5)) + ylim(c(0,1)) +
  geom_vline(aes(xintercept=0),color="grey", linetype="dashed", size=0.2)

cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_KO <- cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb %>% dplyr::filter(grepl('KO', col))
plot_cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_KO <- ggplot(cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_KO, aes(x=value, color=label, linetype=id)) +
  geom_density()+
  scale_color_manual(values=c("#00D100", "black","orange","blue"))+
  scale_linetype_manual(values = c(rep("longdash",1),rep("dashed",1),rep("dotted",1))) + theme_classic() + xlim(c(-4,5)) + ylim(c(0,1)) +
  geom_vline(aes(xintercept=0),color="grey", linetype="dashed", size=0.2)

plot_cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_WT + plot_cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_KO + plot_layout(ncol = 2)

# cutnrunKO and rnase DE genes overlap
rnaseqde_cutnrunKO_de_anno_genes_list <- list(cko_up=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.up_anno_nearest$ensID_Gene),
                                              cko_down=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.down_anno_nearest$ensID_Gene),
                                              rkd_shC1_up = unique(set_rnaseqde_uniq_shC1_up$Gene),
                                              rkd_shS_up = unique(set_rnaseqde_uniq_shS_up$Gene),
                                              rkd_shC1_shS_up = unique(set_rnaseqde_both_shC1_shS_up$Gene),
                                              rkd_none_up = unique(set_rnaseqde_none_up$Gene))

upset(fromList(rnaseqde_cutnrunKO_de_anno_genes_list), order.by = "freq",empty.intersections = "on")

rnaseqde_cutnrunKO_de_anno_genes_list2 <- list(cko_up=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.up_anno_nearest$ensID_Gene),
                                              cko_down=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.down_anno_nearest$ensID_Gene),
                                              rkd_shC1_up = unique(set_rnaseqde_uniq_shC1_up$Gene),
                                              rkd_shS_up = unique(set_rnaseqde_uniq_shS_up$Gene),
                                              rkd_shC1_shS_up = unique(set_rnaseqde_both_shC1_shS_up$Gene))
venn(rnaseqde_cutnrunKO_de_anno_genes_list2)

rnaseqde_cutnrunKO_non_de_anno_genes_list3 <- list(cko_all=unique(cutnrunKOde_KO_WT_re_anno_nearest$ensID_Gene),
                                                   cko_up=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.up_anno_nearest$ensID_Gene),
                                                   cko_down=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.down_anno_nearest$ensID_Gene),
                                                   rkd_shC1_up = unique(set_rnaseqde_uniq_shC1_up$Gene),
                                                   rkd_shS_up = unique(set_rnaseqde_uniq_shS_up$Gene),
                                                   rkd_shC1_shS_up = unique(set_rnaseqde_both_shC1_shS_up$Gene),
                                                   rkd_none_up = unique(set_rnaseqde_none_up$Gene))
upset(fromList(rnaseqde_cutnrunKO_non_de_anno_genes_list3), order.by = "freq",empty.intersections = "on")

rnaseqde_cutnrunKO_non_de_anno_genes_list4 <- list(cko_all=unique(cutnrunKOde_KO_WT_re_anno_nearest$ensID_Gene),
                                                   cko_up=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.up_anno_nearest$ensID_Gene),
                                                   cko_down=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.down_anno_nearest$ensID_Gene),
                                                   rkd_shC1_shS_up = unique(set_rnaseqde_both_shC1_shS_up$Gene))
venn(rnaseqde_cutnrunKO_non_de_anno_genes_list4)
upset(fromList(rnaseqde_cutnrunKO_non_de_anno_genes_list4), order.by = "freq",empty.intersections = "on")


cutnrunKO_de_anno_genes_list1 <- list(cko_all=unique(cutnrunKOde_KO_WT_re_anno_nearest$ensID_Gene),
                                                   cko_up=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.up_anno_nearest$ensID_Gene),
                                                   cko_down=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.down_anno_nearest$ensID_Gene))

rnaseqde_none_anno_genes_list1 <- list(rkd_shC1_up = unique(set_rnaseqde_uniq_shC1_up$Gene),
                                  rkd_shS_up = unique(set_rnaseqde_uniq_shS_up$Gene),
                                  rkd_shC1_shS_up = unique(set_rnaseqde_both_shC1_shS_up$Gene),
                                  rkd_none_up = unique(set_rnaseqde_none_up$Gene))

rnaseqde_cutnrunKO_non_de_anno_genes_int_matrix <- matrix()
for (i in seq_along(cutnrunKO_de_anno_genes_list1)){
  print(length(cutnrunKO_de_anno_genes_list1[[i]]))
  for (j in seq_along(rnaseqde_none_anno_genes_list1)){
    print(length(rnaseqde_none_anno_genes_list1[[j]]))
    rnaseqde_cutnrunKO_non_de_anno_genes_int_matrix[paste0(names(cutnrunKO_de_anno_genes_list1)[[i]],"_",names(rnaseqde_none_anno_genes_list1)[[j]])] <- length(intersect(cutnrunKO_de_anno_genes_list1[[i]], rnaseqde_none_anno_genes_list1[[j]]))
  }
}

for (i in seq_along(rnaseqde_none_anno_genes_list1)){
  print(length(rnaseqde_none_anno_genes_list1[[i]]))
  rnaseqde_cutnrunKO_non_de_anno_genes_int_matrix[paste0("Total","_",names(rnaseqde_none_anno_genes_list1)[[i]])] <- length(rnaseqde_none_anno_genes_list1[[i]])
}

rnaseqde_cutnrunKO_non_de_anno_genes_int_df <- data.frame(Overlap=rnaseqde_cutnrunKO_non_de_anno_genes_int_matrix)
rnaseqde_cutnrunKO_non_de_anno_genes_int_df <- na.omit(rnaseqde_cutnrunKO_non_de_anno_genes_int_df)
rnaseqde_cutnrunKO_non_de_anno_genes_int_df["cko_datatype"] <- unlist(lapply(strsplit(rownames(rnaseqde_cutnrunKO_non_de_anno_genes_int_df), "_rkd"), function(x) x[1]))
rnaseqde_cutnrunKO_non_de_anno_genes_int_df["rkd_datatype"] <- unlist(lapply(strsplit(rownames(rnaseqde_cutnrunKO_non_de_anno_genes_int_df), "_rkd"), function(x) paste0("rkd",x[2])))
rnaseqde_cutnrunKO_non_de_anno_genes_int_df$cko_datatype <- factor(rnaseqde_cutnrunKO_non_de_anno_genes_int_df$cko_datatype, levels = c("Total", "cko_all", "cko_up", "cko_down"))
library(ggbreak)
ggplot(rnaseqde_cutnrunKO_non_de_anno_genes_int_df, aes(x=rkd_datatype, y=Overlap, fill=as.factor(cko_datatype))) + 
  geom_bar(stat = "identity",position=position_dodge()) +
  scale_fill_manual(values = c("#999999","#3CB371","#E69F00", "#56B4E9") )+ theme_minimal()


# # Further subset to genes near differential binding regions
# deseq2_cutnrunKOde_KO_WT_0.05_fc.up_anno_nearest_genes <- data.frame(Gene=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.up_anno_nearest$ensID_Gene), DBR="UP")
# deseq2_cutnrunKOde_KO_WT_0.05_fc.down_anno_nearest_genes <- data.frame(Gene=unique(deseq2_cutnrunKOde_KO_WT_0.05_fc.down_anno_nearest$ensID_Gene), DBR="DOWN")
# deseq2_cutnrunKOde_KO_WT_0.05_fc.up_down_anno_nearest_genes <- rbind.data.frame(deseq2_cutnrunKOde_KO_WT_0.05_fc.up_anno_nearest_genes, deseq2_cutnrunKOde_KO_WT_0.05_fc.down_anno_nearest_genes)
# 
# 
# cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_DE <- merge(cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb, deseq2_cutnrunKOde_KO_WT_0.05_fc.up_down_anno_nearest_genes, by.x="row", by.y="Gene", all.x=T,all.y=T)
# cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_DE["label_DBR"] <- paste0(cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_DE$label, "_", cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_DE$DBR)
# 
# 
# cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_DE_WT <- cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_DE %>% dplyr::filter(grepl('WT', col))
# 
# plot_cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_DE_WT <- ggplot(cutnrunKO_genesupstream_precountdata_myCPM_log_diff_set_2kb_DE_WT, aes(x=value, color=label_DBR, linetype=id)) +
#   geom_density()+
#   scale_color_manual(values=c(rep(c("#08A045", "#6BBF59", "#21D375"),1), rep(c("grey","black","#B0AA7E"),1),rep(c("#FFC8C4", "#FFB0C4", "#FF97C4"),1),rep(c("#5603AD","blue","cyan"),1)))+
#   scale_linetype_manual(values = c(rep("solid",1))) + theme_classic() + xlim(c(-4,5)) + ylim(c(0,1.2)) +
#   geom_vline(aes(xintercept=0),color="grey", linetype="dashed", size=0.2)
# 

# for cutntag gene upstream
cp -r /mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_ko/outfolder/bowtie2/mergedLibrary/gene_*_freeM.txt .

cutntag_genesupstream_path <- "/mnt/home3/reid/av638/cutntag/iva_lab_oct23/outfolder/bowtie2/mergedLibrary/"
cutntag_genesupstream_cov_files <- list.files(cutntag_genesupstream_path, pattern = "_freeM.txt_coverage_se.bed")
cutntag_genesupstream_coverage_list <- list()

for (gene_cov_file in cutntag_genesupstream_cov_files){
  print(gene_cov_file)
  gene_cov_file_tag <- gsub("_freeM.txt_coverage_se.bed","",gsub(".mLb.clN.sorted.bam_gene","", gene_cov_file))
  print(gene_cov_file_tag)
  cutntag_genesupstream_coverage_list[[gene_cov_file_tag]] <- data.frame(fread(paste0(cutntag_genesupstream_path, gene_cov_file)))
}

cutntag_genesupstream_coverage_df <- do.call(cbind.data.frame, cutntag_genesupstream_coverage_list)

cutntag_genesupstream_precountdata <- cutntag_genesupstream_coverage_df[,c(6,7, seq(8,dim(cutntag_genesupstream_coverage_df)[2],11))]
colnames(cutntag_genesupstream_precountdata) <- c("ensID","Gene",colnames(cutntag_genesupstream_precountdata)[3:dim(cutntag_genesupstream_precountdata)[2]])
rownames(cutntag_genesupstream_precountdata) <- paste0(cutntag_genesupstream_precountdata$ensID, "%",cutntag_genesupstream_precountdata$Gene)
cutntag_genesupstream_precountdata <- cutntag_genesupstream_precountdata[,-c(1,2)]
colnames(cutntag_genesupstream_precountdata) <- gsub(".V8","",colnames(cutntag_genesupstream_precountdata))


# use the default cpm function becuase it is a gene-count matrix so for ~60k gene we cannot expect the bias will occur fro any given sample and we do not need give separate lib size.
cutntag_genesupstream_precountdata_myCPM <- data.frame(edgeR::cpm(cutntag_genesupstream_precountdata))
cutntag_genesupstream_precountdata_myCPM_log <- log(cutntag_genesupstream_precountdata_myCPM + 1, 2)

# subtract sample - IgG

# Assuming your data is stored in a data frame named 'df'

library(dplyr)

# Columns to subtract
cutntag_cutntag_columns <- colnames(cutntag_genesupstream_precountdata_myCPM_log)[c(1:9,13:15)]

cutntag_igg_columns <- rep(colnames(cutntag_genesupstream_precountdata_myCPM_log)[10:12],4)


# Subtract H3K27me3 from IgG for each corresponding condition
for (i in seq_along(cutntag_cutntag_columns)) {
  print(paste0(cutntag_cutntag_columns[i], "  ",cutntag_igg_columns[i]))
  cutntag_genesupstream_precountdata_myCPM_log[paste0("diff_",cutntag_cutntag_columns[i], "_",cutntag_igg_columns[i])] <- cutntag_genesupstream_precountdata_myCPM_log[cutntag_cutntag_columns[i]] - cutntag_genesupstream_precountdata_myCPM_log[cutntag_igg_columns[i]]
}

# select only diff columns 
cutntag_genesupstream_precountdata_myCPM_log_diff <- cutntag_genesupstream_precountdata_myCPM_log %>% dplyr::select(contains('diff'))
cutntag_genesupstream_precountdata_myCPM_log_diff_st <- data.frame(stack(as.matrix(cutntag_genesupstream_precountdata_myCPM_log_diff)))
cutntag_genesupstream_precountdata_myCPM_log_diff_st$col <- as.character(cutntag_genesupstream_precountdata_myCPM_log_diff_st$col)
cutntag_genesupstream_precountdata_myCPM_log_diff_set <- merge(cutntag_genesupstream_precountdata_myCPM_log_diff_st, set_rnaseqde_merged_shC1_shS_up, by.x="row", by.y="Gene", all.x=F, all.y=T)
cutntag_genesupstream_precountdata_myCPM_log_diff_set["id"] <- paste0(unlist(lapply(strsplit(cutntag_genesupstream_precountdata_myCPM_log_diff_set$col,"_"), function(x) paste0(x[2:4],collapse = "_"))))

cutntag_genesupstream_precountdata_myCPM_log_diff_set_1kb <- cutntag_genesupstream_precountdata_myCPM_log_diff_set %>% dplyr::filter(grepl('1kb', col))
cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb <- cutntag_genesupstream_precountdata_myCPM_log_diff_set %>% dplyr::filter(grepl('2kb', col))
cutntag_genesupstream_precountdata_myCPM_log_diff_set_3kb <- cutntag_genesupstream_precountdata_myCPM_log_diff_set %>% dplyr::filter(grepl('3kb', col))

cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_H1_4 <- cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb %>% dplyr::filter(grepl('H1_4', id))


plot_cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_H1_4 <- ggplot(cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_H1_4, aes(x=value, color=label, linetype=id)) +
  geom_density()+
  scale_color_manual(values=c("#00D100", "black","orange","blue"))+
  scale_linetype_manual(values = c(rep("solid",1),rep("longdash",1),rep("solid",1),rep("dotted",1))) + theme_classic() + xlim(c(-8,8)) + ylim(c(0,1)) + 
  geom_vline(aes(xintercept=0),color="grey", linetype="dashed", size=0.2)


cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_panH1 <- cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb %>% dplyr::filter(grepl('panH1', id))


plot_cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_panH1 <- ggplot(cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_panH1, aes(x=value, color=label, linetype=id)) +
  geom_density()+
  scale_color_manual(values=c("#00D100", "black","orange","blue"))+
  scale_linetype_manual(values = c(rep("solid",1),rep("longdash",1),rep("solid",1),rep("dotted",1))) + theme_classic() + xlim(c(-8,8)) + ylim(c(0,1)) + 
  geom_vline(aes(xintercept=0),color="grey", linetype="dashed", size=0.2)

cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_K27me3 <- cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb %>% dplyr::filter(grepl('K27me3', id))


plot_cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_K27me3 <- ggplot(cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_K27me3, aes(x=value, color=label, linetype=id)) +
  geom_density()+
  scale_color_manual(values=c("#00D100", "black","orange","blue"))+
  scale_linetype_manual(values = c(rep("solid",1),rep("longdash",1),rep("solid",1),rep("dotted",1))) + theme_classic() + xlim(c(-8,8)) + ylim(c(0,1)) + 
  geom_vline(aes(xintercept=0),color="grey", linetype="dashed", size=0.2)
plot_cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_panH1 + 
  plot_cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_H1_4 +
  plot_cutntag_genesupstream_precountdata_myCPM_log_diff_set_2kb_K27me3 + plot_layout(ncol = 3)
# for atacseq dars
# Over regions extended from gene coordinates  (1Kb , 2kb, 3Kb)
#Annotation to regions , Note only promoter need to be taken not gene body
# Over genes extended with promoter (1500bp)
#Annotation to regions
gene_v41_human_out <- data.frame(fread("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/gene_v41_human_out.txt", header = F))
#1kb
gene_v41_human_onlyTSS1kb <- gene_v41_human_out
# extend TSS code
gene_v41_human_onlyTSS1kb_plus <- gene_v41_human_onlyTSS1kb[which(gene_v41_human_onlyTSS1kb$V4 == "+"),]
# replace column 3 by column 2
gene_v41_human_onlyTSS1kb_plus$V3 <- gene_v41_human_onlyTSS1kb_plus$V2
gene_v41_human_onlyTSS1kb_plus$V2 <- gene_v41_human_onlyTSS1kb_plus$V2 - 1000

gene_v41_human_onlyTSS1kb_minus <- gene_v41_human_onlyTSS1kb[which(gene_v41_human_onlyTSS1kb$V4 == "-"),]
# replace column 2 by column 3
gene_v41_human_onlyTSS1kb_minus$V2 <- gene_v41_human_onlyTSS1kb_minus$V3
gene_v41_human_onlyTSS1kb_minus$V3 <- gene_v41_human_onlyTSS1kb_minus$V3 + 1000

# once written should write again untill run the above extend TSS code again, becuase re-subtraction should be avoided 
write.table(gene_v41_human_onlyTSS1kb_plus, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/gene_v41_human_onlyTSS1kb_plus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)
write.table(gene_v41_human_onlyTSS1kb_minus, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/gene_v41_human_onlyTSS1kb_minus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)

cat gene_v41_human_onlyTSS1kb_plus.txt gene_v41_human_onlyTSS1kb_minus.txt | sort -k1,1 -k2,2n > gene_v41_human_onlyTSS1kb_plus_minus.txt
grep "chrM" gene_v41_human_onlyTSS1kb_plus_minus.txt -v > gene_v41_human_onlyTSS1kb_plus_minus_freeM.txt


#2kb
gene_v41_human_onlyTSS2kb <- gene_v41_human_out
# extend TSS code
gene_v41_human_onlyTSS2kb_plus <- gene_v41_human_onlyTSS2kb[which(gene_v41_human_onlyTSS2kb$V4 == "+"),]
# replace column 3 by column 2
gene_v41_human_onlyTSS2kb_plus$V3 <- gene_v41_human_onlyTSS2kb_plus$V2
gene_v41_human_onlyTSS2kb_plus$V2 <- gene_v41_human_onlyTSS2kb_plus$V2 - 2000

gene_v41_human_onlyTSS2kb_minus <- gene_v41_human_onlyTSS2kb[which(gene_v41_human_onlyTSS2kb$V4 == "-"),]
# replace column 2 by column 3
gene_v41_human_onlyTSS2kb_minus$V2 <- gene_v41_human_onlyTSS2kb_minus$V3
gene_v41_human_onlyTSS2kb_minus$V3 <- gene_v41_human_onlyTSS2kb_minus$V3 + 2000

# once written should write again untill run the above extend TSS code again, becuase re-subtraction should be avoided 
write.table(gene_v41_human_onlyTSS2kb_plus, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/gene_v41_human_onlyTSS2kb_plus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)
write.table(gene_v41_human_onlyTSS2kb_minus, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/gene_v41_human_onlyTSS2kb_minus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)

cat gene_v41_human_onlyTSS2kb_plus.txt gene_v41_human_onlyTSS2kb_minus.txt | sort -k1,1 -k2,2n > gene_v41_human_onlyTSS2kb_plus_minus.txt
grep "chrM" gene_v41_human_onlyTSS2kb_plus_minus.txt -v > gene_v41_human_onlyTSS2kb_plus_minus_freeM.txt

#3kb

gene_v41_human_onlyTSS3kb <- gene_v41_human_out
# extend TSS code
gene_v41_human_onlyTSS3kb_plus <- gene_v41_human_onlyTSS3kb[which(gene_v41_human_onlyTSS3kb$V4 == "+"),]
# replace column 3 by column 2
gene_v41_human_onlyTSS3kb_plus$V3 <- gene_v41_human_onlyTSS3kb_plus$V2
gene_v41_human_onlyTSS3kb_plus$V2 <- gene_v41_human_onlyTSS3kb_plus$V2 - 3000

gene_v41_human_onlyTSS3kb_minus <- gene_v41_human_onlyTSS3kb[which(gene_v41_human_onlyTSS3kb$V4 == "-"),]
# replace column 2 by column 3
gene_v41_human_onlyTSS3kb_minus$V2 <- gene_v41_human_onlyTSS3kb_minus$V3
gene_v41_human_onlyTSS3kb_minus$V3 <- gene_v41_human_onlyTSS3kb_minus$V3 + 3000

# once written should write again untill run the above extend TSS code again, becuase re-subtraction should be avoided 
write.table(gene_v41_human_onlyTSS3kb_plus, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/gene_v41_human_onlyTSS3kb_plus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)
write.table(gene_v41_human_onlyTSS3kb_minus, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/gene_v41_human_onlyTSS3kb_minus.txt", sep="\t", append=F, quote=F, row.names = F, col.names = F)

cat gene_v41_human_onlyTSS3kb_plus.txt gene_v41_human_onlyTSS3kb_minus.txt | sort -k1,1 -k2,2n > gene_v41_human_onlyTSS3kb_plus_minus.txt
grep "chrM" gene_v41_human_onlyTSS3kb_plus_minus.txt -v > gene_v41_human_onlyTSS3kb_plus_minus_freeM.txt


# rename as the names are too big
mv gene_v41_human_onlyTSS3kb_plus_minus_freeM.txt geneTSS_3kb_freeM.txt
mv gene_v41_human_onlyTSS2kb_plus_minus_freeM.txt geneTSS_2kb_freeM.txt
mv gene_v41_human_onlyTSS1kb_plus_minus_freeM.txt geneTSS_1kb_freeM.txt


/mnt/home3/reid/av638/ENCODE/slurm_sub_re.py -j job_gu ./quantify_over_genes_upstream_regions.sh
# for vars in *.name.bam
# do
# echo ${vars}
# for genes_upstream in geneTSS*kb*freeM.txt
# do
# echo ${genes_upstream}
# bedtools coverage -a ${genes_upstream} -b ${vars}_chr.bedpe.bed | sort -k1,1 -k2,2n > ${vars}_${genes_upstream}_coverage_pe.bed
# done
# 
# done

atacseqkd_genesTSS_path <- "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/"
atacseqkd_genesTSS_cov_files <- list.files(atacseqkd_genesTSS_path, pattern = "_freeM.txt_coverage_pe.bed")
atacseqkd_genesTSS_coverage_list <- list()

for (gene_cov_file in atacseqkd_genesTSS_cov_files){
  print(gene_cov_file)
  gene_cov_file_tag <- gsub("_freeM.txt_coverage_pe.bed","",gsub(".mLb.clN.sorted.bam.name.bam_geneTSS","", gene_cov_file))
  print(gene_cov_file_tag)
  atacseqkd_genesTSS_coverage_list[[gene_cov_file_tag]] <- data.frame(fread(paste0(atacseqkd_genesTSS_path, gene_cov_file)))
}

atacseqkd_genesTSS_coverage_df <- do.call(cbind.data.frame, atacseqkd_genesTSS_coverage_list)

atacseqkd_genesTSS_precountdata <- atacseqkd_genesTSS_coverage_df[,c(6,7, seq(8,dim(atacseqkd_genesTSS_coverage_df)[2],11))]
colnames(atacseqkd_genesTSS_precountdata) <- c("ensID","Gene",colnames(atacseqkd_genesTSS_precountdata)[3:dim(atacseqkd_genesTSS_precountdata)[2]])
rownames(atacseqkd_genesTSS_precountdata) <- paste0(atacseqkd_genesTSS_precountdata$ensID, "%",atacseqkd_genesTSS_precountdata$Gene)
atacseqkd_genesTSS_precountdata <- atacseqkd_genesTSS_precountdata[,-c(1,2)]
colnames(atacseqkd_genesTSS_precountdata) <- gsub(".V8","",colnames(atacseqkd_genesTSS_precountdata))


# use the default cpm function becuase it is a gene-count matrix so for ~60k gene we cannot expect the bias will occur fro any given sample and we do not need give separate lib size.
atacseqkd_genesTSS_precountdata_myCPM <- data.frame(edgeR::cpm(atacseqkd_genesTSS_precountdata))
atacseqkd_genesTSS_precountdata_myCPM_log <- log(atacseqkd_genesTSS_precountdata_myCPM + 1, 2)

# subtract sample - IgG

# Assuming your data is stored in a data frame named 'df'

library(dplyr)

# Columns to subtract
atacseqkd_sample_columns <- colnames(atacseqkd_genesTSS_precountdata_myCPM_log)[7:18]

atacseqkd_control_columns <- rep(colnames(atacseqkd_genesTSS_precountdata_myCPM_log)[1:6], 2)


# Subtract H3K27me3 from IgG for each corresponding condition
for (i in seq_along(atacseqkd_sample_columns)) {
  print(paste0(atacseqkd_sample_columns[i], "  ",atacseqkd_control_columns[i]))
  atacseqkd_genesTSS_precountdata_myCPM_log[paste0("diff_",atacseqkd_sample_columns[i], "_",atacseqkd_control_columns[i])] <- atacseqkd_genesTSS_precountdata_myCPM_log[atacseqkd_sample_columns[i]] - atacseqkd_genesTSS_precountdata_myCPM_log[atacseqkd_control_columns[i]]
}

# select only diff columns 
atacseqkd_genesTSS_precountdata_myCPM_log_diff <- atacseqkd_genesTSS_precountdata_myCPM_log %>% dplyr::select(contains('diff'))
atacseqkd_genesTSS_precountdata_myCPM_log_diff_st <- data.frame(stack(as.matrix(atacseqkd_genesTSS_precountdata_myCPM_log_diff)))
atacseqkd_genesTSS_precountdata_myCPM_log_diff_st$col <- as.character(atacseqkd_genesTSS_precountdata_myCPM_log_diff_st$col)
atacseqkd_genesTSS_precountdata_myCPM_log_diff_set <- merge(atacseqkd_genesTSS_precountdata_myCPM_log_diff_st, set_rnaseqde_merged_shC1_shS_up, by.x="row", by.y="Gene", all.x=F, all.y=T)
atacseqkd_genesTSS_precountdata_myCPM_log_diff_set["id"] <- paste0(unlist(lapply(strsplit(atacseqkd_genesTSS_precountdata_myCPM_log_diff_set$col,"_"), function(x) paste0(x[2:4],collapse = "_"))))

atacseqkd_genesTSS_precountdata_myCPM_log_diff_set_1kb <- atacseqkd_genesTSS_precountdata_myCPM_log_diff_set %>% dplyr::filter(grepl('1kb', col))
atacseqkd_genesTSS_precountdata_myCPM_log_diff_set_2kb <- atacseqkd_genesTSS_precountdata_myCPM_log_diff_set %>% dplyr::filter(grepl('2kb', col))
atacseqkd_genesTSS_precountdata_myCPM_log_diff_set_3kb <- atacseqkd_genesTSS_precountdata_myCPM_log_diff_set %>% dplyr::filter(grepl('3kb', col))

ggplot(atacseqkd_genesTSS_precountdata_myCPM_log_diff_set_2kb, aes(x=value, color=label, linetype=id)) +
  geom_density()+
  scale_color_manual(values=c("#00D100", "black","orange","blue"))+
  scale_linetype_manual(values = c(rep("solid",1),rep("solid",1),rep("dotted",1),rep("dotted",1))) + theme_classic() + xlim(c(-2,3)) +
  geom_vline(aes(xintercept=0),color="grey", linetype="dashed", size=0.2)


################################################ END OF ANALYSIS ###########################################################

