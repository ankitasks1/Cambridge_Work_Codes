################################################
###          Seqwalker               ###
################################################
# In a README file mention that this functions works best with the list to store the output

#' deseq2 load rdata
#' @param path Location of the file
#' @param rdata RData object
#'
#' @return dds object
#' @examples
#' dds <- load_rdata(rnaseqkd_path, "outfolder/star_salmon/deseq2_qc/deseq2.dds.RData")
#' @export load_rdata

load_rdata <- function(path, rdata) {
  load(paste0(path, rdata))
  return(dds)
}

#' import genome file
#' @param path Location of the file
#' @param file Gene file in a tab-separated format in order chr,start, end, strand, ensembl ID, gene symbol
#'
#' @return Gene file in a specified format
#' @examples gene <- read_genefile(rnaseqkd_path, "gene_gencode_human__gencode_out_chr.txt")
#' @export read_genefile
read_genefile <- function(path, file){
  out <- read.table(paste0(path, file), header = F, stringsAsFactors = F)
  colnames(out) <- c("chr","start", "end", "strand", "type", "ensID", "gene")
  out["ens_gene"] <- paste0(out$ensID, "%", out$gene)
  return(out)
}

#' deseq2 import raw counts
#' @param dds DESeq2 dds object
#'
#' @return raw counts from dds object
#' @examples counts <- dds_counts(rnaseqkd_list$dds)
#' @export dds_counts
dds_counts <- function(dds){
  counts <- DESeq2::counts(dds, normalized=FALSE)
  return(counts)
}

#' deseq2 filter samples
#' @param countmatrix Gene expression raw counts matrix
#' @param samplefilter Sample to remove if not to be apart of analysis
#' @param columnsremove Index of those samples on column of datagrame
#'
#' @return Filtered raw count matrix selected for desired samples)
#' @examples filtered_count_matrix <- sample_filter(rnaseqkd_list$counts, samplefilter = TRUE, c(7:9))
#' @export sample_filter
sample_filter <- function(countmatrix, samplefilter, columnsremove){
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

#' deseq2 reassign gene names
#' @param countmatrix count matrix
#' @param genedf dataframe of gene coordinates and names
#' @param columnskeep columns you want to keep from the count matrix
#'
#' @return filtered count matrix
#' @examples count_mat <- reassign_genenames(rnaseqkd_list$samples_selected, rnaseqkd_list$gene, c(2:10))
#' @export reassign_genenames
reassign_genenames <- function(countmatrix, genedf, columnskeep){
  genecounts <- merge(countmatrix, genedf ,by.x = "ensID", by.y ="ensID", all.y =T)
  rownames(genecounts) <- genecounts$ens_gene
  genecounts <- genecounts[,columnskeep]
  return(genecounts)
}

#' deseq2 coldata
#' @param dds DESEq2 dds object
#' @param samplefilter If any samples to filter otherwise skip
#' @param rowsremove Index of samples to filter
#'
#' @return coldata/ meta data df
#' @examples coldata <- make_coldata_dds(rnaseqkd_list$dds, samplefilter = TRUE, c(7:9))
#' @export make_coldata_dds
make_coldata_dds <- function(dds, samplefilter, rowsremove){
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

#' deseq2 prepare deseq2 matrix and perform differential analysis
#' @param counts Count matrix (gene are rows and samples are columns and values are unnormalized)
#' @param coldata Meta data
#' @param fdr p.adj-value threshold
#' @param fc fold change threshold
#'
#' @return list (ALL, DE, UP and DOWN) of differentially expressed genes based on a given threshold
#' @examples ddsde <- deseq2(rnaseqkd_list$processed_counts, rnaseqkd_list$coldata, 0.05, 2)
#' @export deseq2_de
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

#' deseq2 plot venn diagram
#' @param list list of samples (n < = 4) with vector of genes symbols
#' @param color vector of colors you want to give eg. c("black", "green", "red") for 3 samples
#'
#' @return venn diagram
#' @examples venn <- plot_venn(list_rnaseqkd_de_shC1_c1509_shS_c1509, c("#0073C2FF", "#EFC000FF","#0073C2FF", "#EFC000FF"))
#' @export plot_venn
plot_venn <- function(list, color){
  ggvenn_list_venn <- ggvenn(list, fill_color = color, stroke_size = 0.5, set_name_size = 4)
  return(ggvenn_list_venn)
}

#' edger create coldata/design, normalise, fit, mds plot
#' @param counts Count matrix (gene are rows and samples are columns and values are unnormalized)
#' @param coldata Meta data
#'
#' @return prepare list of object required by edgeR
#' @examples dge <- edger_create(atacseqkd_edger_list$processed_counts, atacseqkd_edger_list$coldata)
#' @export
edger_create <- function(data, coldata){
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

#' edger differential analysis
#' @param lmfit linear models fitted object prepared while performing edger_create()
#' @param design design matrices for specified variables
#' @param fdr p.adj-value threshold
#' @param fc fold change threshold
#'
#' @return list (ALL, DE, UP and DOWN) of differentially expressed genes based on a given threshold
#' @examples dar_analysis <- fn_edger_de(atacseqkd_edger_list$dge_obj$lmfit, atacseqkd_edger_list$dge_obj$design, 0.05, 2)
#' @export edger_de
edger_de <- function(lmfit, design, fdr, fc){
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


#' diffbind prepare object
#' @param samplesheet Samplesheet formatted as per Diffbind (see Diffbind on bioconductor for detail)
#'
#' @return samplesheet
#' @examples ss <- diffbind_prep(diffbind_atacseqkd_samplesheet)
#' @export diffbind_prep
diffbind_prep <- function(samplesheet){
  diffblist <- list()
  # Creating diffbind object
  obj <- dba(sampleSheet=samplesheet,  scoreCol=5, minOverlap=1)
  diffblist[["dba_obj"]] <- obj
  # Plot of number of peaks that appear in atleast one, two , three so on ... until total samples together.
  ovp_rate <- dba.overlap(obj,mode=DBA_OLAP_RATE)
  diffblist[["ovp_rate"]] <- ovp_rate
  return(diffblist)
}

#' diffbind get counts
#' @param obj diffbind object
#' @param summit expected summit for peaks
#'
#' @return count matrix derived from diffbind
#' @examples counts <- diffbind_count(atacseqkd_diffbind_list$prep$dba_obj, 100)
#' @export diffbind_count
diffbind_count <- function(obj, summit){
  diffblist <- list()
  # Calculate a binding matrix with scores based on read counts,  summits=100 for ATAC-Seq # https://www.biostars.org/p/9493721/
  diffblist[["dba_obj"]] <- dba.count(obj, minOverlap=1, summits=100, score = DBA_SCORE_READS)
  diffblist[["info"]]  <- dba.show(diffblist$dba_obj)
  return(diffblist)
}

#' diffbind normalize counts and perform DE analysis
#' @param obj diffbind object
#'
#' @return Normalize and prepare contrast
#' @examples norm <- diffbind_norm(atacseqkd_diffbind_list$counts$dba_obj)
#' @export diffbind_norm
diffbind_norm <- function(obj){
  diffblist <- list()
  message("Performing normalisation...")
  diffblist[["norm"]] <- dba.normalize(obj, method=DBA_ALL_METHODS)
  message("Performing Diffbind analysis...")
  message("Performing Contrasts...")
  diffblist$norm <- dba.contrast(diffblist$norm, minMembers = 2,categories=DBA_CONDITION)
  diffblist$norm <- dba.analyze(diffblist$norm, method=DBA_ALL_METHODS)
  return(diffblist)
}


#' diffbind extract differential analysis results
#' @param normobj diffbind normalized object
#' @param fdr p.adj-value threshold
#' @param fc fold change threshold
#' @param method edgeR, deseq2
#'
#' @return differentially analysis outputs by diffbind
#' @examples 
#' db_edgeR <- diffbind_de(atacseqkd_diffbind_list$norm$norm, 0.05, 2, "edgeR")
#' db_deseq2 <- diffbind_de(atacseqkd_diffbind_list$norm$norm, 0.05, 2, "deseq2")
#'
#' @export  diffbind_de
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


#' limma create object
#' @param data Count matrix
#' @param coldata Meta data
#'
#' @return Output a list of objects required for limma
#' @examples dge <- limma_create(atacseqkd_limma_list$processed_counts, atacseqkd_limma_list$coldata)
#' @export limma_create
limma_create <- function(data, coldata){
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

#' limma differential analysis
#' @param lmfit linear fitted model object
#' @param contrasts contarst for samples (use makeContrasts before running this function)
#' @param fdr p.adj-value threshold
#' @param fc fold change threshold
#'
#' @return list containing differential analysis results
#' @examples delist <- limma_de(atacseqkd_limma_list$dge_obj$fit, atacseqkd_limma_list$contrasts, 0.05, 2)
#' @export limma_de
limma_de <- function(lmfit, contrasts, fdr, fc){
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

#' diffbind annotation to genes and nearest distance filter for all packages
#' @param assaytype eg. atacseq, rnaseq
#' @param software internal software to be used by Diffbind eg. edger or deseq2
#' @param contrasts contrasts list
#' @param genefile gene file name
#' @param path location of files
#' @param outformat eg. ".txt"
#' @param fdr p.adj-value threshold
#' @param fc fold change threshold
#'
#' @return list of annotated differential regions
#' @examples annotation <- diffbind_gene_anno("atacseqkd", "diffbind",atacseqkd_diffbind_list$dar_analysis$deseq2$contrasts, atacseqkd_diffbind_list$gene, "/mnt/home3/atacseq/outfolder", "txt", 0.05, 2)
#' @export diffbind_gene_anno
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

#' diffbind annotation to genes and nearest distance filter for all packages on diffbind derived counts
#' @param assaytype eg. atacseq, rnaseq
#' @param software internal software to be used by Diffbind eg. edger or deseq2
#' @param contrasts contrasts list
#' @param genefile gene file name
#' @param path location of files
#' @param outformat eg. ".txt"
#' @param fdr p.adj-value threshold
#' @param fc fold change threshold
#'
#' @return list of annotated differential regions
#' @examples annotation <- diffbind_and_packages_gene_anno("atacseqkd", "diffbind_deseq2", atacseqkd_diffbind_deseq2_list$dar_analysis$de_analysis, atacseqkd_diffbind_deseq2_list$gene, "/mnt/home3/outfolder", "txt", 0.05, 2)
#' @export diffbind_and_packages_gene_anno
diffbind_and_packages_gene_anno <- function(assaytype, software, contrasts, genefile, path, outformat, fdr, fc){
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

#' diffbind consensus gene annotation
#' @description consensus peaks from nextflow output but can be applied to other sets
#' @param assaytype eg. atacseq, rnaseq
#' @param software internal software to be used by Diffbind eg. edger or deseq2
#' @param contrasts contrasts list
#' @param genefile gene file name
#' @param path location of files
#' @param outformat eg. ".txt"
#' @param fdr p.adj-value threshold
#' @param fc fold change threshold
#'
#' @return list of annotated differential regions consensus
#' @examples annotation <- consensus_gene_anno("atacseqkd", "deseq2", atacseqkd_deseq2_list$dar_analysis$de_analysis, atacseqkd_deseq2_list$consensus, atacseqkd_deseq2_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder", "txt", 0.05, 2)
#' @export consensus_gene_anno
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

#' diffbind consensus gene annotation for normalized count matrix (consensus peaks from nextflow output but can be applied to other sets)
#' @description consensus peaks from nextflow output but can be applied to other sets
#' @param assaytype eg. atacseq, rnaseq
#' @param software internal software to be used by Diffbind eg. edger or deseq2
#' @param contrasts contrasts list
#' @param genefile gene file name
#' @param path location of files
#' @param outformat eg. ".txt"
#'
#' @return list of annotated differential regions consensus
#' @examples annotation <- consensus_gene_count_anno("cutnrunko", "deseq2", cutnrunkd_deseq2_list$normalized_counts, cutnrunkd_deseq2_list$consensus, cutnrunkd_deseq2_list$gene, "/mnt/home3/reid/av638/cutnrun/iva_lab_oct23/cutnrun_k27_kd/outfolder/", "txt")
#' @export consensus_gene_count_anno
consensus_gene_count_anno <- function(assaytype, software, normcountmatrix, consensus_bed, genefile, path, outformat){
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

#' chipseeker run
#' @param peaks_path location of peaks files (eg. .bed , .narrowPeak, .broadPeak)
#' @param assaytype eg. atacseq
#' @param txdb eg. TxDb.Hsapiens.UCSC.hg38.knownGene (not quoted)
#' @param org_db eg. org.Hs.eg.db
#' @param peakformat eg. ".bed", all peak files should have same format
#'
#' @return list containing objects, granges, plots from chipseeker
#' @examples list <- run_chipseeker("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams/endoH1peaks/", "cutntag", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "*.broadPeak")
#' @export run_chipseeker
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


#' chipseeker run and annotation
#' @param path location of peaks files (eg. .bed , .narrowPeak, .broadPeak)
#' @param assaytype eg. atacseq
#' @param software eg. deseq2
#' @param contrasts contrasts list
#' @param txdb eg. TxDb.Hsapiens.UCSC.hg38.knownGene (not quoted)
#' @param org_db eg. org.Hs.eg.db
#' @param regionformat eg. ".bed", all peak files should have same format
#' @param fdr p.adj-value threshold
#' @param fc fold change threshold
#'
#' @return list containing objects, granges, plots from chipseeker
#' @examples annotation <- chipseeker_region_anno("atacseqkd", "deseq2", atacseqkd_deseq2_list$dar_analysis$de_analysis, atacseqkd_deseq2_list$consensus[,c(1:4)], "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/", TxDb.Hsapiens.UCSC.hg38.knownGene, "org.Hs.eg.db", "region_bed.txt", 0.05, 2)
#' @export chipseeker_region_anno
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


#' chipseeker pie plot
#' @param gr_anno granges annotation
#'
#' @return values for customising pie chart
#' @examples pie_val <- replot_pie_chipseeker(cutntagrmwt_quantify_list$chipseeker$plots$gr_anno)
#' @export replot_pie_chipseeker
replot_pie_chipseeker <- function(gr_anno){
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
# names=c("chr", "start", "end", "feature_id")
# quantify per feature in a given interval using bedtools coverage
# peaks_path =
# bamfile = "\\.name\\.bam$"
# bedfile = "_chr.bedpe"
# sites_file = "atac_shControl_merged.txt"

#' featurecounts quantify per feature in a given interval
#' @param peaks_path location of peaks files
#' @param assaytype eg. "atacseq"
#' @param bamfiles eg. "\\.bam$"
#' @param sites_files eg. "\\_peaks_id.bed$"
#' @param sites_type eg. "histone_marks"
#' @param sites_column_to_rearrange This is suitable for 12-bed format as suggested in Footnotes
#' eg. c(12,4,1:3,7). If you have specific format of peaks files please rearrange appropriately
#' but makes sure you give correct index
#' @param pairedend boolean {TRUE, FALSE}
#' @param refgenome eg. "hg38
#' @param delim eg. "_", the separator you want to use for assigning new id
#' @param merge_sites_files boolean {TRUE, FALSE} keep it FALSE
#'
#' @return It will out put the quantification of bam files withing a given peak(s)/region(s) of interest
#' @examples featurecounts_out <- quantify_featurecounts("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/", "atacseqkd","\\.bam$", "\\_peaks_id.bed$","histone_marks", c(12,4,1:3,7), pairedend=TRUE, "hg38", "_",merge_sites_files=FALSE)
#' @export quantify_featurecounts
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

#' featurecounts quantify for multiple features  subtract
#' @param assay_bam_path location of bam files
#' @param feature_path location of features files (eg. hg38_5kb.txt, gene_promoters.txt
#' @param assay eg. cutnandrun
#' @param bam_extension eg .bam
#' @param feature_list list of features to pick from the feature_path, all features shoudl have a specified arrangement of columns
#' @param columnstorearrange  column rearrangement required
#' @param ref eg. hg38
#' @param control eg. IgG/Control
#' @param pe boolean {TRUE, FALSE}
#' @param diff eg. "single"
#' @param pairs list of sampleID, sample and respective control
#'
#' @return quantification in dataframe format for multiple sampples for multiple feature of interest
#' @examples featurecounts_out <- quantify_featurecounts_multifeature("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "atacseqkd", ".mLb.clN.sorted.bam", atacseqkd_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "shControl",pe=TRUE, diff="multi", atacseqkd_quantify_list$featurecounts$pairs)
#' @export quantify_featurecounts_multifeature
quantify_featurecounts_multifeature <- function(assay_bam_path, feature_path, assay, bam_extension, feature_list, columnstorearrange, ref, control,pe=FALSE, diff="single", pairs){
  storelist <- list()
  for (i in names(feature_list)){
    message(paste0("copying ",feature_list[[i]][2], " to ", assay_bam_path, "..."))
    system(paste0("cp ", feature_path, feature_list[[i]][2], " ", assay_bam_path))
    message("performing featurecounts...")
    storelist[[i]] <- quantify_featurecounts(assay_bam_path, assay, paste0("\\",bam_extension,"$"), paste0("\\",feature_list[[i]][1],"$"), i, columnstorearrange, pairedend=pe, ref, "%",merge_sites_files=FALSE)
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

#' featurecounts quantify for multiple feature using featurecounts_div
#' @param assay_bam_path location of bam files
#' @param feature_path location of features files (eg. hg38_5kb.txt, gene_promoters.txt
#' @param assay eg. cutnandrun
#' @param bam_extension eg .bam
#' @param feature_list list of features to pick from the feature_path, all features shoudl have a specified arrangement of columns
#' @param columnstorearrange  column rearrangement required
#' @param ref eg. hg38
#' @param control eg. IgG/Control
#' @param pe boolean {TRUE, FALSE}
#' @param diff eg. "single"
#' @param pairs list of sampleID, sample and respective control
#'
#' @return quantification in dataframe format for multiple sampples for multiple feature of interest
#' @examples
#' featuresmatrixdiv <- quantify_featurecounts_multifeature_div("/mnt/home3/reid/av638/cutntag/iva_lab_feb2024/bams/", "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration/", "cutntagrmwt", "_sorted.bam",
#' cutntagrmwt_quantify_list$featurecounts$featureslist, c(4,5,1:3,6), "hg38", "IgG",pe=TRUE, diff="multi",
#' cutntagrmwt_quantify_list$featurecounts$pairs)
#' @export quantify_featurecounts_multifeature_div
quantify_featurecounts_multifeature_div <- function(assay_bam_path, feature_path, assay, bam_extension, feature_list, columnstorearrange, ref, control,pe=FALSE, diff="single", pairs){
  storelist <- list()
  for (i in names(feature_list)){
    message(paste0("copying ",feature_list[[i]][2], " to ", assay_bam_path, "..."))
    system(paste0("cp ", feature_path, feature_list[[i]][2], " ", assay_bam_path))
    message("performing featurecounts...")
    storelist[[i]] <- quantify_featurecounts(assay_bam_path, assay, paste0("\\",bam_extension,"$"), paste0("\\",feature_list[[i]][1],"$"), i, columnstorearrange, pairedend=pe, ref, "%",merge_sites_files=FALSE)
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

#' Rsamtools quantify reads in a given bam files
#' @param peaks_path location of files
#' @param assaytype eg. cutandrun
#' @param bamfiles extension , eg. "\\.bam$"
#'
#' @return list of files with reads quantified
#' @examples bams <- quantify_bams("/mnt/home3/reid/av638/atacseq/iva_lab_gencode/boutfolder/bowtie2/merged_library/", "atacseqkd","\\.bam$")
#' @export quantify_bams
quantify_bams <- function(peaks_path, assaytype, bamfiles){
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

#' Further process the featurecounts quantified matrix
#' @param fearturecountmatrix Count matrix
#' @param bamcountsreads Dataframe of total reads in bam files
#' @param labels Labels of features short name eg. H3K4me3_1 K4me3
#'
#' @return Summed counts of reads in a given features
#' @examples summed <- summedreads_per_feature(atacseqkd_quantify_list$featurecounts$histonemarks$histone_marks_countmatrix, atacseqkd_quantify_list$bams$total_reads, atacseqkd_labelpath)
#' @export summedreads_per_feature
summedreads_per_feature <- function(fearturecountmatrix, bamcountsreads, labels){
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


#' aggregate value per feature
#' @param annomatrix dataframe
#' @param columnstoaggregate column to use for aggreagtion, eg. values
#' @param aggregatebycolumn column will be by use for group by process eg. gene
#' @param operation eg. mean (not under the quotes)
#'
#' @return Dataframe of aggregated values by operation of interest
#' @examples agg <- aggregate_feature(atacseqkd_quantify_list$featurecounts$histonemarks$summed$merged_diff_table_afr_st, c(3), "id", mean)
#' @export aggregate_feature
aggregate_feature <- function(annomatrix, columnstoaggregate, aggregatebycolumn, operation){
  storelist <- list()
  message("aggregating...")
  storelist[["aggregated"]] <- stats::aggregate(annomatrix[,columnstoaggregate], by=list(annomatrix[[aggregatebycolumn]]), operation)
  return(storelist)
}

#' pca prcomp using counts
#' @param counts count matrix
#' @param file_extension bam file extension eg. ".bam"
#' @param attribute pca
#'
#' @return list of pca related objects and pca plot
#' @examples pca <- make_pca(cutntagwkd_quantify_list$featurecounts$featuresmatrix$bins5kb$bins5kb_countmatrix$hg38_5kb.txt$cpm, "_R1.target.markdup.sorted.bam", "pca")
#' @export make_pca
make_pca <- function(counts, file_extension, attribute){
  storelist <- list()
  message("filtering the counts for zero rowsums...")
  storelist[["norm_counts_filtered"]] <- counts[rowSums(counts) > 0,]
  colnames(storelist$norm_counts_filtered) <- gsub(file_extension,"", colnames(storelist$norm_counts_filtered))
  message("generating pca")
  storelist[["prcomp"]] <- prcomp(t(storelist$norm_counts_filtered), center = TRUE, scale. = TRUE)
  storelist[[attribute]] <- factoextra::fviz_pca_ind(storelist$prcomp, geom = c("point", "text"), repel = TRUE)
  return(storelist)
}

#' integrate multiple omics data on featurecounts
#' @param inputlist list of dataframes from different samples /results
#' @param sharedcolumn eg. ensID_gene
#'
#' @return Make a consesnus based on given column presence in all samples eg. ensID_gene
#' @examples merge <- integrate_featurecounts(data_integration_list$ensID_gene$list, "ensID_Gene")
#' @export integrate_featurecounts
integrate_featurecounts <- function(inputlist, sharedcolumn){
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

#' integrate multiple omics data on bedtools coverage output
#' @param omics_data_path_lists list of multiple outputs of bedtools coverage
#' @param binsize eg. "5kb
#' @param columns_pos eg. c(1:3)
#' @param columns_samples eg. vector of position of coverage column across all samples, can be extracted by seq(2,30,5) or manually c(2,7,12,17,22,27)
#'
#' @return dataframe for bedtools coverage performed on multiple samples
#' @examples bin5kb <- integrate_bedtools(bin_path_combined, "5kb", c(1:3), seq(4,175,7))
#' @export integrate_bedtools
integrate_bedtools <- function(omics_data_path_lists, binsize, columns_pos, columns_samples){
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


#' Upset plot
#' @param input_list list of samples with vector of features eg. genes
#'
#' @return list of all possible intersections
#' @examples upset <- Upsetout(atacseqkd_venn_packages_list_CRAMP1)
#' @export Upsetout
Upsetout <- function(input_list){
  if(missing(input_list)){
    stop("No input list is provided. Please provide an appropriate input list")
  }
  if(is.null(input_list)){
    stop("Input list is empty. Please add data your list")
  }else{
    message("Input list is correctly provided. performing interactions")
    # Combine all elements in all the the vectors
    all_elements <- unique(as.character(unlist(input_list)))
    # Create the binary list
    binary_list <- lapply(input_list, function(x) as.numeric(all_elements%in%x))
    #Bind the binary matrix
    binary_df <- data.frame(do.call(cbind.data.frame, binary_list))
    rownames(binary_df) <- all_elements
    # make a dataframe which adress intersection of all elemnts and respective samples
    for (names_i in colnames(binary_df)){
      binary_df_names_i <- binary_df[,names_i]
      #Insert name of column into the table elements
      binary_df[names_i] <- ifelse(binary_df_names_i >= 1, names_i, 0)

    }
    binary_df$combos <- apply(binary_df[ , colnames(binary_df)] , 1 , paste , collapse = "_" )

    # Prepare the max combination required
    combos_for_venn <- expand.grid(rep(list(0:1),length(names(input_list))))
    colnames(combos_for_venn) <- names(input_list)
    for (names_i in names(combos_for_venn)){
      combos_for_venn_names_i <- combos_for_venn[,names_i]
      #print(head(temp50_names_i))
      combos_for_venn[names_i] <- ifelse(combos_for_venn_names_i >= 1, names_i, 0)
    }
    combos_for_venn$combos <- apply(combos_for_venn[ , colnames(combos_for_venn)] , 1 , paste , collapse = "_" )
    # Extract combinations as output list
    finalout_upsetlist <- list()
    for (combos in combos_for_venn$combos){
      finalout_upsetlist[[combos]] <- binary_df[which(binary_df$combos == combos),]
    }
    # print(finalout_upsetlist)
    message("Combinations are ready. Explore the output list assigned by you")
    return(finalout_upsetlist)
  }
}

#' ggmaplot maplot
#' @param df dataframe of output from deseq2 or other ways
#' @param columnbyorder since columns we renamed again the exact order need to be specified
#' @param fdr p.adj-value threshold
#' @param fc fold change threshold
#'
#' @return MA plot
#' @examples ma <- maplot(atacseqkd_diffbind_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$all, c(2,1,5,7), 0.05, 2)
#' @export maplot
maplot <- function(df, columnbyorder, fdr, fc){
  plotlist <- list()
  message("rearranging and plotting ma...")
  df <- df[,columnbyorder]
  colnames(df) <- c("baseMean","log2FoldChange","padj","feature")
  plotlist[["maplot"]] <- ggmaplot(df, fdr = fdr, fc = fc, size = 0.3, palette = c("#D22B2B", "#1465AC", "darkgray"),
           genenames = df$feature, legend = "top", top = 20,
           font.label = c("bold", 5), font.legend = "bold", font.main = "bold", ggtheme = ggplot2::theme_minimal())
}

#' ggplot maplot
#' @param df dataframe of output from deseq2 or other ways
#' @param columnbyorder since columns we renamed again the exact order need to be specified
#' @param fdr p.adj-value threshold
#' @param fc fold change threshold
#'
#' @return MA plot
#' @examples ma <- maplot_general(atacseqkd_limma_list$dar_analysis$dars$coef_shSUZ12_shControl$all, c(2,1,5,7), 0.05, 1)
#' @export maplot_general
maplot_general <- function(df, columnbyorder, fdr, logfc){
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

#' deseq2 maplot using ggplot2 using shirnkage
#' @param df dataframe of output from deseq2 or other ways
#' @param columnbyorder since columns we renamed again the exact order need to be specified
#' @param fdr p.adj-value threshold
#' @param fc fold change threshold
#'
#' @return MA plot
#' @examples ma <- maplot_general_shrinkage(rnaseqkd_list$deseq2$de_analysis$shC1_c1509$all, c(1,3,6,7), 0.05, 1)
#' @export maplot_general_shrinkage
maplot_general_shrinkage <- function(df, columnbyorder, fdr, logfc){
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

#' base rearrange groups
#' @param grouplist list of path to specific genes(df) eg. cluster1.txt , cluster2.txt (prepared by list.files(integration_path, pattern="cluster*", full.names = T))
#'
#' @return filtered for random genes which you dont want in your data
#' @examples groups <- rearrange_groups(selected_rkd_gene_groups)
#' @export rearrange_groups
rearrange_groups <- function(grouplist){
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

#' ggplot meta plots
#' @param featurematrix count matrix
#' @param columnstoplot column with value to pick for plotting
#' @param multi FALSE
#' @param value column name with values to pick eg. "value",
#' @param color column name with col to pick  eg. "col"
#' @param linetype column name with col to pick  eg. "col"
#' @param rep eg. 1
#' @param splitside delimiter separator side 1,2,3,4,5,..... eg. 3
#' @param colorbyid FALSE
#'
#' @return Meta plots / Average plots / Density plot
#' @examples meta <- meta_plot(atacseqkd_quantify_list$featurecounts$featuresmatrix, c(7:10), rep=1, splitside=3)
#' @export meta_plot
meta_plot <- function(featurematrix, columnstoplot, multi=FALSE, value="value",color="col", linetype="col", rep, splitside=3, colorbyid=FALSE){
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

#' bedtools annotation
#' @param peakfileslist list of peak file
#' @param genefile datframe of gene file
#' @param path location of peaks files
#' @param assaytype eg. "atacseq"
#' @param software eg. "deseq2"
#' @param outformat eg. ".txt"
#'
#' @return distance calculated by peakfile and gene using bedtools closest  with -d option
#' @examples anno <- annotate_peaks_to_genes(public_list$peakfilelist, public_list$gene, "/mnt/home3/reid/av638/atacseq/iva_lab_gencode/integration","chipseq", "encode", "txt")
#' @export annotate_peaks_to_genes
annotate_peaks_to_genes <- function(peakfileslist, genefile, path, assaytype, software, outformat){
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

#' gprofiler go term analysis
#' @param gene_list list of mutiple samples with genes of vectors
#' @param organism eg. "hg38"
#'
#' @return list of outputs related to gprofiler and gost
#' @examples go <- go_term_analysis(list_rnaseqkd_de_shC1_c1509_shS_c1509, "hsapiens")
#' @export go_term_analysis
go_term_analysis <- function(gene_list, organism){
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

#' ggplot scatter
#' @param df dataframe
#'
#' @return scatterplot
#' @examples sct <- scatterplot(cutntagrmwt_correlation_scatter_divlog2)
#' @export scatterplot
scatterplot <- function(df){
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

#' enrichment_counter
#' @param list1 peaklist1
#' @param list2 peaklist2
#'
#' @return dataframe of overlap between two peak files
#' @examples df <- enrichment_counter(histonemarks_list1, atacseqkd_dars_list1)
#' @export enrichment_counter
enrichment_counter <- function(list1, list2){
  # Create an empty dataframe to store the results
  result_df <- data.frame(matrix(ncol = 7, nrow = length(list1) * length(list2)))  # Assuming you want to store 7 values for each combination
  colnames(result_df) <- c("File1", "File2", "Coverage_File2", "Coverage_File1", "Intersect_Count", "Intersect_Ratio_File2", "Intersect_Ratio_File1")
  index <- 1
  for (i in list1){
    for (j in list2){
      # print(i)
      in_coverage <- as.numeric(fread(cmd=paste0("bedtools intersect -a ", paste0(i), " -b ", paste0(j, " | wc -l"))))
      cov2 = as.numeric(length(readLines(file(j))))
      cov1 = as.numeric(length(readLines(file(i))))
      print(c(basename(j), basename(i), cov2, cov1, in_coverage, in_coverage/cov2, in_coverage/cov1))
      result_df[index, 1] <- basename(j)
      result_df[index, 2] <- basename(i)
      result_df[index, 3] <- cov2
      result_df[index, 4] <- cov1
      result_df[index, 5] <- in_coverage
      result_df[index, 6] <- in_coverage/cov2
      result_df[index, 7] <- in_coverage/cov1
      
      # Increment the index
      index <- index + 1
    }
  }
  return(result_df)
}
