setwd("/mnt/home3/reid/av638/methylseq/gurdonlab/outfolder/methylkit")
#Reading the methylation call files
library(methylKit)
# Import deduplicated, sorted, BAM files
# Store bam files in list
bam_files_list <- as.list(list.files(path = "/mnt/home3/reid/av638/methylseq/gurdonlab/outfolder/bismark_deduplicated/",
                                     pattern = "\\.deduplicated.sorted.bam$",
                                     full.names = TRUE))


# List of sample IDs
sample_ids_list <- lapply(str_split(bam_files_list,"/"), function(x) gsub("_R1_val_1_bismark_bt2_pe.deduplicated.sorted.bam","",x[11]))


# Set minimum CpG coverage desired
# Used in processBismarkAln function
min_coverage <- 5

# Set minimum quality of reads
# Used in processBismarkAln function
min_qual <- 20


# Set minimum methylation percentage difference between groups
# Used in getMethylDiff function; 25 is the default value.
dml_diffs <- 25


#Import from BAM
# Get methylation stats for CpGs with at least min_coverage coverage
meth_stats <- processBismarkAln(location = bam_files_list,
                                sample.id = sample_ids_list,
                                assembly = "GRCm39.fa ",
                                save.folder = NULL, 
                                save.context = c("CpG"), 
                                read.context = "CpG",
                                mincov = min_coverage,
                                minqual = min_qual, 
                                phred64 = FALSE, #Your data is phred33,
                                treatment = rep(0, length(bam_files_list)),
                                save.db = FALSE)
# File count
nFiles <- length(bam_files_list)

#Import from CpG report
#First convert .cov.gz file to .CpG.report.txt.gz using convert_cov_to_CpGreport.sh (custom)
# #test
# #<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
# 
# test_report_list <- as.list(list.files(path = "/mnt/home3/reid/av638/methylseq/gurdonlab/outfolder/bismark_methylation_calls/methylation_coverage",
#                                 pattern = "\\.cov.CpG_test_report.txt.gz",
#                                 full.names = TRUE))
# 
# 
# for (content1 in test_report_list){
#   #print(content1)
#   temp1 <- gsub("/mnt/home3/reid/av638/methylseq/gurdonlab/outfolder/bismark_methylation_calls/methylation_coverage/", "", content1)
#   temp1 <- gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.CpG_test_report.txt.gz", "", temp1)
#   print(temp1)
#   test_report <- fread(content1)
#   colnames(test_report) <- c("chr", "base", "strand", "methylated", "unmethylated", "Ccontext", "trinu_context")
#   test_report <- data.frame(test_report)
#   test_report["chrBase"] <- paste0(test_report$chr,
#                                    ".",
#                                    test_report$base)
# 
#   test_report["coverage"] <- test_report$methylated + test_report$unmethylated
#   test_report["freqC"] <- (test_report$methylated  *  100) / (test_report$methylated + test_report$unmethylated)
#   test_report["freqT"] <- (test_report$unmethylated  *  100) / (test_report$methylated + test_report$unmethylated)
#   test_report <- test_report[,c(8,1,2,3,9,10,11)]
#   write.table(test_report, paste0("/mnt/home3/reid/av638/methylseq/gurdonlab/outfolder/bismark_methylation_calls/methylation_coverage/",temp1,".CpG_test_report.txt"), sep="\t", quote = F, append=F, row.names = F, col.names = T)
# }
# #Create a list of myCpG_report.txt
# test_list <- lapply(test_report_list, function(x) gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov","",(gsub(".gz","",x))))
# test.ids <- lapply(test_list, function(x) gsub("/mnt/home3/reid/av638/methylseq/gurdonlab/outfolder/bismark_methylation_calls/methylation_coverage/","",(gsub(".CpG_test_report.txt","",x))))
# 
# # List of sample IDs
# sample_ids_reports_list <- lapply(str_split(CpG.reports.list,"/"), function(x) gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.CpG_report.txt.gz","",x[11]))
# 
# 
# # read the files to a methylRawList object: myobj
# myobj=methRead(test_list,
#                sample.id=test.ids,
#                assembly="GRCm39",
#                treatment=rep(0, length(test_list)),
#                context="CpG",
#                mincov = 5
# )

#---Now run the script on real data----#
#samples
#<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

CpG.reports.list <- as.list(list.files(path = "/mnt/home3/reid/av638/methylseq/gurdonlab/outfolder/bismark_methylation_calls/methylation_coverage",
                                       pattern = "\\.cov.CpG_report.txt.gz",
                                       full.names = TRUE))

# File count
nFiles <- length(CpG.reports.list)

for (report_content in CpG.reports.list){
  #print(report_content)
  temp_content <- gsub("/mnt/home3/reid/av638/methylseq/gurdonlab/outfolder/bismark_methylation_calls/methylation_coverage/", "", report_content)
  temp_content <- gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.CpG_report.txt.gz", "", temp_content)
  print(temp_content)
  report <- fread(report_content)
  colnames(report) <- c("chr", "base", "strand", "methylated", "unmethylated", "Ccontext", "trinu_context")
  report <- data.frame(report)
  report["chrBase"] <- paste0(report$chr,
                              ".",
                              report$base)
  
  report["coverage"] <- report$methylated + report$unmethylated
  report["freqC"] <- (report$methylated  *  100) / (report$methylated + report$unmethylated)
  report["freqT"] <- (report$unmethylated  *  100) / (report$methylated + report$unmethylated)
  report <- report[,c(8,1,2,3,9,10,11)]
  write.table(report, paste0("/mnt/home3/reid/av638/methylseq/gurdonlab/outfolder/bismark_methylation_calls/methylation_coverage/",temp_content,".CpG_report.txt"), sep="\t", quote = F, append=F, row.names = F, col.names = T)
  rm(report)
  rm(temp_content)
}
#Create a list of myCpG_report.txt
mysamples_list <- lapply(CpG.reports.list, function(x) gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov","",(gsub(".gz","",x))))
samples.id <- lapply(mysamples_list, function(x) gsub("/mnt/home3/reid/av638/methylseq/gurdonlab/outfolder/bismark_methylation_calls/methylation_coverage/","",(gsub(".CpG_report.txt","",x))))



# read the files to a methylRawList object: myobj
myobjDB=methRead(mysamples_list,
               sample.id=samples.id,
               assembly="GRCm39",
               treatment=rep(0, length(mysamples_list)),
               context="CpG",
               mincov = 5,
               dbtype = "tabix",
               dbdir = "methylDB")

print(myobjDB[[1]]@dbpath)

# Generate and save histograms showing Percent CpG Methylation
png("getMethylationStats.png", height = 1000, width = 1500)
par(mfrow=c(4,4))
for (i in 1:nFiles) {
  print(getSampleID(myobjDB[[i]]))
  getMethylationStats(myobjDB[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
}
dev.off()

# Generate and save histograms showing CpG Methylation Coverage
png("getCoverageStats.png", height = 1000, width = 1500)
par(mfrow=c(4,4))
for (i in 1:nFiles) {
  print(getSampleID(myobjDB[[i]]))
  getCoverageStats(myobjDB[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
}
dev.off()


#Filtering samples based on read coverage
filtered.myobjDB <- filterByCoverage(myobjDB, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

#Normalize coverage
norm.filt.objDB <- normalizeCoverage(filtered.myobjDB,method="median")


#Comparative analysis
#Merging samples
#Note: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting
merge.norm.filt.objDB <- unite(norm.filt.objDB, destrand=FALSE)
head(merge.norm.filt.objDB)

#By default, unite function produces bases/regions covered in all samples. 
#That requirement can be relaxed using “min.per.group” option in unite function.
#This creates a methylBase object, where only CpGs covered with at least 1 sample per group will be returned
# there were two groups defined by the treatment vector, 
# given during the creation of myobj: treatment=c(1,1,0,0)
#Not doing it meth.min=unite(myobj,min.per.group=1L)

#Sample Correlation
png("getCorrelation.png", height = 2000, width = 2000)
getCorrelation(merge.norm.filt.objDB,plot=TRUE)
dev.off()

# Clustering dendrogram
png("dendrogram_path.png", height = 600, width = 1000) 
clusterSamples(merge.norm.filt.objDB, dist="correlation", method="ward", plot=TRUE)
dev.off()

# Run a PCA analysis on percent methylation for all samples
png("pca_path.png", height = 1000, width = 1000) 
PCASamples(merge.norm.filt.objDB)
dev.off()

#Run the PCA analysis and plot variances against PC number in a screeplot
png("scree_path.png", height = 1000, width = 1000)
PCASamples(merge.norm.filt.objDB, screeplot = TRUE)
dev.off()


merge.norm.filt.objDB_data <- getData(merge.norm.filt.objDB)
colnames(merge.norm.filt.objDB_data) <- c(colnames(merge.norm.filt.objDB_data)[1:4], 
                                          paste0(unlist(lapply(getSampleID(merge.norm.filt.objDB), function(x) rep(x,3))), "_", 
                                                 colnames(getData(merge.norm.filt.objDB)[5:52])))

#For some situations, it might be desirable to summarize methylation information over tiling windows rather than doing base-pair resolution analysis.
#Tiling windows analysis
myobjDB.low <- methRead(mysamples_list,
                 sample.id=samples.id,
                 assembly="GRCm39",
                 treatment=rep(0, length(mysamples_list)),
                 context="CpG",
                 mincov = 2)

tiles = tileMethylCounts(myobjDB.low,win.size=1000,step.size=1000,cov.bases = 10)
head(tiles[[1]],3)


#Finding differentially methylated bases or regions
myDiff=calculateDiffMeth(merge.norm.filt.objDB)


# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
#
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
#
#
# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

#visualize the distribution of hypo/hyper-methylated bases/regions per chromosome
diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)



#Finding differentially methylated bases using multiple-cores
myDiff=calculateDiffMeth(meth,mc.cores=2)

#Annotating differentially methylated bases or regions
library(genomation)

# read the gene BED file
gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
                                            package = "methylKit"))
# annotate differentially methylated CpGs with 
# promoter/exon/intron using annotation data
#
annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)


# Read the CpG island annotation and annotate our differentially methylated bases/regions with them.
# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
                                     package = "methylKit"),
                         feature.flank.name=c("CpGi","shores"))
#
# convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

# read the CpG island annotation and annotate our differentially methylated bases/regions with them.
# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
                                     package = "methylKit"),
                         feature.flank.name=c("CpGi","shores"))
#
# convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

#Regional analysis
#summarize methylation information over a set of defined regions such as promoters or CpG islands.
promoters=regionCounts(myobj,gene.obj$promoters)

head(promoters[[1]])

#get the distance to TSS and nearest gene name
diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

# target.row is the row number in myDiff25p
head(getAssociationWithTSS(diffAnn))

#It is also desirable to get percentage/number of differentially methylated regions that overlap with intron/exon/promoters
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

#plot the percentage of differentially methylated bases overlapping with exon/intron/promoters
plotTargetAnnotation(diffAnn,precedence=TRUE,
                     main="differential methylation annotation")
#percentage of differentially methylated bases are on CpG islands, CpG island shores and other regions.
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
                     main="differential methylation annotation")
#percentage of intron/exon/promoters that overlap with differentially methylated bases.
getFeatsWithTargetsStats(diffAnn,percentage=TRUE)

#In case want to read back or share tabix file so created with methread, follow these steps:
#Loading methylDB objects from tabix files
#baseDB.obj <- makeMethylDB(methylBase.obj,"my/path")
#mydbpath <- getDBPath(baseDB.obj)
#rm(baseDB.obj)

#methylKit:::checkTabixHeader(mydbpath)
#readMethylDB(mydbpath)

#----------------------------  END OF ANALYSIS -------------------------#
#Full tested scripts are present in /mnt/home3/reid/av638/tutorial/methylkit/spare_scripts.R
