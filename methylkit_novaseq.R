setwd("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit")
#Reading the methylation call files
library(methylKit)


cov_folder <- "/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/bismark_methylation_calls/methylation_coverage/"

#Import from CpG report
#First convert .cov.gz file to .CpG.report.txt.gz using convert_cov_to_CpGreport.sh (custom)
# #<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
# tested as above
#---Now run the script on real data----#
#samples
#<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>

CpG.reports.list <- as.list(list.files(path = cov_folder,
                                       pattern = "\\.cov.CpG_report.txt.gz",
                                       full.names = TRUE))

# File count
nFiles <- length(CpG.reports.list)

# for (report_content in CpG.reports.list){
#   #print(report_content)
#   temp_content <- gsub(paste0(cov_folder,"/"), "", report_content)
#   temp_content <- gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.CpG_report.txt.gz", "", temp_content)
#   print(temp_content)
#   report <- fread(report_content)
#   colnames(report) <- c("chr", "base", "strand", "methylated", "unmethylated", "Ccontext", "trinu_context")
#   report <- data.frame(report)
#   report["chrBase"] <- paste0(report$chr,
#                               ".",
#                               report$base)
#   
#   report["coverage"] <- report$methylated + report$unmethylated
#   report["freqC"] <- (report$methylated  *  100) / (report$methylated + report$unmethylated)
#   report["freqT"] <- (report$unmethylated  *  100) / (report$methylated + report$unmethylated)
#   report <- report[,c(8,1,2,3,9,10,11)]
#   write.table(report, paste0(cov_folder,temp_content,".CpG_report.txt"), sep="\t", quote = F, append=F, row.names = F, col.names = T)
#   rm(report)
#   rm(temp_content)
# }
# #Create a list of myCpG_report.txt
mysamples_list <- lapply(CpG.reports.list, function(x) gsub("_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov","",(gsub(".gz","",x))))
samples.id <- lapply(mysamples_list, function(x) gsub(paste0(cov_folder,"/"),"",(gsub(".CpG_report.txt","",x))))
# 
# 
# # 
# # # read the files to a methylRawList object: myobj
# myobjDB=methRead(mysamples_list,
#                sample.id=samples.id,
#                assembly="GRCm39",
#                treatment=rep(0, length(mysamples_list)),
#                context="CpG",
#                mincov = 3)
# # 
# 
# # 
# # Generate and save histograms showing Percent CpG Methylation
# png("getMethylationStats.png", height = 1000, width = 1500)
# par(mfrow=c(4,4))
# for (i in 1:nFiles) {
#   print(getSampleID(myobjDB[[i]]))
#   getMethylationStats(myobjDB[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
# }
# dev.off()
# # 
# # # Generate and save histograms showing CpG Methylation Coverage
# # png("getCoverageStats.png", height = 1000, width = 1500)
# # par(mfrow=c(4,4))
# # for (i in 1:nFiles) {
# #   print(getSampleID(myobjDB[[i]]))
# #   getCoverageStats(myobjDB[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
# # }
# # dev.off()
# # 
# # 
# # #Filtering samples based on read coverage
# filtered.myobjDB <- filterByCoverage(myobjDB, lo.count=0, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
# # #Normalize coverage
# # norm.filt.objDB <- normalizeCoverage(filtered.myobjDB,method="median")
# # 
# # 
# # #Comparative analysis
# # #Merging samples
# # #Note: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting
# merge.filt.objDB <- unite(filtered.myobjDB, destrand=FALSE)
# head(merge.filt.objDB)
# merge.filt.objDB.df <- getData(merge.filt.objDB)
# merge.filt.objDB.df <- merge.filt.objDB.df[order(merge.filt.objDB.df$chr, merge.filt.objDB.df$start),]
# write.table(merge.filt.objDB.df, 'merge.filt.objDB.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
# grep chr merge.filt.objDB.df.txt > merge.filt.objDB.df_chr.txt
# # 
# #By default, unite function produces bases/regions covered in all samples. 
# #That requirement can be relaxed using “min.per.group” option in unite function.
# #This creates a methylBase object, where only CpGs covered with at least 1 sample per group will be returned
# # there were two groups defined by the treatment vector, 
# # given during the creation of myobj: treatment=c(1,1,0,0)
# #Not doing it meth.min=unite(myobj,min.per.group=1L)
# 
# #Sample Correlation
# png("getCorrelation.png", height = 2000, width = 2000)
# getCorrelation(merge.norm.filt.objDB,plot=TRUE)
# dev.off()
# 
# # Clustering dendrogram
# png("dendrogram_path.png", height = 600, width = 1000) 
# clusterSamples(merge.norm.filt.objDB, dist="correlation", method="ward", plot=TRUE)
# dev.off()
# 
# # Run a PCA analysis on percent methylation for all samples
# png("pca_path.png", height = 1000, width = 1000) 
# PCASamples(merge.norm.filt.objDB)
# dev.off()
# 
# #Run the PCA analysis and plot variances against PC number in a screeplot
# png("scree_path.png", height = 1000, width = 1000)
# PCASamples(merge.norm.filt.objDB, screeplot = TRUE)
# dev.off()
# 
# 
# merge.norm.filt.objDB_data <- getData(merge.norm.filt.objDB)
# colnames(merge.norm.filt.objDB_data) <- c(colnames(merge.norm.filt.objDB_data)[1:4], 
#                                           paste0(unlist(lapply(getSampleID(merge.norm.filt.objDB), function(x) rep(x,3))), "_", 
#                                                  colnames(getData(merge.norm.filt.objDB)[5:52])))
# 
# head(merge.norm.filt.objDB_data)
# 

# #DMRs can be visualized using showOneDMR
# showOneDMR(dmregions.CompareA4[1:6,1], dss.BSobj)
# 
# #use metilene for de novo DMR identification
# # Use pathfindR for network analysis
# 
# 
# #Prepare DSS object
# perc.input.list <- list()
# for (i in seq(5,ncol(merge.norm.filt.objDB_data),3)) {
#   #print(class(as.numeric(i)))
#   #print(head(merge.norm.filt.objDB_data[,c(1:2, i, (i + 1))]))
#   temp3.df <- merge.norm.filt.objDB_data[,c(1:2, i, (i + 1))]
#   colnames(temp3.df) <- c("chr", "pos", "N", "X")
#   temp3.df["perc"] <- temp3.df$X / temp3.df$N
#   temp3.id <- unlist(lapply(str_split(colnames(merge.norm.filt.objDB_data[,c(1:2, i, (i + 1))])[3], "_"), function(x) x[1]))
#   print(temp3.id)
#   perc.input.list[[temp3.id]] <- temp3.df
#   #assign(paste0(temp3.id, "_dss.df"), temp3.df)
#   rm(temp3.df)
#   rm(temp3.id)
# }
# 
# perc.input.list.df <- do.call(cbind, perc.input.list)
# perc.input.list.df.filt <- perc.input.list.df[,grepl(".perc", colnames(perc.input.list.df))]
# EDASeq::plotPCA(as.matrix(perc.input.list.df.filt[,c(1,2,3,6,7)]))
# #For some situations, it might be desirable to summarize methylation information over tiling windows rather than doing base-pair resolution analysis.
# #Tiling windows analysis
# myobjDB.low <- methRead(mysamples_list,
#                  sample.id=samples.id,
#                  assembly="GRCm39",
#                  treatment=rep(0, length(mysamples_list)),
#                  context="CpG",
#                  mincov = 3)
# 
# tiles = tileMethylCounts(myobjDB.low,win.size=1000,step.size=1000,cov.bases = 10)
# head(tiles[[1]],3)

# read the files to a methylRawList object: myobj
myobjDB <- methRead(mysamples_list,
                    sample.id=samples.id,
                    assembly="GRCm39",
                    treatment=rep(0,16),
                    context="CpG",
                    mincov = 5)


png("myobjDB_getMethylationStats.png", height = 1000, width = 1500)
par(mfrow=c(4,4))
for (i in 1:nFiles) {
  print(getSampleID(myobjDB[[i]]))
  getMethylationStats(myobjDB[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
}
dev.off()

#Individual CpGs list need to be prepared 
# Generate and save histograms showing CpG Methylation Coverage
png("myobjDB_getCoverageStats.png", height = 1000, width = 1500)
par(mfrow=c(4,4))
for (i in 1:nFiles) {
  print(getSampleID(myobjDB[[i]]))
  getCoverageStats(myobjDB[[i]], plot = TRUE, both.strands = FALSE) #Get %CpG methylation information
}
dev.off()



#Filtering samples based on read coverage
filtered.myobjDB <- filterByCoverage(myobjDB, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)


#Comparative analysis
#Merging samples
#Note: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting
merge.filt.objDB <- unite(filtered.myobjDB, destrand=FALSE)
head(merge.filt.objDB)
rm(filtered.myobjDB)
rm(myobjDB)
gc()

#Sample Correlation
png("merge.filt.objDB_getCorrelation.png", height = 2000, width = 2000)
getCorrelation(merge.filt.objDB,plot=TRUE)
dev.off()

# Clustering dendrogram
png("merge.filt.objDB_dendrogram_path.png", height = 600, width = 1000)
clusterSamples(merge.filt.objDB, dist="correlation", method="ward.D", plot=TRUE)
dev.off()

# Run a PCA analysis on percent methylation for all samples
png("pca_path.png", height = 1000, width = 1000)
PCASamples(merge.filt.objDB)
dev.off()

#CompareA4

# read the files to a methylRawList object: myobj
myobjDB.CompareA4 <- methRead(mysamples_list[c(6:10)],
                              sample.id=samples.id[c(6:10)],
                              assembly="GRCm39",
                              treatment=c(0,0,1,1,1),
                              context="CpG",
                              mincov = 5)


#Filtering samples based on read coverage
filtered.myobjDB.CompareA4 <- filterByCoverage(myobjDB.CompareA4, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)


#Comparative analysis
#Merging samples
#Note: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting
merge.filt.objDB.CompareA4 <- unite(filtered.myobjDB.CompareA4, destrand=FALSE)
head(merge.filt.objDB.CompareA4)


#Finding differentially methylated bases or regions
myDiff.CompareA4=calculateDiffMeth(merge.filt.objDB.CompareA4,mc.cores=2)
myDiff25p.hyper.CompareA4=getMethylDiff(myDiff.CompareA4,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareA4=getMethylDiff(myDiff.CompareA4,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareA4=getMethylDiff(myDiff.CompareA4,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareA4,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)
myDiff.CompareA4_df <- getData(myDiff.CompareA4)
rm(filtered.myobjDB.CompareA4)
rm(myobjDB.CompareA4)
rm(merge.filt.objDB.CompareA4)
gc()

#CompareB1

# read the files to a methylRawList object: myobj
myobjDB.CompareB1 <- methRead(mysamples_list[c(11,15,16,12,13)],
                              sample.id=samples.id[c(11,15,16,12,13)],
                              assembly="GRCm39",
                              treatment=c(0,0,0,1,1),
                              context="CpG",
                              mincov = 5)


#Filtering samples based on read coverage
filtered.myobjDB.CompareB1 <- filterByCoverage(myobjDB.CompareB1, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

#Comparative analysis
#Merging samples
#Note: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting
merge.filt.objDB.CompareB1 <- unite(filtered.myobjDB.CompareB1, destrand=FALSE)
head(merge.filt.objDB.CompareB1)


#Finding differentially methylated bases or regions
myDiff.CompareB1=calculateDiffMeth(merge.filt.objDB.CompareB1,mc.cores=2)
myDiff25p.hyper.CompareB1=getMethylDiff(myDiff.CompareB1,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareB1=getMethylDiff(myDiff.CompareB1,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareB1=getMethylDiff(myDiff.CompareB1,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareB1,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)
rm(filtered.myobjDB.CompareB1)
rm(myobjDB.CompareB1)
rm(merge.filt.objDB.CompareB1)
gc()


#Compare B2

# read the files to a methylRawList object: myobj
myobjDB.CompareB2 <- methRead(mysamples_list[c(11,15,16,1,2)],
                 sample.id=samples.id[c(11,15,16,1,2)],
                 assembly="GRCm39",
                 treatment=c(0,0,0,1,1),
                 context="CpG",
                 mincov = 5)


#Filtering samples based on read coverage
filtered.myobjDB.CompareB2 <- filterByCoverage(myobjDB.CompareB2, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
rm(myobjDB.CompareB2)
gc()

#Comparative analysis
#Merging samples
#Note: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting
merge.filt.objDB.CompareB2 <- unite(filtered.myobjDB.CompareB2, destrand=FALSE)
head(merge.filt.objDB.CompareB2)


#Finding differentially methylated bases or regions
myDiff.CompareB2=calculateDiffMeth(merge.filt.objDB.CompareB2,mc.cores=2)
myDiff25p.hyper.CompareB2=getMethylDiff(myDiff.CompareB2,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareB2=getMethylDiff(myDiff.CompareB2,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareB2=getMethylDiff(myDiff.CompareB2,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareB2,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
rm(filtered.myobjDB.CompareB2)
rm(myobjDB.CompareB2)
rm(merge.filt.objDB.CompareB2)
gc()

#CompareB3

# read the files to a methylRawList object: myobj
myobjDB.CompareB3 <- methRead(mysamples_list[c(11,15,16,4,5)],
                              sample.id=samples.id[c(11,15,16,4,5)],
                              assembly="GRCm39",
                              treatment=c(0,0,0,1,1),
                              context="CpG",
                              mincov = 5)


#Filtering samples based on read coverage
filtered.myobjDB.CompareB3 <- filterByCoverage(myobjDB.CompareB3, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)


#Comparative analysis
#Merging samples
#Note: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting
merge.filt.objDB.CompareB3 <- unite(filtered.myobjDB.CompareB3, destrand=FALSE)
head(merge.filt.objDB.CompareB3)


#Finding differentially methylated bases or regions
myDiff.CompareB3=calculateDiffMeth(merge.filt.objDB.CompareB3,mc.cores=2)
myDiff25p.hyper.CompareB3=getMethylDiff(myDiff.CompareB3,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareB3=getMethylDiff(myDiff.CompareB3,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareB3=getMethylDiff(myDiff.CompareB3,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareB3,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)
rm(filtered.myobjDB.CompareB3)
rm(myobjDB.CompareB3)
rm(merge.filt.objDB.CompareB3)
gc()

#CompareB4

# read the files to a methylRawList object: myobj
myobjDB.CompareB4 <- methRead(mysamples_list[c(11,15,16,8,9,10)],
                              sample.id=samples.id[c(11,15,16,8,9,10)],
                              assembly="GRCm39",
                              treatment=c(0,0,0,1,1,1),
                              context="CpG",
                              mincov = 5)


#Filtering samples based on read coverage
filtered.myobjDB.CompareB4 <- filterByCoverage(myobjDB.CompareB4, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)


#Comparative analysis
#Merging samples
#Note: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting
merge.filt.objDB.CompareB4 <- unite(filtered.myobjDB.CompareB4, destrand=FALSE)
head(merge.filt.objDB.CompareB4)


#Finding differentially methylated bases or regions
myDiff.CompareB4=calculateDiffMeth(merge.filt.objDB.CompareB4,mc.cores=2)
myDiff25p.hyper.CompareB4=getMethylDiff(myDiff.CompareB4,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareB4=getMethylDiff(myDiff.CompareB4,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareB4=getMethylDiff(myDiff.CompareB4,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareB4,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)
rm(filtered.myobjDB.CompareB4)
rm(myobjDB.CompareB4)
rm(merge.filt.objDB.CompareB4)
gc()

#CompareB5

# read the files to a methylRawList object: myobj
myobjDB.CompareB5 <- methRead(mysamples_list[c(11,15,16,6,7)],
                              sample.id=samples.id[c(11,15,16,6,7)],
                              assembly="GRCm39",
                              treatment=c(0,0,0,1,1),
                              context="CpG",
                              mincov = 5)


#Filtering samples based on read coverage
filtered.myobjDB.CompareB5 <- filterByCoverage(myobjDB.CompareB5, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
rm(myobjDB.CompareB5)
gc()

#Comparative analysis
#Merging samples
#Note: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting
merge.filt.objDB.CompareB5 <- unite(filtered.myobjDB.CompareB5, destrand=FALSE)
head(merge.filt.objDB.CompareB5)


#Finding differentially methylated bases or regions
myDiff.CompareB5=calculateDiffMeth(merge.filt.objDB.CompareB5,mc.cores=2)
myDiff25p.hyper.CompareB5=getMethylDiff(myDiff.CompareB5,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareB5=getMethylDiff(myDiff.CompareB5,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareB5=getMethylDiff(myDiff.CompareB5,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareB5,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)
rm(filtered.myobjDB.CompareB5)
rm(myobjDB.CompareB5)
rm(merge.filt.objDB.CompareB5)
gc()


#CompareC1

# read the files to a methylRawList object: myobj
myobjDB.CompareC1 <- methRead(mysamples_list[c(12,13,8,9,10)],
                              sample.id=samples.id[c(12,13,8,9,10)],
                              assembly="GRCm39",
                              treatment=c(0,0,1,1,1),
                              context="CpG",
                              mincov = 5)


#Filtering samples based on read coverage
filtered.myobjDB.CompareC1 <- filterByCoverage(myobjDB.CompareC1, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
rm(myobjDB.CompareC1)
gc()

#Comparative analysis
#Merging samples
#Note: In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting
merge.filt.objDB.CompareC1 <- unite(filtered.myobjDB.CompareC1, destrand=FALSE)
head(merge.filt.objDB.CompareC1)


#Finding differentially methylated bases or regions
myDiff.CompareC1=calculateDiffMeth(merge.filt.objDB.CompareC1,mc.cores=2)
myDiff25p.hyper.CompareC1=getMethylDiff(myDiff.CompareC1,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareC1=getMethylDiff(myDiff.CompareC1,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareC1=getMethylDiff(myDiff.CompareC1,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareC1,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)
rm(filtered.myobjDB.CompareC1)
rm(myobjDB.CompareC1)
rm(merge.filt.objDB.CompareC1)
gc()

getwd()
#qvalue distribution
summary(myDiff.CompareB1$qvalue)
summary(myDiff.CompareB2$qvalue)
summary(myDiff.CompareB3$qvalue)
summary(myDiff.CompareB4$qvalue)
summary(myDiff.CompareB5$qvalue)

#Export data
myDiff25p.CompareB1.df <- getData(myDiff25p.CompareB1)
write.table(myDiff25p.CompareB1.df, 'myDiff25p.CompareB1.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

myDiff25p.CompareB2.df <- getData(myDiff25p.CompareB2)
write.table(myDiff25p.CompareB2.df, 'myDiff25p.CompareB2.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

myDiff25p.CompareB3.df <- getData(myDiff25p.CompareB3)
write.table(myDiff25p.CompareB3.df, 'myDiff25p.CompareB3.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

myDiff25p.CompareB4.df <- getData(myDiff25p.CompareB4)
write.table(myDiff25p.CompareB4.df, 'myDiff25p.CompareB4.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

myDiff25p.CompareB5.df <- getData(myDiff25p.CompareB5)
write.table(myDiff25p.CompareB5.df, 'myDiff25p.CompareB5.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

myDiff25p.CompareA4.df <- getData(myDiff25p.CompareA4)
write.table(myDiff25p.CompareA4.df, 'myDiff25p.CompareA4.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)


# cut off filter only meth diff 
methdiff.cutoff <- 25
qval <- 0.01

for (comparison in c("CompareB1", "CompareB2", "CompareB3", "CompareB4","CompareB5", "CompareA4", "CompareC1")){
  print(paste0("myDiff.",comparison))
  temp6_sampleid <- paste0("myDiff.",comparison)
  temp6_sampledf <- get(temp6_sampleid)
  print(paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison))
  assign(paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison),getMethylDiff(temp6_sampledf,difference=methdiff.cutoff,qvalue=qval))
  temp6.diffmeth <- get(paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison))
  write.table(temp6.diffmeth, paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison,".df.txt"), quote = F, sep="\t", append = F, row.names = F, col.names = F)
}

#Individual Tiles list need to be prepared 
#coverage can be decreased to 3

# #Tiling windows analysis
myobjDB.low <- methRead(mysamples_list,
                 sample.id=samples.id,
                 assembly="GRCm39",
                 treatment=rep(0, length(mysamples_list)),
                 context="CpG",
                 mincov = 3)

tiles = tileMethylCounts(myobjDB.low,win.size=1000,step.size=1000,cov.bases = 5)
head(tiles[[1]],3)
myobjDB.low.tiles <- unite(tiles, destrand=FALSE)
png("myobjDB.low.tiles_PCA.png", height = 480, width = 480)
PCASamples(myobjDB.low.tiles)
dev.off()
myobjDB.low.tiles_reorg =reorganize(myobjDB.low.tiles,sample.ids=c("S3","S4","S5","S8","S9","S19","S20","S25","S26","S28"),treatment=rep(0, length(c("S3","S4","S5","S8","S9","S19","S20","S25","S26","S28"))))
png("myobjDB.low.tiles_reorg_PCA.png", height = 480, width = 480)
PCASamples(myobjDB.low.tiles_reorg)
dev.off()



#CompareTileA4

# read the files to a methylRawList object: myobj
myobjDB.low.CompareA4 <- methRead(mysamples_list[c(6:10)],
                                  sample.id=samples.id[c(6:10)],
                                  assembly="GRCm39",
                                  treatment=c(0,0,1,1,1),
                                  context="CpG",
                                  mincov = 3)


merge.myobjDB.low.CompareA4  <- unite(myobjDB.low.CompareA4, destrand=FALSE)

merge.myobjDB.low.CompareTileA4 = tileMethylCounts(merge.myobjDB.low.CompareA4,win.size=1000,step.size=1000,cov.bases = 5)


#Finding differentially methylated bases or regions
myDiff.CompareTileA4=calculateDiffMeth(merge.myobjDB.low.CompareTileA4,mc.cores=2)
myDiff25p.hyper.CompareTileA4=getMethylDiff(myDiff.CompareTileA4,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareTileA4=getMethylDiff(myDiff.CompareTileA4,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareTileA4=getMethylDiff(myDiff.CompareTileA4,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareTileA4,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)

rm(myobjDB.low.CompareA4)
gc()


#CompareTileB1

# read the files to a methylRawList object: myobj
myobjDB.low.CompareB1 <- methRead(mysamples_list[c(11,15,16,12,13)],
                                  sample.id=samples.id[c(11,15,16,12,13)],
                                  assembly="GRCm39",
                                  treatment=c(0,0,0,1,1),
                                  context="CpG",
                                  mincov = 3)

merge.myobjDB.low.CompareB1  <- unite(myobjDB.low.CompareB1, destrand=FALSE)

merge.myobjDB.low.CompareTileB1 = tileMethylCounts(merge.myobjDB.low.CompareB1,win.size=1000,step.size=1000,cov.bases = 5)


#Finding differentially methylated bases or regions
myDiff.CompareTileB1=calculateDiffMeth(merge.myobjDB.low.CompareTileB1,mc.cores=2)
myDiff25p.hyper.CompareTileB1=getMethylDiff(myDiff.CompareTileB1,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareTileB1=getMethylDiff(myDiff.CompareTileB1,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareTileB1=getMethylDiff(myDiff.CompareTileB1,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareTileB1,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)

rm(myobjDB.low.CompareB1)
gc()

#CompareTile B2

# read the files to a methylRawList object: myobj
myobjDB.low.CompareB2 <- methRead(mysamples_list[c(11,15,16,1,2)],
                                  sample.id=samples.id[c(11,15,16,1,2)],
                                  assembly="GRCm39",
                                  treatment=c(0,0,0,1,1),
                                  context="CpG",
                                  mincov = 3)
merge.myobjDB.low.CompareB2  <- unite(myobjDB.low.CompareB2, destrand=FALSE)

merge.myobjDB.low.CompareTileB2 = tileMethylCounts(merge.myobjDB.low.CompareB2,win.size=1000,step.size=1000,cov.bases = 5)


#Finding differentially methylated bases or regions
myDiff.CompareTileB2=calculateDiffMeth(merge.myobjDB.low.CompareTileB2,mc.cores=2)
myDiff25p.hyper.CompareTileB2=getMethylDiff(myDiff.CompareTileB2,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareTileB2=getMethylDiff(myDiff.CompareTileB2,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareTileB2=getMethylDiff(myDiff.CompareTileB2,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareTileB2,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)

rm(myobjDB.low.CompareB2)
gc()


#CompareTileB3

# read the files to a methylRawList object: myobj
myobjDB.low.CompareB3 <- methRead(mysamples_list[c(11,15,16,4,5)],
                                  sample.id=samples.id[c(11,15,16,4,5)],
                                  assembly="GRCm39",
                                  treatment=c(0,0,0,1,1),
                                  context="CpG",
                                  mincov = 3)


merge.myobjDB.low.CompareB3  <- unite(myobjDB.low.CompareB3, destrand=FALSE)

merge.myobjDB.low.CompareTileB3 = tileMethylCounts(merge.myobjDB.low.CompareB3,win.size=1000,step.size=1000,cov.bases = 5)


#Finding differentially methylated bases or regions
myDiff.CompareTileB3=calculateDiffMeth(merge.myobjDB.low.CompareTileB3,mc.cores=2)
myDiff25p.hyper.CompareTileB3=getMethylDiff(myDiff.CompareTileB3,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareTileB3=getMethylDiff(myDiff.CompareTileB3,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareTileB3=getMethylDiff(myDiff.CompareTileB3,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareTileB3,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)

rm(myobjDB.low.CompareB3)
gc()

#CompareTileB4

# read the files to a methylRawList object: myobj
myobjDB.low.CompareB4 <- methRead(mysamples_list[c(11,15,16,8,9,10)],
                                  sample.id=samples.id[c(11,15,16,8,9,10)],
                                  assembly="GRCm39",
                                  treatment=c(0,0,0,1,1,1),
                                  context="CpG",
                                  mincov = 3)


merge.myobjDB.low.CompareB4  <- unite(myobjDB.low.CompareB4, destrand=FALSE)

merge.myobjDB.low.CompareTileB4 = tileMethylCounts(merge.myobjDB.low.CompareB4,win.size=1000,step.size=1000,cov.bases = 5)


#Finding differentially methylated bases or regions
myDiff.CompareTileB4=calculateDiffMeth(merge.myobjDB.low.CompareTileB4,mc.cores=2)
myDiff25p.hyper.CompareTileB4=getMethylDiff(myDiff.CompareTileB4,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareTileB4=getMethylDiff(myDiff.CompareTileB4,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareTileB4=getMethylDiff(myDiff.CompareTileB4,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareTileB4,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)

rm(myobjDB.low.CompareB4)
gc()


#CompareTileB5

# read the files to a methylRawList object: myobj
myobjDB.low.CompareB5 <- methRead(mysamples_list[c(11,15,16,6,7)],
                                  sample.id=samples.id[c(11,15,16,6,7)],
                                  assembly="GRCm39",
                                  treatment=c(0,0,0,1,1),
                                  context="CpG",
                                  mincov = 3)


merge.myobjDB.low.CompareB5  <- unite(myobjDB.low.CompareB5, destrand=FALSE)

merge.myobjDB.low.CompareTileB5 = tileMethylCounts(merge.myobjDB.low.CompareB5,win.size=1000,step.size=1000,cov.bases = 5)


#Finding differentially methylated bases or regions
myDiff.CompareTileB5=calculateDiffMeth(merge.myobjDB.low.CompareTileB5,mc.cores=2)
myDiff25p.hyper.CompareTileB5=getMethylDiff(myDiff.CompareTileB5,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareTileB5=getMethylDiff(myDiff.CompareTileB5,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareTileB5=getMethylDiff(myDiff.CompareTileB5,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareTileB5,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)

rm(myobjDB.low.CompareB5)
gc()

#CompareTileC1

# read the files to a methylRawList object: myobj
myobjDB.low.CompareC1 <- methRead(mysamples_list[c(12,13,8,9,10)],
                                  sample.id=samples.id[c(12,13,8,9,10)],
                                  assembly="GRCm39",
                                  treatment=c(0,0,1,1,1),
                                  context="CpG",
                                  mincov = 3)


merge.myobjDB.low.CompareC1  <- unite(myobjDB.low.CompareC1, destrand=FALSE)

merge.myobjDB.low.CompareTileC1 = tileMethylCounts(merge.myobjDB.low.CompareC1,win.size=1000,step.size=1000,cov.bases = 5)


#Finding differentially methylated bases or regions
myDiff.CompareTileC1=calculateDiffMeth(merge.myobjDB.low.CompareTileC1,mc.cores=2)
myDiff25p.hyper.CompareTileC1=getMethylDiff(myDiff.CompareTileC1,difference=25,qvalue=0.01,type="hyper")
myDiff25p.hypo.CompareTileC1=getMethylDiff(myDiff.CompareTileC1,difference=25,qvalue=0.01,type="hypo")
myDiff25p.CompareTileC1=getMethylDiff(myDiff.CompareTileC1,difference=25,qvalue=0.01)
diffMethPerChr(myDiff.CompareTileC1,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)

rm(myobjDB.low.CompareC1)
gc()


#qvalue distribution
summary(myDiff.CompareTileB1$qvalue)
summary(myDiff.CompareTileB2$qvalue)
summary(myDiff.CompareTileB3$qvalue)
summary(myDiff.CompareTileB4$qvalue)
summary(myDiff.CompareTileB5$qvalue)
dim(myDiff25p.CompareTileB1)
dim(merge.filt.objDB.CompareTileB1)

#Export data
myDiff25p.CompareTileB1.df <- getData(myDiff25p.CompareTileB1)
write.table(myDiff25p.CompareTileB1.df, 'myDiff25p.CompareTileB1.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

myDiff25p.CompareTileB2.df <- getData(myDiff25p.CompareTileB2)
write.table(myDiff25p.CompareTileB2.df, 'myDiff25p.CompareTileB2.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

myDiff25p.CompareTileB3.df <- getData(myDiff25p.CompareTileB3)
write.table(myDiff25p.CompareTileB3.df, 'myDiff25p.CompareTileB3.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

myDiff25p.CompareTileB4.df <- getData(myDiff25p.CompareTileB4)
write.table(myDiff25p.CompareTileB4.df, 'myDiff25p.CompareTileB4.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

myDiff25p.CompareTileB5.df <- getData(myDiff25p.CompareTileB5)
write.table(myDiff25p.CompareTileB5.df, 'myDiff25p.CompareTileB5.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

myDiff25p.CompareTileA4.df <- getData(myDiff25p.CompareTileA4)
write.table(myDiff25p.CompareTileA4.df, 'myDiff25p.CompareTileA4.df.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)



methdiff.cutoff <- 25
qval <- 0.01
for (comparison in c("CompareTileB1", "CompareTileB2", "CompareTileB3", "CompareTileB4","CompareTileB5", "CompareTileA4")){
  print(paste0("myDiff.",comparison))
  temp6_sampleid <- paste0("myDiff.",comparison)
  temp6_sampledf <- get(temp6_sampleid)
  print(paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison))
  assign(paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison),getMethylDiff(temp6_sampledf,difference=methdiff.cutoff,qvalue=qval))
  temp6.diffmeth <- get(paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison))
  write.table(temp6.diffmeth, paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison,".df.txt"), quote = F, sep="\t", append = F, row.names = F, col.names = F)
}

methdiff.cutoff <- 25
qval <- 0.01
for (comparison in c("CompareTileC1")){
  print(paste0("myDiff.",comparison))
  temp6_sampleid <- paste0("myDiff.",comparison)
  temp6_sampledf <- get(temp6_sampleid)
  print(paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison))
  assign(paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison),getMethylDiff(temp6_sampledf,difference=methdiff.cutoff,qvalue=qval))
  temp6.diffmeth <- get(paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison))
  write.table(temp6.diffmeth, paste0("myDiff",methdiff.cutoff,"p.", qval, ".",comparison,".df.txt"), quote = F, sep="\t", append = F, row.names = F, col.names = F)
}

ls *.df.txt -1 | awk '{print "grep chr "$1" > "$1"_chr.txt"}' > get_chr_files.sh
chmod +x get_chr_files.sh
./get_chr_files.sh
# #Prepare DSS object
library(DSS)
require(bsseq)
library(stringr)
dss.input.list <- list()
merge.filt.objDB_df <- getData(merge.filt.objDB)
colnames(merge.filt.objDB_df) <- c(colnames(merge.filt.objDB_df)[1:4], 
                                           paste0(unlist(lapply(getSampleID(merge.filt.objDB), function(x) rep(x,3))), "_", 
                                                  colnames(getData(merge.filt.objDB)[5:52])))
 
head(merge.filt.objDB_df)
 

for (i in seq(5,ncol(merge.filt.objDB_df),3)) {
  #print(class(as.numeric(i)))
  #print(head(merge.norm.filt.objDB_data[,c(1:2, i, (i + 1))]))
  temp2.df <- merge.filt.objDB_df[,c(1:2, i, (i + 1))]
  colnames(temp2.df) <- c("chr", "pos", "N", "X")
  temp2.id <- unlist(lapply(str_split(colnames(merge.filt.objDB_df[,c(1:2, i, (i + 1))])[3], "_"), function(x) x[1]))
  print(temp2.id)
  dss.input.list[[temp2.id]] <- temp2.df
  #assign(paste0(temp2.id, "_dss.df"), temp2.df)
  rm(temp2.df)
  rm(temp2.id)
}

dss.BSobj = makeBSseqData(dss.input.list, names(dss.input.list))
dss.BSobj
dmlTest.sm.CompareA4 <- DMLtest(dss.BSobj, group1=c("S19", "S20"), group2=c("S25", "S26","S28"),
                     smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.CompareA4 = callDML(dmlTest.sm.CompareA4, p.threshold=0.001)
head(dmloci.CompareA4)

dmregions.CompareA4 = callDMR(dmlTest.sm.CompareA4, p.threshold=0.01)
head(dmregions.CompareA4)

dmlTest.sm.CompareB1 <- DMLtest(dss.BSobj, group1=c("S3","S8", "S9"), group2=c("S4", "S5"),
                                smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.CompareB1 = callDML(dmlTest.sm.CompareB1, p.threshold=0.001)
head(dmloci.CompareB1)
dim(dmloci.CompareB1)
dmregions.CompareB1 = callDMR(dmlTest.sm.CompareB1, p.threshold=0.01)
head(dmregions.CompareB1)


dmlTest.sm.CompareB2 <- DMLtest(dss.BSobj, group1=c("S3", "S8", "S9"), group2=c("S10", "S11"),
                                smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.CompareB2 = callDML(dmlTest.sm.CompareB2, p.threshold=0.001)
head(dmloci.CompareB2)

dmregions.CompareB2 = callDMR(dmlTest.sm.CompareB2, p.threshold=0.01)
head(dmregions.CompareB2)


dmlTest.sm.CompareB3 <- DMLtest(dss.BSobj, group1=c("S3","S8", "S9"), group2=c("S17", "S18"),
                                smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.CompareB3 = callDML(dmlTest.sm.CompareB3, p.threshold=0.001)
head(dmloci.CompareB3)

dmregions.CompareB3 = callDMR(dmlTest.sm.CompareB3, p.threshold=0.01)
head(dmregions.CompareB3)

dmlTest.sm.CompareB4 <- DMLtest(dss.BSobj, group1=c("S3", "S8", "S9"), group2=c("S25", "S26","S28"),
                                smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.CompareB4 = callDML(dmlTest.sm.CompareB4, p.threshold=0.001)
head(dmloci.CompareB4)

dmregions.CompareB4 = callDMR(dmlTest.sm.CompareB4, p.threshold=0.01)
head(dmregions.CompareB4)


dmlTest.sm.CompareB5 <- DMLtest(dss.BSobj, group1=c("S3", "S8", "S9"), group2=c("S19", "S20"),
                                smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.CompareB5 = callDML(dmlTest.sm.CompareB5, p.threshold=0.001)
head(dmloci.CompareB5)

dmregions.CompareB5 = callDMR(dmlTest.sm.CompareB5, p.threshold=0.01)
head(dmregions.CompareB5)

dim(dmloci.CompareB1)
dim(dmloci.CompareB2)
dim(dmloci.CompareB3)
dim(dmloci.CompareB4)
dim(dmloci.CompareB5)
dim(dmregions.CompareB1)
dim(dmregions.CompareB2)
dim(dmregions.CompareB3)
dim(dmregions.CompareB4)
dim(dmregions.CompareB5)
write.table(dmregions.CompareB1, 'dmregions.CompareB1.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmregions.CompareB2, 'dmregions.CompareB2.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmregions.CompareB3, 'dmregions.CompareB3.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmregions.CompareB4, 'dmregions.CompareB4.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmregions.CompareB5, 'dmregions.CompareB5.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

write.table(dmloci.CompareB1, 'dmloci.CompareB1.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmloci.CompareB2, 'dmloci.CompareB2.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmloci.CompareB3, 'dmloci.CompareB3.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmloci.CompareB4, 'dmloci.CompareB4.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmloci.CompareB5, 'dmloci.CompareB5.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

rm(merge.filt.objDB)
rm(merge.filt.objDB.CompareA4)
rm(merge.filt.objDB.CompareB1)
rm(merge.filt.objDB.CompareB2)
rm(merge.filt.objDB.CompareB3)
rm(merge.filt.objDB.CompareB4)
rm(merge.filt.objDB.CompareB5)
rm(merge.myobjDB.low.CompareTileA4)
rm(merge.myobjDB.low.CompareTileB1)
rm(merge.myobjDB.low.CompareTileB2)
rm(merge.myobjDB.low.CompareTileB3)
rm(merge.myobjDB.low.CompareTileB4)
rm(merge.myobjDB.low.CompareTileB5)
rm(merge.myobjDB.low.CompareA4)
rm(merge.myobjDB.low.CompareB1)
rm(merge.myobjDB.low.CompareB2)
rm(merge.myobjDB.low.CompareB3)
rm(merge.myobjDB.low.CompareB4)
rm(merge.myobjDB.low.CompareB5)



#Tiles based DSS analysis
dss.tiles.input.list <- list()
meth.tiles.input.list <- list()

myobjDB.low.tiles_df <- getData(myobjDB.low.tiles)
colnames(myobjDB.low.tiles_df) <- c(colnames(myobjDB.low.tiles_df)[1:4], 
                                    paste0(unlist(lapply(getSampleID(myobjDB.low.tiles), function(x) rep(x,3))), "_", 
                                           colnames(getData(myobjDB.low.tiles)[5:52])))

head(myobjDB.low.tiles_df)


for (i in seq(5,ncol(myobjDB.low.tiles_df),3)) {
  #print(class(as.numeric(i)))
  #print(head(merge.norm.filt.objDB_data[,c(1:2, i, (i + 1))]))
  temp4.df <- myobjDB.low.tiles_df[,c(1:2, i, (i + 1))]
  colnames(temp4.df) <- c("chr", "pos", "N", "X")
  temp4.id <- unlist(lapply(str_split(colnames(myobjDB.low.tiles_df[,c(1:2, i, (i + 1))])[3], "_"), function(x) x[1]))
  print(temp4.id)
  dss.tiles.input.list[[temp4.id]] <- temp4.df
  temp4.df["meth"] <- temp4.df$X/temp4.df$N
  meth.tiles.input.list[[temp4.id]]  <- temp4.df
  #assign(paste0(temp4.id, "_dss.tiles.df"), temp4.df)
  rm(temp4.df)
  rm(temp4.id)
}

dss.tiles.BSobj = makeBSseqData(dss.tiles.input.list, names(dss.tiles.input.list))
dss.tiles.BSobj

meth.tiles.input.list_df <- do.call(cbind.data.frame, meth.tiles.input.list)
rownames(meth.tiles.input.list_df) <- paste0(meth.tiles.input.list_df$S10.chr,"%",meth.tiles.input.list_df$S10.pos)
meth.tiles.input.list_df <- meth.tiles.input.list_df[,c(seq(5,80,5))]
meth.tiles.input.list_df_st <- stack(as.matrix(meth.tiles.input.list_df))
head(meth.tiles.input.list_df_st)
boxplot(meth.tiles.input.list_df)
#boxplot_meth.tiles.input.list_df.pdf

dss.tiles_reorg.input.list <- list()
meth.tiles_reorg.input.list <- list()

myobjDB.low.tiles_reorg_df <- getData(myobjDB.low.tiles_reorg)
colnames(myobjDB.low.tiles_reorg_df) <- c(colnames(myobjDB.low.tiles_reorg_df)[1:4], 
                                          paste0(unlist(lapply(getSampleID(myobjDB.low.tiles_reorg), function(x) rep(x,3))), "_", 
                                                 colnames(getData(myobjDB.low.tiles_reorg)[5:34])))

head(myobjDB.low.tiles_reorg_df)

for (i in seq(5,ncol(myobjDB.low.tiles_reorg_df),3)) {
  #print(class(as.numeric(i)))
  #print(head(merge.norm.filt.objDB_data[,c(1:2, i, (i + 1))]))
  temp5.df <- myobjDB.low.tiles_reorg_df[,c(1:2, i, (i + 1))]
  colnames(temp5.df) <- c("chr", "pos", "N", "X")
  temp5.id <- unlist(lapply(str_split(colnames(myobjDB.low.tiles_reorg_df[,c(1:2, i, (i + 1))])[3], "_"), function(x) x[1]))
  print(temp5.id)
  dss.tiles_reorg.input.list[[temp5.id]] <- temp5.df
  temp5.df["meth"] <- temp5.df$X/temp5.df$N
  meth.tiles_reorg.input.list[[temp5.id]]  <- temp5.df
  #assign(paste0(temp5.id, "_dss.tiles_reorg.df"), temp5.df)
  rm(temp5.df)
  rm(temp5.id)
}


meth.tiles_reorg.input.list_df <- do.call(cbind.data.frame, meth.tiles_reorg.input.list)
#only first 2 column use for colnmaes
rownames(meth.tiles_reorg.input.list_df) <- paste0(meth.tiles_reorg.input.list_df$S3.chr,"%",meth.tiles_reorg.input.list_df$S3.pos)
meth.tiles_reorg.input.list_df <- meth.tiles_reorg.input.list_df[,c(seq(5,50,5))]
boxplot(meth.tiles_reorg.input.list_df)
#boxplot_meth.tiles_reorg.input.list_df.pdf

dmlTest.sm.tiles.CompareA4 <- DMLtest(dss.tiles.BSobj, group1=c("S19", "S20"), group2=c("S25", "S26","S28"),
                                      smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.tiles.CompareA4 = callDML(dmlTest.sm.tiles.CompareA4, p.threshold=0.001)
head(dmloci.tiles.CompareA4)

dmregions.tiles.CompareA4 = callDMR(dmlTest.sm.tiles.CompareA4, p.threshold=0.01)
head(dmregions.tiles.CompareA4)

dmlTest.sm.tiles.CompareB1 <- DMLtest(dss.tiles.BSobj, group1=c("S3","S8", "S9"), group2=c("S4", "S5"),
                                      smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.tiles.CompareB1 = callDML(dmlTest.sm.tiles.CompareB1, p.threshold=0.001)
head(dmloci.tiles.CompareB1)
dim(dmloci.tiles.CompareB1)
dmregions.tiles.CompareB1 = callDMR(dmlTest.sm.tiles.CompareB1, p.threshold=0.01)
head(dmregions.tiles.CompareB1)


dmlTest.sm.tiles.CompareB2 <- DMLtest(dss.tiles.BSobj, group1=c("S3", "S8", "S9"), group2=c("S10", "S11"),
                                      smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.tiles.CompareB2 = callDML(dmlTest.sm.tiles.CompareB2, p.threshold=0.001)
head(dmloci.tiles.CompareB2)

dmregions.tiles.CompareB2 = callDMR(dmlTest.sm.tiles.CompareB2, p.threshold=0.01)
head(dmregions.tiles.CompareB2)


dmlTest.sm.tiles.CompareB3 <- DMLtest(dss.tiles.BSobj, group1=c("S3","S8", "S9"), group2=c("S17", "S18"),
                                      smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.tiles.CompareB3 = callDML(dmlTest.sm.tiles.CompareB3, p.threshold=0.001)
head(dmloci.tiles.CompareB3)

dmregions.tiles.CompareB3 = callDMR(dmlTest.sm.tiles.CompareB3, p.threshold=0.01)
head(dmregions.tiles.CompareB3)

dmlTest.sm.tiles.CompareB4 <- DMLtest(dss.tiles.BSobj, group1=c("S3", "S8", "S9"), group2=c("S25", "S26","S28"),
                                      smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.tiles.CompareB4 = callDML(dmlTest.sm.tiles.CompareB4, p.threshold=0.001)
head(dmloci.tiles.CompareB4)

dmregions.tiles.CompareB4 = callDMR(dmlTest.sm.tiles.CompareB4, p.threshold=0.01)
head(dmregions.tiles.CompareB4)


dmlTest.sm.tiles.CompareB5 <- DMLtest(dss.tiles.BSobj, group1=c("S3", "S8", "S9"), group2=c("S19", "S20"),
                                      smoothing=TRUE)
#Smoothing is recommended for WGBS data
dmloci.tiles.CompareB5 = callDML(dmlTest.sm.tiles.CompareB5, p.threshold=0.001)
head(dmloci.tiles.CompareB5)

dmregions.tiles.CompareB5 = callDMR(dmlTest.sm.tiles.CompareB5, p.threshold=0.01)
head(dmregions.tiles.CompareB5)

dim(dmloci.tiles.CompareB1)
dim(dmloci.tiles.CompareB2)
dim(dmloci.tiles.CompareB3)
dim(dmloci.tiles.CompareB4)
dim(dmloci.tiles.CompareB5)
dim(dmregions.tiles.CompareB1)
dim(dmregions.tiles.CompareB2)
dim(dmregions.tiles.CompareB3)
dim(dmregions.tiles.CompareB4)
dim(dmregions.tiles.CompareB5)
write.table(dmregions.tiles.CompareB1, 'dmregions.tiles.CompareB1.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmregions.tiles.CompareB2, 'dmregions.tiles.CompareB2.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmregions.tiles.CompareB3, 'dmregions.tiles.CompareB3.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmregions.tiles.CompareB4, 'dmregions.tiles.CompareB4.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmregions.tiles.CompareB5, 'dmregions.tiles.CompareB5.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

write.table(dmloci.tiles.CompareB1, 'dmloci.tiles.CompareB1.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmloci.tiles.CompareB2, 'dmloci.tiles.CompareB2.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmloci.tiles.CompareB3, 'dmloci.tiles.CompareB3.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmloci.tiles.CompareB4, 'dmloci.tiles.CompareB4.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)
write.table(dmloci.tiles.CompareB5, 'dmloci.tiles.CompareB5.txt', quote = F, sep="\t", append = F, row.names = F, col.names = F)

#Import RNA-Seq data
gene_gencode_mouse_grcm39 <- read.table("/mnt/home3/reid/av638/ref_av/gencode/mouse/gene_gencode_mouse_grcm39_out.txt")
colnames(gene_gencode_mouse_grcm39) <- c("chr", "start", "end", "strand", "type", "ensID", "GeneName")
#Read RNA-Seq csv files
rnaseq <-  list.files(path = "/mnt/home3/reid/av638/methylseq/novaseq_em_seq/Re_Differentially_expressed_genes_comparisons_of_all_groups/",pattern = "\\.csv$")
rnaseqlist <- list()
for (files in rnaseq){
  print(files)
  files_df <- read.csv(paste0("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/Re_Differentially_expressed_genes_comparisons_of_all_groups/", files))
  files_df_gene <- merge(files_df, gene_gencode_mouse_grcm39, by="GeneName")
  files_df_gene <- files_df_gene[,c((length(files_df) + 1) :length(files_df_gene),1:length(files_df))]
  files_df_gene["Comparison"] <- gsub(".csv","",gsub("Results_gc_length_corrected_","",files))
  files.id <- gsub(".csv", "", files)
  files_df_gene_de <- files_df_gene %>% filter((logFC > 1) | (logFC < -1) & (FDR < 0.05))
  print(dim(files_df_gene_de))
  rnaseqlist[[files.id]] <- files_df_gene
  rnaseqlist[[paste0(files.id,"de")]] <- files_df_gene_de
  write.table(files_df_gene_de,paste0("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/Re_Differentially_expressed_genes_comparisons_of_all_groups/",files, "_degene.txt"), quote = F, sep="\t", append = F, row.names = F, col.names = F)
}

#hector 32 genes
hector_32_dereg_genes <- read.table("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/hector_32_dereg_genes.txt")
colnames(hector_32_dereg_genes) <- "GeneName"
hector_32_dereg_genes_pos <- merge(gene_gencode_mouse_grcm39, hector_32_dereg_genes, by="GeneName")
hector_32_dereg_genes_pos <- hector_32_dereg_genes_pos[,c(2:7,1)]
write.table(hector_32_dereg_genes_pos, "/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/hector_32_dereg_genes_pos.txt", quote = F, sep="\t", append = F, row.names = F, col.names = F)

#Annotation dmr to nearest gene

loadGenome.pl -name GRCm39 -fasta GRCm39.primary_assembly.genome.fa -gtf GRCm39.gencode.vM30.primary_assembly.annotation.gtf  -org mouse -force -gid
for var in *.df.txt
do 
base=$(basename ${var} .df.txt)
echo ${base}
annotatePeaks.pl ${base}.df.txt GRCm39 -gid > ${base}.df.txt_ann.txt
done
save as annotation.sh

#Annotation
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
dmr_files <- list.files(path="/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit", pattern="*.df_chr.txt")
dmr_files_list = list()
for (file in dmr_files){
  #print(file)
  filenames <- gsub(".df_chr.txt","",file)
  print(filenames)
  print(file.size(paste0("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/", file)))
  if (file.size(paste0("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/", file)) > 0){
    print ("data present")
    dmr_files_list[[filenames]] <- paste0("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/", file)
  }
}
library(org.Mm.eg.db)
peakAnnoList <- lapply(dmr_files_list, annotatePeak, TxDb=txdb,
                       tssRegion=c(-1000, 1000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
#plotAnnoBar_peakAnnoList.pdf
dmr_tagMatrixList <- lapply(dmr_files_list, getTagMatrix, windows=promoter)
plotAvgProf(dmr_tagMatrixList$myDiff10p.0.01.CompareTileA4, xlim=c(-1000, 1000))
plotAvgProf(dmr_tagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=500, facet="row")

library(ChIPpeakAnno)
dmr_adj_files <- list.files(path="/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit", pattern="*.df_chr.adj.txt")
dmr_adj_files <- grep("Tile", dmr_adj_files, value = TRUE)   


dmr_adj_files_list = list()
for (file in dmr_adj_files){
  #print(file)
  filenames <- gsub(".df_chr.adj.txt","",file)
  print(filenames)
  print(file.size(paste0("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/", file)))
  if (file.size(paste0("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/", file)) > 0){
    print ("data present")
    dmr_adj_files_list[[filenames]] <- paste0("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/", file)
  }
}
dmr_adj_bed <- dmr_adj_files_list$myDiff25p.CompareTileB1
dmr_adj_bed_gr <- toGRanges(dmr_adj_bed,  header=FALSE) 


library(TxDb.Mmusculus.UCSC.mm39.refGene)
annoData <- toGRanges(TxDb.Mmusculus.UCSC.mm39.refGene)
annoData[1:2]

seqlevelsStyle(dmr_adj_bed_gr) <- seqlevelsStyle(annoData)

anno1 <- annotatePeakInBatch(dmr_adj_bed_gr, AnnotationData=annoData, 
                            output="overlapping", 
                            FeatureLocForDistance="TSS",
                            bindingRegion=c(-1, 1))
anno1$symbol <- xget(anno1$feature, org.Mm.egSYMBOL)
head(anno1)
write.table(data.frame(anno1), "anno1.txt", quote = F, sep = "\t", append = F, col.names = T, row.names = F)

total_dmr_called <- read.table("total_dmr_called.txt")
ggplot(data=total_dmr_called, aes(x=V1, y=V2)) +
  geom_bar(stat="identity") +coord_flip() + theme_classic() + ggtitle("barplot of DMRs")


#barplot_total_dmr_called.pdf
library(dplyr)
total_dmr_called_sel <- total_dmr_called %>% filter(grepl('A4|B4|B1|C1', V1))
total_dmr_called_sel <- total_dmr_called_sel %>% filter(grepl('Tile', V1))
ggplot(data=total_dmr_called_sel, aes(x=V1, y=V2)) +
  geom_bar(stat="identity") +coord_flip() + theme_classic() + ggtitle("barplot of DMRs")

#barplot_total_dmr_called_sel.pdf
detach(package:dplyr)

library(data.table)
gene_match_A4_dmr <- fread("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/A4/A4_overlapped_hector_gene.txt")

gene_intersect_A4_dmr <- fread("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/A4/A4_intersected_hector_gene.txt")

intersect(unique(gene_match_A4_dmr$V16), unique(gene_intersect_A4_dmr$V7))


gene_match_B4_dmr <- fread("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/B4/B4_overlapped_hector_gene.txt")

gene_intersect_B4_dmr <- fread("/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/B4/B4_intersected_hector_gene.txt")

intersect(unique(gene_match_B4_dmr$V16), unique(gene_intersect_B4_dmr$V7))
# 
# ## one can also try import from rtracklayer
# gene.bed.file = "/mnt/home3/reid/av638/methylseq/novaseq_em_seq/NVS090/data/outfolder/methylkit/gene_gencode_mouse_grcm39_out_chr.txt"
# gene.bed.file.gr <- toGRanges(gene.bed.file, format="BED", header=FALSE)
# oldmr_adj_bed_gr <- findOverlapsOfPeaks(dmr_adj_bed_gr, gene.bed.file.gr)
# #Annotating differentially methylated bases or regions
# library(genomation)
# 
# # read the gene BED file
# gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt",
#                                             package = "methylKit"))
# 
# 
# 
# exons = gffToGRanges(gtf.file, filter = "exon")
# CDS = gffToGRanges(gtf.file, filter = "CDS")
# transcript = gffToGRanges(gtf.file, filter = "transcript")
# genes = gffToGRanges(gtf.file, filter = "gene")
# 
# featuretargets = list(exons = exons, CDS = CDS, transcript=transcript, genes=genes)
# sm.featuretargets = ScoreMatrixList(targets = featuretargets, windows = genes, bin.num = NULL)

# annotate differentially methylated CpGs with
# promoter/exon/intron using annotation data
# #
# annotateWithGeneParts(as(myDiff25p.CompareTileB1,"GRanges"), featuretargets)
# 
# 
# # Read the CpG island annotation and annotate our differentially methylated bases/regions with them.
# # read the shores and flanking regions and name the flanks as shores
# # and CpG islands as CpGi
# cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt",
#                                      package = "methylKit"),
#                          feature.flank.name=c("CpGi","shores"))
# #
# # convert methylDiff object to GRanges and annotate
# diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
#                                     cpg.obj$CpGi,cpg.obj$shores,
#                                     feature.name="CpGi",flank.name="shores")
# 
# # read the CpG island annotation and annotate our differentially methylated bases/regions with them.
# # read the shores and flanking regions and name the flanks as shores
# # and CpG islands as CpGi
# cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt",
#                                      package = "methylKit"),
#                          feature.flank.name=c("CpGi","shores"))
# #
# # convert methylDiff object to GRanges and annotate
# diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
#                                     cpg.obj$CpGi,cpg.obj$shores,
#                                     feature.name="CpGi",flank.name="shores")
# 
# #Regional analysis
# #summarize methylation information over a set of defined regions such as promoters or CpG islands.
# promoters=regionCounts(myobj,gene.obj$promoters)
# 
# head(promoters[[1]])
# 
# #get the distance to TSS and nearest gene name
# diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
# 
# # target.row is the row number in myDiff25p
# head(getAssociationWithTSS(diffAnn))
# 
# #It is also desirable to get percentage/number of differentially methylated regions that overlap with intron/exon/promoters
# getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
# 
# #plot the percentage of differentially methylated bases overlapping with exon/intron/promoters
# plotTargetAnnotation(diffAnn,precedence=TRUE,
#                      main="differential methylation annotation")
# #percentage of differentially methylated bases are on CpG islands, CpG island shores and other regions.
# plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
#                      main="differential methylation annotation")
# #percentage of intron/exon/promoters that overlap with differentially methylated bases.
# getFeatsWithTargetsStats(diffAnn,percentage=TRUE)
# 
# #In case want to read back or share tabix file so created with methread, follow these steps:
# #Loading methylDB objects from tabix files
# #baseDB.obj <- makeMethylDB(methylBase.obj,"my/path")
# #mydbpath <- getDBPath(baseDB.obj)
# #rm(baseDB.obj)
# 
# #methylKit:::checkTabixHeader(mydbpath)
# #readMethylDB(mydbpath)
# 
# 
# 
# # #----------------------------  END OF ANALYSIS -------------------------#
# # #Full tested scripts are present in /mnt/home3/reid/av638/tutorial/methylkit/spare_scripts.R
