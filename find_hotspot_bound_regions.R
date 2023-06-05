setwd("/Users/ankitverma/Documents/ENCODE/K562")
merged_encode_peak_count_500bp <- read.table("MERGED_final_cat_peaking.sorted_hg38_500bp_justbound_selected.bed", header = F, stringsAsFactors = F)
Count_merged_encode_peak_count_500bp <- count(merged_encode_peak_count_500bp, "V2")
library(plyr)
dim(Count_merged_encode_peak_count_500bp)

merged_encode_peak_count_500bp_clustered <- ddply(merged_encode_peak_count_500bp, .(V2), summarize, Proteins = toString(V1))
dim(merged_encode_peak_count_500bp)
head(merged_encode_peak_count_500bp)
dim(merged_encode_peak_count_500bp_clustered)
merged_encode_peak_count_500bp_clustered_count <- merge(merged_encode_peak_count_500bp_clustered, Count_merged_encode_peak_count_500bp, by="V2")
merged_encode_peak_count_500bp_clustered_count_sort <- merged_encode_peak_count_500bp_clustered_count[order(-merged_encode_peak_count_500bp_clustered_count$freq),]
library(splitstackshape)
merged_encode_peak_count_500bp_clustered_count_sort_pos <- cSplit(merged_encode_peak_count_500bp_clustered_count_sort, "V2", "%")
merged_encode_peak_count_500bp_clustered_count_sort_pos <- merged_encode_peak_count_500bp_clustered_count_sort_pos[,c(3,4,5,1,2)]
write.table(merged_encode_peak_count_500bp_clustered_count_sort_pos, "/Users/ankitverma/Documents/ENCODE/K562/merged_encode_peak_count_500bp_clustered_count_sort_pos.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
hist(merged_encode_peak_count_500bp_clustered_count_sort_pos$freq, breaks = 100)


merged_encode_peak_count_200bp <- read.table("MERGED_final_cat_peaking.sorted_hg38_200bp_justbound_selected.bed", header = F, stringsAsFactors = F)
Count_merged_encode_peak_count_200bp <- count(merged_encode_peak_count_200bp, "V2")
library(plyr)
dim(Count_merged_encode_peak_count_200bp)

merged_encode_peak_count_200bp_clustered <- ddply(merged_encode_peak_count_200bp, .(V2), summarize, Proteins = toString(V1))
dim(merged_encode_peak_count_200bp)
head(merged_encode_peak_count_200bp)
dim(merged_encode_peak_count_200bp_clustered)
merged_encode_peak_count_200bp_clustered_count <- merge(merged_encode_peak_count_200bp_clustered, Count_merged_encode_peak_count_200bp, by="V2")
merged_encode_peak_count_200bp_clustered_count_sort <- merged_encode_peak_count_200bp_clustered_count[order(-merged_encode_peak_count_200bp_clustered_count$freq),]
library(splitstackshape)
merged_encode_peak_count_200bp_clustered_count_sort_pos <- cSplit(merged_encode_peak_count_200bp_clustered_count_sort, "V2", "%")
merged_encode_peak_count_200bp_clustered_count_sort_pos <- merged_encode_peak_count_200bp_clustered_count_sort_pos[,c(3,4,5,1,2)]
write.table(merged_encode_peak_count_200bp_clustered_count_sort_pos, "/Users/ankitverma/Documents/ENCODE/K562/merged_encode_peak_count_200bp_clustered_count_sort_pos.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
hist(merged_encode_peak_count_200bp_clustered_count_sort_pos$freq, breaks = 100)



merged_encode_peak_count_1000bp <- read.table("MERGED_final_cat_peaking.sorted_hg38_1000bp_justbound_selected.bed", header = F, stringsAsFactors = F)
Count_merged_encode_peak_count_1000bp <- count(merged_encode_peak_count_1000bp, "V2")
library(plyr)
dim(Count_merged_encode_peak_count_1000bp)

merged_encode_peak_count_1000bp_clustered <- ddply(merged_encode_peak_count_1000bp, .(V2), summarize, Proteins = toString(V1))
dim(merged_encode_peak_count_1000bp)
head(merged_encode_peak_count_1000bp)
dim(merged_encode_peak_count_1000bp_clustered)
merged_encode_peak_count_1000bp_clustered_count <- merge(merged_encode_peak_count_1000bp_clustered, Count_merged_encode_peak_count_1000bp, by="V2")
merged_encode_peak_count_1000bp_clustered_count_sort <- merged_encode_peak_count_1000bp_clustered_count[order(-merged_encode_peak_count_1000bp_clustered_count$freq),]
library(splitstackshape)
merged_encode_peak_count_1000bp_clustered_count_sort_pos <- cSplit(merged_encode_peak_count_1000bp_clustered_count_sort, "V2", "%")
merged_encode_peak_count_1000bp_clustered_count_sort_pos <- merged_encode_peak_count_1000bp_clustered_count_sort_pos[,c(3,4,5,1,2)]
write.table(merged_encode_peak_count_1000bp_clustered_count_sort_pos, "/Users/ankitverma/Documents/ENCODE/K562/merged_encode_peak_count_1000bp_clustered_count_sort_pos.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
hist(merged_encode_peak_count_1000bp_clustered_count_sort_pos$freq, breaks = 100)
merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50 <- merged_encode_peak_count_1000bp_clustered_count_sort_pos[which(merged_encode_peak_count_1000bp_clustered_count_sort_pos$freq <= 50),]

write.table(merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50, "/Users/ankitverma/Documents/ENCODE/K562/merged_encode_peak_count_1000bp_clustered_count_sort_pos_cutoff50.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
save.image("~/Documents/ENCODE/K562/summits/hotspotregions/hsr_workspace.RData")
