#sequence : http://ftp.ensembl.org/pub/release-107/fasta/drosophila_melanogaster/pep/
#transmembrane.txt downloaded from ensembl
#adjusted for column name space
setwd("/Users/ankitverma/Documents/peptide_data")
transmembrane <- read.delim('/Users/ankitverma/Documents/peptide_data/transmembrane.txt', header = T, stringsAsFactors = F)
transmembrane <- transmembrane[,-c(10)]
write.table(transmembrane, "/Users/ankitverma/Documents/peptide_data/transmembrane_re.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)

dim(transmembrane)
transmembrane["row"] <- paste0(transmembrane$Gene_stable_ID,"%",
                              transmembrane$Transcript_stable_ID,"%",
                              transmembrane$Transmembrane_helices,"%",
                              transmembrane$Transmembrane_helices_start,"%",
                              transmembrane$Transmembrane_helices_end,"%",
                              transmembrane$Protein_stable_ID,"%",
                              transmembrane$Gene_name,"%",
                              transmembrane$FlyBase_gene_ID,"%",
                              transmembrane$FlyBase_gene_name_ID)

library(dplyr)
transmembrane_uniq <- transmembrane %>% distinct()
dim(transmembrane_uniq)
detach(package:dplyr)

library(plyr)
transmembrane_me <- ddply(transmembrane_uniq, .(Protein_stable_ID), summarize, 
                          Gene_name = toString(unique(Gene_name)), 
                          FlyBase_gene_name_ID = toString(unique(FlyBase_gene_name_ID)),
                          Transmembrane_helices = toString(unique(Transmembrane_helices)),
                          Transmembrane_helices_start = toString(unique(Transmembrane_helices_start)),
                          Transmembrane_helices_end = toString(unique(Transmembrane_helices_end)),
                          Gene_stable_ID = toString(unique(Gene_stable_ID)))



dim(transmembrane_me)
head(transmembrane_me)



count_transmembrane <- plyr::count(transmembrane_uniq, "Protein_stable_ID")

#Assign count of TM helix to proteins
transmembrane_me_count <- merge(transmembrane_me, count_transmembrane, by="Protein_stable_ID")
dim(transmembrane_me_count)
#remove blank TM helix rows
transmembrane_me_count <- transmembrane_me_count[which(transmembrane_me_count$Transmembrane_helices != ""),]
dim(transmembrane_me_count)
head(transmembrane_me_count)
#Remove undesired spaces
transmembrane_me_count$Transmembrane_helices <- gsub(" ", "", transmembrane_me_count$Transmembrane_helices)
transmembrane_me_count$Transmembrane_helices_start <- gsub(" ", "", transmembrane_me_count$Transmembrane_helices_start)
transmembrane_me_count$Transmembrane_helices_end <- gsub(" ", "", transmembrane_me_count$Transmembrane_helices_end)
transmembrane_me_count$Protein_stable_ID <- gsub(" ", "", transmembrane_me_count$Protein_stable_ID)


transmembrane_me_count["Transmembrane_helices_start"]  <- do.call(rbind.data.frame,lapply(str_split(transmembrane_me_count$Transmembrane_helices_start, pattern =","), function(x) paste0(sort(as.numeric(x)),collapse = ",")))
transmembrane_me_count["Transmembrane_helices_end"] <- do.call(rbind.data.frame,lapply(str_split(transmembrane_me_count$Transmembrane_helices_end, pattern =","), function(x) paste0(sort(as.numeric(x)),collapse = ",")))
head(transmembrane_me_count)
transmembrane_me_count <- transmembrane_me_count[order(-transmembrane_me_count$freq),]

write.table(transmembrane_me_count, "/Users/ankitverma/Documents/peptide_data/transmembrane_me_count.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)


#Get the fasta and overlap with transmembrane file
#python /Users/ankitverma/Documents/peptide_data/signal_peptide_fasta.py
signal_over_fasta_output <- read.table("/Users/ankitverma/Documents/peptide_data/signal_over_fasta_output.txt", header = F, stringsAsFactors = F)
colnames(signal_over_fasta_output) <- c("signal", "num", "sig_start", "sig_end", "Protein_stable_ID", "gene_id", "symbol", "seq")
signal_over_fasta_output["seq_len"] <- nchar(signal_over_fasta_output$seq)

#Get only longest isoform
Re_signal_over_fasta_output <- signal_over_fasta_output[,c(6,5,9)]
library(dplyr)
Re_signal_over_fasta_output <- Re_signal_over_fasta_output %>% distinct()
detach(package:dplyr)
Re_signal_over_fasta_output["P_seq_len"] <- paste0(Re_signal_over_fasta_output$seq_len, "_", Re_signal_over_fasta_output$Protein_stable_ID)

library(plyr)
Re_signal_over_fasta_output_collapse <- ddply(Re_signal_over_fasta_output, .(gene_id), summarize, 
                                      P_seq_len = toString((P_seq_len)))

Re_signal_over_fasta_output_collapse$P_seq_len <- gsub(" ", "", Re_signal_over_fasta_output_collapse$P_seq_len)
Re_signal_over_fasta_output_collapse["P_seq_len_sort"]  <- do.call(rbind.data.frame,lapply(str_split(Re_signal_over_fasta_output_collapse$P_seq_len, pattern =","), function(x) paste0(str_sort(x, numeric = T,decreasing = T),collapse = ",")))
library(splitstackshape)
Re_signal_over_fasta_output_collapse_col <- cSplit(Re_signal_over_fasta_output_collapse, "P_seq_len_sort", ",")
Re_signal_over_fasta_output_collapse_col <- data.frame(Re_signal_over_fasta_output_collapse_col)

single_isoform_ids <- data.frame(pID =unique(unlist(lapply(str_split(Re_signal_over_fasta_output_collapse_col$P_seq_len_sort_01, pattern = "_"), function(x) x[2]))))

signal_over_fasta_transmembrane <- merge(signal_over_fasta_output, transmembrane_me_count, by.x="Protein_stable_ID", by.y="Protein_stable_ID")

#Assign to gene expression data (all EC type)
aEC1 <- read.table("/Users/ankitverma/Documents/peptide_data/aEC1.txt", header = F, stringsAsFactors = F)
colnames(aEC1) <- c("Gene_stable_ID", "aEC1gene_exp", "aEC1readcount", "aEC1id")
aEC2 <- read.table("/Users/ankitverma/Documents/peptide_data/aEC2.txt", header = F, stringsAsFactors = F)
colnames(aEC2) <- c("Gene_stable_ID", "aEC2gene_exp", "aEC2readcount", "aEC2id")
aEC3 <- read.table("/Users/ankitverma/Documents/peptide_data/aEC3.txt", header = F, stringsAsFactors = F)
colnames(aEC3) <- c("Gene_stable_ID", "aEC3gene_exp", "aEC3readcount", "aEC3id")
aEC4 <- read.table("/Users/ankitverma/Documents/peptide_data/aEC4.txt", header = F, stringsAsFactors = F)
colnames(aEC4) <- c("Gene_stable_ID", "aEC4gene_exp", "aEC4readcount", "aEC4id")
pEC1 <- read.table("/Users/ankitverma/Documents/peptide_data/pEC1.txt", header = F, stringsAsFactors = F)
colnames(pEC1) <- c("Gene_stable_ID", "pEC1gene_exp", "pEC1readcount", "pEC1id")
pEC2 <- read.table("/Users/ankitverma/Documents/peptide_data/pEC2.txt", header = F, stringsAsFactors = F)
colnames(pEC2) <- c("Gene_stable_ID", "pEC2gene_exp", "pEC2readcount", "pEC2id")
pEC3 <- read.table("/Users/ankitverma/Documents/peptide_data/pEC3.txt", header = F, stringsAsFactors = F)
colnames(pEC3) <- c("Gene_stable_ID", "pEC3gene_exp", "pEC3readcount", "pEC3id")
dEC <- read.table("/Users/ankitverma/Documents/peptide_data/dEC.txt", header = F, stringsAsFactors = F)
colnames(dEC) <- c("Gene_stable_ID", "dECgene_exp", "dECreadcount", "dECid")
mEC <- read.table("/Users/ankitverma/Documents/peptide_data/mEC.txt", header = F, stringsAsFactors = F)
colnames(mEC) <- c("Gene_stable_ID", "mECgene_exp", "mECreadcount", "mECid")

#Make a master list of all genes expressed in atleast one cell
sc_expression_master  <- data.frame(Gene_stable_ID= unique(c(aEC1$Gene_stable_ID, aEC2$Gene_stable_ID, aEC3$Gene_stable_ID, aEC4$Gene_stable_ID, pEC1$Gene_stable_ID, pEC2$Gene_stable_ID, pEC3$Gene_stable_ID, dEC$Gene_stable_ID, mEC$Gene_stable_ID)))

#Get the matrix of expression pattern in each cell type
sc_expression_master_aEC1 <- merge(sc_expression_master, aEC1, by = "Gene_stable_ID", all.x = T)
sc_expression_master_aEC1_aEC2 <- merge(sc_expression_master_aEC1, aEC2, by = "Gene_stable_ID", all.x = T)
sc_expression_master_aEC1_aEC2_aEC3 <- merge(sc_expression_master_aEC1_aEC2, aEC3, by = "Gene_stable_ID", all.x = T)
sc_expression_master_aEC1_aEC2_aEC3_aEC4 <- merge(sc_expression_master_aEC1_aEC2_aEC3, aEC4, by = "Gene_stable_ID", all.x = T)
sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1 <- merge(sc_expression_master_aEC1_aEC2_aEC3_aEC4, pEC1, by = "Gene_stable_ID", all.x = T)
sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1_pEC2 <- merge(sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1, pEC2, by = "Gene_stable_ID", all.x = T)
sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1_pEC2_pEC3 <- merge(sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1_pEC2, pEC3, by = "Gene_stable_ID", all.x = T)
sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1_pEC2_pEC3_dEC <- merge(sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1_pEC2_pEC3, dEC, by = "Gene_stable_ID", all.x = T)
sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1_pEC2_pEC3_dEC_mEC <- merge(sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1_pEC2_pEC3_dEC, mEC, by = "Gene_stable_ID", all.x = T)
sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1_pEC2_pEC3_dEC_mEC <- sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1_pEC2_pEC3_dEC_mEC[,c(1,2,5,8,11,14,17,20,23,26)]
head(sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1_pEC2_pEC3_dEC_mEC)

sc_expression <- sc_expression_master_aEC1_aEC2_aEC3_aEC4_pEC1_pEC2_pEC3_dEC_mEC
write.table(sc_expression, "/Users/ankitverma/Documents/peptide_data/sc_expression_data.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)

#Get the expression of genes with signal hits
signal_over_fasta_transmembrane_exp <- merge(signal_over_fasta_transmembrane, sc_expression, by.x="Gene_stable_ID", by.y="Gene_stable_ID")
signal_over_fasta_transmembrane_exp <- signal_over_fasta_transmembrane_exp[order(-signal_over_fasta_transmembrane_exp$freq),]
head(signal_over_fasta_transmembrane_exp)
write.table(signal_over_fasta_transmembrane_exp, "/Users/ankitverma/Documents/peptide_data/signal_over_fasta_transmembrane_exp.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)


signal_over_fasta_transmembrane_exp_re <- signal_over_fasta_transmembrane_exp[,c(4,1,8,2,9,10,3,5,6,13:25)]
dim(signal_over_fasta_transmembrane_exp_re)
signal_over_fasta_transmembrane_exp_re["Sum_expression"] <-  rowSums(signal_over_fasta_transmembrane_exp_re[,c(14:22)], na.rm = T)
signal_over_fasta_transmembrane_exp_re["Num_expression"] <- rowSums(!is.na(signal_over_fasta_transmembrane_exp_re[,c(14:22)]))
write.table(signal_over_fasta_transmembrane_exp_re, "/Users/ankitverma/Documents/peptide_data/signal_over_fasta_transmembrane_exp_re.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)

#Extract only the longest isoform from dataframe
signal_over_fasta_transmembrane_exp_re_single_iso <- merge(signal_over_fasta_transmembrane_exp_re, single_isoform_ids, by.x = "Protein_stable_ID", by.y = "pID")

library(plyr)
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed <- ddply(signal_over_fasta_transmembrane_exp_re_single_iso, .(Gene_stable_ID), summarize,
                                                                     symbol = toString(unique(symbol)),
                                                                     Protein_stable_ID = toString(unique(Protein_stable_ID)),
                                                                     seq = toString(unique(seq)),
                                                                     seq_len = toString(unique(seq_len)),
                                                                     signal = toString(unique(signal)),
                                                                     sig_start = toString(unique(sig_start)),
                                                                     sig_end = toString(unique(sig_end)),
                                                                     Transmembrane_helices = toString(unique(Transmembrane_helices)),
                                                                     Transmembrane_helices_start = toString(unique(Transmembrane_helices_start)),
                                                                     Transmembrane_helices_end = toString(unique(Transmembrane_helices_end)),
                                                                     freq = toString(unique(freq)),
                                                                     aEC1gene_exp = toString(unique(aEC1gene_exp)), 
                                                                     aEC2gene_exp = toString(unique(aEC2gene_exp)),
                                                                     aEC3gene_exp = toString(unique(aEC3gene_exp)),
                                                                     aEC4gene_exp = toString(unique(aEC4gene_exp)),
                                                                     pEC1gene_exp = toString(unique(pEC1gene_exp)),
                                                                     pEC2gene_exp = toString(unique(pEC2gene_exp)),
                                                                     pEC3gene_exp = toString(unique(pEC3gene_exp)),
                                                                     dECgene_exp = toString(unique(dECgene_exp)),
                                                                     mECgene_exp = toString(unique(mECgene_exp)),
                                                                     Sum_expression = unique(Sum_expression),
                                                                     Num_expression = unique(Num_expression))

signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Sum_expression <- as.numeric(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Sum_expression)
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Num_expression <- as.numeric(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Num_expression)


signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Protein_stable_ID <- gsub(" ", "", signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Protein_stable_ID)
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$signal <- gsub(" ", "", signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$signal)
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$sig_start <- gsub(" ", "", signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$sig_start)
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$sig_end <- gsub(" ", "", signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$sig_end)
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$freq <- gsub(" ", "", signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$freq)
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Transmembrane_helices_start <- gsub(" ", "", signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Transmembrane_helices_start)
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Transmembrane_helices_end <- gsub(" ", "", signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Transmembrane_helices_end)

signal_over_fasta_transmembrane_exp_re_single_iso_collapsed["Transmembrane_helices_start"]  <- do.call(rbind.data.frame,lapply(str_split(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Transmembrane_helices_start, pattern =","), function(x) paste0(sort(as.numeric(x)),collapse = ",")))
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed["Transmembrane_helices_end"] <- do.call(rbind.data.frame,lapply(str_split(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Transmembrane_helices_end, pattern =","), function(x) paste0(sort(as.numeric(x)),collapse = ",")))


signal_over_fasta_transmembrane_exp_re_single_iso_collapsed <-  data.frame(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed)
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed <- signal_over_fasta_transmembrane_exp_re_single_iso_collapsed[order(-signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Num_expression),]


signal_over_fasta_transmembrane_exp_re_single_iso_collapsed["signalcol"] <- sapply(1:nrow(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed), function(x) paste0(unlist(strsplit(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$sig_start[x], split = ",",  fixed = TRUE)),
                                                                                                                                                                           "_" ,
                                                                                                                                                                           unlist(strsplit(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$signal[x], split = ",",  fixed = TRUE)),
                                                                                                                                                                           "_" ,
                                                                                                                                                                           unlist(strsplit(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$sig_end[x], split = ",",  fixed = TRUE)),
                                                                                                                                                                           collapse = ","))

signal_over_fasta_transmembrane_exp_re_single_iso_collapsed["signalcol"]  <- do.call(rbind.data.frame,lapply(str_split(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$signalcol, pattern =","), function(x) paste0(str_sort(x, numeric = T,decreasing = F),collapse = ",")))
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed["Max_TM_end"] <- do.call(rbind.data.frame,lapply(str_split(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Transmembrane_helices_end, pattern =","), function(x) max(as.numeric(x))))


signal_over_fasta_transmembrane_exp_re_single_iso_collapsed_Max_TM_end_list <-  lapply(str_split(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$Max_TM_end, pattern =","), function(x) as.numeric(x))
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed_sig_start_list <-  lapply(str_split(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$sig_start, pattern=","), function(x) as.numeric(x))
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed["sub_TM_sig_start"] <- do.call(rbind.data.frame, lapply(mapply('-', signal_over_fasta_transmembrane_exp_re_single_iso_collapsed_sig_start_list, signal_over_fasta_transmembrane_exp_re_single_iso_collapsed_Max_TM_end_list, SIMPLIFY = FALSE), function(x) paste0(x, collapse = ",")))
signal_over_fasta_transmembrane_exp_re_single_iso_collapsed["Min_Distance_TM_end_to_sig_start"] <- do.call(rbind.data.frame,lapply(str_split(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed$sub_TM_sig_start, pattern =","), function(x) max(as.numeric(x))))
write.table(signal_over_fasta_transmembrane_exp_re_single_iso_collapsed, "/Users/ankitverma/Documents/peptide_data/signal_over_fasta_transmembrane_exp_re_single_iso_collapsed.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)

#---------------Python---------#

python_signal_over_fasta_tm_exp_output <- read.table("/Users/ankitverma/Documents/peptide_data/withpython/signal_over_fasta_tm_exp_output.txt", header = T, stringsAsFactors = F)
length(unique(signal_over_fasta_transmembrane_exp$Gene_stable_ID))
length(unique(python_signal_over_fasta_tm_exp_output$Gene_Stable_ID))
length(intersect(unique(signal_over_fasta_transmembrane_exp$Gene_stable_ID), unique(python_signal_over_fasta_tm_exp_output$Gene_Stable_ID)))


length(unique(signal_over_fasta_transmembrane_exp$Protein_stable_ID))
length(unique(python_signal_over_fasta_tm_exp_output$Protein_Stable_ID))
length(intersect(unique(signal_over_fasta_transmembrane_exp$Protein_stable_ID), unique(python_signal_over_fasta_tm_exp_output$Protein_Stable_ID)))

#Python and R results matched
