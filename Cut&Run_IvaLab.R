setwd("/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2")
path = getwd()
peak_counts <- read.table("/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/peak_counts.txt")
peak_counts <- peak_counts[order(peak_counts$V2),]
write.table(peak_counts, "peak_counts.txt", sep ="\t", append = F, quote = F, row.names = F, col.names = F)

# Barplot
ggplot(peak_counts, aes(x=V2, y=V1)) + 
  geom_bar(stat = "identity") +
  coord_flip()+theme_bw()
#Barplot_peak_counts.svg
CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene <- read.delim("CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene.txt", header = T)
CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene <- data.frame(CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene)
dim(CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene)
CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene_re <- CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene[,c(2:5,1,6:19)]
write.table(CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene_re, "CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene_re.txt", sep ="\t", append = F, quote = F, row.names = F, col.names = F)
grep chr CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene_re.txt > CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene_re_chr.txt
bedtools intersect -a CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene_re_chr.txt -b ~/motif_analysis/motif_GATA_hg38.bed > CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene_GATA.txt
bedtools intersect -a CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene_re_chr.txt -b ~/motif_analysis/motif_GATA_hg38.bed -v  > CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_gene_noGATA.txt


CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif <- read.delim("CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif.txt", header = T)
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif <- data.frame(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif)
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif)
head(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif)
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif[,c(2:4,1,6,5,7:24)] #column 4 should be ID
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re, "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re.txt", sep ="\t", append = F, quote = F, row.names = F, col.names = F)

CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_annotated_gene <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re[which(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re$Gene.Name != ""),]
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_Notannotated_gene <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re[which(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re$Gene.Name == ""),]

CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_annotated_gene_symbol <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_annotated_gene$Gene.Name
#From wikipedia (https://en.wikipedia.org/wiki/Histone) take genes and save as histone_genes_alias.txt
#upload 
library(xlsx)
histone_genes <- read.xlsx("/mnt/beegfs6/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/SynGO_id_convert_2022-10-16 12;32/idmap.xlsx", 1)
write.table(data.frame(gene=histone_genes$symbol), "all_histone_genes.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)

### Get promoter regions


#Annotation to regions
# library(annotatr)
# annots_basic_gelements = c('hg38_cpgs','hg38_basicgenes', 'hg38_genes_intergenic','hg38_genes_cds')
# 
# annotations_b_hg38_el = build_annotations(genome = 'hg38', annotations = annots_basic_gelements)
# 
# 
# annotations_b_hg38_eldf <- data.frame(annotations_b_hg38_el)
# annotations_b_hg38_eldf <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$start > 0),]
# annotations_b_hg38_eldf <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$end > 0),]
# 
# head(annotations_b_hg38_eldf)
# dim(annotations_b_hg38_eldf)
# write.table(data.frame(annotations_b_hg38_eldf), "/mnt/home3/reid/av638/chipseq/iva_lab_gencode/bams/annotations_b_hg38_eldf.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)

#Get Histone and Non-histone genes promoters
annotations_b_hg38_eldf_promoter <- read.table("/mnt/home3/reid/av638/chipseq/iva_lab_gencode/bams/annotations_b_hg38_eldf_promoter.txt")
colnames(annotations_b_hg38_eldf_promoter) <- c("chr", "start", "end", "size", "strand", "promoterid", "ensTid", "id2", "gene", "type")
annotations_b_hg38_eldf_promoter_histone <- merge(histone_genes, annotations_b_hg38_eldf_promoter, by.x="symbol",by.y = "gene")
dim(annotations_b_hg38_eldf_promoter_histone)
head(annotations_b_hg38_eldf_promoter_histone)
write.table(annotations_b_hg38_eldf_promoter_histone[,c(8:12,1:7,13:16)], "annotations_b_hg38_eldf_promoter_histone.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)

prom_Non_histone_genes <- setdiff(annotations_b_hg38_eldf_promoter$gene, histone_genes$symbol)
prom_Non_histone_genes <- gsub("H4C16","NA",prom_Non_histone_genes)
prom_Non_histone_genes <- data.frame(gene = prom_Non_histone_genes)

annotations_b_hg38_eldf_promoter_Non_histone <- merge(prom_Non_histone_genes, annotations_b_hg38_eldf_promoter, by = "gene")
dim(annotations_b_hg38_eldf_promoter_Non_histone)
head(annotations_b_hg38_eldf_promoter_Non_histone)
write.table(annotations_b_hg38_eldf_promoter_Non_histone[,c(2:6,1,7:10)], "annotations_b_hg38_eldf_promoter_Non_histone.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)


#Get CRAMP1 bound histone genes
CRAMP1_bound_histone_genes <- data.frame(symbol =intersect(histone_genes$symbol, CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_annotated_gene_symbol))
CRAMP1_bound_other_genes <- data.frame(symbol =setdiff(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_annotated_gene_symbol, histone_genes$symbol))
write.table(CRAMP1_bound_histone_genes, "CRAMP1_bound_histone_genes.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
write.table(CRAMP1_bound_other_genes, "CRAMP1_bound_other_genes.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifs.pl CRAMP1_bound_histone_genes.txt promoter_hg38.txt motifResults_CRAMP1_bound_histone_genes_promoter/ -len 10

#Gene bed ##
gene_v41_human_out <- read.table("/mnt/home3/reid/av638/ref_av/gencode/gene_v41_human_out.txt")
colnames(gene_v41_human_out) <- c("chr", "start", "end", "strand", "type", "ensID", "gene")
gene_v41_CRAMP1_bound_histone <- merge(CRAMP1_bound_histone_genes, gene_v41_human_out, by.x="symbol",by.y = "gene")
write.table(gene_v41_CRAMP1_bound_histone[,c(2:6,1,7)], "gene_v41_CRAMP1_bound_histone.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)


#Gene promoterrs
annotations_b_hg38_eldf_promoter_CRAMP1_bound_histone <- merge(CRAMP1_bound_histone_genes, annotations_b_hg38_eldf_promoter, by.x="symbol",by.y = "gene")
dim(annotations_b_hg38_eldf_promoter_CRAMP1_bound_histone)
head(annotations_b_hg38_eldf_promoter_CRAMP1_bound_histone)
write.table(annotations_b_hg38_eldf_promoter_CRAMP1_bound_histone[,c(2:6,1,7:10)], "annotations_b_hg38_eldf_promoter_CRAMP1_bound_histone.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)

#### Promoter ####
bedtools getfasta -fi ./../../../../GRCh38.p13.genome.fa -bed annotations_b_hg38_eldf_promoter_CRAMP1_bound_histone.txt -fo annotations_b_hg38_eldf_promoter_CRAMP1_bound_histone.fa
bedtools getfasta -fi ./../../../../GRCh38.p13.genome.fa -bed /mnt/home3/reid/av638/ref_av/gencode/gene_v41_human_out.txt -fo gene_v41_human_out.fa
findMotifs.pl annotations_b_hg38_eldf_promoter_CRAMP1_bound_histone.fa fasta motifResults_promoter_CRAMP1_bound_histone/ -fasta gene_v41_human_out.fa

#Get peaks near histone and non-histone genes

CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_histone <- merge(histone_genes, CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif, by.x="symbol", by.y = "Gene.Name")
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_histone[,c(9:13,1:8)], "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_histone.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)

# bedtools getfasta -fi ./../../../../GRCh38.p13.genome.fa -bed CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_histone.txt -fo CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_histone.fa
# findMotifs.pl CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_histone.fa fasta motifResults_peak_CRAMP1_bound_histone/ 

# CRAMP1 bound
# Note: use findMotifsGenome.pl find de novo and known motifs in regions in the genome becuase findMotifs.pl find de novo and known motifs in a gene list
## All genes
#---CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re.txt
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re
## Histone genes
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_histone <- merge(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re, histone_genes, by.x = "Gene.Name",by.y="symbol")
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_histone)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_histone[,c(2:6,1,7:30)], "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_histone.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_histone.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_histone
## Non-histone genes

Non_histone_genes <- data.frame(symbol = setdiff(gene_v41_human_out$gene, histone_genes$symbol))
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_Non_histone <- merge(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re, Non_histone_genes, by.x = "Gene.Name",by.y="symbol")
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_Non_histone)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_Non_histone[,c(2:6,1,7:24)], "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_Non_histone.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_Non_histone.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_Non_histone
#Check other genomic features
unique(unlist(lapply(str_split(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re$Annotation, " "), function(x) x[1])))
## Exonic
library(data.table)
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_exon <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re[CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re$Annotation %like% "exon", ]
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_exon)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_exon, "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_exon.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_exon.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_exon

## Intronic
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intron <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re[CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re$Annotation %like% "intron", ]
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intron)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intron, "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intron.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intron.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intron

## Intergenic
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intergenic <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re[CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re$Annotation %like% "Intergenic", ]
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intergenic)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intergenic, "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intergenic.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intergenic.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_intergenic

## Promoter-TSS
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_promoter_TSS <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re[CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re$Annotation %like% "promoter-TSS", ]
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_promoter_TSS)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_promoter_TSS, "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_promoter_TSS.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_promoter_TSS.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_promoter_TSS


#Gene expression
deseq2_res_C1KO_R7508_0.05 <- read.table("/mnt/beegfs6/home3/reid/av638/rnaseq/iva_lab_gencode/outfolder/deseq2_res_C1KO_R7508_0.05.txt")
deseq2_res_C1KO_R7508_0.05_symbol <- data.frame(symbol = unlist(lapply(str_split(rownames(deseq2_res_C1KO_R7508_0.05), "%"), function(x) x[2])))
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO <- merge(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re, deseq2_res_C1KO_R7508_0.05_symbol, by.x = "Gene.Name",by.y="symbol")
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO_histone <- merge(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO, histone_genes, by.x = "Gene.Name",by.y="symbol")
library(gprofiler2)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO <- gost(query=unique(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO$Gene.Name), organism="hsapiens")
gostplot(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted <- gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO$result[order(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO$result$p_value),]
head(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted)
#GO:BP
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp <- gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted[which(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted$source == "GO:BP"),]
head(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp)
dim(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar <- gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp[,c(11,3)]
head(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar_top <- head(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar,10)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar_top$term_name)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar_top$p_value <- -log10(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar_top$p_value)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar_top <- data.frame(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#5d93c4ff", 
          fill = "#5d93c4ff" ,
          sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          sort.val = "asc",
          rotate = TRUE, 
          position = position_dodge(),
          ggtheme = theme_bw())

ggsave("gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re_de_CRAMP1_KO.res.sorted_gobp_bar_top.svg", width=12, height=8, units="cm", dpi=96)


#Clash peak of CRAMP1 and histone marks from chip-seq and cutnrun data
# /mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/peak_histone_cramp_re.sh
# /mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/peak_histone_V5cramp_re.sh

#Now do everything with the peak quality > 100
#score > 100
#awk '{if($5 > 100) print $0}' CRAMP1_NTC_IgG_NTC_peaks.narrowPeak | grep chr > CRAMP1_NTC_IgG_NTC_peaks.narrowPeak_s100_chr.txt
#awk '{if($5 > 100) print $0}' V5-CRAMP1_V5_control_peaks.narrowPeak | grep chr > V5-CRAMP1_V5_control_peaks.narrowPeak_s100_chr.txt
#/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/peak_histone_cramp_re_s100.sh
#/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/peak_histone_V5cramp_re_s100.sh

# CRAMP1 bound
# Note: use findMotifsGenome.pl find de novo and known motifs in regions in the genome becuase findMotifs.pl find de novo and known motifs in a gene list
## All genes
awk '{if($5 > 100) print $0}' CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re.txt  > CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re.txt
#---CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re.txt
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re
## Histone genes
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re <- read.delim("CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re.txt", header = F, stringsAsFactors = F)
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re <- data.frame(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re)
colnames(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re) <- colnames(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_re)
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_histone <- merge(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re, histone_genes, by.x = "Gene.Name",by.y="symbol")
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_histone)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_histone[,c(2:6,1,7:30)], "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_histone.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_histone.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_histone
## Non-histone genes
gene_v41_human_out <- read.table("/mnt/home3/reid/av638/ref_av/gencode/gene_v41_human_out.txt")
colnames(gene_v41_human_out) <- c("chr", "start", "end", "strand", "type", "ensID", "gene")
Non_histone_genes <- data.frame(symbol = setdiff(gene_v41_human_out$gene, histone_genes$symbol))
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_Non_histone <- merge(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re, Non_histone_genes, by.x = "Gene.Name",by.y="symbol")
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_Non_histone)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_Non_histone[,c(2:6,1,7:24)], "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_Non_histone.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_Non_histone.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_Non_histone
#Check other genomic features
unique(unlist(lapply(str_split(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re$Annotation, " "), function(x) x[1])))
## Exonic
library(data.table)
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_exon <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re[CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re$Annotation %like% "exon", ]
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_exon)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_exon, "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_exon.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_exon.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_exon

## Intronic
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intron <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re[CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re$Annotation %like% "intron", ]
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intron)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intron, "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intron.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intron.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intron

## Intergenic
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intergenic <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re[CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re$Annotation %like% "Intergenic", ]
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intergenic)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intergenic, "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intergenic.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intergenic.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_intergenic

## Promoter-TSS
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_promoter_TSS <- CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re[CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re$Annotation %like% "promoter-TSS", ]
dim(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_promoter_TSS)
write.table(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_promoter_TSS, "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_promoter_TSS.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
findMotifsGenome.pl ./../CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_promoter_TSS.txt GRCh38 findMotifsGenome_CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_promoter_TSS


#Gene expression
deseq2_res_C1KO_R7508_0.05 <- read.table("/mnt/beegfs6/home3/reid/av638/rnaseq/iva_lab_gencode/outfolder/deseq2_res_C1KO_R7508_0.05.txt")
deseq2_res_C1KO_R7508_0.05_symbol <- data.frame(symbol = unlist(lapply(str_split(rownames(deseq2_res_C1KO_R7508_0.05), "%"), function(x) x[2])))
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO <- merge(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re, deseq2_res_C1KO_R7508_0.05_symbol, by.x = "Gene.Name",by.y="symbol")
CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO_histone <- merge(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO, histone_genes, by.x = "Gene.Name",by.y="symbol")
library(gprofiler2)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO <- gost(query=unique(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO$Gene.Name), organism="hsapiens")
gostplot(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO)
write.table(unique(CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO$Gene.Name), "CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.Gene.Name.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted <- gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO$result[order(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO$result$p_value),]
head(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted)
#GO:BP
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp <- gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted[which(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted$source == "GO:BP"),]
head(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp)
dim(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar <- gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp[,c(11,3)]
head(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar_top <- head(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar,10)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar_top$term_name)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar_top$p_value <- -log10(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar_top$p_value)
gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar_top <- data.frame(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#5d93c4ff", 
          fill = "#5d93c4ff" ,
          sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          sort.val = "asc",
          rotate = TRUE, 
          position = position_dodge(),
          ggtheme = theme_bw())

ggsave("gost.CRAMP1_NTC_IgG_NTC_peaks_gene_sel_motif_highQ_re_de_CRAMP1_KO.res.sorted_gobp_bar_top.svg", width=12, height=8, units="cm", dpi=96)


#Get histone marks overlap on cutnrun data

#ENCODE data needed some adjustment for peak iD, so added
awk '{print $1"\t"$2"\t"$3"\t""peaks_"NR"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' JUN_K562_GRCh38.bed > JUN_K562_GRCh38_peaks.bed
awk '{print $1"\t"$2"\t"$3"\t""peaks_"NR"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' SETDB1_K562_GRCh38.bed > SETDB1_K562_GRCh38_peaks.bed
awk '{print $1"\t"$2"\t"$3"\t""peaks_"NR"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' NFYA_K562_GRCh38.bed > NFYA_K562_GRCh38_peaks.bed
awk '{print $1"\t"$2"\t"$3"\t""peaks_"NR"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' NFYB_K562_GRCh38.bed > NFYB_K562_GRCh38_peaks.bed
awk '{print $1"\t"$2"\t"$3"\t""peaks_"NR"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' NFYC_HepG2_GRCh38.bed > NFYC_HepG2_GRCh38_peaks.bed
awk '{print $1"\t"$2"\t"$3"\t""peaks_"NR"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' SUZ12_MCF7_GRCh38.bed > SUZ12_MCF7_GRCh38_peaks.bed
awk '{print $1"\t"$2"\t"$3"\t""peaks_"NR"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' GATA2_K562_GRCh38.bed > GATA2_K562_GRCh38_peaks.bed
awk '{print $1"\t"$2"\t"$3"\t""peaks_"NR"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' GATA1_K562_GRCh38.bed > GATA1_K562_GRCh38_peaks.bed
#Run the below python script to get proportion of peak overlapped with histone to the total peaks present in the CUT&RUN sample
#python peak_counter.py > countout.txt
cutnrun_peaks_histone_peaks_counts <- read.table("/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/countswithhistone/countout.txt", header = F, stringsAsFactors = F)
head(cutnrun_peaks_histone_peaks_counts)
cutnrun_peaks_histone_peaks_counts <- cutnrun_peaks_histone_peaks_counts[,c(1,2,5)]
colnames(cutnrun_peaks_histone_peaks_counts) <- c("peak", "histone", "score")
cutnrun_peaks_histone_peaks_counts_tab <- tidyr::pivot_wider(cutnrun_peaks_histone_peaks_counts, names_from = histone, values_from = score)
cutnrun_peaks_histone_peaks_counts_df <- as.data.frame(cutnrun_peaks_histone_peaks_counts_tab)
rownames(cutnrun_peaks_histone_peaks_counts_df) <- cutnrun_peaks_histone_peaks_counts_df$peak
cutnrun_peaks_histone_peaks_counts_df <- cutnrun_peaks_histone_peaks_counts_df[,-1]
rownames(cutnrun_peaks_histone_peaks_counts_df) <- gsub("_peaks.narrowPeak", "", rownames(cutnrun_peaks_histone_peaks_counts_df))
colnames(cutnrun_peaks_histone_peaks_counts_df) <- gsub("_peaks.broadPeak.chr.txt", "", colnames(cutnrun_peaks_histone_peaks_counts_df))
head(cutnrun_peaks_histone_peaks_counts_df)
breaksListd <- seq(0, 1, by = 0.01)
pheatmap(cutnrun_peaks_histone_peaks_counts_df, na_col = "grey",color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListd)))



#python peak_counter.py > countout.txt
ENCODE_countout<- read.table("/mnt/home3/reid/av638/ENCODE/countout.txt", header = F, stringsAsFactors = F)
head(ENCODE_countout)
ENCODE_countout <- ENCODE_countout[,c(1,2,5)]
colnames(ENCODE_countout) <- c("peak", "histone", "score")
ENCODE_countout_tab <- tidyr::pivot_wider(ENCODE_countout, names_from = histone, values_from = score)
ENCODE_countout_df <- as.data.frame(ENCODE_countout_tab)
rownames(ENCODE_countout_df) <- ENCODE_countout_df$peak
ENCODE_countout_df <- ENCODE_countout_df[,-1]
head(ENCODE_countout_df)
rownames(ENCODE_countout_df) <- gsub("_peaks.narrowPeak", "", rownames(ENCODE_countout_df))
rownames(ENCODE_countout_df) <- gsub("_peaks_id.bed", "", rownames(ENCODE_countout_df))

colnames(ENCODE_countout_df) <- gsub("_peaks.broadPeak.chr.txt", "", colnames(ENCODE_countout_df))
head(ENCODE_countout_df)

ENCODE_countout_df_filt <- ENCODE_countout_df
ENCODE_countout_df_filt[is.na(ENCODE_countout_df_filt)] <- 0
breaksListd <- seq(0, 1, by = 0.01)
pheatmap(ENCODE_countout_df_filt, na_col = "grey",color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListd)))

EDASeq::plotPCA(as.matrix(t(ENCODE_countout_df_filt)))


#----------------------------------------------------------#
narrowpeakpath <- dir(path = "/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/analysis_nearest", pattern = ".narrowPeak_chr.bed")
library(ChIPpeakAnno)
library(openxlsx)
narrowpeakpathlist <- list()
for (peaks in narrowpeakpath){
  narrowpeakpathlist[[peaks]] <- paste0("/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/analysis_nearest/",peaks)
}
peak_df_list <- list()
peak_df_dis_list <- list()
peak_histonelist <- list()
for (peaksi in narrowpeakpathlist){
  if (file.info(peaksi)$size > 0){
    print(sapply(strsplit(peaksi,"/"), tail, 1))
    narrowpeak_peaks <- peaksi
    narrowpeak_peaksgr <- toGRanges(narrowpeak_peaks, format="BED", header=FALSE) 
    narrowpeak_peaksgr_anno <- annotatePeak(narrowpeak_peaksgr, tssRegion = c(-1000, 1000), annoDb = "org.Hs.eg.db")
    narrowpeak_peaksgr_anno_df <- as.data.frame(narrowpeak_peaksgr_anno)
    peak_df_list[[gsub("_peaks.narrowPeak_chr.bed","", sapply(strsplit(peaksi,"/"), tail, 1))]] <- narrowpeak_peaksgr_anno_df
    narrowpeak_peaksgr_anno_df_dis <- narrowpeak_peaksgr_anno_df[which(narrowpeak_peaksgr_anno_df$distanceToTSS > -100 & narrowpeak_peaksgr_anno_df$distanceToTSS < 100),]
    peak_df_dis_list[[gsub("_peaks.narrowPeak_chr.bed","", sapply(strsplit(peaksi,"/"), tail, 1))]] <- narrowpeak_peaksgr_anno_df_dis
    narrowpeak_peaksgr_anno_histone <- merge(narrowpeak_peaksgr_anno_df, histone_genes_all_matched, by.x="SYMBOL", by.y="gene_symbol")
    peak_histonelist[[gsub("_peaks.narrowPeak_chr.bed","", sapply(strsplit(peaksi,"/"), tail, 1))]] <- narrowpeak_peaksgr_anno_histone
  }
}
write.xlsx(peak_df_list, file = "/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/peak_df_list.xlsx")
write.xlsx(peak_df_dis_list, file = "/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/peak_df_dis_list.xlsx")
write.xlsx(peak_histonelist, file = "/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/peak_histonelist.xlsx")

### Functional analysis
#Full list
goterm_analysislist_peak <- list()
for (peaktemp in names(peak_df_list)){
  #Get genes
  # print(peaktemp)
  peaktempdf <- peak_df_list[[peaktemp]]
  # print(head(peaktempdf, 1))
  print(peaktemp)
  assign(paste0("gost.cutnrun.",peaktemp),gost(query=unique(peaktempdf$SYMBOL), organism="hsapiens"))
  goterm_analysislist_peak[[peaktemp]] <- get(paste0("gost.cutnrun.",peaktemp))$result
}
write.xlsx(goterm_analysislist_peak, file = "/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/goterm_analysislist_cutnrunpeak.xlsx")

#Full list
goterm_analysis_dis_list_peak <- list()
for (peakdistemp in names(peak_df_dis_list)){
  #Get genes
  # print(peakdistemp)
  if (length(peak_df_dis_list[[peakdistemp]]$seqnames) > 0){
    peakdistempdf <- peak_df_dis_list[[peakdistemp]]
    # print(head(peakdistempdf, 1))
    print(peakdistemp)
    assign(paste0("gost.cutnrun.dis.",peakdistemp),gost(query=unique(peakdistempdf$SYMBOL), organism="hsapiens"))
    goterm_analysis_dis_list_peak[[peakdistemp]] <- get(paste0("gost.cutnrun.dis.",peakdistemp))$result
  }
}
write.xlsx(goterm_analysis_dis_list_peak, file = "/mnt/home3/reid/av638/cutnrun/iva_lab_gencode/outfolder/bwa/mergedLibrary/macs2/goterm_analysislist_cutnrunpeak.xlsx")


# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# annoDataTxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# CRAMP1_NTC_IgG_NTC_peaksgr_anno <- annotatePeak(CRAMP1_NTC_IgG_NTC_peaksgr, tssRegion=c(-0, 0),TxDb=annoDataTxDb, annoDb="org.Hs.eg.db", overlap = "all")

#------------
library(ChIPseeker)


#Prepare annotation data
library(EnsDb.Hsapiens.v86) 
## create annotation file from EnsDb or TxDb
annoData <- toGRanges(EnsDb.Hsapiens.v86, feature="gene")
annoData[1:2]


genomicElementDistribution(CRAMP1_NTC_IgG_NTC_peaksgr, 
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))


CRAMP1_NTC_IgG_NTC_peaksgr_anno <- annotatePeakInBatch(CRAMP1_NTC_IgG_NTC_peaksgr, 
                                     AnnotationData=annoData,
                                     output="nearestBiDirectionalPromoters",
                                     bindingRegion=c(-2000, 500))
library(org.Hs.eg.db)
mart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")
CRAMP1_NTC_IgG_NTC_peaksgr_anno <- addGeneIDs(CRAMP1_NTC_IgG_NTC_peaksgr_anno, orgAnn="org.Hs.eg.db",
                        feature_id_type="gene_symbol", mart = mart,
                        IDs2Add=c("symbol"))
head(CRAMP1_NTC_IgG_NTC_peaksgr_anno)


#cut & run/homer
peaks_path = "/mnt/home3/reid/av638/cutnrun/new_iva_lab_gencode/outfolder/bwa/mergedLibrary/peak_ann/"
cutnrunpeaks_list <- list()
cutnrunpeaks_genes_list <- list()
cutnrunpeaks_gost_list <- list()
cutnrunpeaks_bar_list <- list()

cutnrunpeaks <- list.files(path=peaks_path,pattern = "*_gene.txt")
for (tempeaks in cutnrunpeaks){
  print(paste0("Sample in Analysis: ", tempeaks))
  tempeaksdf <- fread(paste0(peaks_path,tempeaks))
  tempeaksdf <- data.frame(tempeaksdf)
  cutnrunpeaks_list[[tempeaks]] <- tempeaksdf
  #Genes closer than 10Kb
  tempeaksdf <- tempeaksdf[which(tempeaksdf$Distance.to.TSS < 10000),]
  print(dim(tempeaksdf))
  tempeaksdf_genes <- unique(tempeaksdf$Nearest.Ensembl)
  tempeaksdf_genes <- tempeaksdf_genes[tempeaksdf_genes != '']
  print(length(tempeaksdf_genes))
  if (length(tempeaksdf_genes) >= 1){
    #gProfiler analysis
    cutnrunpeaks_genes_list[[tempeaks]] <- tempeaksdf_genes
    if(is.null(tempeaksdf_gost <- gost(query=tempeaksdf_genes, organism="hsapiens"))) next
    tempeaksdf_gost.sorted <- tempeaksdf_gost$result[order(tempeaksdf_gost$result$p_value),]
    tempeaksdf_gost.sorted["term_name_collapse"] <- paste0(tempeaksdf_gost.sorted$term_id,"_",tempeaksdf_gost.sorted$source,"_" ,gsub(" ", ".", tempeaksdf_gost.sorted$term_name))
    tempeaksdf_gost.sorted_gobp <- tempeaksdf_gost.sorted[which(tempeaksdf_gost.sorted$source == "GO:BP"),]
    tempeaksdf_gost.sorted_gobp_bar <- tempeaksdf_gost.sorted_gobp[,c(11,3)]
    tempeaksdf_gost.sorted_gobp_bar_top <- head(tempeaksdf_gost.sorted_gobp_bar,10)
    tempeaksdf_gost.sorted_gobp_bar_top$p_value <- -log10(tempeaksdf_gost.sorted_gobp_bar_top$p_value)
    cutnrunpeaks_gost_list[[tempeaks]] <- tempeaksdf_gost
    
    if (length(tempeaksdf_gost) >= 1){
      plotBar_tempeaks <- ggbarplot(tempeaksdf_gost.sorted_gobp_bar_top, x = "term_name", y = "p_value", color = "#5d93c4ff", fill = "#5d93c4ff" , sort.by.groups = FALSE,x.text.angle = 90, ylab = "-log10(p.value)", xlab = "Biological Process", legend.title = gsub("gost.","",tempeaks), lab.size = 9, sort.val = "asc", rotate = TRUE,  position = position_dodge(),ggtheme = theme_bw())
      cutnrunpeaks_bar_list[[tempeaks]] <- plotBar_tempeaks
    }
  }
}

pdf(file = "cutnrunpeaks_bar_list_new.pdf", height = 20, width = 20)
ggpubr::ggarrange(plotlist = cutnrunpeaks_bar_list, ncol=5, nrow=4, common.legend = F, labels=names(cutnrunpeaks_bar_list), vjust = 1,hjust=-0.5,font.label = list(size = 8, color = "black", face = "bold", family = NULL))
dev.off()


#cut & run/chipseeker
peaks_path = "/mnt/home3/reid/av638/cutnrun/new_iva_lab_gencode/outfolder/bwa/mergedLibrary/peak_ann/"
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
cutnrunpeaksfile <- list.files(path=peaks_path,pattern = "*_chr.txt")
cutnrunpeaksfile_list_gr <- list()
cutnrunpeaksfile_genes_list <- list()
cutnrunpeaksfile_gost_list <- list()
cutnrunpeaksfile_bar_list <- list()
for (peaksfile in cutnrunpeaksfile){
  print(peaksfile)
  peaks_temp <- read.table(paste0(peaks_path,peaksfile), header = F)
  peaks_tempgr <- makeGRangesFromDataFrame(peaks_temp,keep.extra.columns=TRUE, seqnames.field=c("V1"),start.field=c("V2"),end.field=c("V3"))
  cutnrunpeaksfile_list_gr[[peaksfile]] <- peaks_tempgr
  peaks_tempgr_anno <- annotatePeak(peaks_tempgr, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
  peaks_tempgr_anno_df <- as.data.frame(peaks_tempgr_anno)
  peaks_tempgr_anno_df <- peaks_tempgr_anno_df[which(peaks_tempgr_anno_df$distanceToTSS < 10000 | peaks_tempgr_anno_df$distanceToTSS > -10000 ),]
  peaks_tempgr_anno_df_genes <- unique(peaks_tempgr_anno_df$ENSEMBL)
  peaks_tempgr_anno_df_genes <- peaks_tempgr_anno_df_genes[!is.na(peaks_tempgr_anno_df_genes)]

  if (length(peaks_tempgr_anno_df_genes) >= 1){
    #gProfiler analysis
    cutnrunpeaksfile_genes_list[[peaksfile]] <- peaks_tempgr_anno_df_genes
    if(is.null(peaks_tempdf_gost <- gost(query=peaks_tempgr_anno_df_genes, organism="hsapiens"))) next
    peaks_tempdf_gost.sorted <- peaks_tempdf_gost$result[order(peaks_tempdf_gost$result$p_value),]
    peaks_tempdf_gost.sorted["term_name_collapse"] <- paste0(peaks_tempdf_gost.sorted$term_id,"_",peaks_tempdf_gost.sorted$source,"_" ,gsub(" ", ".", peaks_tempdf_gost.sorted$term_name))
    peaks_tempdf_gost.sorted_gobp <- peaks_tempdf_gost.sorted[which(peaks_tempdf_gost.sorted$source == "GO:BP"),]
    peaks_tempdf_gost.sorted_gobp_bar <- peaks_tempdf_gost.sorted_gobp[,c(11,3)]
    peaks_tempdf_gost.sorted_gobp_bar_top <- head(peaks_tempdf_gost.sorted_gobp_bar,10)
    peaks_tempdf_gost.sorted_gobp_bar_top$p_value <- -log10(peaks_tempdf_gost.sorted_gobp_bar_top$p_value)
    cutnrunpeaksfile_gost_list[[peaksfile]] <- peaks_tempdf_gost
    
    if (length(peaks_tempdf_gost) >= 1){
      plotBar_peaks_temp <- ggbarplot(peaks_tempdf_gost.sorted_gobp_bar_top, x = "term_name", y = "p_value", color = "#5d93c4ff", fill = "#5d93c4ff" , sort.by.groups = FALSE,x.text.angle = 90, ylab = "-log10(p.value)", xlab = "Biological Process", legend.title = gsub("gost.","",peaks_temp), lab.size = 9, sort.val = "asc", rotate = TRUE,  position = position_dodge(),ggtheme = theme_bw())
      cutnrunpeaksfile_bar_list[[peaksfile]] <- plotBar_peaks_temp
    }
  }
}

pdf(file = "cutnrunpeaks_chipseeker_bar_list_new.pdf", height = 20, width = 20)
ggpubr::ggarrange(plotlist = cutnrunpeaksfile_bar_list, ncol=5, nrow=4, common.legend = F, labels=names(cutnrunpeaks_bar_list), vjust = 1,hjust=-0.5,font.label = list(size = 8, color = "black", face = "bold", family = NULL))
dev.off()


### Annotating states to genomic features
##plotAnnoPie(temp_gr.txdb)
cutnrunpeaksAnno_list_gr <- lapply(cutnrunpeaksfile_list_gr, annotatePeak, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(cutnrunpeaksAnno_list_gr)
plotDistToTSS(cutnrunpeaksAnno_list_gr)

#get promoters
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
cutnruntagMatrixList <- lapply(cutnrunpeaksfile_list_gr, getTagMatrix, windows=promoter)
plotAvgProf(cutnruntagMatrixList, xlim=c(-1000, 1000))
plotAvgProf(cutnruntagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=500, facet="row")
tagHeatmap(cutnruntagMatrixList, xlim=c(-1000, 1000))
