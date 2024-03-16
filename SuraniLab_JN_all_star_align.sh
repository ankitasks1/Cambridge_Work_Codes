echo "Analysis Started"
#echo "Analysing C8"
#cd C8
#gunzip *
#~/tools_av/STAR/STAR-2.7.10b/bin/Linux_x86_64/STAR --genomeDir /mnt/home3/reid/av638/scrnaseq/surani_lab/star_gencode_ref_v43/ --outFileNamePrefix run_SITTC8_star_ --soloFeatures Gene Velocyto --soloType CB_UMI_Simple --soloCBwhitelist /mnt/home3/reid/av638/tools_av/yard/apps/cellranger-7.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt --soloUMIlen 12 --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --runThreadN 3 --readFilesIn SITTC8_S4_L001_R2_001.fastq SITTC8_S4_L001_R1_001.fastq
#cd ..

echo "Analysing D8"
cd D8
#gunzip *
~/tools_av/STAR/STAR-2.7.10b/bin/Linux_x86_64/STAR --genomeDir /mnt/home3/reid/av638/scrnaseq/surani_lab/star_gencode_ref_v43/ --outFileNamePrefix run_SITTD8_star_ --soloFeatures Gene Velocyto --soloType CB_UMI_Simple --soloCBwhitelist /mnt/home3/reid/av638/tools_av/yard/apps/cellranger-7.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt --soloUMIlen 12 --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --runThreadN 3 --readFilesIn SITTD8_S4_L001_R2_001.fastq SITTD8_S4_L001_R1_001.fastq
cd ..

echo "Analysing E1"
cd E1
#gunzip *
~/tools_av/STAR/STAR-2.7.10b/bin/Linux_x86_64/STAR --genomeDir /mnt/home3/reid/av638/scrnaseq/surani_lab/star_gencode_ref_v43/ --outFileNamePrefix run_SITTE1_star_ --soloFeatures Gene Velocyto --soloType CB_UMI_Simple --soloCBwhitelist /mnt/home3/reid/av638/tools_av/yard/apps/cellranger-7.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt --soloUMIlen 12 --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --runThreadN 3 --readFilesIn SITTE1_S4_L001_R2_001.fastq SITTE1_S4_L001_R1_001.fastq
cd ..

echo "Analysing F1"
cd F1
#gunzip *
~/tools_av/STAR/STAR-2.7.10b/bin/Linux_x86_64/STAR --genomeDir /mnt/home3/reid/av638/scrnaseq/surani_lab/star_gencode_ref_v43/ --outFileNamePrefix run_SITTF1_star_ --soloFeatures Gene Velocyto --soloType CB_UMI_Simple --soloCBwhitelist /mnt/home3/reid/av638/tools_av/yard/apps/cellranger-7.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt --soloUMIlen 12 --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --runThreadN 3 --readFilesIn SITTF1_S4_L001_R2_001.fastq SITTF1_S4_L001_R1_001.fastq
cd ..


echo "Analysing G1"
cd G1
gunzip *
~/tools_av/STAR/STAR-2.7.10b/bin/Linux_x86_64/STAR --genomeDir /mnt/home3/reid/av638/scrnaseq/surani_lab/star_gencode_ref_v43/ --outFileNamePrefix run_SITTG1_star_ --soloFeatures Gene Velocyto --soloType CB_UMI_Simple --soloCBwhitelist /mnt/home3/reid/av638/tools_av/yard/apps/cellranger-7.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt --soloUMIlen 12 --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --runThreadN 3 --readFilesIn SITTG1_S4_L001_R2_001.fastq SITTG1_S4_L001_R1_001.fastq
cd ..

echo "Analysing G9"
cd G9
gunzip *
~/tools_av/STAR/STAR-2.7.10b/bin/Linux_x86_64/STAR --genomeDir /mnt/home3/reid/av638/scrnaseq/surani_lab/star_gencode_ref_v43/ --outFileNamePrefix run_SITTG9_star_ --soloFeatures Gene Velocyto --soloType CB_UMI_Simple --soloCBwhitelist /mnt/home3/reid/av638/tools_av/yard/apps/cellranger-7.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt --soloUMIlen 12 --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --runThreadN 3 --readFilesIn SITTG9_S4_L001_R2_001.fastq SITTG9_S4_L001_R1_001.fastq
cd ..

echo "Analysing H10"
cd H10
gunzip *
~/tools_av/STAR/STAR-2.7.10b/bin/Linux_x86_64/STAR --genomeDir /mnt/home3/reid/av638/scrnaseq/surani_lab/star_gencode_ref_v43/ --outFileNamePrefix run_SITTH10_star_ --soloFeatures Gene Velocyto --soloType CB_UMI_Simple --soloCBwhitelist /mnt/home3/reid/av638/tools_av/yard/apps/cellranger-7.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt --soloUMIlen 12 --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --runThreadN 3 --readFilesIn SITTH10_S4_L001_R2_001.fastq SITTH10_S4_L001_R1_001.fastq
cd ..
echo "Analysis Finished"
