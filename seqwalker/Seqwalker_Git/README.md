# <u> Seqwalker </u>


### Installation
<pre>
library(devtools)
install_github("ankitasks1/Seqwalker", force=TRUE)
library(Seqwalker)
</pre>
  
### Dependencies: Install and load following packages latest version

<pre>
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
</pre>

</n>
<b> <i> This R package offer a comprehensive suite of functions for analyzing various types of omics data, particularly focusing on differential expression/binding analysis and genomic feature annotation. Here's a brief description of what each function does and the overall purpose of the package: </i> </b>

1. **Data Loading and Preprocessing**: Functions like `load_rdata` and `read_genefile` facilitate loading raw data files and gene annotation files, respectively, into R for further analysis.

2. **Data Filtering and Manipulation**: Functions such as `sample_filter` and `reassign_genenames` help filter samples based on certain criteria and reassign gene names, respectively, to ensure data quality and compatibility with downstream analyses.

3. **Differential Expression Analysis**: The package offers functions like `deseq2_de` and `limma_de` for performing differential expression analysis, which helps identify genes that are differentially expressed between different experimental conditions.

4. **Visualization**: Functions like `plot_venn` and various plotting functions (e.g., `ggmaplot`, `scatterplot`) provide visualization tools for exploring and presenting the results of differential expression analysis and other omics data analyses.

5. **Integration and Annotation**: Functions such as `integrate_featurecounts`, `integrate_bedtools`, and `annotate_peaks_to_genes` facilitate the integration of multiple omics datasets and provide annotation for genomic features, which can enhance the interpretability and utility of the analyzed data.

**Full details about each function and parameters are below:**

## Description about each function

### load_rdata
| Parameter | Description |
|-----------|-------------|
| path      | Location of the file |
| rdata     | RData object |

### read_genefile
| Parameter | Description |
|-----------|-------------|
| path      | Location of the file |
| file      | Gene file in tab-separated format |

### dds_counts
| Parameter | Description |
|-----------|-------------|
| dds       | DESeq2 dds object |

### sample_filter
| Parameter      | Description |
|----------------|-------------|
| countmatrix    | Gene expression raw counts matrix |
| samplefilter   | Flag indicating whether to filter samples |
| columnsremove  | Index of samples to remove |

### reassign_genenames
| Parameter      | Description |
|----------------|-------------|
| countmatrix    | Count matrix |
| genedf         | Dataframe of gene coordinates and names |
| columnskeep    | Columns to keep from the count matrix |

### make_coldata_dds
| Parameter      | Description |
|----------------|-------------|
| dds            | DESeq2 dds object |
| samplefilter   | If any samples to filter otherwise skip |
| rowsremove     | Index of samples to filter |

### deseq2_de
| Parameter | Description |
|-----------|-------------|
| counts    | Count matrix (genes are rows and samples are columns with unnormalized values) |
| coldata   | Metadata |
| fdr       | p.adj-value threshold |
| fc        | Fold change threshold |

### plot_venn
| Parameter | Description |
|-----------|-------------|
| list      | List of samples (n <= 4) with vector of gene symbols |
| color     | Vector of colors to assign |

Certainly! Below are tables summarizing each function:

### load_rdata
| Parameter | Description |
|-----------|-------------|
| path      | Location of the file |
| rdata     | RData object |

### read_genefile
| Parameter | Description |
|-----------|-------------|
| path      | Location of the file |
| file      | Gene file in tab-separated format |

### dds_counts
| Parameter | Description |
|-----------|-------------|
| dds       | DESeq2 dds object |

### sample_filter
| Parameter      | Description |
|----------------|-------------|
| countmatrix    | Gene expression raw counts matrix |
| samplefilter   | Flag indicating whether to filter samples |
| columnsremove  | Index of samples to remove |

### reassign_genenames
| Parameter      | Description |
|----------------|-------------|
| countmatrix    | Count matrix |
| genedf         | Dataframe of gene coordinates and names |
| columnskeep    | Columns to keep from the count matrix |

### make_coldata_dds
| Parameter      | Description |
|----------------|-------------|
| dds            | DESeq2 dds object |
| samplefilter   | If any samples to filter otherwise skip |
| rowsremove     | Index of samples to filter |

### plot_venn
| Parameter | Description |
|-----------|-------------|
| list      | List of samples (n <= 4) with vector of gene symbols |
| color     | Vector of colors to assign |

### edger_create
| Parameter | Description |
|-----------|-------------|
| counts    | Count matrix (genes are rows and samples are columns with unnormalized values) |
| coldata   | Metadata |

### edger_de
| Parameter | Description |
|-----------|-------------|
| lmfit     | Linear models fitted object prepared while performing edger_create() |
| design    | Design matrices for specified variables |
| fdr       | p.adj-value threshold |
| fc        | Fold change threshold |

### diffbind_prep
| Parameter    | Description |
|--------------|-------------|
| samplesheet  | Samplesheet formatted as per Diffbind |

### diffbind_count
| Parameter | Description |
|-----------|-------------|
| obj       | Diffbind object |
| summit    | Expected summit for peaks |

### diffbind_norm
| Parameter | Description |
|-----------|-------------|
| obj       | Diffbind object |

### diffbind_de
| Parameter  | Description |
|------------|-------------|
| normobj    | Diffbind normalized object |
| fdr        | p.adj-value threshold |
| fc         | Fold change threshold |
| method     | "edgeR" or "deseq2" |

### limma_create
| Parameter | Description |
|-----------|-------------|
| data      | Count matrix |
| coldata   | Metadata |

### limma_de
| Parameter  | Description |
|------------|-------------|
| lmfit      | Linear fitted model object |
| contrasts  | Contrast for samples |
| fdr        | p.adj-value threshold |
| fc         | Fold change threshold |

### diffbind_gene_anno
| Parameter  | Description |
|------------|-------------|
| assaytype  | Assay type (e.g., "atacseq", "rnaseq") |
| software   | Internal software to be used by Diffbind (e.g., "edger" or "deseq2") |
| contrasts  | Contrasts list |
| genefile   | Gene file name |
| path       | Location of files |
| outformat  | Output format (e.g., ".txt") |
| fdr        | p.adj-value threshold |
| fc         | Fold change threshold |


### diffbind_and_packages_gene_anno
| Parameter  | Description |
|------------|-------------|
| assaytype  | Assay type (e.g., "atacseq", "rnaseq") |
| software   | Internal software to be used by Diffbind (e.g., "edger" or "deseq2") |
| contrasts  | Contrasts list |
| genefile   | Gene file name |
| path       | Location of files |
| outformat  | Output format (e.g., ".txt") |
| fdr        | p.adj-value threshold |
| fc         | Fold change threshold |

### consensus_gene_anno
| Parameter  | Description |
|------------|-------------|
| assaytype  | Assay type (e.g., "atacseq", "rnaseq") |
| software   | Internal software to be used by Diffbind (e.g., "edger" or "deseq2") |
| contrasts  | Contrasts list |
| genefile   | Gene file name |
| path       | Location of files |
| outformat  | Output format (e.g., ".txt") |
| fdr        | p.adj-value threshold |
| fc         | Fold change threshold |

### consensus_gene_count_anno
| Parameter  | Description |
|------------|-------------|
| assaytype  | Assay type (e.g., "atacseq", "rnaseq") |
| software   | Internal software to be used by Diffbind (e.g., "edger" or "deseq2") |
| contrasts  | Contrasts list |
| genefile   | Gene file name |
| path       | Location of files |
| outformat  | Output format (e.g., ".txt") |

### run_chipseeker
| Parameter       | Description |
|-----------------|-------------|
| peaks_path      | Location of peaks files (e.g., .bed, .narrowPeak, .broadPeak) |
| assaytype       | Assay type (e.g., "atacseq") |
| txdb            | TxDb object |
| org_db          | org.db object |
| peakformat      | Peak file format (e.g., ".bed") |

### chipseeker_region_anno
| Parameter       | Description |
|-----------------|-------------|
| path            | Location of peaks files (e.g., .bed, .narrowPeak, .broadPeak) |
| assaytype       | Assay type (e.g., "atacseq") |
| software        | Internal software to be used by Diffbind (e.g., "deseq2") |
| contrasts       | Contrasts list |
| txdb            | TxDb object |
| org_db          | org.db object |
| regionformat    | Region file format (e.g., ".bed") |
| fdr             | p.adj-value threshold |
| fc              | Fold change threshold |

### replot_pie_chipseeker: chipseeker pie plot
| Parameter | Description |
|-----------|-------------|
| gr_anno   | Granges annotation object |

### quantify_featurecounts: featurecounts quantify per feature in a given interval
| Parameter             | Description |
|-----------------------|-------------|
| peaks_path            | Location of peaks files |
| assaytype             | Assay type (e.g., "atacseq") |
| bamfiles              | BAM file extension |
| sites_files           | Site files extension |
| sites_type            | Sites type |
| sites_column_to_rearrange | Columns to rearrange |
| pairedend             | Paired-end flag |
| refgenome             | Reference genome |
| delim                 | Delimiter |

### quantify_featurecounts_multifeature
| Parameter            | Description |
|----------------------|-------------|
| assay_bam_path      | Location of BAM files |
| feature_path        | Location of feature files |
| assay               | Assay type (e.g., "cutnandrun") |
| bam_extension       | BAM file extension |
| feature_list        | List of features |
| columnstorearrange  | Columns to rearrange |
| ref                 | Reference genome |
| control             | Control |
| pe                  | Paired-end flag |
| diff                | Differential mode |
| pairs               | List of sample IDs and controls |

### quantify_featurecounts_multifeature_div
| Parameter           | Description |
|---------------------|-------------|
| assay_bam_path      | Location of BAM files |
| feature_path        | Location of feature files |
| assay               | Assay type (e.g., "cutnandrun") |
| bam_extension       | BAM file extension |
| feature_list        | List of features |
| columnstorearrange  | Columns to rearrange |
| ref                 | Reference genome |
| control             | Control |
| pe                  | Paired-end flag |
| diff                | Differential mode |
| pairs               | List of sample IDs and controls |

### quantify_bams: Rsamtools quantify reads in a given bam files
| Parameter          | Description |
|--------------------|-------------|
| peaks_path         | Location of peak files |
| assaytype          | Assay type (e.g., "cutandrun") |
| bamfiles           | BAM file extension |

### summedreads_per_feature: Further process the featurecounts quantified matrix
| Parameter           | Description |
|---------------------|-------------|
| fearturecountmatrix | Count matrix |
| bamcountsreads      | DataFrame of total reads in BAM files |
| labels              | Labels of features short name |

### aggregate_feature: aggregate value per feature
| Parameter            | Description |
|----------------------|-------------|
| annomatrix           | DataFrame |
| columnstoaggregate   | Columns to aggregate |
| aggregatebycolumn    | Column for group by |
| operation            | Operation (e.g., "mean") |

### make_pca: pca prcomp using counts
| Parameter         | Description |
|-------------------|-------------|
| counts            | Count matrix |
| file_extension    | BAM file extension |
| attribute         | PCA attribute |

### integrate_featurecounts: integrate multiple omics data on featurecounts
| Parameter     | Description |
|---------------|-------------|
| inputlist     | List of dataframes from different samples/results |
| sharedcolumn  | Shared column |

Of course, continuing with the tables for the remaining functions:

### integrate_bedtools: integrate multiple omics data on bedtools coverage output
| Parameter           | Description |
|---------------------|-------------|
| omics_data_path_lists | List of multiple outputs of bedtools coverage |
| binsize             | Bin size (e.g., "5kb") |
| columns_pos         | Columns positions |
| columns_samples     | Positions of coverage column across all samples |
| refgenome           | Reference genome |
| delim               | Delimiter |

### upset_out: Upset plot
| Parameter  | Description |
|------------|-------------|
| input_list | List of samples with vector of features (e.g., genes) |

### maplot: ggmaplot maplot
| Parameter        | Description |
|------------------|-------------|
| df               | DataFrame of output from DESeq2 or other methods |
| columnbyorder    | Columns order |
| fdr              | p.adj-value threshold |
| fc               | Fold change threshold |

### maplot_general: ggplot maplot
| Parameter        | Description |
|------------------|-------------|
| df               | DataFrame of output from DESeq2 or other methods |
| columnbyorder    | Columns order |
| fdr              | p.adj-value threshold |
| fc               | Fold change threshold |

### maplot_general_shrinkage: deseq2 maplot using ggplot2 using shrinkage
| Parameter        | Description |
|------------------|-------------|
| df               | DataFrame of output from DESeq2 or other methods |
| columnbyorder    | Columns order |
| fdr              | p.adj-value threshold |
| fc               | Fold change threshold |

### rearrange_groups: base rearrange groups
| Parameter       | Description |
|-----------------|-------------|
| grouplist       | List of paths to specific genes (DataFrame), e.g., cluster1.txt, cluster2.txt |

### meta_plot: ggplot meta plots
| Parameter        | Description |
|------------------|-------------|
| featurematrix    | Count matrix |
| columnstoplot    | Columns with value to pick for plotting |
| multi            | Multi-flag |
| value            | Column name with values to pick |
| color            | Column name with color |
| linetype         | Column name with line type |
| rep              | Repetition |
| splitside        | Delimiter separator side |
| colorbyid        | Color by ID flag |

### annotate_peaks_to_genes: bedtools annotation
| Parameter       | Description |
|-----------------|-------------|
| peakfileslist   | List of peak files |
| genefile        | DataFrame of gene file |
| path            | Location of peaks files |
| assaytype       | Assay type (e.g., "atacseq") |
| software        | Internal software to be used by Diffbind (e.g., "deseq2") |
| outformat       | Output format (e.g., ".txt") |

### go_term_analysis: gprofiler go term analysis
| Parameter       | Description |
|-----------------|-------------|
| gene_list       | List of multiple samples with genes as vectors |
| organism        | Organism (e.g., "hsapiens") |

### scatterplot: ggplot scatter
| Parameter       | Description |
|-----------------|-------------|
| df              | DataFrame |

### integrate multiple omics data on bedtools coverage output
| Parameter           | Description |
|---------------------|-------------|
| omics_data_path_lists | List of multiple outputs of bedtools coverage |
| binsize             | Bin size (e.g., "5kb") |
| columns_pos         | Columns positions |
| columns_samples     | Positions of coverage column across all samples |
  
### deseq2 maplot using ggplot2 using shrinkage
| Parameter   | Description |
|-------------|-------------|
| df          | DataFrame of output from DESeq2 or other methods |
| columnbyorder | Columns to specify the order |
| fdr         | p.adj-value threshold |
| fc          | Fold change threshold |

### ggplot meta plots
| Parameter         | Description |
|-------------------|-------------|
| featurematrix     | Count matrix |
| columnstoplot     | Columns to plot |
| multi             | Multi-flag |
| value             | Column name with values |
| color             | Column name with colors |
| linetype          | Column name with linetypes |
| rep               | Repetition |
| splitside         | Delimiter separator side |
| colorbyid         | Color by ID flag |


Overall, Seqwalker provides a comprehensive set of tools for analyzing and interpreting omics data, particularly focusing on gene expression analysis and genomic feature annotation. By offering functions for data loading, preprocessing, differential expression analysis, visualization, and integration, the package aims to support researchers in gaining insights into biological processes and identifying potential biomarkers or targets for further investigation.
