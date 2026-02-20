# Multiomics

## Step 1: Add Gene Symbols to the Gene Count Matrix

### Code File: Step1_GeneCountMatrix_AddGeneSymbol.R

Input: gexp_counts.csv

Output: gexp_counts_symbol.csv
 
## Step 2: Assess RNA-seq Sample Groups / Conduct Differential Expression Analysis

### Code File: Step2_RNAseq_Sample_Analysis.R

Input: gexp_counts_symbol.csv from step 1

Output Files
> [Sample Group]__DEG.csv
>
> [Sample Group]__UpDEG.csv
>
> [Sample Group]__DownDEG.csv
>
> gexp_counts_symbol_normalized.csv
>
> PCA_Plot.jpg
>
> Euclidean_Distance_Heatmap.jpg
>
> Pearson_Correlation_Heatmap.jpg

## Step 3: Categorize Differentially Expressed Genes

### Code File: Step3_DEG_Categorization.R

Input Files
> [Sample Group]__UpDEG.csv from step 2
> 
> [Sample Group]__DownDEG.csv from step 2

Output Files
> UpDEG_category.csv
>
> DownDEG_category.csv

## Step 4: Pre-process Public HiC Dataset

### Code File: Step4_Public_HiC_Preprocessing.R

Input Files
> [Sample].boundaries files from public HiC dataset in hg19
> 
> hg19ToHg38.over.chain from UCSC genome browser database

Output File: [Sample]_hg38_TADs.csv

## Step 5: Pre-process ChIP-seq

### Code File: Step5_ChIPseq_Preprocessing.R

Input Files

> [Sample ChIP-seq].broadPeak files from ChIP-seq experiment
>
> [Public ChIP-seq Samples].bed files from ENCODE project

Output Files

> [Sample ChIP-seq] Filtered.csv
>
> [Public ChIP-seq Samples] Filtered.csv

## Step 6: Multiomics

### Code File: Step6_Multiomics.py

Input Files

> UpDEG_category.csv from step 3
>
> [Sample]_hg38_TADs.csv file from Step 4
>
> [Sample ChIP-seq] Filtered.csv from step 5

Output Files

> [Sample]_HiC_RNAseq.csv
>
> [Sample]_HiC_RNAseq_ChIPseq_UpDEG_Positive.csv
>
> [Sample]_HiC_RNAseq_ChIPseq.csv
>
> [Sample]_HiC_RNAseq_ChIPseq_UpDEG_Positive.csv
>
> [Sample]_HiC_RNAseq_ChIPseq_Condensed.csv
>
> [Sample]_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks.csv

## Step 7: Cytoband Addition

### Code File: Step7_Multiomics_Cytoband.R

Input File: [Sample]_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks.csv from step 6

Output Files

> [Sample]_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks_Cytobands.csv
>
> [Sample]_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks_Cytobands_SRA_Peak_Specific.csv

## Step 8: ChIP-seq Peak Overlap Assessment

### Code File: Step8_Multiomics_Overlaps.R

Input Files
>
> [Sample].broadPeak files from ChIP-seq
>
> UpDEG_category.csv from step 3
>
> DownDEG_category.csv from step 3
>
> [Sample]_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks_Cytobands.csv from step 7

Output Files
>
> ChIP_seq_[Sample]_[Sample]_overlap.csv
>
> ChIP_seq_[Sample]_[Sample]_overlap_percentage.csv
>
> ChIPseq_RNAseq_Overlap.csv
>
> ChIPseq_RNAseq_Overlap_Filtered.csv
>
> ChIPseq_RNAseq_Overlap_Filtered_Sorted.csv
>
> ChIPseq_RNAseq_Overlap_Filtered_Sorted_Condensed.csv
>
> [Sample]_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks_Cytobands_Venn Diagram.jpg

## Step 9: TAD Bar Graphs

### Code File: Step9_Multiomics_Bar_Graphs.R

Input Files
> 
> UpDEG_category.csv from step 3
>
> DownDEG_category.csv from step 3
>
> [Sample]_HiC_Boundaries.csv from public dataset

Output Files
>
> UpDEG_category.bed
>
> DownDEG_category.bed
>
> [Sample]_HiC_Boundaries.bed
>
> [Sample]_TAD_Chr[#].jpg

## Step 10: Gene Ontology

### Code File: Step10_Gene_Ontology.R

Input Files
> 
> UpDEG_category.csv from step 3
>
> DownDEG_category.csv from step 3
>
> [Sample] Specific UpDEG Panther GO Analysis.xml from https://geneontology.org/
>
> [Sample] Specific DownDEG Panther GO Analysis.xml from https://geneontology.org/

Output Files
>
> [Sample]_Specific_UpDEG_category.csv
>
> [Sample]_Specific_DownDEG_category.csv
>
> [Sample] Specific UpDEG Panther GO Analysis.csv
>
> [Sample] Specific DownDEG Panther GO Analysis.csv
>
> [Sample] Specific UpDEG Panther GO Analysis Merged.csv
>
> [Sample] Specific DownDEG Panther GO Analysis Merged.csv
>
> [Sample] Specific UpDEG Panther GO Analysis Merged ChIPseq.csv

## Step 11: Peak Counts per 10 Mb Region

### Code File: Step11_ChIPseq_Heatmaps.R

Input Files
> 
> [Sample].broadPeak files from ChIP-seq dataset
>
> [Sample].bed files from ChIP-seq public datasets

Output Files
>
> [Sample]_10Mb.csv files

## Step 12: Visualize Genomic Annotation

### Code File: Step12_Visualize_Genomic_Annotation.R

Input Files
> 
> [Sample].broadPeak files from ChIP-seq dataset

Output Files
>
> [Sample]_Visualize_Genomic_Annotation_hg38.jpg files

## Step 13: Create ChIP-seq TSS Plots

### Code File: Step13_ChIPseq_TSS_Plot.R

Input Files
> 
> [Sample].broadPeak files from ChIP-seq dataset

Output Files
>
> [Sample]_TSS_Plot.jpg files

## Step 14: Create RNA-seq plots

### Code File: Step14_RNAseq_Plots.R

Input Files
> 
> gexp_counts_symbol_normalized.csv from step 2
>
> [Sample Group]__DEG.csv from step 2
>
> UpDEG_category.csv from step 3

Output Files
>
> [Sample]_Normalized_GeneCount_BoxPlot.jpg
>
> [Sample]_log2FC.jpg
>
> [Sample]_VolcanoPlot.jpg

## Step 15: Enhancer analysis

### Code File: Step15_Enhancer.R

Input Files
>
> UpDEG_category.csv from step 3
>
> [Sample ChIP-seq] Filtered.csv from step 5
>
> GeneHancer_AnnotSV_elements_v5.25.txt from GeneCards
>
> GeneHancer_AnnotSV_gene_association_scores_v5.25.txt from GeneCards
>
> GeneHancer_TFBSs_v5.25.txt from GeneCards

Output Files
>
> [Sample]_Enhancer.csv
>
> [Sample]_Enhancer_Gene.csv
>
> [Sample]_Enhancer_Gene_UpDEG.csv
>
> UpDEG_category_merged.csv
>
> [Sample]_Enhancer_Gene_UpDEG_TF.csv
>
> [Sample]_Enhancer_Gene_UpDEG_TF_Coor.csv
>
> [Sample]_Enhancer_Gene_UpDEG_TF_Histone.csv
>
> ChIPseq_Enhancer_Gene_UpDEG_TF Frequency by Category.csv

## Step 16: Analyze Overlap between .bed files

### Code File: Step16_ChIPseq_Overlap.R

Input Files
> 
> .bed files from multiple sources
>
> [Sample ChIP-seq] Filtered.csv from step 5
>
> [Public ChIP-seq Samples] Filtered.csv from step 5

Output Files
>
> trimmed .bed files
>
> [Sample ChIP-seq] Filtered_trim.bed
>
> [Public ChIP-seq Samples] Filtered_trim.bed
>
> [.bed #1]_overlap_[.bed #2].csv
>
> [.bed #1]_overlap_[.bed #2]_summary.csv
>
> [.bed #1]_overlap_[.bed #2]_Venn_Diagram.jpg
>
> [.bed #1]_overlap_[.bed #2]_bar_Chart.jpg
>
> [.bed #1]_overlap_[.bed #2]_Correlation.jpg

## Step 17: UpDEG preidction from Histone Marks

### Code File

> ### Step16_Linux_ChIPseq_Overlap.R
>
> ### Step16_Windows_ChIPseq_Overlap.R

Input Files
> 
> [Sample]_Enhancer_Gene_UpDEG_TF_Histone.csv from step 15
> > [ChIP Peak]_overlap_[UpDEG].bed hand-generated from the csv file
>
> [Sample ChIP-seq] Filtered.csv from step 5
>
> [Public ChIP-seq Samples] Filtered.csv from step 5

Output Files
>
> [ChIP Peak]_overlap_[UpDEG]_trim.bed
>
> [Sample ChIP-seq] Filtered_trim.csv
>
> [Public ChIP-seq Samples] Filtered_trim.csv
>
> [.bed #1]_overlap_[.bed #2].csv
>
> [.bed #1]_overlap_[.bed #2]_summary.csv
>
> [.bed #1]_overlap_[.bed #2]_bar_Chart.jpg

## Step 18: Protein binding preidction from Histone Marks

> ### Step18_1_ML_Prediction_hg38_bin.py
>
> ### Step18_2_ML_Prediction_pos_neg_sets.py
>
> ### Step18_3_ML_Prediction_histone_bigwig.py
>
> ### Step18_4_ML_Prediction_Model_Training.py
>
> ### Step18_5_ML_Prediction_Test_Result_bed.py

Input Files
> 
> [Sample ChIP-seq] Filtered_trim.csv from step 17
>
> hg38.chrom.sizes from UCSC
>
> ENCSR636HFF GRCh38_unified_blacklist.bed from ENCODE
>
> Homo_sapiens.GRCh38.dna.primary_assembly.fa from ENCODE
>
> [Public ChIP-seq Samples].bigWing from ENCODE

Output Files
>
> hg38_bin.csv
>
> hg38_bin_positive.csv
>
> hg38_bin_negative.csv
>
> hg38_bin_positive_histone.csv
>
> hg38_bin_negative_histone.csv
>
> feature_importance_scores.csv
>
> test_confusion_matrix.csv
>
> test_performance_metrics.csv
>
> test_predictions.csv
>
> validation_f1_score.csv
>
> training_validation_metrics.csv
>
> hg38_bin_feature_matrix_histone_test.csv
>
> hg38_bin_feature_matrix_histone_validation.csv
>
> hg38_bin_feature_matrix_histone_training.csv
>
> confusion_matrix_heatmap.jpg
>
> precision_recall_curve.jpg
>
> roc_curve.jpg
 
## MISC: Conversion from .bedGraph to .bigWig

### Code File: Misc_bedGraph_to_bigWig.R

Input File: [Sample].bedGraph

Output File: [Sample].bigWig
