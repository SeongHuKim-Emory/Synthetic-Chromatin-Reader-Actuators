# --- 0. Setup: Load Libraries and Define Paths ---

# Check if BiocManager is installed, install if not
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# List of required Bioconductor packages
bioc_packages <- c("DESeq2", "vsn", "AnnotationDbi", "org.Hs.eg.db") # Added annotation packages

# Install missing Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# List of required CRAN packages
cran_packages <- c("tidyverse", "pheatmap", "RColorBrewer", "ggplot2") # Added ggplot2 explicitly

# Install missing CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

# Load libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(vsn) # For VST transformation used in PCA/Heatmaps
library(ggplot2) # For PCA plot

# --- Define File Paths ---
# Adjust base_path if your script is not in the same directory as the 'public' folder
base_path <- "../public/RNAseq/"
input_file <- file.path(base_path, "gexp_counts_symbol.csv")
output_dir <- base_path

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Output file paths
# Data files
pca_plot_data_file <- file.path(output_dir, "PCA_Plot_Data.csv")
euclidean_heatmap_data_file <- file.path(output_dir, "Euclidean_Distance_Matrix.csv")
pearson_heatmap_data_file <- file.path(output_dir, "Pearson_Correlation_Matrix.csv")
normalized_counts_file <- file.path(output_dir, "gexp_counts_symbol_normalized.csv")

# Plot files (NEW)
pca_plot_file <- file.path(output_dir, "PCA_Plot.jpg")
euclidean_heatmap_plot_file <- file.path(output_dir, "Euclidean_Distance_Heatmap.jpg")
pearson_heatmap_plot_file <- file.path(output_dir, "Pearson_Correlation_Heatmap.jpg")


deg_output_files <- list(
  "MCF7_021_10_vs_0" = file.path(output_dir, "MCF7_021_10_DEG.csv"),
  "MCF7_021_24_vs_0" = file.path(output_dir, "MCF7_021_24_DEG.csv"),
  "MCF7_021_48_vs_0" = file.path(output_dir, "MCF7_021_48_DEG.csv"),
  "MCF7_023_24_vs_0" = file.path(output_dir, "MCF7_023_24_DEG.csv"),
  "MCF7_025_24_vs_0" = file.path(output_dir, "MCF7_025_24_DEG.csv")
)

updeg_output_files <- list(
  "MCF7_021_10_vs_0" = file.path(output_dir, "MCF7_021_10_UpDEG.csv"),
  "MCF7_021_24_vs_0" = file.path(output_dir, "MCF7_021_24_UpDEG.csv"),
  "MCF7_021_48_vs_0" = file.path(output_dir, "MCF7_021_48_UpDEG.csv"),
  "MCF7_023_24_vs_0" = file.path(output_dir, "MCF7_023_24_UpDEG.csv"),
  "MCF7_025_24_vs_0" = file.path(output_dir, "MCF7_025_24_UpDEG.csv")
)

downdeg_output_files <- list(
  "MCF7_021_10_vs_0" = file.path(output_dir, "MCF7_021_10_DownDEG.csv"),
  "MCF7_021_24_vs_0" = file.path(output_dir, "MCF7_021_24_DownDEG.csv"),
  "MCF7_021_48_vs_0" = file.path(output_dir, "MCF7_021_48_DownDEG.csv"),
  "MCF7_023_24_vs_0" = file.path(output_dir, "MCF7_023_24_DownDEG.csv"),
  "MCF7_025_24_vs_0" = file.path(output_dir, "MCF7_025_24_DownDEG.csv")
)

# --- 1. Import Data ---

message("1. Importing count data...")
count_data_full <- readr::read_csv(input_file, show_col_types = FALSE)

# Separate annotation (transcript_id, symbol) from counts
gene_annotation <- count_data_full %>% dplyr::select(transcript_id, symbol)
count_matrix <- count_data_full %>%
  dplyr::select(-transcript_id, -symbol) %>%
  as.matrix()

# Set transcript_id as rownames for the count matrix
rownames(count_matrix) <- count_data_full$transcript_id

# Ensure counts are integers
count_matrix <- round(count_matrix)
mode(count_matrix) <- "integer"

message("Count data dimensions: ", nrow(count_matrix), " genes, ", ncol(count_matrix), " samples.")

# --- Prepare Metadata (colData) ---
message("Preparing sample metadata (colData)...")
sample_names <- colnames(count_matrix)

# Parse sample names to create metadata columns
colData <- data.frame(row.names = sample_names)
colData$sample_id <- sample_names
colData$cell_line <- sapply(strsplit(sample_names, "_"), `[`, 1)
colData$treatment_id <- sapply(strsplit(sample_names, "_"), `[`, 2)
colData$timepoint <- as.factor(sapply(strsplit(sample_names, "_"), `[`, 3)) # Factor for DESeq2 design
colData$replicate <- sapply(strsplit(sample_names, "_"), `[`, 4)

# Create a combined 'condition' factor for DESeq2 design formula
# This represents the unique group for each sample (e.g., MCF7_021_0, MCF7_021_10)
colData$condition <- factor(paste(colData$cell_line, colData$treatment_id, colData$timepoint, sep = "_"))

# Check if rownames of colData match colnames of count_matrix
if (!all(rownames(colData) == colnames(count_matrix))) {
  stop("Mismatch between sample names in count matrix columns and metadata rows.")
}

print("Sample Metadata (colData):")
print(head(colData))

# --- 2. Create DESeqDataSet Object and Run DESeq2 ---

message("Creating DESeqDataSet object...")
# Design formula using the combined condition factor
design_formula <- ~ condition

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = design_formula)

# --- Optional: Pre-filtering lowly expressed genes ---
# Keep genes that have at least 10 reads total across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
message("Removed lowly expressed genes. Genes remaining: ", nrow(dds))

message("Running DESeq2 analysis (normalization, dispersion estimation, GLM fitting)...")
# This single command runs the core DESeq2 steps
dds <- DESeq(dds)
message("DESeq2 analysis complete.")

# --- 3. Normalization and Quality Control ---

# Get normalized counts (median-of-ratios method)
normalized_counts <- counts(dds, normalized=TRUE)

# --- Variance Stabilizing Transformation (VST) for QC ---
# VST is generally recommended for PCA and heatmaps with DESeq2 data
# Use blind=FALSE if many genes are expected to be DE, but blind=TRUE is standard for general QC
message("Performing Variance Stabilizing Transformation (VST)...")
vsd <- vst(dds, blind=TRUE) # Use blind=TRUE for unbiased QC

# --- PCA Analysis ---
message("Performing PCA analysis...")
# Use DESeq2's plotPCA function, retrieve data
pcaData <- plotPCA(vsd, intgroup=c("condition", "timepoint", "treatment_id"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Save PCA data (CSV)
message("Saving PCA data to: ", pca_plot_data_file)
pcaData_to_save <- pcaData %>%
  rownames_to_column("sample_id") %>% # Keep sample names
  dplyr::select(sample_id, PC1, PC2, condition, timepoint, treatment_id, group) # Ensure all relevant columns are present
readr::write_csv(pcaData_to_save, pca_plot_data_file)

# Generate and save the PCA plot itself (JPG) - MODIFIED
message("Generating and saving PCA plot...") # Added message
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=timepoint)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  ggtitle("PCA Plot of Samples") +
  theme_bw()
ggsave(filename = pca_plot_file, plot = pca_plot, width = 8, height = 6, device = "jpeg") # MODIFIED path and device
message("PCA plot saved to: ", pca_plot_file) # MODIFIED message


# --- Sample Distance Heatmaps ---
message("Calculating and saving sample distance/correlation matrices...")

# Euclidean Distance
sampleDists_euc <- dist(t(assay(vsd))) # Transpose VST matrix for sample distances
sampleDistMatrix_euc <- as.matrix(sampleDists_euc)
# rownames(sampleDistMatrix_euc) <- vsd$condition # REMOVED THIS LINE
# colnames(sampleDistMatrix_euc) <- NULL         # REMOVED THIS LINE
colors_euc <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Create the annotation dataframe once
# Ensure column names exist in colData(vsd)
annotation_cols_to_select <- intersect(c("condition", "timepoint"), colnames(colData(vsd)))
if(length(annotation_cols_to_select) == 0) {
  warning("Annotation columns 'condition' or 'timepoint' not found in colData.")
  annotation_df <- NULL # Or handle as appropriate
} else {
  annotation_df <- as.data.frame(colData(vsd)[, annotation_cols_to_select, drop=FALSE])
}


# Save Euclidean Distance Matrix (CSV)
message("Saving Euclidean distance matrix data to: ", euclidean_heatmap_data_file)
# Ensure row names are kept in the CSV if needed, default is FALSE for write.csv unless row.names=TRUE
write.csv(sampleDistMatrix_euc, euclidean_heatmap_data_file, row.names = TRUE) # Added row.names = TRUE explicitly

# Generate and save Euclidean heatmap plot (JPG) - CORRECTED
message("Generating and saving Euclidean distance heatmap...")
if (!is.null(annotation_df)) { # Only add annotation if available
  pheatmap(sampleDistMatrix_euc,
           clustering_distance_rows=sampleDists_euc,
           clustering_distance_cols=sampleDists_euc,
           col=colors_euc,
           main="Sample Euclidean Distances (VST Data)",
           annotation_col = annotation_df, # Use the prepared annotation dataframe
           filename = euclidean_heatmap_plot_file,
           width = 8, height = 7)
} else {
  pheatmap(sampleDistMatrix_euc,
           clustering_distance_rows=sampleDists_euc,
           clustering_distance_cols=sampleDists_euc,
           col=colors_euc,
           main="Sample Euclidean Distances (VST Data)",
           # No annotation_col argument
           filename = euclidean_heatmap_plot_file,
           width = 8, height = 7)
  message("Annotation columns not found, heatmap generated without column annotations.")
}
message("Euclidean distance heatmap saved to: ", euclidean_heatmap_plot_file)


# Pearson Correlation
sampleCor_pearson <- cor(assay(vsd), method="pearson") # Use VST matrix directly for correlation
sampleCorDist_pearson <- as.dist(1 - sampleCor_pearson) # Convert correlation to distance
colors_cor <- colorRampPalette( brewer.pal(9, "YlOrRd") )(255) # Different color scheme for correlation

# Save Pearson Correlation Matrix (CSV)
message("Saving Pearson correlation matrix data to: ", pearson_heatmap_data_file)
write.csv(sampleCor_pearson, pearson_heatmap_data_file, row.names = TRUE) # Added row.names = TRUE explicitly

# Generate and save Pearson Correlation heatmap plot (JPG) - Apply same annotation logic
message("Generating and saving Pearson correlation heatmap...")
if (!is.null(annotation_df)) { # Only add annotation if available
  pheatmap(sampleCor_pearson,
           clustering_distance_rows = sampleCorDist_pearson,
           clustering_distance_cols = sampleCorDist_pearson,
           col = colors_cor,
           main = "Sample Pearson Correlation (VST Data)",
           annotation_col = annotation_df, # Use the prepared annotation dataframe
           filename = pearson_heatmap_plot_file,
           width = 8, height = 7)
} else {
  pheatmap(sampleCor_pearson,
           clustering_distance_rows = sampleCorDist_pearson,
           clustering_distance_cols = sampleCorDist_pearson,
           col = colors_cor,
           main = "Sample Pearson Correlation (VST Data)",
           # No annotation_col argument
           filename = pearson_heatmap_plot_file,
           width = 8, height = 7)
  message("Annotation columns not found, heatmap generated without column annotations.")
}
message("Pearson correlation heatmap saved to: ", pearson_heatmap_plot_file)


# --- 4. Output Normalized Counts ---

message("Preparing and saving normalized count matrix...")
# Convert matrix to dataframe, add transcript_id as a column
normalized_counts_df <- as.data.frame(normalized_counts) %>%
  rownames_to_column("transcript_id")

# Merge with original gene annotation (symbol)
# Need to ensure gene_annotation only has unique transcript_ids if counts were filtered
gene_annotation_filtered <- gene_annotation %>% dplyr::filter(transcript_id %in% normalized_counts_df$transcript_id)
normalized_counts_annotated <- dplyr::left_join(gene_annotation_filtered, normalized_counts_df, by = "transcript_id")

# Save the annotated normalized counts
message("Saving normalized counts to: ", normalized_counts_file)
readr::write_csv(normalized_counts_annotated, normalized_counts_file)

# --- 5. Differential Expression Analysis ---

message("Performing differential expression analysis for specified comparisons...")

# Define the comparisons based on the 'condition' factor levels
# Format: results(dds, contrast=c("factor_name", "numerator_level", "denominator_level"))
comparisons <- list(
  "MCF7_021_10_vs_0" = c("condition", "MCF7_021_10", "MCF7_021_0"),
  "MCF7_021_24_vs_0" = c("condition", "MCF7_021_24", "MCF7_021_0"),
  "MCF7_021_48_vs_0" = c("condition", "MCF7_021_48", "MCF7_021_0"),
  "MCF7_023_24_vs_0" = c("condition", "MCF7_023_24", "MCF7_023_0"),
  "MCF7_025_24_vs_0" = c("condition", "MCF7_025_24", "MCF7_025_0")
)

# Loop through comparisons, extract results, filter, and save
deg_results_list <- list() # Store results for filtering

# Ensure gene_annotation is ready for merging within the loop
gene_annotation_for_merge <- gene_annotation %>% distinct(transcript_id, .keep_all = TRUE)

for (comp_name in names(comparisons)) {
  message("  Processing comparison: ", comp_name)
  contrast_vec <- comparisons[[comp_name]]
  
  # Extract results
  res <- results(dds, contrast=contrast_vec, alpha=0.05) # Use alpha=0.05 for independent filtering
  
  # Optional: Shrink log2 fold changes (good for visualization and ranking, less crucial for basic lists)
  # Requires apeglm package: BiocManager::install("apeglm")
  # try({ # Use try in case apeglm is not installed or shrinkage fails
  #   library(apeglm)
  #   res <- lfcShrink(dds, contrast=contrast_vec, type="apeglm", res=res) # Pass original res object
  # }, silent = TRUE)
  
  # Order results by adjusted p-value
  res_ordered <- res[order(res$padj),]
  
  # Convert to dataframe and add annotations
  res_df <- as.data.frame(res_ordered) %>%
    rownames_to_column("transcript_id") %>%
    dplyr::left_join(gene_annotation_for_merge, by = "transcript_id") %>% # Use distinct annotation
    dplyr::select(transcript_id, symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) # Reorder columns
  
  # Store for later filtering
  deg_results_list[[comp_name]] <- res_df
  
  # Save the full DEG results table
  output_file <- deg_output_files[[comp_name]]
  message("    Saving full DEG results to: ", output_file)
  readr::write_csv(res_df, output_file)
}


# --- 6. Output UpDEGs and DownDEGs ---

message("Filtering and saving Up-regulated and Down-regulated DEGs...")

for (comp_name in names(deg_results_list)) {
  message("  Filtering DEGs for: ", comp_name)
  res_df <- deg_results_list[[comp_name]]
  
  # Define criteria
  padj_threshold <- 0.05
  lfc_threshold_up <- log2(1.5)
  lfc_threshold_down <- log2(1/1.5)
  
  # Filter UpDEGs
  up_degs <- res_df %>%
    dplyr::filter(padj <= padj_threshold & log2FoldChange >= lfc_threshold_up & !is.na(padj)) %>%
    dplyr::arrange(padj, desc(log2FoldChange)) # Sort by significance then fold change
  
  output_file_up <- updeg_output_files[[comp_name]]
  message("    Saving UpDEGs (padj <= ", padj_threshold, ", LFC >= ", lfc_threshold_up, ") to: ", output_file_up)
  readr::write_csv(up_degs, output_file_up)
  
  # Filter DownDEGs
  down_degs <- res_df %>%
    dplyr::filter(padj <= padj_threshold & log2FoldChange <= lfc_threshold_down & !is.na(padj)) %>%
    dplyr::arrange(padj, log2FoldChange) # Sort by significance then fold change
  
  output_file_down <- downdeg_output_files[[comp_name]]
  message("    Saving DownDEGs (padj <= ", padj_threshold, ", LFC <= ", lfc_threshold_down, ") to: ", output_file_down)
  readr::write_csv(down_degs, output_file_down)
}

message("--- Analysis Script Finished ---")