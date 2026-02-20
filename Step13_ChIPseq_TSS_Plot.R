setwd("C:/Users/ski2959/Desktop/NIH R21/Codes/src")

# --- 1. SETUP: INSTALL AND LOAD NECESSARY PACKAGES ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "GenomicRanges")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

cran_packages <- c("data.table", "ggplot2")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(data.table)
library(ggplot2)

cat("All necessary packages are installed and loaded.\n\n")

# --- 2. DEFINE FILE PATHS AND STANDARD CHROMOSOMES ---
# Note: Ensure these relative paths are correct for your directory structure.
file_paths <- c(
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2.broadPeak",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep1.broadPeak",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_0_5_ug_mL_Rep2.broadPeak",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_0_5_ug_mL_Rep1.broadPeak",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_0_1_ug_mL_Rep2.broadPeak",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_0_1_ug_mL_Rep1.broadPeak",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_0_0_ug_mL_Rep2.broadPeak",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_0_0_ug_mL_Rep1.broadPeak",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2.broadPeak",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep1.broadPeak",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_0_5_ug_mL_Rep2.broadPeak",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_0_5_ug_mL_Rep1.broadPeak"
)

standard_chroms <- c(paste0("chr", 1:22), "chrX", "chrY")

# --- 3. PROCESS EACH FILE: IMPORT, FILTER, AND PLOT ---
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

cat("Starting to process", length(file_paths), "files...\n")

for (file_path in file_paths) {
  sample_name <- basename(file_path)
  sample_name <- sub(".broadPeak", "", sample_name)
  
  cat("\n--- Processing:", sample_name, "---\n")
  
  # Import peak file (first three columns: chrom, start, end)
  peak_data <- tryCatch({
    fread(file_path, select = 1:3, col.names = c("chrom", "start", "end"))
  }, error = function(e) {
    cat("Error: Cannot read file ->", file_path, "\nSkipping this file.\n")
    return(NULL)
  })
  
  if (is.null(peak_data)) next
  
  cat("Successfully imported", nrow(peak_data), "total peaks.\n")
  
  # Filter for standard chromosomes
  filtered_peaks <- peak_data[chrom %in% standard_chroms]
  num_filtered_peaks <- nrow(filtered_peaks)
  cat("Found", num_filtered_peaks, "peaks on standard chromosomes.\n")
  
  if (num_filtered_peaks < 1) {
    cat("Skipping plot generation: no peaks left after filtering.\n")
    next
  }
  
  # Convert to GRanges
  peak_granges <- makeGRangesFromDataFrame(filtered_peaks,
                                           keep.extra.columns = FALSE,
                                           seqnames.field = "chrom",
                                           start.field = "start",
                                           end.field = "end")
  
  # Set up promoter regions
  promoter_regions <- getPromoters(TxDb = txdb, upstream = 10000, downstream = 10000)
  
  # Harmonize seqlevels (chromosome names) between peaks and the annotation database
  common_seqlevels <- intersect(seqlevels(peak_granges), seqlevels(promoter_regions))
  peak_granges <- keepSeqlevels(peak_granges, common_seqlevels, pruning.mode = "coarse")
  promoter_regions <- keepSeqlevels(promoter_regions, common_seqlevels, pruning.mode = "coarse")
  
  if (length(peak_granges) < 1) {
    cat("Skipping plot: no peaks remained after harmonizing chromosomes with the TxDb.\n")
    next
  }
  
  # Get tag matrix
  tagMatrix <- tryCatch({
    getTagMatrix(peak_granges, windows = promoter_regions)
  }, error = function(e) {
    cat("Error during getTagMatrix for", sample_name, ":", e$message, "\nSkipping this file.\n")
    return(NULL)
  })
  
  # === REVISED CHECK ===
  # Added a more robust check to ensure the output is a valid, non-empty matrix.
  # This prevents errors if getTagMatrix returns NULL or an empty/malformed object.
  if (!is.matrix(tagMatrix) || nrow(tagMatrix) < 1) {
    cat("Skipping plot for", sample_name, "as no valid tag matrix was generated (likely no peaks in promoter regions).\n")
    next
  }
  
  cat("Generated tag matrix with", nrow(tagMatrix), "rows.\n")
  
  # Define plot limits based on the actual matrix dimensions
  tag_ncol <- ncol(tagMatrix)
  xlim_vals <- c(-ceiling((tag_ncol - 1) / 2), floor((tag_ncol - 1) / 2))
  cat("Using xlims: [", xlim_vals[1], ",", xlim_vals[2], "] for plotting.\n")
  
  # Plot the average profile
  tss_plot <- plotAvgProf(tagMatrix, xlim = xlim_vals,
                          xlab = "Genomic Region (5'->3')",
                          ylab = "Read Count Frequency") +
    ggtitle(paste("TSS Enrichment for", sub("SHK_2024_06_11_", "", sample_name))) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Save the plot
  output_filename <- file.path("../public/ChIPseq", paste0(sample_name, "_TSS_Plot.jpg"))
  ggsave(output_filename, plot = tss_plot, width = 8, height = 6, dpi = 300, device = "jpeg")
  cat("Successfully saved TSS plot to:", output_filename, "\n")
}

cat("\nProcessing complete.\n")
