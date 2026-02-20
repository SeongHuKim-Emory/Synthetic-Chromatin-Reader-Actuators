setwd("C:/Users/ski2959/Desktop/NIH R21/Codes/src")

################################################################################

# R Script for Processing and Exporting RNA-Seq DEG Categories

# This script reads a CSV file containing differentially expressed genes (DEGs),
# processes the data to create a unified category, cleans up gene symbols,
# selects specific columns, and then exports the data into separate BED-like
# files for each category.

# Step 0: Install and load necessary packages
# If you don't have the 'tidyverse' package installed, uncomment the next line and run it once.
# install.packages("tidyverse")

# Load the required libraries for data manipulation.
library(dplyr)
library(tidyr)

# --- Step 1: Import the CSV file ---
# We use read.csv to load the data into a data frame.
# 'check.names = FALSE' is important to prevent R from changing column names
# that contain spaces or special characters (e.g., "Category - SRA specific").
cat("Step 1: Importing data from UpDEG_category.csv...\n")
tryCatch({
  updeg_data <- read.csv("../public/RNAseq/UpDEG_category.csv", check.names = FALSE)
}, error = function(e) {
  stop("Error: The file '../public/RNAseq/UpDEG_category.csv' was not found. Please ensure the file path is correct.")
})


# --- Step 2: Duplicate the table ---
# We create a new data frame to work with, preserving the original data.
cat("Step 2: Duplicating the initial table...\n")
processed_data <- updeg_data


# --- Step 3: Create a new unified 'Category' column ---
# This step merges the 'Category - SRA specific' and 'Category - other' columns.
# First, we replace any empty strings "" with NA (Not Available) to make them easier to handle.
processed_data$`Category - SRA specific`[processed_data$`Category - SRA specific` == ""] <- NA
processed_data$`Category - other`[processed_data$`Category - other` == ""] <- NA

# We then use `mutate` to add a new column named 'Category'.
# `coalesce` is a handy function that takes the first non-NA value from the given columns for each row.
cat("Step 3: Merging category columns into a single 'Category' column...\n")
processed_data <- processed_data %>%
  mutate(Category = coalesce(`Category - SRA specific`, `Category - other`))


# --- Step 4: Fill in missing 'Symbol' values ---
# If a gene's 'Symbol' is missing (empty or NA), we fill it with its 'transcript_id'.
# `ifelse` checks the condition for each row and populates the 'Symbol' column accordingly.
cat("Step 4: Filling empty 'Symbol' values with 'transcript_id'...\n")
processed_data <- processed_data %>%
  mutate(Symbol = ifelse(is.na(Symbol) | Symbol == "", transcript_id, Symbol))


# --- Step 5: Select and reorder the final columns ---
# We select only the columns needed for the final output and discard the rest.
# The columns are also arranged in the desired order.
cat("Step 5: Selecting and reordering final columns...\n")
final_data <- processed_data %>%
  select(Symbol, hg38_chr, hg38_start, hg38_end, Category)


# --- Step 6: Sort the table by the 'Category' column ---
# This step is not strictly required for splitting the file but is good practice
# and was requested. It sorts the entire dataset based on the category names.
cat("Step 6: Sorting the final table by 'Category'...\n")
final_data <- final_data %>%
  arrange(Category)

# --- Step 7: Export data into separate files for each category ---
cat("Step 7: Exporting data into separate files per category...\n")

# Define the output directory
output_dir <- "../public/RNAseq/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Get a list of unique categories
categories <- unique(final_data$Category)

# Loop through each category
for (current_category in categories) {
  # Subset the data for the current category
  subset_data <- final_data %>%
    filter(Category == current_category)
  
  # --- NEW: Add "chr" prefix to the chromosome column ---
  # This converts the numeric chromosome to the required character format (e.g., 8 -> "chr8").
  output_data <- subset_data %>%
    mutate(hg38_chr = paste0("chr", hg38_chr)) %>%
    select(hg38_chr, hg38_start, hg38_end, Symbol)
  
  # Define the output filename
  filename <- file.path(output_dir, paste0("UpDEG_", current_category, ".bed"))
  
  # Use write.table to export the reordered data
  write.table(output_data,
              file = filename,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  cat(paste("  - Successfully exported:", filename, "\n"))
}

cat("\nProcessing complete. All files have been exported.\n")

################################################################################

#
# ChIP-seq Peak Binning and Counting Script
#
# This script reads ChIP-seq peak data from .broadPeak and .bed files,
# divides the human genome (hg38) into 10 megabase (Mb) bins, counts the
# number of peaks falling into each bin, and exports the results to CSV files.
#

# --- 1. SETUP ---

# Install and load the data.table package if not already installed.
# data.table is used for its high performance on large datasets.
if (!require("data.table")) {
  install.packages("data.table")
}
library(data.table)

# Prevent R from using scientific notation for large numbers (e.g., 1.1e+08).
# This ensures chromosome positions are written in full.
options(scipen = 999)

# Define the hg38 chromosome lengths for chr1-22 and chrX.
# This also defines the desired output order for chromosomes.
hg38_lengths <- data.table(
  chr = c(paste0("chr", 1:22), "chrX"),
  length = c(
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
    58617616, 64444167, 46709983, 50818468, 156040895
  )
)

# Define the size of each genomic bin.
bin_size <- 10000000

# --- 2. DEFINE THE CORE PROCESSING FUNCTION ---

#' Process a ChIP-seq file to count peaks in genomic bins.
#'
#' @param input_path character Path to the input .broadPeak or .bed file.
#' @param output_path character Path for the output .csv file.
#' @param hg38_ref data.table Reference data.table with chromosome names and lengths.
#' @param bin_width integer The width of each genomic bin in base pairs.
#'
#' @return The function writes a CSV file to the output_path and prints a status message.
process_chip_file <- function(input_path, output_path, hg38_ref, bin_width) {
  
  # Print a message to indicate which file is being processed.
  cat("Processing:", basename(input_path), "\n")
  
  # --- Step 2a: Read and Prepare Peak Data ---
  # Read the first three columns (chr, start, end) from the input file.
  peaks_dt <- fread(input_path,
                    select = 1:3,
                    col.names = c("chr", "start", "end"))
  
  # Convert 0-based start coordinates from BED/broadPeak to 1-based for R.
  peaks_dt[, start := start + 1]
  
  # Standardize chromosome names to have the 'chr' prefix.
  peaks_dt[!(chr %like% "^chr"), chr := paste0("chr", chr)]
  
  # Filter to keep only the standard chromosomes (1-22 and X).
  peaks_dt <- peaks_dt[chr %in% hg38_ref$chr]
  
  # Check if any peaks remain after filtering.
  if (nrow(peaks_dt) == 0) {
    cat("Warning: No peaks found on standard chromosomes (chr1-22, chrX) for this file.\n")
  }
  
  # --- Step 2b: Generate Genomic Bins ---
  # Create a list to hold the bin data.tables for each chromosome.
  all_bins_list <- list()
  
  # Loop through each chromosome in our reference table.
  for (i in 1:nrow(hg38_ref)) {
    chrom_name <- hg38_ref$chr[i]
    chrom_length <- hg38_ref$length[i]
    
    # Generate the starting positions for each bin on the chromosome (1-based).
    bin_starts <- seq(from = 1, to = chrom_length, by = bin_width)
    
    # Generate the ending positions for each bin.
    bin_ends <- bin_starts + bin_width - 1
    
    # Ensure the last bin does not extend past the end of the chromosome.
    bin_ends[length(bin_ends)] <- chrom_length
    
    # Create a data.table for the bins of the current chromosome.
    chrom_bins <- data.table(
      chr = chrom_name,
      start = bin_starts,
      end = bin_ends
    )
    
    # Add the generated bins to our list.
    all_bins_list[[i]] <- chrom_bins
  }
  
  # Combine the list of data.tables into one large data.table of all genome bins.
  genome_bins_dt <- rbindlist(all_bins_list)
  
  # --- Step 2c: Count Overlapping Peaks ---
  # This is the core step where peaks are assigned to bins.
  # We use foverlaps, a high-performance function for finding interval overlaps.
  
  # Set keys on both data.tables for the overlap join.
  setkey(peaks_dt, chr, start, end)
  setkey(genome_bins_dt, chr, start, end)
  
  # Find all overlaps between peaks (x) and bins (y).
  # The resulting table contains the columns of the bins, plus overlap info.
  overlaps <- foverlaps(peaks_dt, genome_bins_dt, type = "any", nomatch = 0L)
  
  # Count how many peaks fall into each unique bin.
  # We group by the bin's coordinates ('start', 'end').
  counts_dt <- overlaps[, .N, by = .(chr, start, end)]
  
  # Rename the new count column 'N' to 'peak_count' for clarity.
  setnames(counts_dt, "N", "peak_count")
  
  # --- Step 2d: Merge Counts with Bins to Include Zero-Count Bins ---
  # The 'counts_dt' table only has bins with at least one peak.
  # We need to merge it back with the full 'genome_bins_dt' to include
  # bins that have a count of zero.
  
  # Perform a left join from the complete set of bins to the counts table.
  final_dt <- merge(genome_bins_dt, counts_dt, by = c("chr", "start", "end"), all.x = TRUE)
  
  # For bins that had no overlapping peaks, the 'peak_count' will be NA.
  # Replace these NA values with 0.
  final_dt[is.na(peak_count), peak_count := 0]
  
  # --- Step 2e: Sort and Export the Final Table ---
  
  # **REVISION: Enforce Chromosome Order**
  # Define the desired chromosome order based on the reference table.
  chromosome_order <- hg38_ref$chr
  
  # Convert the 'chr' column to a factor with the specified order.
  # This allows for correct sorting.
  final_dt[, chr := factor(chr, levels = chromosome_order)]
  
  # Sort the data table by the specified chromosome order, and then by start position.
  setorder(final_dt, chr, start)
  
  # Ensure the output directory exists.
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write the final data.table to the specified CSV file.
  fwrite(final_dt, file = output_path)
  
  cat("Finished. Output saved to:", output_path, "\n\n")
}

# --- 3. RUN THE PROCESSING FOR ALL FILES ---

# Create a list of file paths to process.
files_to_process <- list(
  list(
    input = "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2.broadPeak",
    output = "../public/Multiomics/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2_10Mb.csv"
  ),
  list(
    input = "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2.broadPeak",
    output = "../public/Multiomics/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2_10Mb.csv"
  ),
  list(
    input = "../public/Public Dataset/ChIPseq/H3K27me3/GSE96363 ENCFF669NUD hg38 MCF7 H3K27me3 ChIPseq.bed",
    output = "../public/Multiomics/GSE96363 ENCFF669NUD hg38 MCF7 H3K27me3 ChIPseq_10Mb.csv"
  ),
  list(
    input = "../public/Public Dataset/ChIPseq/H3K4me3/GSE96506 ENCFF625APE hg38 MCF7 H3K4me3 ChIPseq.bed",
    output = "../public/Multiomics/GSE96506 ENCFF625APE hg38 MCF7 H3K4me3 ChIPseq_10Mb.csv"
  ),
  list(
    input = "../public/RNAseq/UpDEG_category.bed",
    output = "../public/Multiomics/UpDEG_category_10Mb.csv"
  ),
  list(
    input = "../public/RNAseq/UpDEG_Delayed.bed",
    output = "../public/Multiomics/UpDEG_Delayed_10Mb.csv"
  ),
  list(
    input = "../public/RNAseq/UpDEG_Early sustained.bed",
    output = "../public/Multiomics/UpDEG_Early sustained_10Mb.csv"
  ),
  list(
    input = "../public/RNAseq/UpDEG_Early transient.bed",
    output = "../public/Multiomics/UpDEG_Early transient_10Mb.csv"
  ),
  list(
    input = "../public/RNAseq/UpDEG_Late sustained.bed",
    output = "../public/Multiomics/UpDEG_Late sustained_10Mb.csv"
  ),
  list(
    input = "../public/RNAseq/UpDEG_Late transient.bed",
    output = "../public/Multiomics/UpDEG_Late transient_10Mb.csv"
  ),
  list(
    input = "../public/RNAseq/UpDEG_ns.bed",
    output = "../public/Multiomics/UpDEG_ns_10Mb.csv"
  ),
  list(
    input = "../public/RNAseq/UpDEG_PCD-RFP & SRA.bed",
    output = "../public/Multiomics/UpDEG_PCD-RFP & SRA_10Mb.csv"
  ),
  list(
    input = "../public/RNAseq/UpDEG_PCD-RFP.bed",
    output = "../public/Multiomics/UpDEG_PCD-RFP_10Mb.csv"
  ),
  list(
    input = "../public/RNAseq/UpDEG_Stochastic.bed",
    output = "../public/Multiomics/UpDEG_Stochastic_10Mb.csv"
  )
)

# Loop through the list and call the processing function for each file.
for (file_set in files_to_process) {
  process_chip_file(
    input_path = file_set$input,
    output_path = file_set$output,
    hg38_ref = hg38_lengths,
    bin_width = bin_size
  )
}

cat("All files have been processed.\n")