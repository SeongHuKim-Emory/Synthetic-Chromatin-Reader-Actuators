# Load the data.table library for efficient handling of large datasets.
# If you don't have it installed, run: install.packages("data.table")
library(data.table)

# --- 1. Import Files ---

# Define the file paths for the input data
gene_assoc_file <- "../public/Public Dataset/GeneHancer/GeneHancer_AnnotSV_gene_association_scores_v5.25.txt"
element_file <- "../public/Public Dataset/GeneHancer/GeneHancer_AnnotSV_elements_v5.25.txt"

# Read the two text files into data.table objects for performance
gene_associations <- fread(gene_assoc_file)
elements <- fread(element_file)

# --- 2. Create a table named GeneHancer_Genes ---
# --- 3. Copy content to the new table ---

# The 'gene_associations' data.table will be used as the base for the new table.
# We will create a copy to ensure the original data remains unchanged.
GeneHancer_Genes <- copy(gene_associations)

# --- 4. Rename the 'is_elite' column ---

# Rename the 'is_elite' column to 'is_gene_elite' to distinguish it
# from the 'is_elite' column in the elements file.
setnames(GeneHancer_Genes, old = "is_elite", new = "is_gene_elite")

# --- 5 & 6. Add and Populate Columns by Merging Data ---

# Instead of adding empty columns and then looping through to populate them,
# a more efficient and standard approach in R is to perform a merge (or join) operation.
# This single step achieves the goal of finding matching rows by 'GHid' and
# copying the corresponding values from the elements table.

# We perform a 'left join' to ensure all records from GeneHancer_Genes are kept.
# The 'by = "GHid"' argument specifies the common column to match rows.
GeneHancer_Genes <- merge(GeneHancer_Genes, elements, by = "GHid", all.x = TRUE)

# The result is the final table with all the requested columns populated correctly.
# Note: The merge automatically handles the addition and population of columns:
# 'chr', 'element_start', 'element_end', 'is_elite', 'regulatory_element_type', 'enhancer_score'

# --- 7. Export the final table ---

# Define the path for the output file
output_file_path <- "../public/Public Dataset/GeneHancer/GeneHancer_Genes_Elements.csv"

# Write the final merged table to a CSV file.
# 'fwrite' is a fast way to write data frames to disk.
# 'row.names = FALSE' is used to prevent R from writing an extra column for row numbers.
fwrite(GeneHancer_Genes, file = output_file_path, row.names = FALSE)

# --- Completion Message ---
# Print a message to the console to confirm that the script has finished successfully.
cat("Script finished successfully.\n")
cat("The merged data has been exported to:", output_file_path, "\n")

# You can optionally view the first few rows of the final table by uncommenting the line below:
# print(head(GeneHancer_Genes))

################################################################################

################################################################################

# --- GeneHancer Data Processing Script ---
#
# This script imports GeneHancer data, filters it to retain only entries
# associated with protein-coding genes from the hg38 genome, sorts the data,
# and exports the result to a new CSV file.

# --- 1. Setup: Install and Load Required Packages ---

# The 'dplyr' package is used for efficient data manipulation (sorting and filtering).
# The 'biomaRt' package is used to connect to biological databases like Ensembl.

# Check for BiocManager and install if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of required packages
required_packages <- c("dplyr", "biomaRt")

# Install missing packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg == "biomaRt") {
      BiocManager::install("biomaRt")
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

cat("Packages loaded successfully.\n")


# --- 2. Define File Paths ---

# Set the path for the input and output files.
input_file <- "../public/Public Dataset/GeneHancer/GeneHancer_Genes_Elements.csv"
output_file <- "../public/Public Dataset/GeneHancer/GeneHancer_Protein_Coding_Genes_Elements.csv"


# --- 3. Retrieve Protein-Coding Genes from Ensembl (hg38) ---

cat("Connecting to Ensembl to retrieve the list of protein-coding genes...\n")

# Connect to the Ensembl genes database. The 'hsapiens_gene_ensembl' dataset
# defaults to the latest human genome build, which is GRCh38 (hg38).
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

cat("Fetching gene symbols. This may take a moment...\n")

# Retrieve all gene symbols ('hgnc_symbol') that have the "protein_coding" biotype.
protein_coding_genes <- getBM(
  attributes = 'hgnc_symbol',
  filters = 'biotype',
  values = 'protein_coding',
  mart = ensembl
)

# Extract the symbols into a unique character vector for efficient filtering.
# Also, remove any blank entries that might be returned from the database.
protein_coding_symbols <- unique(protein_coding_genes$hgnc_symbol)
protein_coding_symbols <- protein_coding_symbols[protein_coding_symbols != "" & !is.na(protein_coding_symbols)]

cat("Retrieved", length(protein_coding_symbols), "unique protein-coding gene symbols.\n")


# --- 4. Import and Process GeneHancer Data ---

cat("Importing GeneHancer data from:", input_file, "\n")

# Use a tryCatch block to handle potential errors, like the file not being found.
tryCatch({
  # Step 1: Import the file
  genehancer_data <- read.csv(input_file)
  cat("Successfully imported", nrow(genehancer_data), "rows.\n")
  
  cat("Processing data: Sorting by gene symbol and filtering...\n")
  
  # Use dplyr's pipe operator (%>%) for a clean workflow:
  # Step 2: Sort the entire table by the 'symbol' column alphabetically.
  # Step 3: Filter the table, keeping only rows where the 'symbol' value
  #         is present in our list of protein-coding gene symbols.
  processed_data <- genehancer_data %>%
    arrange(symbol) %>%
    filter(symbol %in% protein_coding_symbols)
  
  cat("Processing complete.\n")
  cat("Original number of rows:", nrow(genehancer_data), "\n")
  cat("Number of rows after filtering for protein-coding genes:", nrow(processed_data), "\n")
  
  # --- 5. Export the Filtered Data ---
  
  # Step 4: Export the updated table.
  cat("Exporting the processed data to:", output_file, "\n")
  write.csv(processed_data, output_file, row.names = FALSE, quote = TRUE)
  cat("Data successfully exported.\n")
  
}, error = function(e) {
  # This code will run if an error occurs in the tryCatch block.
  stop(paste("An error occurred during file processing.\n",
             "Please check that the input file exists at the specified path.\n",
             "Path:", input_file, "\n",
             "Original error message:", e$message))
})

cat("Script finished.\n")

################################################################################

################################################################################

# GeneHancer Data Processing Script
# This script reads the GeneHancer dataset, processes it by gene symbol,
# formats the genomic coordinates, and exports the results into separate
# BED-like files for each gene symbol.

# Load the data.table library for efficient data manipulation.
# If you don't have it installed, run: install.packages("data.table")
library(data.table)

# 1. Import file
# Define the input file path.
input_file <- "../public/Public Dataset/GeneHancer/GeneHancer_Protein_Coding_Genes_Elements.csv"

# Check if the file exists before attempting to read it.
if (!file.exists(input_file)) {
  stop(paste("Error: Input file not found at", input_file))
}

# Read the CSV file into a data.table for performance.
# The fread function is generally faster and more memory-efficient than read.csv.
message("Step 1: Reading the dataset...")
gene_data <- fread(input_file)
message("Dataset read successfully.")

# Get the list of unique gene symbols to iterate over.
unique_symbols <- unique(gene_data$symbol)
message(paste("Found", length(unique_symbols), "unique gene symbols."))

# Define the output directory.
output_dir <- "../public/Public Dataset/GeneHancer/"

# Create the output directory if it doesn't exist to prevent errors.
if (!dir.exists(output_dir)) {
  message(paste("Creating output directory:", output_dir))
  dir.create(output_dir, recursive = TRUE)
}

# 2-5. Process and Export Data for Each Symbol
message("Steps 2-5: Processing and exporting data for each symbol...")

# Loop through each unique symbol to process and write its corresponding data.
for (current_symbol in unique_symbols) {
  
  # 2. For each unique symbol, gather the relevant columns.
  # Filter the data.table to get rows for the current symbol.
  # Select only the 'chr', 'element_start', and 'element_end' columns.
  symbol_subset <- gene_data[symbol == current_symbol, .(chr, element_start, element_end)]
  
  # 3. Add "chr" to the front of the "chr" column values.
  # The `:=` operator in data.table modifies the column in place.
  symbol_subset[, chr := paste0("chr", chr)]
  
  # 5. Define the full path for the output file.
  output_file_path <- file.path(output_dir, paste0(current_symbol, "_Enhancer.bed"))
  
  # 4. Export the symbol's information to the specified path.
  # We use fwrite for fast writing.
  # - `col.names = FALSE` removes the header row from the output file.
  # - `sep = "\t"` sets the delimiter to a tab, making it a standard BED file.
  # - `quote = FALSE` ensures values are not surrounded by quotes.
  fwrite(symbol_subset, file = output_file_path, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  
}

message("Processing complete.")
message(paste("All files have been exported to:", output_dir))

################################################################################

################################################################################

# Title: Optimized ChIP-seq Peak Overlap Significance Analysis using Parallel Permutation Testing
# Description: This script imports ChIP-seq peak data from BED files,
#              conducts a parallelized permutation test to assess the significance of overlap
#              between a primary dataset and several public datasets,
#              prints the results, and visualizes the findings.
#              This version is optimized for speed using the BiocParallel package.
#              All console output is saved to a summary CSV file.

# --- 1. Installation and Loading of Required Packages ---

# Check if BiocManager is installed, and install it if not.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define a list of required Bioconductor packages.
# Added "regioneR" for permutation testing and "BiocParallel" for parallelization.
required_packages <- c("ChIPpeakAnno", "ChIPseeker", "GenomicRanges", "TxDb.Hsapiens.UCSC.hg38.knownGene", "rtracklayer", "regioneR", "BiocParallel")

# Loop through the packages to check for installation and install if missing.
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing", pkg, "..."))
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# --- 1.5. Initialize Report Object ---
# This character vector will store all the console output lines for the final CSV report.
csv_report_lines <- c()


# --- 2. Import BED Files as GRanges Objects ---

# Define file paths
single_bed_file <- "../public/Public Dataset/GeneHancer/HOXB9_Enhancer.bed"

# Assign a name for the single bed file dataset for use in console output
single_bed_file_name <- "HOXB9 Enhancers"

bed_files <- c(
  H3K4me1 = "../public/Public Dataset/ChIPseq/GSE86714 ENCFF991HJA hg38 MCF7 H3K4me1 ChIPseq Filtered_trim.bed",
  H3K4me2 = "../public/Public Dataset/ChIPseq/GSE96439 ENCFF188VRU hg38 MCF7 H3K4me2 ChIPseq Filtered_trim.bed",
  H3K4me3 = "../public/Public Dataset/ChIPseq/GSE96506 ENCFF268RXB hg38 MCF7 H3K4me3 ChIPseq Filtered_trim.bed",
  H3K9ac = "../public/Public Dataset/ChIPseq/GSE95898 ENCFF348DEB hg38 MCF7 H3K9ac ChIPseq Filtered_trim.bed",
  H3K9me3 = "../public/Public Dataset/ChIPseq/GSE96517 ENCFF501UHK hg38 MCF7 H3K9me3 ChIPseq Filtered_trim.bed",
  H3K27ac = "../public/Public Dataset/ChIPseq/GSE96352 ENCFF491LQY hg38 MCF7 H3K27ac ChIPseq Filtered_trim.bed",
  H3K27me3 = "../public/Public Dataset/ChIPseq/GSE96363 ENCFF669NUD hg38 MCF7 H3K27me3 ChIPseq Filtered_trim.bed",
  H3K36me3 = "../public/Public Dataset/ChIPseq/GSE174945 ENCFF195FSD hg38 MCF7 H3K36me3 ChIPseq Filtered_trim.bed",
  H4K20me1 = "../public/Public Dataset/ChIPseq/GSE96283 ENCFF714DEQ hg38 MCF7 H4K20me1 ChIPseq Filtered_trim.bed"
)

# Import the primary BED file
cat("Importing primary ChIP-seq peak file...\n")
single_bed_file_peaks <- rtracklayer::import(single_bed_file, format = "BED")
primary_peak_info <- paste("Number of peaks in", single_bed_file_name, "dataset:", length(single_bed_file_peaks))
cat(primary_peak_info, "\n\n")

# Capture primary file import details for the CSV report
csv_report_lines <- c(csv_report_lines,
                      "Importing primary ChIP-seq peak file...",
                      primary_peak_info,
                      "") # Add a blank line for spacing

# Import the public BED files and store them in a list
bed_files_peaks_list <- lapply(bed_files, function(file) {
  cat("Importing:", basename(file), "\n")
  peaks <- rtracklayer::import(file, format = "BED")
  cat("Number of peaks:", length(peaks), "\n")
  return(peaks)
})
cat("\n--- All files imported successfully. ---\n\n")

# Capture public file import details for the CSV report
for (name in names(bed_files_peaks_list)) {
  file_path <- bed_files[name]
  peak_count <- length(bed_files_peaks_list[[name]])
  import_msg <- paste("Importing:", basename(file_path))
  num_peaks_msg <- paste("Number of peaks:", peak_count)
  csv_report_lines <- c(csv_report_lines, import_msg, num_peaks_msg)
}
csv_report_lines <- c(csv_report_lines, "", "--- All files imported successfully. ---", "")


# --- 3. Perform Permutation Test in Parallel ---

# Set a seed for reproducibility of the permutation
set.seed(123)

# Setup parallel processing environment using BiocParallel
# Use n-1 cores to leave one free for system processes.
num_cores <- BiocParallel::bpparam()$workers - 1
if (num_cores < 1) num_cores <- 1 # Ensure at least one core is used

cat("--- Starting Permutation Tests (1000 permutations each) ---\n")
cat("Running tests in parallel on", num_cores, "core(s).\n")
csv_report_lines <- c(csv_report_lines, "--- Starting Permutation Tests (1000 permutations each) ---")

# Use bplapply to run the permutation tests in parallel.
# bplapply distributes the tasks across the registered cores.
permutation_results <- bplapply(bed_files_peaks_list, function(public_peaks) {
  # This function is executed for each element in bed_files_peaks_list on a separate core.
  pt <- permTest(A = single_bed_file_peaks, B = public_peaks, ntimes = 1000,
                 randomize.function = randomizeRegions,
                 evaluate.function = numOverlaps,
                 genome = "hg38",
                 verbose = FALSE) # Verbose is FALSE to keep worker logs clean
  return(pt)
})

cat("\n--- All parallel permutation tests completed. Now generating report. ---\n\n")

# --- 3.5. Process and Print Results Sequentially for Reporting ---

# Now that the computation is done, loop through the results to print them
# in a clean, sequential order and capture the output for the CSV file.
for (bed_files_target in names(permutation_results)) {
  # Define separator and header for console and CSV
  separator <- "================================================================="
  test_header <- paste("Performing test for", single_bed_file_name, "vs.", bed_files_target)
  
  # Print to console
  cat("\n", separator, "\n")
  cat(test_header, "\n")
  cat(separator, "\n")
  
  # Capture for CSV report
  csv_report_lines <- c(csv_report_lines, "", separator, test_header, separator)
  
  # Get the stored result object
  pt_result <- permutation_results[[bed_files_target]]
  
  # Print the full results to the console
  print(pt_result)
  
  # Capture the printed output of the permTest result for the CSV report
  pt_output <- capture.output(print(pt_result))
  csv_report_lines <- c(csv_report_lines, pt_output)
}
cat("\n--- All results processed and logged. ---\n\n")


# --- 4. Export Console Output to CSV File ---

cat("--- Exporting summary report to CSV file... ---\n")

# Construct the dynamic filename
target_names_str <- paste(names(bed_files_peaks_list), collapse = "_")
output_filename <- paste0(single_bed_file_name, "_overlap_", target_names_str, ".csv")
output_path <- file.path("../public/Multiomics/Gene Enhancers", output_filename)

# Ensure the target directory exists
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

# Create a data frame with a single column from the captured lines
output_df <- data.frame(Analysis_Log = csv_report_lines)

# Write the data frame to the specified CSV file
write.csv(output_df, file = output_path, row.names = FALSE, quote = TRUE)

cat("Summary report successfully saved to:", output_path, "\n")


################################################################################

################################################################################

# ChIP-seq Analysis Log Parser
# This script reads a specific log file format from a ChIP-seq analysis,
# extracts key information, and generates a summary CSV file.

# --- Configuration ---
# The user should specify the path to their input log file here.
# This variable, `output_path`, is expected to already exist in the user's environment.
# For demonstration, a placeholder is provided.
# Example: output_path <- "path/to/your/analysis_log.csv"
if (!exists("output_path")) {
  stop("Please define the 'output_path' variable with the path to your log file.")
}


# --- 1. Import and Prepare File ---

# The log file is not a standard CSV, so we read it line by line.
# We also remove potential quotation marks from each line for easier parsing.
tryCatch({
  raw_lines <- readLines(output_path)
  lines <- gsub("\"", "", raw_lines)
}, error = function(e) {
  stop(paste("Error reading file:", output_path, "-", e$message))
})


# --- 2. Extract Data and Generate Table ---

# Initialize a list to store data frames for each result row.
results_list <- list()

# --- Extract Primary File Information ---
# Find the line containing the primary file's name and peak count.
primary_line_index <- grep("Number of peaks in .* dataset:", lines)
if (length(primary_line_index) == 0) {
  stop("Could not find primary ChIP-seq peak file information in the log.")
}
primary_line <- lines[primary_line_index[1]]

# Extract the name and peak count using regular expressions.
primary_name <- sub("Number of peaks in (.*) dataset:.*", "\\1", primary_line)
primary_peaks <- as.numeric(sub(".*: (\\d+)", "\\1", primary_line))

# --- Extract Secondary File Peak Counts ---
# Find all lines with "Number of peaks:", then exclude the primary file's line.
# The remaining counts correspond to the secondary files in the order they were imported.
all_peak_indices <- grep("^Number of peaks: \\d+$", lines)
secondary_peak_indices <- setdiff(all_peak_indices, primary_line_index)

if (length(secondary_peak_indices) > 0) {
  secondary_peak_counts <- as.numeric(gsub("Number of peaks: ", "", lines[secondary_peak_indices]))
} else {
  secondary_peak_counts <- numeric(0)
}

# --- Extract Permutation Test Results ---
# Each test result is in a block starting with "Performing test for...".
test_title_indices <- grep("^Performing test for ", lines)

if (length(test_title_indices) == 0) {
  message("No permutation test results were found in the log file.")
} else {
  # Verify that the number of tests matches the number of secondary files found.
  if (length(test_title_indices) != length(secondary_peak_counts)) {
    warning("Mismatch between the number of secondary files and permutation tests found. The output may be incomplete.")
  }
  
  # Loop through each located test block to extract its data.
  for (i in seq_along(test_title_indices)) {
    start_index <- test_title_indices[i]
    
    # Define the end of the current block (it's right before the next test starts, or the end of the file).
    end_index <- if (i < length(test_title_indices)) test_title_indices[i + 1] - 1 else length(lines)
    
    # Create a subset of lines representing the current test block.
    test_block <- lines[start_index:end_index]
    
    # Extract the name of the secondary file from the test title.
    secondary_name <- sub(".* vs. (.*)", "\\1", test_block[1])
    
    # Extract the P-value, Z-score, and number of overlapping peaks from the block.
    p_value_line <- grep("^P-value:", test_block, value = TRUE)
    p_value <- if (length(p_value_line) > 0) as.numeric(sub("P-value: ", "", p_value_line[1])) else NA
    
    z_score_line <- grep("^Z-score:", test_block, value = TRUE)
    z_score <- if (length(z_score_line) > 0) as.numeric(sub("Z-score: ", "", z_score_line[1])) else NA
    
    overlaps_line <- grep("^Evaluation of the original region set:", test_block, value = TRUE)
    num_overlaps <- if (length(overlaps_line) > 0) as.numeric(sub(".*: (\\d+)", "\\1", overlaps_line[1])) else NA
    
    # Retrieve the corresponding secondary peak count (relies on order).
    secondary_peaks <- if (i <= length(secondary_peak_counts)) secondary_peak_counts[i] else NA
    
    # Assemble the extracted data into a one-row data frame.
    # check.names = FALSE allows column names like "P-value".
    results_list[[i]] <- data.frame(
      Primary_BED_File = primary_name,
      Primary_Num_Peaks = primary_peaks,
      Secondary_BED_File = secondary_name,
      Secondary_Num_Peaks = secondary_peaks,
      Num_Overlap_Peaks = num_overlaps,
      `P-value` = p_value,
      `Z-score` = z_score,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
}

# Combine the list of single-row data frames into the final table.
if (length(results_list) > 0) {
  final_table <- do.call(rbind, results_list)
  
  # --- 3. Export the Table ---
  
  # Construct the output filename by appending "_summary" to the original name.
  output_summary_path <- sub("(\\.[^.]+)$", "_summary\\1", output_path)
  if (output_summary_path == output_path) { # Failsafe for names without extensions
    output_summary_path <- paste0(output_path, "_summary.csv")
  }
  
  # Write the final table to a new CSV file.
  write.csv(final_table, file = output_summary_path, row.names = FALSE, quote = FALSE)
  
  message("Processing complete.")
  message(paste("Summary table exported to:", output_summary_path))
  
} else {
  message("No data was processed to generate a summary file.")
}

################################################################################

################################################################################

# R Script for Peak Overlap Visualization
# This script reads a CSV file with peak overlap data, generates a Venn diagram
# for each entry, and creates a summary bar chart of Z-scores.

# --- 1. SETUP: Install and Load Required Packages ---

# Check if VennDiagram package is installed, if not, install it
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
# The VennDiagram package sometimes requires this logger package
if (!requireNamespace("futile.logger", quietly = TRUE)) {
  install.packages("futile.logger")
}

# Check if ggplot2 package is installed, if not, install it
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# Load the necessary libraries
library(VennDiagram)
library(ggplot2)

# --- 2. DATA IMPORT and PATH DEFINITION ---

# Define the input file path
input_file <- output_summary_path

# Define the output directory
output_dir <- "../public/Multiomics/Gene Enhancers"

# Create the output directory if it doesn't already exist to prevent errors
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read the data from the CSV file
# The check.names=FALSE argument is used to handle column names with special characters like '-'
tryCatch({
  overlap_data <- read.csv(input_file, check.names = FALSE)
  # R might replace '-' with '.' in column names, so we'll refer to P-value as `P-value`
  # using backticks or by its converted name P.value. Let's rename for consistency.
  if("P-value" %in% colnames(overlap_data)) {
    names(overlap_data)[names(overlap_data) == 'P-value'] <- 'P_value'
  }
  if("Z-score" %in% colnames(overlap_data)) {
    names(overlap_data)[names(overlap_data) == 'Z-score'] <- 'Z_score'
  }
  
  print("Successfully imported the data file.")
  
}, error = function(e) {
  stop("Error: Could not read the input file. Please ensure the path is correct: ", input_file)
})

# --- 4. CREATE AND EXPORT SUMMARY BAR CHART ---

# Create a bar chart using ggplot2
# X-axis: Secondary_BED_File
# Y-axis: Z-score
z_score_barchart <- ggplot(overlap_data, aes(x = reorder(Secondary_BED_File, -Z_score), y = Z_score)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(
    aes(label = paste("P =", format(P_value, scientific = TRUE, digits = 2))), 
    vjust = -0.5, # Adjust vertical position to be above the bar
    size = 3.5
  ) +
  labs(
    title = "Peak Overlap Significance (Z-score)",
    x = "Dataset",
    y = "Z-score"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12)
  )

# Define the output file path for the bar chart
# It uses the input file name and appends "_bar_Chart.jpg"
barchart_filename <- sub("\\.csv$", "_bar_Chart.jpg", basename(input_file))
barchart_filepath <- file.path(output_dir, barchart_filename)

# Save the bar chart to a JPG file
ggsave(
  filename = barchart_filepath,
  plot = z_score_barchart,
  device = "jpeg",
  width = 12,
  height = 7,
  dpi = 300
)

print(paste("Generated Bar Chart:", barchart_filepath))
print("--- Script finished successfully! ---")

################################################################################

################################################################################

# Gene Enhancer and Binding Pattern Correlation Analysis
# This script reads gene-specific enhancer data and compares it to a set of 
# general binding pattern data files. It calculates the Pearson correlation 
# between the Z-scores of histone mark overlaps for each gene-pattern pair 
# and exports the results to a new summary file for each gene.

# --- 1. SETUP: Define file paths and constants ---

# Set the base paths for input files
gene_enhancers_path <- "../public/Multiomics/Gene Enhancers/"
multiomics_path <- "../public/Multiomics/"

# Define the output directory. The script will create it if it doesn't exist.
output_dir <- "../public/Multiomics/Gene Enhancers"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define the common suffix for the filenames to help with pattern matching and name extraction.
# This part is consistent across all files of interest.
file_suffix <- "_overlap_H3K4me1_H3K4me2_H3K4me3_H3K9ac_H3K9me3_H3K27ac_H3K27me3_H3K36me3_H4K20me1_summary.csv"

# --- 2. GATHER AND LOAD BINDING PATTERN DATA ---

# Define the names of the binding patterns based on the file prefixes.
binding_pattern_names <- c(
  "SRA",
  "PCD-RFP",
  "SRA at Non-UpDEG Enhancers",
  "PCD-RFP at SRA-Specific-UpDEG Enhancers",
  "PCD-RFP at Non-SRA-Specific-UpDEG Enhancers",
  "PCD-RFP at Non-UpDEG Enhancers"
)

# Construct the full paths for the binding pattern files
binding_pattern_filepaths <- file.path(multiomics_path, paste0(binding_pattern_names, file_suffix))

# Read all binding pattern CSV files into a named list of data frames.
# The list is named with the binding pattern for easy access later.
binding_pattern_data_list <- lapply(binding_pattern_filepaths, read.csv)
names(binding_pattern_data_list) <- binding_pattern_names

cat("Successfully loaded", length(binding_pattern_data_list), "binding pattern files.\n")


# --- 3. IDENTIFY AND PROCESS GENE ENHANCER FILES ---

# List all files in the Gene Enhancers directory that match the specific naming convention.
# The pattern looks for any character (`.`) at the start of the name, followed by the suffix.
gene_enhancer_files <- list.files(
  path = gene_enhancers_path,
  pattern = paste0(".+", " Enhancers", file_suffix),
  full.names = TRUE
)

cat("Found", length(gene_enhancer_files), "gene enhancer files to process.\n\n")

# --- 4. MAIN LOOP: Perform Correlation for Each Gene ---

# Loop through each gene enhancer file found.
for (gene_file_path in gene_enhancer_files) {
  
  # Extract the gene symbol from the filename.
  # It removes the directory path and the common suffix parts.
  gene_filename <- basename(gene_file_path)
  gene_symbol <- gsub(paste0(" Enhancers", file_suffix), "", gene_filename)
  
  cat("Processing gene:", gene_symbol, "...\n")
  
  # Read the current gene's data
  gene_data <- read.csv(gene_file_path)
  
  # Create an empty list to store the results for this gene
  results_list <- list()
  
  # --- 5. Inner Loop: Compare gene against each binding pattern ---
  for (bp_name in names(binding_pattern_data_list)) {
    
    # Retrieve the binding pattern data frame from the list
    bp_data <- binding_pattern_data_list[[bp_name]]
    
    # Ensure both data frames have the Z-score column and are of the same length
    if ("Z.score" %in% names(gene_data) && "Z.score" %in% names(bp_data) && nrow(gene_data) == nrow(bp_data)) {
      
      # Perform the Pearson correlation test between the Z-score vectors
      corr_test_result <- cor.test(gene_data$Z.score, bp_data$Z.score, method = "pearson")
      
      # Store the results in a temporary data frame
      temp_result_df <- data.frame(
        Gene = gene_symbol,
        `Binding Pattern` = bp_name,
        `Pearson Correlation` = corr_test_result$estimate,
        `P-value` = corr_test_result$p.value,
        check.names = FALSE # Prevents R from changing column names with spaces
      )
      
      # Add the result to our list
      results_list[[length(results_list) + 1]] <- temp_result_df
      
    } else {
      warning(paste("Skipping correlation for", gene_symbol, "and", bp_name, "due to data mismatch."))
    }
  }
  
  # --- 6. EXPORT RESULTS FOR THE CURRENT GENE ---
  
  # Combine all the results from the list into a single data frame
  final_results_df <- do.call(rbind, results_list)
  
  # Construct the output filename
  output_filename <- file.path(output_dir, paste0(gene_symbol, "_Binding_Pattern_Correlations.csv"))
  
  # Write the final data frame to a new CSV file
  write.csv(final_results_df, file = output_filename, row.names = FALSE)
  
  cat("  -> Saved correlation results to:", output_filename, "\n")
}

cat("\nAnalysis complete. All files have been processed.\n")
