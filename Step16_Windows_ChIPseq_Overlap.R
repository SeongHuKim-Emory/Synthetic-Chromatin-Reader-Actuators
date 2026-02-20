# R script to import, sort, remove duplicates, and export BED files

# -----------------------------------------------------------------------------
# 1. Define File Paths
# -----------------------------------------------------------------------------
# Create a vector of input file paths.
input_files <- c(
  "../public/Multiomics/025 Peak with SRA UpDEG.bed",
  "../public/Multiomics/025 Peak with Non UpDEG.bed",
  "../public/Multiomics/025 Peak with Non SRA UpDEG.bed",
  "../public/Multiomics/021 Peak with Non UpDEG.bed"
)

# Create a corresponding vector of output file paths by appending "_trim".
output_files <- sub("\\.bed$", "_trim.bed", input_files)


# -----------------------------------------------------------------------------
# 2. Process Each File
# -----------------------------------------------------------------------------
# Loop through the list of files, applying the cleaning steps to each one.
for (i in seq_along(input_files)) {
  
  input_path <- input_files[i]
  output_path <- output_files[i]
  
  # Print a message to the console to track progress.
  cat("Processing:", basename(input_path), "\n")
  
  # Step 1: Import the BED file.
  # read.table is used because BED files are simple, tab-separated text files.
  # 'header = FALSE' specifies that there are no column names in the file.
  # 'sep = "\t"' defines the delimiter as a tab.
  # 'stringsAsFactors = FALSE' is good practice to keep character data as is.
  data <- read.table(input_path, 
                     header = FALSE, 
                     sep = "\t", 
                     stringsAsFactors = FALSE)
  
  # Assign temporary column names for easier handling.
  colnames(data) <- c("V1", "V2", "V3")
  
  # Step 2: Sort the data frame.
  # 'order()' generates the correct row indices for sorting.
  # It first sorts by column V1 (chromosome) and then by V2 (start position)
  # for any ties in V1.
  sorted_data <- data[order(data$V1, data$V2), ]
  
  # Step 3: Remove duplicate rows.
  # 'unique()' returns a new data frame with all identical rows removed.
  unique_data <- unique(sorted_data)
  
  # Report the number of rows before and after processing.
  cat("  - Original row count:", nrow(data), "\n")
  cat("  - Final row count (after removing duplicates):", nrow(unique_data), "\n")
  
  # Step 4: Export the cleaned data to the new file.
  # 'write.table' saves the processed data.
  # - 'quote = FALSE' prevents quotes around the values.
  # - 'sep = "\t"' maintains the tab-separated format.
  # - 'row.names = FALSE' and 'col.names = FALSE' ensure the output is a
  #   standard, headerless BED file.
  write.table(unique_data, 
              file = output_path, 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE)
  
  cat("  - Saved cleaned file to:", basename(output_path), "\n\n")
}

cat("Script finished. All files have been processed.\n")


################################################################################

################################################################################

# R Script for Processing ChIP-seq Data Files

# 1. Define the file paths for the two sets of data.

# First set of files
files_group1 <- c(
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2 Filtered.csv",
  "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2 Filtered.csv"
)

# Second set of files (public datasets)
files_group2 <- c(
  "../public/Public Dataset/ChIPseq/GSE86714 ENCFF991HJA hg38 MCF7 H3K4me1 ChIPseq Filtered.csv",
  "../public/Public Dataset/ChIPseq/GSE96439 ENCFF188VRU hg38 MCF7 H3K4me2 ChIPseq Filtered.csv",
  "../public/Public Dataset/ChIPseq/GSE96506 ENCFF268RXB hg38 MCF7 H3K4me3 ChIPseq Filtered.csv",
  "../public/Public Dataset/ChIPseq/GSE95898 ENCFF348DEB hg38 MCF7 H3K9ac ChIPseq Filtered.csv",
  "../public/Public Dataset/ChIPseq/GSE96517 ENCFF501UHK hg38 MCF7 H3K9me3 ChIPseq Filtered.csv",
  "../public/Public Dataset/ChIPseq/GSE96352 ENCFF491LQY hg38 MCF7 H3K27ac ChIPseq Filtered.csv",
  "../public/Public Dataset/ChIPseq/GSE96363 ENCFF669NUD hg38 MCF7 H3K27me3 ChIPseq Filtered.csv",
  "../public/Public Dataset/ChIPseq/GSE174945 ENCFF195FSD hg38 MCF7 H3K36me3 ChIPseq Filtered.csv",
  "../public/Public Dataset/ChIPseq/GSE96283 ENCFF714DEQ hg38 MCF7 H4K20me1 ChIPseq Filtered.csv"
)

# Define the standard chromosomes to keep
standard_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")

# 2. Process the first group of files
# These files require column selection before filtering and exporting.

cat("Processing Group 1 files...\n")

for (file_path in files_group1) {
  cat("  - Processing:", basename(file_path), "\n")
  
  # Read the CSV file
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Keep only the required columns: "chr", "hg38_Peak_Start", "hg38_Peak_End"
  data_selected <- data[, c("chr", "hg38_Peak_Start", "hg38_Peak_End")]
  
  # Filter to keep only rows with standard chromosomes
  data_filtered <- data_selected[data_selected$chr %in% standard_chromosomes, ]
  
  # Define the output file path by replacing the end of the string
  output_path <- sub("Filtered.csv$", "Filtered_trim.bed", file_path)
  
  # Write the processed data to a tab-delimited file without headers or row names
  write.table(
    data_filtered,
    file = output_path,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  cat("    -> Saved to:", basename(output_path), "\n")
}

cat("Group 1 processing complete.\n\n")

# 3. Process the second group of files (public datasets)
# These files are already in the correct column format.

cat("Processing Group 2 files (Public Datasets)...\n")

for (file_path in files_group2) {
  cat("  - Processing:", basename(file_path), "\n")
  
  # Read the CSV file
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Filter to keep only rows with standard chromosomes
  data_filtered <- data[data$chr %in% standard_chromosomes, ]
  
  # Define the output file path
  output_path <- sub("Filtered.csv$", "Filtered_trim.bed", file_path)
  
  # Write the processed data to a tab-delimited file without headers or row names
  write.table(
    data_filtered,
    file = output_path,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  cat("    -> Saved to:", basename(output_path), "\n")
}

cat("Group 2 processing complete.\n\n")
cat("All tasks finished successfully.\n")


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
output_path <- file.path("../public/Multiomics", output_filename)

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
output_dir <- "../public/Multiomics/"

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


# --- 3. GENERATE VENN DIAGRAMS FOR EACH ROW ---

# Loop through each row of the dataframe
for (i in 1:nrow(overlap_data)) {
  
  # Extract values from the current row
  primary_name <- as.character(overlap_data$Primary_BED_File[i])
  secondary_name <- as.character(overlap_data$Secondary_BED_File[i])
  primary_count <- overlap_data$Primary_Num_Peaks[i]
  secondary_count <- overlap_data$Secondary_Num_Peaks[i]
  overlap_count <- overlap_data$Num_Overlap_Peaks[i]
  
  # Define the output filename for the Venn diagram
  venn_filepath <- file.path(output_dir, paste0(primary_name, "_overlap_", secondary_name, "_Venn_Diagram.jpg"))
  
  # Suppress the log file creation by VennDiagram
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  # Create a new Venn diagram for the current row's data
  venn_plot <- draw.pairwise.venn(
    area1 = primary_count,
    area2 = secondary_count,
    cross.area = overlap_count,
    category = c(primary_name, secondary_name),
    fill = c("skyblue", "lightgreen"),
    alpha = 0.7,
    cex = 1.5,
    cat.cex = 1.5,
    margin = 0.05
  )
  
  # Save the generated Venn diagram to a JPG file
  jpeg(filename = venn_filepath, width = 800, height = 800)
  grid.draw(venn_plot)
  dev.off()
  
  print(paste("Generated Venn Diagram:", venn_filepath))
}


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

# Install and load necessary packages for plotting
# This will install the packages if you don't have them already
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!require("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggplot2)
library(ggrepel)

# 1. Import files

# Define file paths
file_true_binding <- "../public/Multiomics/PCD-RFP chr18_22_X_Y_overlap_H3K4me1_H3K4me2_H3K4me3_H3K9ac_H3K9me3_H3K27ac_H3K27me3_H3K36me3_H4K20me1_summary.csv"
file_pred_binding <- "../public/Multiomics/PCD-RFP Prediction_overlap_H3K4me1_H3K4me2_H3K4me3_H3K9ac_H3K9me3_H3K27ac_H3K27me3_H3K36me3_H4K20me1_summary.csv"

# Read the CSV files into data frames
# The check.names=FALSE argument is used to handle column names with special characters like '-'
data_true <- read.csv(file_true_binding, check.names = FALSE)
data_pred <- read.csv(file_pred_binding, check.names = FALSE)

# 2. Calculate Pearson correlation between the Z-score columns

# Ensure the histone marks are in the same order before comparing Z-scores
# This is a good practice, although in this specific case the files are already ordered.
data_pred <- data_pred[match(data_true$Secondary_BED_File, data_pred$Secondary_BED_File), ]

# Extract Z-score columns
z_scores_true <- data_true$"Z-score"
z_scores_pred <- data_pred$"Z-score"

# Calculate the Pearson correlation coefficient
correlation_result <- cor.test(z_scores_true, z_scores_pred, method = "pearson")
pearson_cor <- correlation_result$estimate

# 3. Make a scatter plot

# Fit a linear model to get the trendline and R-squared value
linear_model <- lm(z_scores_pred ~ z_scores_true)
r_squared <- summary(linear_model)$r.squared

# Create a combined data frame for ggplot
plot_data <- data.frame(
  True_Z_Score = z_scores_true,
  Pred_Z_Score = z_scores_pred,
  Label = data_true$Secondary_BED_File
)

# Create the annotation text for the plot
annotation_text <- paste("Pearson's R =", round(pearson_cor, 3), "\nR² =", round(r_squared, 3))

# Create the scatter plot using ggplot2 and ggrepel for non-overlapping labels
p <- ggplot(plot_data, aes(x = True_Z_Score, y = Pred_Z_Score)) +
  geom_point(color = "blue", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  geom_text_repel(aes(label = Label),
                  box.padding   = 0.4,
                  point.padding = 0.6,
                  segment.color = 'grey50',
                  max.overlaps = Inf) + # Ensure all labels are shown
  labs(
    title = "Correlation of Z-scores: True vs. Predicted Binding",
    x = "True Binding Z-score (chr18_22_X_Y)",
    y = "Predicted Binding Z-score"
  ) +
  theme_bw() + # Use a classic black and white theme
  annotate(
    "text",
    x = -Inf, y = Inf, # Position at the top-left corner
    label = annotation_text,
    hjust = -0.1, vjust = 1.5, # Adjust positioning
    size = 4.5,
    fontface = "bold"
  )

# 4. Export the scatter plot

# Define the output file path
output_file <- "../public/Multiomics/025_Binding_True_vs_Prediction_Correlation.jpg"

# Save the plot using ggsave for better quality control
ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300, units = "in")


# 5. Provide entire code (This script is self-contained)

# Optional: Print the results to the console
print(paste("Pearson Correlation Coefficient:", pearson_cor))
print(paste("R-squared value:", r_squared))
print(paste("Plot saved to:", output_file))

