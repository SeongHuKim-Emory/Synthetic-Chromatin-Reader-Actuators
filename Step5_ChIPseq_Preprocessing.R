setwd("C:/Users/ski2959/Desktop/NIH R21/Codes/src")

################################################################################

################################################################################

# Load necessary libraries
# If you don't have these installed, run: install.packages(c("dplyr", "readr"))
library(dplyr)
library(readr)

# --- Configuration ---

# Define the input broadPeak file path
# Make sure this path is correct relative to where you run the script
# input_file <- "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2 Filtered.broadPeak"
input_file <- "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2.broadPeak"

# Define the output CSV file path
# output_file <- "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2 Filtered.csv"
output_file <- "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2 Filtered.csv"

# --- Processing Steps ---

# 1. Import the broadPeak file
# The file is tab-separated and doesn't have a header row.
# readr::read_tsv automatically handles tab separation.
# col_names = FALSE tells it there's no header.
# show_col_types = FALSE suppresses messages about column types.
cat("Reading input file:", input_file, "\n")
tryCatch({
  peaks_data <- readr::read_tsv(input_file, col_names = FALSE, show_col_types = FALSE)
  
  # Optional: Check the first few rows to confirm import
  # print("Original data head:")
  # print(head(peaks_data))
  
  # 2. Create the new table structure
  cat("Processing data...\n")
  processed_data <- peaks_data %>%
    # Select the required columns using their default names (X1, X2, etc.)
    # and rename them to the desired final names.
    dplyr::select(
      chr = X1,                 # Column 1: Chromosome
      hg38_Peak_Start = X2,     # Column 2: Start position
      hg38_Peak_End = X3,       # Column 3: End position
      Signal_Value = X7         # Column 7: Signal value
    ) %>%
    # Group the data by chromosome. This is crucial for resetting the peak numbering.
    dplyr::group_by(chr) %>%
    # Create the 'Peak_Name' column.
    # paste0() concatenates the chromosome name, an underscore, and the row number within the group.
    # row_number() generates sequential integers starting from 1 for each group (chromosome).
    dplyr::mutate(Peak_Name = paste0(chr, "_", dplyr::row_number())) %>%
    # Ungroup the data frame. It's good practice after grouping operations are done.
    dplyr::ungroup() %>%
    # Reorder the columns to match the desired output format.
    dplyr::select(
      chr,
      Peak_Name,
      hg38_Peak_Start,
      hg38_Peak_End,
      Signal_Value
    )
  
  # Optional: Check the first few rows of the processed data
  # print("Processed data head:")
  # print(head(processed_data))
  # print("Processed data structure:")
  # print(str(processed_data))
  
  
  # 3. Export the table as a CSV file
  # readr::write_csv writes the data frame to a comma-separated file.
  # col_names = TRUE ensures the header row is included in the output file.
  cat("Writing output file:", output_file, "\n")
  readr::write_csv(processed_data, output_file, col_names = TRUE)
  
  cat("Processing complete. Output saved to:", output_file, "\n")
  
}, error = function(e) {
  # Handle potential errors during file reading or processing
  cat("An error occurred:\n")
  print(e$message)
  cat("Please check the input file path and format.\n")
})

################################################################################

################################################################################

# R Script to Process ChIP-seq BED Files

# --- Configuration ---

# Define input file paths
input_files <- c(
  "../public/Public Dataset/ChIPseq/GSE86714 ENCFF991HJA hg38 MCF7 H3K4me1 ChIPseq.bed",
  "../public/Public Dataset/ChIPseq/GSE96439 ENCFF188VRU hg38 MCF7 H3K4me2 ChIPseq.bed",
  "../public/Public Dataset/ChIPseq/GSE96506 ENCFF268RXB hg38 MCF7 H3K4me3 ChIPseq.bed",
  "../public/Public Dataset/ChIPseq/GSE95898 ENCFF348DEB hg38 MCF7 H3K9ac ChIPseq.bed",
  "../public/Public Dataset/ChIPseq/GSE96517 ENCFF501UHK hg38 MCF7 H3K9me3 ChIPseq.bed",
  "../public/Public Dataset/ChIPseq/GSE96352 ENCFF491LQY hg38 MCF7 H3K27ac ChIPseq.bed",
  "../public/Public Dataset/ChIPseq/GSE96363 ENCFF669NUD hg38 MCF7 H3K27me3 ChIPseq.bed",
  "../public/Public Dataset/ChIPseq/GSE174945 ENCFF195FSD hg38 MCF7 H3K36me3 ChIPseq.bed",
  "../public/Public Dataset/ChIPseq/GSE96283 ENCFF714DEQ hg38 MCF7 H4K20me1 ChIPseq.bed"
)

# Define corresponding output file paths
output_files <- c(
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

# --- Processing Function ---

# Function to process a single BED file
process_bed_file <- function(input_path, output_path) {
  # Print status message
  cat("Processing:", basename(input_path), "\n")
  
  # Step 1: Import the BED file
  # BED files are typically tab-separated and don't have headers.
  # We read it into a data frame.
  # Using tryCatch to handle potential file reading errors
  bed_data <- tryCatch({
    read.delim(input_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
  }, error = function(e) {
    cat("Error reading file:", input_path, "\n", conditionMessage(e), "\n")
    return(NULL) # Return NULL if file reading fails
  })
  
  # Check if data reading was successful
  if (is.null(bed_data)) {
    cat("Skipping file due to read error:", basename(input_path), "\n")
    return() # Exit the function for this file
  }
  
  # Check if the data frame has at least 9 columns
  if (ncol(bed_data) < 9) {
    cat("Error: File", basename(input_path), "has fewer than 9 columns. Skipping.\n")
    return()
  }
  
  # Step 2: Filter rows based on the 9th column (V9)
  # Keep rows where the value in the 9th column is >= 1.0
  # Ensure the 9th column is numeric before filtering
  bed_data[, 9] <- as.numeric(bed_data[, 9]) # Convert V9 to numeric
  filtered_data <- bed_data[!is.na(bed_data[, 9]) & bed_data[, 9] >= 1.0, ]
  
  # Check if any data remains after filtering
  if (nrow(filtered_data) == 0) {
    cat("Warning: No rows remaining after filtering for file:", basename(input_path), "\n")
    # Optionally, you might still want to write an empty file or skip writing
    # For now, we'll proceed to create the structure but it will be empty
  }
  
  # Step 3: Create the new table structure
  # Select the first three columns (V1, V2, V3)
  # Rename them as specified
  output_table <- filtered_data[, 1:3]
  colnames(output_table) <- c("chr", "hg38_Peak_Start", "hg38_Peak_End")
  
  # Step 4: Export the new table to a CSV file
  # Use write.csv, ensuring row names are not included in the output
  # Using tryCatch for potential writing errors
  tryCatch({
    write.csv(output_table, file = output_path, row.names = FALSE, quote = FALSE)
    cat("Successfully processed and saved:", basename(output_path), "\n")
  }, error = function(e) {
    cat("Error writing file:", output_path, "\n", conditionMessage(e), "\n")
  })
  
  cat("--------------------\n")
}

# --- Main Execution ---

# Loop through the input files and process each one
# Ensure the number of input and output files match
if (length(input_files) == length(output_files)) {
  for (i in 1:length(input_files)) {
    # Check if the input file exists before processing
    if (file.exists(input_files[i])) {
      process_bed_file(input_files[i], output_files[i])
    } else {
      cat("Error: Input file not found:", input_files[i], "\nSkipping.\n")
      cat("--------------------\n")
    }
  }
} else {
  cat("Error: The number of input files does not match the number of output files.\n")
}

cat("Script finished.\n")
