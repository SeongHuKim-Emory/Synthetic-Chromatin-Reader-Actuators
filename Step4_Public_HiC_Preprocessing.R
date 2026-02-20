# Script 1: Combine Boundary Files

# Load necessary libraries
# data.table is used for efficient file reading (fread) and data manipulation
# stringr is used for string manipulation (extracting chromosome name)
library(data.table)
library(stringr)

# --- 1. Define File Paths ---

# Base directory where the files are located
base_dir <- "../public/Public Dataset/HiC/GSE66733_MCF7_HiC"

# List of chromosomes (1 to 22, and X)
chromosomes <- c(1:22, "X")

# Define the filename pattern using sprintf syntax (%s will be replaced by chromosome)
file_pattern <- "HiCStein-MCF7-WT__hg19__chr%s__C-40000-iced.is1000000.ids240000.insulation.boundaries"

# Generate the full list of file paths
file_paths <- file.path(base_dir, sprintf(file_pattern, chromosomes))

# Optional: Check if files exist (good practice)
# files_exist <- file.exists(file_paths)
# if (!all(files_exist)) {
#   warning("Some input files are missing: ", paste(file_paths[!files_exist], collapse=", "))
#   # Optionally stop execution: stop("Cannot proceed without all input files.")
# }

# --- 2. Process Files and Create Combined Table ---

# Create an empty list to store the data tables from each file
all_data_list <- list()

# Loop through each file path
for (file_path in file_paths) {
  cat("Processing file:", basename(file_path), "\n") # Print progress
  
  tryCatch({
    # Read the file using fread for speed and robustness
    # Specify tab separator and indicate there's a header row
    dt <- fread(file_path, sep = "\t", header = TRUE)
    
    # Check if the data table is empty after reading
    if (nrow(dt) == 0) {
      cat("  File is empty or could not be read correctly, skipping.\n")
      next # Skip to the next file
    }
    
    # Check if required columns exist
    required_cols <- c("header", "start", "end", "insulationScore")
    if (!all(required_cols %in% names(dt))) {
      cat("  Required columns missing in file, skipping. Expected:", paste(required_cols, collapse=", "), " Found:", paste(names(dt), collapse=", "), "\n")
      next # Skip to the next file
    }
    
    
    # Create the 'chr' column
    # Extract 'chr' followed by digits or 'X' from the first 'header' column
    # Example: "boundary.0|hg19|chr1:1640001-1840001" -> "chr1"
    # Using str_extract from the stringr package
    dt[, chr := str_extract(header, "chr[0-9X]+")]
    
    # Create the 'hg19_TAD_Boundary' column
    # Calculate the average of the 'start' and 'end' columns
    # Ensure 'start' and 'end' are numeric (fread usually handles this, but explicit conversion adds safety)
    dt[, hg19_TAD_Boundary := (as.numeric(start) + as.numeric(end)) / 2]
    
    # Select and rename the required columns for the final table
    # Note: The column name in the file is 'insulationScore', renaming to 'Insulation_Score'
    dt_processed <- dt[, .(
      chr = chr,
      hg19_TAD_Boundary = hg19_TAD_Boundary,
      Insulation_Score = insulationScore # Keep original name if preferred: Insulation_Score = insulationScore
    )]
    
    # Add the processed data table to the list
    all_data_list[[basename(file_path)]] <- dt_processed
    
  }, error = function(e) {
    # Handle potential errors during file reading or processing
    cat("  Error processing file:", basename(file_path), "-", e$message, "\n")
  })
}

# Combine all the individual data tables in the list into one large data table
# rbindlist is highly efficient for this task
if (length(all_data_list) > 0) {
  combined_boundaries <- rbindlist(all_data_list)
  
  # Optional: Display the structure and first few rows of the combined table
  cat("\n--- Combined Table Summary ---\n")
  print(str(combined_boundaries))
  print(head(combined_boundaries))
  
  # --- 3. Export the Combined Table ---
  
  # Define the output file path with the prefix
  output_file <- file.path(base_dir, "MCF7_hg19_TAD_Boundaries.csv") # MODIFIED FILENAME
  
  # Write the combined data table to a CSV file
  # Using fwrite for efficiency
  # row.names=FALSE: Do not write row numbers to the file
  # quote=FALSE: Do not put quotes around strings (common for genomic coordinate files)
  fwrite(combined_boundaries, output_file, sep = ",", row.names = FALSE, quote = FALSE)
  
  cat("\nSuccessfully combined data from", length(all_data_list), "files.")
  cat("\nOutput saved to:", output_file, "\n")
  
} else {
  cat("\nNo data was processed. No output file created.\n")
}

# --------------------------------------------------------------------------
# Script 2: Liftover hg19 Boundaries to hg38
# --------------------------------------------------------------------------

# --- Configuration ---
# Define file paths relative to the R working directory
# **MODIFIED**: Use the new prefixed filename for input
tad_file_hg19 <- "../public/Public Dataset/HiC/GSE66733_MCF7_HiC/MCF7_hg19_TAD_Boundaries.csv"
chain_file    <- "../public/Public Dataset/HiC/GSE66733_MCF7_HiC/hg19ToHg38.over.chain"
# **MODIFIED**: Use the new prefixed filename for output
output_file   <- "../public/Public Dataset/HiC/GSE66733_MCF7_HiC/MCF7_hg19_hg38_TAD_Boundaries.csv"

# --- Check and Install Required Packages ---
# This script requires the 'rtracklayer' package from Bioconductor
# and 'dplyr' from CRAN.
message("Checking for required packages...")

# Check for BiocManager, needed to install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager...")
  install.packages("BiocManager")
}

# Check for rtracklayer
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  message("Installing rtracklayer from Bioconductor...")
  BiocManager::install("rtracklayer")
}

# Check for GenomicRanges (often installed as a dependency, but good to check)
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  message("Installing GenomicRanges from Bioconductor...")
  BiocManager::install("GenomicRanges")
}

# Check for dplyr
if (!requireNamespace("dplyr", quietly = TRUE)) {
  message("Installing dplyr from CRAN...")
  install.packages("dplyr")
}

# --- Load Libraries ---
message("Loading required libraries...")
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
})
message("Libraries loaded successfully.")

# --- Step 1: Import hg19 TAD Boundaries ---
message("Importing hg19 TAD boundaries from: ", tad_file_hg19)
tryCatch({
  # Read the CSV file
  tad_data_hg19 <- read.csv(tad_file_hg19, stringsAsFactors = FALSE)
  
  # Ensure correct column types (coordinates should be numeric/integer)
  tad_data_hg19$hg19_TAD_Boundary <- as.numeric(tad_data_hg19$hg19_TAD_Boundary)
  tad_data_hg19$Insulation_Score <- as.numeric(tad_data_hg19$Insulation_Score)
  
  # Check if data was loaded
  if (nrow(tad_data_hg19) == 0) {
    stop("Input file is empty or could not be read correctly.")
  }
  message("Successfully imported ", nrow(tad_data_hg19), " boundaries.")
  
}, error = function(e) {
  # Stop execution if file reading fails
  stop("Error reading hg19 TAD boundaries file '", tad_file_hg19, "': ", e$message)
})

# --- Step 2: Prepare Data for Liftover (Create GRanges object) ---
message("Preparing data for liftover by creating a GRanges object...")
# Liftover works with genomic ranges. We treat each boundary coordinate
# as a 1 base-pair interval for the purpose of liftover.
tryCatch({
  gr_hg19 <- GRanges(
    seqnames = Rle(tad_data_hg19$chr), # Chromosome names
    ranges = IRanges(start = tad_data_hg19$hg19_TAD_Boundary, # Start coordinate
                     end = tad_data_hg19$hg19_TAD_Boundary),   # End coordinate (same as start for a single point)
    strand = "*", # Strand is typically not relevant for TAD boundaries
    # Keep the original insulation score as metadata associated with each range
    mcols = tad_data_hg19[, "Insulation_Score", drop = FALSE]
  )
  # Rename the metadata column for clarity within the GRanges object
  colnames(mcols(gr_hg19)) <- "Insulation_Score"
  
  message("GRanges object created successfully.")
  
}, error = function(e) {
  # Stop execution if GRanges creation fails
  stop("Error creating GRanges object: ", e$message)
})

# --- Step 3: Import Chain File ---
message("Importing chain file: ", chain_file)
tryCatch({
  # Import the liftover chain file
  chain <- import.chain(chain_file)
  message("Chain file imported successfully.")
  
}, error = function(e) {
  # Stop execution if chain file import fails
  stop("Error importing chain file '", chain_file, "': ", e$message)
})

# --- Step 4: Perform Liftover ---
message("Performing liftover from hg19 to hg38...")
# The liftOver function takes the hg19 GRanges and the chain file
# It returns a GRangesList, where each element corresponds to an input range.
# An element might contain zero, one, or more ranges if the input maps ambiguously.
lifted_results <- liftOver(gr_hg19, chain)
message("Liftover complete.")

# --- Step 5: Process Results and Extract hg38 Coordinates ---
message("Processing liftover results...")

# Check if the number of results matches the input length
if (length(lifted_results) != length(gr_hg19)) {
  warning("Length of liftover results (", length(lifted_results),
          ") does not match input length (", length(gr_hg19),
          "). This is unexpected and results might be misaligned.")
}

# We need to extract the hg38 coordinate (start position) from the results.
# If a boundary maps to multiple locations in hg38, we'll take the start
# coordinate of the *first* mapped range.
# If a boundary does not map, the corresponding element in GRangesList is empty,
# and we'll assign NA.
# 'vapply' is efficient for this type of extraction.
hg38_coords <- vapply(lifted_results, function(gr) {
  if (length(gr) > 0) {
    # If mapping exists, return the start coordinate of the first range
    return(start(gr)[1])
  } else {
    # If no mapping found (empty GRanges), return NA
    return(NA_integer_) # Use NA_integer_ for consistency with numeric coords
  }
}, FUN.VALUE = integer(1)) # Specify the expected return type for efficiency

# Report how many boundaries failed to lift over
failed_liftover_count <- sum(is.na(hg38_coords))
if (failed_liftover_count > 0) {
  message(failed_liftover_count, " out of ", length(gr_hg19),
          " boundaries could not be lifted over to hg38 and were assigned NA.")
} else {
  message("All boundaries were successfully lifted over (mapped to at least one location in hg38).")
}

# --- Step 6: Add New Column to Original Data ---
message("Adding hg38 coordinates to the table...")
# Create the final data frame by adding the new hg38 coordinates column
tad_data_final <- tad_data_hg19
tad_data_final$hg38_TAD_Boundary <- hg38_coords

# --- Step 7: Reorder Columns ---
message("Reordering columns to the desired format...")
# Use dplyr's select for easy and readable column reordering
tad_data_final <- tad_data_final %>%
  dplyr::select(chr, hg19_TAD_Boundary, hg38_TAD_Boundary, Insulation_Score) # MODIFIED HERE

# Display the first few rows of the resulting table as a check
message("First few rows of the final table:")
print(head(tad_data_final))

# --- Step 8: Export the Table ---
message("Exporting the final table to: ", output_file)
tryCatch({
  # Check if the output directory exists, create it if not
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    message("Output directory '", output_dir, "' does not exist. Creating it...")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Write the data frame to a CSV file
  # row.names = FALSE: Don't write R row numbers to the file
  # quote = FALSE: Don't put quotes around strings (common for genomic coordinate files)
  # na = "NA": Represent missing values (from failed liftover) as "NA" string
  write.csv(tad_data_final, output_file, row.names = FALSE, quote = FALSE, na = "NA")
  message("Table exported successfully.")
  
}, error = function(e) {
  # Report error if writing fails
  warning("Error exporting the table to '", output_file, "': ", e$message)
})

message("Script finished.")

# --------------------------------------------------------------------------
# Script 3: Create hg38 TAD Table from Boundaries
# --------------------------------------------------------------------------

# --- 1. Load Required Libraries ---
# Install packages if you don't have them:
# install.packages("dplyr")
# install.packages("readr")
library(dplyr)
library(readr)

# --- 2. Define Input and Output File Paths ---
# **MODIFIED**: Use the new prefixed filename for input
input_file <- "../public/Public Dataset/HiC/GSE66733_MCF7_HiC/MCF7_hg19_hg38_TAD_Boundaries.csv"
# **MODIFIED**: Use the new prefixed filename for output
output_file <- "../public/Public Dataset/HiC/GSE66733_MCF7_HiC/MCF7_hg38_TADs.csv"

# --- 3. Define Chromosome Lengths (hg38) ---
# Using standard lengths, but chr1 adjusted based on user example/request
chr_lengths <- c(
  "chr1" = 248956423, # Adjusted based on user example (standard is 248956422)
  "chr2" = 242193529,
  "chr3" = 198295559,
  "chr4" = 190214555,
  "chr5" = 181538259,
  "chr6" = 170805979,
  "chr7" = 159345973,
  "chr8" = 145138636,
  "chr9" = 138394717,
  "chr10" = 133797422,
  "chr11" = 135086622,
  "chr12" = 133275309,
  "chr13" = 114364328,
  "chr14" = 107043718,
  "chr15" = 101991189,
  "chr16" = 90338345,
  "chr17" = 83257441,
  "chr18" = 80373285,
  "chr19" = 58617616,
  "chr20" = 64444167,
  "chr21" = 46709983,
  "chr22" = 50818468,
  "chrX" = 156040895,
  "chrY" = 57227415
)

# --- 4. Import TAD Boundary Data ---
message("Reading input file: ", input_file)
# Use tryCatch for basic error handling during file read
boundaries_df <- tryCatch({
  read_csv(input_file, show_col_types = FALSE) # Suppress column type messages
}, error = function(e) {
  stop("Error reading input file: ", input_file, "\nOriginal error: ", e$message)
})

message("Input data dimensions: ", nrow(boundaries_df), " rows, ", ncol(boundaries_df), " columns.")

# --- 5. Filter and Prepare Data ---
message("Filtering and preparing data...")
boundaries_filtered <- boundaries_df %>%
  # Keep only rows with valid hg38 boundaries
  filter(!is.na(hg38_TAD_Boundary)) %>%
  # Keep only chromosomes for which we have lengths defined
  filter(chr %in% names(chr_lengths)) %>%
  # Ensure the boundary coordinate is numeric
  mutate(hg38_TAD_Boundary = as.numeric(hg38_TAD_Boundary)) %>%
  # Arrange by chromosome and then by boundary position
  arrange(chr, hg38_TAD_Boundary)

message("Filtered data dimensions: ", nrow(boundaries_filtered), " rows.")

if (nrow(boundaries_filtered) == 0) {
  stop("No valid hg38 boundaries found in the input file for the specified chromosomes.")
}

# --- 6. Create TAD Table by Chromosome ---
message("Processing boundaries to create TADs for each chromosome...")
tad_list <- list() # Initialize an empty list to store TAD data frames for each chromosome
processed_chromosomes <- unique(boundaries_filtered$chr)

# Loop through each unique chromosome present in the filtered data
for (current_chr in processed_chromosomes) {
  # Subset data for the current chromosome
  chr_data <- boundaries_filtered %>% filter(chr == current_chr)
  
  # Get the length of the current chromosome
  # Use .subset2 for potentially faster access in a loop
  chr_len <- .subset2(chr_lengths, current_chr) # chr_lengths[[current_chr]]
  
  # Extract the boundary coordinates and insulation scores
  boundaries <- chr_data$hg38_TAD_Boundary
  scores <- chr_data$Insulation_Score
  n_boundaries <- length(boundaries)
  
  # Skip if somehow a chromosome has no boundaries after filtering (shouldn't happen with current logic)
  if (n_boundaries == 0) next
  
  # --- Construct TADs for the current chromosome ---
  # Start coordinates: The first TAD starts at 1, subsequent TADs start at the previous boundary.
  starts <- c(1, boundaries)
  
  # End coordinates: TADs end at the current boundary, the last TAD ends at the chromosome length.
  ends <- c(boundaries, chr_len)
  
  # 5' Insulation Scores: The first TAD has 0, subsequent TADs use the score of the previous boundary.
  score_5_prime <- c(0, scores)
  
  # 3' Insulation Scores: TADs use the score of the current boundary, the last TAD has 0.
  score_3_prime <- c(scores, 0)
  
  # Combine into a data frame for this chromosome
  # Use backticks (`) around column names with special characters like '
  # Set check.names = FALSE to prevent R from modifying column names (e.g., changing ' to .)
  chr_tads <- data.frame(
    chr = current_chr,
    hg38_TAD_Start = starts,
    hg38_TAD_End = ends,
    `Insulation_Score_5'` = score_5_prime,
    `Insulation_Score_3'` = score_3_prime,
    check.names = FALSE
  )
  
  # Add the data frame for the current chromosome to the list
  tad_list[[current_chr]] <- chr_tads
}

# --- 7. Combine Results from All Chromosomes ---
message("Combining TAD data from all processed chromosomes...")
final_tads_df <- bind_rows(tad_list)

# Optional: Ensure correct column order just in case bind_rows changed it
final_tads_df <- final_tads_df %>%
  dplyr::select(chr, hg38_TAD_Start, hg38_TAD_End, `Insulation_Score_5'`, `Insulation_Score_3'`) # MODIFIED HERE

# --- 8. Export the Final TAD Table ---
message("Exporting the final TAD table to: ", output_file)
# Optional: Create the output directory if it doesn't exist
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  message("Creating output directory: ", output_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Use tryCatch for basic error handling during file write
tryCatch({
  write_csv(final_tads_df, output_file)
}, error = function(e) {
  stop("Error writing output file: ", output_file, "\nOriginal error: ", e$message)
})

# --- 9. Confirmation Message and Preview ---
message("\nSuccessfully created hg38 TAD table.")
message("Output saved to: ", output_file)
message("\nFirst few rows of the output table:")
print(head(final_tads_df))
message("\nLast few rows of the output table:")
print(tail(final_tads_df))
message("\nDimensions of the output table: ", nrow(final_tads_df), " rows, ", ncol(final_tads_df), " columns.")


# --------------------------------------------------------------------------
# Script 4: Add TAD Number Column
# --------------------------------------------------------------------------

# Load the dplyr library for data manipulation
# If you don't have it installed, run: install.packages("dplyr")
library(dplyr)

# --- 1. Import the CSV file ---
# Define the input file path
# **MODIFIED**: Use the new prefixed filename for input (the one created by Script 3)
input_file_path <- "../public/Public Dataset/HiC/GSE66733_MCF7_HiC/MCF7_hg38_TADs.csv"

# Read the CSV file into a data frame
# 'check.names = FALSE' prevents R from changing column names like 'Insulation_Score_5''
tad_data <- read.csv(input_file_path, stringsAsFactors = FALSE, check.names = FALSE)

# Print the first few rows to verify import (optional)
# print("Original Data Head:")
# print(head(tad_data))

# --- 2. Add the TAD_Num column ---
# Group the data by chromosome ('chr')
# Then, for each group, create a new column 'TAD_Num' with sequential row numbers starting from 1
tad_data_numbered <- tad_data %>%
  group_by(chr) %>%
  mutate(TAD_Num = row_number()) %>%
  ungroup() # Ungrouping is good practice after grouping operations

# --- 3. Reorder columns ---
# Select columns in the desired order
tad_data_final <- tad_data_numbered %>%
  dplyr::select(chr, TAD_Num, everything()) # MODIFIED HERE, 'everything()' includes all other columns

# Print the first few rows of the modified data (optional)
# print("Modified Data Head:")
# print(head(tad_data_final))

# --- 4. Export the modified table ---
# Define the output file path (same as input in this case, overwriting it with the numbered version)
# **MODIFIED**: Use the new prefixed filename for output (consistent with input)
output_file_path <- "../public/Public Dataset/HiC/GSE66733_MCF7_HiC/MCF7_hg38_TADs.csv"

# Write the modified data frame back to the CSV file
# 'row.names = FALSE' prevents writing the default R row numbers to the file
# 'quote = TRUE' is the default for write.csv and generally safe
write.csv(tad_data_final, file = output_file_path, row.names = FALSE, quote = TRUE)

# Print a confirmation message (optional)
print(paste("Successfully processed data and saved to:", output_file_path))
