# R Script for Finding Overlaps between ChIP-seq Peaks and GeneHancer Elements

# Step 0: Install and load the necessary library
# data.table is used for its efficiency in reading large files (fread)
# and performing fast data manipulation, especially the foverlaps function for genomic intervals.
if (!require("data.table")) install.packages("data.table")
library(data.table)

# --- Step 1: Import Files ---
# Define the paths to the input files.
chipseq_021_path <- "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2 Filtered.csv"
chipseq_025_path <- "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2 Filtered.csv"
genehancer_path <- "../public/Public Dataset/GeneHancer/GeneHancer_AnnotSV_elements_v5.25.txt"

# Use fread for fast and robust file reading.
chipseq_021 <- fread(chipseq_021_path)
chipseq_025 <- fread(chipseq_025_path)
genehancer <- fread(genehancer_path)

# --- Step 5: Create Leniency Variables ---
# These values will be used to extend the start and end coordinates of the ChIP-seq peaks.
Start_Leniency <- 0
End_Leniency <- 0

# --- Pre-processing for Step 6 ---
# To compare the datasets, we must standardize the chromosome names.
# The ChIP-seq files use "chr1", "chr2", etc., while GeneHancer uses "1", "2".
# We will create a new 'chr_key' column in GeneHancer for joining and remove the original
# numeric 'chr' column to prevent a column name collision during the overlap operation.
genehancer[, chr_key := paste0("chr", chr)]
genehancer[, chr := NULL] # Remove the original 'chr' column to avoid conflict.

# Set a key on the GeneHancer data.table. This is required for foverlaps and makes it very fast.
# The key consists of the chromosome and the start/end coordinates of the elements.
setkey(genehancer, chr_key, element_start, element_end)


# --- Step 6: Overlap Finding Function ---
# This function encapsulates the logic for finding overlaps to avoid repeating code.
find_enhancer_overlaps <- function(chipseq_dt, genehancer_dt, start_len, end_len) {
  
  # Select only the necessary columns from the ChIP-seq data.
  # This avoids carrying over extra columns like 'Signal_Value' which aren't needed
  # and prevents potential column name conflicts during the overlap operation, resolving the error.
  peaks <- chipseq_dt[, .(chr, Peak_Name, hg38_Peak_Start, hg38_Peak_End)]
  
  # Apply the leniency values to the peak coordinates to create a search interval.
  peaks[, lenient_start := hg38_Peak_Start - start_len]
  peaks[, lenient_end := hg38_Peak_End + end_len]
  
  # Use foverlaps to find all GeneHancer elements that overlap with the lenient peak intervals.
  # - by.x: Specifies the columns in the 'peaks' table (chr, start, end) for the overlap.
  # - by.y: Specifies the columns in the 'genehancer' table (these are its key columns).
  # - type="any": Finds any overlap, no matter how small.
  # - mult="all": Ensures that if one peak overlaps multiple elements, all are reported (duplicating the peak row).
  # - nomatch=NA: This performs a left-join, keeping all original peaks and filling enhancer columns with NA if no overlap is found.
  overlaps <- foverlaps(
    peaks,
    genehancer_dt,
    by.x = c("chr", "lenient_start", "lenient_end"),
    by.y = c("chr_key", "element_start", "element_end"),
    type = "any",
    mult = "all",
    nomatch = NA
  )
  
  # --- Steps 3 & 4: Select and Format Final Table ---
  # Select the required columns in the specified order to create the final output table.
  # The columns from the original peak file are kept, and the columns from the matched
  # GeneHancer file are appended.
  final_table <- overlaps[, .(
    chr,
    Peak_Name,
    hg38_Peak_Start,
    hg38_Peak_End,
    element_start,
    element_end,
    GHid,
    is_elite,
    regulatory_element_type,
    enhancer_score
  )]
  
  return(final_table)
}


# --- Execute Analysis ---
# Step 2: Create the two new tables by calling the processing function for each ChIP-seq file.
`021_ChIPseq_Enhancer` <- find_enhancer_overlaps(chipseq_021, genehancer, Start_Leniency, End_Leniency)
`025_ChIPseq_Enhancer` <- find_enhancer_overlaps(chipseq_025, genehancer, Start_Leniency, End_Leniency)


# --- Step 7: Export the Resulting Tables ---
# Define output directory and file paths.
output_dir <- "../public/Multiomics"
output_path_021 <- file.path(output_dir, "021_ChIPseq_Enhancer.csv")
output_path_025 <- file.path(output_dir, "025_ChIPseq_Enhancer.csv")

# Create the directory if it doesn't already exist.
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Use fwrite for fast CSV writing. row.names = FALSE is the default and desired behavior.
fwrite(`021_ChIPseq_Enhancer`, output_path_021)
fwrite(`025_ChIPseq_Enhancer`, output_path_025)

# Print a confirmation message to the console.
cat("Analysis complete.\n")
cat("Exported file:", output_path_021, "\n")
cat("Exported file:", output_path_025, "\n")

################################################################################

################################################################################

# Load the dplyr library for data manipulation.
# If you don't have it installed, run: install.packages("dplyr")
library(dplyr)

# 1. Import files
# Define file paths
path_chipseq_021 <- "../public/Multiomics/021_ChIPseq_Enhancer.csv"
path_chipseq_025 <- "../public/Multiomics/025_ChIPseq_Enhancer.csv"
path_genehancer <- "../public/Public Dataset/GeneHancer/GeneHancer_AnnotSV_gene_association_scores_v5.25.txt"

# Read the CSV and TXT files into data frames
chipseq_021 <- read.csv(path_chipseq_021)
chipseq_025 <- read.csv(path_chipseq_025)
genehancer_assoc <- read.delim(path_genehancer, header = TRUE, sep = "\t")

# 2. & 3. Create new tables by copying content
# The data is already loaded into new tables, fulfilling these steps.
# We will work with these data frames directly.
chipseq_021_enhancer_gene <- chipseq_021
chipseq_025_enhancer_gene <- chipseq_025

# 4. Add new columns
# This step is implicitly handled by the join operation later. 
# When we perform a left join, the columns from the right table ('symbol', 'combined_score', etc.)
# are automatically added to the left table. Rows without a match will have NA in these new columns.

# To avoid column name conflicts during the merge, we rename the 'is_elite' column
# in the GeneHancer table to 'is_elite_gene'.
genehancer_assoc <- genehancer_assoc %>%
  rename(is_elite_gene = is_elite)

# 5. Merge the data for both tables
# The left_join function will perform the logic described:
# - It keeps all rows from the ChIP-seq tables.
# - It finds all matching rows in the GeneHancer table based on 'GHid'.
# - It duplicates the ChIP-seq rows to match the number of corresponding GeneHancer entries.
# - It fills in the new columns ('symbol', 'combined_score', 'is_elite_gene') with data from GeneHancer.
# - Rows in the ChIP-seq table with no GHid or no match will be kept, with NA in the new columns.

final_021 <- left_join(chipseq_021_enhancer_gene, genehancer_assoc, by = "GHid")
final_025 <- left_join(chipseq_025_enhancer_gene, genehancer_assoc, by = "GHid")

# 6. Export the final tables to CSV files
# Define output paths
output_path_021 <- "../public/Multiomics/021_ChIPseq_Enhancer_Gene.csv"
output_path_025 <- "../public/Multiomics/025_ChIPseq_Enhancer_Gene.csv"

# Write the data frames to new CSV files.
# row.names = FALSE prevents R from writing row numbers.
# na = "" ensures that missing values are stored as empty fields.
write.csv(final_021, output_path_021, row.names = FALSE, na = "")
write.csv(final_025, output_path_025, row.names = FALSE, na = "")

# 7. Provide confirmation message
print("Processing complete.")
print(paste("File saved to:", output_path_021))
print(paste("File saved to:", output_path_025))

################################################################################

################################################################################

# Load the 'tidyverse' library, which contains helpful packages for data manipulation
# like 'dplyr' for data wrangling and 'readr' for reading/writing files.
# If you don't have it installed, run: install.packages("tidyverse")
library(tidyverse)

# 1. Import the CSV file
# Define the path to your input file.
input_file_path <- "../public/RNAseq/UpDEG_category.csv"

# Read the data from the CSV file into a data frame.
# We use read_csv from the tidyverse, as it has smart defaults.
# col_types = cols(.default = "c") reads all columns as character type initially
# to prevent any data type interpretation issues.
updeg_category <- read_csv(input_file_path, col_types = cols(.default = "c"))

# 2. Create a new table and name it UpDEG_category_merged
# 3. Copy the entire content from the original table
# The new data frame `updeg_merged` is a copy of the original.
updeg_merged <- updeg_category

# 4. Merge the 'Category - SRA specific' and 'Category - other' columns
# The unite() function from the tidyr package is perfect for this task.
# - It creates a new column named "Category".
# - It combines the two specified category columns.
# - 'na.rm = TRUE' ensures that if one of the columns has an NA (missing value),
#   the value from the other column is used.
# - 'sep = ""' specifies that no separator should be placed between the merged values.
# - The original columns are removed by default.
updeg_merged <- updeg_merged %>%
  unite(
    col = "Category",
    c("Category - SRA specific", "Category - other"),
    na.rm = TRUE,
    sep = ""
  )

# 5. For each row, remove the version number from the "transcript_id" column
# We use mutate() combined with str_remove() to modify the 'transcript_id' column.
# The regular expression "\\..*" matches a literal dot (.) followed by any character (.),
# zero or more times (*). This effectively removes the version suffix from each ID.
updeg_merged <- updeg_merged %>%
  mutate(transcript_id = str_remove(transcript_id, "\\..*"))

# 6. If the "Symbol" column value is empty, copy the "transcript_id"
# We use mutate() again with an if_else() condition.
# - It checks if the 'Symbol' column is NA (missing) or an empty string ("").
# - If TRUE, it replaces the value with the content of the 'transcript_id' column for that row.
# - If FALSE, it keeps the existing 'Symbol' value.
updeg_merged <- updeg_merged %>%
  mutate(Symbol = if_else(is.na(Symbol) | Symbol == "", transcript_id, Symbol))

# 7. Export the 'UpDEG_category_merged' table to a new CSV file
# Define the path for the output file.
output_file_path <- "../public/RNAseq/UpDEG_category_merged.csv"

# Write the final, cleaned data frame to a new CSV file.
# - 'na = ""' ensures that any missing values are written as empty fields rather than "NA".
write_csv(updeg_merged, output_file_path, na = "")

# (Optional) Print a confirmation message to the console to indicate completion.
print(paste("Processing complete. The merged file has been saved to:", output_file_path))

# (Optional) Display the first few rows of the final processed data frame.
print("Preview of the final data:")
print(head(updeg_merged))

################################################################################

################################################################################

# R script for importing, merging, and exporting RNA-seq and ChIP-seq data
# This script uses the 'dplyr' package for efficient data manipulation.

# --- Initial Setup ---

# Check if dplyr is installed, and install it if it's not.
# dplyr provides powerful and easy-to-use data manipulation functions.
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

# Load the dplyr library for use in this session.
library(dplyr)

# --- Step 1: Import Files ---

# Print a message to indicate the start of the file import process.
print("Step 1: Importing data files...")

# Define the file paths for clarity and easier modification.
rnaseq_file <- "../public/RNAseq/UpDEG_category_merged.csv"
chipseq_021_file <- "../public/Multiomics/021_ChIPseq_Enhancer_Gene.csv"
chipseq_025_file <- "../public/Multiomics/025_ChIPseq_Enhancer_Gene.csv"

# Read the RNA-seq data containing gene categories.
updeg_data <- read.csv(rnaseq_file)

# Read the two ChIP-seq datasets.
chip_021_data <- read.csv(chipseq_021_file)
chip_025_data <- read.csv(chipseq_025_file)

print("Data import complete.")

# --- Steps 2-5: Prepare and Merge Data ---
# These steps are combined for efficiency. A left_join is used to merge
# data from the RNA-seq file into the ChIP-seq files based on matching
# gene symbols. This single operation creates the new tables, adds the
# required columns, and populates them with data where a match is found.

print("Step 2-5: Merging RNA-seq categories into ChIP-seq data...")

# To avoid adding unnecessary columns, we select only the ones needed for the merge
# from the up-regulated differentially expressed genes (UpDEG) data.
updeg_subset <- select(updeg_data, Symbol, hg38_start, hg38_end, Category)

# Perform a left join on the 021 dataset.
# All rows from `chip_021_data` are kept.
# Where `chip_021_data$symbol` matches `updeg_subset$Symbol`, the corresponding
# hg38_start, hg38_end, and Category values are added.
# If no match is found, these new columns will have NA values.
chip_021_merged <- left_join(chip_021_data, updeg_subset, by = c("symbol" = "Symbol"))

# Perform the same left join on the 025 dataset.
chip_025_merged <- left_join(chip_025_data, updeg_subset, by = c("symbol" = "Symbol"))

print("Data merging complete.")
print("Created tables: 'chip_021_merged' and 'chip_025_merged'")

# --- Step 6 & 7: Export the Merged Tables to CSV ---

print("Step 6 & 7: Exporting merged tables to new CSV files...")

# Define the output file paths.
output_021_file <- "../public/Multiomics/021_ChIPseq_Enhancer_Gene_UpDEG.csv"
output_025_file <- "../public/Multiomics/025_ChIPseq_Enhancer_Gene_UpDEG.csv"

# Write the merged 021 data to a new CSV file.
# `row.names = FALSE` prevents R from writing its internal row numbers.
# `na = ""` ensures that missing values (NAs) are written as empty cells,
# matching the requested output format.
write.csv(chip_021_merged, output_021_file, row.names = FALSE, na = "")

# Write the merged 025 data to its new CSV file with the same settings.
write.csv(chip_025_merged, output_025_file, row.names = FALSE, na = "")

print("Export complete.")
print(paste("File saved to:", output_021_file))
print(paste("File saved to:", output_025_file))

print("Script finished successfully.")

################################################################################

################################################################################

# R script for processing GeneHancer TFBS data

# Install and load the dplyr package if you haven't already.
# dplyr is a powerful package for data manipulation and makes the aggregation step cleaner.
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

# 1. Import the file
# Define the path to the input file.
input_file_path <- "../public/Public Dataset/GeneHancer/GeneHancer_TFBSs_v5.25.txt"

# Read the tab-delimited text file into a data frame.
# 'header = TRUE' specifies that the first row contains column names.
# 'stringsAsFactors = FALSE' ensures that character columns are not converted to factors.
genehancer_data <- read.delim(input_file_path, header = TRUE, stringsAsFactors = FALSE)

# 2. & 3. Create a new table and duplicate the content
# This creates a new data frame to work with, leaving the original data untouched.
GeneHancer_TFBSs_v5.25_MCF7 <- genehancer_data

# 4. Filter for rows where the 'tissues' column contains "MCF-7"
# The grepl() function searches for the specified pattern ("MCF-7") within the 'tissues' column
# and returns a logical vector (TRUE/FALSE), which is then used to select the matching rows.
GeneHancer_TFBSs_v5.25_MCF7 <- GeneHancer_TFBSs_v5.25_MCF7[grepl("MCF-7", GeneHancer_TFBSs_v5.25_MCF7$tissues, fixed = TRUE), ]

# 5. Change the 'tissues' column value for all rows to "MCF7"
# This standardizes the tissue name for the filtered dataset.
GeneHancer_TFBSs_v5.25_MCF7$tissues <- "MCF7"

# 6. Group by 'GHid' and combine 'TF' values
# We use the dplyr package for efficient and readable data aggregation.
GeneHancer_TFBSs_v5.25_MCF7 <- GeneHancer_TFBSs_v5.25_MCF7 %>%
  # Group the data frame by the 'GHid' column.
  group_by(GHid) %>%
  # For each group, summarize the data.
  summarise(
    # Combine all 'TF' values within the group into a single string, separated by a semicolon.
    TF = paste(TF, collapse = ";"),
    # Since all 'tissues' values are now "MCF7", we can just take the first one for the summary row.
    tissues = first(tissues)
  ) %>%
  # It's good practice to ungroup after a summarise operation.
  ungroup() %>%
  # Ensure the final column order matches the desired output.
  select(GHid, TF, tissues)

# 7. Export the resulting table to a CSV file
# Define the name for the output file.
output_file_path <- "../public/Public Dataset/GeneHancer/GeneHancer_TFBSs_v5.25_MCF7.csv"

# Write the data frame to a CSV file.
# 'row.names = FALSE' prevents R from writing its internal row numbers to the file.
write.csv(GeneHancer_TFBSs_v5.25_MCF7, file = output_file_path, row.names = FALSE)

# Optional: Print a message to the console to confirm completion.
cat("Processing complete. The output file has been saved as:", output_file_path, "\n")

################################################################################

################################################################################

# R Script for Merging ChIP-seq Data with Transcription Factor Information

# Suppress warnings for a cleaner output, can be removed for debugging
# options(warn = -1)

# Announce the start of the script
cat("Starting the data processing script...\n")

# --- 1. Import Files ---
# This section reads the three required CSV files into data frames.
# It includes error handling to stop the script if a file cannot be found.
tryCatch({
  cat("Importing GeneHancer TFBS data...\n")
  genehancer_tfbs <- read.csv("../public/Public Dataset/GeneHancer/GeneHancer_TFBSs_v5.25_MCF7.csv")
  
  cat("Importing ChIP-seq dataset 021...\n")
  chipseq_021 <- read.csv("../public/Multiomics/021_ChIPseq_Enhancer_Gene_UpDEG.csv")
  
  cat("Importing ChIP-seq dataset 025...\n")
  chipseq_025 <- read.csv("../public/Multiomics/025_ChIPseq_Enhancer_Gene_UpDEG.csv")
  
  cat("All files imported successfully.\n\n")
}, error = function(e) {
  stop("Fatal Error: Could not read one or more input files. Please check paths and file integrity.\n", e)
})


# --- 2. Create a new table for the 021 dataset ---
# This creates a new data frame by copying the original 021 data.
cat("Creating new table for dataset 021...\n")
chipseq_021_tf <- chipseq_021


# --- 3. Create a new table for the 025 dataset ---
# This creates a new data frame by copying the original 025 data.
cat("Creating new table for dataset 025...\n")
chipseq_025_tf <- chipseq_025


# --- 4. Add a new 'TF' column ---
# A new column named "TF" is appended to both new data frames.
# It is initialized with NA (Not Available) values.
cat("Adding 'TF' column to both new tables...\n")
chipseq_021_tf$TF <- NA
chipseq_025_tf$TF <- NA


# --- 5. Populate the 'TF' Column ---
# This section finds matching 'GHid' values between the ChIP-seq tables
# and the GeneHancer TFBS table. Where a match is found, the corresponding
# transcription factor (TF) is copied over.
# The `match()` function is an efficient way to do this without loops.
cat("Populating 'TF' column for the 021 table based on GHid...\n")
# Find the row indices in `genehancer_tfbs` that correspond to each `GHid` in `chipseq_021_tf`
matched_indices_021 <- match(chipseq_021_tf$GHid, genehancer_tfbs$GHid)
# Use these indices to retrieve the TF values. Non-matches will result in NA.
chipseq_021_tf$TF <- genehancer_tfbs$TF[matched_indices_021]

cat("Populating 'TF' column for the 025 table based on GHid...\n")
# Repeat the same process for the second table.
matched_indices_025 <- match(chipseq_025_tf$GHid, genehancer_tfbs$GHid)
chipseq_025_tf$TF <- genehancer_tfbs$TF[matched_indices_025]
cat("'TF' columns populated successfully.\n\n")


# --- 6. Export the modified tables to new CSV files ---
# The resulting data frames are written to new CSV files in the target directory.
# `row.names = FALSE` prevents R from writing an extra column for row numbers.
# `na = ""` ensures that missing values are written as empty cells instead of "NA".
tryCatch({
  cat("Exporting 021_ChIPseq_Enhancer_Gene_UpDEG_TF.csv...\n")
  write.csv(chipseq_021_tf, 
            file = "../public/Multiomics/021_ChIPseq_Enhancer_Gene_UpDEG_TF.csv", 
            row.names = FALSE, 
            na = "")
  
  cat("Exporting 025_ChIPseq_Enhancer_Gene_UpDEG_TF.csv...\n")
  write.csv(chipseq_025_tf, 
            file = "../public/Multiomics/025_ChIPseq_Enhancer_Gene_UpDEG_TF.csv", 
            row.names = FALSE, 
            na = "")
  
  cat("All files exported successfully.\n")
}, error = function(e) {
  stop("Fatal Error: Could not write one or more output files. Please check directory permissions.\n", e)
})

cat("\nScript finished.\n")

################################################################################

################################################################################

# R Script for Transcription Factor Frequency Analysis

# --- 1. SETUP ---
# Define file paths for input and output
input_file <- "../public/Multiomics/SHK 2025_09_29 025_ChIPseq_Enhancer_Gene_UpDEG_TF by Category.csv"
output_file <- "../public/Multiomics/SHK 2025_09_29 025_ChIPseq_Enhancer_Gene_UpDEG_TF Frequency by Category.csv"

# --- 2. IMPORT DATA ---
# Read the CSV file. 
# 'check.names = FALSE' prevents R from changing column names (e.g., "SRA-specific" to "SRA.specific").
# 'na.strings = ""' ensures that empty cells are treated as missing values (NA).
cat("Reading data from:", input_file, "\n")
tryCatch({
  data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
}, error = function(e) {
  stop("Error reading the input file. Please ensure the path is correct and the file exists.\n", e)
})

# --- 3. PROCESS DATA ---
# Function to count TF occurrences in a given column
get_tf_counts <- function(column_data, count_column_name) {
  # Remove any NA (blank) entries from the column
  tfs_in_column <- column_data[!is.na(column_data)]
  
  # Split the semicolon-delimited strings into a single vector of TFs
  all_tfs <- unlist(strsplit(tfs_in_column, ";"))
  
  # Return NULL if no TFs are found after splitting
  if (length(all_tfs) == 0) {
    return(NULL)
  }
  
  # Calculate the frequency of each unique TF
  tf_freq <- as.data.frame(table(all_tfs))
  
  # Rename columns for clarity
  colnames(tf_freq) <- c("TF", count_column_name)
  
  return(tf_freq)
}

# Apply the function to each of the three specified columns
cat("Calculating TF frequencies for each category...\n")
sra_specific_counts <- get_tf_counts(data$'SRA-specific', "SRA_specific_Count")
non_specific_counts <- get_tf_counts(data$'Non-Specific', "Non_Specific_Count")
non_updeg_counts <- get_tf_counts(data$'Non-UpDEG', "Non_UpDEG_Count")

# Create a list of the resulting data frames for merging
count_list <- list(sra_specific_counts, non_specific_counts, non_updeg_counts)
# Filter out any NULL results (from empty columns)
count_list <- count_list[!sapply(count_list, is.null)]

# --- 4. MERGE RESULTS ---
# Merge all the count data frames together using a full outer join
# This ensures that all TFs from all columns are included in the final result
cat("Merging results...\n")
if (length(count_list) > 0) {
  # Use Reduce to iteratively merge all data frames in the list
  final_counts <- Reduce(function(x, y) merge(x, y, by = "TF", all = TRUE), count_list)
  
  # Replace any NA values (where a TF was not present in a category) with 0
  final_counts[is.na(final_counts)] <- 0
} else {
  # Create an empty data frame with expected columns if no TFs were found at all
  final_counts <- data.frame(TF = character(), SRA_specific_Count = integer(), 
                             Non_Specific_Count = integer(), Non_UpDEG_Count = integer())
  warning("No transcription factors were found in any of the specified columns.")
}


# --- 5. EXPORT RESULTS ---
# Write the final combined data frame to a new CSV file
# 'row.names = FALSE' prevents writing an extra column for row numbers
cat("Exporting results to:", output_file, "\n")
tryCatch({
  write.csv(final_counts, output_file, row.names = FALSE)
  cat("Analysis complete. Results saved successfully.\n")
}, error = function(e) {
  stop("Error writing the output file. Please check permissions.\n", e)
})

################################################################################

################################################################################

# Title: R script to fill missing gene coordinates in ChIP-seq data
# Description: This script reads two CSV files, fetches gene coordinates for hg38 from Ensembl,
#              updates rows with missing coordinates, and saves the results to new CSV files.

# --- 1. SETUP: INSTALL AND LOAD LIBRARIES ---

# Check for BiocManager, install if not present. BiocManager is used to install Bioconductor packages.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required packages if they are not already installed.
# We need 'biomaRt' for accessing Ensembl data, 'dplyr' for data manipulation, and 'readr' for fast CSV reading/writing.
required_packages <- c("biomaRt", "dplyr", "readr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "biomaRt") {
      BiocManager::install("biomaRt")
    } else {
      install.packages(pkg)
    }
  }
}

# Load the libraries into the R session
library(biomaRt)
library(dplyr)
library(readr)

cat("Libraries loaded successfully.\n")

# --- 2. FETCH GENE COORDINATES FROM ENSEMBL (HG38) ---

# This step can take a moment as it connects to the Ensembl database over the internet.
cat("Connecting to Ensembl to fetch hg38 gene coordinates...\n")
# Select the Ensembl 'genes' BioMart database and the human ('hsapiens_gene_ensembl') dataset.
# We use a specific Ensembl version for reproducibility (GRCh38.p13 corresponds to Ensembl version 105).
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 115)

# Retrieve a dataframe of genes with their HGNC symbol, chromosome name, start, and end positions.
gene_coords <- getBM(
  attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
  mart = ensembl
)
cat("Gene coordinates downloaded.\n")

# --- 3. PRE-PROCESS THE ENSEMBL COORDINATES DATA ---

# Clean the downloaded data for easier merging.
gene_coords_processed <- gene_coords %>%
  # Remove rows that don't have a gene symbol, as we can't use them for mapping.
  filter(hgnc_symbol != "") %>%
  # Standardize chromosome names to match the input files (e.g., from '1' to 'chr1').
  mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
  # Handle cases where one gene symbol might map to multiple locations (e.g., different patches).
  # We will group by symbol and chromosome, and take the outermost start and end positions to create one canonical entry.
  group_by(hgnc_symbol, chromosome_name) %>%
  summarise(
    start_position = min(start_position),
    end_position = max(end_position),
    .groups = 'drop' # Drop the grouping structure after summarising.
  )

cat("Gene coordinate data has been processed.\n")

# --- 4. DEFINE A FUNCTION TO PROCESS EACH FILE ---

# This function encapsulates the logic to read, update, and write a single file.
process_chipseq_file <- function(input_path, output_path, gene_data) {
  cat("Processing file:", input_path, "...\n")
  
  # Read the input CSV. We specify column types to prevent parsing errors, especially with empty columns.
  # Start by reading all columns as character to avoid misinterpretation.
  chip_data <- read_csv(input_path, col_types = cols(.default = "c")) %>%
    # Then, convert columns to their correct data types.
    type_convert()
  
  # Left join the ChIP-seq data with our processed gene coordinates data.
  # The join is performed on both 'symbol'/'hgnc_symbol' and 'chr'/'chromosome_name' for accuracy.
  updated_data <- chip_data %>%
    left_join(gene_data, by = c("symbol" = "hgnc_symbol", "chr" = "chromosome_name"))
  
  # Fill in missing coordinates.
  # The coalesce() function is used to take the first non-missing value.
  # If 'hg38_start' has a value, it's kept. If it's missing (NA), the value from 'start_position' (from Ensembl) is used.
  updated_data <- updated_data %>%
    mutate(
      hg38_start = coalesce(hg38_start, start_position),
      hg38_end = coalesce(hg38_end, end_position)
    ) %>%
    # Remove the temporary columns that were added by the join.
    # Explicitly use dplyr::select to avoid conflicts with other packages (like BiocGenerics).
    dplyr::select(-start_position, -end_position)
  
  # Write the fully updated dataframe to the specified output CSV file.
  # Use na = "" to ensure NA values are written as empty cells.
  write_csv(updated_data, output_path, na = "")
  
  cat("Successfully updated and saved to:", output_path, "\n")
}

# --- 5. DEFINE FILE PATHS AND PROCESS FILES ---

# Define the paths for the input and output files.
file1_in <- "../public/Multiomics/021_ChIPseq_Enhancer_Gene_UpDEG_TF.csv"
file2_in <- "../public/Multiomics/025_ChIPseq_Enhancer_Gene_UpDEG_TF.csv"
file1_out <- "../public/Multiomics/021_ChIPseq_Enhancer_Gene_UpDEG_TF_Coor.csv"
file2_out <- "../public/Multiomics/025_ChIPseq_Enhancer_Gene_UpDEG_TF_Coor.csv"

# Call the processing function for each file.
process_chipseq_file(file1_in, file1_out, gene_coords_processed)
process_chipseq_file(file2_in, file2_out, gene_coords_processed)

# --- 6. SCRIPT COMPLETION ---

cat("\nScript finished successfully. All files have been processed.\n")

################################################################################

################################################################################

# R script for annotating genomic regions with histone mark overlap counts
#
# This script will:
# 1. Import enhancer/gene coordinate files and ChIP-seq peak files for various histone marks.
# 2. Create new data frames for annotation.
# 3. Add new columns to store the overlap counts.
# 4. Calculate the number of overlapping histone mark peaks for three types of genomic regions:
#    - "Peak" regions (from hg38_Peak_Start/End)
#    - "Element" regions (from element_start/end)
#    - "UpDEG" regions (from hg38_start/end)
# 5. Export the fully annotated data frames to new CSV files.

# --- 1. Import Files ---

# Install and load the 'readr' package for faster file reading if you don't have it.
# if (!requireNamespace("readr", quietly = TRUE)) {
#   install.packages("readr")
# }
# library(readr)

print("Starting the script...")

# Define file paths for clarity and easy modification
main_files <- c(
  "021" = "../public/Multiomics/021_ChIPseq_Enhancer_Gene_UpDEG_TF_Coor.csv",
  "025" = "../public/Multiomics/025_ChIPseq_Enhancer_Gene_UpDEG_TF_Coor.csv"
)

histone_files <- c(
  "H3K4me1"  = "../public/Public Dataset/ChIPseq/GSE86714 ENCFF991HJA hg38 MCF7 H3K4me1 ChIPseq Filtered.csv",
  "H3K4me2"  = "../public/Public Dataset/ChIPseq/GSE96439 ENCFF188VRU hg38 MCF7 H3K4me2 ChIPseq Filtered.csv",
  "H3K4me3"  = "../public/Public Dataset/ChIPseq/GSE96506 ENCFF268RXB hg38 MCF7 H3K4me3 ChIPseq Filtered.csv",
  "H3K9me3"  = "../public/Public Dataset/ChIPseq/GSE96517 ENCFF501UHK hg38 MCF7 H3K9me3 ChIPseq Filtered.csv",
  "H3K27ac"  = "../public/Public Dataset/ChIPseq/GSE96352 ENCFF491LQY hg38 MCF7 H3K27ac ChIPseq Filtered.csv",
  "H3K27me3" = "../public/Public Dataset/ChIPseq/GSE96363 ENCFF669NUD hg38 MCF7 H3K27me3 ChIPseq Filtered.csv",
  "H3K36me3" = "../public/Public Dataset/ChIPseq/GSE174945 ENCFF195FSD hg38 MCF7 H3K36me3 ChIPseq Filtered.csv"
)

# Load main data files
print("Loading main coordinate files...")
df_021_coor <- read.csv(main_files["021"])
df_025_coor <- read.csv(main_files["025"])

# Load all histone ChIP-seq data into a named list for easy access
print("Loading histone ChIP-seq files...")
histone_data_list <- lapply(histone_files, read.csv)
print("All files loaded successfully.")

# --- 2 & 3. Create New Tables for Histone Annotation ---
df_021_histone <- df_021_coor
df_025_histone <- df_025_coor

# --- 4. Add New Columns ---
# Define the names for the 21 new columns
new_column_names <- c(
  "Peak_H3K4me1", "Peak_H3K4me2", "Peak_H3K4me3", "Peak_H3K9me3", "Peak_H3K27ac", "Peak_H3K27me3", "Peak_H3K36me3",
  "Element_H3K4me1", "Element_H3K4me2", "Element_H3K4me3", "Element_H3K9me3", "Element_H3K27ac", "Element_H3K27me3", "Element_H3K36me3",
  "UpDEG_H3K4me1", "UpDEG_H3K4me2", "UpDEG_H3K4me3", "UpDEG_H3K9me3", "UpDEG_H3K27ac", "UpDEG_H3K27me3", "UpDEG_H3K36me3"
)

# Add the new columns to both data frames and initialize with NA
df_021_histone[new_column_names] <- NA
df_025_histone[new_column_names] <- NA

print("New columns added to the tables.")

# --- 5, 6, & 7. Calculate Overlaps for All Regions and Histone Marks ---

# Helper function to count overlaps for a single genomic region against a histone data frame
count_overlaps <- function(region_chr, region_start, region_end, histone_df) {
  # If coordinates are missing (NA), no overlap can be calculated.
  if (is.na(region_chr) || is.na(region_start) || is.na(region_end)) {
    return(NA_integer_)
  }
  
  # Filter the histone data to only include the relevant chromosome
  histone_subset <- histone_df[histone_df$chr == region_chr & !is.na(histone_df$chr), ]
  
  # If no histone peaks exist on that chromosome, there are no overlaps.
  if (nrow(histone_subset) == 0) {
    return(0)
  }
  
  # The condition for overlap is: (StartA < EndB) AND (EndA > StartB)
  overlap_condition <- (histone_subset$hg38_Peak_Start < region_end) & (histone_subset$hg38_Peak_End > region_start)
  
  # Return the sum of TRUE values, which is the count of overlapping rows.
  return(sum(overlap_condition, na.rm = TRUE))
}

# Main function to process an entire data frame
process_genomic_data <- function(main_df, all_histone_data) {
  # Define the regions to process and their corresponding coordinate columns
  regions_config <- list(
    Peak    = list(start = "hg38_Peak_Start", end = "hg38_Peak_End"),
    Element = list(start = "element_start",   end = "element_end"),
    UpDEG   = list(start = "hg38_start",      end = "hg38_end")
  )
  
  # Iterate over each region type (Peak, Element, UpDEG)
  for (region_name in names(regions_config)) {
    start_col <- regions_config[[region_name]]$start
    end_col <- regions_config[[region_name]]$end
    
    # Iterate over each histone mark (H3K4me1, H3K4me2, etc.)
    for (mark_name in names(all_histone_data)) {
      # Construct the full column name, e.g., "Peak_H3K4me1"
      output_col_name <- paste(region_name, mark_name, sep = "_")
      current_histone_df <- all_histone_data[[mark_name]]
      
      print(paste("   Processing column:", output_col_name))
      
      # Use mapply to apply the count_overlaps function row-wise
      overlap_counts <- mapply(
        count_overlaps,
        main_df$chr,
        main_df[[start_col]],
        main_df[[end_col]],
        MoreArgs = list(histone_df = current_histone_df)
      )
      
      main_df[[output_col_name]] <- overlap_counts
    }
  }
  return(main_df)
}

# Process both data frames
print("--- Starting annotation for 021 data frame ---")
df_021_processed <- process_genomic_data(df_021_histone, histone_data_list)

print("--- Starting annotation for 025 data frame ---")
df_025_processed <- process_genomic_data(df_025_histone, histone_data_list)

# --- 8 & 9. Export the final annotated tables ---
output_files <- c(
  "021" = "../public/Multiomics/021_ChIPseq_Enhancer_Gene_UpDEG_TF_Histone.csv",
  "025" = "../public/Multiomics/025_ChIPseq_Enhancer_Gene_UpDEG_TF_Histone.csv"
)

print("Exporting annotated files...")

# Write the data frames to CSV files. `na = ""` ensures NAs are written as empty cells.
write.csv(df_021_processed, file = output_files["021"], row.names = FALSE, na = "")
write.csv(df_025_processed, file = output_files["025"], row.names = FALSE, na = "")

print("--- Script finished successfully! ---")
print(paste("File 1 exported to:", output_files["021"]))
print(paste("File 2 exported to:", output_files["025"]))
