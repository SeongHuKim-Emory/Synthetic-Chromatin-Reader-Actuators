setwd("C:/Users/ski2959/Desktop/NIH R21/Codes/src")

################################################################################

################################################################################

# --- Prerequisites ---

# 1. Check if BiocManager is installed, install if not
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 2. Check if GenomicRanges is installed, install if not
if (!requireNamespace("GenomicRanges", quietly = TRUE))
  BiocManager::install("GenomicRanges")

# 3. Check if rtracklayer is installed, install if not
if (!requireNamespace("rtracklayer", quietly = TRUE))
  BiocManager::install("rtracklayer")

################################################################################

################################################################################

# --- Load Libraries ---
library(GenomicRanges)
library(rtracklayer)
library(methods) # Sometimes needed for GRanges operations

# --- Define File Paths ---
input_file <- "../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks.csv"
output_file <- "../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks_Cytobands.csv"

# --- Step 1: Import Data ---
# Use check.names = FALSE to keep original column names with special characters like '
cat("Reading input file:", input_file, "\n")
tryCatch({
  data <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat("Successfully read", nrow(data), "rows and", ncol(data), "columns.\n")
  # Keep original column names for later reordering
  original_colnames <- names(data)
}, error = function(e) {
  stop("Error reading input file '", input_file, "': ", e$message)
})

# --- Step 2: Prepare Genomic Ranges for TADs ---
cat("Creating GRanges object for TADs...\n")
# Ensure column names match exactly (case-sensitive, including special chars)
if (!all(c("chr", "hg38_TAD_Start", "hg38_TAD_End") %in% names(data))) {
  stop("Input data must contain columns 'chr', 'hg38_TAD_Start', and 'hg38_TAD_End'. Found columns: ", paste(names(data), collapse=", "))
}

# Ensure start coordinates are at least 1 (GRanges is 1-based)
data$hg38_TAD_Start <- pmax(1, data$hg38_TAD_Start)

tryCatch({
  tads_gr <- GRanges(
    seqnames = data$chr,
    ranges = IRanges(start = data$hg38_TAD_Start, end = data$hg38_TAD_End),
    # Add original row index to map results back easily
    original_index = 1:nrow(data) 
  )
  cat("GRanges object created for TADs.\n")
}, error = function(e) {
  stop("Error creating GRanges object from input data: ", e$message)
})

# --- Step 3: Get hg38 Cytoband Data from UCSC ---
cat("Fetching hg38 cytoband data from UCSC...\n")
tryCatch({
  session <- browserSession("UCSC")
  genome(session) <- "hg38"
  
  # --- MODIFIED LINE ---
  # Explicitly query for the table named "cytoBandIdeo" instead of specifying it as a track
  query <- ucscTableQuery(session, table = "cytoBandIdeo") 
  # --- END MODIFIED LINE ---
  
  cytobands_df <- getTable(query)
  cat("Cytoband data fetched successfully.\n")
  
  # Check if data was actually returned
  if (nrow(cytobands_df) == 0) {
    stop("UCSC query returned 0 rows for cytoBandIdeo table.")
  }
  
  # Convert to GRanges object
  # UCSC uses 0-based start coordinates, GRanges uses 1-based.
  cytobands_gr <- GRanges(
    seqnames = cytobands_df$chrom,
    ranges = IRanges(start = cytobands_df$chromStart + 1, end = cytobands_df$chromEnd),
    band = cytobands_df$name,
    stain = cytobands_df$gieStain
  )
  cat("Cytoband data converted to GRanges.\n")
  
  # Ensure seqlevels style matches (e.g., "chr1" vs "1") - tads_gr uses "chr" prefix
  # Make sure tads_gr exists and has seqlevels before trying to modify it
  if (exists("tads_gr") && length(seqlevels(tads_gr)) > 0) {
    seqlevelsStyle(tads_gr) <- "UCSC"
  } else {
    warning("tads_gr object not found or has no seqlevels before attempting style conversion.")
    # If tads_gr wasn't created yet, this step might need adjustment later
  }
  seqlevelsStyle(cytobands_gr) <- "UCSC"
  
}, error = function(e) {
  # Provide a more informative error message including the original error
  stop("Error fetching or processing cytoband data from UCSC: ", e$message)
})



# --- Step 4: Map TADs to Cytobands ---
cat("Mapping TADs to cytobands...\n")

# Find which cytoband overlaps the START of each TAD
tad_starts_gr <- GRanges(seqnames(tads_gr), IRanges(start(tads_gr), start(tads_gr)))
start_overlaps_idx <- findOverlaps(tad_starts_gr, cytobands_gr, select = "first")

# Find which cytoband overlaps the END of each TAD
tad_ends_gr <- GRanges(seqnames(tads_gr), IRanges(end(tads_gr), end(tads_gr)))
end_overlaps_idx <- findOverlaps(tad_ends_gr, cytobands_gr, select = "first")

# Initialize the cytoband result column
cytoband_results <- character(length(tads_gr))

# Process each TAD
for (i in 1:length(tads_gr)) {
  chrom <- as.character(seqnames(tads_gr[i]))
  start_idx <- start_overlaps_idx[i]
  end_idx <- end_overlaps_idx[i]
  
  start_band <- if (!is.na(start_idx)) mcols(cytobands_gr[start_idx])$band else NA
  end_band <- if (!is.na(end_idx)) mcols(cytobands_gr[end_idx])$band else NA
  
  if (is.na(start_band) && is.na(end_band)) {
    # If neither start nor end overlaps (unlikely for valid TADs covering genome regions)
    # Check for any overlap as a fallback
    any_overlap_idx <- findOverlaps(tads_gr[i], cytobands_gr, select = "arbitrary")
    if(!is.na(any_overlap_idx)) {
      bands <- mcols(cytobands_gr[findOverlaps(tads_gr[i], cytobands_gr)])$band
      unique_bands <- unique(bands)
      sorted_bands <- sort(unique_bands)
      if (length(sorted_bands) == 1) {
        cytoband_results[i] <- paste0(chrom, sorted_bands[1])
      } else {
        # Use first and last based on alphabetical/numerical sort as fallback
        cytoband_results[i] <- paste0(chrom, sorted_bands[1], "-", sorted_bands[length(sorted_bands)])
      }
    } else {
      cytoband_results[i] <- NA_character_ # No overlap found
    }
  } else if (is.na(start_band)) {
    # Overlaps end but not start (e.g., TAD starts exactly at band boundary?)
    cytoband_results[i] <- paste0(chrom, end_band)
  } else if (is.na(end_band)) {
    # Overlaps start but not end
    cytoband_results[i] <- paste0(chrom, start_band)
  } else if (start_band == end_band) {
    # TAD is fully contained within a single cytoband
    cytoband_results[i] <- paste0(chrom, start_band)
  } else {
    # TAD spans multiple cytobands
    # Use the band names corresponding to the start and end overlaps
    # Assume the order from UCSC cytoBandIdeo file (and thus start_idx/end_idx) reflects genomic order
    if (start_idx <= end_idx) {
      cytoband_results[i] <- paste0(chrom, start_band, "-", end_band)
    } else {
      # Should be rare (e.g., crossing centromere, or issue in ordering)
      # Follow example format, assuming start_band comes first conceptually if indices are swapped
      cytoband_results[i] <- paste0(chrom, start_band, "-", end_band) 
      # Or potentially swap them if index order implies reversed band order:
      # cytoband_results[i] <- paste0(chrom, end_band, "-", start_band) 
      # Sticking to start-end convention for now.
    }
  }
}
cat("Cytoband mapping complete.\n")

# --- Step 5: Add Cytoband Column to Original Data ---
data$cytoband <- cytoband_results

# --- Step 6: Reorder Columns ---
cat("Reordering columns...\n")
# Define the desired start of the column order
desired_start_order <- c("chr", "TAD_Num", "cytoband", "hg38_TAD_Start", "hg38_TAD_End")

# Get the remaining columns in their original relative order
# Use the saved original_colnames
other_colnames <- original_colnames[!(original_colnames %in% desired_start_order)]

# Combine to get the final order
final_col_order <- c(desired_start_order, other_colnames)

# Check if all expected columns exist in the data frame 'data' before reordering
if (!all(final_col_order %in% names(data))) {
  warning("Mismatch between expected columns and columns present in data after processing. Some columns might be missing or renamed.")
  # Use only columns that are actually present
  final_col_order <- final_col_order[final_col_order %in% names(data)]
  # Ensure all columns from 'data' are included
  missing_from_final <- names(data)[!(names(data) %in% final_col_order)]
  if(length(missing_from_final) > 0) {
    final_col_order <- c(final_col_order, missing_from_final)
    warning("Added columns missing from defined order: ", paste(missing_from_final, collapse=", "))
  }
}

# Apply the new column order
data_final <- data[, final_col_order, drop = FALSE]
cat("Columns reordered.\n")

# --- Step 7: Export Updated Table ---
cat("Writing updated data to:", output_file, "\n")
tryCatch({
  write.csv(data_final, output_file, row.names = FALSE, quote = FALSE, na = "")
  cat("Successfully wrote updated data.\n")
}, error = function(e) {
  stop("Error writing output file '", output_file, "': ", e$message)
})

cat("Script finished successfully.\n")

################################################################################

################################################################################

# R code to import, filter, and export multiomics data

# 1. Import the data file
# Define the input file path
input_file_path <- "../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks_Cytobands.csv"

# Read the CSV file into a data frame
# Added check.names = FALSE to prevent R from altering column names (e.g., prepending 'X' to '021_Peak_Count')
# stringsAsFactors = FALSE is generally recommended practice.
tryCatch({
  multiomics_data <- read.csv(input_file_path, stringsAsFactors = FALSE, check.names = FALSE)
  print(paste("Successfully imported data from:", input_file_path))
  print("Dimensions of imported data:")
  print(dim(multiomics_data))
  print("First few rows of imported data:")
  print(head(multiomics_data))
  print("Column names as R sees them:") # Added for debugging
  print(colnames(multiomics_data))      # Added for debugging
}, error = function(e) {
  stop(paste("Error importing file:", input_file_path, "\nCheck the file path and permissions.\nOriginal error:", e$message))
})


# 2. Filter the data
# Ensure multiomics_data was loaded before proceeding
if (!exists("multiomics_data")) {
  stop("Data import failed. Cannot proceed with filtering.")
}

# Remove rows where the '021_Peak_Count' column is equal to 0
# Check if '021_Peak_Count' column exists
if (!'021_Peak_Count' %in% colnames(multiomics_data)) {
  # If it's still not found, print available columns to help diagnose
  print("Available column names after import:")
  print(colnames(multiomics_data))
  stop("Error: '021_Peak_Count' column not found in the imported data even after check.names=FALSE. Please verify the exact column name in the CSV header and the output above.")
}

# Perform the filtering
# Using $ operator with backticks for column names that might be non-standard,
# though with check.names=FALSE, direct access should work if the name is exact.
# For safety, explicitly referencing the column as a string with `[[...]]` is often more robust.
filtered_data <- multiomics_data[multiomics_data[['021_Peak_Count']] != 0, ]

# Optional: Check the dimensions before and after filtering
print("Dimensions before filtering:")
print(dim(multiomics_data))
print("Dimensions after filtering (removing rows where 021_Peak_Count == 0):")
print(dim(filtered_data))


# 3. Export the updated data frame
# Define the output file path
output_file_path <- "../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks_Cytobands_SRA_Peak_Specific.csv"

# Write the filtered data frame to a new CSV file
# row.names = FALSE prevents R from writing the data frame row numbers as a column
tryCatch({
  write.csv(filtered_data, output_file_path, row.names = FALSE)
  print(paste("Successfully exported filtered data to:", output_file_path))
}, error = function(e) {
  stop(paste("Error exporting file:", output_file_path, "\nCheck the file path and permissions.\nOriginal error:", e$message))
})

# End of script
