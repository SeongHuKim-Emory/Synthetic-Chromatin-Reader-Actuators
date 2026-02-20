setwd("C:/Users/ski2959/Desktop/NIH R21/Codes/src")

################################################################################
# Generate the following files
# ../public/RNAseq/SRA_Specific_UpDEG_category.csv
# ../public/RNAseq/SRA_Specific_DownDEG_category.csv
################################################################################

# Load the dplyr package for data manipulation
# If you haven't installed it yet, run: install.packages("dplyr")
library(dplyr)

# Define file paths
up_deg_file_path <- "../public/RNAseq/UpDEG_category.csv"
down_deg_file_path <- "../public/RNAseq/DownDEG_category.csv"

output_up_deg_file_path <- "../public/RNAseq/SRA_Specific_UpDEG_category.csv"
output_down_deg_file_path <- "../public/RNAseq/SRA_Specific_DownDEG_category.csv"

# --- Process UpDEG_category.csv ---
message("Processing UpDEG_category.csv...")

# 1. Import UpDEG_category.csv
# Using check.names = FALSE to keep original column names with spaces
up_deg_data_raw <- read.csv(up_deg_file_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))

# 2. Create a new table and duplicate specified columns
# The column "Category - SRA specific" is the 12th column in the list of columns to keep.
# The original CSV has 13 columns, we are selecting the first 12.
columns_to_select_up <- c(
  "transcript_id", "Symbol", "hg38_chr", "hg38_start", "hg38_end", "hg38_strand",
  "023_24hr_Up", "025_24hr_Up", "021_10hr_Up", "021_24hr_Up", "021_48hr_Up",
  "Category - SRA specific"
)
# Explicitly use dplyr::select and dplyr::all_of
up_deg_table <- up_deg_data_raw %>%
  dplyr::select(dplyr::all_of(columns_to_select_up))

message("Selected columns for UpDEG table. Current dimensions: ", paste(dim(up_deg_table), collapse = "x"))
# print(head(up_deg_table, 3)) # For debugging

# 4. Find rows with empty "Category - SRA specific" column value then delete that entire rows
# We filter for rows where "Category - SRA specific" is NOT NA and not an empty string.
# read.csv with na.strings = c("", "NA") converts empty strings to NA.
# So we only need to check for !is.na()
# Explicitly use dplyr::filter
up_deg_table_filtered <- up_deg_table %>%
  dplyr::filter(!is.na(`Category - SRA specific`))

message("Filtered UpDEG table. Current dimensions: ", paste(dim(up_deg_table_filtered), collapse = "x"))
# print(head(up_deg_table_filtered, 3)) # For debugging

# 5. Remove decimal dot and sub-decimal digits from the "transcript_id" column
# Explicitly use dplyr::mutate
up_deg_table_final <- up_deg_table_filtered %>%
  dplyr::mutate(transcript_id = sub("\\..*", "", transcript_id))

message("Modified transcript_id for UpDEG table.")
# print(head(up_deg_table_final, 3)) # For debugging

# 6. Export the processed table
write.csv(up_deg_table_final, output_up_deg_file_path, row.names = FALSE, na = "")
message("Exported processed UpDEG data to: ", output_up_deg_file_path)


# --- Process DownDEG_category.csv ---
message("\nProcessing DownDEG_category.csv...")

# 1. Import DownDEG_category.csv
down_deg_data_raw <- read.csv(down_deg_file_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))

# 2. Create a new table and duplicate specified columns
columns_to_select_down <- c(
  "transcript_id", "Symbol", "hg38_chr", "hg38_start", "hg38_end", "hg38_strand",
  "023_24hr_Down", "025_24hr_Down", "021_10hr_Down", "021_24hr_Down", "021_48hr_Down",
  "Category - SRA specific"
)
# Explicitly use dplyr::select and dplyr::all_of
down_deg_table <- down_deg_data_raw %>%
  dplyr::select(dplyr::all_of(columns_to_select_down))

message("Selected columns for DownDEG table. Current dimensions: ", paste(dim(down_deg_table), collapse = "x"))
# print(head(down_deg_table, 3)) # For debugging

# 4. Find rows with empty "Category - SRA specific" column value then delete that entire rows
# Explicitly use dplyr::filter
down_deg_table_filtered <- down_deg_table %>%
  dplyr::filter(!is.na(`Category - SRA specific`))

message("Filtered DownDEG table. Current dimensions: ", paste(dim(down_deg_table_filtered), collapse = "x"))
# print(head(down_deg_table_filtered, 3)) # For debugging

# 5. Remove decimal dot and sub-decimal digits from the "transcript_id" column
# Explicitly use dplyr::mutate
down_deg_table_final <- down_deg_table_filtered %>%
  dplyr::mutate(transcript_id = sub("\\..*", "", transcript_id))

message("Modified transcript_id for DownDEG table.")
# print(head(down_deg_table_final, 3)) # For debugging

# 6. Export the processed table
write.csv(down_deg_table_final, output_down_deg_file_path, row.names = FALSE, na = "")
message("Exported processed DownDEG data to: ", output_down_deg_file_path)

message("\nScript finished successfully!")


################################################################################
# Generates following files
# ../public/RNAseq/SRA Specific UpDEG Panther GO Analysis.csv
# ../public/RNAseq/SRA Specific DownDEG Panther GO Analysis.csv
################################################################################

# R Script to Process Panther GO Analysis XML Files and Convert to CSV

# Function to process a single Panther GO XML file and save as CSV
process_panther_xml <- function(input_xml_path, output_csv_path) {
  # Ensure xml2 package is available and loaded
  if (!requireNamespace("xml2", quietly = TRUE)) {
    install.packages("xml2")
  }
  library(xml2)
  
  # Announce processing
  message(paste("Processing XML file:", input_xml_path))
  
  # Read and parse the XML file
  # Adding a tryCatch for robustness in file reading, with fallback for encoding issues
  doc <- tryCatch({
    # Attempt 1: Use read_xml's default behavior (should respect XML declaration or auto-detect)
    message(paste("Attempting to read XML with default/declared encoding:", input_xml_path))
    read_xml(input_xml_path)
  }, error = function(e) {
    warning(paste("Initial XML read attempt failed for", input_xml_path, ":", e$message,
                  "\nThis might be due to an encoding issue or malformed XML. Attempting to read as UTF-8 as a fallback."))
    tryCatch({
      # Attempt 2: Explicitly read file lines as UTF-8 and then parse the resulting string.
      # This can sometimes help if the declared encoding (e.g., UTF-16) is problematic or the file has BOM issues.
      message(paste("Attempting to read XML by forcing UTF-8 encoding for:", input_xml_path))
      xml_content_as_string <- paste(readLines(input_xml_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
      
      # It's good practice to check if any content was read
      if (nchar(xml_content_as_string) == 0) {
        stop("File appears to be empty or could not be read as UTF-8.")
      }
      read_xml(xml_content_as_string)
    }, error = function(e_fallback) {
      # If UTF-8 fallback also fails, then stop with a more comprehensive error message
      stop(paste("Error reading XML file:", input_xml_path,
                 "\nInitial attempt (default/declared encoding) failed with:", e$message,
                 "\nFallback attempt (forcing UTF-8) also failed with:", e_fallback$message,
                 "\nPlease verify the XML file's integrity and encoding."))
    })
  })
  
  # Find all 'result' nodes in the XML
  # The XPath ".//result" searches for 'result' nodes anywhere in the document
  result_nodes <- xml_find_all(doc, ".//result")
  
  # Check if any result nodes were found
  if (length(result_nodes) == 0) {
    message(paste("No <result> nodes found in", input_xml_path, ". An empty CSV will be created if an output path is specified."))
    # Create an empty data frame if no results
    final_df <- data.frame()
  } else {
    # List to store data for each column (term)
    # Each element of the list will be named by the term label and contain a vector of mapped IDs
    columns_data_list <- list()
    
    # Iterate through each 'result' node to extract data
    for (i in seq_along(result_nodes)) {
      node <- result_nodes[[i]]
      
      # Extract term label (this will become a column header)
      # XPath ".//term/label" finds the 'label' node within the 'term' node
      label_node <- xml_find_first(node, ".//term/label")
      
      term_label <- "" # Initialize term_label
      if (!is.na(label_node) && length(xml_text(label_node)) > 0 && nzchar(trimws(xml_text(label_node)))) {
        term_label <- trimws(xml_text(label_node))
      } else {
        # Fallback if label is missing, empty, or only whitespace
        term_label <- paste0("Unnamed_Term_", i)
        warning(paste("Missing or empty label for result node", i, "in", input_xml_path, ". Using fallback name:", term_label))
      }
      
      # Extract mapped IDs for this term
      # XPath ".//input_list/mapped_id_list/mapped_id" finds all 'mapped_id' nodes
      mapped_id_nodes <- xml_find_all(node, ".//input_list/mapped_id_list/mapped_id")
      mapped_ids <- xml_text(mapped_id_nodes)
      
      # Handle potential duplicate term labels by making them unique if necessary
      # This is important because list elements and data frame columns must have unique names
      original_term_label <- term_label
      count <- 1
      while (term_label %in% names(columns_data_list)) {
        term_label <- paste0(original_term_label, "_", count)
        count <- count + 1
        if (count == 2) { # Only warn on the first time it's made unique
          warning(paste0("Duplicate term label '", original_term_label, "' found. Renaming to '", term_label, "'."))
        }
      }
      columns_data_list[[term_label]] <- mapped_ids
    }
    
    # If columns_data_list is populated, proceed to create the data frame
    if (length(columns_data_list) > 0) {
      # Determine the maximum number of rows needed
      # This is the maximum number of mapped IDs found for any single term
      max_rows <- max(sapply(columns_data_list, length))
      
      # If all terms had empty mapped_id_lists, max_rows would be 0 (or -Inf if no terms, handled earlier)
      # Ensure max_rows is at least 0 if there are columns but no data.
      if (is.infinite(max_rows) || max_rows < 0) max_rows <- 0
      
      
      # Pad shorter lists of mapped IDs with NA
      # This makes all columns of equal length, as required for a data frame
      padded_columns_data <- lapply(columns_data_list, function(ids_vector) {
        if (length(ids_vector) < max_rows) {
          c(ids_vector, rep(NA, max_rows - length(ids_vector)))
        } else {
          ids_vector
        }
      })
      
      # Create the data frame from the list of padded columns
      final_df <- as.data.frame(padded_columns_data, stringsAsFactors = FALSE)
      # The names of the list `padded_columns_data` automatically become the column names of `final_df`
    } else {
      # This case might occur if result_nodes were found but all failed to yield labels or data
      final_df <- data.frame() 
      message(paste("No data could be extracted into columns for", input_xml_path))
    }
  }
  
  # Ensure the output directory exists before writing the file
  output_dir <- dirname(output_csv_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    message(paste("Created output directory:", output_dir))
  }
  
  # Write the data frame to a CSV file
  # `row.names = FALSE` prevents writing R's default row numbers
  # `na = ""` represents NA values as empty strings in the CSV output
  tryCatch({
    write.csv(final_df, file = output_csv_path, row.names = FALSE, na = "")
    message(paste("Successfully created CSV:", output_csv_path))
  }, error = function(e) {
    stop(paste("Error writing CSV file:", output_csv_path, "\nOriginal error:", e$message))
  })
}

# --- Main script execution ---

# Define file paths for Up-regulated Differentially Expressed Genes (UpDEG)
updeg_xml_file <- "../public/RNAseq/SRA Specific UpDEG Panther GO Analysis.xml"
updeg_csv_file <- "../public/RNAseq/SRA Specific UpDEG Panther GO Analysis.csv"

# Define file paths for Down-regulated Differentially Expressed Genes (DownDEG)
downdeg_xml_file <- "../public/RNAseq/SRA Specific DownDEG Panther GO Analysis.xml"
downdeg_csv_file <- "../public/RNAseq/SRA Specific DownDEG Panther GO Analysis.csv"

# Process the UpDEG file
process_panther_xml(input_xml_path = updeg_xml_file, output_csv_path = updeg_csv_file)

# Process the DownDEG file
process_panther_xml(input_xml_path = downdeg_xml_file, output_csv_path = downdeg_csv_file)

message("All processing complete.")


################################################################################
# Generates following files
# ../public/RNAseq/SRA Specific UpDEG Panther GO Analysis Merged.csv
# ../public/RNAseq/SRA Specific DownDEG Panther GO Analysis Merged.csv
################################################################################

# R Script for Merging RNAseq Category Data with Panther GO Analysis

# --- 1. Import Files ---

# Define base path for input files
base_path <- "../public/RNAseq/"

# File names
file_up_category <- paste0(base_path, "SRA_Specific_UpDEG_category.csv")
file_down_category <- paste0(base_path, "SRA_Specific_DownDEG_category.csv")
file_up_go <- paste0(base_path, "SRA Specific UpDEG Panther GO Analysis.csv")
file_down_go <- paste0(base_path, "SRA Specific DownDEG Panther GO Analysis.csv")

# Output file names
output_file_up_merged <- paste0(base_path, "SRA Specific UpDEG Panther GO Analysis Merged.csv")
output_file_down_merged <- paste0(base_path, "SRA Specific DownDEG Panther GO Analysis Merged.csv")

# Read the CSV files
# It's good practice to use stringsAsFactors = FALSE, though default in R >= 4.0.0
# The fill=TRUE argument for read.csv is default and handles varying numbers of items per column in GO files by padding with NAs.
message("Importing files...")
tryCatch({
  up_deg_category <- read.csv(file_up_category, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  message(paste("Successfully read:", file_up_category, "- Rows:", nrow(up_deg_category), "Cols:", ncol(up_deg_category)))
  
  down_deg_category <- read.csv(file_down_category, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  message(paste("Successfully read:", file_down_category, "- Rows:", nrow(down_deg_category), "Cols:", ncol(down_deg_category)))
  
  # For GO analysis files, each column is a GO term and rows are transcript IDs.
  # read.csv handles varying numbers of entries by padding with NAs.
  up_go_analysis <- read.csv(file_up_go, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  message(paste("Successfully read:", file_up_go, "- Rows:", nrow(up_go_analysis), "Cols:", ncol(up_go_analysis)))
  
  down_go_analysis <- read.csv(file_down_go, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  message(paste("Successfully read:", file_down_go, "- Rows:", nrow(down_go_analysis), "Cols:", ncol(down_go_analysis)))
  
}, error = function(e) {
  stop(paste("Error reading one or more CSV files. Please check paths and file integrity.\nOriginal error:", e$message))
})


# --- 2. Create new tables (copies of category data) ---
message("\nCreating initial merged tables...")
up_deg_merged <- up_deg_category
down_deg_merged <- down_deg_category

# --- Helper function to process and merge GO terms ---
process_go_terms <- function(merged_df, go_analysis_df, id_column_name = "transcript_id") {
  
  # Get the transcript IDs from the main data frame
  transcript_ids_in_merged_df <- merged_df[[id_column_name]]
  
  # Get the GO term names (column headers from GO analysis file)
  go_term_names <- colnames(go_analysis_df)
  
  if (length(go_term_names) == 0) {
    message("Warning: GO analysis data frame has no columns (no GO terms found).")
    return(merged_df)
  }
  
  message(paste("Processing", length(go_term_names), "GO terms..."))
  
  # Iterate over each GO term
  for (go_term in go_term_names) {
    # Get the list of transcript IDs associated with the current GO term
    # Remove NAs and any empty strings that might have been read
    transcript_ids_for_current_go_term <- go_analysis_df[[go_term]]
    transcript_ids_for_current_go_term <- transcript_ids_for_current_go_term[!is.na(transcript_ids_for_current_go_term) & transcript_ids_for_current_go_term != ""]
    
    # Add a new column for this GO term to the merged_df, initialized to 0
    # Then, mark 1 if the transcript_id from merged_df is found in the list for the current GO term
    merged_df[[go_term]] <- as.integer(transcript_ids_in_merged_df %in% transcript_ids_for_current_go_term)
  }
  
  return(merged_df)
}

# --- 4 & 5. Add columns from GO analysis and populate them ---

# Process Up-regulated DEGs
message("\nProcessing Up-regulated DEGs GO terms...")
if (nrow(up_deg_merged) > 0 && ncol(up_go_analysis) > 0) {
  up_deg_merged <- process_go_terms(up_deg_merged, up_go_analysis, "transcript_id")
  message("Finished processing Up-regulated DEGs.")
  message(paste("UpDEG merged table dimensions - Rows:", nrow(up_deg_merged), "Cols:", ncol(up_deg_merged)))
} else {
  message("Skipping Up-regulated DEGs processing due to empty category data or empty GO analysis data.")
  # If GO analysis was empty, still add empty columns if any were expected
  if (ncol(up_go_analysis) > 0) {
    for (go_term in colnames(up_go_analysis)) {
      if (!go_term %in% colnames(up_deg_merged)) { # Corrected this line
        up_deg_merged[[go_term]] <- 0
      }
    }
  }
}


# Process Down-regulated DEGs
message("\nProcessing Down-regulated DEGs GO terms...")
if (nrow(down_deg_merged) > 0 && ncol(down_go_analysis) > 0) {
  down_deg_merged <- process_go_terms(down_deg_merged, down_go_analysis, "transcript_id")
  message("Finished processing Down-regulated DEGs.")
  message(paste("DownDEG merged table dimensions - Rows:", nrow(down_deg_merged), "Cols:", ncol(down_deg_merged)))
} else {
  message("Skipping Down-regulated DEGs processing due to empty category data or empty GO analysis data.")
  # If GO analysis was empty, still add empty columns if any were expected
  if (ncol(down_go_analysis) > 0) {
    for (go_term in colnames(down_go_analysis)) {
      if (!go_term %in% colnames(down_deg_merged)) { # Corrected this line
        down_deg_merged[[go_term]] <- 0
      }
    }
  }
}

# Display a sample of the updated tables (optional)
message("\nSample of UpDEG merged table (first 6 rows, first 15 columns if available):")
print(head(up_deg_merged[, 1:min(15, ncol(up_deg_merged))]))

message("\nSample of DownDEG merged table (first 6 rows, first 15 columns if available):")
print(head(down_deg_merged[, 1:min(15, ncol(down_deg_merged))]))


# --- 6. Export the tables ---
message("\nExporting merged tables...")
tryCatch({
  write.csv(up_deg_merged, file = output_file_up_merged, row.names = FALSE, quote = TRUE)
  message(paste("Successfully exported UpDEG merged data to:", output_file_up_merged))
  
  write.csv(down_deg_merged, file = output_file_down_merged, row.names = FALSE, quote = TRUE)
  message(paste("Successfully exported DownDEG merged data to:", output_file_down_merged))
  
}, error = function(e) {
  stop(paste("Error exporting one or more CSV files.\nOriginal error:", e$message))
})

message("\nScript finished successfully.")

################################################################################
# Generates following files
# ../public/RNAseq/SRA Specific UpDEG Panther GO Analysis Merged ChIPseq.csv
################################################################################

# Load necessary library
# If you don't have tidyverse installed, uncomment the next line and run it once:
# install.packages("tidyverse")
library(tidyverse) # This loads dplyr, tidyr, readr, etc.

# --- 1. Import files ---
# Define file paths
go_analysis_file <- "../public/RNAseq/SRA Specific UpDEG Panther GO Analysis Merged.csv"
chip_seq_file <- "../public/Multiomics/ChIPseq_025_RNAseq_Overlap_Filtered_Sorted_Condensed.csv"
output_file <- "../public/Multiomics/SRA Specific UpDEG Panther GO Analysis Merged ChIPseq.csv"

# Check if input files exist
if (!file.exists(go_analysis_file)) {
  stop(paste("Error: GO analysis file not found at", go_analysis_file))
}
if (!file.exists(chip_seq_file)) {
  stop(paste("Error: ChIP-seq file not found at", chip_seq_file))
}

# Read the CSV files
# Using read_csv for better type handling and speed
# Adding tryCatch for more robust error handling during file reading
go_analysis_data <- tryCatch({
  readr::read_csv(go_analysis_file, show_col_types = FALSE) # Explicitly use readr::read_csv
}, error = function(e) {
  stop(paste("Error reading GO analysis file:", go_analysis_file, "\nOriginal error:", e$message))
})

chip_seq_data <- tryCatch({
  readr::read_csv(chip_seq_file, show_col_types = FALSE) # Explicitly use readr::read_csv
}, error = function(e) {
  stop(paste("Error reading ChIP-seq file:", chip_seq_file, "\nOriginal error:", e$message))
})

# Display structure of imported data (optional, for verification)
# print("Structure of go_analysis_data:")
# print(str(go_analysis_data))
# print("First few rows of go_analysis_data:")
# print(head(go_analysis_data))

# print("Structure of chip_seq_data:")
# print(str(chip_seq_data))
# print("First few rows of chip_seq_data:")
# print(head(chip_seq_data))

# --- 2. Make a new table and duplicate contents ---
processed_chip_seq <- chip_seq_data

# --- 3. Filter rows: remove if "Category - SRA specific" is empty ---
# Note: read_csv might interpret empty strings as NA by default.
# We filter out rows where `Category - SRA specific` is NA or an empty string.
# Column name in CSV is "Category - SRA specific", R might convert it to "Category..SRA.specific"
# Using backticks if column names have spaces or special characters after import.
# Let's check actual column names after import
# print("Column names in processed_chip_seq before filtering:")
# print(colnames(processed_chip_seq))

# Assuming the column name is "Category - SRA specific" as per the CSV header.
# If read_csv converts it, adjust accordingly. `dplyr` handles spaces well with backticks.
if ("Category - SRA specific" %in% colnames(processed_chip_seq)) {
  processed_chip_seq <- processed_chip_seq %>%
    dplyr::filter(!is.na(`Category - SRA specific`) & `Category - SRA specific` != "") # Explicitly use dplyr::filter
} else if ("Category...SRA.specific" %in% colnames(processed_chip_seq)) {
  # Fallback if R converted the name (less likely with read_csv and backticks)
  processed_chip_seq <- processed_chip_seq %>%
    dplyr::filter(!is.na(Category...SRA.specific) & Category...SRA.specific != "") # Explicitly use dplyr::filter
  # Rename for consistency if it was converted
  # processed_chip_seq <- processed_chip_seq %>% dplyr::rename(`Category - SRA specific` = Category...SRA.specific)
} else {
  warning("Column 'Category - SRA specific' not found for filtering. Please check column names.")
}


# --- 4. Clean "transcript_id": remove decimal and values below ---
# Ensure the column exists before trying to modify it
if ("transcript_id" %in% colnames(processed_chip_seq)) {
  processed_chip_seq <- processed_chip_seq %>%
    dplyr::mutate(transcript_id = sub("\\..*", "", transcript_id)) # Explicitly use dplyr::mutate
} else {
  warning("Column 'transcript_id' not found in processed_chip_seq for cleaning.")
}

# --- 5. Remove the "Category - other" column ---
# Ensure the column exists before trying to remove it
if ("Category - other" %in% colnames(processed_chip_seq)) {
  processed_chip_seq <- processed_chip_seq %>%
    dplyr::select(-`Category - other`) # Explicitly use dplyr::select
} else {
  warning("Column 'Category - other' not found for removal. Please check column names.")
}


# --- 6. Define GO term columns to be added (and potentially merged from go_analysis_data) ---
# These are the columns we expect to get from go_analysis_data
go_term_columns_to_add <- c(
  "regulation.of.membrane.potential", "regulation.of.trans.synaptic.signaling",
  "cell.cell.signaling", "modulation.of.chemical.synaptic.transmission",
  "metal.ion.transport", "system.development", "cell.communication", "signaling",
  "monoatomic.cation.transport", "biological.regulation", "regulation.of.heart.rate",
  "positive.regulation.of.sodium.ion.transport", "regulation.of.biological.process",
  "trans.synaptic.signaling", "cardiac.conduction", "multicellular.organismal.process",
  "monoatomic.ion.transport", "nervous.system.development",
  "multicellular.organism.development", "potassium.ion.transport",
  "regulation.of.sodium.ion.transmembrane.transport",
  "regulation.of.multicellular.organismal.process",
  "inorganic.cation.transmembrane.transport", "synaptic.signaling",
  "regulation.of.developmental.process",
  "regulation.of.heart.rate.by.cardiac.conduction",
  "anterograde.trans.synaptic.signaling", "monoatomic.cation.transmembrane.transport"
)

# Prepare go_analysis_data: select transcript_id and the GO term columns
# Also ensure its transcript_id is clean (though it seems to be already)
if ("transcript_id" %in% colnames(go_analysis_data)) {
  go_analysis_data_subset <- go_analysis_data %>%
    dplyr::mutate(transcript_id = sub("\\..*", "", transcript_id)) %>% # Clean just in case, explicitly use dplyr::mutate
    dplyr::select(transcript_id, dplyr::any_of(go_term_columns_to_add)) # any_of avoids error if a col is missing, explicitly use dplyr::select and dplyr::any_of
} else {
  stop("Column 'transcript_id' not found in go_analysis_data.")
}

# Check which GO term columns are actually present in go_analysis_data_subset
actual_go_columns_in_subset <- colnames(go_analysis_data_subset)[colnames(go_analysis_data_subset) != "transcript_id"]
missing_go_cols <- setdiff(go_term_columns_to_add, actual_go_columns_in_subset)
if (length(missing_go_cols) > 0) {
  warning(paste("The following GO term columns were requested but not found in go_analysis_data:", paste(missing_go_cols, collapse=", ")))
}


# --- 7. Merge GO term data into processed_chip_seq ---
# Perform a left join. Rows in processed_chip_seq will be kept.
# Matching rows from go_analysis_data_subset will add GO term data.
# Non-matching rows will have NA for these new columns.
processed_chip_seq <- processed_chip_seq %>%
  dplyr::left_join(go_analysis_data_subset, by = "transcript_id") # Explicitly use dplyr::left_join

# For any GO term columns that were joined (or should exist), replace NAs with 0
# This applies to columns that were successfully joined from go_analysis_data_subset
# and also ensures any explicitly listed GO columns that might not have been in go_analysis_data
# (and thus wouldn't be created by the join) are added and set to 0 if they don't exist.

# First, ensure all go_term_columns_to_add actually exist in processed_chip_seq after the join,
# adding them with NA if they don't (e.g., if they were missing from go_analysis_data_subset)
for (col_name in go_term_columns_to_add) {
  if (!col_name %in% colnames(processed_chip_seq)) {
    processed_chip_seq[[col_name]] <- NA_integer_ # Or NA_real_ if they can be non-integer
  }
}

# Now, replace NAs with 0 in all specified GO term columns
# Explicitly use dplyr::mutate, dplyr::across, dplyr::any_of, and tidyr::replace_na
processed_chip_seq <- processed_chip_seq %>%
  dplyr::mutate(dplyr::across(dplyr::any_of(go_term_columns_to_add), ~tidyr::replace_na(., 0)))


# --- 8. Export the table ---
# Create directory if it doesn't exist
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

tryCatch({
  readr::write_csv(processed_chip_seq, output_file) # Explicitly use readr::write_csv
  print(paste("Successfully exported the processed table to:", output_file))
}, error = function(e) {
  stop(paste("Error writing output file:", output_file, "\nOriginal error:", e$message))
})

# Display structure and first few rows of the final table (optional)
# print("Structure of the final processed_chip_seq table:")
# print(str(processed_chip_seq))
# print("First few rows of the final processed_chip_seq table:")
# print(head(processed_chip_seq))
# print("Column names in the final table:")
# print(colnames(processed_chip_seq))


################################################################################

################################################################################
