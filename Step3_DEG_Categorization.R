# Load necessary libraries
# Ensure tidyverse is loaded first to minimize masking issues,
# although explicit calls (dplyr::) are the safest approach.
library(tidyverse) # Includes readr, dplyr, purrr, stringr, etc.
library(biomaRt)   # Added for fetching coordinates

# --- Define File Paths ---

# Input file paths
up_files <- c(
  "../public/RNAseq/MCF7_021_10_UpDEG.csv",
  "../public/RNAseq/MCF7_021_24_UpDEG.csv",
  "../public/RNAseq/MCF7_021_48_UpDEG.csv",
  "../public/RNAseq/MCF7_023_24_UpDEG.csv",
  "../public/RNAseq/MCF7_025_24_UpDEG.csv"
)

down_files <- c(
  "../public/RNAseq/MCF7_021_10_DownDEG.csv",
  "../public/RNAseq/MCF7_021_24_DownDEG.csv",
  "../public/RNAseq/MCF7_021_48_DownDEG.csv",
  "../public/RNAseq/MCF7_023_24_DownDEG.csv",
  "../public/RNAseq/MCF7_025_24_DownDEG.csv"
)

# Output file paths
output_up_file <- "../public/RNAseq/UpDEG_category.csv" # Updated name for clarity
output_down_file <- "../public/RNAseq/DownDEG_category.csv" # Updated name for clarity

# --- Function to generate column names from file paths ---
# Example input: "../public/RNAseq/MCF7_021_10_UpDEG.csv"
# Example output: "021_10hr_Up"
generate_col_name <- function(filepath) {
  filename <- basename(filepath)
  # Extract parts: MCF7_XXX_YY_DirectionDEG.csv
  # Use tryCatch for robustness in case pattern doesn't match
  parts <- tryCatch({
    stringr::str_match(filename, "MCF7_(\\d+)_(\\d+)_(\\w+)DEG\\.csv")
  }, error = function(e) {
    warning("Could not parse filename pattern: ", filename, " - Error: ", e$message)
    return(matrix(NA, nrow = 1, ncol = 4)) # Return NA matrix on error
  })
  
  # Check if matching worked
  if (any(is.na(parts[1, 2:4]))) {
    # Try matching if "DEG" is missing before .csv (more robust)
    parts <- tryCatch({
      stringr::str_match(filename, "MCF7_(\\d+)_(\\d+)_(\\w+)\\.csv")
    }, error = function(e) {
      warning("Could not parse filename pattern (alternative): ", filename, " - Error: ", e$message)
      return(matrix(NA, nrow = 1, ncol = 4))
    })
    if (any(is.na(parts[1, 2:4]))) {
      stop("Filename does not match expected pattern: ", filename)
    }
  }
  
  condition_code <- parts[, 2]
  time_point <- parts[, 3]
  direction <- parts[, 4] # Up or Down
  
  # Format time
  time_hr <- paste0(time_point, "hr")
  
  # Combine parts
  col_name <- paste(condition_code, time_hr, direction, sep = "_")
  return(col_name)
}


# --- Generic Function to Read and Process DEG Files ---
process_deg_files <- function(file_list, direction_suffix) {
  cat(sprintf("Processing %s-regulated files...\n", direction_suffix))
  
  data_list <- purrr::map(file_list, ~{
    # Check if file exists before trying to read
    if (!file.exists(.x)) {
      warning("File not found: ", .x, ". Skipping.")
      return(NULL) # Return NULL if file doesn't exist
    }
    col_name <- generate_col_name(.x)
    cat("Reading:", .x, " -> Column:", col_name, "\n")
    
    # Use tryCatch for reading robustness
    df <- tryCatch({
      readr::read_csv(.x, show_col_types = FALSE)
    }, error = function(e) {
      warning("Error reading file: ", .x, " - Error: ", e$message)
      return(NULL) # Return NULL on read error
    })
    
    # Proceed only if reading was successful and df is not NULL
    if (!is.null(df) && nrow(df) > 0) {
      # Ensure required columns exist
      if (!all(c("transcript_id", "symbol") %in% names(df))) {
        warning("File ", .x, " is missing required columns 'transcript_id' or 'symbol'. Skipping.")
        return(NULL)
      }
      df %>%
        # Use dplyr::select explicitly to avoid masking conflicts
        dplyr::select(transcript_id, symbol) %>%
        # Use dplyr::mutate explicitly
        dplyr::mutate(!!col_name := 1) # Create indicator column with value 1
    } else if (!is.null(df) && nrow(df) == 0) {
      warning("File is empty: ", .x)
      # Return an empty tibble with expected columns but specific indicator
      tibble(transcript_id = character(),
             symbol = character(),
             !!col_name := numeric()) # Use numeric() for type consistency
    } else {
      # Handle case where read_csv failed and returned NULL
      NULL
    }
  })
  
  # Filter out NULL elements (from file not found, read errors, or missing columns) before reducing
  data_list <- purrr::compact(data_list)
  
  # Check if list is empty after filtering
  if (length(data_list) == 0) {
    stop(sprintf("No valid %sDEG files were read successfully. Cannot proceed.", direction_suffix))
  }
  
  # Define expected columns based on input files
  expected_cols <- purrr::map_chr(file_list, generate_col_name)
  
  # Combine all data frames using full join, replacing NAs with 0s
  combined_data <- purrr::reduce(data_list,
                                 dplyr::full_join,
                                 by = c("transcript_id", "symbol")) %>%
    # Use dplyr::mutate and dplyr::across explicitly
    dplyr::mutate(dplyr::across(where(is.numeric), ~replace_na(., 0)))
  
  # Ensure all expected columns exist, adding them with 0 if missing
  # This handles cases where a file was empty or skipped entirely
  missing_cols <- setdiff(expected_cols, names(combined_data))
  if(length(missing_cols) > 0) {
    for(col in missing_cols) {
      combined_data <- combined_data %>% dplyr::mutate(!!col := 0)
    }
    warning("The following expected columns were missing and added with 0s: ", paste(missing_cols, collapse=", "))
  }
  
  
  return(combined_data)
}

# --- Process Up-regulated Genes ---
combined_up_data <- process_deg_files(up_files, "Up")

# --- Process Down-regulated Genes ---
combined_down_data <- process_deg_files(down_files, "Down")


# --- Fetch Coordinates using biomaRt ---
cat("\nFetching coordinates from Ensembl (biomaRt)...\n")

# Get unique IDs from both datasets
all_transcript_ids <- unique(c(combined_up_data$transcript_id, combined_down_data$transcript_id))

# Extract gene ID part (assuming format ENSGXXXX.version)
# If your IDs are different (e.g., already just ENSGXXX), adjust this
gene_ids_no_version <- unique(sub("\\.\\d+$", "", all_transcript_ids))

coords_processed <- NULL # Initialize

if (length(gene_ids_no_version) > 0) {
  # Set up biomaRt connection to Ensembl GRCh38 (hg38)
  ensembl <- NULL
  tryCatch({
    # Using https is generally preferred
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
    # For GRCh37/hg19 use: host = "https://grch37.ensembl.org"
  }, error = function(e) {
    warning("Could not connect to Ensembl biomart (https://www.ensembl.org): ", e$message,
            "\nAttempting http connection...")
    # Try http as a fallback
    tryCatch({
      ensembl <<- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://www.ensembl.org")
    }, error = function(e_http) {
      warning("Could not connect to Ensembl biomart (http://www.ensembl.org either): ", e_http$message)
    })
  })
  
  
  if (!is.null(ensembl)) {
    cat("Connected to Ensembl. Querying for coordinates...\n")
    coords <- NULL
    tryCatch({
      coords <- getBM(attributes = c('ensembl_gene_id',    # The ID we are filtering by
                                     'chromosome_name',    # e.g., 1, 2, X
                                     'start_position',     # Gene start
                                     'end_position',       # Gene end
                                     'strand'),            # 1 or -1
                      filters = 'ensembl_gene_id',
                      values = gene_ids_no_version,
                      mart = ensembl)
    }, error = function(e) {
      warning("biomaRt query failed: ", e$message)
      coords <<- NULL # Ensure coords is NULL on error
    })
    
    if (!is.null(coords) && nrow(coords) > 0) {
      cat("Processing retrieved coordinates...\n")
      # Process coordinates: rename, translate strand, handle potential duplicates
      coords_processed <- coords %>%
        dplyr::as_tibble() %>% # Ensure it's a tibble
        dplyr::mutate(
          hg38_strand = ifelse(strand == 1, "forward strand", "reverse strand"),
          # Prepend "chr" if desired, although Ensembl usually doesn't use it for primary chromosomes
          # hg38_chr = ifelse(grepl("^[0-9XYM]+$", chromosome_name), paste0("chr", chromosome_name), chromosome_name)
          hg38_chr = chromosome_name # Keep original chromosome name for now
        ) %>%
        dplyr::select(
          gene_id_no_version = ensembl_gene_id, # Keep for joining
          hg38_chr,
          hg38_start = start_position,
          hg38_end = end_position,
          hg38_strand
        ) %>%
        # Handle cases where one gene ID might return multiple locations (e.g., patches)
        # Here, we just keep the first one encountered for simplicity.
        dplyr::distinct(gene_id_no_version, .keep_all = TRUE)
      
      # Report how many IDs were successfully mapped
      found_ids <- length(intersect(gene_ids_no_version, coords_processed$gene_id_no_version))
      cat(sprintf("Successfully retrieved coordinates for %d out of %d unique gene IDs.\n",
                  found_ids, length(gene_ids_no_version)))
      if (found_ids < length(gene_ids_no_version)) {
        missing_ids_sample <- setdiff(gene_ids_no_version, coords_processed$gene_id_no_version)
        cat(sprintf("Example missing IDs: %s\n",
                    paste(head(missing_ids_sample, 5), collapse=", ")))
      }
      
    } else if (!is.null(coords) && nrow(coords) == 0) {
      warning("biomaRt query returned 0 results for the provided gene IDs.")
    } else {
      # Error occurred during getBM or connection was NULL
      warning("Coordinate retrieval failed. Coordinate columns will be NA.")
    }
  } else {
    warning("Could not establish connection to biomaRt. Coordinate columns will be NA.")
  }
} else {
  warning("No valid gene IDs found to query for coordinates.")
}


# --- Prepare data for Joining Coordinates ---

# Add the gene ID without version number for joining
# Also ensure symbol is character to avoid potential issues later
combined_up_data <- combined_up_data %>%
  dplyr::mutate(gene_id_no_version = sub("\\.\\d+$", "", transcript_id),
                symbol = as.character(symbol)) # Ensure symbol is character

combined_down_data <- combined_down_data %>%
  dplyr::mutate(gene_id_no_version = sub("\\.\\d+$", "", transcript_id),
                symbol = as.character(symbol)) # Ensure symbol is character


# --- Join Coordinates if Available ---
if (!is.null(coords_processed) && nrow(coords_processed) > 0) {
  cat("Joining coordinates to UpDEG data...\n")
  combined_up_data <- combined_up_data %>%
    dplyr::left_join(coords_processed, by = "gene_id_no_version")
  
  cat("Joining coordinates to DownDEG data...\n")
  combined_down_data <- combined_down_data %>%
    dplyr::left_join(coords_processed, by = "gene_id_no_version")
} else {
  cat("Skipping coordinate joining as no coordinates were retrieved.\n")
  # Add empty coordinate columns if they don't exist to maintain structure
  coord_cols_to_add <- c("hg38_chr", "hg38_start", "hg38_end", "hg38_strand")
  for (col in coord_cols_to_add) {
    if (!(col %in% names(combined_up_data))) {
      combined_up_data <- combined_up_data %>% dplyr::mutate(!!col := NA)
    }
    if (!(col %in% names(combined_down_data))) {
      combined_down_data <- combined_down_data %>% dplyr::mutate(!!col := NA)
    }
  }
  # Ensure correct types if added as NA logical
  combined_up_data <- combined_up_data %>%
    dplyr::mutate(hg38_chr = as.character(hg38_chr),
                  hg38_start = as.integer(hg38_start),
                  hg38_end = as.integer(hg38_end),
                  hg38_strand = as.character(hg38_strand))
  combined_down_data <- combined_down_data %>%
    dplyr::mutate(hg38_chr = as.character(hg38_chr),
                  hg38_start = as.integer(hg38_start),
                  hg38_end = as.integer(hg38_end),
                  hg38_strand = as.character(hg38_strand))
}


# --- Apply Categorization Rules (UpDEGs) ---
cat("Applying categorization rules to UpDEGs...\n")
categorized_up_data <- combined_up_data %>%
  dplyr::mutate(
    # Initialize category columns
    `Category - SRA specific` = "",
    `Category - other` = "",
    
    # Apply rules using case_when (prioritized order)
    # Ensure column names used here EXACTLY match those generated earlier
    `Category - other` = dplyr::case_when(
      `023_24hr_Up` == 1 ~ "ns",
      # Check for 025=1 AND any 021 is non-zero (sum > 0)
      `025_24hr_Up` == 1 & (`021_10hr_Up` + `021_24hr_Up` + `021_48hr_Up` > 0) ~ "PCD-RFP & SRA",
      # Check for 025=1 AND all 021 are zero (sum = 0)
      `025_24hr_Up` == 1 & (`021_10hr_Up` + `021_24hr_Up` + `021_48hr_Up` == 0) ~ "PCD-RFP",
      TRUE ~ `Category - other` # Keep existing value (or default "") if no 023/025 match
    ),
    
    # Apply SRA specific rules ONLY if 023 and 025 conditions didn't apply
    `Category - SRA specific` = dplyr::case_when(
      `023_24hr_Up` == 0 & `025_24hr_Up` == 0 & `021_10hr_Up` == 1 & `021_24hr_Up` == 0 & `021_48hr_Up` == 0 ~ "Early transient",
      `023_24hr_Up` == 0 & `025_24hr_Up` == 0 & `021_10hr_Up` == 1 & `021_24hr_Up` == 1 & `021_48hr_Up` == 0 ~ "Early transient",
      `023_24hr_Up` == 0 & `025_24hr_Up` == 0 & `021_10hr_Up` == 1 & `021_24hr_Up` == 1 & `021_48hr_Up` == 1 ~ "Early sustained",
      `023_24hr_Up` == 0 & `025_24hr_Up` == 0 & `021_10hr_Up` == 1 & `021_24hr_Up` == 0 & `021_48hr_Up` == 1 ~ "Stochastic",
      `023_24hr_Up` == 0 & `025_24hr_Up` == 0 & `021_10hr_Up` == 0 & `021_24hr_Up` == 1 & `021_48hr_Up` == 1 ~ "Late sustained",
      `023_24hr_Up` == 0 & `025_24hr_Up` == 0 & `021_10hr_Up` == 0 & `021_24hr_Up` == 1 & `021_48hr_Up` == 0 ~ "Late transient",
      `023_24hr_Up` == 0 & `025_24hr_Up` == 0 & `021_10hr_Up` == 0 & `021_24hr_Up` == 0 & `021_48hr_Up` == 1 ~ "Delayed",
      TRUE ~ `Category - SRA specific` # Keep existing value (or default "") otherwise
    )
  ) %>%
  # Select and reorder columns for the final output (use dplyr::select)
  dplyr::select(
    transcript_id,
    Symbol = symbol, # Rename symbol column
    hg38_chr,        # Added
    hg38_start,      # Added
    hg38_end,        # Added
    hg38_strand,     # Added
    `023_24hr_Up`,
    `025_24hr_Up`,
    `021_10hr_Up`,
    `021_24hr_Up`,
    `021_48hr_Up`,
    `Category - SRA specific`,
    `Category - other`,
    -gene_id_no_version # Remove the temporary join key
  )

# --- Apply Categorization Rules (DownDEGs) ---
cat("Applying categorization rules to DownDEGs...\n")
categorized_down_data <- combined_down_data %>%
  dplyr::mutate(
    # Initialize category columns
    `Category - SRA specific` = "",
    `Category - other` = "",
    
    # Apply rules using case_when (prioritized order)
    `Category - other` = dplyr::case_when(
      `023_24hr_Down` == 1 ~ "ns",
      # Check for 025=1 AND any 021 is non-zero (sum > 0)
      `025_24hr_Down` == 1 & (`021_10hr_Down` + `021_24hr_Down` + `021_48hr_Down` > 0) ~ "PCD-RFP & SRA",
      # Check for 025=1 AND all 021 are zero (sum = 0)
      `025_24hr_Down` == 1 & (`021_10hr_Down` + `021_24hr_Down` + `021_48hr_Down` == 0) ~ "PCD-RFP",
      TRUE ~ `Category - other` # Keep existing value (or default "") if no 023/025 match
    ),
    
    # Apply SRA specific rules ONLY if 023 and 025 conditions didn't apply
    `Category - SRA specific` = dplyr::case_when(
      `023_24hr_Down` == 0 & `025_24hr_Down` == 0 & `021_10hr_Down` == 1 & `021_24hr_Down` == 0 & `021_48hr_Down` == 0 ~ "Early transient",
      `023_24hr_Down` == 0 & `025_24hr_Down` == 0 & `021_10hr_Down` == 1 & `021_24hr_Down` == 1 & `021_48hr_Down` == 0 ~ "Early transient",
      `023_24hr_Down` == 0 & `025_24hr_Down` == 0 & `021_10hr_Down` == 1 & `021_24hr_Down` == 1 & `021_48hr_Down` == 1 ~ "Early sustained",
      `023_24hr_Down` == 0 & `025_24hr_Down` == 0 & `021_10hr_Down` == 1 & `021_24hr_Down` == 0 & `021_48hr_Down` == 1 ~ "Stochastic",
      `023_24hr_Down` == 0 & `025_24hr_Down` == 0 & `021_10hr_Down` == 0 & `021_24hr_Down` == 1 & `021_48hr_Down` == 1 ~ "Late sustained",
      `023_24hr_Down` == 0 & `025_24hr_Down` == 0 & `021_10hr_Down` == 0 & `021_24hr_Down` == 1 & `021_48hr_Down` == 0 ~ "Late transient",
      `023_24hr_Down` == 0 & `025_24hr_Down` == 0 & `021_10hr_Down` == 0 & `021_24hr_Down` == 0 & `021_48hr_Down` == 1 ~ "Delayed",
      TRUE ~ `Category - SRA specific` # Keep existing value (or default "") otherwise
    )
  ) %>%
  # Select and reorder columns for the final output (use dplyr::select)
  dplyr::select(
    transcript_id,
    Symbol = symbol, # Rename symbol column
    hg38_chr,        # Added
    hg38_start,      # Added
    hg38_end,        # Added
    hg38_strand,     # Added
    `023_24hr_Down`,
    `025_24hr_Down`,
    `021_10hr_Down`,
    `021_24hr_Down`,
    `021_48hr_Down`,
    `Category - SRA specific`,
    `Category - other`,
    -gene_id_no_version # Remove the temporary join key
  )

# --- Export Results ---

# Export the categorized UpDEG data
cat("\nWriting UpDEG output to:", output_up_file, "\n")
# Create directory if it doesn't exist
output_dir_up <- dirname(output_up_file)
if (!dir.exists(output_dir_up)) {
  dir.create(output_dir_up, recursive = TRUE)
}
# Replace NA values in character columns with empty strings for cleaner CSV output
categorized_up_data_out <- categorized_up_data %>%
  dplyr::mutate(dplyr::across(where(is.character), ~replace_na(., "")))
readr::write_csv(categorized_up_data_out, output_up_file, na = "") # Write remaining NAs (numeric) as empty string


# Export the categorized DownDEG data
cat("Writing DownDEG output to:", output_down_file, "\n")
# Create directory if it doesn't exist
output_dir_down <- dirname(output_down_file)
if (!dir.exists(output_dir_down)) {
  dir.create(output_dir_down, recursive = TRUE)
}
# Replace NA values in character columns with empty strings for cleaner CSV output
categorized_down_data_out <- categorized_down_data %>%
  dplyr::mutate(dplyr::across(where(is.character), ~replace_na(., "")))
readr::write_csv(categorized_down_data_out, output_down_file, na = "") # Write remaining NAs (numeric) as empty string

cat("\nProcessing complete.\n")