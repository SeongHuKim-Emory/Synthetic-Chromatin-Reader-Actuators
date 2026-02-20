# --- 1. Setup: Install and Load Required Packages ---

# Check if BiocManager is installed, install if not
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# List of required packages
required_packages <- c("readr",      # For efficient CSV reading/writing
                       "dplyr",      # For data manipulation verbs (select, mutate)
                       "stringr",    # For string manipulation (removing version)
                       "AnnotationDbi", # Core annotation functions
                       "org.Hs.eg.db")  # Human gene annotations (change if not human)

# Check and install missing CRAN/Bioconductor packages
for (pkg in required_packages) {
  if (pkg %in% c("AnnotationDbi", "org.Hs.eg.db")) { # Bioconductor packages
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg)
  } else { # CRAN packages
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg)
  }
  # Load the package after ensuring it's installed
  library(pkg, character.only = TRUE)
}

# --- 2. Define File Paths ---
input_file <- "../public/RNAseq/gexp_counts.csv"
output_file <- "../public/RNAseq/gexp_counts_symbol.csv"

# --- 3. Import the Count Data ---
message("Reading input file: ", input_file)
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
counts_data <- readr::read_csv(input_file, show_col_types = FALSE)

# --- 4. Prepare Gene IDs for Mapping ---
# Extract the transcript IDs (which look like Ensembl Gene IDs with versions)
ensembl_ids_version <- counts_data$transcript_id

# Remove the version suffix (e.g., ".14") to get base Ensembl Gene IDs
# Assumes the format is ENSGXXXXXXXXXXX.version
ensembl_ids_base <- stringr::str_remove(ensembl_ids_version, "\\.\\d+$")
# Alternative using base R:
# ensembl_ids_base <- gsub("\\.\\d+$", "", ensembl_ids_version)

message("Extracted ", length(unique(ensembl_ids_base)), " unique base Ensembl IDs for mapping.")

# --- 5. Map Ensembl IDs to Gene Symbols ---
# Use AnnotationDbi::mapIds
# We map from ENSEMBL (base ID) to SYMBOL using the human annotation database
# `multiVals = "first"`: If an Ensembl ID maps to multiple symbols, take the first one.
#   Other options: "list" (returns a list), "filter" (removes multi-matches),
#   "CharacterList" (returns a list-like object), "asNA" (returns NA for multi-matches)
message("Mapping Ensembl IDs to Gene Symbols using org.Hs.eg.db...")
gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys = ensembl_ids_base,
                                      column = "SYMBOL",
                                      keytype = "ENSEMBL",
                                      multiVals = "first")

# `mapIds` returns a named vector where names are the input keys (ensembl_ids_base)
# We need to match these back to the original data frame order

# --- 6. Add the Symbol Column to the Data Frame ---
# Create a new column 'symbol'
# For each row, find its base Ensembl ID and look up the corresponding symbol
counts_data_symbol <- counts_data %>%
  dplyr::mutate(
    ensembl_base = stringr::str_remove(transcript_id, "\\.\\d+$"), # Re-calculate base ID for matching
    symbol = gene_symbols[ensembl_base], # Look up symbol using the named vector
    .after = transcript_id # Place the new columns right after transcript_id
  ) %>%
  dplyr::select(-ensembl_base) # Remove the temporary base ID column

# Optional: Report how many IDs were successfully mapped
mapped_count <- sum(!is.na(counts_data_symbol$symbol))
total_count <- nrow(counts_data_symbol)
message("Successfully mapped ", mapped_count, " out of ", total_count, " IDs (",
        round(100 * mapped_count / total_count, 1), "%).")
unmapped_count <- total_count - mapped_count
if (unmapped_count > 0) {
  message(unmapped_count, " IDs could not be mapped and have NA in the symbol column.")
}

# Ensure the column order is exactly: transcript_id, symbol, then counts...
# dplyr::mutate with .after should handle this, but double-check with select if needed:
# counts_data_symbol <- counts_data_symbol %>%
#   dplyr::select(transcript_id, symbol, dplyr::everything())

# --- 7. Export the Modified Data Frame ---
message("Writing output file: ", output_file)
readr::write_csv(counts_data_symbol, output_file)

message("Script finished successfully.")

# --- Optional: Session Info for Reproducibility ---
# sessionInfo()