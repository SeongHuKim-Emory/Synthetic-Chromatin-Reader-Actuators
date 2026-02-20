# R script for processing TAD insulation scores and RNA-seq categories
# Version: 4.4.0

# --- 0. Setup ---
# Set file paths
tad_file_path <- "../public/Multiomics/Step26_TAD_Insulation_Score/MCF7_TAD_Analysis_Insulation_0_4.csv"
deg_file_path <- "../public/RNAseq/UpDEG_category.csv"
output_file_path <- "../public/Multiomics/Step26_TAD_Insulation_Score/MCF7_TAD_Analysis_Insulation_0_4_Category.csv"

# Define the categories to be added
categories_to_add <- c(
  "Early transient", 
  "Early sustained", 
  "Late transient", 
  "Late sustained", 
  "Delayed", 
  "Stochastic"
)

cat("Starting script...\n")

# --- 1. Import Files ---
cat("Importing files...\n")

# Import TAD data
# check.names = FALSE prevents R from changing column names (e.g., "Insulation_Score_5'" to "Insulation_Score_5.")
tad_data <- read.csv(
  tad_file_path, 
  stringsAsFactors = FALSE, 
  check.names = FALSE
)

# Import DEG category data
# check.names = FALSE is important here for "Category - SRA specific"
deg_data <- read.csv(
  deg_file_path, 
  stringsAsFactors = FALSE, 
  check.names = FALSE
)

cat("Files imported successfully.\n")

# --- 2. Create New Table ---
cat("Creating new table 'MCF7_TAD_Analysis_Insulation_0_4_Category'...\n")
MCF7_TAD_Analysis_Insulation_0_4_Category <- tad_data

# --- 3. Add New Category Columns ---
cat("Adding new category columns and initializing to 0...\n")
for (category_name in categories_to_add) {
  # Use [[]] notation to create new columns with names containing spaces
  MCF7_TAD_Analysis_Insulation_0_4_Category[[category_name]] <- 0
}

# --- 4. Populate Category Counts ---
cat("Populating category counts based on transcript IDs...\n")

# Create a fast lookup map (named vector) for transcript_id -> category
# This is much faster than searching the deg_data data.frame inside the loop

# Filter DEG data to only include relevant transcripts and categories
deg_data_filtered <- deg_data[
  !is.na(deg_data$transcript_id) & 
    !is.na(deg_data[["Category - SRA specific"]]) & 
    deg_data[["Category - SRA specific"]] %in% categories_to_add, 
]

# Create the lookup map
# setNames(values, keys)
transcript_to_category_map <- setNames(
  deg_data_filtered[["Category - SRA specific"]], 
  deg_data_filtered$transcript_id
)

# Iterate through each row of the new TAD table
for (i in 1:nrow(MCF7_TAD_Analysis_Insulation_0_4_Category)) {
  
  # Get the string of transcript IDs for the current row
  ids_string <- MCF7_TAD_Analysis_Insulation_0_4_Category[i, "transcript_ids"]
  
  # Skip if the cell is NA or empty
  if (is.na(ids_string) || ids_string == "") {
    next
  }
  
  # Split the string by spaces and get unique transcript IDs
  # This handles the requirement to "ignore duplicate individual transcript_ids"
  transcript_ids_list <- strsplit(ids_string, " ")[[1]]
  unique_transcript_ids <- unique(transcript_ids_list)
  
  # Process each unique transcript ID
  for (t_id in unique_transcript_ids) {
    
    # Skip if the ID is an empty string (e.g., from double spaces)
    if (t_id == "") {
      next
    }
    
    # Look up the category using the fast map
    category <- transcript_to_category_map[[t_id]]
    
    # If a valid category was found (not NULL)
    if (!is.null(category)) {
      # Increment the count in the corresponding category column for row 'i'
      MCF7_TAD_Analysis_Insulation_0_4_Category[i, category] <- 
        MCF7_TAD_Analysis_Insulation_0_4_Category[i, category] + 1
    }
  }
}
cat("Category counts populated.\n")

# --- 5. Clean Up Duplicates and Recalculate Counts ---
cat("Cleaning duplicate 'transcript_ids' and 'Symbols', and recalculating 'transcript_Count'...\n")

# Define a helper function to process space-delimited strings
process_string_duplicates <- function(s) {
  # Return NA or 0 if input is NA
  if (is.na(s)) {
    return(NA)
  }
  # Return "" or 0 if input is empty
  if (s == "") {
    return("")
  }
  
  # Split, find unique, filter out empty strings, and paste back together
  s_split <- strsplit(s, " ")[[1]]
  s_unique <- unique(s_split)
  s_filtered <- s_unique[s_unique != ""] # Remove empty strings
  
  return(paste(s_filtered, collapse = " "))
}

# Apply the function to the 'transcript_ids' and 'Symbols' columns
MCF7_TAD_Analysis_Insulation_0_4_Category$transcript_ids <- sapply(
  MCF7_TAD_Analysis_Insulation_0_4_Category$transcript_ids,
  process_string_duplicates,
  USE.NAMES = FALSE
)

MCF7_TAD_Analysis_Insulation_0_4_Category$Symbols <- sapply(
  MCF7_TAD_Analysis_Insulation_0_4_Category$Symbols,
  process_string_duplicates,
  USE.NAMES = FALSE
)

# Define a helper function to count unique, non-empty transcripts
count_unique_transcripts <- function(s) {
  # Return 0 if input is NA or empty
  if (is.na(s) || s == "") {
    return(0)
  }
  
  # Split the (already unique) string and count non-empty elements
  s_split <- strsplit(s, " ")[[1]]
  s_filtered <- s_split[s_split != ""]
  
  return(length(s_filtered))
}

# Apply the counting function to the 'transcript_ids' column to update 'transcript_Count'
MCF7_TAD_Analysis_Insulation_0_4_Category$transcript_Count <- sapply(
  MCF7_TAD_Analysis_Insulation_0_4_Category$transcript_ids,
  count_unique_transcripts,
  USE.NAMES = FALSE
)

cat("Cleanup and recounting complete.\n")

# --- 6. Export the Final Table ---
cat(paste("Exporting final table to", output_file_path, "...\n"))

write.csv(
  MCF7_TAD_Analysis_Insulation_0_4_Category, 
  output_file_path, 
  row.names = FALSE, # Do not save the R row indices
  quote = TRUE       # Quote strings (standard CSV practice)
)

cat("Script finished successfully.\n")