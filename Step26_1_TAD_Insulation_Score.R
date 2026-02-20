setwd("X:/haynes-lab/Kim_SH/Projects/NIH R21/Codes/src")

# R script for merging TADs based on insulation score thresholds
# R Version: 4.4.0
# Required Package: data.table

# --- 1. Setup ---

# Install data.table if not already installed
if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

# Define file paths
input_file <- "../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks_Cytobands.csv"
output_dir <- "../public/Multiomics/Step26_TAD_Insulation_Score"

# Create the output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. Load Data ---

cat("Loading data from:", input_file, "\n")
# Use fread for fast CSV reading
tryCatch({
  dt_orig <- fread(input_file)
}, error = function(e) {
  stop("Error reading file: ", e$message, "\nCheck if the file path is correct:\n", input_file)
})

# Define columns to remove as requested
cols_to_remove <- c(
  "021_Peak_Names", "025_Peak_Names", "H3K4me1_Peak_Count",
  "H3K4me2_Peak_Count", "H3K4me3_Peak_Count", "H3K9me3_Peak_Count",
  "H3K27ac_Peak_Count", "H3K27me3_Peak_Count", "H3K36me3_Peak_Count"
)

# --- 3. Define Thresholds ---

# Generate the sequence of insulation thresholds
insulation_thresholds <- seq(0, 1, by = 0.05)

cat("Starting processing for", length(insulation_thresholds), "thresholds...\n")

# --- 4. Main Processing Loop ---

for (threshold in insulation_thresholds) {
  
  cat("Processing threshold:", threshold, "\n")
  
  # a. Make a copy of the original data for this iteration
  dt_proc <- copy(dt_orig)
  
  # b. Remove the specified columns
  # data.table way to remove columns by name
  dt_proc[, (cols_to_remove) := NULL]
  
  # c. Perform the merging logic
  # We group by chromosome ('chr') because merges should not cross chromosomes.
  # For each chromosome group (.SD), we perform the merge logic.
  
  dt_merged <- dt_proc[, {
    
    # .SD is the data.table subset for the current 'chr'
    # .N is the number of rows in this subset
    if (.N == 0) {
      # Handle empty chromosome group if it ever occurs
      .SD
    } else {
      
      # Get the insulation scores that separate rows
      # The score at row 'i' separates row 'i' from 'i+1'
      # We don't need the score from the very last row (.N)
      scores <- .SD$`Insulation_Score_3'`
      
      # Find the indices where a merge *stops* (i.e., the boundary is met or exceeded)
      # These are the rows *before* the boundary.
      break_indices <- which(scores[-.N] >= threshold)
      
      # Define the start and end rows for each *new* merged group
      group_starts <- c(1, break_indices + 1)
      group_ends <- c(break_indices, .N)
      
      # Create a list to hold the new merged rows (as data.tables)
      merged_rows_list <- vector("list", length(group_starts))
      
      for (k in 1:length(group_starts)) {
        start_idx <- group_starts[k]
        end_idx <- group_ends[k]
        
        # Get all rows that will be part of this new merged row
        rows_to_merge <- .SD[start_idx:end_idx]
        
        # Get the first and last row for easy access
        first_row <- rows_to_merge[1]
        last_row <- rows_to_merge[.N]
        
        # --- Create the new aggregated values ---
        
        # 1. TAD_Num: "57" or "57 to 60"
        tad_num_str <- if (start_idx == end_idx) {
          as.character(first_row$TAD_Num)
        } else {
          paste(first_row$TAD_Num, "to", last_row$TAD_Num)
        }
        
        # 2. cytoband: "chr17q21.32" or "chr17q21.32-chr17q21.33"
        #    Note: The example logic is complex. This implementation pastes
        #    the first and last cytoband names if they are different.
        first_cyto <- first_row$cytoband
        last_cyto <- last_row$cytoband
        cytoband_str <- if (first_cyto == last_cyto) {
          first_cyto
        } else {
          paste(first_cyto, last_cyto, sep = "-")
        }
        
        # 3. String concatenation (transcript_ids, Symbols)
        #    Filter out empty strings ("") before pasting
        transcripts <- paste(
          rows_to_merge$transcript_ids[rows_to_merge$transcript_ids != ""], 
          collapse = " "
        )
        symbols <- paste(
          rows_to_merge$Symbols[rows_to_merge$Symbols != ""], 
          collapse = " "
        )
        
        # 4. Numeric Summation
        sum_transcript_Count <- sum(rows_to_merge$transcript_Count, na.rm = TRUE)
        sum_021_Peak_Count <- sum(rows_to_merge$`021_Peak_Count`, na.rm = TRUE)
        sum_025_Peak_Count <- sum(rows_to_merge$`025_Peak_Count`, na.rm = TRUE)
        
        # 5. Build the new row as a data.table
        merged_rows_list[[k]] <- data.table(
          chr = first_row$chr,
          TAD_Num = tad_num_str,
          cytoband = cytoband_str,
          hg38_TAD_Start = first_row$hg38_TAD_Start,
          hg38_TAD_End = last_row$hg38_TAD_End,
          `Insulation_Score_5'` = first_row$`Insulation_Score_5'`,
          `Insulation_Score_3'` = last_row$`Insulation_Score_3'`,
          transcript_ids = transcripts,
          Symbols = symbols,
          transcript_Count = sum_transcript_Count,
          `021_Peak_Count` = sum_021_Peak_Count,
          `025_Peak_Count` = sum_025_Peak_Count
        )
      }
      
      # Combine all the new merged rows for this chromosome
      rbindlist(merged_rows_list)
    }
    
  }, by = chr] # End of grouping by 'chr'
  
  # --- 5. Export the Result ---
  
  # Format threshold for file name (e.g., 0.3 -> 0_3, 0.05 -> 0_05)
  threshold_str_file <- gsub("\\.", "_", as.character(threshold))
  
  # Construct table name and file path
  table_name <- paste0("MCF7_TAD_Analysis_Insulation_", threshold_str_file)
  output_file <- file.path(output_dir, paste0(table_name, ".csv"))
  
  # Write to CSV using fwrite
  # Using quote = TRUE to match the quoted text in the example
  fwrite(dt_merged, file = output_file, row.names = FALSE, quote = TRUE)
  
} # End of main loop

cat("Processing complete. Files are saved in:", output_dir, "\n")

################################################################################

################################################################################

