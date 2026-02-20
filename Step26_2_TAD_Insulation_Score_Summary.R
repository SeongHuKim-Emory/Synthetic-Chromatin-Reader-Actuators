# R script for TAD Insulation Score Analysis
# Compatible with R version 4.4.0

# Set working directory to the script's location (optional, good practice)
# If running in RStudio, this can be set automatically.
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- 1. Import Files ---

# Define the base path and file components
base_path <- "../public/Multiomics/Step26_TAD_Insulation_Score/"

# Suffixes for the files, corresponding to thresholds 0, 0.1, ..., 0.9, 1
file_suffixes <- c("0", "0_1", "0_2", "0_3", "0_4", "0_5", "0_6", "0_7", "0_8", "0_9", "1")

# Corresponding threshold values as strings (for naming the list)
threshold_keys <- c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1")

# Generate the full file names and paths
file_names <- paste0("MCF7_TAD_Analysis_Insulation_", file_suffixes, ".csv")
file_paths <- file.path(base_path, file_names)

# Import all CSV files into a list
# We use read.csv, which is suitable for the provided structure.
# We set stringsAsFactors = FALSE for modern R compatibility.
data_list <- lapply(file_paths, function(path) {
  tryCatch({
    read.csv(path, stringsAsFactors = FALSE)
  }, error = function(e) {
    warning("Could not read file: ", path, ". Error: ", e$message)
    return(NULL) # Return NULL if file is missing or corrupt
  })
})

# Name the list elements using the threshold keys for easy access
names(data_list) <- threshold_keys

cat("Step 1: Successfully imported and processed", length(data_list), "files.\n")

# --- 2. Create Insulation_Score_Summary Table ---

# Initialize the summary data.frame with the specified columns and default values.
# R will recycle the 0 and 0.0 values for all 11 rows.
Insulation_Score_Summary <- data.frame(
  Insulation_threshold = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
  Total_TADs = 0,
  Mean_TAD_Size_bp = 0.0,
  SRA_TADs = 0,
  PCD_RFP_TADs = 0,
  UpDEG_TADs = 0,
  SRA_PCD_RFP_TADs = 0,
  SRA_UpDEG_TADs = 0,
  PCD_RFP_UpDEG_TADs = 0,
  SRA_PCD_RFP_UpDEG_TADs = 0
)

cat("Step 2: Created initial Insulation_Score_Summary table.\n")

# --- 3. Populate the Summary Table ---

# Loop through each row of the summary table
for (i in 1:nrow(Insulation_Score_Summary)) {
  
  # Get the insulation threshold for the current row
  threshold_val <- Insulation_Score_Summary$Insulation_threshold[i]
  
  # Convert the numeric threshold (e.g., 0.1) to its string key (e.g., "0.1")
  threshold_key <- as.character(threshold_val)
  
  # Retrieve the corresponding data frame from our list
  current_data <- data_list[[threshold_key]]
  
  # Check if data was loaded successfully and is not empty
  if (!is.null(current_data) && nrow(current_data) > 0) {
    
    # NOTE: R's read.csv (with default check.names=TRUE) converts column names
    # that start with numbers. "021_Peak_Count" becomes "X021_Peak_Count"
    # and "025_Peak_Count" becomes "X025_Peak_Count".
    
    # Calculate Total_TADs (number of rows)
    Insulation_Score_Summary$Total_TADs[i] <- nrow(current_data)
    
    # Calculate Mean_TAD_Size_bp
    tad_sizes <- current_data$hg38_TAD_End - current_data$hg38_TAD_Start
    Insulation_Score_Summary$Mean_TAD_Size_bp[i] <- mean(tad_sizes, na.rm = TRUE)
    
    # Calculate SRA_TADs (021_Peak_Count > 0)
    Insulation_Score_Summary$SRA_TADs[i] <- sum(current_data$X021_Peak_Count > 0, na.rm = TRUE)
    
    # Calculate PCD_RFP_TADs (025_Peak_Count > 0)
    Insulation_Score_Summary$PCD_RFP_TADs[i] <- sum(current_data$X025_Peak_Count > 0, na.rm = TRUE)
    
    # Calculate UpDEG_TADs (transcript_Count > 0)
    Insulation_Score_Summary$UpDEG_TADs[i] <- sum(current_data$transcript_Count > 0, na.rm = TRUE)
    
    # Calculate SRA_PCD_RFP_TADs (021 > 0 AND 025 > 0)
    Insulation_Score_Summary$SRA_PCD_RFP_TADs[i] <- sum(
      current_data$X021_Peak_Count > 0 & current_data$X025_Peak_Count > 0,
      na.rm = TRUE
    )
    
    # Calculate SRA_UpDEG_TADs (021 > 0 AND transcript_Count > 0)
    Insulation_Score_Summary$SRA_UpDEG_TADs[i] <- sum(
      current_data$X021_Peak_Count > 0 & current_data$transcript_Count > 0,
      na.rm = TRUE
    )
    
    # Calculate PCD_RFP_UpDEG_TADs (025 > 0 AND transcript_Count > 0)
    Insulation_Score_Summary$PCD_RFP_UpDEG_TADs[i] <- sum(
      current_data$X025_Peak_Count > 0 & current_data$transcript_Count > 0,
      na.rm = TRUE
    )
    
    # Calculate SRA_PCD_RFP_UpDEG_TADs (All three > 0)
    Insulation_Score_Summary$SRA_PCD_RFP_UpDEG_TADs[i] <- sum(
      current_data$X021_Peak_Count > 0 & 
        current_data$X025_Peak_Count > 0 &
        current_data$transcript_Count > 0,
      na.rm = TRUE
    )
    
  } else {
    # If file was missing or empty, corresponding row will keep its 0s
    cat("  - No data found for threshold", threshold_key, ". Skipping calculations.\n")
  }
}

cat("Step 3: Successfully populated the summary table.\n")

# Display the first few rows of the result
cat("\n--- Summary Table (Head) ---\n")
print(head(Insulation_Score_Summary))
cat("-----------------------------\n\n")

# --- 4. Export Insulation_Score_Summary Table ---

output_path <- "../public/Multiomics/Step26_TAD_Insulation_Score/Insulation_Score_Summary.csv"

# Create directory if it doesn't exist (recursive = TRUE creates parent dirs)
output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Export the data.frame to a CSV file
# row.names = FALSE prevents R from writing an extra column for row numbers
write.csv(Insulation_Score_Summary, file = output_path, row.names = FALSE)

cat("Step 4: Successfully exported summary table to:\n", output_path, "\n")
cat("\nScript finished.\n")