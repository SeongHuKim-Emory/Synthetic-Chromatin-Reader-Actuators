setwd("C:/Users/ski2959/Desktop/NIH R21/Codes/src")

################################################################################
# Generate the following files
# ../public/RNAseq/UpDEG_category.bed
# ../public/RNAseq/DownDEG_category.bed
################################################################################

# R Script for Converting RNA-seq data to BED format

# Function to process a single CSV file and convert it to BED format
process_deg_to_bed <- function(input_csv_path, output_bed_path) {
  # 1. Import the CSV file
  # Assuming the CSV has a header and uses commas as separators.
  # stringsAsFactors = FALSE is important to keep character columns as characters.
  # na.strings = c("", "NA") will treat empty strings as NA, which helps in checking.
  tryCatch({
    deg_data <- read.csv(input_csv_path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
    cat("Successfully read:", input_csv_path, "\n")
    cat("Number of rows initially:", nrow(deg_data), "\n")
  }, error = function(e) {
    stop("Error reading file ", input_csv_path, ": ", e$message)
  })
  
  # Check if essential columns exist
  required_cols <- c("hg38_chr", "hg38_start", "hg38_end", "Symbol", "transcript_id")
  if (!all(required_cols %in% names(deg_data))) {
    missing_cols <- required_cols[!required_cols %in% names(deg_data)]
    stop("Missing required columns in ", input_csv_path, ": ", paste(missing_cols, collapse=", "))
  }
  
  # 2. Create the new table (data frame)
  # Filter out rows where 'hg38_chr' is NA or an empty string
  # Note: na.strings in read.csv already converted empty strings in hg38_chr to NA
  deg_data_filtered <- deg_data[!is.na(deg_data$hg38_chr), ]
  cat("Number of rows after filtering for 'hg38_chr':", nrow(deg_data_filtered), "\n")
  
  if (nrow(deg_data_filtered) == 0) {
    cat("No valid rows found in", input_csv_path, "after filtering for 'hg38_chr'. Output file will be empty or not created.\n")
    # Optionally, create an empty file or skip writing
    if (!is.null(output_bed_path)) {
      write.table(data.frame(), file = output_bed_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      cat("Empty BED file created:", output_bed_path, "\n")
    }
    return() # Exit the function if no data to process
  }
  
  # 3. Fill in the new table columns
  # Column 1: hg38_chr
  bed_col1_chr <- deg_data_filtered$hg38_chr
  
  # Column 2: hg38_start
  # Ensure start positions are integers
  bed_col2_start <- as.integer(deg_data_filtered$hg38_start)
  if(any(is.na(bed_col2_start))) {
    cat("Warning: Some 'hg38_start' values were NA or could not be converted to integer in", input_csv_path, "\n")
    # Handle rows with NA start if necessary, e.g., by removing them or setting a default
    # For now, we'll proceed, and write.table will write NA as "NA"
  }
  
  
  # Column 3: hg38_end
  # Ensure end positions are integers
  bed_col3_end <- as.integer(deg_data_filtered$hg38_end)
  if(any(is.na(bed_col3_end))) {
    cat("Warning: Some 'hg38_end' values were NA or could not be converted to integer in", input_csv_path, "\n")
  }
  
  # Column 4: Symbol (or transcript_id if Symbol is NA/empty)
  # is.na() checks for NA. `nchar(x) == 0` can check for empty strings if they weren't converted to NA.
  # Since we used na.strings in read.csv, NA check is sufficient.
  bed_col4_name <- ifelse(is.na(deg_data_filtered$Symbol),
                          deg_data_filtered$transcript_id,
                          deg_data_filtered$Symbol)
  
  # Ensure no NA values in the name column after substitution, if transcript_id could also be NA
  if(any(is.na(bed_col4_name))) {
    cat("Warning: Some 'name' values (Symbol/transcript_id) are NA in", input_csv_path, "\n")
    # Decide how to handle: e.g., replace NA with a placeholder like "unknown"
    # bed_col4_name[is.na(bed_col4_name)] <- "unknown" 
  }
  
  
  # Combine into a data frame
  bed_output_df <- data.frame(
    chr = bed_col1_chr,
    start = bed_col2_start,
    end = bed_col3_end,
    name = bed_col4_name,
    stringsAsFactors = FALSE
  )
  
  # Remove rows where start or end might have become NA due to conversion issues,
  # or if they were NA in the original data and not filtered.
  # BED format requires start and end to be integers.
  bed_output_df <- bed_output_df[!is.na(bed_output_df$start) & !is.na(bed_output_df$end), ]
  cat("Number of rows after ensuring start/end are not NA:", nrow(bed_output_df), "\n")
  
  
  # 4. Export the new table
  # write.table is used for flexibility, sep="\t" for tab-delimited.
  # row.names=FALSE and col.names=FALSE as per BED format requirements.
  # quote=FALSE ensures values are not quoted.
  if (nrow(bed_output_df) > 0) {
    tryCatch({
      write.table(bed_output_df, file = output_bed_path, sep = "\t",
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
      cat("Successfully wrote BED file:", output_bed_path, "\n")
    }, error = function(e) {
      stop("Error writing file ", output_bed_path, ": ", e$message)
    })
  } else {
    cat("No data to write for", output_bed_path, "after all processing steps.\n")
    # Optionally, create an empty file if it wasn't created earlier
    if (!file.exists(output_bed_path)) {
      write.table(data.frame(), file = output_bed_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      cat("Empty BED file created:", output_bed_path, "\n")
    }
  }
}

# --- Main script execution ---
output_dir_rnaseq <- "../public/RNAseq/" # Renamed to avoid conflict with other output_dir
if (!dir.exists(output_dir_rnaseq)) {
  dir.create(output_dir_rnaseq, recursive = TRUE)
  cat("Created directory:", output_dir_rnaseq, "\n")
}


input_up_path <- "../public/RNAseq/UpDEG_category.csv"
output_up_path <- file.path(output_dir_rnaseq, "UpDEG_category.bed") # Use file.path for robustness

input_down_path <- "../public/RNAseq/DownDEG_category.csv"
output_down_path <- file.path(output_dir_rnaseq, "DownDEG_category.bed")

# Process the UpREG file
cat("\nProcessing UpDEG file...\n")
process_deg_to_bed(input_up_path, output_up_path)

# Process the DownREG file
cat("\nProcessing DownDEG file...\n")
process_deg_to_bed(input_down_path, output_down_path)

cat("\nScript finished processing DEG files.\n") # Clarified message

################################################################################
# Generate following files
# MCF7_hg38_TAD_Boundaries.bed
################################################################################

# R script for processing TAD boundary data

# 1. Import files
# Define the input file path
input_file_tad <- "../public/Public Dataset/HiC/GSE66733_MCF7_HiC/MCF7_hg19_hg38_TAD_Boundaries.csv" # Renamed to avoid conflict

# Read the CSV file into a data frame
# The header is present, and the separator is a comma
tryCatch({
  tad_data <- read.csv(input_file_tad, header = TRUE, stringsAsFactors = FALSE, na.strings = "NA")
  print("Successfully imported the CSV file for TAD boundaries.") # Clarified message
  print("First few rows of the imported TAD data:")
  print(head(tad_data))
}, error = function(e) {
  stop("Error reading the TAD CSV file: ", e$message)
})

# 2. Remove rows if the value of the third column (hg38_TAD_Boundary) is "NA"
# The third column in R will be tad_data[,3] or by its name tad_data$hg38_TAD_Boundary
# We need to ensure that the column name is correct as per the CSV header.
# From the example, the column names are: chr, hg19_TAD_Boundary, hg38_TAD_Boundary, Insulation_Score

# Check if the column 'hg38_TAD_Boundary' exists
if (!"hg38_TAD_Boundary" %in% colnames(tad_data)) {
  stop("The column 'hg38_TAD_Boundary' was not found in the TAD input file. Please check the CSV header.")
}

# Filter out rows where hg38_TAD_Boundary is NA
# The na.strings = "NA" in read.csv should have already converted "NA" strings to actual NA values.
filtered_tad_data <- tad_data[!is.na(tad_data$hg38_TAD_Boundary), ]

print("First few rows of the TAD data after removing NA values in 'hg38_TAD_Boundary':")
print(head(filtered_tad_data))
print(paste("Number of rows in TAD data before filtering:", nrow(tad_data)))
print(paste("Number of rows in TAD data after filtering:", nrow(filtered_tad_data)))

# 3. Create a new table
# 1st column: chromosome number (from 'chr' column)
# 2nd column: [hg38_TAD_Boundary] value - 20000
# 3rd column: [hg38_TAD_Boundary] value + 20000

# Ensure the 'hg38_TAD_Boundary' column is numeric for calculations
if (!is.numeric(filtered_tad_data$hg38_TAD_Boundary)) {
  # Attempt to convert to numeric, coercing errors to NA (though NAs should be removed)
  filtered_tad_data$hg38_TAD_Boundary <- as.numeric(as.character(filtered_tad_data$hg38_TAD_Boundary))
  # Re-check for NAs introduced by coercion, if any, and remove them
  # This step is a safeguard; ideally, the column is already numeric or cleanly convertible.
  if(any(is.na(filtered_tad_data$hg38_TAD_Boundary))) {
    warning("NAs were introduced during conversion of 'hg38_TAD_Boundary' to numeric. Removing these rows.")
    filtered_tad_data <- filtered_tad_data[!is.na(filtered_tad_data$hg38_TAD_Boundary), ]
  }
}

# Create the new table
# The problem statement implies the chromosome identifiers (e.g., "chr1", "chrX") should be kept as is.
new_table_tad <- data.frame( # Renamed to avoid conflict
  chromosome = filtered_tad_data$chr,
  start_coord = filtered_tad_data$hg38_TAD_Boundary - 20000,
  end_coord = filtered_tad_data$hg38_TAD_Boundary + 20000,
  stringsAsFactors = FALSE # Keep chromosome as character
)

print("First few rows of the new TAD table:")
print(head(new_table_tad))

# 4. Export the new table
# Define the output file path
output_file_tad_bed <- "../public/Public Dataset/HiC/GSE66733_MCF7_HiC/MCF7_hg38_TAD_Boundaries.bed" # Renamed to avoid conflict

# Write the table to a tab-delimited file without row names and column headers (typical for BED files)
tryCatch({
  write.table(new_table_tad,
              file = output_file_tad_bed,
              sep = "\t",           # Tab delimited
              row.names = FALSE,  # Do not write row names
              col.names = FALSE,  # Do not write column names (standard for BED)
              quote = FALSE)      # Do not quote strings
  print(paste("Successfully exported the new TAD table to:", output_file_tad_bed))
}, error = function(e) {
  stop("Error writing the TAD BED file: ", e$message)
})

print("Script finished processing TAD boundary data.") # Clarified message

################################################################################
# Generate bar graphs
################################################################################

# R Script for Generating Multi-track Bar Charts from Genomic Data

# --- 1. Load Required Libraries ---
# Ensure these packages are installed: install.packages(c("readr", "dplyr", "ggplot2", "tidyr", "stringr"))
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# --- Define Fixed Y-axis Limits for Specific Tracks ---
# These limits will be applied to the respective tracks across all chromosome plots.
track_y_axis_limits <- list(
  "transcript_Count"    = list(min = 0, max = 5),
  "021_Peak_Count"      = list(min = 0, max = 100),
  "025_Peak_Count"      = list(min = 0, max = 350),
  "H3K27me3_Peak_Count" = list(min = 0, max = 650),
  "H3K27ac_Peak_Count"  = list(min = 0, max = 250),
  "H3K4me3_Peak_Count"  = list(min = 0, max = 150)
)

# --- Define Bar Colors for Specific Tracks ---
# These colors will be applied to the respective tracks.
base_track_colors <- list(
  "transcript_Count"    = "coral1",
  "021_Peak_Count"      = "darkorange",
  "025_Peak_Count"      = "indianred",
  "H3K27me3_Peak_Count" = "steelblue",
  "H3K27ac_Peak_Count"  = "darkred",
  "H3K4me3_Peak_Count"  = "palevioletred4"
)

# --- 2. Define File Paths ---
input_file_multiomics <- "../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks_Cytobands.csv" # Renamed
output_dir_multiomics <- "../public/Multiomics/" # Renamed


# --- 3. Create Output Directory if it doesn't exist ---
if (!dir.exists(output_dir_multiomics)) {
  dir.create(output_dir_multiomics, recursive = TRUE, showWarnings = FALSE)
  message(paste("Created output directory:", output_dir_multiomics))
}

# --- 4. Read the Input Data ---
# read_csv from readr handles column names with special characters or leading numbers well.
message(paste("Reading multiomics data from:", input_file_multiomics))
tryCatch({
  full_data <- readr::read_csv(input_file_multiomics, show_col_types = FALSE)
}, error = function(e) {
  stop(paste("Error reading multiomics CSV file:", e$message, "\nPlease ensure the file path is correct, the file exists, and R has permission to access it."))
})

# --- 5. Identify Columns for Plotting ---
# These are the columns that will form the different tracks in the bar chart.
track_columns <- c(
  "transcript_Count",
  "021_Peak_Count",  
  "025_Peak_Count",  
  "H3K27me3_Peak_Count",
  "H3K27ac_Peak_Count",
  "H3K4me3_Peak_Count"
)

# Verify that all track columns exist in the loaded data
actual_colnames <- colnames(full_data)

# --- 6. Get Unique Chromosome Names ---
# Ensure the 'chr' column exists
if (!"chr" %in% actual_colnames) {
  stop("The required 'chr' column was not found in the multiomics input CSV file. Please check the CSV header.")
}
unique_chromosomes <- unique(full_data$chr)
if (length(unique_chromosomes) == 0 || all(is.na(unique_chromosomes))) {
  stop("No chromosome identifiers found in the 'chr' column of multiomics data or they are all NA.")
}
message(paste("Found chromosomes in multiomics data:", paste(unique_chromosomes, collapse=", ")))

# --- 7. Loop Through Each Chromosome to Generate Plots ---
for (current_chr_name in unique_chromosomes) {
  message(paste("Processing chromosome:", current_chr_name))
  
  # Filter data for the current chromosome
  chr_data <- full_data %>%
    dplyr::filter(chr == current_chr_name) # Using dplyr::filter for clarity, though not strictly necessary if no conflict
  
  if (nrow(chr_data) == 0) {
    message(paste("No data found for chromosome:", current_chr_name, "- Skipping."))
    next
  }
  
  # Ensure the 'TAD_Num' column exists for the current chromosome's data
  if (!"TAD_Num" %in% names(chr_data)) {
    warning(paste("'TAD_Num' column not found for chromosome:", current_chr_name, ". This chromosome will be skipped or plotting might fail."))
    next  
  }
  
  # Create a mapping from desired track_columns to actual column names in chr_data
  cols_for_pivot <- c()
  final_track_order_names <- c()  
  current_plot_limits_list <- list() # To store data.frames for y-axis limits
  
  # Map from actual column name found in data (pname_found) to original base name (col_base_name)
  # This is crucial for applying colors based on the original track names.
  pname_to_basename_map <- list()  
  
  for (col_base_name in track_columns) {
    potential_names <- c(col_base_name,  
                         make.names(col_base_name), 
                         paste0("X", col_base_name)) 
    
    found_match = FALSE
    pname_found = "" # Store the matched column name
    for (pname in unique(potential_names)) {  
      if (pname %in% names(chr_data)) {
        cols_for_pivot <- c(cols_for_pivot, pname)
        final_track_order_names <- c(final_track_order_names, pname)
        pname_found <- pname # Store the name as it was found in the data
        
        # Store the mapping from the actual data column name to the original base name
        pname_to_basename_map[[pname_found]] <- col_base_name
        
        found_match = TRUE
        break  
      }
    }
    if (!found_match) {
      warning(paste("Column '", col_base_name, "' or a variant of it was not found in the data for chromosome ", current_chr_name, ". It will be skipped.", sep=""))
    } else {
      # If a match was found, check if this base name has defined y-axis limits
      if (col_base_name %in% names(track_y_axis_limits)) {
        limits <- track_y_axis_limits[[col_base_name]]
        # Use pname_found (actual name in data) for the 'Track' identifier in the limits data frame
        current_plot_limits_list[[length(current_plot_limits_list) + 1]] <-
          data.frame(Track = pname_found,  
                     min_y_val = limits$min,  
                     max_y_val = limits$max,  
                     stringsAsFactors = FALSE)
      }
    }
  }
  
  if(length(cols_for_pivot) == 0){
    message(paste("No data columns found for plotting for chromosome:", current_chr_name, "- Skipping."))
    next
  }
  
  # Combine the list of limit data.frames into one, if any limits were defined for current tracks
  facet_limits_df <- NULL
  if (length(current_plot_limits_list) > 0) {
    facet_limits_df <- do.call(rbind, current_plot_limits_list)
    # Ensure 'Track' column is a factor with the same levels as in chr_data_long for proper matching
    # This will use final_track_order_names which dictates the facet order
    facet_limits_df$Track <- factor(facet_limits_df$Track, levels = final_track_order_names)
  }
  
  # Select relevant columns and pivot to long format
  # MODIFIED LINE: Added dplyr:: to select
  chr_data_long <- chr_data %>%
    dplyr::select(TAD_Num, all_of(cols_for_pivot)) %>% 
    pivot_longer(
      cols = all_of(cols_for_pivot),
      names_to = "Track",
      values_to = "Count",
      values_drop_na = FALSE  
    )
  
  # Replace NA counts with 0
  chr_data_long$Count[is.na(chr_data_long$Count)] <- 0
  
  # Ensure the 'Track' column is a factor with the specified order
  chr_data_long$Track <- factor(chr_data_long$Track, levels = final_track_order_names)
  
  # --- 8. Generate the Plot ---
  plot_title <- paste("MCF7 SRA TAD Data - Chromosome:", sub("chr", "", current_chr_name))
  
  # Prepare colors for the current plot
  # Uses the mapping (pname_to_basename_map) to get the original base name for color lookup
  current_plot_colors <- sapply(final_track_order_names, function(pname_in_data) {
    original_basename <- pname_to_basename_map[[pname_in_data]]
    color <- base_track_colors[[original_basename]]
    if (is.null(color)) { # If no color defined for this original base name
      # Check if the pname_in_data itself is a key in base_track_colors (fallback for direct match)
      color <- base_track_colors[[pname_in_data]]
      if(is.null(color)) "grey50" else color # Default if still not found
    } else {
      color
    }
  }, USE.NAMES = TRUE)
  # Ensure the names of current_plot_colors match the levels of chr_data_long$Track
  names(current_plot_colors) <- final_track_order_names
  
  
  p <- ggplot(chr_data_long, aes(x = factor(TAD_Num), y = Count)) +
    geom_bar(stat = "identity", aes(fill = Track), width = 0.7) + # MODIFIED: aes(fill = Track)
    scale_fill_manual(values = current_plot_colors, guide = "none") # ADDED: Apply custom colors, hide legend
  
  # MODIFIED: Add geom_blank layers if facet_limits_df is available and has rows
  # Use inherit.aes = FALSE to prevent inheriting x aesthetic from global ggplot call
  if (!is.null(facet_limits_df) && nrow(facet_limits_df) > 0) {
    p <- p +  
      geom_blank(data = facet_limits_df, aes(y = min_y_val), inherit.aes = FALSE) +
      geom_blank(data = facet_limits_df, aes(y = max_y_val), inherit.aes = FALSE)
  }
  
  p <- p + facet_wrap(~ Track, ncol = 1, scales = "free_y", strip.position = "left") +
    labs(
      title = plot_title,
      x = paste("TAD Number on Chromosome", sub("chr", "", current_chr_name)),
      y = NULL  # Y-axis label is effectively handled by strip text
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),  
      axis.title.x = element_text(size = 20, margin = margin(t = 10)),
      strip.background = element_blank(),  
      strip.placement = "outside",       
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size=20),  
      panel.spacing.y = unit(2, "lines"),  
      plot.margin = margin(1, 1, 1, 1, "cm")
      # legend.position = "none" # Alternative way to hide legend, if not using guide="none" in scale
    )
  
  # --- 9. Save the Plot ---
  chr_identifier <- sub("chr", "", current_chr_name)  
  output_filename <- paste0(output_dir_multiomics, "MCF7_SRA_TAD_Chr", chr_identifier, ".jpg") # Used renamed output_dir
  
  # num_tads <- length(unique(chr_data_long$TAD_Num)) # No longer needed for dynamic sizing
  # num_tracks <- length(unique(chr_data_long$Track)) # No longer needed for dynamic sizing
  
  # Set fixed plot dimensions
  plot_width <- 42 
  plot_height <- 17 
  
  message(paste("Saving plot:", output_filename, "with dimensions", plot_width, "x", plot_height, "inches"))
  
  ggsave(
    filename = output_filename,
    plot = p,
    width = plot_width,
    height = plot_height,
    units = "in",
    dpi = 300,
    device = "jpeg"
  )
  message(paste("Saved plot:", output_filename))
}

message("Multiomics bar graph processing complete.") # Clarified message


################################################################################

################################################################################
