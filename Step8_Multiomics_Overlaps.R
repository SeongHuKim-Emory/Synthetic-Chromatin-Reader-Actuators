# Set the working directory (ensure this path is correct for your system)
# setwd("C:/Users/ski2959/Desktop/NIH R21/Codes/src") # Uncomment and set if needed

################################################################################

# 1. Import files
file1_path <- "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2.broadPeak"
file2_path <- "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2.broadPeak"
peak_cols <- c("chr", "start", "end")

# Read DBN021 peaks
tryCatch({
  dbn021_peaks <- read.table(file1_path,
                             sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                             col.names = c(peak_cols, rep("NULL", 6)),
                             colClasses = c("character", "integer", "integer", rep("NULL", 6)))
  if (ncol(dbn021_peaks) > length(peak_cols)) {
    dbn021_peaks <- dbn021_peaks[, 1:length(peak_cols)]
  }
  colnames(dbn021_peaks) <- peak_cols
}, error = function(e) {
  stop(paste("Error reading file:", file1_path, "-", e$message))
})

# Read DBN025 peaks
tryCatch({
  dbn025_peaks <- read.table(file2_path,
                             sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                             col.names = c(peak_cols, rep("NULL", 6)),
                             colClasses = c("character", "integer", "integer", rep("NULL", 6)))
  if (ncol(dbn025_peaks) > length(peak_cols)) {
    dbn025_peaks <- dbn025_peaks[, 1:length(peak_cols)]
  }
  colnames(dbn025_peaks) <- peak_cols
}, error = function(e) {
  stop(paste("Error reading file:", file2_path, "-", e$message))
})

# 2. Determine peak overlaps (DBN021 primary, find ALL DBN025 overlaps)
Overlap_Gap <- 0 # Only standard overlaps if 0
all_found_rows_list1 <- list()

for (i in 1:nrow(dbn021_peaks)) {
  peak021_chr   <- dbn021_peaks[i, "chr"]
  peak021_start <- dbn021_peaks[i, "start"]
  peak021_end   <- dbn021_peaks[i, "end"]
  
  partners_for_this_peak021 <- list()
  
  for (j in 1:nrow(dbn025_peaks)) {
    peak025_chr   <- dbn025_peaks[j, "chr"]
    peak025_start <- dbn025_peaks[j, "start"]
    peak025_end   <- dbn025_peaks[j, "end"]
    
    if (peak021_chr == peak025_chr) {
      is_standard_overlap <- (max(peak021_start, peak025_start) <= min(peak021_end, peak025_end))
      # With Overlap_Gap = 0, gap conditions are effectively false
      # is_left_gap_overlap <- FALSE
      # if (peak021_end < peak025_start) {
      #   if ((peak025_start - peak021_end) < Overlap_Gap) is_left_gap_overlap <- TRUE
      # }
      # is_right_gap_overlap <- FALSE
      # if (peak025_end < peak021_start) {
      #   if ((peak021_start - peak025_end) < Overlap_Gap) is_right_gap_overlap <- TRUE
      # }
      
      if (is_standard_overlap) { # || is_left_gap_overlap || is_right_gap_overlap) {
        partners_for_this_peak021 <- append(partners_for_this_peak021, list(
          data.frame(
            `021_Chr`   = peak021_chr,
            `021_Start` = peak021_start,
            `021_End`   = peak021_end,
            `025_Chr`   = peak025_chr,
            `025_Start` = peak025_start,
            `025_End`   = peak025_end,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        ))
      }
    }
  }
  
  if (length(partners_for_this_peak021) > 0) {
    all_found_rows_list1 <- c(all_found_rows_list1, partners_for_this_peak021)
  } else {
    all_found_rows_list1 <- append(all_found_rows_list1, list(
      data.frame(
        `021_Chr`   = peak021_chr,
        `021_Start` = peak021_start,
        `021_End`   = peak021_end,
        `025_Chr`   = NA_character_,
        `025_Start` = NA_integer_,
        `025_End`   = NA_integer_,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    ))
  }
}
results_df1 <- do.call(rbind, all_found_rows_list1)

# 3. Filter and Export results table for DBN021 vs DBN025
standard_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
results_df1_filtered <- results_df1[results_df1$"021_Chr" %in% standard_chromosomes, ]

output_file_path1 <- "../public/ChIPseq/ChIP_seq_021_025_overlap.csv"
tryCatch({
  write.csv(results_df1_filtered, file = output_file_path1, row.names = FALSE, quote = FALSE, na = "NA")
  print(paste("Filtered results (021 vs 025) successfully exported to:", output_file_path1))
  print(paste("Original number of DBN021 peaks considered (all chrom):", nrow(dbn021_peaks)))
  print(paste("Number of rows in intermediate file (021 vs 025, std chrom):", nrow(results_df1_filtered)))
}, error = function(e) {
  stop(paste("Error writing CSV file:", output_file_path1, "-", e$message))
})

################################################################################
# Script Part 2: Calculate percentage for DBN021 vs DBN025
################################################################################
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

input_file1_perc <- "../public/ChIPseq/ChIP_seq_021_025_overlap.csv"
output_file1_perc <- "../public/ChIPseq/ChIP_seq_021_025_overlap_percentage.csv"

tryCatch({
  chip_data1 <- read.csv(input_file1_perc, stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA")
  cat(paste("Successfully imported data from:", input_file1_perc, "for 021 vs 025 percentage calculation\n"))
}, error = function(e) {
  stop(paste("Error importing file:", input_file1_perc, "\n", e$message))
})

# Ensure `021_Chr` is not NA (should be guaranteed by previous filtering, but good practice)
chip_data1 <- chip_data1[!is.na(chip_data1$`021_Chr`) & chip_data1$`021_Chr` != "", ]

overlap_summary1 <- chip_data1 %>%
  group_by(`021_Chr`) %>%
  summarise(
    # Total unique DBN021 peaks for this chromosome in the input
    Num_021_Peaks_Total_In_Chr = n_distinct(paste(`021_Start`, `021_End`)),
    # Total number of pairwise overlaps on this chromosome
    Num_Total_Overlap_Pairs = sum(!is.na(`025_Chr`)),
    # Number of unique DBN021 peaks on this chromosome that participate in any overlap
    Num_021_Contributing_To_Pairs = n_distinct(paste(`021_Start`[!is.na(`025_Chr`)], `021_End`[!is.na(`025_Chr`)])),
    .groups = 'drop'
  ) %>%
  mutate(
    Overlap_Percentage = ifelse(Num_021_Peaks_Total_In_Chr > 0, round((Num_021_Contributing_To_Pairs / Num_021_Peaks_Total_In_Chr) * 100), 0)
  )

results_table1 <- overlap_summary1 %>%
  dplyr::select( # Explicitly use dplyr::select
    Chromosome = `021_Chr`,
    Overlap_Percentage,
    Num_021_Peaks = Num_021_Peaks_Total_In_Chr,
    Num_Overlapping_025_Peaks = Num_Total_Overlap_Pairs # This column now stores N_pairs
  )
results_table1$Overlap_Percentage <- as.integer(results_table1$Overlap_Percentage)

cat("\nCalculated overlap percentage table (021 vs 025):\n")
print(results_table1)

tryCatch({
  write.csv(results_table1, output_file1_perc, row.names = FALSE, quote = TRUE)
  cat(paste("\nSuccessfully exported results (021 vs 025 percentage) to:", output_file1_perc, "\n"))
}, error = function(e) {
  stop(paste("Error exporting file:", output_file1_perc, "\n", e$message))
})


################################################################################
# Script Part 3: DBN025 peaks as primary, find all overlapping DBN021 peaks
################################################################################
# DBN021 and DBN025 peaks are already loaded.

all_found_rows_list2 <- list()

for (i in 1:nrow(dbn025_peaks)) {
  peak025_chr   <- dbn025_peaks[i, "chr"]
  peak025_start <- dbn025_peaks[i, "start"]
  peak025_end   <- dbn025_peaks[i, "end"]
  
  partners_for_this_peak025 <- list()
  
  for (j in 1:nrow(dbn021_peaks)) {
    peak021_chr   <- dbn021_peaks[j, "chr"]
    peak021_start <- dbn021_peaks[j, "start"]
    peak021_end   <- dbn021_peaks[j, "end"]
    
    if (peak025_chr == peak021_chr) {
      is_standard_overlap <- (max(peak025_start, peak021_start) <= min(peak025_end, peak021_end))
      # With Overlap_Gap = 0, gap conditions are effectively false
      
      if (is_standard_overlap) {
        partners_for_this_peak025 <- append(partners_for_this_peak025, list(
          data.frame(
            `025_Chr`   = peak025_chr,
            `025_Start` = peak025_start,
            `025_End`   = peak025_end,
            `021_Chr`   = peak021_chr,
            `021_Start` = peak021_start,
            `021_End`   = peak021_end,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
        ))
      }
    }
  }
  
  if (length(partners_for_this_peak025) > 0) {
    all_found_rows_list2 <- c(all_found_rows_list2, partners_for_this_peak025)
  } else {
    all_found_rows_list2 <- append(all_found_rows_list2, list(
      data.frame(
        `025_Chr`   = peak025_chr,
        `025_Start` = peak025_start,
        `025_End`   = peak025_end,
        `021_Chr`   = NA_character_,
        `021_Start` = NA_integer_,
        `021_End`   = NA_integer_,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    ))
  }
}
results_df2 <- do.call(rbind, all_found_rows_list2)

# Filter and Export results table for DBN025 vs DBN021
results_df2_filtered <- results_df2[results_df2$"025_Chr" %in% standard_chromosomes, ]

output_file_path2 <- "../public/ChIPseq/ChIP_seq_025_021_overlap.csv"
tryCatch({
  write.csv(results_df2_filtered, file = output_file_path2, row.names = FALSE, quote = FALSE, na = "NA")
  print(paste("Filtered results (025 vs 021) successfully exported to:", output_file_path2))
  print(paste("Original number of DBN025 peaks considered (all chrom):", nrow(dbn025_peaks)))
  print(paste("Number of rows in intermediate file (025 vs 021, std chrom):", nrow(results_df2_filtered)))
}, error = function(e) {
  stop(paste("Error writing CSV file:", output_file_path2, "-", e$message))
})


################################################################################
# Script Part 4: Calculate percentage for DBN025 vs DBN021
################################################################################
# dplyr is already loaded

input_file2_perc <- "../public/ChIPseq/ChIP_seq_025_021_overlap.csv"
output_file2_perc <- "../public/ChIPseq/ChIP_seq_025_021_overlap_percentage.csv"

tryCatch({
  chip_data2 <- read.csv(input_file2_perc, stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA")
  cat(paste("Successfully imported data from:", input_file2_perc, "for 025 vs 021 percentage calculation\n"))
}, error = function(e) {
  stop(paste("Error importing file:", input_file2_perc, "\n", e$message))
})

chip_data2 <- chip_data2[!is.na(chip_data2$`025_Chr`) & chip_data2$`025_Chr` != "", ]

overlap_summary2 <- chip_data2 %>%
  group_by(`025_Chr`) %>%
  summarise(
    Num_025_Peaks_Total_In_Chr = n_distinct(paste(`025_Start`, `025_End`)),
    Num_Total_Overlap_Pairs = sum(!is.na(`021_Chr`)),
    Num_025_Contributing_To_Pairs = n_distinct(paste(`025_Start`[!is.na(`021_Chr`)], `025_End`[!is.na(`021_Chr`)])),
    .groups = 'drop'
  ) %>%
  mutate(
    Overlap_Percentage = ifelse(Num_025_Peaks_Total_In_Chr > 0, round((Num_025_Contributing_To_Pairs / Num_025_Peaks_Total_In_Chr) * 100), 0)
  )

results_table2 <- overlap_summary2 %>%
  dplyr::select( # Explicitly use dplyr::select
    Chromosome = `025_Chr`,
    Overlap_Percentage,
    Num_025_Peaks = Num_025_Peaks_Total_In_Chr,
    Num_Overlapping_021_Peaks = Num_Total_Overlap_Pairs # This column now stores N_pairs
  )
results_table2$Overlap_Percentage <- as.integer(results_table2$Overlap_Percentage)

cat("\nCalculated overlap percentage table (025 vs 021):\n")
print(results_table2)

tryCatch({
  write.csv(results_table2, output_file2_perc, row.names = FALSE, quote = TRUE)
  cat(paste("\nSuccessfully exported results (025 vs 021 percentage) to:", output_file2_perc, "\n"))
}, error = function(e) {
  stop(paste("Error exporting file:", output_file2_perc, "\n", e$message))
})


################################################################################

# R Script for processing two sets of ChIP-seq data against RNA-seq data

# 0. Define global variable for distance
ChIP_UpDEG_Distance <- 353000

# Load necessary libraries at the beginning
# Install dplyr if you haven't already: install.packages("dplyr")
# Note: dplyr is already loaded earlier, but this ensures it's loaded if this part is run separately.
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)


################################################################################
# PROCESSING PIPELINE 1: ChIP_seq_025_021_overlap.csv
################################################################################

print("# STARTING PROCESSING PIPELINE 1: ChIP_seq_025_021_overlap.csv                 #")

# Script 1.1: Overlap detection for ChIP_seq_025_021_overlap.csv
# ------------------------------------------------------------------------------
print("--- Script 1.1: Starting Overlap Detection for ChIP_seq_025_021_overlap.csv ---")
# 1. Import files
chip_file_path_025 <- "../public/ChIPseq/ChIP_seq_025_021_overlap.csv"
rna_file_path <- "../public/RNAseq/UpDEG_category.csv" # Common RNA-seq file

# Print message for user
print(paste("Attempting to read ChIP-seq data (025) from:", chip_file_path_025))
print(paste("Attempting to read RNA-seq data from:", rna_file_path))

# Check if files exist
if (!file.exists(chip_file_path_025)) {
  stop(paste("ChIP-seq file (025) not found at:", chip_file_path_025,
             "Please check the path. Current working directory is:", getwd()))
}
if (!file.exists(rna_file_path)) {
  stop(paste("RNA-seq file not found at:", rna_file_path,
             "Please check the path. Current working directory is:", getwd()))
}

# Read the CSV files
chip_data_025 <- read.csv(chip_file_path_025, stringsAsFactors = FALSE, na.strings = "NA", check.names = FALSE)
rna_data <- read.csv(rna_file_path, stringsAsFactors = FALSE, na.strings = "NA", check.names = FALSE)

print("Files for pipeline 1 (025) imported successfully.")
print("First few rows of ChIP data (025):")
print(head(chip_data_025))
print("First few rows of RNA data (common):")
print(head(rna_data))

# Ensure coordinate columns are numeric.
cols_to_numeric_chip <- c("025_Start", "025_End", "021_Start", "021_End")
for (col_name in cols_to_numeric_chip) {
  if (col_name %in% colnames(chip_data_025)) {
    chip_data_025[[col_name]] <- as.numeric(as.character(chip_data_025[[col_name]]))
  } else {
    print(paste("Warning: Column", col_name, "not found in ChIP data (025)."))
  }
}

cols_to_numeric_rna <- c("hg38_start", "hg38_end", "023_24hr_Up", "025_24hr_Up",
                         "021_10hr_Up", "021_24hr_Up", "021_48hr_Up")
for (col_name in cols_to_numeric_rna) {
  if (col_name %in% colnames(rna_data)) {
    rna_data[[col_name]] <- as.numeric(as.character(rna_data[[col_name]]))
  } else {
    print(paste("Warning: Column", col_name, "not found in RNA data."))
  }
}
print("Numeric conversions applied for pipeline 1 (025).")

# 2. Prepare for output
output_column_names <- c(
  "025_Chr", "025_Start", "025_End", "021_Chr", "021_Start", "021_End",
  "transcript_id", "Symbol", "hg38_chr", "hg38_start", "hg38_end", "hg38_strand",
  "023_24hr_Up", "025_24hr_Up", "021_10hr_Up", "021_24hr_Up", "021_48hr_Up",
  "Category - SRA specific", "Category - other"
)
results_list_025 <- list()

print("Starting overlap detection for pipeline 1 (025)...")
# 3. Find overlaps
for (i in 1:nrow(chip_data_025)) {
  chip_row <- chip_data_025[i, , drop = FALSE]
  
  # Extract ChIP peak info (PRIMARY: 025 for this pipeline)
  chip_chr_val <- chip_row[["025_Chr"]]
  chip_start_val <- chip_row[["025_Start"]]
  chip_end_val <- chip_row[["025_End"]]
  
  can_process_chip_peak <- !is.na(chip_chr_val) && !is.na(chip_start_val) && !is.na(chip_end_val)
  found_overlap_for_this_chip_peak <- FALSE
  
  if (can_process_chip_peak) {
    cleaned_chip_chr <- sub("^chr", "", chip_chr_val)
    
    for (j in 1:nrow(rna_data)) {
      rna_row <- rna_data[j, , drop = FALSE]
      
      rna_chr_val <- rna_row[["hg38_chr"]]
      rna_start_val <- rna_row[["hg38_start"]]
      rna_end_val <- rna_row[["hg38_end"]]
      
      if (is.na(rna_chr_val) || is.na(rna_start_val) || is.na(rna_end_val)) {
        next
      }
      
      cleaned_rna_chr <- as.character(rna_chr_val) # Ensure character comparison
      is_overlap <- FALSE
      if (cleaned_chip_chr == cleaned_rna_chr) {
        extended_chip_start <- chip_start_val - ChIP_UpDEG_Distance
        extended_chip_end <- chip_end_val + ChIP_UpDEG_Distance
        if (extended_chip_start <= rna_end_val && extended_chip_end >= rna_start_val) {
          is_overlap <- TRUE
        }
      }
      
      if (is_overlap) {
        current_result_df <- data.frame(matrix(NA, nrow = 1, ncol = length(output_column_names)))
        colnames(current_result_df) <- output_column_names
        
        for(col_n in colnames(chip_row)){
          if(col_n %in% output_column_names){
            current_result_df[1, col_n] <- chip_row[[col_n]]
          }
        }
        for(col_n in colnames(rna_row)){
          if(col_n %in% output_column_names){
            current_result_df[1, col_n] <- rna_row[[col_n]]
          }
        }
        results_list_025[[length(results_list_025) + 1]] <- current_result_df
        found_overlap_for_this_chip_peak <- TRUE
      }
    }
  }
  
  if (!found_overlap_for_this_chip_peak) {
    current_result_df <- data.frame(matrix(NA, nrow = 1, ncol = length(output_column_names)))
    colnames(current_result_df) <- output_column_names
    for(col_n in colnames(chip_row)){
      if(col_n %in% output_column_names){
        current_result_df[1, col_n] <- chip_row[[col_n]]
      }
    }
    results_list_025[[length(results_list_025) + 1]] <- current_result_df
  }
  
  if (i %% 100 == 0) {
    print(paste("Processed", i, "ChIP peaks (025) out of", nrow(chip_data_025)))
  }
}
print("Overlap detection finished for pipeline 1 (025).")

if (length(results_list_025) > 0) {
  final_output_table_025 <- do.call(rbind, results_list_025)
  all_numeric_cols_in_output <- unique(c(names(which(sapply(chip_data_025, is.numeric))), 
                                         names(which(sapply(rna_data, is.numeric))),
                                         cols_to_numeric_chip, cols_to_numeric_rna)) # Ensure defined numeric cols are included
  all_numeric_cols_in_output <- intersect(all_numeric_cols_in_output, colnames(final_output_table_025))
  
  
  for(col in all_numeric_cols_in_output){
    if(col %in% colnames(final_output_table_025)){
      suppressWarnings(final_output_table_025[[col]] <- as.numeric(as.character(final_output_table_025[[col]])))
    }
  }
  print("Final table assembled for pipeline 1 (025) and numeric types re-applied.")
  print("First few rows of the final output table (025):")
  print(head(final_output_table_025))
} else {
  final_output_table_025 <- data.frame(matrix(ncol = length(output_column_names), nrow = 0))
  colnames(final_output_table_025) <- output_column_names
  print("No results to combine for pipeline 1 (025). An empty table with specified columns will be written.")
}

output_dir <- "../public/Multiomics/"
if (!dir.exists(output_dir)) {
  print(paste("Output directory", output_dir, "does not exist. Attempting to create it."))
  dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)
}
output_file_path_025_overlap <- file.path(output_dir, "ChIPseq_025_RNAseq_Overlap.csv")
write.csv(final_output_table_025, output_file_path_025_overlap, row.names = FALSE, na = "NA")
print(paste("Output table for pipeline 1 (025) saved to:", output_file_path_025_overlap))
print("--- Script 1.1: Finished Overlap Detection for ChIP_seq_025_021_overlap.csv ---")
print("") # newline for separation

# Script 1.2: Filtering for ChIPseq_025_RNAseq_Overlap.csv
# ------------------------------------------------------------------------------
print("--- Script 1.2: Starting Filtering for ChIPseq_025_RNAseq_Overlap.csv ---")
input_file_path_025_filter <- "../public/Multiomics/ChIPseq_025_RNAseq_Overlap.csv"
if (!file.exists(input_file_path_025_filter)) {
  stop(paste("Error: Input file for filtering (025) not found at", input_file_path_025_filter))
}
tryCatch({
  chip_rna_data_025 <- read.csv(input_file_path_025_filter, stringsAsFactors = FALSE, na.strings = "NA", check.names = FALSE)
  print(paste("Successfully imported data for filtering (025) from:", input_file_path_025_filter))
  print("Original column names (025):")
  print(colnames(chip_rna_data_025))
}, error = function(e) {
  stop(paste("Error reading CSV file for filtering (025):", e$message))
})

chip_rna_data_copy_025 <- chip_rna_data_025
if (!"transcript_id" %in% colnames(chip_rna_data_copy_025)) {
  stop("Error: 'transcript_id' column not found in 025 data. Please check column names.")
}
filtered_data_025 <- chip_rna_data_copy_025[!is.na(chip_rna_data_copy_025$transcript_id), ]
print(paste("Rows before filtering NA in 'transcript_id' (025):", nrow(chip_rna_data_copy_025)))
print(paste("Rows after filtering NA in 'transcript_id' (025):", nrow(filtered_data_025)))

output_file_path_025_filtered <- "../public/Multiomics/ChIPseq_025_RNAseq_Overlap_Filtered.csv"
output_dir_025_filtered <- dirname(output_file_path_025_filtered)
if (!dir.exists(output_dir_025_filtered)) {
  dir.create(output_dir_025_filtered, recursive = TRUE)
}
tryCatch({
  write.csv(filtered_data_025, file = output_file_path_025_filtered, row.names = FALSE, na = "NA")
  print(paste("Successfully exported filtered data (025) to:", output_file_path_025_filtered))
}, error = function(e) {
  stop(paste("Error writing filtered CSV file (025):", e$message))
})
print("--- Script 1.2: Finished Filtering for ChIPseq_025_RNAseq_Overlap.csv ---")
print("") # newline

# Script 1.3: Sorting for ChIPseq_025_RNAseq_Overlap_Filtered.csv
# ------------------------------------------------------------------------------
print("--- Script 1.3: Starting Sorting for ChIPseq_025_RNAseq_Overlap_Filtered.csv ---")
input_file_path_025_sort <- "../public/Multiomics/ChIPseq_025_RNAseq_Overlap_Filtered.csv"
output_file_path_025_sorted <- "../public/Multiomics/ChIPseq_025_RNAseq_Overlap_Filtered_Sorted.csv"

tryCatch({
  original_data_025_sort <- read.csv(input_file_path_025_sort, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  print("Successfully imported filtered data (025) for sorting.")
  data_to_sort_025 <- original_data_025_sort
  if ("transcript_id" %in% colnames(data_to_sort_025)) {
    sorted_data_025 <- data_to_sort_025[order(data_to_sort_025$transcript_id), ]
    print("Data (025) sorted successfully by 'transcript_id'.")
    tryCatch({
      write.csv(sorted_data_025, file = output_file_path_025_sorted, row.names = FALSE, quote = TRUE, na = "NA")
      print(paste("Sorted data (025) successfully exported to:", output_file_path_025_sorted))
    }, error = function(e_export) {
      print(paste("Error exporting the sorted CSV file (025):", e_export$message))
    })
  } else {
    print("Error: 'transcript_id' column not found in 025 data for sorting.")
  }
}, error = function(e_import) {
  print(paste("Error importing CSV file for sorting (025):", e_import$message))
})
print("--- Script 1.3: Finished Sorting for ChIPseq_025_RNAseq_Overlap_Filtered.csv ---")
print("") # newline

# Script 1.4: Condensing for ChIPseq_025_RNAseq_Overlap_Filtered_Sorted.csv
# ------------------------------------------------------------------------------
print("--- Script 1.4: Starting Condensing for ChIPseq_025_RNAseq_Overlap_Filtered_Sorted.csv ---")
input_file_025_condense <- "../public/Multiomics/ChIPseq_025_RNAseq_Overlap_Filtered_Sorted.csv"
tryCatch({
  data_025_condense <- read.csv(input_file_025_condense, stringsAsFactors = FALSE, na.strings = "NA", check.names = FALSE)
  cat("Successfully imported sorted data (025) for condensing:", input_file_025_condense, "\n")
}, error = function(e) {
  stop(paste("Error importing file for condensing (025):", e$message))
})

condensed_data_025 <- data_025_condense %>%
  group_by(transcript_id) %>%
  summarise(
    hg38_chr = first(hg38_chr),
    hg38_start = first(hg38_start),
    hg38_end = first(hg38_end),
    Symbol = first(Symbol),
    `Category - SRA specific` = first(`Category - SRA specific`),
    `Category - other` = first(`Category - other`),
    Num_025_Peak = n(), # Primary peak for this pipeline
    Num_021_Peak = sum(!is.na(`021_Chr`)), # Secondary peak, use backticks
    .groups = 'drop' # Added to avoid grouped tbl warning with dplyr >= 1.0.0
  ) %>%
  ungroup() # Ensure it's fully ungrouped if summarise didn't do it (older dplyr)

condensed_data_025 <- condensed_data_025 %>%
  dplyr::select( # Explicitly use dplyr::select
    hg38_chr, hg38_start, hg38_end, transcript_id, Symbol,
    `Category - SRA specific`, `Category - other`,
    Num_025_Peak, Num_021_Peak
  )
cat("Successfully condensed data (025) based on transcript_id.\n")

condensed_data_025 <- condensed_data_025 %>%
  mutate(
    chr_numeric = suppressWarnings(as.integer(gsub("chr", "", hg38_chr))),
    chr_char = gsub("[0-9]", "", hg38_chr) # Extracts non-numeric part like X, Y, M
  ) %>%
  arrange(chr_numeric, chr_char, hg38_start) %>%
  dplyr::select(-chr_numeric, -chr_char) # Explicitly use dplyr::select
cat("Successfully sorted the condensed data (025).\n")

output_file_025_condensed <- "../public/Multiomics/ChIPseq_025_RNAseq_Overlap_Filtered_Sorted_Condensed.csv"
tryCatch({
  write.csv(condensed_data_025, output_file_025_condensed, row.names = FALSE, quote = TRUE, na = "NA")
  cat("Successfully exported condensed data (025) to:", output_file_025_condensed, "\n")
}, error = function(e) {
  stop(paste("Error exporting condensed file (025):", e$message))
})
cat("\nFirst few rows of the condensed and sorted table (025):\n")
print(head(condensed_data_025))
cat("\nDimensions of the condensed table (025) (rows, columns):\n")
print(dim(condensed_data_025))
print("--- Script 1.4: Finished Condensing for ChIPseq_025_RNAseq_Overlap_Filtered_Sorted.csv ---")
print("# FINISHED PROCESSING PIPELINE 1: ChIP_seq_025_021_overlap.csv                #")
print("\n\n") # Extra newlines for separation


################################################################################
# PROCESSING PIPELINE 2: ChIP_seq_021_025_overlap.csv
################################################################################

print("# STARTING PROCESSING PIPELINE 2: ChIP_seq_021_025_overlap.csv                 #")

# Script 2.1: Overlap detection for ChIP_seq_021_025_overlap.csv
# ------------------------------------------------------------------------------
print("--- Script 2.1: Starting Overlap Detection for ChIP_seq_021_025_overlap.csv ---")
# 1. Import files
chip_file_path_021 <- "../public/ChIPseq/ChIP_seq_021_025_overlap.csv"
# rna_file_path is already defined and rna_data is already loaded

# Print message for user
print(paste("Attempting to read ChIP-seq data (021) from:", chip_file_path_021))
# RNA data already loaded

# Check if files exist
if (!file.exists(chip_file_path_021)) {
  stop(paste("ChIP-seq file (021) not found at:", chip_file_path_021,
             "Please check the path. Current working directory is:", getwd()))
}

# Read the CSV files
chip_data_021 <- read.csv(chip_file_path_021, stringsAsFactors = FALSE, na.strings = "NA", check.names = FALSE)
# rna_data is already loaded and processed (numeric conversions done)

print("Files for pipeline 2 (021) imported successfully.")
print("First few rows of ChIP data (021):")
print(head(chip_data_021))
# RNA data already shown

# Ensure coordinate columns are numeric for the new ChIP data.
# cols_to_numeric_chip is already defined c("025_Start", "025_End", "021_Start", "021_End")
for (col_name in cols_to_numeric_chip) {
  if (col_name %in% colnames(chip_data_021)) {
    chip_data_021[[col_name]] <- as.numeric(as.character(chip_data_021[[col_name]]))
  } else {
    print(paste("Warning: Column", col_name, "not found in ChIP data (021)."))
  }
}
# Numeric conversion for rna_data already done.
print("Numeric conversions applied for ChIP data in pipeline 2 (021).")

# 2. Prepare for output
# output_column_names is already defined
results_list_021 <- list()

print("Starting overlap detection for pipeline 2 (021)...")
# 3. Find overlaps
for (i in 1:nrow(chip_data_021)) {
  chip_row <- chip_data_021[i, , drop = FALSE]
  
  # Extract ChIP peak info (PRIMARY: 021 for this pipeline)
  chip_chr_val <- chip_row[["021_Chr"]] # Changed from "025_Chr"
  chip_start_val <- chip_row[["021_Start"]] # Changed from "025_Start"
  chip_end_val <- chip_row[["021_End"]]   # Changed from "025_End"
  
  can_process_chip_peak <- !is.na(chip_chr_val) && !is.na(chip_start_val) && !is.na(chip_end_val)
  found_overlap_for_this_chip_peak <- FALSE
  
  if (can_process_chip_peak) {
    cleaned_chip_chr <- sub("^chr", "", chip_chr_val)
    
    for (j in 1:nrow(rna_data)) {
      rna_row <- rna_data[j, , drop = FALSE]
      
      rna_chr_val <- rna_row[["hg38_chr"]]
      rna_start_val <- rna_row[["hg38_start"]]
      rna_end_val <- rna_row[["hg38_end"]]
      
      if (is.na(rna_chr_val) || is.na(rna_start_val) || is.na(rna_end_val)) {
        next
      }
      
      cleaned_rna_chr <- as.character(rna_chr_val) # Ensure character comparison
      is_overlap <- FALSE
      if (cleaned_chip_chr == cleaned_rna_chr) {
        extended_chip_start <- chip_start_val - ChIP_UpDEG_Distance
        extended_chip_end <- chip_end_val + ChIP_UpDEG_Distance
        if (extended_chip_start <= rna_end_val && extended_chip_end >= rna_start_val) {
          is_overlap <- TRUE
        }
      }
      
      if (is_overlap) {
        current_result_df <- data.frame(matrix(NA, nrow = 1, ncol = length(output_column_names)))
        colnames(current_result_df) <- output_column_names
        
        # Fill ChIP part (from chip_row of ChIP_seq_021_025_overlap.csv)
        for(col_n in colnames(chip_row)){
          if(col_n %in% output_column_names){
            current_result_df[1, col_n] <- chip_row[[col_n]]
          }
        }
        # Fill RNA part
        for(col_n in colnames(rna_row)){
          if(col_n %in% output_column_names){
            current_result_df[1, col_n] <- rna_row[[col_n]]
          }
        }
        results_list_021[[length(results_list_021) + 1]] <- current_result_df
        found_overlap_for_this_chip_peak <- TRUE
      }
    }
  }
  
  if (!found_overlap_for_this_chip_peak) {
    current_result_df <- data.frame(matrix(NA, nrow = 1, ncol = length(output_column_names)))
    colnames(current_result_df) <- output_column_names
    # Fill ChIP part (from chip_row of ChIP_seq_021_025_overlap.csv)
    for(col_n in colnames(chip_row)){
      if(col_n %in% output_column_names){
        current_result_df[1, col_n] <- chip_row[[col_n]]
      }
    }
    # RNA part remains NA
    results_list_021[[length(results_list_021) + 1]] <- current_result_df
  }
  
  if (i %% 100 == 0) {
    print(paste("Processed", i, "ChIP peaks (021) out of", nrow(chip_data_021)))
  }
}
print("Overlap detection finished for pipeline 2 (021).")

if (length(results_list_021) > 0) {
  final_output_table_021 <- do.call(rbind, results_list_021)
  # all_numeric_cols_in_output is already defined and should cover relevant columns
  # Re-check which columns should be numeric based on the inputs for this specific table
  all_numeric_cols_in_output_021 <- unique(c(names(which(sapply(chip_data_021, is.numeric))), 
                                             names(which(sapply(rna_data, is.numeric))),
                                             cols_to_numeric_chip, cols_to_numeric_rna))
  all_numeric_cols_in_output_021 <- intersect(all_numeric_cols_in_output_021, colnames(final_output_table_021))
  
  
  for(col in all_numeric_cols_in_output_021){ # Use the specific set for 021
    if(col %in% colnames(final_output_table_021)){
      suppressWarnings(final_output_table_021[[col]] <- as.numeric(as.character(final_output_table_021[[col]])))
    }
  }
  print("Final table assembled for pipeline 2 (021) and numeric types re-applied.")
  print("First few rows of the final output table (021):")
  print(head(final_output_table_021))
} else {
  final_output_table_021 <- data.frame(matrix(ncol = length(output_column_names), nrow = 0))
  colnames(final_output_table_021) <- output_column_names
  print("No results to combine for pipeline 2 (021). An empty table with specified columns will be written.")
}

# output_dir is already defined
output_file_path_021_overlap <- file.path(output_dir, "ChIPseq_021_RNAseq_Overlap.csv") # New output file name
write.csv(final_output_table_021, output_file_path_021_overlap, row.names = FALSE, na = "NA")
print(paste("Output table for pipeline 2 (021) saved to:", output_file_path_021_overlap))
print("--- Script 2.1: Finished Overlap Detection for ChIP_seq_021_025_overlap.csv ---")
print("") # newline

# Script 2.2: Filtering for ChIPseq_021_RNAseq_Overlap.csv
# ------------------------------------------------------------------------------
print("--- Script 2.2: Starting Filtering for ChIPseq_021_RNAseq_Overlap.csv ---")
input_file_path_021_filter <- "../public/Multiomics/ChIPseq_021_RNAseq_Overlap.csv" # New input file name
if (!file.exists(input_file_path_021_filter)) {
  stop(paste("Error: Input file for filtering (021) not found at", input_file_path_021_filter))
}
tryCatch({
  chip_rna_data_021 <- read.csv(input_file_path_021_filter, stringsAsFactors = FALSE, na.strings = "NA", check.names = FALSE)
  print(paste("Successfully imported data for filtering (021) from:", input_file_path_021_filter))
  print("Original column names (021):")
  print(colnames(chip_rna_data_021))
}, error = function(e) {
  stop(paste("Error reading CSV file for filtering (021):", e$message))
})

chip_rna_data_copy_021 <- chip_rna_data_021
if (!"transcript_id" %in% colnames(chip_rna_data_copy_021)) {
  stop("Error: 'transcript_id' column not found in 021 data. Please check column names.")
}
filtered_data_021 <- chip_rna_data_copy_021[!is.na(chip_rna_data_copy_021$transcript_id), ]
print(paste("Rows before filtering NA in 'transcript_id' (021):", nrow(chip_rna_data_copy_021)))
print(paste("Rows after filtering NA in 'transcript_id' (021):", nrow(filtered_data_021)))

output_file_path_021_filtered <- "../public/Multiomics/ChIPseq_021_RNAseq_Overlap_Filtered.csv" # New output file name
output_dir_021_filtered <- dirname(output_file_path_021_filtered)
if (!dir.exists(output_dir_021_filtered)) {
  dir.create(output_dir_021_filtered, recursive = TRUE)
}
tryCatch({
  write.csv(filtered_data_021, file = output_file_path_021_filtered, row.names = FALSE, na = "NA")
  print(paste("Successfully exported filtered data (021) to:", output_file_path_021_filtered))
}, error = function(e) {
  stop(paste("Error writing filtered CSV file (021):", e$message))
})
print("--- Script 2.2: Finished Filtering for ChIPseq_021_RNAseq_Overlap.csv ---")
print("") # newline

# Script 2.3: Sorting for ChIPseq_021_RNAseq_Overlap_Filtered.csv
# ------------------------------------------------------------------------------
print("--- Script 2.3: Starting Sorting for ChIPseq_021_RNAseq_Overlap_Filtered.csv ---")
input_file_path_021_sort <- "../public/Multiomics/ChIPseq_021_RNAseq_Overlap_Filtered.csv" # New input file name
output_file_path_021_sorted <- "../public/Multiomics/ChIPseq_021_RNAseq_Overlap_Filtered_Sorted.csv" # New output file name

tryCatch({
  original_data_021_sort <- read.csv(input_file_path_021_sort, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  print("Successfully imported filtered data (021) for sorting.")
  data_to_sort_021 <- original_data_021_sort
  if ("transcript_id" %in% colnames(data_to_sort_021)) {
    sorted_data_021 <- data_to_sort_021[order(data_to_sort_021$transcript_id), ]
    print("Data (021) sorted successfully by 'transcript_id'.")
    tryCatch({
      write.csv(sorted_data_021, file = output_file_path_021_sorted, row.names = FALSE, quote = TRUE, na = "NA")
      print(paste("Sorted data (021) successfully exported to:", output_file_path_021_sorted))
    }, error = function(e_export) {
      print(paste("Error exporting the sorted CSV file (021):", e_export$message))
    })
  } else {
    print("Error: 'transcript_id' column not found in 021 data for sorting.")
  }
}, error = function(e_import) {
  print(paste("Error importing CSV file for sorting (021):", e_import$message))
})
print("--- Script 2.3: Finished Sorting for ChIPseq_021_RNAseq_Overlap_Filtered.csv ---")
print("") # newline

# Script 2.4: Condensing for ChIPseq_021_RNAseq_Overlap_Filtered_Sorted.csv
# ------------------------------------------------------------------------------
print("--- Script 2.4: Starting Condensing for ChIPseq_021_RNAseq_Overlap_Filtered_Sorted.csv ---")
input_file_021_condense <- "../public/Multiomics/ChIPseq_021_RNAseq_Overlap_Filtered_Sorted.csv" # New input file name
tryCatch({
  data_021_condense <- read.csv(input_file_021_condense, stringsAsFactors = FALSE, na.strings = "NA", check.names = FALSE)
  cat("Successfully imported sorted data (021) for condensing:", input_file_021_condense, "\n")
}, error = function(e) {
  stop(paste("Error importing file for condensing (021):", e$message))
})

condensed_data_021 <- data_021_condense %>%
  group_by(transcript_id) %>%
  summarise(
    hg38_chr = first(hg38_chr),
    hg38_start = first(hg38_start),
    hg38_end = first(hg38_end),
    Symbol = first(Symbol),
    `Category - SRA specific` = first(`Category - SRA specific`),
    `Category - other` = first(`Category - other`),
    Num_025_Peak = sum(!is.na(`025_Chr`)), # Secondary peak for this pipeline, use backticks
    Num_021_Peak = n(),                   # Primary peak for this pipeline
    .groups = 'drop' # Added to avoid grouped tbl warning with dplyr >= 1.0.0
  ) %>%
  ungroup() # Ensure it's fully ungrouped

condensed_data_021 <- condensed_data_021 %>%
  dplyr::select( # Explicitly use dplyr::select. Ensure order matches the original script's select for 025, if Num_025_Peak and Num_021_Peak are swapped, adjust here
    hg38_chr, hg38_start, hg38_end, transcript_id, Symbol,
    `Category - SRA specific`, `Category - other`,
    Num_025_Peak, Num_021_Peak # Order maintained as per original output structure (Num_025_Peak then Num_021_Peak)
  )
cat("Successfully condensed data (021) based on transcript_id.\n")

condensed_data_021 <- condensed_data_021 %>%
  mutate(
    chr_numeric = suppressWarnings(as.integer(gsub("chr", "", hg38_chr))),
    chr_char = gsub("[0-9]", "", hg38_chr) # Extracts non-numeric part like X, Y, M
  ) %>%
  arrange(chr_numeric, chr_char, hg38_start) %>%
  dplyr::select(-chr_numeric, -chr_char) # Explicitly use dplyr::select
cat("Successfully sorted the condensed data (021).\n")

output_file_021_condensed <- "../public/Multiomics/ChIPseq_021_RNAseq_Overlap_Filtered_Sorted_Condensed.csv" # New output file name
tryCatch({
  write.csv(condensed_data_021, output_file_021_condensed, row.names = FALSE, quote = TRUE, na = "NA")
  cat("Successfully exported condensed data (021) to:", output_file_021_condensed, "\n")
}, error = function(e) {
  stop(paste("Error exporting condensed file (021):", e$message))
})
cat("\nFirst few rows of the condensed and sorted table (021):\n")
print(head(condensed_data_021))
cat("\nDimensions of the condensed table (021) (rows, columns):\n")
print(dim(condensed_data_021))
print("--- Script 2.4: Finished Condensing for ChIPseq_021_RNAseq_Overlap_Filtered_Sorted.csv ---")
print("# FINISHED PROCESSING PIPELINE 2: ChIP_seq_021_025_overlap.csv                 #")

print("\n\nAll script executions completed.")

################################################################################

################################################################################

# Title: Venn Diagram Generation for Multiomics Data
# Description: This script reads a CSV file with peak and transcript counts,
#              and generates a three-way, area-proportional Euler diagram.
# Author: Gemini
# Date: 2025-10-15

# --- 1. Installation and Library Loading ---

# Check if the 'eulerr' package is installed. If not, install it.
# This package is excellent for creating area-proportional Euler diagrams.
if (!requireNamespace("eulerr", quietly = TRUE)) {
  install.packages("eulerr")
}

# Load the eulerr library.
library(eulerr)

# --- 2. File Paths and Data Import ---

# Define the input and output file paths for clarity and easy modification.
inputFile <- "../public/Multiomics/SHK 2025_09_24 MCF7_SRA_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks_Cytobands_Venn Diagram.csv"
outputFile <- "../public/Multiomics/SHK 2025_09_24 MCF7_SRA_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks_Cytobands_Venn Diagram.jpg"

# Read the data from the CSV file into a data frame.
# The 'check.names=TRUE' default in read.csv will convert column names
# that start with numbers to valid R names (e.g., "X021_Peak_Count").
cat("Reading data from:", inputFile, "\n")
data <- read.csv(inputFile)

# --- 3. Data Preparation ---

# To create the diagram, we first need a logical matrix or data frame
# where each column represents a set and each row an element.
# The value is TRUE if the element is in the set, and FALSE otherwise.
logical_data <- data.frame(
  UpDEG = data$transcript_Count > 0,
  SRA = data$X021_Peak_Count > 0,
  `PCD-RFP` = data$X025_Peak_Count > 0
)

# --- 4. Euler Diagram Generation and Export ---

# Calculate the Euler diagram fit. The euler() function takes the logical
# data frame and computes the ideal positions and sizes for the shapes.
fit <- euler(logical_data)

# To export a high-quality JPG, we open a jpeg graphics device first.
jpeg(filename = outputFile, width = 800, height = 800, quality = 100, res = 120)

# Now, plot the fitted Euler diagram. The plot() function for eulerr objects
# has several options for customization.
plot(fit,
     # --- Styling Options ---
     # Display the counts within the regions.
     quantities = TRUE,
     
     # Set the labels for the sets.
     labels = c("UpDEG", "SRA", "PCD-RFP"),
     
     # Set the fill colors with transparency.
     fills = list(
       fill = c("gray", "red", "blue"),
       alpha = 0.5
     ),
     
     # Set the font size for labels and quantities.
     cex = 1.5
)

# Close the jpeg graphics device to finalize writing the file.
dev.off()

# Print a confirmation message to the console.
cat("Successfully generated and saved Euler diagram to:", outputFile, "\n")

