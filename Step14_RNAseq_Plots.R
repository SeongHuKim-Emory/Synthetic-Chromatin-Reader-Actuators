setwd("C:/Users/ski2959/Desktop/NIH R21/Codes/src")

################################################################################

################################################################################

# R Script for RNA-seq Data Processing and Averaging

# --- 1. Import Files ---

# This section imports all the necessary CSV files.
# The primary file for this script is 'gexp_counts_symbol_normalized.csv'.
# The other files are imported as requested, although they are not used in the subsequent calculations.

cat("Importing data files...\n")

# Import the main gene expression counts table
gexp_counts <- read.csv("../public/RNAseq/gexp_counts_symbol_normalized.csv")

# Import the DEG category file
updeg_category <- read.csv("../public/RNAseq/UpDEG_category.csv")

# Import the individual DEG files
mcf7_021_10_updeg <- read.csv("../public/RNAseq/MCF7_021_10_UpDEG.csv")
mcf7_021_24_updeg <- read.csv("../public/RNAseq/MCF7_021_24_UpDEG.csv")
mcf7_021_48_updeg <- read.csv("../public/RNAseq/MCF7_021_48_UpDEG.csv")
mcf7_023_24_updeg <- read.csv("../public/RNAseq/MCF7_023_24_UpDEG.csv")
mcf7_025_24_updeg <- read.csv("../public/RNAseq/MCF7_025_24_UpDEG.csv")

cat("All files imported successfully.\n")


# --- 2 & 3. Create New Table and Calculate Mean of Sample Groups ---

# This section creates a new data frame to store the results.
# It then calculates the row-wise average for each sample group and populates the new data frame.

cat("Calculating mean expression for sample groups...\n")

# Initialize the new data frame, starting with the identifier columns from the original table.
# This fulfills the "Duplicate" step, as we are creating the new table structure.
gexp_counts_symbol_normalized_mean <- gexp_counts[, c("transcript_id", "symbol")]

# Define the unique prefixes for each sample group to identify their replicate columns.
sample_group_prefixes <- c(
  "MCF7_023_0",
  "MCF7_023_24",
  "MCF7_025_0",
  "MCF7_025_24",
  "MCF7_021_0",
  "MCF7_021_10",
  "MCF7_021_24",
  "MCF7_021_48"
)

# Loop through each defined sample group prefix.
for (prefix in sample_group_prefixes) {
  
  # Find all column names in the original data that start with the current prefix.
  # The `grep` function with `^` searches for the pattern at the beginning of the string.
  replicate_cols <- grep(paste0("^", prefix), names(gexp_counts), value = TRUE)
  
  # Check if any columns were found for the prefix
  if (length(replicate_cols) > 0) {
    # Calculate the mean across the identified replicate columns for each row (gene).
    # `rowMeans` is a highly efficient function for this task.
    mean_values <- rowMeans(gexp_counts[, replicate_cols], na.rm = TRUE)
    
    # Create a new, descriptive column name for the calculated means (e.g., "MCF7_023_0_Mean").
    new_col_name <- paste0(prefix, "_Mean")
    
    # Add the new column with the mean values to our results data frame.
    gexp_counts_symbol_normalized_mean[[new_col_name]] <- mean_values
  }
}

cat("Mean calculations are complete.\n")


# --- 4. Export the Result Table ---

# This section saves the newly created data frame with the averaged values to a CSV file.

cat("Exporting the final table...\n")

# Write the data frame to the specified file path.
# `row.names = FALSE` prevents R from writing its internal row numbers into the file.
write.csv(gexp_counts_symbol_normalized_mean, 
          file = "../public/RNAseq/gexp_counts_symbol_normalized_mean.csv", 
          row.names = FALSE)

cat("Process finished. The file 'gexp_counts_symbol_normalized_mean.csv' has been created.\n")

################################################################################

################################################################################

# RNAseq Data Merging and Annotation Script
# This script imports gene category and expression data, merges them based on
# transcript ID, and exports the combined dataset.

# Install and load the 'dplyr' package if you don't have it
# dplyr is a powerful library for data manipulation.
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

# --- Step 1: Import all required files ---
# Define file paths for clarity and easy modification
updeg_category_file <- "../public/RNAseq/UpDEG_category.csv"
gexp_counts_file <- "../public/RNAseq/gexp_counts_symbol_normalized_mean.csv"
output_file <- "../public/RNAseq/UpDEG_category_NormCount.csv"

# Reading the primary data files
message("Reading input files...")
updeg_category <- read.csv(updeg_category_file)
gexp_counts <- read.csv(gexp_counts_file)

# The user also requested importing the following files.
# While they are not used in the final merge step described, we load them as requested.
mcf7_021_10_deg <- read.csv("../public/RNAseq/MCF7_021_10_UpDEG.csv")
mcf7_021_24_deg <- read.csv("../public/RNAseq/MCF7_021_24_UpDEG.csv")
mcf7_021_48_deg <- read.csv("../public/RNAseq/MCF7_021_48_UpDEG.csv")
mcf7_023_24_deg <- read.csv("../public/RNAseq/MCF7_023_24_UpDEG.csv")
mcf7_025_24_deg <- read.csv("../public/RNAseq/MCF7_025_24_UpDEG.csv")
message("All files imported successfully.")

# --- Step 2: Duplicate the category table ---
# Create a new table that will hold the merged data.
# This keeps the original data frame unchanged.
UpDEG_category_NormCount_log2FC <- updeg_category

# --- Steps 3 & 4: Add and populate normalized count columns ---
# We will merge the duplicated table with the gene expression data.
# A 'left_join' is used to ensure all original genes from the category list are kept.
# The merge is performed by matching the 'transcript_id' column in both tables.
# This single operation effectively adds and fills the new columns.

message("Merging data frames based on 'transcript_id'...")

# To prevent a duplicate 'symbol' column in our final table, we will select
# all columns from the gexp_counts table EXCEPT for its 'symbol' column before joining.
gexp_counts_subset <- select(gexp_counts, -symbol)

# Perform the left join
final_data <- left_join(
  UpDEG_category_NormCount_log2FC,
  gexp_counts_subset,
  by = "transcript_id"
)

message("Data merged successfully.")

# --- Step 5: Export the final table to a CSV file ---
# The resulting data frame is saved to the specified output path.
# 'row.names = FALSE' prevents R from writing an extra column for row numbers.
message(paste("Exporting final table to:", output_file))
write.csv(
  final_data,
  file = output_file,
  row.names = FALSE
)
message("Script finished successfully.")

# You can uncomment the line below to see a preview of the first few rows
# of the final, merged data frame in your R console.
# print(head(final_data))

################################################################################

################################################################################

# RNA-seq Data Processing Script
# This script imports normalized count data and several differential expression
# (DEG) files, merges them based on transcript_id, and exports the combined data.

# Install and load the dplyr package if it's not already installed.
# dplyr is used for efficient data manipulation, particularly the 'left_join' function.
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

# 1. Import Files
# Define file paths for clarity and easier management.
main_file_path <- "../public/RNAseq/UpDEG_category_NormCount.csv"
deg_files_paths <- c(
  "023_24" = "../public/RNAseq/MCF7_023_24_DEG.csv",
  "025_24" = "../public/RNAseq/MCF7_025_24_DEG.csv",
  "021_10" = "../public/RNAseq/MCF7_021_10_DEG.csv",
  "021_24" = "../public/RNAseq/MCF7_021_24_DEG.csv",
  "021_48" = "../public/RNAseq/MCF7_021_48_DEG.csv"
)

# Read the main data file.
updeg_data <- read.csv(main_file_path)

# Read all the DEG files into a list of data frames.
deg_data_list <- lapply(deg_files_paths, read.csv)

# 2. Duplicate the main data frame
# This creates a new table where we will add the log2FoldChange values.
updeg_data_log2fc <- updeg_data

# 3 & 5. Add and Populate log2FC Columns (Combined for efficiency)
# We will iterate through each DEG file, select the necessary columns,
# rename the log2FoldChange column, and then merge it into our main data frame.

for (name in names(deg_data_list)) {
  # Create the new column name for the log2FC value
  log2fc_col_name <- paste0(name, "_log2FC")
  
  # Select only the transcript_id and log2FoldChange from the current DEG file.
  # Rename 'log2FoldChange' to its specific name (e.g., '023_24_log2FC').
  # The `!!` (bang-bang) operator is used to evaluate the variable `log2fc_col_name`
  # as a column name within the select function.
  deg_subset <- deg_data_list[[name]] %>%
    select(transcript_id, !!log2fc_col_name := log2FoldChange)
  
  # Join the subset with our main data frame by 'transcript_id'.
  # This adds the new log2FC column and fills it with the corresponding values.
  updeg_data_log2fc <- left_join(updeg_data_log2fc, deg_subset, by = "transcript_id")
}


# 4. Remove "X" from column names
# The `gsub` function finds all occurrences of "X" at the beginning of a column name
# that is followed by a number and replaces it with an empty string.
colnames(updeg_data_log2fc) <- gsub("^X([0-9])", "\\1", colnames(updeg_data_log2fc))


# 6. Export the final data frame to a new CSV file
# `row.names = FALSE` prevents R from writing the default row index as a column.
output_file_path <- "../public/RNAseq/UpDEG_category_NormCount_log2FC.csv"
write.csv(updeg_data_log2fc, output_file_path, row.names = FALSE)

# Print a confirmation message to the console.
cat("Processing complete. The file has been saved to:", output_file_path, "\n")

# Display the first few rows and column names of the final data frame to verify.
cat("\nFirst 6 rows of the final data frame:\n")
print(head(updeg_data_log2fc))

cat("\nColumn names of the final data frame:\n")
print(colnames(updeg_data_log2fc))

################################################################################

################################################################################

# RNA Sequence Analysis and Visualization Script
# This script imports gene expression data, sorts it by specific categories,
# and generates box and whisker plots for each category, saving them as JPG files.

# --- 1. Load Required Libraries ---
# Ensure you have these packages installed. If not, run:
# install.packages(c("dplyr", "tidyr", "ggplot2"))

library(dplyr)    # For data manipulation (filter, arrange, mutate)
library(tidyr)    # For reshaping data (pivot_longer)
library(ggplot2)  # For creating plots

# --- 2. Define File Paths ---
# Set the input file path and the output directory for the plots.
input_file <- "../public/RNAseq/UpDEG_category_NormCount_log2FC.csv"
output_dir <- "../public/RNAseq/"

# Create the output directory if it doesn't already exist to prevent errors.
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. Import and Preprocess Data ---
# Read the CSV file into a data frame.
# 'check.names = FALSE' is used to prevent R from changing column names
# that have special characters (like "...").
rna_data <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

# --- 4. Sort Data by Categories ---
# As requested, the data is sorted first by "Category...SRA.specific" and then
# by "Category...other". We define the desired order for each column.
sra_specific_order <- c("Delayed", "Early sustained", "Early transient", "Late sustained", "Late transient", "Stochastic")
other_order <- c("ns", "PCD-RFP", "PCD-RFP & SRA")

# The sorting is done by converting the category columns to factors with the
# defined levels and then using arrange().
rna_data_sorted <- rna_data %>%
  mutate(
    `Category...SRA.specific` = factor(`Category...SRA.specific`, levels = sra_specific_order),
    `Category...other` = factor(`Category...other`, levels = other_order)
  ) %>%
  arrange(`Category...SRA.specific`, `Category...other`)

# --- 5. Consolidate Categories and Reshape Data for Plotting ---
# Combine the two category columns into a single 'Category' column for easier processing.
# If `Category...SRA.specific` is empty, the value from `Category...other` is used.
data_with_single_category <- rna_data_sorted %>%
  mutate(Category = ifelse(!is.na(`Category...SRA.specific`) & `Category...SRA.specific` != "",
                           as.character(`Category...SRA.specific`),
                           as.character(`Category...other`))) %>%
  # Filter out any rows that do not belong to one of the target categories.
  filter(!is.na(Category) & Category != "")

# Reshape the data from a wide format to a long format. This is the standard
# approach for creating plots with multiple groups in ggplot2.
mean_columns <- c(
  "MCF7_023_0_Mean", "MCF7_023_24_Mean", "MCF7_025_0_Mean",
  "MCF7_025_24_Mean", "MCF7_021_0_Mean", "MCF7_021_10_Mean",
  "MCF7_021_24_Mean", "MCF7_021_48_Mean"
)

long_format_data <- data_with_single_category %>%
  pivot_longer(
    cols = all_of(mean_columns),
    names_to = "Sample",
    values_to = "Normalized_Count"
  )

# --- 6. Generate and Export Plots for Each Category ---
# Get the list of unique categories to loop through.
unique_categories <- unique(long_format_data$Category)

# Loop through each category, create a plot, and save it.
for (current_category in unique_categories) {
  
  # Filter the data for the current category.
  category_subset <- long_format_data %>%
    filter(Category == current_category)
  
  # Create the box and whisker plot using ggplot2.
  # geom_jitter() is added to display the individual data points.
  plot <- ggplot(category_subset, aes(x = Sample, y = Normalized_Count, fill = Sample)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "black") +
    scale_y_log10() + # Set the y-axis to a log10 scale
    labs(
      title = paste("Normalized Gene Counts for:", current_category),
      subtitle = "Boxplot with individual data points",
      x = "Sample Condition",
      y = "Normalized Gene Count (Mean, log10 scale)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none" # The x-axis is already descriptive.
    )
  
  # Create a file-safe name for the category (e.g., replacing spaces and '&').
  safe_filename <- gsub(" & ", "_and_", current_category)
  safe_filename <- gsub(" ", "_", safe_filename)
  output_filename <- paste0(safe_filename, "_Normalized_GeneCount_BoxPlot.jpg")
  output_filepath <- file.path(output_dir, output_filename)
  
  # Save the plot to the specified directory.
  ggsave(
    output_filepath,
    plot = plot,
    device = "jpeg",
    width = 12,
    height = 9,
    units = "in",
    dpi = 300
  )
  
  # Print a confirmation message to the console.
  cat(paste("Successfully exported plot for category '", current_category, "' to:\n", output_filepath, "\n\n", sep = ""))
}

cat("--- All plots have been generated successfully. ---\n")

################################################################################

################################################################################

# Title: RNA-seq Data Analysis and Visualization
# Description: This script imports RNA-seq data, sorts it by category,
# and generates violin plots of log2 fold change for each category.

# --- 1. SETUP: Load necessary packages ---
# Ensure you have these packages installed: install.packages(c("dplyr", "tidyr", "ggplot2"))
library(dplyr)
library(tidyr)
library(ggplot2)

cat("Step 1: Packages loaded.\n")

# --- 2. DATA IMPORT AND PREPARATION ---
data_file_path <- "../public/RNAseq/UpDEG_category_NormCount_log2FC.csv"

# Check if the file exists before attempting to read it
if (!file.exists(data_file_path)) {
  stop(paste("Error: The data file was not found at the specified path:", data_file_path))
}

# Read the CSV file
rna_data <- read.csv(data_file_path, header = TRUE, stringsAsFactors = FALSE)
cat(paste("Step 2: Successfully imported data from", data_file_path, "\n"))

# Clean column names to make them syntactically valid for R
# This changes names like '023_24_log2FC' to 'X023_24_log2FC'
names(rna_data) <- make.names(names(rna_data))

# Combine the two category columns into a single, comprehensive 'Category' column.
# This uses the value from 'Category...SRA.specific' if available, otherwise it uses 'Category...other'.
rna_data <- rna_data %>%
  mutate(
    Category = ifelse(
      Category...SRA.specific != "" & !is.na(Category...SRA.specific),
      Category...SRA.specific,
      Category...other
    )
  ) %>%
  # Remove any rows where the final Category is blank or NA
  filter(Category != "" & !is.na(Category))

cat("Step 3: Cleaned column names and combined category columns.\n")


# --- 3. SORTING DATA BY CATEGORY ---
# Define the specific order for the categories as requested
sra_specific_order <- c("Delayed", "Early sustained", "Early transient", "Late sustained", "Late transient", "Stochastic")
other_order <- c("ns", "PCD-RFP", "PCD-RFP & SRA")
full_category_order <- c(sra_specific_order, other_order)

# Convert the 'Category' column to an ordered factor based on the defined order
rna_data$Category <- factor(rna_data$Category, levels = full_category_order)

# Sort the entire data frame based on the ordered 'Category' factor
rna_data_sorted <- rna_data %>%
  arrange(Category)

cat("Step 4: Sorted data frame by the specified category order.\n")


# --- 4. GENERATE AND EXPORT PLOTS ---
# Define the output directory and create it if it doesn't exist
output_directory <- "../public/RNAseq/"
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

# Get the list of unique categories present in the data (will be in the correct factor order)
unique_categories_in_data <- levels(rna_data_sorted$Category)[levels(rna_data_sorted$Category) %in% unique(rna_data_sorted$Category)]

cat("Step 5: Starting plot generation for each category...\n\n")

# !!! STEP 1: DEFINE YOUR CUSTOM ORDER HERE !!!
# Assuming your conditions are named "0", "1", "2", "3", "4" after cleaning.
# Change these strings if your condition names are different.
condition_order <- c("023_24_log2FC", "025_24_log2FC", "021_10_log2FC", "021_24_log2FC", "021_48_log2FC") 

# Loop through each unique category to create and save a plot
for (current_category in unique_categories_in_data) {
  
  # Filter the data for the current category
  category_subset <- rna_data_sorted %>%
    filter(Category == current_category)
  
  # Reshape the data from a wide format to a long format suitable for ggplot2
  plot_data_long <- category_subset %>%
    pivot_longer(
      cols = starts_with("X") & ends_with("log2FC"),
      names_to = "Condition",
      values_to = "log2FC"
    ) %>%
    # Clean up the condition names for better axis labels
    mutate(Condition = gsub("^X|\\.log2FC$", "", Condition) %>% gsub("\\.", "_", .)) %>%
    
    # !!! STEP 2: APPLY THE ORDER BY CONVERTING TO A FACTOR !!!
    mutate(Condition = factor(Condition, levels = condition_order))
  
  # Proceed only if there's data to plot for the category
  if (nrow(plot_data_long) > 0) {
    # Define y-axis breaks with a step of 1, and include the custom 0.58 tick
    y_range <- range(plot_data_long$log2FC, na.rm = TRUE)
    y_axis_breaks <- sort(unique(c(seq(floor(y_range[1]), ceiling(y_range[2]), by = 1), 0.58)))
    
    # Create the violin plot, with an overlaid boxplot, and individual data points as jitter
    p <- ggplot(plot_data_long, aes(x = Condition, y = log2FC)) +
      geom_violin(
        trim = FALSE, # Show the full density plot without trimming tails
        fill = "#4e79a7",  # A nice blue color
        color = "#2c4660",
        alpha = 0.7
      ) +
      geom_jitter(
        width = 0.2,       # Control the horizontal spread of points
        alpha = 0.6,       # Make points semi-transparent
        color = "black",
        size = 1.5
      ) +
      geom_boxplot(
        width = 0.1,
        fill = "white",
        outlier.shape = NA # We use jitter to show all points, so hide default outliers
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = 0.58, linetype = "dashed", color = "red") +
      scale_y_continuous(breaks = y_axis_breaks) +
      labs(
        title = paste("log2 Fold Change for Category:", current_category),
        subtitle = "Each point represents a transcript",
        x = "Experimental Condition",
        y = "log2 Fold Change"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()
      )
    
    # Sanitize the category name to create a valid and safe file name
    safe_filename <- gsub("[^a-zA-Z0-9_.-]", "_", current_category)
    output_path <- file.path(output_directory, paste0(safe_filename, "_log2FC.jpg"))
    
    # Save the plot as a JPEG file
    ggsave(
      filename = output_path,
      plot = p,
      device = "jpeg",
      width = 10,
      height = 8,
      dpi = 300
    )
    
    cat(paste0("-> Successfully created plot: '", output_path, "'\n"))
    
  } else {
    cat(paste0("-> Skipping category '", current_category, "' as it contained no data.\n"))
  }
}

cat("\nStep 6: All plots have been generated and saved successfully.\n")



################################################################################

################################################################################

# R Script to Generate Volcano Plots from DEG Analysis Files

# --- 1. Setup: Install and Load Required Packages ---

# This section checks if the required packages are installed and, if not, installs them.
required_packages <- c("ggplot2", "dplyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")

# Load the packages
library(ggplot2)
library(dplyr)
cat("Successfully loaded required packages: ggplot2, dplyr\n")


# --- 2. Define File Paths and Parameters ---

# List of input CSV files containing differential expression data
input_files <- c(
  "../public/RNAseq/MCF7_023_24_DEG.csv",
  "../public/RNAseq/MCF7_025_24_DEG.csv",
  "../public/RNAseq/MCF7_021_10_DEG.csv",
  "../public/RNAseq/MCF7_021_24_DEG.csv",
  "../public/RNAseq/MCF7_021_48_DEG.csv"
)

# Define the significance thresholds
padj_threshold <- 0.05
# The log2FoldChange cutoff is calculated from a 1.5-fold change
log2fc_threshold <- log2(1.5)


# --- 3. Main Loop: Process Each File and Generate Plot ---

# This loop iterates over each file path specified in the 'input_files' vector.
for (file_path in input_files) {
  
  # Check if the file exists before trying to read it
  if (!file.exists(file_path)) {
    cat("Warning: File not found at", file_path, "- Skipping.\n")
    next # Skip to the next file in the loop
  }
  
  cat("Processing file:", file_path, "\n")
  
  # --- Data Preparation ---
  
  # Read the CSV data into a data frame
  deg_data <- read.csv(file_path)
  
  # Prepare the data for plotting by adding a column to classify genes
  deg_data <- deg_data %>%
    filter(!is.na(padj)) %>%
    mutate(
      status = case_when(
        log2FoldChange > log2fc_threshold & padj < padj_threshold  ~ "Upregulated",
        log2FoldChange < -log2fc_threshold & padj < padj_threshold ~ "Downregulated",
        TRUE                                                      ~ "Not Significant"
      ),
      status = factor(status, levels = c("Upregulated", "Downregulated", "Not Significant"))
    )
  
  # --- Gene Counting ---
  
  upregulated_count <- sum(deg_data$status == "Upregulated")
  downregulated_count <- sum(deg_data$status == "Downregulated")
  
  cat("  - Upregulated genes:", upregulated_count, "\n")
  cat("  - Downregulated genes:", downregulated_count, "\n")
  
  # --- Plot Generation (ggplot2) ---
  
  plot_title <- basename(file_path) %>% sub("_DEG.csv", "", .)
  x_axis_limit <- max(abs(deg_data$log2FoldChange), na.rm = TRUE) * 1.05
  y_max <- max(-log10(deg_data$padj), na.rm = TRUE)
  
  volcano_plot <- ggplot(
    data = deg_data, 
    aes(x = log2FoldChange, y = -log10(padj))
  ) +
    geom_point(aes(color = status), alpha = 0.6, size = 1.8) +
    scale_color_manual(
      name = "Gene Regulation",
      values = c("Upregulated" = "#E64B35", "Downregulated" = "#3C5488", "Not Significant" = "grey")
    ) +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
    
    # --- FIX 1: Add clip = 'off' to prevent labels from being cut off ---
    coord_cartesian(xlim = c(-x_axis_limit, x_axis_limit), clip = 'off') +
    
    annotate("text", x = x_axis_limit, y = -log10(padj_threshold), 
             label = paste("padj =", padj_threshold), vjust = -0.5, hjust = 1.0, size = 4) +
    annotate("text", x = log2fc_threshold, y = y_max * 0.95, 
             label = paste0(round(2^log2fc_threshold, 1), "-fold"), hjust = -0.2, size = 4) +
    
    # --- FIX 2: Add vjust > 1 to push the text down from the top edge ---
    annotate(
      "text", 
      x = x_axis_limit * 0.6, 
      y = y_max, # Anchor to the top...
      label = paste("Up:", upregulated_count), 
      color = "#E64B35", 
      size = 5,
      fontface = "bold",
      hjust = 0.5,
      vjust = 1.5  # ...but push it down slightly
    ) +
    annotate(
      "text", 
      x = -x_axis_limit * 0.6, 
      y = y_max, # Anchor to the top...
      label = paste("Down:", downregulated_count), 
      color = "#3C5488", 
      size = 5,
      fontface = "bold",
      hjust = 0.5,
      vjust = 1.5 # ...but push it down slightly
    ) +
    
    labs(
      title = paste("Volcano Plot for", plot_title),
      x = bquote(~Log[2]~ "Fold Change"),
      y = bquote(~-Log[10]~ "(Adjusted P-value)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # --- Export Plot ---
  
  output_path <- sub("_DEG.csv", "_DEG_VolcanoPlot.jpg", file_path)
  ggsave(
    filename = output_path,
    plot = volcano_plot,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300
  )
  
  cat("  - Plot saved to:", output_path, "\n\n")
}

cat("--- Script finished ---\n")