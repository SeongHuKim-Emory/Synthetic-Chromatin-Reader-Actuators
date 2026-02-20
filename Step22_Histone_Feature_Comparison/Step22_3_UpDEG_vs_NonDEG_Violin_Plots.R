# R code to compare features between SRA-Specific UpDEG and Non-DEG genes with box plots

# Step 0: Install and load necessary packages
# If you don't have these packages installed, uncomment the following lines to install them.
# install.packages("ggplot2")
# install.packages("ggpubr")
# install.packages("here")

# The original script used setwd, but also loaded 'here'.
# We will follow the file path logic from the original script.
# setwd("X:/haynes-lab/Kim_SH/Projects/NIH R21/Codes/src/Step21_Feature_Comparison")

library(ggplot2)
library(ggpubr)
library(here) # Using the 'here' package for robust path management

# Step 1: Import files using the 'here' package
# This assumes your R session is running within a project where 'public' is a subfolder.
# Note: Using here() makes the setwd() call redundant if the R session/project root is correct.
# We'll use the relative paths as provided.
updeg_file_path <- "../../public/Multiomics/Step22_Histone_Feature_Comparison/SRA_Specific_UpDEG_with_GO_CCLE_Select_Features.csv"
nondeg_file_path <- "../../public/Multiomics/Step22_Histone_Feature_Comparison/Non_DEG_with_GO_CCLE_Select_Features.csv"

# Handle potential 'here' package pathing if the above relative paths fail
# This checks if 'here' is available and the files exist relative to the project root
# If you are not running this in an RStudio Project, the manual paths above will be used.
if (require(here) & all(file.exists(here(updeg_file_path)), file.exists(here(nondeg_file_path)))) {
  updeg_data <- read.csv(here(updeg_file_path), header = TRUE, check.names = TRUE)
  nondeg_data <- read.csv(here(nondeg_file_path), header = TRUE, check.names = TRUE)
  print("Loaded data using 'here' package paths.")
} else {
  print("Loading data using relative paths. Make sure your working directory is correct.")
  # This setwd is from the original script. Uncomment if you need to set it manually.
  # setwd("X:/haynes-lab/Kim_SH/Projects/NIH R21/Codes/src/Step21_Feature_Comparison")
  updeg_data <- read.csv(updeg_file_path, header = TRUE, check.names = TRUE)
  nondeg_data <- read.csv(nondeg_file_path, header = TRUE, check.names = TRUE)
}


# Step 2: Define feature list and loop to create plots
feature_columns <- c(
  "Max_DSPP_025", "Mean_DSPP_025", "Max_DSPP_H3K4me1", "Mean_DSPP_H3K4me1", "Max_DSPP_H3K4me2", "Mean_DSPP_H3K4me2",
  "Max_DSPP_H3K4me3", "Mean_DSPP_H3K4me3", "Max_DSPP_H3K9ac", "Mean_DSPP_H3K9ac", "Max_DSPP_H3K9me3", "Mean_DSPP_H3K9me3",
  "Max_DSPP_H3K27ac", "Mean_DSPP_H3K27ac", "Max_DSPP_H3K27me3", "Mean_DSPP_H3K27me3", "Max_DSPP_H3K36me3", "Mean_DSPP_H3K36me3",
  "Max_DSPP_H4K20me1", "Mean_DSPP_H4K20me1", "Max_CP_025", "Mean_CP_025", "Max_CP_H3K4me1", "Mean_CP_H3K4me1",
  "Max_CP_H3K4me2", "Mean_CP_H3K4me2", "Max_CP_H3K4me3", "Mean_CP_H3K4me3", "Max_CP_H3K9ac", "Mean_CP_H3K9ac",
  "Max_CP_H3K9me3", "Mean_CP_H3K9me3", "Max_CP_H3K27ac", "Mean_CP_H3K27ac", "Max_CP_H3K27me3", "Mean_CP_H3K27me3",
  "Max_CP_H3K36me3", "Mean_CP_H3K36me3", "Max_CP_H4K20me1", "Mean_CP_H4K20me1", "Max_USPP_025", "Mean_USPP_025",
  "Max_USPP_H3K4me1", "Mean_USPP_H3K4me1", "Max_USPP_H3K4me2", "Mean_USPP_H3K4me2", "Max_USPP_H3K4me3", "Mean_USPP_H3K4me3",
  "Max_USPP_H3K9ac", "Mean_USPP_H3K9ac", "Max_USPP_H3K9me3", "Mean_USPP_H3K9me3", "Max_USPP_H3K27ac", "Mean_USPP_H3K27ac",
  "Max_USPP_H3K27me3", "Mean_USPP_H3K27me3", "Max_USPP_H3K36me3", "Mean_USPP_H3K36me3", "Max_USPP_H4K20me1", "Mean_USPP_H4K20me1",
  "Max_Enhancer_025", "Mean_Enhancer_025", "Max_Enhancer_H3K4me1", "Mean_Enhancer_H3K4me1", "Max_Enhancer_H3K4me2", "Mean_Enhancer_H3K4me2",
  "Max_Enhancer_H3K4me3", "Mean_Enhancer_H3K4me3", "Max_Enhancer_H3K9ac", "Mean_Enhancer_H3K9ac", "Max_Enhancer_H3K9me3", "Mean_Enhancer_H3K9me3",
  "Max_Enhancer_H3K27ac", "Mean_Enhancer_H3K27ac", "Max_Enhancer_H3K27me3", "Mean_Enhancer_H3K27me3", "Max_Enhancer_H3K36me3", "Mean_Enhancer_H3K36me3",
  "Max_Enhancer_H4K20me1", "Mean_Enhancer_H4K20me1"
)
# Note: R's read.csv automatically converts special characters like ':' in column names
# to '.', so "GO:0099050" becomes "GO.0099050". The 'check.names=TRUE' default handles this.

# Create a directory to save plots if it doesn't exist
output_dir <- "../../public/Multiomics/Step22_Histone_Feature_Comparison"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Initialize a data frame to store p-values
p_values_summary <- data.frame(Feature = character(), P_value = numeric(), stringsAsFactors = FALSE)

# --- Define custom colors ---
box_colors <- c("SRA-Specific UpDEG" = "red", "Non-DEG" = "blue")

# Loop through each feature to generate a plot and calculate statistics
print("Generating individual feature plots...")
for (feature in feature_columns) {
  
  # Prepare data for plotting
  combined_data <- rbind(
    data.frame(group = "SRA-Specific UpDEG", value = updeg_data[[feature]]),
    data.frame(group = "Non-DEG", value = nondeg_data[[feature]])
  )
  
  # Add a small constant to values to prevent issues with log10(0)
  combined_data$value <- combined_data$value + 0.001
  
  # Perform t-test to get p-value
  ttest_result <- t.test(value ~ group, data = combined_data)
  
  # Store the feature and its p-value
  p_values_summary <- rbind(p_values_summary, data.frame(Feature = feature, P_value = ttest_result$p.value))
  
  # Create the box plot
  box_plot <- ggplot(combined_data, aes(x = group, y = value, fill = group)) +
    
    # --- MODIFICATION: Made box width thinner (0.3) ---
    geom_boxplot(width = 0.3, alpha = 0.7, outlier.shape = NA) + # outlier.shape = NA hides outliers
    
    # --- MODIFICATION: Removed geom_jitter ---
    # geom_jitter(shape = 16, position = position_jitter(width = 0.2), size = 0.8, alpha = 0.6) +
    
    # Apply custom colors
    scale_fill_manual(values = box_colors) +
    
    # Apply log10 scale to the y-axis
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(
      title = paste("Box Plot Comparison of", feature), # Updated title
      subtitle = "SRA-Specific UpDEG vs. Non-DEG",
      x = "Group",
      y = "Value (log10 scale)"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    
    # --- MODIFICATION: Updated stat_compare_means to add p-value bar ---
    # By adding the 'comparisons' argument, ggpubr will draw the bracket.
    stat_compare_means(
      method = "t.test",
      comparisons = list(c("SRA-Specific UpDEG", "Non-DEG")),
      label = "p.format", # Use "p.format" for the value, or "p.signif" for stars
      size = 4
    )
  
  # Step 3: Export each box plot
  # Updated file name to reflect box plot
  plot_file_path <- file.path(output_dir, paste0(feature, "_UpDEG_vs_NonDEG_Boxplot.jpg"))
  ggsave(plot_file_path, plot = box_plot, width = 8, height = 6, dpi = 300)
}
print("Individual plots complete.")

# Step 4: Export the p-values table
p_values_file_path <- file.path(output_dir, "UpDEG_vs_NonDEG_p_values.csv")
write.csv(p_values_summary, p_values_file_path, row.names = FALSE)
print(paste("P-value summary table saved to:", p_values_file_path))


# ==============================================================================
# Step 5: Generate Grouped Box Plots (NEW SECTION)
# ==============================================================================

print("Starting generation of grouped box plots...")

# Define the 8 groups based on their prefixes
groups_list <- c(
  "Max_CP", "Max_DSPP", "Max_Enhancer", "Max_USPP",
  "Mean_CP", "Mean_DSPP", "Mean_Enhancer", "Mean_USPP"
)

for (group_prefix in groups_list) {
  group_name <- paste0(group_prefix, "_Group")
  
  # Find all feature columns that belong to this group
  group_features <- grep(paste0("^", group_prefix, "_"), feature_columns, value = TRUE)
  
  if (length(group_features) == 0) {
    print(paste("No features found for group:", group_prefix, ". Skipping."))
    next
  }
  
  # --- a. Prepare long-form data for this group ---
  long_data_list <- list()
  for (feat in group_features) {
    # Add data for UpDEG
    long_data_list[[length(long_data_list) + 1]] <- data.frame(
      feature_name = feat,
      group = "SRA-Specific UpDEG",
      value = updeg_data[[feat]]
    )
    # Add data for Non-DEG
    long_data_list[[length(long_data_list) + 1]] <- data.frame(
      feature_name = feat,
      group = "Non-DEG",
      value = nondeg_data[[feat]]
    )
  }
  combined_long_data <- do.call(rbind, long_data_list)
  
  # Add small constant for log10 scale
  combined_long_data$value <- combined_long_data$value + 0.001
  
  # --- MODIFICATION: Create short labels for x-axis ---
  # 1. Define the new labels by removing prefix
  new_labels <- gsub(paste0(group_prefix, "_"), "", group_features)
  # 2. Substitute "025" with "PCD-RFP"
  new_labels[new_labels == "025"] <- "PCD-RFP"
  
  # 3. Add short_label column to main data
  combined_long_data$short_label <- gsub(paste0(group_prefix, "_"), "", combined_long_data$feature_name)
  combined_long_data$short_label[combined_long_data$short_label == "025"] <- "PCD-RFP"
  
  # 4. Factor the new short labels to maintain order
  combined_long_data$short_label <- factor(combined_long_data$short_label, levels = new_labels)
  
  # Factor original feature_name (still useful for joining)
  combined_long_data$feature_name <- factor(combined_long_data$feature_name, levels = group_features)
  
  
  # --- b. Prepare p-value data for brackets ---
  
  # Get p-values for just this group's features
  pvals_subset <- subset(p_values_summary, Feature %in% group_features)
  
  # Find the maximum y-value for the entire group to position brackets
  overall_max_y <- max(combined_long_data$value, na.rm = TRUE)
  bracket_y_pos <- 10^(log10(overall_max_y) + 0.3) # Position brackets 0.3 OOMs above max
  
  # Format the data for stat_pvalue_manual
  stat_data <- pvals_subset
  stat_data$group1 <- "SRA-Specific UpDEG"
  stat_data$group2 <- "Non-DEG"
  
  # --- MODIFICATION: Format p-value string (p=0.118 or p<0.001) ---
  raw_p_values <- scales::pvalue(stat_data$P_value, accuracy = 0.001)
  stat_data$label <- ifelse(startsWith(raw_p_values, "<"), 
                            paste0("p", raw_p_values), 
                            paste0("p=", raw_p_values))
  
  # Set y-position for brackets
  stat_data$y.position <- bracket_y_pos
  # Rename 'Feature' to 'feature_name' to match the plot's x-axis aesthetic
  colnames(stat_data)[colnames(stat_data) == 'Feature'] <- 'feature_name'
  
  # --- MODIFICATION: Add short_label to p-value data ---
  stat_data$short_label <- gsub(paste0(group_prefix, "_"), "", stat_data$feature_name)
  stat_data$short_label[stat_data$short_label == "025"] <- "PCD-RFP"
  # We need to get the factor levels from the data prep step
  stat_data$short_label <- factor(stat_data$short_label, levels = new_labels)
  
  
  # --- c. Create the grouped plot ---
  plot_title <- paste("Box Plot Comparison for", group_name)
  
  # --- MODIFICATION: Use 'short_label' for x-axis ---
  grouped_plot <- ggplot(combined_long_data, aes(x = short_label, y = value, fill = group)) +
    
    # Add dodged boxplots
    geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA) +
    
    # Apply custom colors
    scale_fill_manual(values = box_colors) +
    
    # Apply log10 scale
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    
    # Add the p-value brackets
    # --- MODIFICATION: Use 'short_label' for x mapping ---
    ggpubr::stat_pvalue_manual(
      stat_data,
      x = "short_label", # Map to the x-axis variable
      y.position = "y.position",
      label = "label",
      position = position_dodge(0.8) # Dodge brackets to match boxes
    ) +
    
    # Labels and titles
    labs(
      title = plot_title,
      subtitle = "SRA-Specific UpDEG vs. Non-DEG",
      x = "Feature",
      y = "Value (log10 scale)",
      fill = "Group" # Add a legend title
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      # Rotate x-axis labels for readability
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom" # Move legend to the bottom
    )
  
  # --- d. Save the grouped plot ---
  grouped_plot_file_path <- file.path(output_dir, paste0(group_name, "_Grouped_Boxplot.jpg"))
  ggsave(grouped_plot_file_path, plot = grouped_plot, width = 12, height = 8, dpi = 300)
  
  print(paste("Saved grouped plot:", grouped_plot_file_path))
}


print("Processing complete. All individual plots, grouped plots, and the p-value summary table have been generated and saved.")

