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
updeg_file_path <- "../../public/Multiomics/Step21_Feature_Comparison/SRA_Specific_UpDEG_with_GO_CCLE_Final_Features.csv"
nondeg_file_path <- "../../public/Multiomics/Step21_Feature_Comparison/Non_DEG_with_GO_CCLE_Final_Features.csv"

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
  "Mean_DSPP_025", "Mean_DSPP_CEBPB", "Max_DSPP_CHD1", "Mean_DSPP_CHD1", "Max_DSPP_CLOCK", "Max_DSPP_CTCF",
  "Max_DSPP_E4F1", "Mean_DSPP_EGR1", "Mean_DSPP_EP300", "Mean_DSPP_FOS", "Max_DSPP_GATA3", "Max_DSPP_H3K4me3",
  "Max_DSPP_H3K9ac", "Max_DSPP_H3K9me3", "Max_DSPP_H3K27ac", "Mean_DSPP_H3K27ac", "Max_DSPP_H3K27me3",
  "Mean_DSPP_H3K27me3", "Max_DSPP_H3K36me3", "Mean_DSPP_H3K36me3", "Max_DSPP_H4K20me1", "Mean_DSPP_HCFC1",
  "Mean_DSPP_HDAC2", "Max_DSPP_HDGF", "Max_DSPP_HSF1", "Max_DSPP_JUND", "Mean_DSPP_LARP7", "Max_DSPP_MAX",
  "Max_DSPP_MAZ", "Mean_DSPP_PAX8", "Mean_DSPP_PPP1R10", "Max_DSPP_RAD21", "Max_DSPP_RAD51", "Max_DSPP_SIX4",
  "Mean_DSPP_SREBF1", "Max_DSPP_SRF", "Mean_DSPP_TARDBP", "Max_DSPP_TRIM22", "Mean_DSPP_ZBTB11",
  "Mean_DSPP_ZBTB40", "Mean_DSPP_ZFX", "Max_DSPP_ZKSCAN1", "Mean_DSPP_ZKSCAN1", "Max_DSPP_ZNF217",
  "Max_DSPP_ZNF512B", "Max_DSPP_ZNF574", "Mean_DSPP_ZNF574", "Max_DSPP_ZNF687", "Mean_DSPP_ZNF687",
  "Mean_CP_025", "Mean_CP_CEBPB", "Max_CP_COPS2", "Max_CP_CUX1", "Max_CP_EGR1", "Mean_CP_ELK1",
  "Max_CP_FOSL2", "Mean_CP_FOXA1", "Mean_CP_FOXK2", "Max_CP_GABPA", "Max_CP_H3K9ac", "Mean_CP_H3K27ac",
  "Max_CP_H3K36me3", "Mean_CP_HDGF", "Max_CP_HES1", "Mean_CP_MBD2", "Max_CP_NFXL1", "Mean_CP_POLR2A",
  "Max_CP_RAD51", "Mean_CP_RAD51", "Max_CP_RFX1", "Mean_CP_RFX5", "Mean_CP_SNIP1", "Mean_CP_TCF7L2",
  "Max_CP_TEAD4", "Mean_CP_ZFX", "Max_CP_ZHX2", "Max_CP_ZNF8", "Mean_CP_ZNF8", "Max_CP_ZNF24",
  "Mean_CP_ZNF512B", "Mean_USPP_ATF7", "Max_USPP_CREB1", "Max_USPP_DPF2", "Max_USPP_EGR1", "Mean_USPP_ELF1",
  "Max_USPP_ELK1", "Mean_USPP_EP300", "Max_USPP_FOXK2", "Mean_USPP_FOXM1", "Max_USPP_H3K4me1",
  "Mean_USPP_H3K4me2", "Mean_USPP_H3K27me3", "Mean_USPP_H3K36me3", "Max_USPP_HDAC2", "Mean_USPP_JUN",
  "Mean_USPP_JUND", "Max_USPP_MAX", "Mean_USPP_MAZ", "Max_USPP_MBD2", "Mean_USPP_MBD2", "Max_USPP_MTA3",
  "Max_USPP_NCOA3", "Mean_USPP_NFRKB", "Max_USPP_PPP1R10", "Max_USPP_RAD51", "Max_USPP_REST",
  "Max_USPP_SMARCA5", "Max_USPP_SNIP1", "Max_USPP_SREBF1", "Max_USPP_TEAD4", "Max_USPP_ZBTB11",
  "Mean_USPP_ZBTB33", "Mean_USPP_ZFX", "Mean_USPP_ZNF8", "Mean_USPP_ZNF574", "Mean_USPP_ZNF579",
  "Max_USPP_ZNF687", "Mean_USPP_ZNF687", "Max_Enhancer_ARID3A", "Max_Enhancer_BMI1", "Mean_Enhancer_BMI1",
  "Mean_Enhancer_COPS2", "Mean_Enhancer_DPF2", "Mean_Enhancer_E2F8", "Mean_Enhancer_FOS",
  "Max_Enhancer_H3K27me3", "Mean_Enhancer_H3K27me3", "Max_Enhancer_H3K36me3", "Max_Enhancer_H4K20me1",
  "Mean_Enhancer_H4K20me1", "Mean_Enhancer_HDGF", "Mean_Enhancer_HSF1", "Max_Enhancer_LARP7",
  "Mean_Enhancer_LARP7", "Mean_Enhancer_MBD2", "Max_Enhancer_MTA2", "Max_Enhancer_NEUROD1",
  "Mean_Enhancer_NFRKB", "Max_Enhancer_PAX8", "Mean_Enhancer_TEAD4", "Max_Enhancer_ZFX", "Mean_Enhancer_ZNF8",
  "Max_Enhancer_ZNF207", "Max_Enhancer_ZNF592", "Mean_Enhancer_ZNF592", "TPM", "GO.0099050", "GO.0005788",
  "GO.2000741", "MS_Value"
)
# Note: R's read.csv automatically converts special characters like ':' in column names
# to '.', so "GO:0099050" becomes "GO.0099050". The 'check.names=TRUE' default handles this.

# Create a directory to save plots if it doesn't exist
output_dir <- "../../public/Multiomics/Step21_Feature_Comparison"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Initialize a data frame to store p-values
p_values_summary <- data.frame(Feature = character(), P_value = numeric(), stringsAsFactors = FALSE)

# Loop through each feature to generate a plot and calculate statistics
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
  
  # Define custom colors
  # --- MODIFICATION: Updated box colors ---
  box_colors <- c("SRA-Specific UpDEG" = "red", "Non-DEG" = "blue") 
  
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

# Step 4: Export the p-values table
p_values_file_path <- file.path(output_dir, "UpDEG_vs_NonDEG_p_values.csv")
write.csv(p_values_summary, p_values_file_path, row.names = FALSE)

print("Processing complete. All box plots and the p-value summary table have been generated and saved.")

