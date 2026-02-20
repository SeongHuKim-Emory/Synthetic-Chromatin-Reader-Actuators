# R Script for ChIP-seq Peak Annotation and Visualization using ChIPseeker

# --- 1. Installation of Required Packages ---
# This section automatically checks for and installs any missing packages.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of required Bioconductor and CRAN packages
required_packages <- c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "ggupset", "ggimage")

# Loop through the packages, check if they are installed, and install if not
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    print(paste("Installing package:", pkg))
    # Use BiocManager for Bioconductor packages and install.packages for CRAN packages
    if (pkg %in% c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}


# --- 2. Load Libraries ---
# Load the necessary R packages for the analysis.
print("Loading libraries...")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggupset)
library(ggimage)


# --- 3. Define File Paths ---
# Define the paths for the input peak files and the desired output figures.
# Ensure these relative paths are correct based on your working directory.

# Input files
file1_path <- "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2.broadPeak"
file2_path <- "../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2.broadPeak"

# Output files
output1_path <- "../public/ChIPseq/SRA_Visualize_Genomic_Annotation_hg38.jpg"
output2_path <- "../public/ChIPseq/PCD_RFP_Visualize_Genomic_Annotation_hg38.jpg"


# --- 4. Process and Visualize the First File (DBN021) ---

# Print a message to the console to track progress.
print(paste("Processing file:", file1_path))

# Read the peak data from the broadPeak file.
peak1 <- readPeakFile(file1_path)

# Annotate the peaks using the hg38 human genome database.
# tssRegion defines the promoter region around the Transcription Start Site (TSS).
peakAnno1 <- annotatePeak(peak1, tssRegion=c(-3000, 3000),
                          TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                          annoDb="org.Hs.eg.db")

# Export the visualization for the first file.
# We open a JPEG device to save the plot directly to a file.
# Set physical dimensions in inches and resolution (res) for high quality.
jpeg(output1_path, width = 10, height = 8, units = "in", res = 300)

# Generate the upset plot with a Venn pie chart overlay, as requested.
# This visualization shows the distribution and overlap of genomic features.
upsetplot(peakAnno1, vennpie=TRUE)

# Close the JPEG device to finalize the file saving.
dev.off()

# Print a confirmation message.
print(paste("Figure saved to:", output1_path))


# --- 5. Process and Visualize the Second File (DBN025) ---

# Print a message to the console to track progress.
print(paste("Processing file:", file2_path))

# Read the peak data from the second broadPeak file.
peak2 <- readPeakFile(file2_path)

# Annotate the peaks, similar to the first file.
peakAnno2 <- annotatePeak(peak2, tssRegion=c(-3000, 3000),
                          TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                          annoDb="org.Hs.eg.db")

# Export the visualization for the second file.
# Set physical dimensions in inches and resolution (res) for high quality.
jpeg(output2_path, width = 10, height = 8, units = "in", res = 300)

# Generate the upset plot with the Venn pie chart overlay.
upsetplot(peakAnno2, vennpie=TRUE)

# Close the JPEG device.
dev.off()

# Print a final confirmation message.
print(paste("Figure saved to:", output2_path))

print("Analysis complete.")
