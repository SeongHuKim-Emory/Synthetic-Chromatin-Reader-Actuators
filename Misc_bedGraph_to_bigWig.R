setwd("C:/Users/ski2959/Desktop/NIH R21/Codes/src")

################################################################################

################################################################################

# R Script to Convert .bedGraph to .bigWig

# --- 1. Install and Load rtracklayer ---
# Check if rtracklayer is installed, if not, provide installation instructions.
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  message("The 'rtracklayer' package is not installed.")
  message("To install it, run the following commands in your R console:")
  message('if (!requireNamespace("BiocManager", quietly = TRUE))')
  message('    install.packages("BiocManager")')
  message('BiocManager::install("rtracklayer")')
  stop("Please install 'rtracklayer' and try again.", call. = FALSE)
}
library(rtracklayer)

# --- 2. Define File Paths ---
# Input .bedGraph file
bedGraph_file <- "../public/Public Dataset/ChIPseq/H3K27me3/GSE139199 GSM4133800 hg38 MCF7 H3K27me3 ChIPseq.bedGraph"

# Output .bigWig file
bigWig_file <- "../public/Public Dataset/ChIPseq/H3K27me3/GSE139199 GSM4133800 hg38 MCF7 H3K27me3 ChIPseq.bigWig"

# Genome assembly (from filename, e.g., hg38)
genome_assembly <- "hg38" # This is important for fetching chromosome lengths

# --- 3. Check if input file exists ---
if (!file.exists(bedGraph_file)) {
  stop(paste("Input file not found:", bedGraph_file), call. = FALSE)
}

# --- 4. Import .bedGraph file ---
message(paste("Importing .bedGraph file:", bedGraph_file))
tryCatch({
  # Import the bedGraph file. The 'score' column will be automatically recognized.
  gr_data <- import(bedGraph_file, format = "bedGraph")
  message("Successfully imported .bedGraph file.")
}, error = function(e) {
  stop(paste("Error importing .bedGraph file:", e$message), call. = FALSE)
})

# --- 5. Fetch Chromosome Information ---
# .bigWig files require chromosome lengths (seqinfo).
# We will try to fetch this information from UCSC for the specified genome.
message(paste("Fetching chromosome information for genome:", genome_assembly))
tryCatch({
  # Attempt to get seqinfo for the specified genome
  # This requires an internet connection.
  # Sometimes UCSC might be slow or temporarily unavailable.
  seq_info <- Seqinfo(genome = genome_assembly)
  
  # If Seqinfo(genome=genome_assembly) doesn't directly fetch lengths for all chromosomes
  # present in your data, you might need to fetch them more explicitly or provide a .chrom.sizes file.
  # For common assemblies, this usually works.
  # Let's ensure the seqlevels in our data match those in seq_info
  # and assign the seqinfo to our GRanges object.
  
  # Check if all seqlevels in gr_data are present in seq_info
  missing_seqlevels <- setdiff(seqlevels(gr_data), seqlevels(seq_info))
  if (length(missing_seqlevels) > 0) {
    warning(paste("The following chromosomes from your bedGraph are not found in the fetched UCSC seqinfo for", genome_assembly, ":", paste(missing_seqlevels, collapse=", ")))
    warning("This might lead to issues during bigWig conversion or an incomplete bigWig file.")
    warning("Consider providing a custom .chrom.sizes file or ensuring chromosome names match UCSC (e.g., 'chr1' vs '1').")
  }
  
  # Filter seq_info to only include chromosomes present in the data,
  # and ensure the order matches. This helps prevent errors if the
  # bedGraph contains fewer chromosomes than the full genome assembly.
  common_seqlevels <- intersect(seqlevels(gr_data), seqlevels(seq_info))
  if (length(common_seqlevels) == 0) {
    stop("No common chromosomes found between your data and the fetched genome information. Please check chromosome naming (e.g. 'chr1' vs '1') and the genome assembly string.")
  }
  
  # Reorder seq_info to match gr_data's seqlevels if necessary, and subset
  # This step is crucial: seqinfo must be compatible with the GRanges object.
  seqlevels(gr_data, pruning.mode="coarse") <- common_seqlevels # Keep only common seqlevels in gr_data
  current_seq_info <- seq_info[common_seqlevels] # Subset and reorder seq_info
  
  # Assign the fetched and filtered sequence information to the GRanges object
  seqinfo(gr_data) <- current_seq_info
  
  message("Successfully fetched and assigned chromosome information.")
  
  # Verify that seqlengths are now present
  if(any(is.na(seqlengths(gr_data)))){
    warning("Some chromosomes in your data still have NA seqlengths after fetching info. This may cause issues.")
    message("Chromosomes with NA lengths: ", paste(seqlevels(gr_data)[is.na(seqlengths(gr_data))], collapse=", "))
  }
  
}, error = function(e) {
  warning(paste("Could not automatically fetch chromosome information for", genome_assembly, "from UCSC:", e$message))
  warning("Exporting to .bigWig might fail or produce an invalid file without correct chromosome lengths.")
  warning("You might need to create a Seqinfo object manually with chromosome lengths, e.g., from a .chrom.sizes file.")
  # Example for manual Seqinfo (if UCSC fetch fails):
  # chrom_sizes_file <- "path/to/your/hg38.chrom.sizes"
  # if(file.exists(chrom_sizes_file)){
  #   cs_df <- read.table(chrom_sizes_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  #   si <- Seqinfo(seqnames=cs_df$V1, seqlengths=cs_df$V2, genome=genome_assembly)
  #   seqinfo(gr_data) <- si[seqlevels(gr_data)] # ensure it matches data
  # } else {
  #   message("Chromosome sizes file not found. Proceeding without full seqinfo.")
  # }
})


# --- 6. Export to .bigWig file ---
# The GRanges object must have a 'score' metadata column for export to bigWig.
# import() from bedGraph should automatically create this from the 4th column.
if (!"score" %in% names(mcols(gr_data))) {
  stop("The imported GRanges object does not have a 'score' metadata column. This is required for .bigWig conversion.", call. = FALSE)
}

message(paste("Exporting to .bigWig file:", bigWig_file))
tryCatch({
  export(gr_data, bigWig_file, format = "bigWig")
  message("Successfully exported to .bigWig file.")
  message(paste("Output file located at:", bigWig_file))
}, error = function(e) {
  stop(paste("Error exporting to .bigWig file:", e$message), call. = FALSE)
})

# --- Script End ---

################################################################################

################################################################################