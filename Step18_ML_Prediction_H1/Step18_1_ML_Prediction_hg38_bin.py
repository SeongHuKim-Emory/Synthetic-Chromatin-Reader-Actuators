# chipseq_binning.py
#
# This script performs the following steps:
# 1. Fetches hg38 chromosome sizes from the UCSC Genome Browser.
# 2. Bins the standard chromosomes (chr1-22, X, Y) into 200 bp intervals.
# 3. Reads a Protein-of-Interest (PoI) ChIP-seq BED file and nine histone mark ChIP-seq BED files.
# 4. For each 200 bp bin, it checks for overlaps with the peaks in each ChIP-seq file.
# 5. Creates a final table where each row is a bin and columns indicate the presence (1) or absence (0) of an overlap for the PoI and each histone mark.
# 6. Exports the resulting table to a CSV file.
#
# Required libraries:
# pip install pandas pyranges requests tqdm

import pandas as pd
import pyranges as pr
import requests
import os
from io import StringIO
from tqdm import tqdm


def get_hg38_chrom_sizes():
    """
    Fetches hg38 chromosome sizes from UCSC and filters for standard chromosomes.

    Returns:
        dict: A dictionary mapping standard chromosome names to their sizes.
    """
    print("Fetching hg38 chromosome sizes from UCSC...")
    url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes

        chrom_sizes_df = pd.read_csv(StringIO(response.text), sep='\t', header=None, names=['Chromosome', 'Size'])

        # Define standard chromosomes
        standard_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

        # Filter for standard chromosomes and convert to dictionary
        filtered_chrom_sizes = chrom_sizes_df[chrom_sizes_df['Chromosome'].isin(standard_chroms)]

        print(f"Successfully fetched and filtered for {len(filtered_chrom_sizes)} standard chromosomes.")
        return dict(zip(filtered_chrom_sizes.Chromosome, filtered_chrom_sizes.Size))

    except requests.exceptions.RequestException as e:
        print(f"Error fetching chromosome sizes: {e}")
        return None


def create_genomic_bins(chrom_sizes, bin_size=200):
    """
    Generates 200 bp genomic bins for the given chromosomes.

    Args:
        chrom_sizes (dict): Dictionary of chromosome names and their sizes.
        bin_size (int): The size of each genomic bin.

    Returns:
        pd.DataFrame: A DataFrame with 'chr', 'start', and 'end' columns for each bin.
    """
    print(f"Generating {bin_size} bp genomic bins...")
    bins = []
    for chrom, size in tqdm(chrom_sizes.items(), desc="Creating bins for chromosomes"):
        for start in range(0, size, bin_size):
            end = min(start + bin_size, size)
            bins.append({'chr': chrom, 'start': start, 'end': end})

    print(f"Generated a total of {len(bins)} bins.")
    return pd.DataFrame(bins)


def process_chip_seq_files(bins_df, chip_files):
    """
    Processes ChIP-seq files to find overlaps with genomic bins.

    Args:
        bins_df (pd.DataFrame): The DataFrame of genomic bins.
        chip_files (dict): A dictionary mapping column names to ChIP-seq BED file paths.

    Returns:
        pd.DataFrame: The bins DataFrame with added columns indicating overlaps.
    """
    # Create a copy to avoid modifying the original DataFrame in place
    annotated_bins = bins_df.copy()

    # Rename columns for pyranges compatibility
    bins_pr = pr.PyRanges(annotated_bins.rename(columns={'chr': 'Chromosome', 'start': 'Start', 'end': 'End'}))

    for column_name, file_path in tqdm(chip_files.items(), desc="Processing ChIP-seq files"):

        if not os.path.exists(file_path):
            print(f"  WARNING: File not found: {file_path}. Skipping and filling column '{column_name}' with 0.")
            annotated_bins[column_name] = 0
            continue

        try:
            # Read BED file, names argument handles standard BED3 format
            peaks_pr = pr.read_bed(file_path, as_df=True)
            if peaks_pr.empty:
                print(f"  WARNING: File {file_path} is empty. Skipping.")
                annotated_bins[column_name] = 0
                continue

            # Use only the first 3 columns and ensure correct naming for PyRanges
            peaks_pr = peaks_pr.iloc[:, [0, 1, 2]]
            peaks_pr.columns = ['Chromosome', 'Start', 'End']
            peaks_pr = pr.PyRanges(peaks_pr)

            # Find which bins overlap with the peaks
            # The result of join is a PyRanges object containing only the bins that overlap
            overlapping_bins = bins_pr.join(peaks_pr).df

            # Create a unique identifier for each bin (chromosome, start)
            # This is a robust way to mark overlaps
            overlapping_indices = pd.MultiIndex.from_frame(overlapping_bins[['Chromosome', 'Start']])
            original_indices = pd.MultiIndex.from_frame(annotated_bins[['chr', 'start']])

            # Mark bins that have an overlap with 1
            annotated_bins[column_name] = original_indices.isin(overlapping_indices).astype(int)

        except Exception as e:
            print(f"  ERROR: Could not process file {file_path}. Error: {e}")
            print(f"  Skipping and filling column '{column_name}' with 0.")
            annotated_bins[column_name] = 0

    return annotated_bins


def main():
    """
    Main function to execute the full data processing pipeline.
    """
    # --- 1. Define File Paths ---
    poi_file = '../../public/Public Dataset/ChIPseq/Other Cell Lines/H1/H1 CBX8 ChIPseq GSE123227 ENCFF483UZG hg38.bed'

    histone_files = {
        'H3K4me1': '../../public/Public Dataset/ChIPseq/GSE86714 ENCFF991HJA hg38 MCF7 H3K4me1 ChIPseq Filtered_trim.bed',
        'H3K4me2': '../../public/Public Dataset/ChIPseq/GSE96439 ENCFF188VRU hg38 MCF7 H3K4me2 ChIPseq Filtered_trim.bed',
        'H3K4me3': '../../public/Public Dataset/ChIPseq/GSE96506 ENCFF268RXB hg38 MCF7 H3K4me3 ChIPseq Filtered_trim.bed',
        'H3K9ac': '../../public/Public Dataset/ChIPseq/GSE95898 ENCFF348DEB hg38 MCF7 H3K9ac ChIPseq Filtered_trim.bed',
        'H3K9me3': '../../public/Public Dataset/ChIPseq/GSE96517 ENCFF501UHK hg38 MCF7 H3K9me3 ChIPseq Filtered_trim.bed',
        'H3K27ac': '../../public/Public Dataset/ChIPseq/GSE96352 ENCFF491LQY hg38 MCF7 H3K27ac ChIPseq Filtered_trim.bed',
        'H3K27me3': '../../public/Public Dataset/ChIPseq/GSE96363 ENCFF669NUD hg38 MCF7 H3K27me3 ChIPseq Filtered_trim.bed',
        'H3K36me3': '../../public/Public Dataset/ChIPseq/GSE174945 ENCFF195FSD hg38 MCF7 H3K36me3 ChIPseq Filtered_trim.bed',
        'H4K20me1': '../../public/Public Dataset/ChIPseq/GSE96283 ENCFF714DEQ hg38 MCF7 H4K20me1 ChIPseq Filtered_trim.bed'
    }

    # Combine all files into one dictionary for processing
    all_chip_files = {'PoI': poi_file, **histone_files}

    output_path = '../../public/Multiomics/Step18_ML_Prediction_H1/'
    output_file = os.path.join(output_path, 'hg38_bin.csv')

    # --- 2. Get Chromosome Sizes ---
    chrom_sizes = get_hg38_chrom_sizes()
    if not chrom_sizes:
        print("Could not retrieve chromosome sizes. Exiting.")
        return

    # --- 3. Create Genomic Bins ---
    hg38_bin_df = create_genomic_bins(chrom_sizes, bin_size=200)

    # --- 4. Process all ChIP-seq files ---
    final_df = process_chip_seq_files(hg38_bin_df, all_chip_files)

    # --- 5. Export the results ---
    print(f"Exporting the final table to {output_file}...")
    try:
        # Ensure the output directory exists
        os.makedirs(output_path, exist_ok=True)
        final_df.to_csv(output_file, index=False)
        print("Export complete.")
        print("\nFinal table preview:")
        print(final_df.head())
        print(f"\nTotal rows: {len(final_df)}")
        print("\nOverlap summary:")
        print(final_df.iloc[:, 3:].sum())

    except Exception as e:
        print(f"Error exporting file: {e}")


if __name__ == '__main__':
    main()
