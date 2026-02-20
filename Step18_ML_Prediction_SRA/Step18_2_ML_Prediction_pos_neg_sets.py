
import pandas as pd
import numpy as np
from pybedtools import BedTool
import tempfile
import os
from Bio import SeqIO
import random
import requests

def calculate_gc_content(sequence):
    """Calculate GC content manually"""
    if len(sequence) == 0:
        return 50.0
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def fetch_hg38_chromosome_sizes():
    """
    Fetch up-to-date hg38 chromosome sizes from UCSC
    Returns a dictionary with chromosome names as keys and sizes as values
    """

    # UCSC provides the chromosome sizes file
    ucsc_url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes"

    try:
        print("Fetching chromosome sizes from UCSC...")
        response = requests.get(ucsc_url, timeout=30)
        response.raise_for_status()

        # Parse the tab-delimited file
        genome_sizes = {}
        lines = response.text.strip().split('\n')

        for line in lines:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 2:
                    chrom = parts[0]
                    try:
                        size = int(parts[1])
                        # Only keep main chromosomes (chr1-22, chrX, chrY, chrM)
                        if (chrom.startswith('chr') and
                            (chrom[3:].isdigit() or chrom[3:] in ['X', 'Y', 'M'])):
                            genome_sizes[chrom] = size
                    except ValueError:
                        continue

        print(f"Successfully fetched sizes for {len(genome_sizes)} chromosomes")
        return genome_sizes

    except Exception as e:
        print(f"Error fetching from UCSC: {e}")
        print("Falling back to hardcoded values...")
        # Fallback to hardcoded values if network fails
        return {
            'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
            'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
            'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
            'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
            'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
            'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
        }

def find_gc_matched_negatives_fast(positive_df, negative_df, tolerance=5.0, ratio=10):
    """Fast GC matching using binning approach"""

    print(f"Fast GC matching: {len(positive_df)} positive vs {len(negative_df)} negative candidates")

    # Create GC content bins (2% bins)
    gc_bins = np.arange(0, 101, 2)
    positive_df['gc_bin'] = pd.cut(positive_df['gc_content'], bins=gc_bins, labels=False, include_lowest=True)
    negative_df['gc_bin'] = pd.cut(negative_df['gc_content'], bins=gc_bins, labels=False, include_lowest=True)

    matched_negatives = []

    # For each GC bin, sample proportionally
    for bin_idx in positive_df['gc_bin'].unique():
        if pd.isna(bin_idx):
            continue

        pos_in_bin = positive_df[positive_df['gc_bin'] == bin_idx]
        neg_in_bin = negative_df[negative_df['gc_bin'] == bin_idx]

        needed = len(pos_in_bin) * ratio
        available = len(neg_in_bin)

        if available > 0:
            n_sample = min(needed, available)
            sampled = neg_in_bin.sample(n=n_sample, random_state=42)
            matched_negatives.append(sampled)
            print(f"GC bin {bin_idx*2:.0f}-{(bin_idx+1)*2:.0f}%: {len(pos_in_bin)} pos -> {n_sample} neg")

    if matched_negatives:
        return pd.concat(matched_negatives, ignore_index=True)
    else:
        return pd.DataFrame()

def main():
    print("Starting the analysis pipeline...")

    # 1. Import files
    print("Step 1: Importing files...")

    # Import the main ChIP-seq data
    chip_data = pd.read_csv("../../public/Multiomics/Step18_ML_Prediction_SRA/hg38_bin.csv")
    print(f"Loaded ChIP-seq data with {len(chip_data)} rows")

    # Import blacklist regions
    blacklist = pd.read_csv("../../public/Public Dataset/ChIPseq/ENCSR636HFF GRCh38_unified_blacklist.bed",
                           sep='\t', header=None, names=['chr', 'start', 'end'])
    print(f"Loaded blacklist with {len(blacklist)} regions")

    # Import reference genome
    print("Loading reference genome...")
    genome_file = "../../public/Multiomics/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    genome_dict = {}
    for record in SeqIO.parse(genome_file, "fasta"):
        # Only keep main chromosomes (chr1-22, X, Y)
        if record.id in [f'{i}' for i in range(1, 23)] + ['X', 'Y']:
            genome_dict[f'chr{record.id}'] = str(record.seq)
            print(f"Loaded chromosome {record.id}")

    print(f"Loaded {len(genome_dict)} chromosomes")

    # 2. Create the Positive Set
    print("\nStep 2: Creating positive set...")
    positive_set = chip_data[chip_data['PoI'] == 1].copy()
    print(f"Positive set contains {len(positive_set)} regions")

    # 3. Create Negative Set using pybedtools
    print("\nStep 3: Creating negative set...")

    # Create BedTool objects
    positive_bed = BedTool.from_dataframe(positive_set[['chr', 'start', 'end']])
    blacklist_bed = BedTool.from_dataframe(blacklist)

    # Add buffer around positive regions (1kb = 1000bp)
    buffer_size = 1000

    # Fetch up-to-date chromosome sizes from UCSC
    genome_sizes = fetch_hg38_chromosome_sizes()

    # Create temporary genome file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.genome') as f:
        for chrom, size in genome_sizes.items():
            f.write(f"{chrom}\t{size}\n")
        genome_file_path = f.name

    # Add buffer around positive regions
    positive_buffered = positive_bed.slop(b=buffer_size, g=genome_file_path)

    # Combine positive regions (with buffer) and blacklist regions
    excluded_regions = positive_buffered.cat(blacklist_bed, postmerge=False).sort().merge()

    # Generate random regions excluding the excluded regions
    # We want 15x more negative regions than positive for better GC matching
    num_negative = len(positive_set) * 15

    print(f"Generating {num_negative} random regions...")
    random_regions = BedTool().random(l=200, n=num_negative, g=genome_file_path).subtract(excluded_regions)

    # Convert to dataframe and handle variable number of columns
    random_df = random_regions.to_dataframe()
    print(f"Random regions dataframe shape: {random_df.shape}")
    print(f"Random regions columns: {list(random_df.columns)}")

    # Take only the first 3 columns (chr, start, end) regardless of how many columns bedtools returns
    random_df = random_df.iloc[:, :3]  # Select first 3 columns
    random_df.columns = ['chr', 'start', 'end']

    print(f"Generated {len(random_df)} candidate negative regions")

    # Function to calculate GC content for a genomic region
    def get_gc_content(chrom, start, end, genome_dict):
        try:
            if chrom in genome_dict:
                seq = genome_dict[chrom][start:end]
                if len(seq) > 0:
                    return calculate_gc_content(seq)
                else:
                    return 50.0
            else:
                return 50.0  # Default GC content if chromosome not found
        except Exception as e:
            print(f"Error calculating GC content for {chrom}:{start}-{end}: {e}")
            return 50.0

    # Calculate GC content for positive regions
    print("Calculating GC content for positive regions...")
    positive_set['gc_content'] = positive_set.apply(
        lambda row: get_gc_content(row['chr'], row['start'], row['end'], genome_dict),
        axis=1
    )

    # Calculate GC content for negative candidates
    print("Calculating GC content for negative candidates...")
    random_df['gc_content'] = random_df.apply(
        lambda row: get_gc_content(row['chr'], row['start'], row['end'], genome_dict),
        axis=1
    )

    # Match GC content between positive and negative sets using fast approach
    print("Matching GC content between positive and negative sets...")
    negative_set = find_gc_matched_negatives_fast(positive_set, random_df, tolerance=5.0, ratio=10)

    print(f"Final negative set contains {len(negative_set)} regions")
    print(f"Ratio of negative to positive: {len(negative_set) / len(positive_set):.1f}:1")

    # Add histone mark columns to negative set (all zeros)
    histone_marks = ['H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9ac', 'H3K9me3',
                     'H3K27ac', 'H3K27me3', 'H3K36me3', 'H4K20me1']

    for mark in histone_marks:
        if mark not in negative_set.columns:
            negative_set[mark] = 0

    # Add PoI column (0 for negative)
    negative_set['PoI'] = 0

    # Reorder columns to match positive set
    column_order = ['chr', 'start', 'end', 'PoI'] + histone_marks
    negative_set = negative_set[column_order]

    # 4. Export positive set
    print("\nStep 4: Exporting positive set...")
    os.makedirs("../../public/Multiomics", exist_ok=True)
    positive_set_export = positive_set[column_order]  # Remove gc_content column
    positive_set_export.to_csv("../../public/Multiomics/Step18_ML_Prediction_SRA/hg38_bin_positive.csv", index=False)
    print(f"Exported {len(positive_set_export)} positive regions to hg38_bin_positive.csv")

    # 5. Export negative set
    print("Step 5: Exporting negative set...")
    negative_set.to_csv("../../public/Multiomics/Step18_ML_Prediction_SRA/hg38_bin_negative.csv", index=False)
    print(f"Exported {len(negative_set)} negative regions to hg38_bin_negative.csv")

    # Clean up temporary files
    os.unlink(genome_file_path)

    # Print summary statistics
    print("\nSummary:")
    print(f"Positive regions: {len(positive_set_export)}")
    print(f"Negative regions: {len(negative_set)}")
    print(f"Ratio (negative:positive): {len(negative_set) / len(positive_set_export):.1f}:1")

    print("\nGC content statistics:")
    if 'gc_content' in positive_set.columns and 'gc_content' in negative_set.columns:
        print(f"Positive set GC content: {positive_set['gc_content'].mean():.2f}% ± {positive_set['gc_content'].std():.2f}%")
        print(f"Negative set GC content: {negative_set['gc_content'].mean():.2f}% ± {negative_set['gc_content'].std():.2f}%")

    print("\nChromosome sizes used:")
    for chrom, size in sorted(genome_sizes.items(), key=lambda x: (len(x[0]), x[0])):
        print(f"  {chrom}: {size:,} bp")

    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()
