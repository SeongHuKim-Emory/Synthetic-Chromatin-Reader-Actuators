#!/usr/bin/env python3

"""
This script loads two CSV files (UpDEGs and Non-DEGs) containing genomic regions,
and for each row, calculates the max and mean signal from a set of BigWig files
for several specified sub-regions (DSPP, CP, USPP, Distal, Enhancer).

The processing is parallelized using concurrent.futures.ProcessPoolExecutor.

*** Optimization Notes ***
This version is optimized for speed by:
1.  Using an `initializer` function for the ProcessPoolExecutor.
2.  Each worker process opens all BigWig files ONCE and stores the
    file handles in a global variable (OPEN_BW_HANDLES).
3.  The row processing function (`process_row`) accesses these
    already-open file handles, eliminating millions of expensive
    file open/close operations.
4.  Replaced slow `df.apply(axis=1)` with fast `df.itertuples()`.
5.  Removed pre-allocation of columns, as it's faster to build
    a new DataFrame from results and join it.

Python Version: 3.12.3
Required Libraries: pandas, pyBigWig, numpy, concurrent.futures, os, tqdm
"""

import pandas as pd
import pyBigWig
import numpy as np
import os
import concurrent.futures
import time
from tqdm import tqdm  # <-- ADDED FOR PROGRESS BAR

# =============================================================================
# Configuration
# =============================================================================

print("Step 1: Initializing Configuration...")

# A. Input CSV Files
UpDEG_CSV = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Pro.csv'
Non_DEG_CSV = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Pro.csv'

# B. BigWig File Configuration
# I. Primary .bigWig file
bw_file_025 = {
    "025": "../../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2_ppois.bigWig"
}

# II. Other histone mark and CTCF .bigWig files
other_bw_files_list = [
    "../../public/Public Dataset/ChIPseq/MCF7 ARID3A ChIPseq GSE91595 ENCFF907DWT hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ATF7 ChIPseq GSE106025 ENCFF703ZRP hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 BMI1 ChIPseq GSE105933 ENCFF977TJK hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 CEBPB ChIPseq GSM1010889 ENCFF464JOC hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 CHD1 ChIPseq GSE91698 ENCFF923MDX hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 CLOCK ChIPseq GSE127438 ENCFF045RXD hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 COPS2 ChIPseq GSE105356 ENCFF102ABZ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 CREB1 ChIPseq GSE105525 ENCFF954WHX hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 CSDE1 ChIPseq GSE105291 ENCFF722GUE hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 CTBP1 ChIPseq GSE91938 ENCFF922FKY hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 CTCF ChIPseq GSE123219 ENCFF507CRU hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 CUX1 ChIPseq GSE91415 ENCFF663YXL hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 DDX20 ChIPseq GSE105517 ENCFF895IXP hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 DPF2 ChIPseq GSE91598 ENCFF142KHL hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 E2F8 ChIPseq GSE127615 ENCFF919YCR hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 E4F1 ChIPseq GSE127582 ENCFF555HAI hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 EGR1 ChIPseq GSM1010844 ENCFF487DKA hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ELF1 ChIPseq GSE105503 ENCFF751TAA hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ELK1 ChIPseq GSE91713 ENCFF045XJM hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 EP300 ChIPseq GSM1010800 ENCFF708NMR hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ESRRA ChIPseq GSE92187 ENCFF059ZSC hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 FOS ChIPseq GSE105734 ENCFF327LTA hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 FOSL2 ChIPseq GSM1010768 ENCFF704HAX hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 FOXA1 ChIPseq GSE105305 ENCFF512UGW hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 FOXK2 ChIPseq GSE91785 ENCFF477VZD hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 FOXM1 ChIPseq GSM1010769 ENCFF405RBL hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 GABPA ChIPseq GSM1010864 ENCFF494GAR hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 GATA3 ChIPseq GSE127656 ENCFF971AZB hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 GATAD2B ChIPseq GSE91723 ENCFF560WGB hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 GTF2F1 ChIPseq GSE91869 ENCFF736KUY hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 H3K4me1 ChIPseq GSE86714 ENCFF763NCP hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 H3K4me2 ChIPseq GSE96439 ENCFF442RRY hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 H3K4me3 ChIPseq GSE96506 ENCFF163MXP hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 H3K9ac ChIPseq GSE95898 ENCFF327XJC hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 H3K9me3 ChIPseq GSE96517 ENCFF481DZL hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 H3K27ac ChIPseq GSE96352 ENCFF138YNG hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 H3K27me3 ChIPseq GSE96363 ENCFF163QKN hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 H3K36me3 ChIPseq GSE174945 ENCFF910BRP hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 H4K20me1 ChIPseq GSE96283 ENCFF366GLZ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 HCFC1 ChIPseq GSE91992 ENCFF103ABY hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 HDAC2 ChIPseq GSM1010825 ENCFF178WAQ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 HDGF ChIPseq GSE91567 ENCFF178LPH hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 HES1 ChIPseq GSE105287 ENCFF509HTJ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 HSF1 ChIPseq GSE91444 ENCFF505ZBM hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 JUN ChIPseq GSE91550 ENCFF730TVS hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 JUND ChIPseq GSM1010892 ENCFF990FGN hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 LARP7 ChIPseq GSE105864 ENCFF148SFF hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 MAFK ChIPseq GSE127366 ENCFF174HAP hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 MAX ChIPseq GSM1010863 ENCFF315FIO hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 MAZ ChIPseq GSE91633 ENCFF297WGG hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 MBD2 ChIPseq GSE96480 ENCFF517UQD hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 MLLT1 ChIPseq GSE91754 ENCFF654ISN hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 MNT ChIPseq GSE91968 ENCFF623IRB hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 MTA1 ChIPseq GSE91687 ENCFF883HBF hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 MTA2 ChIPseq GSE91864 ENCFF557XGU hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 MTA3 ChIPseq GSE91727 ENCFF486PVE hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 NBN ChIPseq GSE91899 ENCFF242QKO hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 NCOA3 ChIPseq GSE105681 ENCFF770MLI hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 NEUROD1 ChIPseq GSE127346 ENCFF501HJF hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 NFIB ChIPseq GSE105751 ENCFF631HNZ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 NFRKB ChIPseq GSE105388 ENCFF411WFM hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 NFXL1 ChIPseq GSE105498 ENCFF537TBX hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 NONO ChIPseq GSE92159 ENCFF875CIS hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 NR2F2 ChIPseq GSM1010837 ENCFF678MPN hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 NRF1 ChIPseq GSE91522 ENCFF113MFU hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 PAX8 ChIPseq GSE105903 ENCFF775POG hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 PKNOX1 ChIPseq GSE92210 ENCFF737NUU hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 PML ChIPseq GSM1010838 ENCFF911FVQ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 POLR2A ChIPseq GSM822295 ENCFF827YIP hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 PPP1R10 ChIPseq GSE105558 ENCFF167JMR hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 RAD21 ChIPseq GSM1010791 ENCFF775EKJ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 RAD51 ChIPseq GSE105597 ENCFF522WFI hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 RCOR1 ChIPseq GSE91726 ENCFF536KZN hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 REST ChIPseq GSM1010891 ENCFF068TEY hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 RFX1 ChIPseq GSE91448 ENCFF310ZPD hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 RFX5 ChIPseq GSE105874 ENCFF085XIV hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 SIN3A ChIPseq GSE91789 ENCFF807QKK hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 SIX4 ChIPseq GSE91630 ENCFF032MSQ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 SMARCA5 ChIPseq GSE105722 ENCFF297JNO hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 SMARCE1 ChIPseq GSE105546 ENCFF375LRP hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 SNIP1 ChIPseq GSE105188 ENCFF231LYV hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 SP1 ChIPseq GSE92014 ENCFF367HVY hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 SREBF1 ChIPseq GSE91561 ENCFF791DFW hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 SRF ChIPseq GSM1010839 ENCFF534AJF hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 SUZ12 ChIPseq GSE105981 ENCFF598JRC hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 TAF1 ChIPseq GSM1010811 ENCFF418WWV hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 TARDBP ChIPseq GSE105812 ENCFF065RQM hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 TCF7L2 ChIPseq GSM816438 ENCFF711CKB hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 TCF12 ChIPseq GSM1010861 ENCFF031XVC hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 TEAD4 ChIPseq GSM1010860 ENCFF977ZGZ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 TOE1 ChIPseq GSE105831 ENCFF318YIK hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 TRIM22 ChIPseq GSE127607 ENCFF100OCC hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 YBX1 ChIPseq GSE92114 ENCFF033VON hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZBTB1 ChIPseq GSE105482 ENCFF062ZHZ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZBTB7B ChIPseq GSE105418 ENCFF073COQ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZBTB11 ChIPseq GSE95952 ENCFF510SQD hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZBTB33 ChIPseq GSE91596 ENCFF506CTQ hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZBTB40 ChIPseq GSE91661 ENCFF539ZHK hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZFX ChIPseq GSE105562 ENCFF482RYD hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZHX2 ChIPseq GSE96441 ENCFF724GPW hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZKSCAN1 ChIPseq GSE91769 ENCFF239BWH hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF8 ChIPseq GSE127545 ENCFF550FAH hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF24 ChIPseq GSE105531 ENCFF505ORH hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF207 ChIPseq GSE91475 ENCFF668FJE hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF217 ChIPseq GSE105662 ENCFF491CXP hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF217 ChIPseq GSE127625 ENCFF343WRL hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF444 ChIPseq GSE127367 ENCFF101ION hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF507 ChIPseq GSE127652 ENCFF659LTF hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF512B ChIPseq GSE105987 ENCFF413BLF hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF574 ChIPseqGSE127638 ENCFF914VHT hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF579 ChIPseq GSE105224 ENCFF913XOW hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF592 ChIPseq GSE91425 ENCFF950SDY hg38.bigWig",
    "../../public/Public Dataset/ChIPseq/MCF7 ZNF687 ChIPseq GSE92150 ENCFF910XUA hg38.bigWig",
    "../../public/Public Dataset/DNase_seq/MCF7 DNase_seq GSM1008565 ENCFF137FTO hg38.bigWig",
    "../../public/Public Dataset/ATACseq/MCF7 ATACseq GSE169929 ENCFF976UNK hg38.bigWig"
]

# III. Dynamically build the final bigwig_files dictionary
bigwig_files = {}
bigwig_files.update(bw_file_025)

for file_path in other_bw_files_list:
    try:
        base_filename = os.path.basename(file_path)
        # Split by space, the mark name is the second word (index 1)
        mark_name = base_filename.split(' ')[1]
        bigwig_files[mark_name] = file_path
    except IndexError:
        print(f"Warning: Could not parse mark name from filename: {file_path}")
    except Exception as e:
        print(f"Error processing filename {file_path}: {e}")

print(f"Successfully configured {len(bigwig_files)} BigWig files:")
for mark, path in bigwig_files.items():
    print(f"  - {mark}: {os.path.basename(path)}")

# C. Genomic Region Configuration
regions = [
    ("DSPP", "DSPP_start", "DSPP_end"),
    ("CP", "CP_start", "CP_end"),
    ("USPP", "USPP_start", "USPP_end"),
    ("Distal", "Distal_start", "Distal_end"),
    ("Enhancer", "Enhancer_start", "Enhancer_end")
]

print(f"Configured {len(regions)} regions to process: {[r[0] for r in regions]}")

# D. Output File Configuration
Output_UpDEG_CSV = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone.csv'
Output_Non_DEG_CSV = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone.csv'

# =============================================================================
# Helper Functions for Parallel Processing
# =============================================================================

# These will be 'global' within each worker process
OPEN_BW_HANDLES = {}
REGIONS_CONFIG = []


def init_worker(bw_files_from_main, regions_config_from_main):
    """
    Initializer function for each worker process.
    Opens all BigWig files ONCE and stores them in a global variable.
    """
    global OPEN_BW_HANDLES, REGIONS_CONFIG
    REGIONS_CONFIG = regions_config_from_main

    # Suppress worker-level print statements for a cleaner progress bar
    # print(f"Worker {os.getpid()}: Initializing and opening {len(bw_files_from_main)} BigWig files...")
    for mark_name, file_path in bw_files_from_main.items():
        try:
            OPEN_BW_HANDLES[mark_name] = pyBigWig.open(file_path)
        except Exception as e:
            # Log error if a file can't be opened, and skip it
            print(f"WORKER_INIT_ERROR (pid {os.getpid()}): Failed to open {file_path}. Error: {e}")
    # print(f"Worker {os.getpid()}: Initialization complete. {len(OPEN_BW_HANDLES)} files open.")


def process_row(row):
    """
    Processes a single row (as a namedtuple from itertuples).
    Calculates max and mean for all configured regions using the
    pre-opened BigWig file handles from the worker's global scope.
    Returns a dictionary of the new values.
    """
    # 'row' is a namedtuple, access attributes by name
    try:
        chr_val = row.chr
    except AttributeError:
        # Handle rows without a 'chr' column
        return {}  # Return empty dict, will be filtered out

    if not chr_val or pd.isna(chr_val):
        return {}  # Cannot process without a chromosome

    results = {}

    for region_prefix, start_col, end_col in REGIONS_CONFIG:
        try:
            # Get region coordinates using getattr for dynamic column names
            region_start = getattr(row, start_col)
            region_end = getattr(row, end_col)

            # Ensure coordinates are valid
            if pd.isna(region_start) or pd.isna(region_end):
                continue  # Skip this region if coordinates are missing

            region_start = int(region_start)
            region_end = int(region_end)

            if region_start >= region_end:
                continue  # Skip if region is invalid (start >= end)

        except (ValueError, TypeError, AttributeError):
            # Handle cases where coords are not numeric or cols are missing
            continue  # Skip this region

        # Inside region loop, iterate through each PRE-OPENED BigWig file
        for mark_name, bw_handle in OPEN_BW_HANDLES.items():
            max_val = 0.0
            mean_val = 0.0

            try:
                # Check if chromosome exists in the BigWig file
                if bw_handle.chroms().get(chr_val) is None:
                    continue  # Skip if chromosome is not found

                # Check if range is valid for the chromosome
                chrom_len = bw_handle.chroms().get(chr_val)
                if region_start < 0 or region_end > chrom_len or region_start > chrom_len:
                    continue  # Skip invalid ranges

                # Calculate max value
                max_stats = bw_handle.stats(chr_val, region_start, region_end, type="max")
                if max_stats and max_stats[0] is not None:
                    max_val = max_stats[0]

                # Calculate mean value
                mean_stats = bw_handle.stats(chr_val, region_start, region_end, type="mean")
                if mean_stats and mean_stats[0] is not None:
                    mean_val = mean_stats[0]

            except Exception as e:
                # Log error but continue processing
                if "no data" not in str(e):
                    # Suppress most warnings to keep progress bar clean
                    # You can re-enable this for debugging
                    pass
                    # print(
                    #     f"Warning (pid {os.getpid()}): Error processing {mark_name} for {chr_val}:{region_start}-{region_end}. Error: {e}")

            # Store results in the dictionary
            results[f"Max_{region_prefix}_{mark_name}"] = max_val
            results[f"Mean_{region_prefix}_{mark_name}"] = mean_val

    return results


def process_dataframe_chunk(df_chunk):
    """
    Applies the 'process_row' function to each row in a DataFrame chunk
    using the highly efficient `itertuples()`.

    This function is mapped to parallel processes.

    Returns a DataFrame containing the *original* columns plus the
    *newly processed* columns.
    """
    # itertuples() is much faster than apply(axis=1)
    # It yields namedtuples, which 'process_row' is designed to accept
    results_list = [process_row(row) for row in df_chunk.itertuples()]

    # Create a new DataFrame from the list of result dictionaries
    # Use the original chunk's index to ensure correct alignment
    new_cols_df = pd.DataFrame(results_list, index=df_chunk.index)

    # Join the new columns back to the original chunk
    # This is much faster than pre-allocating and modifying
    return df_chunk.join(new_cols_df)


# =============================================================================
# Main Execution
# =============================================================================

def main():
    start_time = time.time()

    # --- Step 2: Table Initialization ---
    print("Step 2: Loading input CSV files...")
    try:
        df_updeg = pd.read_csv(UpDEG_CSV)
        df_non_deg = pd.read_csv(Non_DEG_CSV)
        print(f"  Loaded UpDEG data: {df_updeg.shape[0]} rows")
        print(f"  Loaded Non_DEG data: {df_non_deg.shape[0]} rows")
    except FileNotFoundError as e:
        print(f"Error: Input file not found. {e}")
        return
    except Exception as e:
        print(f"Error loading CSV files: {e}")
        return

    dataframes_to_process = [df_updeg, df_non_deg]
    df_names = ["UpDEG", "Non_DEG"]

    # --- Step 3: (REMOVED) ---
    # We no longer pre-allocate columns. It's faster to build
    # the new columns and join them at the end of chunk processing.
    print("Step 3: Skipping column pre-allocation (optimization).")

    # --- Step 4: Process BigWig Data (Parallelized) ---
    print("Step 4: Starting parallel processing of BigWig data...")

    # Determine number of workers
    num_workers = os.cpu_count()
    print(f"  Using {num_workers} processes.")

    processed_dataframes = []

    # Pass the file dictionary and region config to the initializer
    init_args = (bigwig_files, regions)

    with concurrent.futures.ProcessPoolExecutor(
            max_workers=num_workers,
            initializer=init_worker,
            initargs=init_args) as executor:

        for df, name in zip(dataframes_to_process, df_names):
            if df.empty:
                print(f"  DataFrame {name} is empty, skipping.")
                processed_dataframes.append(df)
                continue

            print(f"  Submitting {name} DataFrame ({df.shape[0]} rows) for processing...")
            df_start_time = time.time()

            # Split DataFrame into chunks for parallel processing
            # Using 2-4 chunks per worker is often a good balance
            num_chunks = min(num_workers * 4, df.shape[0])
            if num_chunks == 0:
                processed_dataframes.append(df)  # Handle tiny dataframes
                continue

            # np.array_split is efficient for splitting DataFrames
            df_chunks = np.array_split(df, num_chunks)

            # Map the processing function to the chunks
            try:
                # --- MODIFICATION: Added tqdm for progress bar ---
                # Wrap the executor.map iterator with tqdm
                # total=num_chunks sets the max value for the progress bar
                # desc=... adds a label (e.g., "  Processing UpDEG chunks")
                # list() collects all results, ensuring all chunks are done
                results = list(tqdm(
                    executor.map(process_dataframe_chunk, df_chunks),
                    total=num_chunks,
                    desc=f"  Processing {name} chunks",
                    unit="chunk"
                ))
                # --- END MODIFICATION ---

                # Combine results back into a single DataFrame
                processed_df = pd.concat(results)
                processed_dataframes.append(processed_df)

                df_end_time = time.time()
                print(f"  Finished processing {name} DataFrame in {df_end_time - df_start_time:.2f} seconds.")

            except Exception as e:
                print(f"Error during parallel processing of {name}: {e}")
                # Append original df to avoid crashing the script
                processed_dataframes.append(df)

    # Ensure we have the correct dataframes, even if one was empty/failed
    df_updeg_processed = processed_dataframes[0] if len(processed_dataframes) > 0 else pd.DataFrame()
    df_non_deg_processed = processed_dataframes[1] if len(processed_dataframes) > 1 else pd.DataFrame()

    # --- Step 5: Export Results ---
    print("Step 5: Exporting processed DataFrames to CSV...")

    try:
        # Ensure output directories exist
        os.makedirs(os.path.dirname(Output_UpDEG_CSV), exist_ok=True)
        os.makedirs(os.path.dirname(Output_Non_DEG_CSV), exist_ok=True)

        # Export UpDEG
        if not df_updeg_processed.empty:
            df_updeg_processed.to_csv(Output_UpDEG_CSV, index=False)
            print(f"  Successfully saved processed UpDEG data to: {Output_UpDEG_CSV}")
        else:
            print("  UpDEG DataFrame is empty, nothing to save.")

        # Export Non_DEG
        if not df_non_deg_processed.empty:
            df_non_deg_processed.to_csv(Output_Non_DEG_CSV, index=False)
            print(f"  Successfully saved processed Non_DEG data to: {Output_Non_DEG_CSV}")
        else:
            print("  Non_DEG DataFrame is empty, nothing to save.")

    except Exception as e:
        print(f"Error exporting CSV files: {e}")

    end_time = time.time()
    print(f"\nTotal execution time: {end_time - start_time:.2f} seconds.")


# This guard is crucial for multiprocessing to work correctly
if __name__ == "__main__":
    main()