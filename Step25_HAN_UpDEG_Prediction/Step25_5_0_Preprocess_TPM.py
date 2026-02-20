#!/usr/bin/env python3
"""
This script preprocesses an RNA-seq TPM file using Polars.

It performs the following steps:
1. Imports the specified CSV file.
2. Creates a table (DataFrame) from it.
3. Removes rows where 'gene_symbol' is 'N/A'.
4. Removes rows where 'gene_id' contains an underscore ('_').
5. Saves the cleaned data to a new CSV file.

This script is designed to be efficient and will automatically utilize
multiple CPU cores (like your 36 CPUs) for filtering and I/O.
"""

import polars as pl
import os
import sys

def preprocess_tpm_data(input_path, output_path):
    """
    Reads, filters, and saves TPM data using Polars.
    """
    try:
        # --- Step 1 & 2: Import file and create table (DataFrame) ---
        # Polars' read_csv is multi-threaded and very fast.
        print(f"Step 1/2: Reading file into 'TPM' table from: {input_path}")
        
        # Check if file exists before reading
        if not os.path.exists(input_path):
            print(f"Error: Input file not found at {input_path}", file=sys.stderr)
            return

        tpm_df = pl.read_csv(input_path)
        print(f"  > Original shape (rows, cols): {tpm_df.shape}")

        # --- Step 3 & 4: Filter the TPM table ---
        # We can chain both filtering conditions together for maximum
        # efficiency. Polars will run these operations in parallel.
        print("Step 3/4: Applying filters...")
        
        # Condition 3: 'gene_symbol' is not "N/A"
        condition_3 = pl.col("gene_symbol") != "N/A"
        
        # Condition 4: 'gene_id' does not contain "_"
        condition_4 = ~pl.col("gene_id").str.contains("_")

        # Apply both conditions at once
        tpm_processed_df = tpm_df.filter(condition_3 & condition_4)
        
        print(f"  > Shape after filtering (rows, cols): {tpm_processed_df.shape}")

        # --- Step 5: Save the processed TPM table ---
        # Ensure the output directory exists
        output_dir = os.path.dirname(output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"  > Created output directory: {output_dir}")

        print(f"Step 5: Saving processed file to: {output_path}")
        tpm_processed_df.write_csv(output_path)
        
        print("\nPreprocessing complete.")

    except pl.exceptions.ComputeError as e:
        print(f"An error occurred during Polars computation: {e}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Could not find the input file at {input_path}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)

def main():
    # Define file paths
    # Note: Using os.path.join for better cross-platform compatibility
    base_dir = os.path.join("..", "..", "public", "Multiomics", "Step19_ML_UpDEG_Prediction")
    
    input_file = os.path.join(
        base_dir,
        "GSE175204 ENCFF721BRA hg38 MCF7 total RNAseq ENSG Only Symbol and TPM.csv"
    )
    
    output_file = os.path.join(
        base_dir,
        "MCF7_RNAseq_TPM_Preprocessed_GSE175204_ENCFF721BRA_hg38.csv"
    )
    
    # Run the preprocessing function
    preprocess_tpm_data(input_file, output_file)

if __name__ == "__main__":
    # Polars automatically detects and uses available CPU cores.
    # No special setup is needed for your 36 CPUs.
    main()