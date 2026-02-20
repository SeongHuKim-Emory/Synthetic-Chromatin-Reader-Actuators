#!/usr/bin/env python3
#
# Python 3.12.3
#
# This script processes gene node data from two input files, merges them,
# adds classification labels, and formats the data into a final 'Nodes.csv'
# file for use in a Graph Neural Network (GNN).
#
# This script uses the pandas library, which is the industry standard for
# this type of tabular data manipulation. While 25 CPUs and CUDA were
# mentioned, they are not applicable or efficient for this workflow.
# This I/O-bound task (reading/writing CSVs) is most efficiently handled
# by a direct pandas implementation.
#

import pandas as pd
from pathlib import Path
import sys
import time

def main():
    """
    Main function to execute the full data processing pipeline.
    """
    print("Starting data processing pipeline...")
    start_time = time.time()

    # Define file paths using pathlib for cross-platform compatibility
    base_input_dir = Path("../../public/Multiomics/Step19_ML_UpDEG_Prediction")
    output_dir = Path("../../public/Multiomics/Step24_GNN_UpDEG_Prediction")

    # Input file paths
    updeg_file = base_input_dir / "GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Pro.csv"
    nondeg_file = base_input_dir / "GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Pro.csv"

    # Ensure output directory exists (Step 10, 11, 12, 17)
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Ensured output directory exists: {output_dir}")
    except Exception as e:
        print(f"Error creating output directory: {e}", file=sys.stderr)
        sys.exit(1)

    # Columns to be removed (Step 4)
    cols_to_drop = [
        'element_start', 'element_end', 'gene_chr', 'gene_start', 'gene_end',
        'gene_strand', '025_on_element', 'Element_TSS_Distance', 'DSPP_start',
        'DSPP_end', 'CP_start', 'CP_end', 'USPP_start', 'USPP_end',
        'Distal_start', 'Distal_end', 'Enhancer_start', 'Enhancer_end'
    ]

    try:
        # Step 1 & 2: Import UpDEG file
        print(f"Loading UpDEG data from {updeg_file}...")
        SRA_Specific_UpDEG_Gene_Nodes = pd.read_csv(updeg_file)
        print(f"Loaded {len(SRA_Specific_UpDEG_Gene_Nodes)} UpDEG records.")

        # Step 1 & 3: Import Non-DEG file
        print(f"Loading Non-DEG data from {nondeg_file}...")
        Non_DEG_Gene_Nodes = pd.read_csv(nondeg_file)
        print(f"Loaded {len(Non_DEG_Gene_Nodes)} Non-DEG records.")

    except FileNotFoundError as e:
        print(f"Error: Input file not found. {e}", file=sys.stderr)
        print("Please ensure the file paths are correct.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)

    # --- Start Processing ---

    # Step 4: Remove columns
    # We filter the drop list by columns that actually exist in the DataFrame
    # to prevent KeyErrors if the script is re-run or columns change.
    print("Step 4: Removing specified columns...")
    updeg_cols_exist = [col for col in cols_to_drop if col in SRA_Specific_UpDEG_Gene_Nodes.columns]
    SRA_Specific_UpDEG_Gene_Nodes = SRA_Specific_UpDEG_Gene_Nodes.drop(columns=updeg_cols_exist)

    nondeg_cols_exist = [col for col in cols_to_drop if col in Non_DEG_Gene_Nodes.columns]
    Non_DEG_Gene_Nodes = Non_DEG_Gene_Nodes.drop(columns=nondeg_cols_exist)

    # Step 5: Add 'is_upregulated' = 1 to UpDEG nodes
    print("Step 5: Adding 'is_upregulated' = 1 to UpDEG nodes...")
    SRA_Specific_UpDEG_Gene_Nodes['is_upregulated'] = 1

    # Step 6: Add 'is_upregulated' = 0 to Non-DEG nodes
    print("Step 6: Adding 'is_upregulated' = 0 to Non-DEG nodes...")
    Non_DEG_Gene_Nodes['is_upregulated'] = 0

    # Step 7 & 8: Merge and sort
    print("Step 7 & 8: Merging tables and sorting by 'symbol'...")
    Gene_Nodes = pd.concat(
        [SRA_Specific_UpDEG_Gene_Nodes, Non_DEG_Gene_Nodes],
        ignore_index=True
    )
    Gene_Nodes = Gene_Nodes.sort_values(by='symbol')
    Gene_Nodes = Gene_Nodes.reset_index(drop=True)

    # Step 9: Add 'node_id' column
    print("Step 9: Adding 'node_id' column...")
    Gene_Nodes.insert(0, 'node_id', Gene_Nodes.index)

    # --- Export Intermediate Files ---

    # Step 10: Export SRA_Specific_UpDEG_Gene_Nodes
    path_10 = output_dir / "SRA_Specific_UpDEG_Gene_Nodes.csv"
    print(f"Step 10: Exporting to {path_10}...")
    SRA_Specific_UpDEG_Gene_Nodes.to_csv(path_10, index=False)

    # Step 11: Export Non_DEG_Gene_Nodes
    path_11 = output_dir / "Non_DEG_Gene_Nodes.csv"
    print(f"Step 11: Exporting to {path_11}...")
    Non_DEG_Gene_Nodes.to_csv(path_11, index=False)

    # Step 12: Export Gene_Nodes
    path_12 = output_dir / "Gene_Nodes.csv"
    print(f"Step 12: Exporting to {path_12}...")
    Gene_Nodes.to_csv(path_12, index=False)

    # --- Create Final Nodes Table ---

    # Step 13: Create Nodes table
    print("Step 13: Creating final 'Nodes' table...")
    Nodes = Gene_Nodes.copy()

    # Step 14: Remove 'is_upregulated' and add 'node_type'
    print("Step 14: Modifying 'Nodes' table structure...")
    if 'is_upregulated' in Nodes.columns:
        Nodes = Nodes.drop(columns=['is_upregulated'])
    # Insert 'node_type' after 'symbol'
    Nodes.insert(Nodes.columns.get_loc('chr'), 'node_type', 'gene')

    # Step 15: Add 'start' and 'end' columns from 'TSS'
    print("Step 15: Adding 'start' and 'end' columns...")
    Nodes['start'] = Nodes['TSS']
    Nodes['end'] = Nodes['TSS']

    # Step 16: Remove 'TSS' and rename 'symbol' to 'identifier'
    print("Step 16: Finalizing columns...")
    if 'TSS' in Nodes.columns:
        Nodes = Nodes.drop(columns=['TSS'])
    Nodes = Nodes.rename(columns={'symbol': 'identifier'})

    # Step 17: Export the final Nodes table
    path_17 = output_dir / "Nodes.csv"
    print(f"Step 17: Exporting final 'Nodes' table to {path_17}...")
    Nodes.to_csv(path_17, index=False)

    # --- End of Processing ---
    end_time = time.time()
    print(f"\nProcessing complete.")
    print(f"Total records in final 'Nodes' table: {len(Nodes)}")
    print(f"Total execution time: {end_time - start_time:.2f} seconds.")


if __name__ == "__main__":
    # Step 18: Provide entire code (this is it)
    main()
