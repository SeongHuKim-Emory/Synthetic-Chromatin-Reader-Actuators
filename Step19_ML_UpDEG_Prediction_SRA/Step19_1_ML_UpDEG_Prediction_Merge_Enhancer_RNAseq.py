import pandas as pd
import os

# Define file paths
gexp_counts_path = '../../public/RNAseq/gexp_counts_symbol.csv'
genehancer_path = '../../public/Public Dataset/GeneHancer/GeneHancer_Genes_Elements.csv'
output_dir = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/'

# --- Step 1: Create output directory if it doesn't exist ---
# This ensures that the script won't fail if the destination folder is missing.
os.makedirs(output_dir, exist_ok=True)
print(f"Output directory '{output_dir}' is ready.")

# --- Step 2: Import gexp_counts_symbol.csv and process transcript_id ---
print(f"Reading gene expression data from '{gexp_counts_path}'...")
try:
    gexp_counts_df = pd.read_csv(gexp_counts_path)

    # Create a copy to avoid modifying the original DataFrame in memory
    gexp_counts_symbol_no_id_version = gexp_counts_df.copy()

    # Remove decimal and subsequent values from the 'transcript_id' column
    # For example, 'ENSG00000000003.14' becomes 'ENSG00000000003'
    gexp_counts_symbol_no_id_version['transcript_id'] = \
    gexp_counts_symbol_no_id_version['transcript_id'].str.split('.').str[0]
    print("Successfully processed 'transcript_id' column.")

    # --- Step 3 & 4: Export the modified gene expression table ---
    output_gexp_path = os.path.join(output_dir, 'gexp_counts_symbol_no_id_version.csv')
    gexp_counts_symbol_no_id_version.to_csv(output_gexp_path, index=False)
    print(f"Exported processed gene expression data to '{output_gexp_path}'.")

    # --- Step 5: Filter GeneHancer data based on the processed gene expression data ---
    print(f"Reading GeneHancer data from '{genehancer_path}'...")
    genehancer_df = pd.read_csv(genehancer_path)

    # Create a set of all valid identifiers from both 'transcript_id' and 'symbol' columns
    # Using a set provides a fast lookup for the filtering step.
    valid_ids = set(gexp_counts_symbol_no_id_version['transcript_id']).union(
        set(gexp_counts_symbol_no_id_version['symbol']))

    # Keep rows in GeneHancer data where the 'symbol' exists in our set of valid_ids
    initial_rows = len(genehancer_df)
    GeneHancer_Genes_Elements_RNAseq = genehancer_df[genehancer_df['symbol'].isin(valid_ids)].copy()
    final_rows = len(GeneHancer_Genes_Elements_RNAseq)

    print(f"Filtered GeneHancer data. Kept {final_rows} rows out of {initial_rows}.")

    # --- Step 6 & 7: Export the filtered GeneHancer table ---
    output_genehancer_path = os.path.join(output_dir, 'GeneHancer_Genes_Elements_RNAseq.csv')
    GeneHancer_Genes_Elements_RNAseq.to_csv(output_genehancer_path, index=False)
    print(f"Exported filtered GeneHancer data to '{output_genehancer_path}'.")

    print("\nScript finished successfully.")

except FileNotFoundError as e:
    print(f"Error: {e}. Please ensure the input file paths are correct.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
