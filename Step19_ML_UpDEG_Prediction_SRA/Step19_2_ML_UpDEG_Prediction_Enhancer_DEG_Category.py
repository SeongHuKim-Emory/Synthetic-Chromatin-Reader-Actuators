import pandas as pd
import os
import numpy as np

# --- 1. Import Files ---

# Define file paths for clarity
genehancer_path = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq.csv'
updeg_path = '../../public/RNAseq/UpDEG_category.csv'
downdeg_path = '../../public/RNAseq/DownDEG_category.csv'

# Define output directory
output_dir = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/'

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

try:
    # Load the datasets into pandas DataFrames
    df_genehancer = pd.read_csv(genehancer_path)
    df_updeg = pd.read_csv(updeg_path)
    df_downdeg = pd.read_csv(downdeg_path)
    print("Successfully imported all CSV files.")
except FileNotFoundError as e:
    print(f"Error: {e}. Please ensure the input file paths are correct.")
    exit()

# --- 2. Clean 'transcript_id' Column ---
# In both DEG tables, remove the version number (decimal and following digits)
# from the 'transcript_id' column.

print("Cleaning 'transcript_id' columns...")
df_updeg['transcript_id'] = df_updeg['transcript_id'].str.split('.').str[0]
df_downdeg['transcript_id'] = df_downdeg['transcript_id'].str.split('.').str[0]


# --- 3. Create Updated DataFrame Copies ---
# These copies will be modified further.
UpDEG_category_update = df_updeg.copy()
DownDEG_category_update = df_downdeg.copy()
print("Created up-to-date copies of DEG tables.")

# --- 4. Fill Empty 'Symbol' Columns ---
# If a row's 'Symbol' is empty (NaN), replace it with the 'transcript_id' value.

print("Filling empty 'Symbol' values...")
# Replace empty strings with NaN to ensure fillna works correctly
UpDEG_category_update['Symbol'].replace('', np.nan, inplace=True)
DownDEG_category_update['Symbol'].replace('', np.nan, inplace=True)

UpDEG_category_update['Symbol'].fillna(UpDEG_category_update['transcript_id'], inplace=True)
DownDEG_category_update['Symbol'].fillna(DownDEG_category_update['transcript_id'], inplace=True)


# --- 5. Export Updated Category Tables ---
updeg_output_path = os.path.join(output_dir, 'UpDEG_category_update.csv')
downdeg_output_path = os.path.join(output_dir, 'DownDEG_category_update.csv')

print(f"Exporting updated category tables to {output_dir}...")
UpDEG_category_update.to_csv(updeg_output_path, index=False)
DownDEG_category_update.to_csv(downdeg_output_path, index=False)


# --- 6-8. Create and Populate SRA-Specific UpDEG Table ---

print("Filtering for SRA-specific Upregulated DEG enhancers...")
# Get a set of symbols from the UpDEG table where 'Category - SRA specific' is not empty.
# Using .dropna() handles both NaN and None values.
sra_specific_symbols = set(UpDEG_category_update.dropna(subset=['Category - SRA specific'])['Symbol'])

# Filter the original GeneHancer table to keep rows where the symbol is in our set.
GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG = df_genehancer[df_genehancer['symbol'].isin(sra_specific_symbols)].copy()


# --- 9. Export SRA-Specific UpDEG Table ---
sra_updeg_output_path = os.path.join(output_dir, 'GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG.csv')

print(f"Exporting SRA-specific UpDEG table to {sra_updeg_output_path}...")
GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG.to_csv(sra_updeg_output_path, index=False)


# --- 10 & 11. Create Non-DEG Table by Copying Original Data ---

print("Creating the Non-DEG table...")
GeneHancer_Genes_Elements_RNAseq_Non_DEG = df_genehancer.copy()


# --- 12 & 13. Filter Non-DEG Table ---
# Get all unique symbols from both updated UpDEG and DownDEG tables.
updeg_symbols = set(UpDEG_category_update['Symbol'])
downdeg_symbols = set(DownDEG_category_update['Symbol'])

# Combine them into a single set of all DEGs
all_deg_symbols = updeg_symbols.union(downdeg_symbols)

# Filter the Non-DEG table to *exclude* any rows where the symbol is in the combined set.
# The '~' operator inverts the boolean mask created by .isin().
print("Removing all DEG-associated enhancers from the Non-DEG table...")
GeneHancer_Genes_Elements_RNAseq_Non_DEG = GeneHancer_Genes_Elements_RNAseq_Non_DEG[
    ~GeneHancer_Genes_Elements_RNAseq_Non_DEG['symbol'].isin(all_deg_symbols)
]


# --- 14. Export Non-DEG Table ---
non_deg_output_path = os.path.join(output_dir, 'GeneHancer_Genes_Elements_RNAseq_Non_DEG.csv')

print(f"Exporting Non-DEG table to {non_deg_output_path}...")
GeneHancer_Genes_Elements_RNAseq_Non_DEG.to_csv(non_deg_output_path, index=False)

print("\nAll tasks completed successfully!")
