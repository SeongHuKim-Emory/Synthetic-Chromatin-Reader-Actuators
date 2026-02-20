# This script performs data cleaning, transformation, and merging of multiomics data.
# It processes RNA-seq data to extract TPM values, maps gene IDs to gene symbols,
# and then merges this expression data with two other histone modification datasets.

# Required libraries: pandas for data manipulation and mygene for gene ID conversion.
# You can install them using pip:
# pip install pandas mygene

import pandas as pd
import mygene
import os


def process_data():
    """
    Main function to execute the entire data processing pipeline.
    """
    # --- File Paths ---
    # Define input and output directories for better organization.
    input_dir = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/'
    rnaseq_input_dir = '../../public/Public Dataset/RNAseq/'
    output_dir = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/'

    # Ensure the output directory exists.
    os.makedirs(output_dir, exist_ok=True)

    # Input file paths
    updeg_histone_file = os.path.join(input_dir,
                                      'GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone.csv')
    nondeg_histone_file = os.path.join(input_dir,
                                       'GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone.csv')
    rnaseq_file = os.path.join(rnaseq_input_dir, 'GSE175204 ENCFF721BRA hg38 MCF7 total RNAseq ENSG Only.tsv')

    # Output file paths
    processed_rnaseq_output_file = os.path.join(output_dir,
                                                'GSE175204 ENCFF721BRA hg38 MCF7 total RNAseq ENSG Only Symbol and TPM.csv')
    updeg_output_file = os.path.join(output_dir,
                                     'GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone_RNA.csv')
    nondeg_output_file = os.path.join(output_dir,
                                      'GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone_RNA.csv')

    print("--- Step 1: Loading initial datasets ---")
    try:
        updeg_histone_df = pd.read_csv(updeg_histone_file)
        nondeg_histone_df = pd.read_csv(nondeg_histone_file)
        rnaseq_df = pd.read_csv(rnaseq_file, sep='\t')
        print("Initial files loaded successfully.")
    except FileNotFoundError as e:
        print(f"Error: {e}. Please make sure the input files are in the correct directories.")
        return

    # --- Steps 2-6: Process RNA-seq data ---
    print("\n--- Steps 2-6: Processing RNA-seq data ---")

    # 2. Keep only "gene_id" and "TPM" columns
    rnaseq_processed = rnaseq_df[['gene_id', 'TPM']].copy()
    print("Filtered RNA-seq data to 'gene_id' and 'TPM' columns.")

    # 3. Add a new "gene_symbol" column
    rnaseq_processed.insert(1, 'gene_symbol', '')
    print("Added 'gene_symbol' column.")

    # 5. Fill in the gene_symbol values using mygene.info service
    print("Fetching gene symbols for Ensembl IDs. This may take a moment...")
    mg = mygene.MyGeneInfo()
    # Remove version numbers from Ensembl IDs for querying
    gene_ids_no_version = rnaseq_processed['gene_id'].str.split('.').str[0]
    gene_info = mg.querymany(gene_ids_no_version, scopes='ensembl.gene', fields='symbol', species='human',
                             as_dataframe=True)

    # Create a mapping from the original gene_id to the new symbol
    # Reset index to make 'query' (the Ensembl ID) a column
    gene_info.reset_index(inplace=True)

    # Handle cases where mygene returns duplicate results for a single query
    # by keeping only the first result. This ensures the index for mapping is unique.
    gene_info.drop_duplicates(subset='query', keep='first', inplace=True)

    symbol_map = gene_info.set_index('query')['symbol']

    # Map symbols back to the dataframe
    rnaseq_processed['gene_symbol'] = gene_ids_no_version.map(symbol_map)
    # Fill any missing symbols with NaN or a placeholder if preferred
    rnaseq_processed['gene_symbol'].fillna('N/A', inplace=True)
    print("Gene symbols have been populated.")

    # 6. Name the table
    print("Renaming processed RNA-seq table in memory.")
    gse_processed_df = rnaseq_processed

    # --- Step 7: Export the processed RNA-seq table ---
    print("\n--- Step 7: Exporting processed RNA-seq data ---")
    gse_processed_df.to_csv(processed_rnaseq_output_file, index=False)
    print(f"Successfully exported processed RNA-seq data to:\n{processed_rnaseq_output_file}")

    # --- Steps 8-11: Create copies of histone dataframes ---
    print("\n--- Steps 8-11: Creating new tables for merging ---")
    # 8-9. Create and populate the UpDEG RNA table
    updeg_rna_df = updeg_histone_df.copy()
    print("Created 'GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone_RNA'.")

    # 10-11. Create and populate the Non-DEG RNA table
    nondeg_rna_df = nondeg_histone_df.copy()
    print("Created 'GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone_RNA'.")

    # --- Step 12: Add 'TPM' column to new tables ---
    print("\n--- Step 12: Adding 'TPM' column to histone tables ---")
    updeg_rna_df['TPM'] = 0.0
    nondeg_rna_df['TPM'] = 0.0
    print("'TPM' column added to both tables.")

    # --- Step 13: Map TPM values based on gene symbol ---
    print("\n--- Step 13: Mapping TPM values to histone tables ---")
    # Create a mapping dictionary from the processed RNA-seq data
    # Drop duplicates to ensure a unique mapping for each gene symbol
    tpm_map = gse_processed_df.drop_duplicates(subset=['gene_symbol']).set_index('gene_symbol')['TPM']

    # Map the TPM values to the new dataframes
    updeg_rna_df['TPM'] = updeg_rna_df['symbol'].map(tpm_map)
    nondeg_rna_df['TPM'] = nondeg_rna_df['symbol'].map(tpm_map)

    # Fill any symbols that didn't have a match with 0
    updeg_rna_df['TPM'].fillna(0, inplace=True)
    nondeg_rna_df['TPM'].fillna(0, inplace=True)
    print("TPM values have been successfully mapped.")

    # --- Steps 14-15: Export the final merged tables ---
    print("\n--- Steps 14-15: Exporting final merged tables ---")

    # 14. Export the UpDEG merged table
    updeg_rna_df.to_csv(updeg_output_file, index=False)
    print(f"Successfully exported UpDEG data with RNA TPM to:\n{updeg_output_file}")

    # 15. Export the Non-DEG merged table
    nondeg_rna_df.to_csv(nondeg_output_file, index=False)
    print(f"Successfully exported Non-DEG data with RNA TPM to:\n{nondeg_output_file}")

    print("\n--- Script finished successfully! ---")


# Execute the main function when the script is run
if __name__ == "__main__":
    process_data()

