import pandas as pd
import requests
import time
import os

# --- Configuration ---
# Define the base directory relative to the script location
BASE_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'public', 'Multiomics', 'Step19_ML_UpDEG_Prediction')

# Define input and output file paths using the base directory
INPUT_FILES = {
    "updeg": os.path.join(BASE_DIR, "GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG.csv"),
    "nondeg": os.path.join(BASE_DIR, "GeneHancer_Genes_Elements_RNAseq_Non_DEG.csv")
}

OUTPUT_FILES = {
    "updeg": os.path.join(BASE_DIR, "GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS.csv"),
    "nondeg": os.path.join(BASE_DIR, "GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS.csv")
}


# --- Optimized Ensembl API Function ---

def get_gene_info_batch(identifiers):
    """
    Fetches gene details for a list of identifiers using efficient batch requests to the Ensembl API.
    Handles both Ensembl IDs and Gene Symbols and processes them in chunks of 1000.
    """
    gene_info_map = {}

    # Separate identifiers into symbols and Ensembl IDs to use the correct API endpoints
    symbols = [i for i in identifiers if isinstance(i, str) and not i.startswith("ENSG")]
    ensembl_ids = [i for i in identifiers if isinstance(i, str) and i.startswith("ENSG")]

    server = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    # Process symbols in batches of 1000 (the API limit)
    for i in range(0, len(symbols), 1000):
        batch = symbols[i:i + 1000]
        ext_symbol = "/lookup/symbol/homo_sapiens"
        try:
            r_symbol = requests.post(server + ext_symbol, headers=headers, json={"symbols": batch}, timeout=60)
            r_symbol.raise_for_status()
            decoded_symbols = r_symbol.json()
            for symbol, data in decoded_symbols.items():
                if data:
                    strand_val = data.get('strand', 0)
                    gene_info_map[symbol] = {
                        "gene_chr": f"chr{data.get('seq_region_name', '')}",
                        "gene_start": data.get('start', ''),
                        "gene_end": data.get('end', ''),
                        "gene_strand": "reverse" if strand_val == -1 else "forward" if strand_val == 1 else ""
                    }
        except requests.exceptions.RequestException as e:
            print(f"Warning: Batch symbol lookup failed for a chunk. Error: {e}")

    # Process Ensembl IDs in batches of 1000
    for i in range(0, len(ensembl_ids), 1000):
        batch = ensembl_ids[i:i + 1000]
        ext_id = "/lookup/id"
        try:
            r_id = requests.post(server + ext_id, headers=headers, json={"ids": batch}, timeout=60)
            r_id.raise_for_status()
            decoded_ids = r_id.json()
            for ensg_id, data in decoded_ids.items():
                if data:
                    strand_val = data.get('strand', 0)
                    gene_info_map[ensg_id] = {
                        "gene_chr": f"chr{data.get('seq_region_name', '')}",
                        "gene_start": data.get('start', ''),
                        "gene_end": data.get('end', ''),
                        "gene_strand": "reverse" if strand_val == -1 else "forward" if strand_val == 1 else ""
                    }
        except requests.exceptions.RequestException as e:
            print(f"Warning: Batch ID lookup failed for a chunk. Error: {e}")

    return gene_info_map


def process_file(input_path, output_path):
    """
    Main processing function to apply all transformation steps to a given file.
    """
    print(f"--- Starting processing for: {os.path.basename(input_path)} ---")

    # Step 1 & 3: Import file and create a new table (DataFrame)
    try:
        df = pd.read_csv(input_path)
        print(f"Step 1: Successfully loaded {len(df)} rows.")
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}")
        return

    # Step 6: Remove specified columns
    columns_to_drop = [
        'GHid', 'combined_score', 'is_gene_elite', 'is_elite',
        'regulatory_element_type', 'enhancer_score'
    ]
    df.drop(columns=columns_to_drop, inplace=True, errors='ignore')
    print("Step 6: Removed specified columns.")

    # Step 7: Add "chr" prefix to the 'chr' column
    df['chr'] = 'chr' + df['chr'].astype(str)
    print("Step 7: Added 'chr' prefix to chromosome numbers.")

    # Step 8: Add new columns for gene information
    new_cols = ['gene_chr', 'gene_start', 'gene_end', 'gene_strand', 'TSS']
    for col in new_cols:
        df[col] = pd.NA
    print("Step 8: Added new columns for gene details and TSS.")

    # Step 9: Fetch gene info in a batch and populate columns
    print("Step 9: Fetching gene information from Ensembl using batch API. This should be much faster...")
    unique_symbols = df['symbol'].dropna().unique()

    # Get the mapping from the batch function
    gene_data_map = get_gene_info_batch(list(unique_symbols))

    # Use the efficient pandas `map` and `apply` methods to populate the new columns
    # This is significantly faster than iterating row by row.
    gene_info_series = df['symbol'].map(gene_data_map)

    # Create a temporary DataFrame from the mapped dictionaries
    gene_info_df = pd.DataFrame(gene_info_series.dropna().tolist(), index=gene_info_series.dropna().index)

    # Update the main DataFrame with the fetched data
    if not gene_info_df.empty:
        df.update(gene_info_df)

    print("  > Gene information populated.")

    # Step 10: Calculate TSS based on gene strand
    # Ensure coordinate columns are numeric before calculation
    df['gene_start'] = pd.to_numeric(df['gene_start'], errors='coerce')
    df['gene_end'] = pd.to_numeric(df['gene_end'], errors='coerce')

    # For 'forward' strand, TSS is the gene_start
    df.loc[df['gene_strand'] == 'forward', 'TSS'] = df['gene_start']
    # For 'reverse' strand, TSS is the gene_end
    df.loc[df['gene_strand'] == 'reverse', 'TSS'] = df['gene_end']
    print("Step 10: Calculated TSS based on gene strand.")

    # Step 11/12: Export the processed DataFrame to a new CSV file
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"Step 11/12: Successfully exported processed data to {output_path}\n")


def main():
    """
    Main function to run the processing for both DEG and Non-DEG files.
    """
    print("Starting GeneHancer Data Processing Pipeline...")
    process_file(INPUT_FILES["updeg"], OUTPUT_FILES["updeg"])
    process_file(INPUT_FILES["nondeg"], OUTPUT_FILES["nondeg"])
    print("Pipeline finished successfully.")


if __name__ == "__main__":
    main()

