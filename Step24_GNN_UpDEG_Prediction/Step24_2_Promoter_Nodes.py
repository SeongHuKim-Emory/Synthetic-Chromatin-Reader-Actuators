import pandas as pd
import os
import sys

# --- Configuration ---

# Note on "25 CPUs": This task is I/O and memory-bound, not CPU-bound
# in a way that parallelizes well. Using pandas's optimized C-backend
# for mapping (as done below) is far more efficient than using
# multiprocessing for this specific lookup task.
# NUM_CPUS = 25 # Not used for reasons stated above.

# Note on "CUDA": CUDA is for GPU-based numerical computation (e.g.,
# deep learning) and is not applicable to this pandas-based ETL task.

# --- File Paths ---
UPDEG_FILE = '../../public/Multiomics/Step19_ML_UpDEG_Prediction/GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Pro.csv'
NONDEG_FILE = '../../public/Multiomics/Step19_ML_UpDEG_Prediction/GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Pro.csv'
NODES_FILE_IN = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/Nodes.csv'

OUTPUT_DIR = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/'
DSPP_NODES_FILE_OUT = os.path.join(OUTPUT_DIR, 'DSPP_Nodes.csv')
CP_NODES_FILE_OUT = os.path.join(OUTPUT_DIR, 'CP_Nodes.csv')
USPP_NODES_FILE_OUT = os.path.join(OUTPUT_DIR, 'USPP_Nodes.csv')
NODES_FILE_OUT = os.path.join(OUTPUT_DIR, 'Nodes.csv')

def create_new_nodes(base_nodes_df, current_main_nodes_df, info_map, node_type, suffix, start_col_name, end_col_name):
    """
    Generates a new nodes DataFrame based on the provided parameters.

    Args:
        base_nodes_df (pd.DataFrame): DataFrame with the original 'identifier' column.
        current_main_nodes_df (pd.DataFrame): The current state of the main nodes table
                                              (used to find max node_id).
        info_map (dict): A pre-computed dictionary mapping gene symbols to their info.
        node_type (str): The string to use for the 'node_type' column.
        suffix (str): The suffix to append to the 'identifier' (e.g., "_DSPP").
        start_col_name (str): The key in info_map to use for 'start' (e.g., "DSPP_start").
        end_col_name (str): The key in info_map to use for 'end' (e.g., "DSPP_end").

    Returns:
        pd.DataFrame: A newly created nodes DataFrame.
    """
    print(f"--- Creating {node_type} nodes ---")

    # Step 3, 11, 19: Create new table
    new_nodes = pd.DataFrame()

    # Step 4, 12, 20: Copy 'identifier' column from the *original* Nodes table
    new_nodes['identifier'] = base_nodes_df['identifier']

    # Step 5, 13, 21: Fill in 'node_id' starting from max(node_id) + 1
    max_id = current_main_nodes_df['node_id'].max()
    new_nodes['node_id'] = range(max_id + 1, max_id + 1 + len(new_nodes))
    print(f"New node IDs range from {max_id + 1} to {max_id + len(new_nodes)}")

    # Step 6, 14, 22: Fill in 'node_type'
    new_nodes['node_type'] = node_type

    # Step 7, 15, 23: Map 'chr', 'start', and 'end' values
    print("Mapping chr, start, and end columns...")
    # Create specific maps from the main info_map for efficient pandas mapping
    chr_map = {k: v.get('chr') for k, v in info_map.items()}
    start_map = {k: v.get(start_col_name) for k, v in info_map.items()}
    end_map = {k: v.get(end_col_name) for k, v in info_map.items()}

    new_nodes['chr'] = new_nodes['identifier'].map(chr_map)
    new_nodes['start'] = new_nodes['identifier'].map(start_map)
    new_nodes['end'] = new_nodes['identifier'].map(end_map)

    # Step 8, 16, 24: Add suffix to 'identifier'
    new_nodes['identifier'] = new_nodes['identifier'] + suffix
    
    # Handle potential NaNs if some identifiers weren't in the map
    original_len = len(new_nodes)
    new_nodes = new_nodes.dropna(subset=['chr', 'start', 'end'])
    if len(new_nodes) < original_len:
        print(f"Warning: Dropped {original_len - len(new_nodes)} rows due to missing info for {node_type}")

    # Convert start/end to integer types as requested (using nullable Int64)
    new_nodes['start'] = new_nodes['start'].astype(pd.Int64Dtype())
    new_nodes['end'] = new_nodes['end'].astype(pd.Int64Dtype())

    # Reorder columns to match request
    new_nodes = new_nodes[['node_id', 'identifier', 'node_type', 'chr', 'start', 'end']]
    
    print(f"Finished creating {node_type} nodes. Shape: {new_nodes.shape}")
    return new_nodes

def main():
    print(f"Script started. Python version: {sys.version}")
    
    # Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # --- Step 1: Import files ---
    print(f"Loading {UPDEG_FILE}...")
    updeg_df = pd.read_csv(UPDEG_FILE)
    print(f"Loading {NONDEG_FILE}...")
    nondeg_df = pd.read_csv(NONDEG_FILE)
    print(f"Loading {NODES_FILE_IN}...")
    nodes_df = pd.read_csv(NODES_FILE_IN)
    print("All source files loaded.")
    print(f"Initial Nodes.csv shape: {nodes_df.shape}")

    # --- Step 2: Make a new table named Nodes (already done by loading) ---
    # We will call our main DataFrame `main_nodes_df`
    main_nodes_df = nodes_df.copy()
    
    # Ensure start/end in the original nodes are integers (using nullable Int64)
    if 'start' in main_nodes_df.columns:
        main_nodes_df['start'] = main_nodes_df['start'].astype(pd.Int64Dtype())
    if 'end' in main_nodes_df.columns:
        main_nodes_df['end'] = main_nodes_df['end'].astype(pd.Int64Dtype())
    
    # Create a static copy of the *original* identifiers for creating the new node types
    base_identifiers_df = main_nodes_df[['identifier']].copy()

    # Combine UpDEG and Non-DEG files and create a fast lookup map
    print("Combining gene info files and creating lookup map...")
    gene_info_df = pd.concat([updeg_df, nondeg_df]).drop_duplicates(subset=['symbol'])
    # Convert to a dictionary for fast lookups
    # Format: {'A1BG': {'chr': 'chr19', 'DSPP_start': 123, 'CP_start': 456, ...}, ...}
    gene_info_map = gene_info_df.set_index('symbol').to_dict('index')
    print(f"Lookup map created with {len(gene_info_map)} unique symbols.")

    # --- Steps 3-10: Process DSPP_Nodes ---
    dspp_nodes_df = create_new_nodes(
        base_nodes_df=base_identifiers_df,
        current_main_nodes_df=main_nodes_df,
        info_map=gene_info_map,
        node_type="promoter_proximal_downstream",
        suffix="_DSPP",
        start_col_name="DSPP_start",
        end_col_name="DSPP_end"
    )
    
    # Step 9: Export DSPP_Nodes
    dspp_nodes_df.to_csv(DSPP_NODES_FILE_OUT, index=False)
    print(f"Exported DSPP_Nodes to {DSPP_NODES_FILE_OUT}")

    # Step 10: Merge DSPP_Nodes into main table
    main_nodes_df = pd.concat([main_nodes_df, dspp_nodes_df], ignore_index=True)
    print(f"Merged DSPP_Nodes. Main table shape: {main_nodes_df.shape}")

    # --- Steps 11-18: Process CP_Nodes ---
    cp_nodes_df = create_new_nodes(
        base_nodes_df=base_identifiers_df,
        current_main_nodes_df=main_nodes_df, # Pass the *updated* main table
        info_map=gene_info_map,
        node_type="promoter_core",
        suffix="_CP",
        start_col_name="CP_start",
        end_col_name="CP_end"
    )
    
    # Step 17: Export CP_Nodes
    cp_nodes_df.to_csv(CP_NODES_FILE_OUT, index=False)
    print(f"Exported CP_Nodes to {CP_NODES_FILE_OUT}")

    # Step 18: Merge CP_Nodes into main table
    main_nodes_df = pd.concat([main_nodes_df, cp_nodes_df], ignore_index=True)
    print(f"Merged CP_Nodes. Main table shape: {main_nodes_df.shape}")

    # --- Steps 19-26: Process USPP_Nodes ---
    uspp_nodes_df = create_new_nodes(
        base_nodes_df=base_identifiers_df,
        current_main_nodes_df=main_nodes_df, # Pass the *updated* main table
        info_map=gene_info_map,
        node_type="promoter_proximal_upstream",
        suffix="_USPP",
        start_col_name="USPP_start",
        end_col_name="USPP_end"
    )

    # Step 25: Export USPP_Nodes
    uspp_nodes_df.to_csv(USPP_NODES_FILE_OUT, index=False)
    print(f"Exported USPP_Nodes to {USPP_NODES_FILE_OUT}")

    # Step 26: Merge USPP_Nodes into main table
    main_nodes_df = pd.concat([main_nodes_df, uspp_nodes_df], ignore_index=True)
    print(f"Merged USPP_Nodes. Main table shape: {main_nodes_df.shape}")

    # --- Step 27: Export final Nodes table ---
    # Ensure final data types are reasonable (start/end should be integers)
    main_nodes_df['node_id'] = main_nodes_df['node_id'].astype(int)
    
    # Per user request, ensure start/end are (nullable) integers
    main_nodes_df['start'] = main_nodes_df['start'].astype(pd.Int64Dtype())
    main_nodes_df['end'] = main_nodes_df['end'].astype(pd.Int64Dtype())

    main_nodes_df.to_csv(NODES_FILE_OUT, index=False)
    print(f"Exported final updated Nodes table to {NODES_FILE_OUT}")
    print("--- Processing Complete ---")

if __name__ == "__main__":
    main()
