import polars as pl
import os
import time

# --- Configuration ---

# Input file path
NODES_INPUT_PATH = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/Nodes.csv'
# Output file path
EDGES_OUTPUT_PATH = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/Edges.csv'

# --- Helper Functions ---

def create_gene_promoter_edges_pl(nodes_lf, gene_nodes_lf, node_type, suffix, edge_flags):
    """
    Generic function to create edges between a promoter type and genes using Polars Lazy API.
    
    Args:
        nodes_lf (pl.LazyFrame): The complete nodes LazyFrame.
        gene_nodes_lf (pl.LazyFrame): A pre-filtered LazyFrame of gene nodes.
        node_type (str): The node_type to filter for (e.g., 'promoter_core').
        suffix (str): The suffix to remove from the identifier (e.g., '_CP').
        edge_flags (dict): A dictionary of the flag columns for this edge type.

    Returns:
        pl.LazyFrame: A LazyFrame representing the new edges.
    """
    print(f"  Defining lazy query for '{node_type}' nodes...")
    
    # Filter for the specific promoter type
    promoter_nodes_lf = nodes_lf.filter(pl.col('node_type') == node_type)
    
    # Derive the base gene identifier by removing the suffix
    # Create edges by joining with the gene nodes
    edges_lf = promoter_nodes_lf.join(
        gene_nodes_lf, 
        left_on=pl.col('identifier').str.replace(f"{suffix}$", ""), 
        right_on='identifier'
    ).select(
        pl.col('node_id').alias('node_id_source'),
        pl.col('node_id_target'),
        *[pl.lit(v, dtype=pl.UInt8).alias(k) for k, v in edge_flags.items()] # Add flags
    )
    
    return edges_lf

def create_enhancer_loop_edges_pl(nodes_lf, edge_flags):
    """
    Creates edges between 'enhancer_loop_..._x' and 'enhancer_loop_..._y' nodes using Polars Lazy API.
    
    Args:
        nodes_lf (pl.LazyFrame): The complete nodes LazyFrame.
        edge_flags (dict): A dictionary of the flag columns for this edge type.

    Returns:
        pl.LazyFrame: A LazyFrame representing the new edges.
    """
    print("  Defining lazy query for enhancer loop edges...")
    
    # Filter for all enhancer nodes
    enhancer_nodes_lf = nodes_lf.filter(pl.col('node_type') == 'enhancer').select(['node_id', 'identifier'])
    
    # Separate '_x' nodes (sources)
    enhancer_x_lf = enhancer_nodes_lf.filter(
        pl.col('identifier').str.ends_with('_x')
    ).rename({'node_id': 'node_id_source'})
    
    # Separate '_y' nodes (targets)
    enhancer_y_lf = enhancer_nodes_lf.filter(
        pl.col('identifier').str.ends_with('_y')
    ).rename({'node_id': 'node_id_target'})
    
    # Join x and y nodes on the base identifier
    edges_lf = enhancer_x_lf.join(
        enhancer_y_lf,
        left_on=pl.col('identifier').str.replace("_x$", ""),
        right_on=pl.col('identifier').str.replace("_y$", "")
    ).select(
        pl.col('node_id_source'),
        pl.col('node_id_target'),
        *[pl.lit(v, dtype=pl.UInt8).alias(k) for k, v in edge_flags.items()] # Add flags
    )
    
    return edges_lf

# --- Main Execution ---

def main():
    """
    Main function to load nodes, generate edges, and save the results using Polars.
    """
    
    print("--- Edge Generation Script Started (Polars Version) ---")
    start_time = time.time()

    # 1. & 3. Import Nodes.csv (Lazy)
    print(f"Scanning nodes from {NODES_INPUT_PATH}...")
    try:
        # Use scan_csv for lazy loading
        nodes_lf = pl.scan_csv(NODES_INPUT_PATH)
    except Exception as e: # Catch potential file not found or read errors
        print(f"Error: Input file not found or could not be read: {NODES_INPUT_PATH}")
        print(f"Details: {e}")
        return
        
    print(f"Nodes scanner created in {time.time() - start_time:.2f} seconds.")

    # 2. Initialize Edges list (will store lazy frames)
    all_edges_lf = []
    
    # --- Pre-processing Step ---
    # Create a lazy frame for gene nodes to be joined
    print("Defining gene lookup frame...")
    prep_start = time.time()
    gene_nodes_lf = nodes_lf.filter(
        pl.col('node_type') == 'gene'
    ).select(
        'identifier', 
        pl.col('node_id').alias('node_id_target')
    ).cache() # Cache this intermediate result as it's used multiple times
    print(f"Gene frame defined in {time.time() - prep_start:.2f} seconds.")

    # 4. Process promoter_proximal_downstream (DSPP)
    print("\nProcessing Gene-DSPP edges...")
    dspp_flags = {'is_gene_DSPP': 1, 'is_gene_CP': 0, 'is_gene_USPP': 0, 'is_loop': 0, 'is_TAD': 0, 'is_overlap': 0}
    dspp_edges_lf = create_gene_promoter_edges_pl(nodes_lf, gene_nodes_lf, 'promoter_proximal_downstream', '_DSPP', dspp_flags)
    all_edges_lf.append(dspp_edges_lf)

    # 5. Process promoter_core (CP)
    print("\nProcessing Gene-CP edges...")
    cp_flags = {'is_gene_DSPP': 0, 'is_gene_CP': 1, 'is_gene_USPP': 0, 'is_loop': 0, 'is_TAD': 0, 'is_overlap': 0}
    cp_edges_lf = create_gene_promoter_edges_pl(nodes_lf, gene_nodes_lf, 'promoter_core', '_CP', cp_flags)
    all_edges_lf.append(cp_edges_lf)
    
    # 6. Process promoter_proximal_upstream (USPP)
    print("\nProcessing Gene-USPP edges...")
    uspp_flags = {'is_gene_DSPP': 0, 'is_gene_CP': 0, 'is_gene_USPP': 1, 'is_loop': 0, 'is_TAD': 0, 'is_overlap': 0}
    uspp_edges_lf = create_gene_promoter_edges_pl(nodes_lf, gene_nodes_lf, 'promoter_proximal_upstream', '_USPP', uspp_flags)
    all_edges_lf.append(uspp_edges_lf)

    # 7. Process enhancer loops
    print("\nProcessing enhancer loop edges...")
    loop_flags = {'is_gene_DSPP': 0, 'is_gene_CP': 0, 'is_gene_USPP': 0, 'is_loop': 1, 'is_TAD': 0, 'is_overlap': 0}
    loop_edges_lf = create_enhancer_loop_edges_pl(nodes_lf, loop_flags)
    all_edges_lf.append(loop_edges_lf)
    
    # --- Finalization Step ---
    print(f"\nTotal edge queries defined: {len(all_edges_lf)}")
    print("Concatenating and executing lazy queries...")
    exec_start = time.time()
    
    if not all_edges_lf:
        print("\nWarning: No edge queries were defined.")
        edges_df = pl.DataFrame(columns=[
            'edge_id', 'node_id_source', 'node_id_target', 'is_gene_DSPP', 
            'is_gene_CP', 'is_gene_USPP', 'is_loop', 'is_TAD', 'is_overlap'
        ],
        schema={
            'edge_id': pl.UInt64, 'node_id_source': pl.UInt64, 'node_id_target': pl.UInt64,
            'is_gene_DSPP': pl.UInt8, 'is_gene_CP': pl.UInt8, 'is_gene_USPP': pl.UInt8,
            'is_loop': pl.UInt8, 'is_TAD': pl.UInt8, 'is_overlap': pl.UInt8
        })
    else:
        # Concatenate all lazy frames vertically
        final_lf = pl.concat(all_edges_lf, how='vertical')
        
        # Add the 'edge_id' column, starting from 0
        final_lf = final_lf.with_row_index('edge_id')
        
        # Reorder columns to match the requested specification
        final_columns = [
            'edge_id', 'node_id_source', 'node_id_target', 'is_gene_DSPP', 
            'is_gene_CP', 'is_gene_USPP', 'is_loop', 'is_TAD', 'is_overlap'
        ]
        
        # Eagerly collect the final result
        edges_df = final_lf.select(final_columns).collect()
        
        print(f"DataFrame computed in {time.time() - exec_start:.2f} seconds.")
        print(f"Total edges found: {len(edges_df)}")

    # 8. Export the Edges table
    print(f"Saving DataFrame to {EDGES_OUTPUT_PATH}...")
    save_start = time.time()
    try:
        # Ensure the output directory exists
        output_dir = os.path.dirname(EDGES_OUTPUT_PATH)
        if output_dir: # Check if not an empty string (e.g., saving to root)
            os.makedirs(output_dir, exist_ok=True)
            
        # Save to CSV
        edges_df.write_csv(EDGES_OUTPUT_PATH)
        print(f"File saved successfully in {time.time() - save_start:.2f} seconds.")
    except Exception as e:
        print(f"An error occurred while saving the file: {e}")

    print(f"\n--- Script finished in {time.time() - start_time:.2f} seconds ---")

if __name__ == "__main__":
    # This block ensures the main() function runs when the script is executed
    main()
