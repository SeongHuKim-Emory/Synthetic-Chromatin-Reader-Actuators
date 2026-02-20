import pandas as pd
import os
import sys

def process_genomic_nodes():
    """
    Reads node and HiC loop files, generates new enhancer nodes from
    the loops, and merges them into the main node list.
    """
    
    # --- File Paths ---
    nodes_in_path = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/Nodes.csv'
    hic_in_path = '../../public/Public Dataset/HiC/MCF7 HiC Loops GSE237722 ENCFF797NNQ hg38.bedpe'
    
    enhancer_out_path = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/Enhancer_Nodes.csv'
    nodes_out_path = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/Nodes.csv' # Output overwrites original
    
    # Ensure output directories exist
    os.makedirs(os.path.dirname(enhancer_out_path), exist_ok=True)

    try:
        # --- Step 1: Import files ---
        print(f"Step 1: Importing files...")
        # 1. Import Nodes.csv
        nodes_df = pd.read_csv(nodes_in_path)
        print(f"  Loaded {len(nodes_df)} nodes from {nodes_in_path}")

        # 1. Import MCF7 HiC Loops...
        # We set comment=None because the header row itself starts with '#'
        hic_df = pd.read_csv(hic_in_path, sep='\t', comment=None)
        print(f"  Loaded {len(hic_df)} HiC loops from {hic_in_path}")

        # --- Step 2: Make a new table named Nodes ---
        # This is implicitly done by loading 'nodes_df'. We will operate on this DataFrame.
        print(f"Step 2: 'Nodes' table created in memory.")

        # --- Step 3: Make a new table named Enhancer_Nodes ---
        # We will build this DataFrame populated with data.
        print(f"Step 3: Preparing 'Enhancer_Nodes' table...")

        # --- Steps 4, 5, 6, 7 (Combined for efficiency) ---
        print(f"Steps 4-7: Populating 'Enhancer_Nodes'...")
        
        num_loops = len(hic_df)
        if num_loops == 0:
            print("  Warning: HiC file is empty. No enhancer nodes will be generated.")
            enhancer_nodes_df = pd.DataFrame(columns=['node_id', 'identifier', 'node_type', 'chr', 'start', 'end'])
        
        else:
            # Create loop numbers (1-based index)
            loop_numbers = range(1, num_loops + 1)
            
            # Create 'x' nodes (from #chr1, x1, x2)
            x_nodes = pd.DataFrame({
                'identifier': [f'enhancer_loop_{i}_x' for i in loop_numbers],
                'chr': hic_df['#chr1'].values,
                'start': hic_df['x1'].values,
                'end': hic_df['x2'].values,
                'loop_num': loop_numbers # Helper for sorting
            })
            
            # Create 'y' nodes (from chr2, y1, y2)
            y_nodes = pd.DataFrame({
                'identifier': [f'enhancer_loop_{i}_y' for i in loop_numbers],
                'chr': hic_df['chr2'].values,
                'start': hic_df['y1'].values,
                'end': hic_df['y2'].values,
                'loop_num': loop_numbers # Helper for sorting
            })
            
            # Combine x and y nodes
            combined_nodes = pd.concat([x_nodes, y_nodes])
            
            # Sort by loop_num (1, 1, 2, 2, ...) and then identifier (x, y)
            combined_nodes = combined_nodes.sort_values(
                by=['loop_num', 'identifier']
            ).reset_index(drop=True)
            
            # Create the final Enhancer_Nodes DataFrame
            enhancer_nodes_df = combined_nodes.drop(columns=['loop_num'])

            # --- Step 4: Fill in node_id ---
            max_node_id = nodes_df['node_id'].max()
            num_new_nodes = len(enhancer_nodes_df)
            new_node_ids = range(max_node_id + 1, max_node_id + 1 + num_new_nodes)
            enhancer_nodes_df['node_id'] = new_node_ids
            
            # --- Step 7: Fill in node_type ---
            enhancer_nodes_df['node_type'] = 'enhancer'
            
            # --- Step 6 (cont.): Ensure integer types ---
            enhancer_nodes_df['start'] = enhancer_nodes_df['start'].astype(int)
            enhancer_nodes_df['end'] = enhancer_nodes_df['end'].astype(int)
            
            # Reorder columns to match specification
            enhancer_nodes_df = enhancer_nodes_df[
                ['node_id', 'identifier', 'node_type', 'chr', 'start', 'end']
            ]

        print(f"  Generated {len(enhancer_nodes_df)} new enhancer nodes.")
        
        # --- Step 8: Export the Enhancer_Nodes table ---
        enhancer_nodes_df.to_csv(enhancer_out_path, index=False)
        print(f"Step 8: Exported 'Enhancer_Nodes' to {enhancer_out_path}")

        # --- Step 9: Merge the rows of Enhancer_Nodes into Nodes ---
        merged_nodes_df = pd.concat([nodes_df, enhancer_nodes_df], ignore_index=True)
        print(f"Step 9: Merged tables. Total nodes: {len(merged_nodes_df)}")

        # --- Step 10: Export the Nodes table ---
        merged_nodes_df.to_csv(nodes_out_path, index=False)
        print(f"Step 10: Exported final 'Nodes' table to {nodes_out_path}")
        
        print("\nProcess completed successfully.")

    except FileNotFoundError as e:
        print(f"Error: Input file not found.")
        print(e, file=sys.stderr)
    except KeyError as e:
        print(f"Error: A required column was not found. Check file formats.")
        print(f"Missing column: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred:")
        print(e, file=sys.stderr)

if __name__ == "__main__":
    # --- Step 11: Provide entire code ---
    process_genomic_nodes()
