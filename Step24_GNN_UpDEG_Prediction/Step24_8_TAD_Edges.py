import polars as pl
from itertools import product
import sys
import os

def main():
    """
    Main function to load genomic nodes and edges, identify new edges
    based on TAD (Topologically Associating Domain) co-localization
    of genes and enhancers, and save the updated edge list.
    Uses polars library.
    """
    
    # --- 1. Define File Paths ---
    # Define the relative paths to the input and output files
    nodes_path = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/Nodes.csv'
    edges_path = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/Edges.csv'
    tad_path = '../../public/Public Dataset/HiC/MCF7 HiC Contact Domains GSE237722 ENCFF164AGX hg38.bedpe'

    print("Starting script...")

    try:
        # --- 2. & 3. Import Files ---
        # Create new tables (DataFrames) from the CSV files
        print(f"Loading nodes from {nodes_path}...")
        nodes_df = pl.read_csv(nodes_path)
        
        print(f"Loading edges from {edges_path}...")
        try:
            edges_df = pl.read_csv(edges_path)
        except pl.exceptions.NoDataError:
            print("Edges file is empty, starting with an empty DataFrame.")
            # Define schema based on description to allow concatenation later
            edges_df = pl.DataFrame({
                'edge_id': pl.Series([], dtype=pl.Int64),
                'node_id_source': pl.Series([], dtype=pl.Int64),
                'node_id_target': pl.Series([], dtype=pl.Int64),
                'is_gene_DSPP': pl.Series([], dtype=pl.Int64),
                'is_gene_CP': pl.Series([], dtype=pl.Int64),
                'is_gene_USPP': pl.Series([], dtype=pl.Int64),
                'is_loop': pl.Series([], dtype=pl.Int64),
                'is_TAD': pl.Series([], dtype=pl.Int64),
                'is_overlap': pl.Series([], dtype=pl.Int64)
            })

        print(f"Loading TAD definitions from {tad_path}...")
        # The .bedpe file is typically tab-separated.
        # polars read_csv handles the header starting with '#' fine.
        tad_df = pl.read_csv(tad_path, separator='\t')

        print("Files loaded successfully.")
        print(f"Original node count: {len(nodes_df)}")
        print(f"Original edge count: {len(edges_df)}")
        print(f"TAD count: {len(tad_df)}")

        # --- 4. Process TADs and Create New Edges ---
        
        # This list will store dictionaries, each representing a new edge row
        new_edges_list = []

        # Find the maximum existing edge_id to start new IDs from
        if not edges_df.is_empty():
            max_edge_id = edges_df['edge_id'].max()
            if max_edge_id is None: # Handle case of empty but schema-defined df
                max_edge_id = -1
        else:
            max_edge_id = -1  # Start from 0 if no edges exist
            
        current_new_edge_id = max_edge_id + 1
        
        # Define the column structure for new edges, matching the original Edges.csv
        # Using a schema dict for polars
        edge_schema = {
            'edge_id': pl.Int64, 
            'node_id_source': pl.Int64, 
            'node_id_target': pl.Int64, 
            'is_gene_DSPP': pl.Int64, 
            'is_gene_CP': pl.Int64, 
            'is_gene_USPP': pl.Int64, 
            'is_loop': pl.Int64, 
            'is_TAD': pl.Int64, 
            'is_overlap': pl.Int64
        }

        # Pre-filter nodes to only include genes and enhancers for efficiency
        # This avoids filtering the entire node list inside the loop
        genes_and_enhancers_df = nodes_df.filter(
            pl.col('node_type').is_in(['gene', 'enhancer'])
        )

        print("Processing TADs to find gene-enhancer pairs...")
        
        # Iterate over each TAD definition
        # Using iter_rows(named=True) as a replacement for pandas iterrows()
        for tad_row in tad_df.iter_rows(named=True):
            # Get the genomic range of the current TAD
            # The .bedpe header has '#chr1' as the column name
            tad_chr = tad_row['#chr1']
            tad_start = tad_row['x1']
            tad_end = tad_row['x2']

            # Find all genes and enhancers located *within* this TAD's range
            # We check:
            # 1. Chromosome matches
            # 2. Node 'start' is greater than or equal to TAD 'start'
            # 3. Node 'end' is less than or equal to TAD 'end'
            nodes_in_tad = genes_and_enhancers_df.filter(
                (pl.col('chr') == tad_chr) &
                (pl.col('start') >= tad_start) &
                (pl.col('end') <= tad_end)
            )

            # Separate the found nodes into two groups
            genes_in_tad = nodes_in_tad.filter(pl.col('node_type') == 'gene')
            enhancers_in_tad = nodes_in_tad.filter(pl.col('node_type') == 'enhancer')

            # If there are both genes AND enhancers in this TAD, create new edges
            if not genes_in_tad.is_empty() and not enhancers_in_tad.is_empty():
                
                # Use itertools.product to create all possible (gene, enhancer) pairs
                # .to_dicts() is the polars equivalent of .to_dict('records')
                for gene, enhancer in product(genes_in_tad.to_dicts(), enhancers_in_tad.to_dicts()):
                    
                    # Create the new edge row as a dictionary
                    new_edge = {
                        'edge_id': current_new_edge_id,
                        'node_id_source': gene['node_id'],    # Gene node_id
                        'node_id_target': enhancer['node_id'], # Enhancer node_id
                        'is_gene_DSPP': 0,
                        'is_gene_CP': 0,
                        'is_gene_USPP': 0,
                        'is_loop': 0,
                        'is_TAD': 1,  # Mark this as a TAD-derived edge
                        'is_overlap': 0
                    }
                    
                    new_edges_list.append(new_edge)
                    current_new_edge_id += 1  # Increment for the next new edge

        print(f"Generated {len(new_edges_list)} new TAD-based edges.")

        # --- 5. Combine and Export Edges ---
        
        if new_edges_list:
            # Convert the list of new edges into a DataFrame
            # Provide the schema for correct types from the start
            new_edges_df = pl.DataFrame(new_edges_list, schema=edge_schema)
            
            # Concatenate the original edges with the new edges
            updated_edges_df = pl.concat([edges_df, new_edges_df])
            
        else:
            print("No new edges were generated. Original edges will be saved.")
            updated_edges_df = edges_df

        # Export the final (updated) Edges table back to the original file path
        print(f"Saving updated edges to {edges_path}...")
        
        # Ensure the directory exists before saving
        os.makedirs(os.path.dirname(edges_path), exist_ok=True)
        
        # write_csv does not write an index by default
        updated_edges_df.write_csv(edges_path)
        
        print(f"Script finished successfully. Total edges: {len(updated_edges_df)}")

    except FileNotFoundError as e:
        print(f"Error: File not found.")
        print(f"Details: {e}", file=sys.stderr)
    except pl.exceptions.NoDataError as e:
        print(f"Error: One of the input files is empty or has no data.")
        print(f"Details: {e}", file=sys.stderr)
    except (pl.exceptions.ColumnNotFoundError, KeyError) as e:
        print(f"Error: A required column is missing from an input file.")
        print(f"Make sure columns (like '#chr1', 'x1', 'x2', 'node_type', 'chr', 'start', 'end') are correct.")
        print(f"Details: {e}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()