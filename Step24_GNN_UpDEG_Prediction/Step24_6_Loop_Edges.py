import polars as pl
from pathlib import Path
import sys

def find_genomic_overlaps(nodes_path, edges_path, output_path):
    """
    Finds genomic overlaps between enhancer nodes and promoter nodes and
    adds new edges for these overlaps using Polars.

    Args:
        nodes_path (Path): Path to the Nodes.csv file.
        edges_path (Path): Path to the Edges.csv file.
        output_path (Path): Path to save the updated Edges.csv file.
    """
    try:
        # 1. Import files & 2/3. Create tables (DataFrames)
        print(f"Loading nodes from: {nodes_path}")
        # Eagerly load the data and ensure correct types for coordinates
        # FIX: Use schema_overrides instead of deprecated dtypes
        nodes_df = pl.read_csv(
            nodes_path,
            schema_overrides={'start': pl.Int64, 'end': pl.Int64}
        )
        
        print(f"Loading edges from: {edges_path}")
        edges_df = pl.read_csv(edges_path)
        
        print(f"Original nodes shape: {nodes_df.shape}")
        print(f"Original edges shape: {edges_df.shape}")

    except FileNotFoundError as e:
        print(f"Error: File not found. {e}", file=sys.stderr)
        return
    except Exception as e:
        print(f"Error loading files: {e}", file=sys.stderr)
        return

    try:
        # 4. Find overlaps
        print("Separating enhancers and promoters...")
        
        # Get all enhancer rows and rename columns for joining
        enhancers = (
            nodes_df
            .filter(pl.col('node_type') == 'enhancer')
            .rename({
                'node_id': 'node_id_enh',
                'start': 'start_enh',
                'end': 'end_enh'
            })
        )
        
        # Get all promoter rows and rename columns for joining
        promoter_types = [
            'promoter_proximal_downstream',
            'promoter_core',
            'promoter_proximal_upstream'
        ]
        promoters = (
            nodes_df
            .filter(pl.col('node_type').is_in(promoter_types))
            .rename({
                'node_id': 'node_id_pro',
                'start': 'start_pro',
                'end': 'end_pro'
            })
            # Select only necessary columns for the join
            .select(['chr', 'node_id_pro', 'start_pro', 'end_pro'])
        )
        
        print(f"Found {len(enhancers)} enhancers and {len(promoters)} promoters.")
        
        if enhancers.is_empty() or promoters.is_empty():
            print("No enhancers or promoters found. No new edges to add.")
            updated_edges_df = edges_df
            
        else:
            # Use polars join to efficiently find all pairs on the same chromosome
            print("Finding potential overlaps by joining on 'chr'...")
            merged_df = enhancers.join(promoters, on='chr')
            
            print(f"Found {len(merged_df)} potential chr pairs.")

            # Apply the genomic overlap condition:
            # (StartA <= EndB) and (EndA >= StartB)
            overlap_condition = (
                (pl.col('start_enh') <= pl.col('end_pro')) &
                (pl.col('end_enh') >= pl.col('start_pro'))
            )
            
            overlaps = merged_df.filter(overlap_condition)
            
            num_overlaps = len(overlaps)
            print(f"Found {num_overlaps} actual genomic overlaps.")
            
            if num_overlaps > 0:
                # Create the new rows for the Edges table
                print("Creating new edges...")
                
                # Get the last edge_id
                # .max().item() safely extracts the single value
                last_edge_id = edges_df['edge_id'].max()
                if last_edge_id is None:
                    last_edge_id = -1 # Start from 0 if edges_df was empty
                
                # Generate new, sequential edge_ids
                new_ids = pl.arange(
                    last_edge_id + 1,
                    last_edge_id + 1 + num_overlaps,
                    eager=True
                )
                
                new_edges_df = (
                    overlaps.select([
                        pl.lit(None).alias('edge_id'), # Placeholder
                        pl.col('node_id_enh').alias('node_id_source'),
                        pl.col('node_id_pro').alias('node_id_target'),
                    ])
                    .with_columns(
                        pl.lit(0).alias('is_gene_DSPP'),
                        pl.lit(0).alias('is_gene_CP'),
                        pl.lit(0).alias('is_gene_USPP'),
                        pl.lit(0).alias('is_loop'),
                        pl.lit(0).alias('is_TAD'),
                        pl.lit(1).alias('is_overlap'),
                        new_ids.alias('edge_id'), # Add the real IDs
                    )
                    # Reorder columns to match the original edges_df
                    .select(edges_df.columns)
                )

                # --- FIX: Cast columns before concat to prevent SchemaError ---
                # Define the integer columns that must match
                int_cols = [
                    'edge_id', 'node_id_source', 'node_id_target', 
                    'is_gene_DSPP', 'is_gene_CP', 'is_gene_USPP', 
                    'is_loop', 'is_TAD', 'is_overlap'
                ]
                
                print("Casting columns to Int64 to ensure schema match...")
                # Cast all specified columns to Int64 in both dataframes
                edges_df = edges_df.with_columns(
                    [pl.col(c).cast(pl.Int64) for c in int_cols if c in edges_df.columns]
                )
                
                new_edges_df = new_edges_df.with_columns(
                    [pl.col(c).cast(pl.Int64) for c in int_cols if c in new_edges_df.columns]
                )
                # --- End of FIX ---

                # Append the new edges to the original edges table
                print("Concatenating DataFrames...")
                updated_edges_df = pl.concat([edges_df, new_edges_df])
                
                # This block is now slightly redundant but ensures final type safety
                updated_edges_df = updated_edges_df.with_columns(
                    [pl.col(c).cast(pl.Int64) for c in int_cols if c in updated_edges_df.columns]
                )
                
            else:
                print("No overlaps found, edges file will not be modified.")
                updated_edges_df = edges_df

        # 5. Export the Edges table
        print(f"Saving updated edges to: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        updated_edges_df.write_csv(output_path)
        
        print(f"Process complete. New edges shape: {updated_edges_df.shape}")

    except Exception as e:
        print(f"An error occurred during overlap processing: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()

def main():
    # Define file paths using pathlib
    base_path = Path('../../public/Multiomics/Step24_GNN_UpDEG_Prediction')
    nodes_file = base_path / 'Nodes.csv'
    edges_file = base_path / 'Edges.csv'
    
    # Run the overlap detection function
    find_genomic_overlaps(
        nodes_path=nodes_file,
        edges_path=edges_file,
        output_path=edges_file  # Overwrites the original file as requested
    )

if __name__ == "__main__":
    main()
