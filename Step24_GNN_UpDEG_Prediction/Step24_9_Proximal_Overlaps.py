import polars as pl
import sys
import os

"""
This script loads node and edge data for a genomic graph using Polars.
Its main purpose is to find new 'proximal overlap' edges between
enhancer nodes and promoter nodes (promoter_core, promoter_proximal_downstream,
promoter_proximal_upstream) that are within a 1kb leniency.
It adds these new edges to the edge list and exports the updated file.
"""

# --- Configuration ---
# Define file paths
BASE_DIR = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction'
NODE_FILE_PATH = os.path.join(BASE_DIR, 'Nodes.csv')
EDGE_FILE_PATH = os.path.join(BASE_DIR, 'Edges.csv')
OUTPUT_FILE_PATH = os.path.join(BASE_DIR, 'Edges.csv') # Overwrites the original

# Define overlap parameters
LENIENCY_KB = 1000 # 1kb as specified
PROMOTER_TYPES = [
    'promoter_proximal_downstream',
    'promoter_core',
    'promoter_proximal_upstream'
]

def load_data(file_path: str) -> pl.DataFrame:
    """Loads a CSV file into a polars DataFrame."""
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}", file=sys.stderr)
        sys.exit(1)
    try:
        return pl.read_csv(file_path)
    except Exception as e:
        print(f"Error loading {file_path} with polars: {e}", file=sys.stderr)
        sys.exit(1)

def find_overlaps(nodes_df: pl.DataFrame, edges_df: pl.DataFrame) -> pl.DataFrame | None:
    """
    Finds proximal overlaps between enhancers and promoters
    and returns a DataFrame of new edges.
    """
    print("Finding enhancer-promoter overlaps...")
    
    # 5. Filter for enhancers
    enhancers = nodes_df.filter(pl.col('node_type') == 'enhancer')
    
    # 5. Filter for promoter types
    promoters = nodes_df.filter(pl.col('node_type').is_in(PROMOTER_TYPES))
    
    print(f"Found {enhancers.height} enhancers and {promoters.height} total promoters.")

    # Ensure 'chr' is string for correct matching (e.g., 'chr1' != 'chr10')
    enhancers = enhancers.with_columns(pl.col('chr').cast(pl.String))
    promoters = promoters.with_columns(pl.col('chr').cast(pl.String))

    # Apply leniency to enhancer ranges for comparison
    enhancers = enhancers.with_columns(
        start_lenient = pl.col('start') - LENIENCY_KB,
        end_lenient = pl.col('end') + LENIENCY_KB
    )

    # Get the current highest edge_id to start incrementing from
    if edges_df.height == 0:
        current_max_edge_id = -1
    else:
        # .item() safely extracts the single value
        current_max_edge_id = edges_df.select(pl.max('edge_id')).item()

    # Join enhancers and promoters on the same chromosome
    overlapping_pairs = enhancers.join(promoters, on='chr', suffix='_promo')
    
    # Filter for actual genomic overlaps
    # Overlap condition: (promoter_start <= enhancer_end) AND (promoter_end >= enhancer_start)
    overlapping_pairs = overlapping_pairs.filter(
        (pl.col('start_promo') <= pl.col('end_lenient')) &
        (pl.col('end_promo') >= pl.col('start_lenient'))
    )
    
    num_new_overlaps = overlapping_pairs.height
    print(f"Found {num_new_overlaps} new proximal overlaps.")
    
    if num_new_overlaps == 0:
        return None
    
    # Create the new edge rows
    # We use .with_row_count() to generate sequential IDs
    new_edges_df = overlapping_pairs.with_row_count(name="row_num").select(
        (pl.col("row_num") + current_max_edge_id + 1).cast(pl.Int64).alias("edge_id"),
        pl.col("node_id").alias("node_id_source"),
        pl.col("node_id_promo").alias("node_id_target"),
        pl.lit(0).cast(pl.Int64).alias("is_gene_DSPP"),
        pl.lit(0).cast(pl.Int64).alias("is_gene_CP"),
        pl.lit(0).cast(pl.Int64).alias("is_gene_USPP"),
        pl.lit(0).cast(pl.Int64).alias("is_loop"),
        pl.lit(0).cast(pl.Int64).alias("is_TAD"),
        pl.lit(0).cast(pl.Int64).alias("is_overlap"),
        pl.lit(1).cast(pl.Int64).alias("is_proximal_overlap")
    )
    
    # Ensure new edges DataFrame has columns in the same order as the original
    # This helps ensure the concat operation is clean, though polars concat
    # can align by name.
    column_order = list(edges_df.columns)
    if 'is_proximal_overlap' not in column_order:
         column_order = list(edges_df.columns) + ['is_proximal_overlap']
         
    new_edges_df = new_edges_df.select(column_order)

    return new_edges_df

def main():
    print("Starting script...")

    # 1. Import files
    print(f"Loading Nodes from {NODE_FILE_PATH}...")
    nodes_df = load_data(NODE_FILE_PATH) # 4. Create Nodes table
    print(f"Loading Edges from {EDGE_FILE_PATH}...")
    edges_df = load_data(EDGE_FILE_PATH) # 2. Create Edges table

    print(f"Original Edges table has {edges_df.height} rows.")

    # 3. Add a new "is_proximal_overlap" column
    # Ensure all existing integer columns and the new one are Int64
    # to prevent type mismatches on concatenation.
    edges_df = edges_df.with_columns(
        pl.col(pl.Int64), # Cast all existing integer columns to Int64
        is_proximal_overlap=pl.lit(0).cast(pl.Int64)
    )
    print("Added 'is_proximal_overlap' column to existing edges.")

    # 5. Find new overlaps
    new_edges_df = find_overlaps(nodes_df, edges_df)

    # Add new edges to the main Edges DataFrame
    if new_edges_df is not None:
        # polars.concat stacks DataFrames vertically
        edges_df = pl.concat([edges_df, new_edges_df])
        print(f"Added {new_edges_df.height} new edges to the table.")

    # 6. Export the Edges table
    try:
        edges_df.write_csv(OUTPUT_FILE_PATH)
        print(f"\nSuccessfully exported updated Edges table ({edges_df.height} total rows) to {OUTPUT_FILE_PATH}")
    except IOError as e:
        print(f"Error writing to file: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during file writing: {e}", file=sys.stderr)
        sys.exit(1)

    print("Script finished.")

if __name__ == "__main__":
    # Ensure polars is installed
    try:
        import polars as pl
    except ImportError:
        print("Error: The 'polars' library is required. Please install it using 'pip install polars'", file=sys.stderr)
        sys.exit(1)
        
    main()

