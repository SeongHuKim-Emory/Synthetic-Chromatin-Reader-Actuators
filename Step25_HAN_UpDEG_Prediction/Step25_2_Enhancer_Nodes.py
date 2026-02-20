import polars as pl
from pathlib import Path
import sys

def main():
    """
    Main function to process gene and enhancer nodes as per the request.
    """
    print(f"--- Script Start ---")
    print(f"Python version: {sys.version.split()[0]}")
    print(f"Polars version: {pl.__version__}")

    # --- Path Definitions ---
    # Define paths relative to the script's location for robustness
    try:
        # This works when running as a script
        SCRIPT_DIR = Path(__file__).parent.resolve()
    except NameError:
        # Fallback for interactive environments (e.g., Jupyter, REPL)
        print("Warning: __file__ not defined. Using current working directory as base.")
        SCRIPT_DIR = Path.cwd()

    # 1. Import files
    NODES_GENE_PATH = SCRIPT_DIR / "../../public/Multiomics/Step25_HAN_UpDEG_Prediction/Nodes_with_gene.csv"
    GENEHANCER_PATH = SCRIPT_DIR / "../../public/Public Dataset/GeneHancer/GeneHancer_AnnotSV_elements_v5.25.txt"
    
    # Define output paths
    OUTPUT_DIR = SCRIPT_DIR / "../../public/Multiomics/Step25_HAN_UpDEG_Prediction"
    ENHANCER_NODES_PATH = OUTPUT_DIR / "Enhancer_Nodes.csv"
    NODES_GENE_ENHANCER_PATH = OUTPUT_DIR / "Nodes_with_gene_enhancer.csv"

    # Ensure output directory exists
    try:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"Error creating output directory {OUTPUT_DIR}: {e}", file=sys.stderr)
        return

    try:
        print(f"\n[Step 1] Loading files...")
        print(f"  > Loading gene nodes: {NODES_GENE_PATH}")
        nodes_gene_df = pl.read_csv(NODES_GENE_PATH)
        print(f"  > Loaded {len(nodes_gene_df):,} rows from Nodes_with_gene.csv")

        print(f"  > Loading GeneHancer elements: {GENEHANCER_PATH}")
        genehancer_df = pl.read_csv(
            GENEHANCER_PATH, 
            separator='\t',
            # Add schema_overrides to correctly read 'chr' as String
            # This fixes the error caused by non-numeric values like 'X', 'Y', 'M'
            schema_overrides={"chr": pl.String}
        )
        print(f"  > Loaded {len(genehancer_df):,} rows from GeneHancer_AnnotSV_elements_v5.25.txt")

    except FileNotFoundError as e:
        print(f"\n[Error] File not found: {e.filename}", file=sys.stderr)
        print("Please ensure the file paths are correct relative to the script's location.", file=sys.stderr)
        return
    except Exception as e:
        print(f"\n[Error] An error occurred during file loading: {e}", file=sys.stderr)
        return

    # [Step 2 & 3] Filter for 'Enhancer' elements and create the base Enhancer_Nodes table
    print("\n[Step 2 & 3] Filtering for 'Enhancer' elements and creating Enhancer_Nodes table...")
    enhancer_nodes_df = genehancer_df.filter(
        pl.col("regulatory_element_type") == "Enhancer"
    ).select(
        pl.lit(0).cast(pl.Int64).alias("node_id"),  # Add placeholder node_id (Step 3)
        pl.col("GHid").alias("identifier"),
        pl.col("regulatory_element_type").alias("node_type"),
        pl.col("chr"),
        pl.col("element_start").alias("start"),
        pl.col("element_end").alias("end"),
        pl.col("is_elite"),
        pl.col("enhancer_score")
    )
    print(f"  > Found {len(enhancer_nodes_df):,} 'Enhancer' elements.")

    # [Step 4] Add "chr" prefix to the 'chr' column
    print("\n[Step 4] Adding 'chr' prefix to 'chr' column...")
    enhancer_nodes_df = enhancer_nodes_df.with_columns(
        # Concatenate the string "chr" with the 'chr' column cast to string
        pl.concat_str([pl.lit("chr"), pl.col("chr").cast(pl.String)]).alias("chr")
    )
    print("  > 'chr' column updated.")

    # [Step 5] Assign new, sequential node_ids
    print("\n[Step 5] Assigning new node_ids starting from 44219 (or max_gene_id + 1)...")
    
    # Find the maximum node_id from the gene nodes file
    max_gene_node_id = nodes_gene_df.select(pl.col("node_id").max()).item()
    start_id = max_gene_node_id + 1
    
    print(f"  > Max gene node_id: {max_gene_node_id}")
    print(f"  > New enhancer node_ids will start from: {start_id}")
    
    # Generate the new range of node_ids
    new_node_ids = pl.Series("node_id", range(start_id, start_id + len(enhancer_nodes_df)))
    
    # Update the node_id column
    enhancer_nodes_df = enhancer_nodes_df.with_columns(new_node_ids)
    print(f"  > Assigned node_ids from {start_id} to {start_id + len(enhancer_nodes_df) - 1}.")

    # [Step 6] Export the Enhancer_Nodes table
    print(f"\n[Step 6] Exporting Enhancer_Nodes to: {ENHANCER_NODES_PATH}")
    try:
        enhancer_nodes_df.write_csv(ENHANCER_NODES_PATH)
        print(f"  > Successfully exported {ENHANCER_NODES_PATH}")
    except Exception as e:
        print(f"  > [Error] Failed to export Enhancer_Nodes: {e}", file=sys.stderr)

    # [Step 7] Create Nodes_with_gene_enhancer table (starts with gene nodes)
    # We will use the already-loaded 'nodes_gene_df' as the base for concatenation
    print("\n[Step 7 & 8] Merging gene nodes and enhancer nodes...")

    # [Step 8] Merge (concatenate) all rows from Enhancer_Nodes
    # First, select only the columns that match the Nodes_with_gene structure
    enhancer_nodes_to_merge = enhancer_nodes_df.select(
        "node_id", "identifier", "node_type", "chr", "start", "end"
    )

    # Concatenate the two dataframes
    nodes_with_gene_enhancer_df = pl.concat(
        [nodes_gene_df, enhancer_nodes_to_merge]
    )
    print(f"  > Original gene nodes: {len(nodes_gene_df):,}")
    print(f"  > Added enhancer nodes: {len(enhancer_nodes_to_merge):,}")
    print(f"  > Total nodes in merged table: {len(nodes_with_gene_enhancer_df):,}")

    # [Step 9] Export Nodes_with_gene_enhancer table
    print(f"\n[Step 9] Exporting Nodes_with_gene_enhancer to: {NODES_GENE_ENHANCER_PATH}")
    try:
        nodes_with_gene_enhancer_df.write_csv(NODES_GENE_ENHANCER_PATH)
        print(f"  > Successfully exported {NODES_GENE_ENHANCER_PATH}")
    except Exception as e:
        print(f"  > [Error] Failed to export Nodes_with_gene_enhancer: {e}", file=sys.stderr)

    print("\n--- Script Finished ---")

if __name__ == "__main__":
    main()
