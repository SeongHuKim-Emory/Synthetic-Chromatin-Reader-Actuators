import polars as pl
import sys

def main():
    """
    Main function to process genomic node and annotation files.
    """
    
    # Define file paths
    nodes_file_path = "../../public/Multiomics/Step25_HAN_UpDEG_Prediction/Nodes_with_gene_enhancer.csv"
    genehancer_file_path = "../../public/Public Dataset/GeneHancer/GeneHancer_AnnotSV_elements_v5.25.txt"
    promoter_output_path = "../../public/Multiomics/Step25_HAN_UpDEG_Prediction/Promoter_Nodes.csv"
    final_output_path = "../../public/Multiomics/Step25_HAN_UpDEG_Prediction/Nodes_with_gene_enhancer_promoter.csv"

    # --- Step 1: Import files ---
    print("Step 1: Importing files...")
    try:
        # Import the existing nodes file
        nodes_df = pl.read_csv(nodes_file_path)
        print(f"Successfully loaded {nodes_file_path}")
        print(f"Initial nodes shape: {nodes_df.shape}")

        # Import the GeneHancer annotation file (it's tab-separated)
        # We add schema_overrides={"chr": pl.Utf8} to handle non-numeric chromosome names like 'X' or 'Y'
        genehancer_df = pl.read_csv(
            genehancer_file_path, 
            separator='\t', 
            schema_overrides={"chr": pl.Utf8}
        )
        print(f"Successfully loaded {genehancer_file_path}")
        print(f"GeneHancer annotations shape: {genehancer_df.shape}")

    except FileNotFoundError as e:
        print(f"Error: File not found. {e}", file=sys.stderr)
        print("Please check your file paths.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred during file loading: {e}", file=sys.stderr)
        sys.exit(1)

    # --- Step 2 & 3: Create Promoter_Nodes table and populate it ---
    print("\nStep 2 & 3: Filtering for 'Promoter' elements and creating Promoter_Nodes table...")
    
    # Filter for rows where regulatory_element_type is exactly "Promoter"
    promoter_nodes_df = genehancer_df.filter(
        pl.col("regulatory_element_type") == "Promoter"
    )

    # Select and rename columns as specified
    promoter_nodes_df = promoter_nodes_df.select(
        pl.lit(0).alias("node_id"),  # Set node_id to 0 for now
        pl.col("GHid").alias("identifier"),
        pl.col("regulatory_element_type").alias("node_type"),
        pl.col("chr"),
        pl.col("element_start").alias("start"),
        pl.col("element_end").alias("end"),
        pl.col("is_elite"),
        pl.col("enhancer_score")
    )
    
    print(f"Found {promoter_nodes_df.height} 'Promoter' elements.")

    # --- Step 4: Add "chr" prefix to the 'chr' column ---
    print("\nStep 4: Adding 'chr' prefix to chromosome values...")
    promoter_nodes_df = promoter_nodes_df.with_columns(
        # Cast 'chr' column to string (Utf8) and concatenate "chr" at the start
        pl.concat_str([pl.lit("chr"), pl.col("chr").cast(pl.Utf8)]).alias("chr")
    )

    # --- Step 5: Assign new node_id values ---
    print("\nStep 5: Assigning new sequential node_ids starting from 360264...")
    start_id = 360264
    
    promoter_nodes_df = promoter_nodes_df.with_columns(
        # Create a new range of IDs starting from start_id
        pl.arange(start_id, start_id + pl.count()).alias("node_id")
    )

    print("Promoter_Nodes table after ID assignment (first 5 rows):")
    print(promoter_nodes_df.head())

    # --- Step 6: Export the Promoter_Nodes table ---
    print(f"\nStep 6: Exporting Promoter_Nodes table to {promoter_output_path}...")
    try:
        promoter_nodes_df.write_csv(promoter_output_path)
        print("Export successful.")
    except Exception as e:
        print(f"Error exporting Promoter_Nodes.csv: {e}", file=sys.stderr)

    # --- Step 7: Make Nodes_with_gene_enhancer_promoter table ---
    # We already have nodes_df in memory, which is the content of 
    # Nodes_with_gene_enhancer.csv (as per your file import).
    # We will use this as the base for our combined table.
    print("\nStep 7: Using loaded Nodes_with_gene_enhancer.csv as base...")
    nodes_with_all_df = nodes_df

    # --- Step 8: Merge (Append) Promoter_Nodes to the main table ---
    print("\nStep 8: Appending new promoter nodes to the main node table...")
    
    # Select only the columns needed for the final combined file
    promoter_nodes_to_append = promoter_nodes_df.select(
        "node_id", "identifier", "node_type", "chr", "start", "end"
    )
    
    # Concatenate the original nodes with the new promoter nodes
    nodes_with_all_df = pl.concat([nodes_with_all_df, promoter_nodes_to_append])

    print(f"New combined table shape: {nodes_with_all_df.shape}")
    print("Combined table (last 5 rows):")
    print(nodes_with_all_df.tail())

    # --- Step 9: Export Nodes_with_gene_enhancer_promoter table ---
    print(f"\nStep 9: Exporting final table to {final_output_path}...")
    try:
        nodes_with_all_df.write_csv(final_output_path)
        print("Export successful.")
    except Exception as e:
        print(f"Error exporting Nodes_with_gene_enhancer_promoter.csv: {e}", file=sys.stderr)

    print("\n--- All steps completed. ---")

if __name__ == "__main__":
    main()
