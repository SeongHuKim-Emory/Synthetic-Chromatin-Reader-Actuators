import polars as pl
import sys
import os

def main():
    """
    Main function to process genomic node and GeneHancer files.
    """
    
    # --- File Paths ---
    # Input files
    nodes_file = '../../public/Multiomics/Step25_HAN_UpDEG_Prediction/Nodes_with_gene_enhancer_promoter.csv'
    genehancer_file = '../../public/Public Dataset/GeneHancer/GeneHancer_AnnotSV_elements_v5.25.txt'
    
    # Output files
    output_pro_en_file = '../../public/Multiomics/Step25_HAN_UpDEG_Prediction/Pro_En_Nodes.csv'
    output_combined_file = '../../public/Multiomics/Step25_HAN_UpDEG_Prediction/Nodes_with_gene_enhancer_promoter_pro_en.csv'

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_pro_en_file), exist_ok=True)

    try:
        # --- Step 1: Import files ---
        print(f"Step 1: Importing files...")
        
        # Import Nodes_with_gene_enhancer_promoter.csv
        nodes_df = pl.read_csv(nodes_file)
        print(f"  Loaded '{nodes_file}': {nodes_df.shape}")
        print(nodes_df.head(3))

        # Import GeneHancer_AnnotSV_elements_v5.25.txt (it's tab-separated)
        genehancer_df = pl.read_csv(
            genehancer_file, 
            separator='\t',
            # Force 'chr' column to be read as String to handle 'X', 'Y' etc.
            schema_overrides={'chr': pl.String}
        )
        print(f"  Loaded '{genehancer_file}': {genehancer_df.shape}")
        print(genehancer_df.head(3))

        # --- Step 2 & 3: Filter GeneHancer for 'Promoter/Enhancer' and create Pro_En_Nodes ---
        print("\nStep 2 & 3: Filtering for 'Promoter/Enhancer' and creating 'Pro_En_Nodes' table...")
        
        pro_en_nodes_df = genehancer_df.filter(
            pl.col('regulatory_element_type') == 'Promoter/Enhancer'
        ).select(
            pl.lit(0).cast(pl.Int64).alias('node_id'),  # Step 3: Start with node_id 0
            pl.col('GHid').alias('identifier'),
            pl.col('regulatory_element_type').alias('node_type'),
            pl.col('chr'),
            pl.col('element_start').alias('start'),
            pl.col('element_end').alias('end'),
            pl.col('is_elite'),
            pl.col('enhancer_score')
        )
        
        print(f"  Found {pro_en_nodes_df.height} 'Promoter/Enhancer' rows.")
        print("  'Pro_En_Nodes' table (before 'chr' and 'node_id' update):")
        print(pro_en_nodes_df.head(3))

        # --- Step 4: Add "chr" in front of the "chr" column value ---
        print("\nStep 4: Adding 'chr' prefix to 'chr' column...")
        
        pro_en_nodes_df = pro_en_nodes_df.with_columns(
            pl.concat_str([
                pl.lit('chr'), 
                pl.col('chr').cast(pl.String)
            ]).alias('chr')
        )
        
        print("  'Pro_En_Nodes' table (after 'chr' update):")
        print(pro_en_nodes_df.head(3))

        # --- Step 5: Update node_id starting from 367313 ---
        print("\nStep 5: Updating 'node_id' to start from 367313...")
        start_id = 367313
        
        # Use pl.arange (lazily) and pl.len() to create the new sequential IDs
        pro_en_nodes_df = pro_en_nodes_df.with_columns(
            pl.arange(start_id, start_id + pl.len()).alias('node_id')
        )
        
        # Ensure final column order
        pro_en_nodes_df = pro_en_nodes_df.select([
            'node_id', 'identifier', 'node_type', 'chr', 'start', 'end', 
            'is_elite', 'enhancer_score'
        ])

        print("  'Pro_En_Nodes' table (after 'node_id' update):")
        print(pro_en_nodes_df.head(3))
        print(pro_en_nodes_df.tail(3))

        # --- Step 6: Export the Pro_En_Nodes table ---
        print(f"\nStep 6: Exporting 'Pro_En_Nodes' to '{output_pro_en_file}'...")
        pro_en_nodes_df.write_csv(output_pro_en_file)
        print("  Export complete.")

        # --- Step 7: Create Nodes_with_gene_enhancer_promoter_pro_en table ---
        # This is just the 'nodes_df' loaded in Step 1, which we will append to.
        print(f"\nStep 7: Using '{nodes_file}' as base for 'Nodes_with_gene_enhancer_promoter_pro_en'.")
        
        # --- Step 8: Merge rows from Pro_En_Nodes to the main nodes table ---
        print("Step 8: Appending new 'Pro_En_Nodes' (subset of columns) to base table...")
        
        # Select only the columns needed for the final combined file
        columns_to_merge = ['node_id', 'identifier', 'node_type', 'chr', 'start', 'end']
        pro_en_nodes_subset = pro_en_nodes_df.select(columns_to_merge)
        
        print(f"  Base table shape: {nodes_df.shape}")
        print(f"  New nodes to append: {pro_en_nodes_subset.shape}")

        # Ensure data types are compatible before concatenation (e.g., node_id)
        # nodes_df node_id is likely int, pro_en_nodes_subset node_id is int.
        # start/end in nodes_df are int, start/end in subset are int. This should be fine.
        
        # Vertically stack the two DataFrames
        nodes_combined_df = pl.concat([nodes_df, pro_en_nodes_subset])
        
        print(f"  Final combined table shape: {nodes_combined_df.shape}")
        print("  Tail of combined table:")
        print(nodes_combined_df.tail(5))

        # --- Step 9: Export Nodes_with_gene_enhancer_promoter_pro_en table ---
        print(f"\nStep 9: Exporting combined table to '{output_combined_file}'...")
        nodes_combined_df.write_csv(output_combined_file)
        print("  Export complete.")
        
        print("\n--- All steps completed successfully. ---")

    except pl.exceptions.ShapeError as e:
        print(f"\nError: DataFrame shape mismatch during operation. {e}", file=sys.stderr)
    except FileNotFoundError as e:
        print(f"\nError: File not found. {e}", file=sys.stderr)
        print("Please check that the input file paths are correct.", file=sys.stderr)
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
