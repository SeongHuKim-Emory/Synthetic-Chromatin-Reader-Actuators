import polars as pl
import sys

# --- Script Configuration ---
# (Paths remain the same)
UPDEG_PATH = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone_RNA.csv'
NON_DEG_PATH = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone_RNA.csv'
GO_ANNOTATION_PATH = '../../public/Public Dataset/GO/goa_human header trimmed.csv'

GO_MATRIX_OUTPUT_PATH = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GO_Matrix.csv'

# FIX: Shorten output file paths to avoid OS path length limits.
UPDEG_GO_OUTPUT_PATH = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/SRA_Specific_UpDEG_with_GO.csv'
NON_DEG_GO_OUTPUT_PATH = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/Non_DEG_with_GO.csv'


# --- Main Script Logic ---

try:
    print(f"Using Polars version: {pl.__version__}")

    # OPTIMIZATION: Use a StringCache context to ensure all
    # categorical columns ('symbol', 'GO_ID') from different
    # files can be compared and joined.
    with pl.StringCache():
        # 1. Import files using Polars for speed
        print("Step 1: Importing CSV files using Polars (multi-threaded)...")
        
        # We set symbol to Categorical for faster joins later.
        updeg_df = pl.read_csv(UPDEG_PATH, dtypes={'symbol': pl.Categorical})
        non_deg_df = pl.read_csv(NON_DEG_PATH, dtypes={'symbol': pl.Categorical})
        
        # For go_df, we set both key columns to Categorical.
        go_df = pl.read_csv(GO_ANNOTATION_PATH, dtypes={'symbol': pl.Categorical, 'GO_ID': pl.Categorical})
        
        print("Files imported successfully.")

        # 2. Create the GO_Matrix (Polars-Optimized)
        print("\nStep 2: Creating the GO_Matrix using Polars 'pivot'...")

        print("Dropping duplicate symbol-GO_ID pairs...")
        # Get unique pairs, as we only need to mark (symbol, GO_ID) as 1, not count them.
        go_unique_pairs = go_df.select(['symbol', 'GO_ID']).unique()
        
        # OPTIMIZATION: Free memory
        del go_df

        print("Creating GO_Matrix with pivot operation...")
        # OPTIMIZATION: Polars' pivot is the multi-threaded equivalent of
        # crosstab() or building a dense matrix.
        # We create a new column of literal 1s to use as the value.
        go_matrix = go_unique_pairs.with_columns(
            pl.lit(1).cast(pl.Int8).alias("value")
        ).pivot(
            index='symbol',
            columns='GO_ID',
            values='value'
        )
        
        # We need the list of all GO ID columns for the fill_null step later
        # The first column is 'symbol', so we skip it.
        all_go_ids = go_matrix.columns[1:]

        print("GO_Matrix created successfully.")

        # 3. Export the GO_Matrix
        print(f"\nStep 3: Exporting GO_Matrix to '{GO_MATRIX_OUTPUT_PATH}' using Polars...")
        
        # OPTIMIZATION: Polars' write_csv is multi-threaded.
        go_matrix.write_csv(GO_MATRIX_OUTPUT_PATH)
        
        print("GO_Matrix exported.")

        # 4. Merge GO data into the main dataframes
        print("\nStep 4: Merging GO data into the UpDEG and Non-DEG tables using Polars...")

        # OPTIMIZATION: Polars' joins are multi-threaded and highly optimized.
        # We can create the lookup table much more efficiently.
        
        print("Finding all unique symbols across both dataframes...")
        unique_symbols_updeg = updeg_df.select('symbol').unique()
        unique_symbols_nondeg = non_deg_df.select('symbol').unique()
        
        # Combine and find unique symbols
        # This pl.concat is where the error occurred, and it is now
        # protected by the StringCache context.
        all_unique_symbols_df = pl.concat(
            [unique_symbols_updeg, unique_symbols_nondeg],
            how='vertical' # Use 'vertical' for clarity
        ).unique()
        
        del unique_symbols_updeg
        del unique_symbols_nondeg

        print(f"Found {all_unique_symbols_df.height} unique symbols in total.")
        
        print("Creating dense GO lookup table (joining and filling nulls)...")
        # Join the unique symbols with the matrix. This creates nulls.
        # Then, use fill_null(0) to replace all nulls with 0.
        # This is the Polars equivalent of the pandas Step 4, but much faster.
        symbols_go_lookup = all_unique_symbols_df.join(
            go_matrix, on='symbol', how='left'
        ).fill_null(0)

        # We can also cast all GO columns to Int8 at once for memory safety
        symbols_go_lookup = symbols_go_lookup.with_columns(
            pl.col(all_go_ids).cast(pl.Int8)
        )
        
        print("Lookup table created successfully (dense, int8).")
        
        # OPTIMIZATION: Free up go_matrix memory
        del go_matrix
        del all_unique_symbols_df

        # ---
        # Now, perform the final joins. These are multi-threaded.
        # ---

        print("Joining dense GO lookup table with UpDEG dataframe (using 'inner' join)...")
        updeg_go_df = updeg_df.join(symbols_go_lookup, on='symbol', how='inner')
        
        print("Joining dense GO lookup table with Non-DEG dataframe (using 'inner' join)...")
        # FIX: Ensure the join uses the original 'non_deg_df' (on the right)
        # to create the new 'non_deg_go_df' (on the left).
        # A typo here (e.g., non_deg_go_df.join) would cause the NameError you saw.
        non_deg_go_df = non_deg_df.join(symbols_go_lookup, on='symbol', how='inner')

        # OPTIMIZATION: Free up memory
        print("Freeing memory from original dataframes...")
        del symbols_go_lookup
        del updeg_df
        del non_deg_df
        
        print("Merge and optimization complete.")

        # 5. Export the final merged tables using Polars
        print(f"\nStep 5: Exporting final tables using Polars...")

        # OPTIMIZATION: Polars' write_csv is multi-threaded.
        print(f"Writing UpDEG table to CSV: {UPDEG_GO_OUTPUT_PATH}")
        updeg_go_df.write_csv(UPDEG_GO_OUTPUT_PATH)
        del updeg_go_df # Free memory
        print(f"-> Saved: {UPDEG_GO_OUTPUT_PATH}")

        print(f"Writing Non-DEG table to CSV: {NON_DEG_GO_OUTPUT_PATH}")
        non_deg_go_df.write_csv(NON_DEG_GO_OUTPUT_PATH)
        del non_deg_go_df # Free memory
        print(f"-> Saved: {NON_DEG_GO_OUTPUT_PATH}")

        print("\nAll tasks completed successfully!")

except ImportError:
    print(f"\nError: The 'polars' library is required for this optimized script.", file=sys.stderr)
    print("Please install it by running: pip install polars", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"\nAn unexpected error occurred: {e}", file=sys.stderr)

