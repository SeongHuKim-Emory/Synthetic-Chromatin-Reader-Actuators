import polars as pl
import os

def create_edges_file():
    """
    Imports node and association files, processes them to create an edges table,
    and exports it to a CSV file.
    """
    print("Script started...")

    # --- 1. Define File Paths ---
    nodes_file = "../../public/Multiomics/Step25_HAN_UpDEG_Prediction/Nodes_with_gene_enhancer_promoter_pro_en.csv"
    assoc_file = "../../public/Public Dataset/GeneHancer/GeneHancer_AnnotSV_gene_association_scores_v5.25.txt"
    output_dir = "../../public/Multiomics/Step25_HAN_UpDEG_Prediction"
    output_file = os.path.join(output_dir, "Edges.csv")

    try:
        # --- 2. Import Files ---
        print(f"Loading nodes file: {nodes_file}")
        # Load the nodes file
        nodes_df = pl.read_csv(nodes_file)

        print(f"Loading associations file: {assoc_file}")
        # Load the GeneHancer associations file (assuming tab-separated)
        assoc_df = pl.read_csv(assoc_file, separator="\t")

        print("Files loaded successfully.")

        # --- 3. Process Data for Joins ---
        
        # Create a DataFrame for gene nodes
        # We will link this to the 'symbol' column from the associations file
        gene_nodes = (
            nodes_df.filter(pl.col("node_type") == "gene")
            .select(
                pl.col("identifier").alias("symbol"),  # Alias to match assoc_df
                pl.col("node_id").alias("node_id_1")   # This will be node_id_1
            )
        )

        # Create a DataFrame for GeneHancer (gh) nodes (Enhancer, Promoter, etc.)
        # We will link this to the 'GHid' column from the associations file
        gh_nodes = (
            nodes_df.filter(pl.col("node_type") != "gene")
            .select(
                pl.col("identifier").alias("GHid"),     # Alias to match assoc_df
                pl.col("node_id").alias("node_id_2"),   # This will be node_id_2
                pl.col("node_type")                     # Keep node_type for one-hot encoding
            )
        )

        # --- 4. Create the Edges Table ---
        print("Creating edges by joining tables...")

        # The logic "if both... exist" is handled by 'inner' joins.
        
        # Start with the associations file
        edges_df = (
            assoc_df
            # Find the corresponding gene row ("node_id_1")
            .join(
                gene_nodes, 
                on="symbol", 
                how="inner"
            )
            # Find the corresponding GeneHancer row ("node_id_2" and "node_type")
            .join(
                gh_nodes, 
                on="GHid", 
                how="inner"
            )
            # --- 5. Add new columns based on conditions ---
            .with_columns(
                # Use expressions to create the one-hot encoded columns
                # .eq() creates a boolean, .cast(pl.UInt8) converts True->1, False->0
                pl.col("node_type").eq(pl.lit("Enhancer")).cast(pl.UInt8).alias("is_gene_Enhancer"),
                pl.col("node_type").eq(pl.lit("Promoter")).cast(pl.UInt8).alias("is_gene_Promoter"),
                pl.col("node_type").eq(pl.lit("Promoter/Enhancer")).cast(pl.UInt8).alias("is_gene_Promoter_Enhancer")
            )
        )

        # --- 6. Select and Order Final Columns ---
        final_columns = [
            "node_id_1",
            "node_id_2",
            "is_gene_Enhancer",
            "is_gene_Promoter",
            "is_gene_Promoter_Enhancer",
            "is_elite",
            "combined_score"
        ]
        
        edges_final_df = edges_df.select(final_columns)

        print(f"Successfully created Edges table with {len(edges_final_df)} rows.")

        # --- 7. Export the Edges Table ---
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        print(f"Exporting Edges table to: {output_file}")
        edges_final_df.write_csv(output_file)
        
        print("Export complete. Script finished successfully.")

    except pl.exceptions.ShapeError as e:
        print(f"Error: Empty DataFrame. This might be due to no matches found during joins.")
        print(f"Details: {e}")
    except pl.exceptions.ComputeError as e:
        print(f"Error during data processing: {e}")
    except FileNotFoundError as e:
        print(f"Error: Input file not found.")
        print(f"Details: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    create_edges_file()