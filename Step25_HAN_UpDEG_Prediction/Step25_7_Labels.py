import polars as pl
import os

# Define file paths
# Use os.path.join for better cross-platform compatibility
# Assumes the script is run from a directory where this relative path is valid
base_path = os.path.join('..', '..', 'public', 'Multiomics', 'Step25_HAN_UpDEG_Prediction')
input_file = os.path.join(base_path, 'Gene_Nodes.csv')
output_file = os.path.join(base_path, 'Labels.csv')

try:
    # 1. Import files
    print(f"Attempting to read file from: {input_file}")
    # Explicitly state the Python version this script is intended for (as requested)
    # Note: polars functionality used here is compatible with Python 3.12.3
    
    gene_nodes_df = pl.read_csv(input_file)
    print("Successfully read Gene_Nodes.csv.")
    print("\nOriginal DataFrame (head):")
    print(gene_nodes_df.head())

    # 2. Create a new table named Labels and copy entire content
    # In polars, this means creating a new DataFrame using .clone()
    labels_df = gene_nodes_df.clone()
    print("\nSuccessfully cloned data to 'Labels' DataFrame.")

    # 3. Remove the symbol, chr, TSS in the Labels table
    columns_to_remove = ['symbol', 'chr', 'TSS']
    
    # Check which of the columns to remove actually exist in the DataFrame
    existing_columns_to_remove = [col for col in columns_to_remove if col in labels_df.columns]
    
    if existing_columns_to_remove:
        # .drop() in polars returns a new DataFrame
        labels_df = labels_df.drop(existing_columns_to_remove)
        print(f"Successfully removed columns: {', '.join(existing_columns_to_remove)}")
    else:
        print(f"Columns to remove ({', '.join(columns_to_remove)}) not found in the DataFrame.")

    print("\nModified 'Labels' DataFrame (head):")
    print(labels_df.head())

    # 4. Export the Labels table
    
    # Ensure the output directory exists
    os.makedirs(base_path, exist_ok=True)
    
    # Export the DataFrame
    # Polars' write_csv does not include an index by default
    labels_df.write_csv(output_file)
    print(f"\nSuccessfully exported 'Labels' table to: {output_file}")

except FileNotFoundError:
    print(f"Error: The file was not found at {input_file}")
    print("Please check the relative path and ensure the file exists.")
except Exception as e:
    print(f"An error occurred: {e}")

