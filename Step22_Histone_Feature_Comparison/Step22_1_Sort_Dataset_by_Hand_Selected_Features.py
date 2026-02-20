import polars as pl
import os

def process_and_filter_data_polars(input_path: str, output_path: str):
    """
    Reads a CSV file using polars, filters its columns based on specified criteria,
    and saves the result to a new CSV file.

    Args:
        input_path (str): The path to the input CSV file.
        output_path (str): The path to save the filtered CSV file.
    """
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")

    # Read the CSV file with polars
    try:
        df = pl.read_csv(input_path)
        print(f"Successfully loaded {input_path} with polars.")
    except Exception as e:
        print(f"Error loading file {input_path}: {e}")
        return

    # --- Column Filtering Logic ---
    
    # Define the suffixes for columns to be retained
    suffixes_to_keep = [
        '_025', '_H3K4me1', '_H3K4me2', '_H3K4me3', '_H3K9ac',
        '_H3K9me3', '_H3K27ac', '_H3K27me3', '_H3K36me3', '_H4K20me1'
    ]
    
    original_column_count = len(df.columns)
    
    # Build the list of columns to keep
    columns_to_keep = ['Symbol']  # Always include the 'Symbol' column

    for col_name in df.columns:
        if col_name == 'Symbol':
            continue
        
        # Check if the column name ends with any of the desired suffixes
        # and does not contain '_Distal_'
        if any(col_name.endswith(suffix) for suffix in suffixes_to_keep) and '_Distal_' not in col_name:
            columns_to_keep.append(col_name)

    # Use the select method to create the new filtered DataFrame
    df_filtered = df.select(columns_to_keep)
    
    print(f"Filtered columns. Kept {len(df_filtered.columns)} out of {original_column_count} columns.")

    # Write the filtered DataFrame to a new CSV file
    try:
        df_filtered.write_csv(output_path)
        print(f"Filtered data successfully saved to {output_path}")
    except Exception as e:
        print(f"Error saving file to {output_path}: {e}")

# --- Main execution ---

# Define relative paths. Adjust if your script's location is different.
base_input_dir = '../../public/Multiomics/Step19_ML_UpDEG_Prediction/'
base_output_dir = '../../public/Multiomics/Step22_Histone_Feature_Comparison/'

# File 1: Process the Up-regulated DEGs file
updeg_input_file = os.path.join(base_input_dir, 'SRA_Specific_UpDEG_with_GO_CCLE.csv')
updeg_output_file = os.path.join(base_output_dir, 'SRA_Specific_UpDEG_with_GO_CCLE_Select_Features.csv')
process_and_filter_data_polars(updeg_input_file, updeg_output_file)

print("-" * 50)

# File 2: Process the Non-DEGs file
non_deg_input_file = os.path.join(base_input_dir, 'Non_DEG_with_GO_CCLE.csv')
non_deg_output_file = os.path.join(base_output_dir, 'Non_DEG_with_GO_CCLE_Select_Features.csv')
process_and_filter_data_polars(non_deg_input_file, non_deg_output_file)

print("\nProcessing complete.")
