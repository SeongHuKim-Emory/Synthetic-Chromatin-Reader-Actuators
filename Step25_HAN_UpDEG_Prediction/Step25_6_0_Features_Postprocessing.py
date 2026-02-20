import polars as pl
import sys

def process_file(file_path: str):
    """
    Reads a CSV file, deduplicates based on 'node_id' by keeping the row
    with the maximum 'TPM' value, and overwrites the original file.
    """
    try:
        # 1. Read the CSV file
        df = pl.read_csv(file_path)
        
        # 2. Identify rows with identical "node_id" values and keep the one with max "TPM"
        # We group by 'node_id'. Within each group, we sort all columns
        # by 'TPM' in descending order (with nulls last) and take the 'first' row.
        # This effectively selects the row with the max 'TPM' for each 'node_id'.
        print(f"  Original shape of {file_path}: {df.shape}")
        
        df_deduplicated = (
            df.group_by("node_id")
            .agg(pl.all().sort_by("TPM", descending=True, nulls_last=True).first())
            .sort("node_id") # Sort by node_id for a clean, ordered output file
        )
        
        # 3. Save the updated DataFrame back to the original file path
        print(f"  New shape of {file_path}: {df_deduplicated.shape}")
        
        df_deduplicated.write_csv(file_path)
        
        print(f"Successfully processed and overwrote {file_path}")

    except pl.exceptions.ComputeError as e:
        print(f"Error processing file {file_path}: {e}", file=sys.stderr)
        print("Please check if the 'node_id' and 'TPM' columns exist.", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}", file=sys.stderr)
        print("Please check the file path.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred with {file_path}: {e}", file=sys.stderr)

def main():
    # Define file paths
    node_features_path = "../../public/Multiomics/Step25_HAN_UpDEG_Prediction/Node_Features.csv"
    features_path = "../../public/Multiomics/Step25_HAN_UpDEG_Prediction/Features.csv"
    
    print("Starting deduplication process...")
    
    # Process Node_Features.csv
    print(f"\nProcessing {node_features_path}...")
    process_file(node_features_path)
    
    # Process Features.csv
    print(f"\nProcessing {features_path}...")
    process_file(features_path)
    
    print("\nDeduplication process finished.")

if __name__ == "__main__":
    main()