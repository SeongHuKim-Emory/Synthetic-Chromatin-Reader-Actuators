import pandas as pd
import os


def process_gene_data(file_path):
    """
    Reads a gene data CSV, filters out rows with empty TSS values,
    formats specific columns to integers, and saves the updated data
    back to the original file.

    Args:
        file_path (str): The path to the input and output CSV file.
    """
    try:
        # Step 1: Import the CSV file into a pandas DataFrame
        print(f"Reading data from {file_path}...")
        df = pd.read_csv(file_path)
        print("Initial data loaded successfully.")
        print(f"Original number of rows: {len(df)}")

        # Step 2: For each row, if "TSS" column value is empty, delete the entire row
        # The .dropna() method with a subset handles this efficiently.
        df.dropna(subset=['TSS'], inplace=True)
        print("Rows with empty 'TSS' values have been removed.")
        print(f"Number of rows after filtering: {len(df)}")

        # Step 3: For gene_start, gene_end, TSS columns, convert them to integers
        # This removes the decimal point and any numbers to the right of it.
        columns_to_convert = ['gene_start', 'gene_end', 'TSS']
        for col in columns_to_convert:
            # Ensure the column is numeric before converting to integer to avoid errors
            df[col] = pd.to_numeric(df[col], errors='coerce').astype('Int64')

        print("Formatted 'gene_start', 'gene_end', and 'TSS' columns to integers.")

        # Step 4: Export the updated table back to the original file path
        df.to_csv(file_path, index=False)
        print(f"Successfully exported the updated data to {file_path}")

    except FileNotFoundError:
        print(f"Error: The file at {file_path} was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    # Define the file path. Using os.path.join for better cross-platform compatibility.
    # The path is relative to the script's location as per the request.
    file_path = os.path.join('../../', 'public', 'Multiomics', 'Step19_ML_UpDEG_Prediction_SRA',
                             'GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS.csv')

    # Create a dummy file for demonstration purposes if it doesn't exist.
    if not os.path.exists(file_path):
        print("Creating a dummy CSV file for demonstration...")
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        dummy_data = {
            'symbol': ['OR4F16', 'ENSG00000278791', 'MIR12136', 'MTATP6P1', 'MTATP8P1', 'MIR12136'],
            'chr': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
            'element_start': [629101, 631681, 631681, 631681, 631681, 633541],
            'element_end': [629896, 632470, 632470, 632470, 632470, 633690],
            'gene_chr': ['chr1', 'chr1', '', 'chr1', 'chr1', ''],
            'gene_start': [685716.0, 632325.0, None, 633696.0, 633535.0, None],
            'gene_end': [686654.0, 632413.0, None, 634376.0, 633741.0, None],
            'gene_strand': ['reverse', 'reverse', '', 'forward', 'forward', ''],
            'TSS': [686654.0, 632413.0, None, 633696.0, 633535.0, None]
        }
        pd.DataFrame(dummy_data).to_csv(file_path, index=False)

    process_gene_data(file_path)
