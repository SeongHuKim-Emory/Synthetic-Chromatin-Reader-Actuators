import pandas as pd
import os


def process_predictions_to_bed(input_path, output_path):
    """
    Reads a CSV file, filters and processes the data, and saves it as a
    tab-delimited BED file without a header.

    Args:
        input_path (str): The path to the input CSV file.
        output_path (str): The path where the output BED file will be saved.
    """
    try:
        # Step 1: Import the CSV file into a pandas DataFrame.
        # Assuming the file has a header row as shown in the example.
        print(f"Reading data from {input_path}...")
        df = pd.read_csv(input_path)

        # Step 2: Remove rows where the "predicted_label" column values are 0.
        print("Filtering rows where 'predicted_label' is not 1...")
        # We keep rows where the predicted_label is 1.
        df_filtered = df[df['predicted_label'] != 0].copy()

        # Step 3: Remove the 'true_label', 'predicted_label', and
        # 'predicted_probability' columns.
        print("Selecting 'chr', 'start', and 'end' columns...")
        df_final = df_filtered[['chr', 'start', 'end']]

        # Ensure the output directory exists.
        output_dir = os.path.dirname(output_path)
        if not os.path.exists(output_dir) and output_dir:
            os.makedirs(output_dir)
            print(f"Created directory: {output_dir}")

        # Step 4, 5, and 6: Remove the title row, change to tab-delimited,
        # and export to the specified file.
        # The to_csv method handles these requirements with its parameters:
        # - sep='\t': Sets the delimiter to a tab.
        # - header=False: Excludes the column headers from the output file.
        # - index=False: Prevents writing the DataFrame index.
        print(f"Exporting the processed data to {output_path}...")
        df_final.to_csv(output_path, sep='\t', header=False, index=False)

        print("\nProcessing complete.")
        print(f"The output has been successfully saved to: {output_path}")

    except FileNotFoundError:
        print(f"Error: The input file was not found at {input_path}")
    except KeyError as e:
        print(f"Error: The CSV file is missing a required column: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


# --- Main execution block ---
if __name__ == "__main__":
    # Define the input and output file paths.
    # The script assumes it is run from a directory where this relative path is valid.
    input_file_path = '../../public/multiomics/test_predictions.csv'
    output_file_path = '../../public/multiomics/Predicted_025_binding.bed'

    # Run the processing function
    process_predictions_to_bed(input_file_path, output_file_path)

