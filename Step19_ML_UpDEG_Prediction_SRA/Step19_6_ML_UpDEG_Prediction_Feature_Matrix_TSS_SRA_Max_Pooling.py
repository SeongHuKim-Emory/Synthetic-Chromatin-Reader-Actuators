import pandas as pd
import os


def filter_max_values_per_symbol(input_filepath, output_filepath):
    """
    Reads a CSV file, filters it to keep only the row with the maximum
    '025_on_element' value for each unique 'symbol'. In case of a tie,
    it keeps the row with the highest 'Element_TSS_Distance'.

    Args:
        input_filepath (str): The path to the input CSV file.
        output_filepath (str): The path where the processed CSV file will be saved.
    """
    try:
        # 1. Import the CSV file into a pandas DataFrame
        df = pd.read_csv(input_filepath)

        # 2. Sort the DataFrame.
        #    - First by 'symbol' (ascending).
        #    - Then by '025_on_element' (descending) to get the max value on top.
        #    - Finally, by 'Element_TSS_Distance' (descending) to handle ties in '025_on_element'.
        df_sorted = df.sort_values(
            by=['symbol', '025_on_element', 'Element_TSS_Distance'],
            ascending=[True, False, False]
        )

        # 3. For each unique symbol, keep only the first row.
        #    Due to the sorting, this row corresponds to the maximum desired values.
        df_max_pooled = df_sorted.drop_duplicates(subset='symbol', keep='first')

        # 4. Create the directory for the output file if it doesn't exist
        output_dir = os.path.dirname(output_filepath)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # 5. Export the updated DataFrame to a new CSV file
        df_max_pooled.to_csv(output_filepath, index=False)

        print(f"Successfully processed '{input_filepath}'")
        print(f"Result saved to '{output_filepath}'")

    except FileNotFoundError:
        print(f"Error: The file '{input_filepath}' was not found.")
    except Exception as e:
        print(f"An error occurred while processing {input_filepath}: {e}")


# --- Main execution ---
if __name__ == "__main__":
    # Define file paths for the UpDEG dataset
    updeg_input_path = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA.csv'
    updeg_output_path = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool.csv'

    # Define file paths for the Non-DEG dataset
    non_deg_input_path = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA.csv'
    non_deg_output_path = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool.csv'

    # Process both files
    print("Starting file processing...")
    filter_max_values_per_symbol(updeg_input_path, updeg_output_path)
    print("-" * 30)
    filter_max_values_per_symbol(non_deg_input_path, non_deg_output_path)
    print("\nProcessing complete.")
