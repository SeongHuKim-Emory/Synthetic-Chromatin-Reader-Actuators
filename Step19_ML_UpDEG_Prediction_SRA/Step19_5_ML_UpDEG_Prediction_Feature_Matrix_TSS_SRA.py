import pandas as pd
import pyBigWig
import numpy as np
import os


def process_genomic_data(input_csv_path, bigwig_reader, output_csv_path):
    """
    Processes a genomic data CSV file by sorting, calculating distances,
    and annotating with max signal values from a bigWig file.

    Args:
        input_csv_path (str): Path to the input CSV file.
        bigwig_reader (pyBigWig): An open pyBigWig file object.
        output_csv_path (str): Path to save the processed output CSV file.
    """
    print(f"--- Processing {os.path.basename(input_csv_path)} ---")

    # Steps 2-5: Read the CSV file into a pandas DataFrame.
    # This effectively creates a new table in memory which is a copy of the original.
    try:
        df = pd.read_csv(input_csv_path)
        print(f"Successfully loaded {input_csv_path}. Shape: {df.shape}")
    except FileNotFoundError:
        print(f"Error: The file {input_csv_path} was not found.")
        return

    # Step 6: Sort rows by the "symbol" column values.
    # We use ignore_index=True to reset the index after sorting.
    df.sort_values(by='symbol', inplace=True, ignore_index=True)
    print("Sorted DataFrame by the 'symbol' column.")

    # Step 7: Add "025_on_element" and "Element_TSS_Distance" columns.
    # They are initialized with NaN (Not a Number), a standard placeholder for missing numeric data.
    df['025_on_element'] = np.nan
    df['Element_TSS_Distance'] = np.nan
    print("Added '025_on_element' and 'Element_TSS_Distance' columns.")

    # Step 8: Calculate "Element_TSS_Distance".
    # This vectorized operation is much faster than looping through rows.
    # It calculates the absolute difference between the TSS and both the element start and end,
    # and then takes the maximum of those two differences for each row.
    df['Element_TSS_Distance'] = np.maximum(
        (df['element_start'] - df['TSS']).abs(),
        (df['element_end'] - df['TSS']).abs()
    )
    print("Calculated 'Element_TSS_Distance' values.")

    # Step 9: Get the maximum value from the bigWig file for each element's genomic range.
    def get_max_from_bigwig(row):
        """Helper function to retrieve max signal for a given row's coordinates."""
        try:
            # Ensure coordinates are standard Python integers, as pyBigWig can be sensitive to numpy types.
            chrom = row['chr']
            start = int(row['element_start'])
            end = int(row['element_end'])

            # Query the bigWig file for the max value in the specified range.
            # The stats method returns a list containing the single max value, e.g., [0.123].
            max_val = bigwig_reader.stats(chrom, start, end, type='max')

            # If a value is found, extract it from the list. Otherwise, return NaN.
            if max_val and max_val[0] is not None:
                return max_val[0]
            else:
                return np.nan
        except (RuntimeError, ValueError) as e:
            # This handles cases where the chromosome name (e.g., 'chrM') from the CSV
            # might not be present in the bigWig file, preventing the script from crashing.
            # print(f"Warning: Could not process row for {chrom}:{start}-{end}. Error: {e}")
            return np.nan

    print("Annotating with max values from bigWig file. This may take a moment...")
    df['025_on_element'] = df.apply(get_max_from_bigwig, axis=1)
    print("Annotation complete.")

    # Steps 10-11: Export the processed DataFrame to a new CSV file.
    # index=False prevents pandas from writing the DataFrame index as a column.
    output_dir = os.path.dirname(output_csv_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    df.to_csv(output_csv_path, index=False)
    print(f"Successfully saved processed data to {output_csv_path}")


if __name__ == "__main__":
    # --- Configuration: Define all file paths ---

    # Input file paths
    up_deg_input_csv = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS.csv'
    non_deg_input_csv = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS.csv'
    bigwig_file_path = '../../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2_ppois.bigWig'

    # Output file paths
    up_deg_output_csv = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA.csv'
    non_deg_output_csv = '../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA/GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA.csv'

    # --- Main Execution ---

    bw = None  # Initialize bigWig reader variable
    try:
        # Step 1: Import the bigWig file
        print(f"Opening bigWig file: {bigwig_file_path}")
        bw = pyBigWig.open(bigwig_file_path)

        # Process the UpDEG file
        process_genomic_data(
            input_csv_path=up_deg_input_csv,
            bigwig_reader=bw,
            output_csv_path=up_deg_output_csv
        )

        print("\n" + "=" * 50 + "\n")

        # Process the Non-DEG file
        process_genomic_data(
            input_csv_path=non_deg_input_csv,
            bigwig_reader=bw,
            output_csv_path=non_deg_output_csv
        )

    except FileNotFoundError:
        print(f"Error: The bigWig file was not found at {bigwig_file_path}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    finally:
        # Ensure the bigWig file is closed properly
        if bw:
            bw.close()
            print("\nbigWig file closed.")
        print("\nScript execution finished.")
