import pandas as pd
import pyBigWig
import os
import numpy as np
import concurrent.futures
from functools import partial
import multiprocessing

# --- Worker Function ---
# This MUST be a top-level function so it can be 'pickled' and sent
# to other processes by the multiprocessing module.

def process_chunk(chunk_df, bigwig_file_path):
    """
    Worker function to process a chunk of the DataFrame against a single bigWig file.
    This function is executed in a separate, parallel process.
    
    Args:
        chunk_df (pd.DataFrame): A piece of the original DataFrame.
        bigwig_file_path (str): The path to the .bigWig file to query.
        
    Returns:
        list: A list of mean values (or 0.0) for this chunk.
    """
    mean_values = []
    try:
        # Each worker process must open its own file handle for the bigWig file.
        # File handles cannot be shared across processes.
        bw = pyBigWig.open(bigwig_file_path)
        
        # Get valid chromosomes for this specific bigWig file
        valid_chroms = set(bw.chroms().keys())

        # Iterate over the rows of this specific chunk
        for row in chunk_df.itertuples(index=False):
            if row.chr in valid_chroms:
                # Ensure coordinates are int (they should be, but good to be safe)
                start, end = int(row.start), int(row.end)
                
                # Fetch the mean signal value
                stats = bw.stats(row.chr, start, end, type='mean')
                
                if stats is not None and stats[0] is not None:
                    mean_values.append(stats[0])
                else:
                    # No data for this interval
                    mean_values.append(0.0)
            else:
                # Chromosome not in this bigWig file
                mean_values.append(0.0)
        
        bw.close()
        return mean_values
        
    except Exception as e:
        print(f"[Worker Error] Error processing {bigwig_file_path}: {e}")
        # On failure, return a list of 0s of the correct length for this chunk
        return [0.0] * len(chunk_df)


def process_histone_data(input_csv_path, output_csv_path, histone_bigwig_map):
    """
    Processes a CSV file of genomic coordinates to fill in specified columns
    with mean signal values from corresponding bigWig files using pyBigWig.
    
    This version uses a ProcessPoolExecutor to parallelize the calculations 
    for each bigWig file across multiple CPU cores.
    """
    print(f"Starting processing for: {os.path.basename(input_csv_path)}")

    try:
        # 1. Import the CSV file into a pandas DataFrame
        df = pd.read_csv(input_csv_path)
        print(f"Successfully loaded {input_csv_path}. Shape: {df.shape}")

        # Ensure coordinate columns are of the correct type
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)

        # Determine the number of workers to use.
        # Using cpu_count() - 1 is often a good balance to leave one core
        # free for the main process and operating system tasks.
        num_workers = max(1, multiprocessing.cpu_count() - 1)
        print(f"Using {num_workers} worker processes.")

        # 2. For each entry in the map, consult the respective .bigWig file
        for column_name, bigwig_file in histone_bigwig_map.items():
            print(f"  Calculating mean signal for {column_name} (in parallel)...")

            if not os.path.exists(bigwig_file):
                print(f"    WARNING: BigWig file not found at {bigwig_file}. Skipping {column_name}.")
                if column_name not in df.columns:
                    df[column_name] = 0.0
                else:
                    print(f"    WARNING: Column {column_name} already exists. Setting to 0.0 as file not found.")
                    df[column_name] = 0.0
                continue

            try:
                # Split the DataFrame into 'num_workers' chunks
                # np.array_split is efficient for this.
                chunks = np.array_split(df, num_workers)

                # Use functools.partial to create a new function that has the
                # 'bigwig_file_path' argument fixed. This is needed because
                # executor.map takes a function and *one* iterable (our chunks).
                worker_func = partial(process_chunk, bigwig_file_path=bigwig_file)

                mean_values = []
                # Create the process pool
                with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
                    
                    # map applies the 'worker_func' to each 'chunk' in parallel
                    # It returns the results in the same order as the input chunks.
                    results_list = executor.map(worker_func, chunks)
                    
                    # The result is a list of lists (one list from each worker).
                    # We flatten it into a single list.
                    for sublist in results_list:
                        mean_values.extend(sublist)
                
                # Assign the collated results to the DataFrame column
                df[column_name] = mean_values
                print(f"    Done with {column_name}.")

            except Exception as e:
                print(f"    ERROR processing {column_name} with {bigwig_file}: {e}")
                print(f"    Setting column {column_name} to 0 as a fallback.")
                df[column_name] = 0.0

        # 3. Export the updated DataFrame to a new CSV file
        df.to_csv(output_csv_path, index=False)
        print(f"Processing complete. Updated data saved to: {output_csv_path}\n")

    except FileNotFoundError:
        print(f"ERROR: Input file not found at {input_csv_path}")
    except Exception as e:
        print(f"An unexpected error occurred while processing {input_csv_path}: {e}")


if __name__ == '__main__':
    # Define paths for the bigWig files
    base_path_chipseq = '../../public/Public Dataset/ChIPseq/'
    
    # This map now includes the 'PoI' column and its corresponding bigWig file
    histone_map = {
        # 'PoI': os.path.join('../../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2_ppois.bigWig'),
        'PoI': os.path.join('../../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2_ppois.bigWig'),
        'H3K4me1': os.path.join(base_path_chipseq, 'MCF7 H3K4me1 ChIPseq GSE86714 ENCFF763NCP hg38.bigWig'),
        'H3K4me2': os.path.join(base_path_chipseq, 'MCF7 H3K4me2 ChIPseq GSE96439 ENCFF442RRY hg38.bigWig'),
        'H3K4me3': os.path.join(base_path_chipseq, 'MCF7 H3K4me3 ChIPseq GSE96506 ENCFF163MXP hg38.bigWig'),
        'H3K9ac': os.path.join(base_path_chipseq, 'MCF7 H3K9ac ChIPseq GSE95898 ENCFF327XJC hg38.bigWig'),
        'H3K9me3': os.path.join(base_path_chipseq, 'MCF7 H3K9me3 ChIPseq GSE96517 ENCFF481DZL hg38.bigWig'),
        'H3K27ac': os.path.join(base_path_chipseq, 'MCF7 H3K27ac ChIPseq GSE96352 ENCFF138YNG hg38.bigWig'),
        'H3K27me3': os.path.join(base_path_chipseq, 'MCF7 H3K27me3 ChIPseq GSE96363 ENCFF163QKN hg38.bigWig'),
        'H3K36me3': os.path.join(base_path_chipseq, 'MCF7 H3K36me3 ChIPseq GSE174945 ENCFF910BRP hg38.bigWig'),
        'H4K20me1': os.path.join(base_path_chipseq, 'MCF7 H4K20me1 ChIPseq GSE96283 ENCFF366GLZ hg38.bigWig')
    }

    # Define paths for the input and output CSV files
    base_path_multiomics = '../../public/Multiomics/'
    positive_input_path = os.path.join(base_path_multiomics, 'hg38_bin_positive.csv')
    negative_input_path = os.path.join(base_path_multiomics, 'hg38_bin_negative.csv')
    positive_output_path = os.path.join(base_path_multiomics, 'Step18_ML_Prediction_Regression/hg38_bin_positive_histone.csv')
    negative_output_path = os.path.join(base_path_multiomics, 'Step18_ML_Prediction_Regression/hg38_bin_negative_histone.csv')

    # Process the positive dataset
    process_histone_data(positive_input_path, positive_output_path, histone_map)

    # Process the negative dataset
    process_histone_data(negative_input_path, negative_output_path, histone_map)

    print("All tasks completed.")
