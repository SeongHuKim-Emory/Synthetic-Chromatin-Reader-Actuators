import pandas as pd
import pyBigWig
import os
import multiprocessing
import numpy as np

def process_dataframe_chunk(args):
    """
    Worker function to process a chunk of the DataFrame.
    This function contains the logic previously inside the main loop.
    
    Args:
        args (tuple): A tuple containing (df_chunk, histone_bigwig_map).
    
    Returns:
        pd.DataFrame: The processed chunk with new histone columns.
    """
    df_chunk, histone_bigwig_map = args
    
    # We work on a copy to avoid SettingWithCopy warnings in parallel contexts
    df_chunk = df_chunk.copy()

    # Iterate through each histone mark map provided
    for histone_mark, bigwig_file in histone_bigwig_map.items():
        if not os.path.exists(bigwig_file):
            # Silent print in workers to avoid console clutter, or keep if debugging needed
            # print(f"    WARNING: BigWig file not found at {bigwig_file}. Skipping {histone_mark}.")
            df_chunk[histone_mark] = 0.0
            continue

        try:
            # Open the bigWig file for reading (Must be done inside the worker)
            bw = pyBigWig.open(bigwig_file)

            # Get the set of valid chromosomes
            valid_chroms = set(bw.chroms().keys())

            mean_values = []
            
            # Use itertuples for efficient row iteration
            for row in df_chunk.itertuples(index=False):
                # Check if the chromosome exists in the bigWig file
                if row.chr in valid_chroms:
                    # Fetch the mean signal value
                    stats = bw.stats(row.chr, row.start, row.end, type='mean')

                    if stats is not None and stats[0] is not None:
                        mean_values.append(stats[0])
                    else:
                        mean_values.append(0.0)
                else:
                    mean_values.append(0.0)

            # Assign values to the chunk
            df_chunk[histone_mark] = mean_values
            bw.close()

        except Exception as e:
            # Fallback logic for errors
            print(f"    ERROR processing {histone_mark} in chunk: {e}")
            df_chunk[histone_mark] = 0.0
            
    return df_chunk

def process_histone_data(input_csv_path, output_csv_path, histone_bigwig_map, n_cores=25):
    """
    Processes a CSV file of genomic coordinates to fill in histone mark columns
    using parallel processing.

    Args:
        input_csv_path (str): The file path for the input CSV.
        output_csv_path (str): The file path for the output CSV.
        histone_bigwig_map (dict): A dictionary mapping histone mark names to their bigWig file paths.
        n_cores (int): Number of CPUs to use for parallel processing.
    """
    print(f"Starting parallel processing ({n_cores} cores) for: {os.path.basename(input_csv_path)}")

    try:
        # 1. Import the CSV file into a pandas DataFrame
        df = pd.read_csv(input_csv_path)
        print(f"Successfully loaded {input_csv_path}. Shape: {df.shape}")

        # Ensure coordinate columns are of the correct type
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)

        # 2. Split the DataFrame into chunks for parallel processing
        # If the dataframe is smaller than n_cores, split into len(df) chunks
        num_chunks = min(n_cores, len(df))
        df_chunks = np.array_split(df, num_chunks)

        # Prepare arguments for the worker function: (chunk, map)
        # We need to pass the map to every worker
        worker_args = [(chunk, histone_bigwig_map) for chunk in df_chunks]

        print(f"Distributing work across {num_chunks} chunks...")

        # 3. Initialize Multiprocessing Pool
        with multiprocessing.Pool(processes=n_cores) as pool:
            # Map the worker function to the chunks
            processed_chunks = pool.map(process_dataframe_chunk, worker_args)

        # 4. Concatenate the processed chunks back into a single DataFrame
        print("Parallel processing complete. Concatenating results...")
        result_df = pd.concat(processed_chunks, ignore_index=True)

        # 5. Export the updated DataFrame to a new CSV file
        result_df.to_csv(output_csv_path, index=False)
        print(f"Processing complete. Updated data saved to: {output_csv_path}\n")

    except FileNotFoundError:
        print(f"ERROR: Input file not found at {input_csv_path}")
    except Exception as e:
        print(f"An unexpected error occurred while processing {input_csv_path}: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    # Define paths for the bigWig files
    base_path_chipseq = '../../public/Public Dataset/ChIPseq/'
    histone_map = {
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
    base_path_multiomics = '../../public/Multiomics/Step18_ML_Prediction_SRA/'
    positive_input_path = os.path.join(base_path_multiomics, 'hg38_bin_positive.csv')
    negative_input_path = os.path.join(base_path_multiomics, 'hg38_bin_negative.csv')
    positive_output_path = os.path.join(base_path_multiomics, 'hg38_bin_positive_histone.csv')
    negative_output_path = os.path.join(base_path_multiomics, 'hg38_bin_negative_histone.csv')

    # CPU Count (Hardcoded to 25 as requested)
    CPUS_TO_USE = 25

    # Process the positive dataset
    process_histone_data(positive_input_path, positive_output_path, histone_map, n_cores=CPUS_TO_USE)

    # Process the negative dataset
    process_histone_data(negative_input_path, negative_output_path, histone_map, n_cores=CPUS_TO_USE)

    print("All tasks completed.")