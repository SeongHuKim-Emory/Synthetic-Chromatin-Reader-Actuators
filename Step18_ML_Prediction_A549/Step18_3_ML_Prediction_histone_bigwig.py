import pandas as pd
import pyBigWig
import os


def process_histone_data(input_csv_path, output_csv_path, histone_bigwig_map):
    """
    Processes a CSV file of genomic coordinates to fill in histone mark columns
    with mean signal values from corresponding bigWig files using pyBigWig.

    Args:
        input_csv_path (str): The file path for the input CSV.
        output_csv_path (str): The file path for the output CSV.
        histone_bigwig_map (dict): A dictionary mapping histone mark names to their bigWig file paths.
    """
    print(f"Starting processing for: {os.path.basename(input_csv_path)}")

    try:
        # 1. Import the CSV file into a pandas DataFrame
        df = pd.read_csv(input_csv_path)
        print(f"Successfully loaded {input_csv_path}. Shape: {df.shape}")

        # Ensure coordinate columns are of the correct type
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)

        # 2. For each histone mark, consult the respective .bigWig file
        for histone_mark, bigwig_file in histone_bigwig_map.items():
            print(f"  Calculating mean signal for {histone_mark}...")

            if not os.path.exists(bigwig_file):
                print(f"    WARNING: BigWig file not found at {bigwig_file}. Skipping {histone_mark}.")
                df[histone_mark] = 0.0
                continue

            try:
                # Open the bigWig file for reading
                bw = pyBigWig.open(bigwig_file)

                # Get the set of valid chromosomes from the bigWig file for quick lookups
                valid_chroms = set(bw.chroms().keys())

                mean_values = []
                # Use itertuples for more efficient row iteration
                for row in df.itertuples(index=False):
                    # Check if the chromosome from the dataframe exists in the bigWig file
                    if row.chr in valid_chroms:
                        # Fetch the mean signal value. It can return None if no data exists for the interval.
                        stats = bw.stats(row.chr, row.start, row.end, type='mean')

                        # Append the value, or 0.0 if it's None
                        if stats is not None and stats[0] is not None:
                            mean_values.append(stats[0])
                        else:
                            mean_values.append(0.0)
                    else:
                        # If chromosome is not present in the bigWig file, append 0.0
                        mean_values.append(0.0)

                # Assign the calculated mean values to the new column
                df[histone_mark] = mean_values
                bw.close()
                print(f"    Done with {histone_mark}.")

            except Exception as e:
                print(f"    ERROR processing {histone_mark} with {bigwig_file}: {e}")
                print(f"    Setting column {histone_mark} to 0 as a fallback.")
                df[histone_mark] = 0.0

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
    base_path_multiomics = '../../public/Multiomics/'
    positive_input_path = os.path.join(base_path_multiomics, 'hg38_bin_positive.csv')
    negative_input_path = os.path.join(base_path_multiomics, 'hg38_bin_negative.csv')
    positive_output_path = os.path.join(base_path_multiomics, 'hg38_bin_positive_histone.csv')
    negative_output_path = os.path.join(base_path_multiomics, 'hg38_bin_negative_histone.csv')

    # Process the positive dataset
    process_histone_data(positive_input_path, positive_output_path, histone_map)

    # Process the negative dataset
    process_histone_data(negative_input_path, negative_output_path, histone_map)

    print("All tasks completed.")
