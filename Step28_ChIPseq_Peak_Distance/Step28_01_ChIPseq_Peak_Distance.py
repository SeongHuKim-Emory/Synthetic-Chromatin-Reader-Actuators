import pandas as pd
import numpy as np
import os

def main():
    # ---------------------------------------------------------
    # 1. Define File Paths
    # ---------------------------------------------------------
    file_path_021 = "../../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2 Filtered_trim.bed"
    file_path_025 = "../../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2 Filtered_trim.bed"
    output_path = "../../public/Multiomics/Step28_ChIPseq_Peak_Distance/Peak_Distance.csv"

    # Check if input files exist to avoid runtime errors
    if not os.path.exists(file_path_021):
        print(f"Error: File not found at {file_path_021}")
        return
    if not os.path.exists(file_path_025):
        print(f"Error: File not found at {file_path_025}")
        return

    # ---------------------------------------------------------
    # 2. Import Files
    # ---------------------------------------------------------
    print("Loading files...")
    # Read files, assuming no header as per the prompt structure
    # Assigning column names for internal processing
    col_names = ['chr', 'start', 'end']
    
    df_021 = pd.read_csv(file_path_021, sep='\t', header=None, names=col_names)
    df_025 = pd.read_csv(file_path_025, sep='\t', header=None, names=col_names)

    # ---------------------------------------------------------
    # 3. Prepare Peak_Distance Table
    # ---------------------------------------------------------
    print("Initializing Peak_Distance table...")
    
    # Copy content from 021 file
    peak_distance = df_021.copy()
    
    # Rename columns as requested
    peak_distance.columns = ['021_chr', '021_start', '021_end']
    
    # Initialize new columns with empty values (NaN)
    peak_distance['025_chr'] = ""
    peak_distance['025_start'] = np.nan
    peak_distance['025_end'] = np.nan
    peak_distance['Distance'] = np.nan

    # ---------------------------------------------------------
    # 4. Find Closest Peaks
    # ---------------------------------------------------------
    print("Calculating distances (this may take a moment)...")

    # Optimization: Group df_025 by chromosome to avoid filtering the whole dataframe in every loop
    df_025_grouped = dict(tuple(df_025.groupby('chr')))

    # Lists to collect results to assign at once (faster than row-by-row DataFrame update)
    res_025_chr = []
    res_025_start = []
    res_025_end = []
    res_distance = []

    for index, row in peak_distance.iterrows():
        chrom = row['021_chr']
        start = row['021_start']
        end = row['021_end']

        # Get subset of 025 peaks on the same chromosome
        if chrom in df_025_grouped:
            targets = df_025_grouped[chrom]
            
            # Vectorized distance calculation
            # Distance logic: max(0, start2 - end1, start1 - end2)
            # If overlap, distance is 0.
            
            # Calculate distance from current peak (row) to all targets
            # targets['start'] - end  -> distance if target is to the right
            # start - targets['end']  -> distance if target is to the left
            
            dist_right = targets['start'] - end
            dist_left = start - targets['end']
            
            # Element-wise max to handle overlap logic
            # The calculation implies: if intervals overlap, one value is negative.
            # If fully overlapping or contained, both are negative.
            # We want the max of (0, dist_right, dist_left).
            
            # Using numpy for speed
            dists = np.maximum(0, np.maximum(dist_right.values, dist_left.values))
            
            # Find the index of the minimum distance
            min_idx_loc = np.argmin(dists)
            min_dist = dists[min_idx_loc]
            
            # Get the actual row from targets
            closest_peak = targets.iloc[min_idx_loc]
            
            res_025_chr.append(closest_peak['chr'])
            res_025_start.append(closest_peak['start'])
            res_025_end.append(closest_peak['end'])
            res_distance.append(min_dist)
            
        else:
            # No peaks found on this chromosome in file 025
            res_025_chr.append(None)
            res_025_start.append(None)
            res_025_end.append(None)
            res_distance.append(None)

    # Assign results back to DataFrame
    peak_distance['025_chr'] = res_025_chr
    peak_distance['025_start'] = res_025_start
    peak_distance['025_end'] = res_025_end
    peak_distance['Distance'] = res_distance

    # Cleaning: Fill NaN in 025_chr if any (converts None to empty string or keeps NaN)
    # Ensure start/end/distance are integers where possible, but keep NaN if missing
    
    # ---------------------------------------------------------
    # 5. Export to CSV
    # ---------------------------------------------------------
    print(f"Exporting to {output_path}...")
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Save to CSV
    peak_distance.to_csv(output_path, index=False)
    print("Done.")

if __name__ == "__main__":
    main()