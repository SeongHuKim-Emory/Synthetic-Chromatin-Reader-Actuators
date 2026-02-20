import pandas as pd
import os

########################################################################################################################

########################################################################################################################

# --- Configuration ---
tad_file = '../public/Public Dataset/HiC/GSE66733_MCF7_HiC/MCF7_hg38_TADs.csv'
rnaseq_file = '../public/RNAseq/UpDEG_category.csv'
output_dir = '../public/Multiomics/'
output_file = os.path.join(output_dir, 'MCF7_SRA_HiC_RNAseq.csv')

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# --- 1. Import the files ---
print(f"Loading TAD file: {tad_file}")
try:
    tads_df = pd.read_csv(tad_file)
    print(f"TADs loaded successfully. Shape: {tads_df.shape}")
    # Ensure coordinate columns are numeric
    tads_df['hg38_TAD_Start'] = pd.to_numeric(tads_df['hg38_TAD_Start'])
    tads_df['hg38_TAD_End'] = pd.to_numeric(tads_df['hg38_TAD_End'])
except FileNotFoundError:
    print(f"Error: TAD file not found at {tad_file}")
    exit()
except Exception as e:
    print(f"Error loading TAD file: {e}")
    exit()


print(f"Loading RNAseq file: {rnaseq_file}")
try:
    rnaseq_df = pd.read_csv(rnaseq_file)
    print(f"RNAseq data loaded successfully. Shape: {rnaseq_df.shape}")
    # Ensure coordinate columns are numeric and handle potential non-numeric values
    rnaseq_df['hg38_start'] = pd.to_numeric(rnaseq_df['hg38_start'], errors='coerce')
    rnaseq_df['hg38_end'] = pd.to_numeric(rnaseq_df['hg38_end'], errors='coerce')
    # Drop rows where coordinates couldn't be converted (or handle differently if needed)
    rnaseq_df.dropna(subset=['hg38_start', 'hg38_end'], inplace=True)
    rnaseq_df['hg38_start'] = rnaseq_df['hg38_start'].astype(int)
    rnaseq_df['hg38_end'] = rnaseq_df['hg38_end'].astype(int)

    # Add 'chr' prefix to chromosome column for matching
    rnaseq_df['chr_match'] = 'chr' + rnaseq_df['hg38_chr'].astype(str)

except FileNotFoundError:
    print(f"Error: RNAseq file not found at {rnaseq_file}")
    exit()
except Exception as e:
    print(f"Error loading RNAseq file: {e}")
    exit()

# --- 2. Merge the information ---

results = []
# Keep track of which TAD indices have been matched with at least one transcript
matched_tad_indices = set()
# Get column names for creating placeholder rows later
rnaseq_cols = rnaseq_df.columns.tolist()
# Create a dictionary with NA values for transcript columns
rna_na_dict = {col: pd.NA for col in rnaseq_cols}

print("Merging data based on genomic coordinates...")
# Iterate through each transcript (gene)
for gene_idx, gene_row in rnaseq_df.iterrows():
    gene_chr = gene_row['chr_match']
    gene_start = gene_row['hg38_start']
    gene_end = gene_row['hg38_end']
    gene_data = gene_row.to_dict() # Get transcript data as a dictionary

    # Find potential matching TADs on the same chromosome
    potential_tads = tads_df[tads_df['chr'] == gene_chr]

    # Check for overlap: gene_start < TAD_end AND gene_end > TAD_start
    # This condition ensures any overlap (contained, crossing start, crossing end, spanning)
    overlapping_tads = potential_tads[
        (potential_tads['hg38_TAD_Start'] < gene_end) &
        (potential_tads['hg38_TAD_End'] > gene_start)
    ]

    if overlapping_tads.empty:
         print(f"Warning: Transcript {gene_row.get('transcript_id', 'N/A')} ({gene_chr}:{gene_start}-{gene_end}) did not overlap any TADs.")
         # Decide how to handle transcripts that don't fit in any TAD.
         # Option 1: Skip them (current behavior implicitly does this if loop below doesn't run)
         # Option 2: Add a row with transcript info and NA for TAD info (more complex)

    # For each TAD this transcript overlaps, create a combined row
    for tad_idx, tad_row in overlapping_tads.iterrows():
        matched_tad_indices.add(tad_idx) # Mark this TAD index as matched
        tad_data = tad_row.to_dict() # Get TAD data as a dictionary

        # Combine TAD data and transcript data
        # The order matters if there are duplicate column names (none expected here)
        combined_row = {**tad_data, **gene_data}
        results.append(combined_row)

# Identify TADs that had no transcripts overlapping them
all_tad_indices = set(tads_df.index)
unmatched_tad_indices = all_tad_indices - matched_tad_indices

print(f"Found {len(results)} transcript-TAD overlaps.")
print(f"Adding {len(unmatched_tad_indices)} TADs with no overlapping transcripts.")

# Add the unmatched TADs with NA values for the transcript columns
for tad_idx in unmatched_tad_indices:
    tad_data = tads_df.loc[tad_idx].to_dict()
    # Combine TAD data with the NA dictionary for transcript columns
    combined_row = {**tad_data, **rna_na_dict}
    results.append(combined_row)

# Create the final DataFrame
final_df = pd.DataFrame(results)

# --- Define final column order ---
# Start with TAD columns, then add RNAseq columns (excluding the temporary 'chr_match')
final_columns = tads_df.columns.tolist() + [col for col in rnaseq_cols if col != 'chr_match']
# Ensure all expected columns are present and handle potential missing ones if logic changes
final_columns = [col for col in final_columns if col in final_df.columns]
# Reorder the DataFrame
final_df = final_df[final_columns]

# Sort results for better readability (optional, e.g., by chromosome and TAD start)
final_df = final_df.sort_values(by=['chr', 'hg38_TAD_Start', 'hg38_start'])


# --- 3. Export the new table ---
print(f"Exporting merged data to: {output_file}")
try:
    final_df.to_csv(output_file, index=False, na_rep='NA') # Use NA for missing values
    print("Export complete.")
except Exception as e:
    print(f"Error exporting file: {e}")

print("Script finished.")

########################################################################################################################

########################################################################################################################

import pandas as pd
import os
import numpy as np  # Retained for pd.NA or np.nan if needed, though pd.NA is preferred

# --- Configuration ---
hic_rnaseq_file = '../public/Multiomics/MCF7_SRA_HiC_RNAseq.csv'
chipseq_file_021 = '../public/ChIPseq/SHK_2024_06_11_MCF7_DBN021_1_0_ug_mL_Rep2 Filtered.csv'
chipseq_file_025 = '../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2 Filtered.csv'  # New file
output_dir = '../public/Multiomics/'
output_file = os.path.join(output_dir, 'MCF7_SRA_HiC_RNAseq_ChIPseq.csv')  # Output filename remains the same

PREFIX_021 = "021_"
PREFIX_025 = "025_"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)


# --- Helper function to load and preprocess ChIP-seq data ---
def load_and_preprocess_chipseq(file_path, file_description):
    """
    Loads a ChIP-seq CSV file, preprocesses it by converting column types,
    handling missing values for key columns, and dropping rows with invalid coordinates.
    """
    print(f"Loading ChIPseq file ({file_description}): {file_path}")
    try:
        df = pd.read_csv(file_path)
        print(f"{file_description} data loaded successfully. Shape: {df.shape}")

        # Check for essential coordinate and chromosome columns
        required_coord_cols = ['hg38_Peak_Start', 'hg38_Peak_End', 'chr']
        for col in required_coord_cols:
            if col not in df.columns:
                print(f"Error: Required column '{col}' not found in {file_description} ({file_path}). Cannot proceed.")
                exit()

        # Convert coordinate columns to numeric, coercing errors
        df['hg38_Peak_Start'] = pd.to_numeric(df['hg38_Peak_Start'], errors='coerce')
        df['hg38_Peak_End'] = pd.to_numeric(df['hg38_Peak_End'], errors='coerce')

        # Handle optional but expected columns (Signal_Value, Peak_Name)
        if 'Signal_Value' in df.columns:
            df['Signal_Value'] = pd.to_numeric(df['Signal_Value'], errors='coerce')
        else:
            print(
                f"Warning: 'Signal_Value' column not found in {file_description} ({file_path}). It will be filled with pd.NA for this source.")
            df['Signal_Value'] = pd.NA

        if 'Peak_Name' not in df.columns:
            print(
                f"Warning: 'Peak_Name' column not found in {file_description} ({file_path}). It will be filled with pd.NA for this source.")
            df['Peak_Name'] = pd.NA

        # Drop rows where essential Peak coordinates or chromosome are missing after coercion
        df.dropna(subset=['chr', 'hg38_Peak_Start', 'hg38_Peak_End'], inplace=True)

        if not df.empty:
            # Convert Peak coordinates to integer type
            df['hg38_Peak_Start'] = df['hg38_Peak_Start'].astype(int)
            df['hg38_Peak_End'] = df['hg38_Peak_End'].astype(int)
        else:
            print(
                f"Warning: {file_description} DataFrame is empty after dropping NA from essential coordinate columns.")
            # Ensure all expected columns exist even if DataFrame is empty, for consistency downstream
            expected_cols = ['chr', 'Peak_Name', 'hg38_Peak_Start', 'hg38_Peak_End', 'Signal_Value']
            for col_to_ensure in expected_cols:
                if col_to_ensure not in df.columns:
                    # Define dtype for empty series to avoid issues
                    dtype = 'object' if col_to_ensure in ['Peak_Name', 'chr'] else 'float64'  # float64 can hold pd.NA
                    df[col_to_ensure] = pd.Series(dtype=dtype)
        return df
    except FileNotFoundError:
        print(f"Error: {file_description} file not found at {file_path}")
        exit()
    except Exception as e:
        print(f"Error loading {file_description} file ({file_path}): {e}")
        exit()


# --- 1. Import and preprocess HiC/RNAseq data ---
print(f"Loading HiC/RNAseq file: {hic_rnaseq_file}")
try:
    hic_rnaseq_df = pd.read_csv(hic_rnaseq_file)
    print(f"HiC/RNAseq data loaded successfully. Shape: {hic_rnaseq_df.shape}")

    # Ensure coordinate columns are numeric, coercing errors to NaN
    hic_rnaseq_df['hg38_TAD_Start'] = pd.to_numeric(hic_rnaseq_df['hg38_TAD_Start'], errors='coerce')
    hic_rnaseq_df['hg38_TAD_End'] = pd.to_numeric(hic_rnaseq_df['hg38_TAD_End'], errors='coerce')

    # Process other coordinate columns if they exist
    if 'hg38_start' in hic_rnaseq_df.columns:  # Gene start
        hic_rnaseq_df['hg38_start'] = pd.to_numeric(hic_rnaseq_df['hg38_start'], errors='coerce')
    if 'hg38_end' in hic_rnaseq_df.columns:  # Gene end
        hic_rnaseq_df['hg38_end'] = pd.to_numeric(hic_rnaseq_df['hg38_end'], errors='coerce')
    if 'hg38_chr' in hic_rnaseq_df.columns:  # Gene chromosome (often numeric like '1' vs 'chr1')
        hic_rnaseq_df['hg38_chr'] = pd.to_numeric(hic_rnaseq_df['hg38_chr'],
                                                  errors='coerce')  # Keep as float to allow NA

    # Drop rows where essential TAD coordinates are missing
    hic_rnaseq_df.dropna(subset=['chr', 'hg38_TAD_Start', 'hg38_TAD_End'], inplace=True)
    if not hic_rnaseq_df.empty:
        # Convert TAD coordinates to integer type after handling NaNs
        hic_rnaseq_df['hg38_TAD_Start'] = hic_rnaseq_df['hg38_TAD_Start'].astype(int)
        hic_rnaseq_df['hg38_TAD_End'] = hic_rnaseq_df['hg38_TAD_End'].astype(int)
    else:
        print("Warning: HiC/RNAseq DataFrame is empty after dropping NA from essential coordinate columns.")

except FileNotFoundError:
    print(f"Error: HiC/RNAseq file not found at {hic_rnaseq_file}")
    exit()
except Exception as e:
    print(f"Error loading HiC/RNAseq file: {e}")
    exit()

# --- 2. Import and preprocess ChIPseq data ---
chipseq_df_021 = load_and_preprocess_chipseq(chipseq_file_021, "ChIPseq 021")
chipseq_df_025 = load_and_preprocess_chipseq(chipseq_file_025, "ChIPseq 025")

# Define columns from ChIP-seq data that will be added (excluding 'chr')
# These are the original names; they will be prefixed later.
# load_and_preprocess_chipseq ensures these columns exist in the dataframes, possibly filled with NA.
chipseq_cols_to_add = ['Peak_Name', 'hg38_Peak_Start', 'hg38_Peak_End', 'Signal_Value']

# --- 3. Merge the information ---
results = []
print("Merging ChIP-seq peaks into TADs based on genomic coordinates...")

if hic_rnaseq_df.empty:
    print("HiC/RNAseq data is empty. No merging will be performed.")
else:
    # Iterate through each TAD
    for tad_idx, tad_row in hic_rnaseq_df.iterrows():
        if tad_idx > 0 and tad_idx % 500 == 0:  # Progress update
            print(f"Processing TAD {tad_idx}/{len(hic_rnaseq_df)}...")

        tad_data = tad_row.to_dict()
        tad_chr = tad_row['chr']
        tad_start = tad_row['hg38_TAD_Start']
        tad_end = tad_row['hg38_TAD_End']

        # --- Process ChIP-seq 021 ---
        processed_peaks_021 = []
        if not chipseq_df_021.empty:
            potential_peaks_021 = chipseq_df_021[chipseq_df_021['chr'] == tad_chr]
            if not potential_peaks_021.empty:
                overlapping_peaks_021_df = potential_peaks_021[
                    (potential_peaks_021['hg38_Peak_Start'] < tad_end) &
                    (potential_peaks_021['hg38_Peak_End'] > tad_start)
                    ]
                if not overlapping_peaks_021_df.empty:
                    for _, peak_row in overlapping_peaks_021_df.iterrows():
                        peak_info = {PREFIX_021 + col: peak_row.get(col, pd.NA) for col in chipseq_cols_to_add}
                        processed_peaks_021.append(peak_info)

        if not processed_peaks_021:  # If no overlaps or chipseq_df_021 was empty
            na_dict_021 = {PREFIX_021 + col: pd.NA for col in chipseq_cols_to_add}
            processed_peaks_021.append(na_dict_021)

        # --- Process ChIP-seq 025 ---
        processed_peaks_025 = []
        if not chipseq_df_025.empty:
            potential_peaks_025 = chipseq_df_025[chipseq_df_025['chr'] == tad_chr]
            if not potential_peaks_025.empty:
                overlapping_peaks_025_df = potential_peaks_025[
                    (potential_peaks_025['hg38_Peak_Start'] < tad_end) &
                    (potential_peaks_025['hg38_Peak_End'] > tad_start)
                    ]
                if not overlapping_peaks_025_df.empty:
                    for _, peak_row in overlapping_peaks_025_df.iterrows():
                        peak_info = {PREFIX_025 + col: peak_row.get(col, pd.NA) for col in chipseq_cols_to_add}
                        processed_peaks_025.append(peak_info)

        if not processed_peaks_025:  # If no overlaps or chipseq_df_025 was empty
            na_dict_025 = {PREFIX_025 + col: pd.NA for col in chipseq_cols_to_add}
            processed_peaks_025.append(na_dict_025)

        # --- Combine TAD data with processed peak data (Cartesian product) ---
        for p_021_data in processed_peaks_021:
            for p_025_data in processed_peaks_025:
                combined_row = {**tad_data, **p_021_data, **p_025_data}
                results.append(combined_row)

print(f"Generated {len(results)} combined rows.")

# --- 4. Create Final DataFrame ---
# Define base column lists
hic_cols = []
if not hic_rnaseq_df.empty:
    hic_cols = hic_rnaseq_df.columns.tolist()
else:  # Fallback if hic_rnaseq_df was empty but we need its schema
    try:
        hic_cols = pd.read_csv(hic_rnaseq_file, nrows=0).columns.tolist()
        print("Used header from HiC/RNAseq file for empty DataFrame schema.")
    except Exception as e:
        print(f"Could not read header from HiC/RNAseq file for schema: {e}")
        # Define a minimal set of columns if schema cannot be read
        hic_cols = ['chr', 'TAD_Num', 'hg38_TAD_Start', 'hg38_TAD_End']

chip_021_final_cols = [PREFIX_021 + col for col in chipseq_cols_to_add]
chip_025_final_cols = [PREFIX_025 + col for col in chipseq_cols_to_add]
all_final_cols = hic_cols + chip_021_final_cols + chip_025_final_cols

if not results:
    print("Warning: No results generated. The output file will contain only headers if possible, or be empty.")
    final_df = pd.DataFrame(columns=all_final_cols)
else:
    final_df = pd.DataFrame(results)
    # Ensure all expected columns are present, reorder, and fill missing ones with pd.NA
    for col in all_final_cols:
        if col not in final_df.columns:
            final_df[col] = pd.NA
    final_df = final_df[all_final_cols]

# --- Sort results ---
if not final_df.empty:
    # For sorting, we need the prefixed peak start columns
    sort_col_021 = PREFIX_021 + 'hg38_Peak_Start'
    sort_col_025 = PREFIX_025 + 'hg38_Peak_Start'

    # Ensure sort columns exist, add them as NA if not (e.g. if a chipseq file was totally empty)
    if sort_col_021 not in final_df.columns: final_df[sort_col_021] = pd.NA
    if sort_col_025 not in final_df.columns: final_df[sort_col_025] = pd.NA

    final_df = final_df.sort_values(
        by=['chr', 'hg38_TAD_Start', sort_col_021, sort_col_025],
        na_position='first'  # Puts NAs at the beginning for peak starts
    )
else:
    print("Final DataFrame is empty, skipping sort.")

# --- 5. Export the new table ---
print(f"Exporting merged data to: {output_file}")
try:
    final_df.to_csv(output_file, index=False, na_rep='NA')  # Use 'NA' for string representation
    print("Export complete.")
except Exception as e:
    print(f"Error exporting file: {e}")

print("Script finished.")

########################################################################################################################

########################################################################################################################

import pandas as pd
import numpy as np  # numpy is often used with pandas, especially for 'NA' like values


def process_multiomics_data(file_path):
    """
    Processes a multiomics CSV file to conditionally replace values in specific columns with 'NA'.

    Args:
        file_path (str): The path to the CSV file.
    """
    try:
        # 1. Import the CSV file
        # Read the CSV file into a pandas DataFrame.
        # We use na_filter=False to prevent pandas from interpreting "NA" strings as actual NaN values initially,
        # so we can reliably check for the string "NA".
        # We also specify dtype=str to read all columns as strings initially,
        # to avoid issues with mixed types or pandas inferring types incorrectly for this specific task.
        df = pd.read_csv(file_path, na_filter=False, dtype=str)
        print(f"Successfully read {len(df)} rows from {file_path}")

        # Define the columns to be set to "NA" if the condition is met
        columns_to_na = [
            'transcript_id', 'Symbol', 'hg38_chr', 'hg38_start', 'hg38_end',
            'hg38_strand', '023_24hr_Up', '025_24hr_Up', '021_10hr_Up',
            '021_24hr_Up', '021_48hr_Up', 'Category - SRA specific'
            # 'Category - other' will also be set to 'NA' as per the example logic,
            # even though it's the condition column.
        ]

        # Column that triggers the NA conversion
        condition_column = "Category - other"

        # 2. For each row, check the value of "Category - other" column
        # We iterate using df.iterrows() for clarity, though for very large files,
        # vectorized operations would be more performant.

        rows_changed_count = 0
        for index, row in df.iterrows():
            # Check if the "Category - other" column value is not "NA" (as a string)
            if row[condition_column] != "NA" and pd.notna(row[condition_column]) and row[
                condition_column].strip() != "":
                # If the condition is met, convert the values of the specified columns to "NA" (as a string)
                for col in columns_to_na:
                    if col in df.columns:  # Ensure column exists
                        df.loc[index, col] = "NA"
                    else:
                        print(f"Warning: Column '{col}' not found in the DataFrame.")
                # Also set the "Category - other" column itself to "NA" as per the example's output
                df.loc[index, condition_column] = "NA"
                rows_changed_count += 1

        print(f"Processed {rows_changed_count} rows based on the condition.")

        # 3. Save the filtered file
        # Save the modified DataFrame back to the same CSV file, overwriting it.
        # We ensure that "NA" is written as the string "NA" and not as an empty field.
        df.to_csv(file_path, index=False, na_rep='NA')
        print(f"Successfully saved the processed data to {file_path}")

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{file_path}' is empty.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    # Define the path to your CSV file
    # IMPORTANT: Adjust this path if your script is not in the same directory level
    # as the parent of 'public'. For example, if your script is in a 'scripts' folder
    # and 'public' is a sibling to 'scripts', the path might be '../public/...'
    # If the script is in the root and 'public' is a subfolder, it might be 'public/...'

    # Assuming the script is located in a way that '../public/' correctly points to the target directory.
    # If your script is, for example, in /some_project/scripts/
    # and the file is in /some_project/public/Multiomics/
    # then the relative path from /some_project/scripts/ would be '../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq.csv'

    # For the provided path "../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq.csv"
    # This implies the script is located one directory level above the 'public' directory's parent.
    # Or, more likely, if your project structure is:
    # project_root/
    #   scripts/ (your_script.py is here)
    #   public/
    #     Multiomics/
    #       MCF7_SRA_HiC_RNAseq_ChIPseq.csv
    # Then the path from your_script.py would be '../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq.csv'

    # If your script is in the same directory as the 'public' folder (e.g. project_root/your_script.py)
    # then the path would be 'public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq.csv'

    # Please verify this path based on your actual file structure.
    csv_file_path = "../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq.csv"

    # Create a dummy CSV for testing if the actual file doesn't exist or for development
    # This part is for testing and can be removed if you have the actual file.
    try:
        pd.read_csv(csv_file_path, nrows=1)
    except FileNotFoundError:
        print(f"Warning: File {csv_file_path} not found. Creating a dummy file for demonstration.")
        # Create a dummy directory structure if it doesn't exist
        import os

        os.makedirs(os.path.dirname(csv_file_path), exist_ok=True)

        dummy_data = {
            'chr': ['chr1', 'chr1', 'chr1', 'chr1'],
            'TAD_Num': [1, 2, 10, 9],
            'hg38_TAD_Start': [1, 1808562, 8399941, 7759941],
            'hg38_TAD_End': [1808562, 2128562, 9639943, 8399941],
            'Insulation_Score_5\'': [0.0, 0.259062031198192, 0.278826252364477, 0.5069133690219],
            'Insulation_Score_3\'': [0.259062031198192, 0.196631710459799, 0.287701186416925, 0.278826252364477],
            'transcript_id': ['ENSG00000269981.1', 'NA', 'ENSG00000171608.15', 'ENSG00000162426.14'],
            'Symbol': ['NA', 'NA', 'PIK3CD', 'SLC45A1'],
            'hg38_chr': [1.0, 'NA', 1.0, 1.0],
            'hg38_start': [137682.0, 'NA', 9629889.0, 8318114.0],
            'hg38_end': [137965.0, 'NA', 9729114.0, 8344167.0],
            'hg38_strand': ['reverse strand', 'NA', 'forward strand', 'forward strand'],
            '023_24hr_Up': [0.0, 'NA', 0.0, 0.0],
            '025_24hr_Up': [1.0, 'NA', 0.0, 1.0],
            '021_10hr_Up': [0.0, 'NA', 1.0, 0.0],
            '021_24hr_Up': [0.0, 'NA', 0.0, 0.0],
            '021_48hr_Up': [0.0, 'NA', 0.0, 0.0],
            'Category - SRA specific': ['NA', 'NA', 'Early transient', 'NA'],
            'Category - other': ['PCD-RFP', 'NA', 'NA', 'PCD-RFP'],  # First and last row will be modified
            '021_Peak_Name': ['NA', 'NA', 'NA', 'NA'],
            '021_hg38_Peak_Start': ['NA', 'NA', 'NA', 'NA'],
            '021_hg38_Peak_End': ['NA', 'NA', 'NA', 'NA'],
            '021_Signal_Value': ['NA', 'NA', 'NA', 'NA'],
            '025_Peak_Name': ['NA', 'NA', 'chr1_1', 'NA'],
            '025_hg38_Peak_Start': ['NA', 'NA', 8753611, 'NA'],
            '025_hg38_Peak_End': ['NA', 'NA', 8754169, 'NA'],
            '025_Signal_Value': ['NA', 'NA', 3.14876, 'NA']
        }
        dummy_df = pd.DataFrame(dummy_data)
        # Ensure all columns from the original header are present, filling missing ones with "NA"
        header = "chr,TAD_Num,hg38_TAD_Start,hg38_TAD_End,Insulation_Score_5',Insulation_Score_3',transcript_id,Symbol,hg38_chr,hg38_start,hg38_end,hg38_strand,023_24hr_Up,025_24hr_Up,021_10hr_Up,021_24hr_Up,021_48hr_Up,Category - SRA specific,Category - other,021_Peak_Name,021_hg38_Peak_Start,021_hg38_Peak_End,021_Signal_Value,025_Peak_Name,025_hg38_Peak_Start,025_hg38_Peak_End,025_Signal_Value".split(
            ',')
        for col in header:
            if col not in dummy_df.columns:
                dummy_df[col] = "NA"
        dummy_df = dummy_df[header]  # Ensure correct column order
        dummy_df.to_csv(csv_file_path, index=False, na_rep='NA')
        print(f"Created and populated dummy file: {csv_file_path}")

    process_multiomics_data(csv_file_path)

    # Optional: Print the first few lines of the modified file to verify
    try:
        df_modified = pd.read_csv(csv_file_path, na_filter=False)  # Read again to check, keep "NA" as string
        print("\nFirst 5 rows of the modified file:")
        print(df_modified.head().to_string())
    except Exception as e:
        print(f"Error reading modified file for verification: {e}")

########################################################################################################################

########################################################################################################################

import pandas as pd
import os

# Define the input and output file paths
# Assuming the script is run from a directory where '../public/Multiomics/' is a valid relative path
# If not, you might need to adjust the base path or use absolute paths.
input_file_path = '../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq.csv'
output_file_path = '../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq_UpDEG_Positive.csv'

try:
    # Read the input CSV file into a pandas DataFrame
    # It's good practice to specify the dtype for columns if known,
    # especially if 'NA' might be misinterpreted as a missing value indicator
    # by default. However, based on the example, 'NA' seems to be a string.
    df = pd.read_csv(input_file_path)

    # Filter the DataFrame to keep only rows where 'transcript_id' is not 'NA'
    # The .ne() method checks for 'not equal to'
    filtered_df = df[df['transcript_id'].ne('NA') & df['transcript_id'].notna()]

    # Ensure the output directory exists, create it if it doesn't
    output_dir = os.path.dirname(output_file_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")

    # Save the filtered DataFrame to a new CSV file
    # index=False prevents pandas from writing the DataFrame index as a column
    filtered_df.to_csv(output_file_path, index=False)

    print(f"Successfully read '{input_file_path}'.")
    print(f"Filtered data saved to '{output_file_path}'.")
    print(f"Original number of rows: {len(df)}")
    print(f"Number of rows after filtering (transcript_id != 'NA'): {len(filtered_df)}")

except FileNotFoundError:
    print(f"Error: Input file not found at '{input_file_path}'. Please check the path.")
except Exception as e:
    print(f"An error occurred: {e}")

########################################################################################################################

########################################################################################################################

# Import necessary library
import pandas as pd
import numpy as np
import os

# Define input and output file paths
# Assuming the script is run from a directory where '../public/Multiomics/' is accessible
input_file = '../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq.csv'
output_dir = '../public/Multiomics/'
output_file = os.path.join(output_dir, 'MCF7_SRA_HiC_RNAseq_ChIPseq_Condensed.csv')

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

print(f"Reading input file: {input_file}")

try:
    # Read the input CSV file into a pandas DataFrame
    df = pd.read_csv(input_file)

    print("Input data loaded successfully. Columns:", df.columns.tolist())
    print("Input data shape:", df.shape)
    # print("First 5 rows of input data:\n", df.head()) # Can be verbose

    # Replace potential string 'NA' values with actual NaN values for easier handling
    # This helps in nunique() and dropna() operations
    df.replace('NA', np.nan, inplace=True)
    # Also replace empty strings if they might occur and should be treated as NA
    df.replace('', np.nan, inplace=True)

    # Define the aggregation logic
    # For columns not specified here, they will be dropped by the groupby().agg() operation
    # unless they are numeric and a default aggregation is applied (not the case here).
    aggregation_functions = {
        'hg38_TAD_Start': 'first',
        'hg38_TAD_End': 'first',
        "Insulation_Score_5'": 'first',
        "Insulation_Score_3'": 'first',
        'transcript_id': lambda x: ' '.join(x.dropna().astype(str).unique()),
        'Symbol': lambda x: ' '.join(x.dropna().astype(str).unique()),
        # Add new peak name aggregations for 021 and 025 datasets
        '021_Peak_Name': lambda x: ' '.join(x.dropna().astype(str).unique()),
        '025_Peak_Name': lambda x: ' '.join(x.dropna().astype(str).unique())
        # Note: Other columns from the input (like 021_hg38_Peak_Start, 021_Signal_Value, etc.)
        # will not be included in the condensed output unless added to aggregation_functions.
        # This script focuses on condensing to TAD-level info with peak names and counts.
    }

    # Identify columns that are actually present in the DataFrame to avoid KeyErrors
    # This is important if the input CSV might sometimes lack optional columns
    # For this script, we expect the columns from the previous script to be present.
    cols_to_aggregate = {k: v for k, v in aggregation_functions.items() if k in df.columns}
    if len(cols_to_aggregate) < len(aggregation_functions):
        missing_cols = set(aggregation_functions.keys()) - set(cols_to_aggregate.keys())
        print(
            f"Warning: The following columns specified for aggregation were not found in the input file and will be skipped: {missing_cols}")

    print("Grouping data by 'chr' and 'TAD_Num' and aggregating...")

    # Group by 'chr' and 'TAD_Num' and apply the aggregation
    # Using observed=True for potential performance improvement if categories are used (default in newer pandas)
    if not df.empty and cols_to_aggregate:
        grouped_df = df.groupby(['chr', 'TAD_Num'], observed=True, as_index=False).agg(cols_to_aggregate)
    elif not df.empty and not cols_to_aggregate:
        print("Warning: No columns to aggregate. Creating a DataFrame with group keys only.")
        grouped_df = df[['chr', 'TAD_Num']].drop_duplicates().reset_index(drop=True)
    else:
        print("Input DataFrame is empty. Creating an empty grouped DataFrame.")
        grouped_df = pd.DataFrame(columns=['chr', 'TAD_Num'] + list(cols_to_aggregate.keys()))

    print("Aggregation complete. Calculating counts...")

    # Calculate unique non-NA transcript count for each group
    if 'transcript_id' in df.columns:  # Ensure column exists before grouping
        transcript_counts = df.groupby(['chr', 'TAD_Num'], observed=True)['transcript_id'].nunique().reset_index(
            name='transcript_Count')
        grouped_df = pd.merge(grouped_df, transcript_counts, on=['chr', 'TAD_Num'], how='left')
    else:
        grouped_df['transcript_Count'] = 0  # Or pd.NA if preferred

    # Calculate unique non-NA peak count for 021 data
    if '021_Peak_Name' in df.columns:
        peak_counts_021 = df.groupby(['chr', 'TAD_Num'], observed=True)['021_Peak_Name'].nunique().reset_index(
            name='021_Peak_Count')
        grouped_df = pd.merge(grouped_df, peak_counts_021, on=['chr', 'TAD_Num'], how='left')
    else:
        grouped_df['021_Peak_Count'] = 0  # Or pd.NA

    # Calculate unique non-NA peak count for 025 data
    if '025_Peak_Name' in df.columns:
        peak_counts_025 = df.groupby(['chr', 'TAD_Num'], observed=True)['025_Peak_Name'].nunique().reset_index(
            name='025_Peak_Count')
        grouped_df = pd.merge(grouped_df, peak_counts_025, on=['chr', 'TAD_Num'], how='left')
    else:
        grouped_df['025_Peak_Count'] = 0  # Or pd.NA

    # condensed_df is now grouped_df after merges
    condensed_df = grouped_df

    print("Counts calculated. Renaming and reordering columns...")

    # Rename columns to match the desired output format
    # The aggregated columns are already named (e.g., '021_Peak_Name' from agg).
    # We rename them to plural form (e.g., '021_Peak_Names').
    rename_map = {
        'transcript_id': 'transcript_ids',  # From aggregation_functions keys
        'Symbol': 'Symbols',  # From aggregation_functions keys
        '021_Peak_Name': '021_Peak_Names',  # From aggregation_functions keys
        '025_Peak_Name': '025_Peak_Names'  # From aggregation_functions keys
    }
    # Filter rename_map for columns that actually exist in condensed_df
    actual_rename_map = {k: v for k, v in rename_map.items() if k in condensed_df.columns}
    condensed_df.rename(columns=actual_rename_map, inplace=True)

    # Define the desired column order
    output_columns = [
        'chr',
        'TAD_Num',
        'hg38_TAD_Start',
        'hg38_TAD_End',
        "Insulation_Score_5'",
        "Insulation_Score_3'",
        'transcript_ids',  # Renamed from transcript_id
        'Symbols',  # Renamed from Symbol
        'transcript_Count',  # Calculated
        '021_Peak_Names',  # Renamed from 021_Peak_Name
        '021_Peak_Count',  # Calculated
        '025_Peak_Names',  # Renamed from 025_Peak_Name
        '025_Peak_Count'  # Calculated
    ]

    # Ensure all output columns exist in condensed_df, add them with NA if missing
    for col in output_columns:
        if col not in condensed_df.columns:
            print(f"Warning: Expected output column '{col}' not found after aggregation/renaming. Adding it as NA.")
            condensed_df[col] = np.nan  # Use np.nan, will be converted to 'NA' string later

    # Reorder columns to the defined order
    condensed_df = condensed_df[output_columns]

    print("Columns processed. Final DataFrame shape:", condensed_df.shape)
    # print("First 5 rows of condensed data:\n", condensed_df.head()) # Can be verbose

    # Export the condensed DataFrame to a new CSV file
    # Replace NaN values back to 'NA' string for the output CSV
    condensed_df.fillna('NA', inplace=True)
    condensed_df.to_csv(output_file, index=False, na_rep='NA')

    print(f"Successfully condensed data and saved to: {output_file}")

except FileNotFoundError:
    print(f"Error: Input file not found at {input_file}")
    exit()  # Exit if input file is critical and not found
except KeyError as e:
    print(f"An error occurred due to a missing column (KeyError): {e}")
    print("This might be due to the input CSV not having expected columns.")
    exit()
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    exit()

########################################################################################################################

########################################################################################################################

import pandas as pd
import os
import time

print("Starting script...")

# --- Define file paths ---
# Main input and output files
main_multiomics_file = '../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq_Condensed.csv'
output_file = '../public/Multiomics/MCF7_SRA_HiC_RNAseq_ChIPseq_Condensed_HistoneMarks.csv'

# New list of 7 ChIP-seq peak files
chipseq_files = {
    'H3K4me1': "../public/Public Dataset/ChIPseq/GSE86714 ENCFF991HJA hg38 MCF7 H3K4me1 ChIPseq Filtered.csv",
    'H3K4me2': "../public/Public Dataset/ChIPseq/GSE96439 ENCFF188VRU hg38 MCF7 H3K4me2 ChIPseq Filtered.csv",
    'H3K4me3': "../public/Public Dataset/ChIPseq/GSE96506 ENCFF268RXB hg38 MCF7 H3K4me3 ChIPseq Filtered.csv",
    'H3K9me3': "../public/Public Dataset/ChIPseq/GSE96517 ENCFF501UHK hg38 MCF7 H3K9me3 ChIPseq Filtered.csv",
    'H3K27ac': "../public/Public Dataset/ChIPseq/GSE96352 ENCFF491LQY hg38 MCF7 H3K27ac ChIPseq Filtered.csv",
    'H3K27me3': "../public/Public Dataset/ChIPseq/GSE96363 ENCFF669NUD hg38 MCF7 H3K27me3 ChIPseq Filtered.csv",
    'H3K36me3': "../public/Public Dataset/ChIPseq/GSE174945 ENCFF195FSD hg38 MCF7 H3K36me3 ChIPseq Filtered.csv"
}


# Create output directory if it doesn't exist
output_dir = os.path.dirname(output_file)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created output directory: {output_dir}")

# --- 1. Import files ---
try:
    print(f"Loading main multiomics file: {main_multiomics_file}")
    main_df = pd.read_csv(main_multiomics_file)
    print(f"Loaded main multiomics file. Shape: {main_df.shape}")

    # Load all ChIP-seq files into a dictionary of dataframes
    chipseq_dfs = {}
    for mark, file_path in chipseq_files.items():
        print(f"Loading {mark} peaks: {file_path}")
        chipseq_dfs[mark] = pd.read_csv(file_path)
        print(f"Loaded {mark} peaks. Shape: {chipseq_dfs[mark].shape}")

except FileNotFoundError as e:
    print(f"Error loading file: {e}")
    print("Please ensure all input files exist at the specified paths.")
    exit()  # Exit the script if files are missing
except Exception as e:
    print(f"An unexpected error occurred during file loading: {e}")
    exit()

# --- 2. Add new columns and populate counts ---

# Add the new peak count columns and initialize with 0
print("\nAdding new peak count columns to the main DataFrame...")
for mark in chipseq_files.keys():
    column_name = f'{mark}_Peak_Count'
    main_df[column_name] = 0
    print(f"  Added column: {column_name}")
print("New columns added.")


# Define a function to process peaks for a given histone mark
def add_peak_counts(main_tad_df, peak_df, count_column_name):
    """
    Adds peak counts to the main TAD DataFrame based on peak locations.

    Args:
        main_tad_df (pd.DataFrame): The DataFrame containing TAD information.
        peak_df (pd.DataFrame): The DataFrame containing peak information
                                (chr, hg38_Peak_Start, hg38_Peak_End).
        count_column_name (str): The name of the column in main_tad_df
                                 to increment.
    """
    print(f"\nProcessing peaks for {count_column_name}...")
    start_time = time.time()
    peaks_processed = 0
    total_peaks = len(peak_df)

    # Ensure correct data types for comparison
    main_tad_df['hg38_TAD_Start'] = main_tad_df['hg38_TAD_Start'].astype(int)
    main_tad_df['hg38_TAD_End'] = main_tad_df['hg38_TAD_End'].astype(int)
    peak_df['hg38_Peak_Start'] = peak_df['hg38_Peak_Start'].astype(int)
    peak_df['hg38_Peak_End'] = peak_df['hg38_Peak_End'].astype(int)

    # Group TADs by chromosome for faster lookup
    tads_by_chrom = main_tad_df.groupby('chr')

    # Iterate through each peak
    for index, peak in peak_df.iterrows():
        peak_chr = peak['chr']
        peak_start = peak['hg38_Peak_Start']
        peak_end = peak['hg38_Peak_End']

        # Find relevant TADs only on the same chromosome
        if peak_chr in tads_by_chrom.groups:
            relevant_tads = tads_by_chrom.get_group(peak_chr)

            # Find TADs that overlap with the peak
            # Overlap condition: TAD start < peak end AND TAD end > peak start
            overlapping_tads_indices = relevant_tads[
                (relevant_tads['hg38_TAD_Start'] < peak_end) &
                (relevant_tads['hg38_TAD_End'] > peak_start)
                ].index

            # Increment the count for each overlapping TAD in the original main_tad_df
            if not overlapping_tads_indices.empty:
                main_tad_df.loc[overlapping_tads_indices, count_column_name] += 1

        peaks_processed += 1
        if peaks_processed % 10000 == 0:  # Print progress update
            print(f"  Processed {peaks_processed}/{total_peaks} peaks for {count_column_name}...")

    end_time = time.time()
    print(f"Finished processing {count_column_name}. Time taken: {end_time - start_time:.2f} seconds.")
    print(f"Summary for {count_column_name}: Non-zero counts in {(main_tad_df[count_column_name] > 0).sum()} TADs.")


# Process peaks for all 7 histone marks
for mark, peak_df in chipseq_dfs.items():
    column_name = f'{mark}_Peak_Count'
    add_peak_counts(main_df, peak_df, column_name)


# --- 3. Export the updated table ---
print(f"\nExporting updated DataFrame to: {output_file}")
try:
    main_df.to_csv(output_file, index=False)
    print("Export successful.")
    print(f"\nFinal DataFrame head:\n{main_df.head()}")
    print(f"Final DataFrame shape: {main_df.shape}")

except Exception as e:
    print(f"An error occurred during file export: {e}")

print("\nScript finished.")