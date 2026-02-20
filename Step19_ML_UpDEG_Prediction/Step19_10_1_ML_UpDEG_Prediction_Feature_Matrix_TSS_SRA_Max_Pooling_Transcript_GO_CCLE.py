import polars as pl
import time
import os

# --- Configuration ---
# (These are the file paths from your *first* script)
ms_file_path = '../../public/Public Dataset/Protein_Array/harmonized_MS_CCLE_Gygi_MCF7.csv'
mapping_file_path = '../../public/Public Dataset/Protein_Array/uniprot_hugo_entrez_id_mapping.csv'
updeg_file_path = '../../public/Multiomics/Step19_ML_UpDEG_Prediction/SRA_Specific_UpDEG_with_GO.csv'
nondeg_file_path = '../../public/Multiomics/Step19_ML_UpDEG_Prediction/Non_DEG_with_GO.csv'

output_ms_processed_file = '../../public/Public Dataset/Protein_Array/harmonized_MS_CCLE_Gygi_MCF7_Processed.csv'
output_updeg_file = '../../public/Multiomics/Step19_ML_UpDEG_Prediction/SRA_Specific_UpDEG_with_GO_CCLE.csv'
output_nondeg_file = '../../public/Multiomics/Step19_ML_UpDEG_Prediction/Non_DEG_with_GO_CCLE.csv'

# --- Main Processing Logic (using Polars) ---

print("Starting data processing (using Polars)...")
start_time = time.time()

# 1. & 2. Load MS data and filter out rows with missing MS_Value
try:
    print(f"Loading MS data from: {ms_file_path}")
    # Use scan_csv (lazy) and drop_nulls, then collect (eager)
    ms_data_pl = pl.scan_csv(ms_file_path).drop_nulls(subset=['MS_Value']).collect()
    print(f"  Shape after loading and removing missing MS_Values: {ms_data_pl.shape}")

except FileNotFoundError:
    print(f"Error: File not found at {ms_file_path}")
    exit()
except Exception as e:
    print(f"An error occurred while reading {ms_file_path}: {e}")
    exit()

# 3. Load Mapping Data
try:
    print(f"Loading mapping data from: {mapping_file_path}")
    # Lazily scan, select only needed columns, use unique() for lazy, then collect
    mapping_pl = (
        pl.scan_csv(mapping_file_path)
        .select(['UniprotID', 'Symbol'])
        .unique(subset=['UniprotID'], keep='first') # <-- FIX: Was .drop_duplicates
        .collect()
    )
    print(f"  Mapping data shape: {mapping_pl.shape}")

except FileNotFoundError:
    print(f"Error: Mapping file not found at {mapping_file_path}")
    exit()
except Exception as e:
    print(f"An error occurred while processing {mapping_file_path}: {e}")
    exit()

# 4. & 5. Map Symbols to MS data using a Polars 'join'
print("Mapping gene symbols to MS data...")
# An 'inner' join is the Polars equivalent of mapping and then dropping NaNs
ms_data_processed_pl = ms_data_pl.join(
    mapping_pl, on='UniprotID', how='inner'
)

# Reorder columns using 'select'
ms_data_processed_pl = ms_data_processed_pl.select(
    ['Symbol', 'UniprotID', 'MS_Value']
)
print(f"  Shape after mapping and dropping unmapped: {ms_data_processed_pl.shape}")

# --- 6. Export the processed MS data ---
try:
    print(f"Saving processed MS data to: {output_ms_processed_file}")
    # Use polars .write_csv
    ms_data_processed_pl.write_csv(output_ms_processed_file)
except Exception as e:
    print(f"Error saving processed MS data: {e}")
    exit()

# --- Prepare MS data for merging with DEG data ---
# We only need the 'Symbol' and 'MS_Value' for the next join
ms_values_to_join_pl = ms_data_processed_pl.select(['Symbol', 'MS_Value'])

# --- 7. & 8. Process SRA_Specific_UpDEG_with_GO ---
print(f"\nProcessing UpDEG file: {updeg_file_path}")
try:
    # Use read_csv (eager) to immediately check columns
    updeg_data_pl = pl.read_csv(updeg_file_path)
    print(f"  Original shape: {updeg_data_pl.shape}")

    # Check if 'symbol' column exists
    if 'symbol' not in updeg_data_pl.columns:
        raise ValueError(f"'symbol' column not found in {updeg_file_path}")

    # Rename 'symbol' to 'Symbol' to match for joining
    updeg_data_pl = updeg_data_pl.rename({'symbol': 'Symbol'})

    # 11 & 12. Perform an inner join on 'Symbol'
    updeg_merged_pl = updeg_data_pl.join(
        ms_values_to_join_pl, on='Symbol', how='inner'
    )

    print(f"  Shape after merging with MS data: {updeg_merged_pl.shape}")

    # 13. Export the final UpDEG DataFrame
    print(f"Saving processed UpDEG data to: {output_updeg_file}")
    updeg_merged_pl.write_csv(output_updeg_file)

except FileNotFoundError:
    print(f"Error: File not found at {updeg_file_path}")
except ValueError as ve:
    print(ve)
except Exception as e:
    print(f"An error occurred while processing {updeg_file_path}: {e}")

# --- 9. & 10. Process Non_DEG_with_GO ---
print(f"\nProcessing Non-DEG file: {nondeg_file_path}")
try:
    # Use read_csv (eager) to immediately check columns
    nondeg_data_pl = pl.read_csv(nondeg_file_path)
    print(f"  Original shape: {nondeg_data_pl.shape}")

    # Check if 'symbol' column exists
    if 'symbol' not in nondeg_data_pl.columns:
        raise ValueError(f"'symbol' column not found in {nondeg_file_path}")

    # Rename 'symbol' to 'Symbol' to match for joining
    nondeg_data_pl = nondeg_data_pl.rename({'symbol': 'Symbol'})

    # 11 & 12. Perform an inner join on 'Symbol'
    nondeg_merged_pl = nondeg_data_pl.join(
        ms_values_to_join_pl, on='Symbol', how='inner'
    )

    print(f"  Shape after merging with MS data: {nondeg_merged_pl.shape}")

    # 14. Export the final Non-DEG DataFrame
    print(f"Saving processed Non-DEG data to: {output_nondeg_file}")
    nondeg_merged_pl.write_csv(output_nondeg_file)

except FileNotFoundError:
    print(f"Error: File not found at {nondeg_file_path}")
except ValueError as ve:
    print(ve)
except Exception as e:
    print(f"An error occurred while processing {nondeg_file_path}: {e}")

end_time = time.time()
print(f"\nProcessing finished in {end_time - start_time:.2f} seconds.")


