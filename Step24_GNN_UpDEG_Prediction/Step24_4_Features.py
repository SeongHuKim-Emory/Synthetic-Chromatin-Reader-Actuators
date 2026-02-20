import polars as pl
import pyBigWig
import os
import sys
from multiprocessing import Pool, current_process
from functools import partial
from tqdm import tqdm

# --- 1. Define File Paths ---

# Dictionary to hold all file paths for easy management
PATHS = {
    'nodes': '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/Nodes.csv',
    'tpm': '../../public/Multiomics/Step19_ML_UpDEG_Prediction/GSE175204 ENCFF721BRA hg38 MCF7 total RNAseq ENSG Only Symbol and TPM.csv',
    
    # BigWig files
    'PCD_RFP': '../../public/ChIPseq/SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2_ppois.bigWig',
    'H3K4me1': '../../public/Public Dataset/ChIPseq/MCF7 H3K4me1 ChIPseq GSE86714 ENCFF763NCP hg38.bigWig',
    'H3K4me2': '../../public/Public Dataset/ChIPseq/MCF7 H3K4me2 ChIPseq GSE96439 ENCFF442RRY hg38.bigWig',
    'H3K4me3': '../../public/Public Dataset/ChIPseq/MCF7 H3K4me3 ChIPseq GSE96506 ENCFF163MXP hg38.bigWig',
    'H3K9ac': '../../public/Public Dataset/ChIPseq/MCF7 H3K9ac ChIPseq GSE95898 ENCFF327XJC hg38.bigWig',
    'H3K9me3': '../../public/Public Dataset/ChIPseq/MCF7 H3K9me3 ChIPseq GSE96517 ENCFF481DZL hg38.bigWig',
    'H3K27ac': '../../public/Public Dataset/ChIPseq/MCF7 H3K27ac ChIPseq GSE96352 ENCFF138YNG hg38.bigWig',
    'H3K27me3': '../../public/Public Dataset/ChIPseq/MCF7 H3K27me3 ChIPseq GSE96363 ENCFF163QKN hg38.bigWig',
    'H3K36me3': '../../public/Public Dataset/ChIPseq/MCF7 H3K36me3 ChIPseq GSE174945 ENCFF910BRP hg38.bigWig',
    'H4K20me1': '../../public/Public Dataset/ChIPseq/MCF7 H4K20me1 ChIPseq GSE96283 ENCFF366GLZ hg38.bigWig',
    
    # Output files
    'node_features_out': '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/Node_Features.csv',
    'features_out': '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/Features.csv'
}

# Define the order of BigWig columns (must match the worker function)
BW_COLUMN_NAMES = [
    'PCD_RFP', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9ac', 
    'H3K9me3', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H4K20me1'
]

# --- Multiprocessing Worker Functions ---

# Global variable to hold BigWig file handles *within each worker process*
# This avoids pickling errors, as file handles cannot be passed between processes.
worker_bw_handles = {}

def init_worker():
    """
    Initializer function for each worker process.
    Opens all necessary BigWig files once per process.
    """
    global worker_bw_handles
    pid = current_process().pid
    # print(f"Initializing worker {pid}...")
    for col_name in BW_COLUMN_NAMES:
        path = PATHS[col_name]
        try:
            worker_bw_handles[col_name] = pyBigWig.open(path)
        except Exception as e:
            print(f"Worker {pid} ERROR: Could not open {path}: {e}")
            worker_bw_handles[col_name] = None

def process_row_bw(row_tuple):
    """
    Processes a single row to extract BigWig features.
    `row_tuple` is expected to be (node_type, chr, start, end)
    """
    global worker_bw_handles
    node_type, chrom, start, end = row_tuple
    
    # 6. If node_type is 'gene', return all zeros
    if node_type == 'gene':
        return [0.0] * len(BW_COLUMN_NAMES)

    features = []
    
    # Ensure coordinates are integers
    try:
        # Assuming 1-based, inclusive coordinates from CSV.
        # pyBigWig stats() needs 0-based, exclusive coordinates.
        # (start-1) converts 1-based start to 0-based start.
        # (end) is used as-is, as 1-based inclusive end is 0-based exclusive end.
        query_start = int(start) - 1
        query_end = int(end)
        
        # Handle single-base-pair regions (e.g., gene TSS)
        if query_start == query_end:
            query_end += 1 # pyBigWig needs at least 1bp region
            
        if query_start < 0: # Safety check
             query_start = 0
             
    except ValueError:
        # Handle non-integer coordinates if any
        return [0.0] * len(BW_COLUMN_NAMES)

    # 6. If not 'gene', extract mean values
    for col_name in BW_COLUMN_NAMES:
        bw = worker_bw_handles.get(col_name)
        
        if bw is None:
            features.append(0.0)
            continue
            
        try:
            # Check if chromosome exists in the BigWig file
            if chrom not in bw.chroms():
                features.append(0.0)
                continue
                
            # Get mean value for the region
            mean_val = bw.stats(chrom, query_start, query_end, type='mean')
            
            if mean_val is None or mean_val[0] is None:
                features.append(0.0)
            else:
                features.append(mean_val[0])
                
        except RuntimeError:
            # This can happen if the region is invalid or other pyBigWig error
            features.append(0.0)
        except Exception as e:
            print(f"Error processing {chrom}:{query_start}-{query_end} for {col_name}: {e}")
            features.append(0.0)
            
    return features

# --- Main Execution ---

def main():
    print("Starting Genomic Feature Engineering with Polars...")
    
    # Set number of CPUs
    NUM_CPUS = 25
    print(f"Using {NUM_CPUS} CPU cores for parallel processing.")

    # --- Steps 1 & 2: Import files and create Node_Features table ---
    print("Step 1/2: Loading Nodes.csv and initializing Node_Features table...")
    try:
        nodes_df = pl.read_csv(PATHS['nodes'])
    except FileNotFoundError:
        print(f"ERROR: Cannot find file {PATHS['nodes']}. Exiting.")
        sys.exit(1)
        
    # Create the base Node_Features table
    node_features_df = nodes_df.clone()
    
    print(f"Loaded {len(node_features_df)} nodes.")

    # --- Step 4: One-hot encode the 'node_type' column ---
    print("Step 4: One-hot encoding 'node_type' column...")
    node_features_df = node_features_df.with_columns(
        pl.when(pl.col('node_type') == 'gene').then(1).otherwise(0).alias('is_gene'),
        pl.when(pl.col('node_type') == 'promoter_proximal_downstream').then(1).otherwise(0).alias('is_promoter_proximal_downstream'),
        pl.when(pl.col('node_type') == 'promoter_core').then(1).otherwise(0).alias('is_promoter_core'),
        pl.when(pl.col('node_type') == 'promoter_proximal_upstream').then(1).otherwise(0).alias('is_promoter_proximal_upstream'),
        pl.when(pl.col('node_type') == 'enhancer').then(1).otherwise(0).alias('is_enhancer')
    )

    # --- Step 5: Add TPM data ---
    print("Step 5: Mapping TPM data to 'gene' nodes...")
    try:
        tpm_df = pl.read_csv(PATHS['tpm']).select(['gene_symbol', 'TPM'])
        
        # --- FIX: Deduplicate TPM data before joining ---
        # The source TPM file may have multiple entries per gene.
        # We sort by TPM descending, so the highest value is first.
        # Then, we take the unique 'gene_symbol', keeping the first (highest TPM) entry.
        print(f"Loaded {len(tpm_df)} TPM entries. Deduplicating...")
        tpm_df = (
            tpm_df
            .sort(by=['gene_symbol', 'TPM'], descending=[False, True])
            .unique(subset=['gene_symbol'], keep='first')
        )
        print(f"Using {len(tpm_df)} unique TPM entries for join.")
        # --- End Fix ---

    except FileNotFoundError:
        print(f"ERROR: Cannot find file {PATHS['tpm']}. Exiting.")
        sys.exit(1)
        
    # Join the TPM data onto the node_features table
    # This will no longer create duplicates because tpm_df is now unique by gene_symbol
    node_features_df = node_features_df.join(
        tpm_df,
        left_on='identifier',
        right_on='gene_symbol',
        how='left'
    )
    
    # Apply logic: TPM only for 'gene' nodes, 0 otherwise. Fill nulls from join.
    node_features_df = node_features_df.with_columns(
        pl.when(pl.col('node_type') == 'gene')
          .then(pl.col('TPM'))
          .otherwise(0.0)
          .fill_null(0.0)
          .alias('TPM')
    )

    # --- Step 6: Add BigWig data using multiprocessing ---
    print("Step 6: Extracting BigWig features (this may take a while)...")
    
    # Create an iterable of row tuples for the workers
    # We only need 'node_type', 'chr', 'start', 'end'
    rows_to_process = node_features_df.select(['node_type', 'chr', 'start', 'end']).iter_rows()
    
    # Calculate a reasonable chunk size
    num_rows = len(node_features_df)
    chunksize = max(1, num_rows // (NUM_CPUS * 4))

    bw_results = []
    
    # Create the multiprocessing Pool
    with Pool(processes=NUM_CPUS, initializer=init_worker) as pool:
        # Use pool.imap to get results in order and show progress with tqdm
        # imap is generally more memory-efficient for large iterables
        print(f"Processing {num_rows} regions with chunksize {chunksize}...")
        bw_results = list(tqdm(
            pool.imap(process_row_bw, rows_to_process, chunksize=chunksize),
            total=num_rows,
            desc="Extracting BigWig Signals"
        ))

    print("BigWig processing complete. Integrating results...")
    
    # Convert list of lists into a Polars DataFrame
    bw_results_df = pl.DataFrame(bw_results, schema=BW_COLUMN_NAMES)
    
    # Concatenate the results horizontally
    node_features_df = pl.concat([node_features_df, bw_results_df], how='horizontal')

    # --- Step 3 & 7: Reorder columns and Export Node_Features ---
    print("Step 3/7: Reordering columns and exporting Node_Features.csv...")
    
    # Define the final column order as requested
    final_node_features_cols = [
        'node_id', 'identifier', 'node_type', 'chr', 'start', 'end',
        'is_gene', 'is_promoter_proximal_downstream', 'is_promoter_core',
        'is_promoter_proximal_upstream', 'is_enhancer', 'TPM'
    ] + BW_COLUMN_NAMES
    
    # Reorder the DataFrame
    node_features_df = node_features_df.select(final_node_features_cols)
    
    # Export
    node_features_df.write_csv(PATHS['node_features_out'], float_precision=6)
    print(f"Successfully saved {PATHS['node_features_out']}")

    # --- Step 8 & 9: Create and slim down Features table ---
    print("Step 8/9: Creating and modifying Features table...")
    
    cols_to_drop = ['identifier', 'node_type', 'chr', 'start', 'end']
    features_df = node_features_df.drop(cols_to_drop)

    # --- Step 10: Export Features table ---
    print("Step 10: Exporting Features.csv...")
    features_df.write_csv(PATHS['features_out'], float_precision=6)
    print(f"Successfully saved {PATHS['features_out']}")
    
    print("\nAll steps completed successfully.")

if __name__ == "__main__":
    # This check is crucial for multiprocessing to work correctly on
    # Windows and macOS (spawn start method)
    main()
