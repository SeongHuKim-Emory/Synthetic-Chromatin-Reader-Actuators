import sys
import polars as pl
import pyBigWig
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import time

# --- 1. Define File Paths ---
BASE_PATH = Path("../../public")
STEP25_PATH = BASE_PATH / "Multiomics/Step25_HAN_UpDEG_Prediction"
STEP19_PATH = BASE_PATH / "Multiomics/Step19_ML_UpDEG_Prediction"
CHIPSEQ_PUBLIC_PATH = BASE_PATH / "Public Dataset/ChIPseq"
CHIPSEQ_SHK_PATH = BASE_PATH / "ChIPseq"

# Input CSVs
NODE_FILE = STEP25_PATH / "Nodes_with_gene_enhancer_promoter_pro_en.csv"
ENHANCER_FILE = STEP25_PATH / "Enhancer_Nodes.csv"
PROMOTER_FILE = STEP25_PATH / "Promoter_Nodes.csv"
PRO_EN_FILE = STEP25_PATH / "Pro_En_Nodes.csv"
TPM_FILE = STEP19_PATH / "MCF7_RNAseq_TPM_Preprocessed_GSE175204_ENCFF721BRA_hg38.csv"

# Input BigWigs
PCD_RFP_FILE = CHIPSEQ_SHK_PATH / "SHK_2024_06_11_MCF7_DBN025_1_0_ug_mL_Rep2_ppois.bigWig"
H3K4me1_FILE = CHIPSEQ_PUBLIC_PATH / "MCF7 H3K4me1 ChIPseq GSE86714 ENCFF763NCP hg38.bigWig"
H3K4me2_FILE = CHIPSEQ_PUBLIC_PATH / "MCF7 H3K4me2 ChIPseq GSE96439 ENCFF442RRY hg38.bigWig"
H3K4me3_FILE = CHIPSEQ_PUBLIC_PATH / "MCF7 H3K4me3 ChIPseq GSE96506 ENCFF163MXP hg38.bigWig"
H3K9ac_FILE = CHIPSEQ_PUBLIC_PATH / "MCF7 H3K9ac ChIPseq GSE95898 ENCFF327XJC hg38.bigWig"
H3K9me3_FILE = CHIPSEQ_PUBLIC_PATH / "MCF7 H3K9me3 ChIPseq GSE96517 ENCFF481DZL hg38.bigWig"
H3K27ac_FILE = CHIPSEQ_PUBLIC_PATH / "MCF7 H3K27ac ChIPseq GSE96352 ENCFF138YNG hg38.bigWig"
H3K27me3_FILE = CHIPSEQ_PUBLIC_PATH / "MCF7 H3K27me3 ChIPseq GSE96363 ENCFF163QKN hg38.bigWig"
H3K36me3_FILE = CHIPSEQ_PUBLIC_PATH / "MCF7 H3K36me3 ChIPseq GSE174945 ENCFF910BRP hg38.bigWig"
H4K20me1_FILE = CHIPSEQ_PUBLIC_PATH / "MCF7 H4K20me1 ChIPseq GSE96283 ENCFF366GLZ hg38.bigWig"

# Output Files
NODE_FEATURES_OUT_FILE = STEP25_PATH / "Node_Features.csv"
FEATURES_OUT_FILE = STEP25_PATH / "Features.csv"

# --- Define BigWig Columns and File Map (Global for worker processes) ---
# pyBigWig requires string paths
BIGWIG_FILES = {
    "PCD_RFP": str(PCD_RFP_FILE),
    "H3K4me1": str(H3K4me1_FILE),
    "H3K4me2": str(H3K4me2_FILE),
    "H3K4me3": str(H3K4me3_FILE),
    "H3K9ac": str(H3K9ac_FILE),
    "H3K9me3": str(H3K9me3_FILE),
    "H3K27ac": str(H3K27ac_FILE),
    "H3K27me3": str(H3K27me3_FILE),
    "H3K36me3": str(H3K36me3_FILE),
    "H4K20me1": str(H4K20me1_FILE)
}
BIGWIG_COLS = list(BIGWIG_FILES.keys())


# --- Worker Function for Parallel BigWig Processing ---
def process_bw_row(row_tuple):
    """
    Processes a single node's genomic region against all BigWig files.
    
    Args:
        row_tuple: A tuple containing (node_id, node_type, chr, start, end)
        
    Returns:
        A tuple: (node_id, [list_of_mean_bw_values])
    """
    node_id, node_type, chrom, start, end = row_tuple
    
    # Skip genes, they will have 0.0 for these features
    if node_type == "gene":
        return (node_id, [0.0] * len(BIGWIG_COLS))
        
    results = []
    # Ensure coordinates are integers for pyBigWig
    try:
        int_start = int(start)
        int_end = int(end)
    except ValueError:
        # Handle potential non-integer coordinates if any
        return (node_id, [0.0] * len(BIGWIG_COLS))

    for col_name in BIGWIG_COLS:  # Iterate in the defined order
        bw_path = BIGWIG_FILES[col_name]
        mean_val = 0.0
        try:
            # Open BigWig file *inside* the worker
            with pyBigWig.open(bw_path) as bw:
                # Get mean value for the region
                val = bw.stats(chrom, int_start, int_end, "mean")
                
                # stats() returns [mean] or [None] or None
                if val is not None and val[0] is not None:
                    mean_val = val[0]
                # else: mean_val remains 0.0
        except RuntimeError:
            # Catch errors like "chromosome not found" or invalid region
            mean_val = 0.0
        except Exception as e:
            # Catch other potential errors
            # print(f"Warning: Error processing {node_id} ({chrom}:{start}-{end}) for {col_name}: {e}", file=sys.stderr)
            mean_val = 0.0
            
        results.append(mean_val)
    
    return (node_id, results)


# --- Main Execution ---
def main():
    start_time = time.time()
    print("Script started...")

    # --- Step 1, 2, 3: Load base nodes and initialize Node_Features table ---
    print(f"Loading base nodes from {NODE_FILE}...")
    try:
        node_features_df = pl.read_csv(NODE_FILE, dtypes={
            "node_id": pl.Int64,
            "identifier": pl.Utf8,
            "node_type": pl.Categorical,
            "chr": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64
        })
    except Exception as e:
        print(f"Error loading {NODE_FILE}: {e}", file=sys.stderr)
        return

    print("Step 1, 2, 3: Initializing Node_Features table... done.")

    # --- Step 4: Populate features from other CSVs using joins ---
    print("Step 4: Loading and joining supplementary CSV data...")
    
    # Load supplementary files
    try:
        tpm_df = pl.read_csv(TPM_FILE).select(["gene_symbol", "TPM"]).rename({"TPM": "TPM_new"})
        
        enhancer_df = pl.read_csv(ENHANCER_FILE).select(["identifier", "is_elite", "enhancer_score"])
        promoter_df = pl.read_csv(PROMOTER_FILE).select(["identifier", "is_elite", "enhancer_score"])
        pro_en_df = pl.read_csv(PRO_EN_FILE).select(["identifier", "is_elite", "enhancer_score"])
        
        # Combine all regulatory features into one DataFrame
        regulatory_df = pl.concat([enhancer_df, promoter_df, pro_en_df], how="vertical") \
                          .rename({"is_elite": "is_elite_new", "enhancer_score": "enhancer_score_new"})
                          
    except Exception as e:
        print(f"Error loading supplementary CSVs: {e}", file=sys.stderr)
        return

    # Join features
    node_features_df = node_features_df.join(tpm_df, left_on="identifier", right_on="gene_symbol", how="left")
    node_features_df = node_features_df.join(regulatory_df, on="identifier", how="left")

    # Create final columns using when/then and fill nulls
    # This combines steps 2, 3, and 4 efficiently
    node_features_df = node_features_df.with_columns([
        # One-hot encode node_type
        pl.when(pl.col("node_type") == "gene").then(1).otherwise(0).alias("is_gene"),
        pl.when(pl.col("node_type") == "Enhancer").then(1).otherwise(0).alias("is_Enhancer"),
        pl.when(pl.col("node_type") == "Promoter").then(1).otherwise(0).alias("is_Promoter"),
        pl.when(pl.col("node_type") == "Promoter/Enhancer").then(1).otherwise(0).alias("is_Promoter_Enhancer"),
        
        # Assign TPM: use new value if gene, else 0. Fill nulls (genes not in TPM file) with 0.
        pl.when(pl.col("node_type") == "gene").then(pl.col("TPM_new")).otherwise(0.0).fill_null(0.0).alias("TPM"),
        
        # Assign is_elite: use new value if not gene, else 0. Fill nulls (reg elements not in files) with 0.
        pl.when(pl.col("node_type") != "gene").then(pl.col("is_elite_new")).otherwise(0).fill_null(0).cast(pl.Int64).alias("is_elite"),
        
        # Assign enhancer_score: use new value if not gene, else 0.0. Fill nulls with 0.0.
        pl.when(pl.col("node_type") != "gene").then(pl.col("enhancer_score_new")).otherwise(0.0).fill_null(0.0).alias("enhancer_score"),

        # Initialize BigWig columns (will be populated in step 5)
        pl.lit(0.0).alias("PCD_RFP"),
        pl.lit(0.0).alias("H3K4me1"),
        pl.lit(0.0).alias("H3K4me2"),
        pl.lit(0.0).alias("H3K4me3"),
        pl.lit(0.0).alias("H3K9ac"),
        pl.lit(0.0).alias("H3K9me3"),
        pl.lit(0.0).alias("H3K27ac"),
        pl.lit(0.0).alias("H3K27me3"),
        pl.lit(0.0).alias("H3K36me3"),
        pl.lit(0.0).alias("H4K20me1")
    ])

    # Select final columns in order and drop temporary helper columns
    final_columns = [
        "node_id", "identifier", "node_type", "chr", "start", "end",
        "is_gene", "is_Enhancer", "is_Promoter", "is_Promoter_Enhancer",
        "is_elite", "enhancer_score", "TPM",
        "PCD_RFP", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9ac", "H3K9me3",
        "H3K27ac", "H3K27me3", "H3K36me3", "H4K20me1"
    ]
    node_features_df = node_features_df.select(final_columns)
    
    print("Step 4: Joining supplementary data... done.")

    # --- Step 5: Process BigWig files in parallel ---
    print(f"Step 5: Starting BigWig processing with 36 workers...")
    
    # Prepare data for parallel processing: list of tuples
    rows_to_process = node_features_df.select(["node_id", "node_type", "chr", "start", "end"]).iter_rows()
    total_rows = len(node_features_df)
    
    bw_results_list = []
    
    with ProcessPoolExecutor(max_workers=36) as executor:
        # map() processes the iterator and returns results in the original order
        bw_results_iter = executor.map(process_bw_row, rows_to_process)
        
        # Collect results and show progress
        for i, result in enumerate(bw_results_iter):
            bw_results_list.append(result)
            if (i + 1) % 20000 == 0 or (i + 1) == total_rows:
                print(f"  ... processed {i+1}/{total_rows} rows for BigWig signals", file=sys.stderr)

    print("BigWig processing complete. Merging results...")

    # Convert results list back to a Polars DataFrame
    # Schema: node_id + all the BigWig columns
    bw_schema = [("node_id", pl.Int64)] + [(col, pl.Float64) for col in BIGWIG_COLS]
    # Unpack the results: (node_id, [val1, val2, ...]) -> (node_id, val1, val2, ...)
    bw_data_for_df = [(res[0], *res[1]) for res in bw_results_list]
    
    bw_results_df = pl.DataFrame(data=bw_data_for_df, schema=bw_schema)
    
    # Update the main DataFrame with the new BigWig values
    # This efficiently replaces the 0.0s with the calculated values
    # FIX: Re-assign the result of update() back to the DataFrame
    node_features_df = node_features_df.update(bw_results_df, on="node_id", how="left")
    
    print("Step 5: BigWig processing and merging... done.")

    # --- Step 6: Export Node_Features table ---
    print(f"Step 6: Exporting Node_Features table to {NODE_FEATURES_OUT_FILE}...")
    try:
        node_features_df.write_csv(NODE_FEATURES_OUT_FILE)
    except Exception as e:
        print(f"Error writing {NODE_FEATURES_OUT_FILE}: {e}", file=sys.stderr)
        
    print("Step 6: Export... done.")

    # --- Step 7, 8, 9: Create and Export Features table ---
    print("Step 7, 8: Creating Features table by removing identifier columns...")
    
    cols_to_drop = ["identifier", "node_type", "chr", "start", "end"]
    features_df = node_features_df.drop(cols_to_drop)
    
    print(f"Step 9: Exporting Features table to {FEATURES_OUT_FILE}...")
    try:
        features_df.write_csv(FEATURES_OUT_FILE)
    except Exception as e:
        print(f"Error writing {FEATURES_OUT_FILE}: {e}", file=sys.stderr)
        
    print("Step 7, 8, 9: Export... done.")

    end_time = time.time()
    print(f"\nScript finished in {end_time - start_time:.2f} seconds.")


if __name__ == "__main__":
    # Set Polars to use many threads for its own operations (like joins)
    # This is separate from the ProcessPoolExecutor
    pl.set_random_seed(42)
    
    main()