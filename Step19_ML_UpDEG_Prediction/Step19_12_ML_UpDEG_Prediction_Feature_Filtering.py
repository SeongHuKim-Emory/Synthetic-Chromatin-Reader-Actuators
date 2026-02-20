# Python 3.12.3
# UpDEG vs Non-DEG classifier with 2-Stage RFE + Optuna
# - Loads two CSVs (positive/negative)
# - *** REMOVES SPECIFIED DISTAL COLUMNS FIRST ***
# - *** REVISION: Combines train/validation sets for CV ***
# - Splits each by chromosome into train/validation/test
# - *** QC Step: Checks for NaN values ***
# - *** STAGE 1: WRAPPER (RFE) - Step 1 to select top 1000 features ***
# - *** STAGE 2: Optuna Tuning (finds best XGB params AND best n_features for RFE Step 2) ***
# - *** REVISION: Uses 10-fold Group (Chromosome-based) CV for tuning evaluation ***
# - *** STAGE 3: Runs RFE - Step 2 (with best n_features) & Trains final XGBoost model ***
# - *** REVISION: Trains final model on all data (no eval_set) using best iteration from CV ***
# - *** STAGE 4: Evaluates final model on Test Set (NOW INCLUDES PLOTS) ***
# - Saves model (.bin and .json), metrics, predictions, confusion matrix, feature importance, and plots.

# ---
# This version is optimized for I/O and processing performance using the 'polars' library.
#
# *** CUDA MODIFICATION ***
# This script will automatically use an NVIDIA GPU with CUDA if xgboost
# is installed with CUDA support and a GPU is detected.
# Otherwise, it will fall back to CPU-based training.
#
# You must install it first: pip install polars xgboost scikit-learn pandas optuna scipy matplotlib
# ---

import os
import json
import warnings
from typing import Dict, List, Tuple
import time
import sys
import threading
import pickle

# --- NEW IMPORTS FOR PLOTTING ---
import matplotlib
matplotlib.use('Agg') # Use non-interactive backend for server-side plotting
import matplotlib.pyplot as plt
# --- END NEW IMPORTS ---

import numpy as np
import pandas as pd
import polars as pl  # <-- Import polars for fast I/O and processing
from scipy.sparse import csr_matrix, hstack  # <-- Import hstack for sparse matrix
from sklearn.metrics import (
    log_loss,
    accuracy_score,
    roc_auc_score,
    average_precision_score,
    f1_score,
    confusion_matrix,
    roc_curve,            # <-- ADDED
    precision_recall_curve # <-- ADDED
)
# --- REVISION: Import GroupKFold instead of StratifiedKFold ---
from sklearn.model_selection import GroupKFold
# --- END REVISION ---

from sklearn.feature_selection import RFE, SelectKBest, f_classif
import xgboost as xgb
import optuna

# ----------------------------
# CUDA Availability Check
# ----------------------------

def check_cuda_availability():
    """Checks if XGBoost can utilize CUDA."""
    print("Checking for CUDA availability...")
    try:
        # Try to train a minimal model on GPU
        X = np.array([[1, 2], [3, 4]])
        y = np.array([0, 1])
        # Use 'hist' with 'device' for modern XGBoost
        clf = xgb.XGBClassifier(tree_method="hist", device="cuda")
        clf.fit(X, y)
        print("  CUDA is available! XGBoost will use the GPU.")
        return True
    except Exception as e:
        # Check if the error is related to CUDA
        if "CUDA" in str(e) or "GPU" in str(e) or "cuda" in str(e):
            print("  WARNING: CUDA check failed.")
            print(f"  Error: {e}")
            print("  XGBoost will run on CPU. This may be very slow.")
        else:
            print(f"  CUDA check failed with an unexpected error: {e}")
            print("  Defaulting to CPU.")
        return False

IS_CUDA_AVAILABLE = check_cuda_availability()

if IS_CUDA_AVAILABLE:
    # Parameters for GPU execution
    DEVICE_PARAMS = {
        "tree_method": "hist",
        "device": "cuda"
    }
    # Parameters for GPU execution inside Optuna trials (potentially fewer resources)
    OPTUNA_DEVICE_PARAMS = {
        "tree_method": "hist",
        "device": "cuda"
    }
else:
    # Fallback parameters for CPU execution
    DEVICE_PARAMS = {
        "n_jobs": 25  # Use 25 CPU cores (as in original script)
    }
    # Fallback for Optuna trials
    OPTUNA_DEVICE_PARAMS = {
        "n_jobs": 4   # Use 4 CPU cores (as in original script)
    }

print("-" * 30)


# ----------------------------
# Configuration and constants
# ----------------------------
print("Script started. Initializing configuration...")
BASE_DIR = "../../public/Multiomics/Step19_ML_UpDEG_Prediction"

# Input files for the positive dataset splits
POS_TRAIN_PATH = os.path.join(
    BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone_Training.csv",
)
POS_VALID_PATH = os.path.join(
    BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone_Validating.csv",
)
POS_TEST_PATH = os.path.join(
    BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone_Testing.csv",
)

# Input files for the negative dataset splits
NEG_TRAIN_PATH = os.path.join(
    BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone_Training.csv",
)
NEG_VALID_PATH = os.path.join(
    BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone_Validating.csv",
)
NEG_TEST_PATH = os.path.join(
    BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone_Testing.csv",
)

# Model and metrics output paths
MODEL_BIN_PATH = os.path.join(BASE_DIR, "UpDEG_Prediction_trained_model.bin")
MODEL_JSON_PATH = os.path.join(BASE_DIR, "UpDEG_Prediction_trained_model.json")
# --- REMOVED: train_validation_metrics.csv is no longer generated ---
# TRAIN_VALID_METRICS_PATH = os.path.join(BASE_DIR, "train_validation_metrics.csv")
TEST_PREDICTIONS_PATH = os.path.join(BASE_DIR, "test_predictions.csv")
TEST_METRICS_PATH = os.path.join(BASE_DIR, "test_metrics.csv")
TEST_CONFUSION_PATH = os.path.join(BASE_DIR, "test_confusion_matrix.csv")
FEATURE_IMPORTANCE_PATH = os.path.join(BASE_DIR, "feature_importance.csv")

# --- MODIFIED FILE PATHS ---
RFE_STEP1_FEATURES_PATH = os.path.join(BASE_DIR, "selected_top_1000_features_rfe.json")
SELECTED_FEATURES_PATH = os.path.join(BASE_DIR, "selected_final_features.json")

# New paths for Optuna tuning
BEST_PARAMS_PATH = os.path.join(BASE_DIR, "best_model_params.json")
TUNING_TRIALS_PATH = os.path.join(BASE_DIR, "tuning_trials.csv")
# --- REMOVED OPTUNA_STUDY_PATH ---

# --- NEW FILE PATHS FOR PLOTS ---
AUC_CURVE_PATH = os.path.join(BASE_DIR, "UpDEG_Prediction_AUC_curve.jpg")
AUCPR_CURVE_PATH = os.path.join(BASE_DIR, "UpDEG_Prediction_AUCPR_curve.jpg")
# --- END NEW FILE PATHS ---


# Chromosome splits for datasets
# --- REVISION: TRAIN_CHRS and VALID_CHRS will be combined ---
TRAIN_CHRS = {f"chr{i}" for i in range(1, 14)}
VALID_CHRS = {"chr14", "chr15", "chr16", "chr17"}
# --- TEST_CHRS are still reserved exclusively for testing ---
TEST_CHRS = {"chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}

RANDOM_STATE = 42
N_TRIALS = 100  # Number of Optuna tuning trials
N_SPLITS_CV = 10 # Number of folds for cross-validation
warnings.filterwarnings("ignore", category=UserWarning)

# --- COLUMNS TO REMOVE ---
COLUMNS_TO_REMOVE = [
    "Max_Distal_025", "Mean_Distal_025", "Max_Distal_ARID3A", "Mean_Distal_ARID3A",
    "Max_Distal_ATF7", "Mean_Distal_ATF7", "Max_Distal_BMI1", "Mean_Distal_BMI1",
    "Max_Distal_CEBPB", "Mean_Distal_CEBPB", "Max_Distal_CHD1", "Mean_Distal_CHD1",
    "Max_Distal_CLOCK", "Mean_Distal_CLOCK", "Max_Distal_COPS2", "Mean_Distal_COPS2",
    "Max_Distal_CREB1", "Mean_Distal_CREB1", "Max_Distal_CSDE1", "Mean_Distal_CSDE1",
    "Max_Distal_CTBP1", "Mean_Distal_CTBP1", "Max_Distal_CTCF", "Mean_Distal_CTCF",
    "Max_Distal_CUX1", "Mean_Distal_CUX1", "Max_Distal_DDX20", "Mean_Distal_DDX20",
    "Max_Distal_DPF2", "Mean_Distal_DPF2", "Max_Distal_E2F8", "Mean_Distal_E2F8",
    "Max_Distal_E4F1", "Mean_Distal_E4F1", "Max_Distal_EGR1", "Mean_Distal_EGR1",
    "Max_Distal_ELF1", "Mean_Distal_ELF1", "Max_Distal_ELK1", "Mean_Distal_ELK1",
    "Max_Distal_EP300", "Mean_Distal_EP300", "Max_Distal_ESRRA", "Mean_Distal_ESRRA",
    "Max_Distal_FOS", "Mean_Distal_FOS", "Max_Distal_FOSL2", "Mean_Distal_FOSL2",
    "Max_Distal_FOXA1", "Mean_Distal_FOXA1", "Max_Distal_FOXK2", "Mean_Distal_FOXK2",
    "Max_Distal_FOXM1", "Mean_Distal_FOXM1", "Max_Distal_GABPA", "Mean_Distal_GABPA",
    "Max_Distal_GATA3", "Mean_Distal_GATA3", "Max_Distal_GATAD2B", "Mean_Distal_GATAD2B",
    "Max_Distal_GTF2F1", "Mean_Distal_GTF2F1", "Max_Distal_H3K4me1", "Mean_Distal_H3K4me1",
    "Max_Distal_H3K4me2", "Mean_Distal_H3K4me2", "Max_Distal_H3K4me3", "Mean_Distal_H3K4me3",
    "Max_Distal_H3K9ac", "Mean_Distal_H3K9ac", "Max_Distal_H3K9me3", "Mean_Distal_H3K9me3",
    "Max_Distal_H3K27ac", "Mean_Distal_H3K27ac", "Max_Distal_H3K27me3", "Mean_Distal_H3K27me3",
    "Max_Distal_H3K36me3", "Mean_Distal_H3K36me3", "Max_Distal_H4K20me1", "Mean_Distal_H4K20me1",
    "Max_Distal_HCFC1", "Mean_Distal_HCFC1", "Max_Distal_HDAC2", "Mean_Distal_HDAC2",
    "Max_Distal_HDGF", "Mean_Distal_HDGF", "Max_Distal_HES1", "Mean_Distal_HES1",
    "Max_Distal_HSF1", "Mean_Distal_HSF1", "Max_Distal_JUN", "Mean_Distal_JUN",
    "Max_Distal_JUND", "Mean_Distal_JUND", "Max_Distal_LARP7", "Mean_Distal_LARP7",
    "Max_Distal_MAFK", "Mean_Distal_MAFK", "Max_Distal_MAX", "Mean_Distal_MAX",
    "Max_Distal_MAZ", "Mean_Distal_MAZ", "Max_Distal_MBD2", "Mean_Distal_MBD2",
    "Max_Distal_MLLT1", "Mean_Distal_MLLT1", "Max_Distal_MNT", "Mean_Distal_MNT",
    "Max_Distal_MTA1", "Mean_Distal_MTA1", "Max_Distal_MTA2", "Mean_Distal_MTA2",
    "Max_Distal_MTA3", "Mean_Distal_MTA3", "Max_Distal_NBN", "Mean_Distal_NBN",
    "Max_Distal_NCOA3", "Mean_Distal_NCOA3", "Max_Distal_NEUROD1", "Mean_Distal_NEUROD1",
    "Max_Distal_NFIB", "Mean_Distal_NFIB", "Max_Distal_NFRKB", "Mean_Distal_NFRKB",
    "Max_Distal_NFXL1", "Mean_Distal_NFXL1", "Max_Distal_NONO", "Mean_Distal_NONO",
    "Max_Distal_NR2F2", "Mean_Distal_NR2F2", "Max_Distal_NRF1", "Mean_Distal_NRF1",
    "Max_Distal_PAX8", "Mean_Distal_PAX8", "Max_Distal_PKNOX1", "Mean_Distal_PKNOX1",
    "Max_Distal_PML", "Mean_Distal_PML", "Max_Distal_POLR2A", "Mean_Distal_POLR2A",
    "Max_Distal_PPP1R10", "Mean_Distal_PPP1R10", "Max_Distal_RAD21", "Mean_Distal_RAD21",
    "Max_Distal_RAD51", "Mean_Distal_RAD51", "Max_Distal_RCOR1", "Mean_Distal_RCOR1",
    "Max_Distal_REST", "Mean_Distal_REST", "Max_Distal_RFX1", "Mean_Distal_RFX1",
    "Max_Distal_RFX5", "Mean_Distal_RFX5", "Max_Distal_SIN3A", "Mean_Distal_SIN3A",
    "Max_Distal_SIX4", "Mean_Distal_SIX4", "Max_Distal_SMARCA5", "Mean_Distal_SMARCA5",
    "Max_Distal_SMARCE1", "Mean_Distal_SMARCE1", "Max_Distal_SNIP1", "Mean_Distal_SNIP1",
    "Max_Distal_SP1", "Mean_Distal_SP1", "Max_Distal_SREBF1", "Mean_Distal_SREBF1",
    "Max_Distal_SRF", "Mean_Distal_SRF", "Max_Distal_SUZ12", "Mean_Distal_SUZ12",
    "Max_Distal_TAF1", "Mean_Distal_TAF1", "Max_Distal_TARDBP", "Mean_Distal_TARDBP",
    "Max_Distal_TCF7L2", "Mean_Distal_TCF7L2", "Max_Distal_TCF12", "Mean_Distal_TCF12",
    "Max_Distal_TEAD4", "Mean_Distal_TEAD4", "Max_Distal_TOE1", "Mean_Distal_TOE1",
    "Max_Distal_TRIM22", "Mean_Distal_TRIM22", "Max_Distal_YBX1", "Mean_Distal_YBX1",
    "Max_Distal_ZBTB1", "Mean_Distal_ZBTB1", "Max_Distal_ZBTB7B", "Mean_Distal_ZBTB7B",
    "Max_Distal_ZBTB11", "Mean_Distal_ZBTB11", "Max_Distal_ZBTB33", "Mean_Distal_ZBTB33",
    "Max_Distal_ZBTB40", "Mean_Distal_ZBTB40", "Max_Distal_ZFX", "Mean_Distal_ZFX",
    "Max_Distal_ZHX2", "Mean_Distal_ZHX2", "Max_Distal_ZKSCAN1", "Mean_Distal_ZKSCAN1",
    "Max_Distal_ZNF8", "Mean_Distal_ZNF8", "Max_Distal_ZNF24", "Mean_Distal_ZNF24",
    "Max_Distal_ZNF207", "Mean_Distal_ZNF207", "Max_Distal_ZNF217", "Mean_Distal_ZNF217",
    "Max_Distal_ZNF444", "Mean_Distal_ZNF444", "Max_Distal_ZNF507", "Mean_Distal_ZNF507",
    "Max_Distal_ZNF512B", "Mean_Distal_ZNF512B", "Max_Distal_ZNF574", "Mean_Distal_ZNF574",
    "Max_Distal_ZNF579", "Mean_Distal_ZNF579", "Max_Distal_ZNF592", "Mean_Distal_ZNF592",
    "Max_Distal_ZNF687", "Mean_Distal_ZNF687", "Max_Distal_DNase_seq", "Mean_Distal_DNase_seq",
    "Max_Distal_ATACseq", "Mean_Distal_ATACseq"
]
# --- END COLUMNS TO REMOVE ---

# --- REVISION: This variable will be defined after loading data ---
FEATURE_COLS = []
# --- REVISION: This global variable will hold chromosome info for CV ---
train_groups = None

print("Configuration loaded.")
print(f"  NOTE: {len(COLUMNS_TO_REMOVE)} columns will be removed from all input files.")

# ----------------------------
# Utility functions
# ----------------------------

def split_by_chr_pl(
    df: pl.DataFrame,
) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """
    Splits a Polars DataFrame into training, validation, and test sets
    based on chromosome using efficient Polars filtering.

    --- NOTE: This function is no longer called in the revised script
    but is kept for potential future use or reference. ---
    """
    # Use Polars native filtering, which is much faster than pandas
    train = df.filter(pl.col("chr").is_in(TRAIN_CHRS))
    valid = df.filter(pl.col("chr").is_in(VALID_CHRS))
    test = df.filter(pl.col("chr").is_in(TEST_CHRS))
    return train, valid, test


def compute_metrics(
    y_true: np.ndarray, y_prob: np.ndarray, y_pred: np.ndarray
) -> Dict[str, float]:
    """Computes a dictionary of classification metrics."""
    metrics = {}
    try:
        metrics["logloss"] = float(log_loss(y_true, y_prob))
    except (ValueError, TypeError):
        metrics["logloss"] = float("nan")

    try:
        metrics["error"] = float(1.0 - accuracy_score(y_true, y_pred))
    except (ValueError, TypeError):
        metrics["error"] = float("nan")

    try:
        metrics["auc"] = float(roc_auc_score(y_true, y_prob))
    except (ValueError, TypeError):
        metrics["auc"] = float("nan")

    try:
        metrics["aucpr"] = float(average_precision_score(y_true, y_prob))
    except (ValueError, TypeError):
        metrics["aucpr"] = float("nan")

    try:
        metrics["f1"] = float(f1_score(y_true, y_pred))
    except (ValueError, TypeError):
        metrics["f1"] = float("nan")
    return metrics


def get_feature_importance(
    booster: xgb.Booster, feature_names: List[str]
) -> pd.DataFrame:
    """Extracts and formats feature importance scores from the model."""
    # Set feature names explicitly to ensure scores are found
    booster.feature_names = feature_names

    imp_weight = booster.get_score(importance_type="weight")
    imp_gain = booster.get_score(importance_type="gain")
    imp_cover = booster.get_score(importance_type="cover")

    rows = [
        {
            "feature": feat,
            "weight": float(imp_weight.get(feat, 0.0)),
            "gain": float(imp_gain.get(feat, 0.0)),
            "cover": float(imp_cover.get(feat, 0.0)),
        }
        for feat in feature_names
    ]
    return pd.DataFrame(rows)


# ----------------------------
# Main execution
# ----------------------------
print(f"Ensuring base directory exists: {BASE_DIR}")
os.makedirs(BASE_DIR, exist_ok=True)


# 1. Load pre-split files (Using polars for speed) and REMOVE SPECIFIED COLUMNS
print("\n" + "=" * 50)
print("STEP 1: Loading pre-split datasets and removing specified columns...")
print("=" * 50)
start_load_time = time.time()

try:
    print(f"  Loading {POS_TRAIN_PATH}...")
    pos_train_pl = pl.read_csv(POS_TRAIN_PATH)
    # Remove columns that exist in the dataframe
    cols_to_drop = [c for c in COLUMNS_TO_REMOVE if c in pos_train_pl.columns]
    if cols_to_drop:
        pos_train_pl = pos_train_pl.drop(cols_to_drop)
        print(f"    Removed {len(cols_to_drop)} columns from POS_TRAIN")

    print(f"  Loading {POS_VALID_PATH}...")
    pos_valid_pl = pl.read_csv(POS_VALID_PATH)
    cols_to_drop = [c for c in COLUMNS_TO_REMOVE if c in pos_valid_pl.columns]
    if cols_to_drop:
        pos_valid_pl = pos_valid_pl.drop(cols_to_drop)
        print(f"    Removed {len(cols_to_drop)} columns from POS_VALID")

    print(f"  Loading {POS_TEST_PATH}...")
    pos_test_pl = pl.read_csv(POS_TEST_PATH)
    cols_to_drop = [c for c in COLUMNS_TO_REMOVE if c in pos_test_pl.columns]
    if cols_to_drop:
        pos_test_pl = pos_test_pl.drop(cols_to_drop)
        print(f"    Removed {len(cols_to_drop)} columns from POS_TEST")

    print(f"  Loading {NEG_TRAIN_PATH}...")
    neg_train_pl = pl.read_csv(NEG_TRAIN_PATH)
    cols_to_drop = [c for c in COLUMNS_TO_REMOVE if c in neg_train_pl.columns]
    if cols_to_drop:
        neg_train_pl = neg_train_pl.drop(cols_to_drop)
        print(f"    Removed {len(cols_to_drop)} columns from NEG_TRAIN")

    print(f"  Loading {NEG_VALID_PATH}...")
    neg_valid_pl = pl.read_csv(NEG_VALID_PATH)
    cols_to_drop = [c for c in COLUMNS_TO_REMOVE if c in neg_valid_pl.columns]
    if cols_to_drop:
        neg_valid_pl = neg_valid_pl.drop(cols_to_drop)
        print(f"    Removed {len(cols_to_drop)} columns from NEG_VALID")

    print(f"  Loading {NEG_TEST_PATH}...")
    neg_test_pl = pl.read_csv(NEG_TEST_PATH)
    cols_to_drop = [c for c in COLUMNS_TO_REMOVE if c in neg_test_pl.columns]
    if cols_to_drop:
        neg_test_pl = neg_test_pl.drop(cols_to_drop)
        print(f"    Removed {len(cols_to_drop)} columns from NEG_TEST")

except pl.exceptions.ComputeError as e:
    print(f"\n--- ERROR: Could not read one of the split files. ---")
    print(f"  File error: {e}")
    print("  Please ensure all six split CSV files exist in the correct directory.")
    exit(1)  # Exit the script if files are missing

print(f"\nAll data splits loaded and columns removed in {time.time() - start_load_time:.2f} seconds.")
print(
    f"  Positive - Train: {pos_train_pl.shape[0]} rows x {pos_train_pl.shape[1]} cols"
)
print(
    f"  Positive - Valid: {pos_valid_pl.shape[0]} rows x {pos_valid_pl.shape[1]} cols"
)
print(
    f"  Positive - Test: {pos_test_pl.shape[0]} rows x {pos_test_pl.shape[1]} cols"
)
print(
    f"  Negative - Train: {neg_train_pl.shape[0]} rows x {neg_train_pl.shape[1]} cols"
)
print(
    f"  Negative - Valid: {neg_valid_pl.shape[0]} rows x {pos_valid_pl.shape[1]} cols"
)
print(
    f"  Negative - Test: {neg_test_pl.shape[0]} rows x {neg_test_pl.shape[1]} cols"
)
print("=" * 50 + "\n")

# -----------------------------------------------------------------
# --- NEW: QC Step for NaN values
# -----------------------------------------------------------------
print("\n" + "=" * 50)
print("STEP 2: QC - Checking for NaN values...")
print("=" * 50)

# Store loaded dataframes and their names for iteration
dataframes_to_check = {
    "POS_TRAIN": pos_train_pl,
    "POS_VALID": pos_valid_pl,
    "POS_TEST": pos_test_pl,
    "NEG_TRAIN": neg_train_pl,
    "NEG_VALID": neg_valid_pl,
    "NEG_TEST": neg_test_pl,
}

total_nans = 0
for name, df in dataframes_to_check.items():
    # Polars `.null_count()` returns a DataFrame of counts.
    # We must sum the counts horizontally to get a single total.
    nan_count = df.null_count().select(pl.sum_horizontal(pl.all())).item()
    if nan_count > 0:
        print(f"  WARNING: Found {nan_count} NaN values in {name} file.")
    total_nans += nan_count

if total_nans == 0:
    print("  QC PASSED: No NaN values found in any input files.")
else:
    print(f"  QC COMPLETE: Found a total of {total_nans} NaN values.")
    print("  (NaN values will be filled with 0 during feature selection/training)")

print("=" * 50 + "\n")
# -----------------------------------------------------------------
# --- END: NEW QC Step
# -----------------------------------------------------------------


# 2. Dynamically define feature columns
print("Defining feature columns...")
# Get columns from the loaded training set.
# This assumes the 'label' column is the LAST column, and features start at index 21,
# matching the logic of the original script.
FEATURE_COLS = pos_train_pl.columns[21:-1]
if not FEATURE_COLS:
    print("--- WARNING ---")
    print("No feature columns were detected based on the logic `columns[21:-1]`.")
    print("Please check your input CSV files. The script may fail.")
print(
    f"Dynamically defined {len(FEATURE_COLS)} feature columns (from index 21 to -1)."
)


# Prepare combined datasets for machine learning
print("\nPreparing combined datasets for training (using Polars)...")
# Concatenate in Polars, which is very fast
# Use how="vertical_relaxed" to automatically handle schema mismatches
# by casting to a common supertype (e.g., Int64 -> Float64)
train_pl_orig = pl.concat([pos_train_pl, neg_train_pl], how="vertical_relaxed")
valid_pl = pl.concat([pos_valid_pl, neg_valid_pl], how="vertical_relaxed")
test_pl = pl.concat([pos_test_pl, neg_test_pl], how="vertical_relaxed")

# --- REVISION: Combine Train and Valid sets for CV ---
print("REVISION: Combining Train and Validation sets for 10-fold CV...")
train_pl = pl.concat([train_pl_orig, valid_pl], how="vertical_relaxed")
# --- `valid_pl` is no longer needed as a separate entity ---
del valid_pl
del train_pl_orig
print(f"  New combined training set size: {train_pl.shape[0]} rows")
# --- END REVISION ---

# --- REVISION: Extract chromosome groups for GroupKFold ---
# We do this here so the global `train_groups` variable is available
# for the `objective` function.
print("Extracting chromosome groups for GroupKFold...")
try:
    # This `train_groups` variable will be accessed by the `objective` function
    train_groups = train_pl.select("chr").to_series()
    print(f"  Found {len(train_groups.unique())} unique chromosome groups for CV.")
except pl.exceptions.ColumnNotFoundError:
    print("\n--- ERROR: 'chr' column not found in training data. ---")
    print("  Cannot perform chromosome-based GroupKFold.")
    print("  Please ensure the 'chr' column exists in your input CSVs.")
    exit(1)
# --- END REVISION ---


# Prepare data for XGBoost
# XGBoost can accept Polars DataFrames/Series directly!
# This avoids a massive .to_pandas() conversion.
print("Creating DMatrix components...")
X_train = train_pl.select(FEATURE_COLS)
y_train = train_pl.select("label").to_series()

# --- REVISION: X_valid/y_valid are removed ---
# X_valid = valid_pl.select(FEATURE_COLS)
# y_valid = valid_pl.select("label").to_series()
# --- END REVISION ---

X_test = test_pl.select(FEATURE_COLS)
y_test = test_pl.select("label").to_series()

# We only need numpy arrays for scikit-learn metrics and scale_pos_weight calc
print("Converting labels to NumPy for metrics...")
y_train_np = y_train.to_numpy().astype(int)
# --- REVISION: y_valid_np removed ---
# y_valid_np = y_valid.to_numpy().astype(int)
y_test_np = y_test.to_numpy().astype(int)

print("ML datasets prepared.")
print(f"  Total Train samples (Combined): {X_train.shape[0]}")
# --- REVISION: Valid samples print removed ---
# print(f"  Total Valid samples: {X_valid.shape[0]}")
print(f"  Total Test samples:  {X_test.shape[0]}")

# 6. Calculate scale_pos_weight for imbalanced data
print("\nCalculating scale_pos_weight for class imbalance (on combined train set)...")
# Use the numpy array for this calculation (now the combined set)
num_pos = int((y_train_np == 1).sum())
num_neg = int((y_train_np == 0).sum())
if num_pos == 0:
    print("ERROR: Training set contains no positive examples.")
    raise ValueError("Training set contains no positive examples.")
scale_pos_weight = float(num_neg) / float(num_pos)
print(f"  Negative samples (train): {num_neg}")
print(f"  Positive samples (train): {num_pos}")
print(f"  Calculated scale_pos_weight: {scale_pos_weight:.4f}")


# -----------------------------------------------------------------
# --- STAGE 1: WRAPPER (RFE) - Step 1 (to 1000 features)
# -----------------------------------------------------------------
print("\n" + "=" * 50)
print("STAGE 1: WRAPPER (RFE) - Step 1 (to 1000 features)")
print(f"  (Starting from {len(FEATURE_COLS)} features)")
print("=" * 50)

# 1. Define a base estimator for RFE.
print("  1. Defining base XGBoost estimator for RFE...")
if IS_CUDA_AVAILABLE:
    print("     (Using GPU-accelerated estimator)")
else:
    print("     (Using CPU-based estimator)")
    
rfe_estimator = xgb.XGBClassifier(
    n_estimators=100,  # Fewer estimators for faster RFE iterations
    max_depth=5,
    learning_rate=0.1,  # Standard learning rate
    colsample_bytree=0.8,
    reg_lambda=1.0,
    random_state=RANDOM_STATE,
    objective="binary:logistic",
    scale_pos_weight=scale_pos_weight,
    **DEVICE_PARAMS # <-- CUDA MODIFICATION: Unpack GPU (or CPU) params
)

# 2. Initialize RFE to select 1000 features
print("  2. Initializing RFE to select 1000 features...")
rfe_step1 = RFE(
    estimator=rfe_estimator,
    n_features_to_select=1000,  # <-- MODIFICATION: Hard-code to 1000
    step=0.1,  # Remove 10% of features per iteration
    verbose=1,
)

# 3. Fit RFE to the training data.
print(f"  3. Fitting RFE on {X_train.shape[0]} samples (combined set) and {X_train.shape[1]} features...")
print("       (This may take a long time with all features)")
start_rfe_time = time.time()

# We fit RFE on the Polars training data (combined set) directly.
# Filling with 0 ensures consistency and prevents errors.
rfe_step1.fit(X_train.fill_null(0), y_train)

print(f"  RFE fitting complete in {time.time() - start_rfe_time:.2f} seconds.")

# 4. Get the list of selected features
print("  4. Extracting list of 1000 features...")
selected_mask_step1 = rfe_step1.support_
STEP1_FEATURE_COLS = [
    col for col, selected in zip(FEATURE_COLS, selected_mask_step1) if selected
]

# 5. Save the selected feature list
print(f"  5. Saving Step 1 feature list ({len(STEP1_FEATURE_COLS)} features) to {RFE_STEP1_FEATURES_PATH}...")
with open(RFE_STEP1_FEATURES_PATH, "w") as f:
    json.dump(STEP1_FEATURE_COLS, f, indent=4)

# 6. Filter all datasets (X_train, X_test) to these features
print("  6. Filtering Train (Combined) and Test sets to Step 1 features...")
# CRITICAL: Overwrite the global FEATURE_COLS list
FEATURE_COLS = STEP1_FEATURE_COLS

X_train = X_train.select(FEATURE_COLS)
# --- REVISION: X_valid removed ---
# X_valid = X_valid.select(FEATURE_COLS)
X_test = X_test.select(FEATURE_COLS)

print("      Datasets successfully filtered.")
print("=" * 50)
print("COMPLETED: STAGE 1: WRAPPER (RFE) - Step 1")
print("=" * 50 + "\n")
# -----------------------------------------------------------------
# --- END: STAGE 1
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# --- STAGE 2: HYPERPARAMETER TUNING (incl. RFE Step 2)
# -----------------------------------------------------------------
print("\n" + "=" * 50)
print(f"STAGE 2: Hyperparameter Tuning ({N_TRIALS} trials)")
print("=" * 50)
print(f"  (Optimizing XGBoost params AND n_features from {len(FEATURE_COLS)} features)")
print(f"  (Using {N_SPLITS_CV}-Fold Group (Chromosome-based) Cross-Validation)")


def objective(trial: optuna.Trial) -> float:
    """Optuna objective function to maximize mean validation AUCPR from CV."""
    
    # Access the global 'train_groups' variable defined outside this function
    global train_groups
    if train_groups is None:
        raise ValueError("train_groups is not set. CV cannot proceed.")

    # --- Define search space ---
    # 1. NEW: Number of features to select in RFE Step 2
    n_features = trial.suggest_int("n_features_to_select", 100, 500, step=50)

    # 2. scale_pos_weight: Center search around the calculated value
    scale_factor = trial.suggest_float("scale_pos_weight_factor", 0.2, 0.3)

    # 3. Core tree parameters
    params = {
        "objective": "binary:logistic",
        "eval_metric": "aucpr", # We optimize for AUCPR
        "scale_pos_weight": scale_pos_weight * scale_factor,
        "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.3, log=True),
        "max_depth": trial.suggest_int("max_depth", 1, 5),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.1, 1.0),
        "reg_lambda": trial.suggest_float("reg_lambda", 1.0, 100.0, log=True),
        "random_state": RANDOM_STATE,
        **OPTUNA_DEVICE_PARAMS # <-- CUDA MODIFICATION: Unpack GPU (or CPU) params
    }

    # --- RFE Step 2 (inside trial) ---
    # Use a lighter estimator for speed inside the trial
    rfe_step2_estimator = xgb.XGBClassifier(
        n_estimators=50, # Lighter model for RFE inside Optuna
        max_depth=5,
        learning_rate=0.1,
        random_state=RANDOM_STATE,
        scale_pos_weight=scale_pos_weight,
        **OPTUNA_DEVICE_PARAMS # <-- CUDA MODIFICATION: Unpack GPU (or CPU) params
    )

    rfe_step2 = RFE(
        estimator=rfe_step2_estimator,
        n_features_to_select=n_features,
        step=0.2, # Larger step for faster RFE
        verbose=0, # Must be silent
    )

    # Fit RFE on the 1000-feature *combined* training data
    rfe_step2.fit(X_train.fill_null(0), y_train)

    # Filter data down to the features selected in *this trial*
    selected_mask_step2 = rfe_step2.support_
    trial_feature_cols = [
        col for col, selected in zip(FEATURE_COLS, selected_mask_step2) if selected
    ]

    X_train_trial = X_train.select(trial_feature_cols)

    # --- REVISION: Model Training (inside trial) using 10-fold Group CV ---
    
    # 1. Create the DMatrix with the RFE-selected features
    dtrain = xgb.DMatrix(X_train_trial.fill_null(0), label=y_train)
    
    # 2. Create the GroupKFold generator
    #    Note: GroupKFold does not shuffle; it just ensures groups are not split.
    #    The "randomness" comes from which chromosomes are in the original data.
    gkf = GroupKFold(n_splits=N_SPLITS_CV)
    #    We pass the full y_train and train_groups, as X_train_trial
    #    only had its columns changed, not its rows.
    folds = list(gkf.split(X_train_trial, y_train, groups=train_groups))
    
    # 3. Run xgb.cv
    try:
        cv_results_df = xgb.cv(
            params=params,          # <-- This now contains GPU/CPU params
            dtrain=dtrain,
            num_boost_round=2000,  # High number, will be stopped by early_stopping
            folds=folds,            # <-- Pass the chromosome-grouped folds
            metrics=["aucpr"],
            early_stopping_rounds=50,
            seed=RANDOM_STATE,
            verbose_eval=False, # Silence the CV output
        )
        
        # 4. Get the best mean AUCPR from the last iteration
        best_aucpr = cv_results_df["test-aucpr-mean"].iloc[-1]
        
        # 5. Save the best number of iterations (boosting rounds)
        best_iteration = len(cv_results_df)
        trial.set_user_attr("best_iteration", best_iteration)

        return float(best_aucpr)

    except Exception as e:
        print(f"  WARNING: Trial {trial.number} failed with error: {e}")
        # Prune the trial if CV fails (e.g., no positive samples in a fold)
        raise optuna.exceptions.TrialPruned()


# 1. Create new Optuna study
# --- MODIFICATION: Always create a new study to reset trials ---
print(f"  1. Creating new Optuna study to maximize mean {N_SPLITS_CV}-fold Group CV AUCPR...")
study = optuna.create_study(
    direction="maximize", 
    sampler=optuna.samplers.TPESampler(seed=RANDOM_STATE),
    pruner=optuna.pruners.MedianPruner(n_warmup_steps=5),
)

# 2. Run the optimization
try:
    print(f"  2. Running {N_TRIALS} tuning trials...")
    study.optimize(
        objective, 
        n_trials=N_TRIALS, 
        n_jobs=4,  # <-- Run fewer trials in parallel (e.g., 4) since each trial now does RFE + CV
        show_progress_bar=True,
        # --- MODIFICATION: Removed callback, no need to save resume-able study ---
    )
    print(f"  Tuning complete. Best Mean CV AUCPR: {study.best_value:.4f}")

except KeyboardInterrupt:
    print("\n  Tuning interrupted. Saving best results so far...")
    pass

# 3. Save the best parameters and all trial results
print(f"  3. Saving best parameters to {BEST_PARAMS_PATH}...")

# --- MODIFICATION: Use .pop() with a default value to prevent KeyError ---
# This handles resuming a study where the best trial is from a previous run
# that did not include these parameters.
best_factor = study.best_params.pop("scale_pos_weight_factor", 1.0)
best_n_features = study.best_params.pop("n_features_to_select", 100) # Default to 100
# --- END MODIFICATION ---

# --- REVISION: Get best_iteration from the trial user_attrs ---
try:
    best_iteration = study.best_trial.user_attrs["best_iteration"]
except (KeyError, AttributeError):
    print("  WARNING: Could not find 'best_iteration' in best trial. Defaulting to 100.")
    best_iteration = 100
# --- END REVISION ---


final_params = study.best_params # 'final_params' now only has XGBoost params
final_params["scale_pos_weight"] = scale_pos_weight * best_factor

# --- CUDA MODIFICATION ---
# Add GPU (or CPU fallback) parameters for the final model
final_params.update(DEVICE_PARAMS)
# --- END MODIFICATION ---

print(f"  Best n_features_to_select found (or defaulted): {best_n_features}")
print(f"  Best n_estimators (boosting rounds) found: {best_iteration}")
print(f"  Best XGBoost params found: {final_params}")

with open(BEST_PARAMS_PATH, "w") as f:
    json.dump(final_params, f, indent=4) # Save all final params

print(f"  4. Saving all trial results to {TUNING_TRIALS_PATH}...")
trials_df_pd = study.trials_dataframe()
trials_pl_df = pl.from_pandas(trials_df_pd)
if 'duration' in trials_pl_df.columns:
    if trials_pl_df['duration'].dtype == pl.Duration:
        trials_pl_df = trials_pl_df.with_columns(
            pl.col('duration').cast(pl.Int64).cast(pl.Utf8)
        )
trials_pl_df.write_csv(TUNING_TRIALS_PATH)

print("=" * 50)
print("COMPLETED: STAGE 2: Hyperparameter Tuning")
print("=" * 50 + "\n")
# -----------------------------------------------------------------
# --- END: STAGE 2
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# --- STAGE 3: RFE - Step 2 & Final Model Training
# -----------------------------------------------------------------
print("\n" + "=" * 50)
print("STAGE 3: RFE - Step 2 & Final Model Training")
print("=" * 50)

# 1. Run RFE Step 2 ONE MORE TIME using the best n_features found
print(f"  1. Running RFE - Step 2 to select best {best_n_features} features...")
# We can use the more robust 'rfe_estimator' from Stage 1 here
# (rfe_estimator already has GPU/CPU params set from DEVICE_PARAMS)
rfe_final = RFE(
    estimator=rfe_estimator, # Use the robust, device-aware estimator
    n_features_to_select=best_n_features,
    step=0.1, # Use a smaller step for the final run
    verbose=1,
)

# Fit on the 1000-feature *combined* training data
rfe_final.fit(X_train.fill_null(0), y_train)

# Get the FINAL list of features
final_mask = rfe_final.support_
FINAL_FEATURE_COLS = [
    col for col, selected in zip(FEATURE_COLS, final_mask) if selected
]

# 2. Save the FINAL feature list
print(f"  2. Saving final feature list ({len(FINAL_FEATURE_COLS)} features) to {SELECTED_FEATURES_PATH}...")
with open(SELECTED_FEATURES_PATH, "w") as f:
    json.dump(FINAL_FEATURE_COLS, f, indent=4)

# 3. Filter all datasets to the FINAL feature list
print("  3. Filtering Train (Combined) and Test sets to final features...")
# CRITICAL: Overwrite global FEATURE_COLS and dataframes
FEATURE_COLS = FINAL_FEATURE_COLS
X_train = X_train.select(FEATURE_COLS)
# --- REVISION: X_valid removed ---
# X_valid = X_valid.select(FEATURE_COLS)
X_test = X_test.select(FEATURE_COLS)
print(f"  Datasets filtered down to {len(FEATURE_COLS)} features.")

# 4. Train the FINAL model
print("\n  4. Initializing FINAL XGBoost classifier...")
print(f"       Using best parameters from tuning: {final_params}")
print(f"       Training for {best_iteration} boosting rounds (from CV).")

# --- REVISION: Use best_iteration, remove early_stopping ---
model = xgb.XGBClassifier(
    n_estimators=best_iteration,  # <-- Use best rounds from CV
    random_state=RANDOM_STATE,
    objective="binary:logistic",
    eval_metric=["logloss", "auc", "aucpr"], # Still useful for .predict_proba
    early_stopping_rounds=None, # <-- No early stopping
    # --- CUDA MODIFICATION ---
    # n_jobs=25, # REMOVED: Replaced by DEVICE_PARAMS inside final_params
    **final_params # Unpack the tuned XGBoost parameters (now includes device or n_jobs)
    # --- END MODIFICATION ---
)
print("Final model initialized.")

print("Starting FINAL model training (on combined data)...")
start_time = time.time()
# --- REVISION: Fit on combined data, no eval_set ---
model.fit(
    X_train.fill_null(0),  # <-- Pass FINAL filtered DataFrame (combined)
    y_train,               # <-- Pass combined labels
    # eval_set=[(X_train.fill_null(0), y_train), (X_valid.fill_null(0), y_valid)], # <-- REMOVED
    verbose=False,
)
print(f"Model training complete in {time.time() - start_time:.2f} seconds.")
# --- END REVISION ---

# --- REVISION: Remove eval results print block ---
# print(f"  Best iteration: {model.best_iteration}") # No longer relevant
# all_results = model.evals_result() # No longer exists
# print("  Best Validation Scores:") # No longer exists
# --- END REVISION ---


# 5. Export the trained model object
print(f"\n  5. Saving trained model (binary) to {MODEL_BIN_PATH}...")
model.save_model(MODEL_BIN_PATH)
print(f"       Saving trained model (JSON) to {MODEL_JSON_PATH}...")
model.save_model(MODEL_JSON_PATH)
print("Model saved.")

# --- REVISION: Remove section 6 (Export training/validation metrics) ---
# This data is no longer generated by model.fit()
print("  (Skipping train/validation metrics history export as CV was used)")
# --- END REVISION ---

# -----------------------------------------------------------------
# --- END: STAGE 3
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# --- STAGE 4: Final Model Evaluation on Test Set
# -----------------------------------------------------------------
print("\n" + "=" * 50)
print("STAGE 4: Final Model Evaluation on Test Set")
print("=" * 50)
print("\nEvaluating model on the test set...")
# We use the 'model' object trained in Stage 3 and 'X_test' filtered in Stage 3
# model.predict_proba will use the GPU if the model was trained on it.
test_prob = model.predict_proba(X_test.fill_null(0))[:, 1]
test_pred = (test_prob >= 0.5).astype(int)
print("Test set evaluation complete.")

# Export predictions
print(f"Saving test predictions to {TEST_PREDICTIONS_PATH}...")
pred_out = pd.DataFrame(
    {
        "pred_proba": test_prob,
        "pred_label": test_pred,
        "true_label": y_test_np,
    }
)
if "symbol" in test_pl.columns:
    pred_out["symbol"] = test_pl.select("symbol").to_series().to_numpy()
if "chr" in test_pl.columns:
    pred_out["chr"] = test_pl.select("chr").to_series().to_numpy()
if "symbol" in pred_out.columns and "chr" in pred_out.columns:
    pred_out = pred_out[
        ["symbol", "chr", "true_label", "pred_label", "pred_proba"]
    ]
pl.DataFrame(pred_out).write_csv(TEST_PREDICTIONS_PATH)
print("Test predictions saved.")

# Export test metrics
print(f"Saving test metrics to {TEST_METRICS_PATH}...")
test_metrics = compute_metrics(y_test_np, test_prob, test_pred)
test_metrics_df = pd.DataFrame([test_metrics])
print("  Test Metrics:")
print(test_metrics_df.to_string())
pl.DataFrame(test_metrics_df).write_csv(TEST_METRICS_PATH)
print("Test metrics saved.")

# Export confusion matrix
print(f"Saving confusion matrix to {TEST_CONFUSION_PATH}...")
tn, fp, fn, tp = confusion_matrix(y_test_np, test_pred).ravel()
confusion_df = pd.DataFrame(
    [{"TP": int(tp), "TN": int(tn), "FP": int(fp), "FN": int(fn)}]
)
print("  Confusion Matrix:")
print(confusion_df.to_string())
pl.DataFrame(confusion_df).write_csv(TEST_CONFUSION_PATH)
print("Confusion matrix saved.")

# --- NEW: Generate and save plots ---
print("\nGenerating and saving evaluation plots...")

try:
    # 1. AUC Curve
    print(f"  Generating AUC curve and saving to {AUC_CURVE_PATH}...")
    fpr, tpr, _ = roc_curve(y_test_np, test_prob)
    auc_score = test_metrics.get("auc", 0.0) # Get score from metrics dict

    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {auc_score:.4f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (AUC) Curve')
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.savefig(AUC_CURVE_PATH, dpi=300, bbox_inches='tight')
    plt.close() # Close the figure to free memory
    print("  AUC curve saved.")

    # 2. AUCPR Curve
    print(f"  Generating AUCPR curve and saving to {AUCPR_CURVE_PATH}...")
    precision, recall, _ = precision_recall_curve(y_test_np, test_prob)
    aucpr_score = test_metrics.get("aucpr", 0.0) # Get score from metrics dict

    plt.figure(figsize=(8, 6))
    plt.plot(recall, precision, color='blue', lw=2, label=f'PR curve (area = {aucpr_score:.4f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('Precision-Recall (AUCPR) Curve')
    plt.legend(loc="lower left")
    plt.grid(True)
    plt.savefig(AUCPR_CURVE_PATH, dpi=300, bbox_inches='tight')
    plt.close() # Close the figure
    print("  AUCPR curve saved.")

except Exception as e:
    print(f"\n--- WARNING: Could not generate plots. ---")
    print(f"  Error: {e}")
    print("  Continuing without plots.")
# --- END: NEW PLOTS ---


# Export feature importance
print(f"\nSaving feature importance to {FEATURE_IMPORTANCE_PATH}...")
booster = model.get_booster()
# FEATURE_COLS was updated in Stage 3 to the final list
booster.feature_names = FEATURE_COLS 
fi_df = get_feature_importance(booster, FEATURE_COLS)
pl.DataFrame(fi_df).write_csv(FEATURE_IMPORTANCE_PATH)
print("Feature importance saved.")

print("\n" + "=" * 50)
print("Script completed successfully!")
print("=" * 50)
