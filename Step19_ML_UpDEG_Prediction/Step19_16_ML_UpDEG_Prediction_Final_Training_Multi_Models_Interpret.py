# Python 3.12.3
# UpDEG vs Non-DEG classifier with 2-Stage RFE + Optuna
# - Loads two CSVs (positive/negative)
# - *** REVISION: NO LONGER REMOVES DISTAL COLUMNS MANUALLY ***
# - *** REVISION: Combines train/validation sets for CV ***
# - Splits each by chromosome into train/validation/test
# - *** QC Step: Checks for NaN values ***
# - *** REVISION: STAGE 1 SKIPPED - Loads pre-selected 1000 features from JSON ***
# - *** REVISION: STAGE 2 SKIPPED - Loads pre-tuned hyperparameters from JSON ***
# - *** REVISION: STAGE 3 Runs RFE - Step 2 (with best n_features) ***
# - *** REVISION: STAGE 3 RFE estimator matched to original tuning RFE to ensure reproducibility ***
# - *** REVISION: STAGE 3 Uses 10-fold Group CV to find best iteration ***
# - *** REVISION: STAGE 3 Trains final model on all data using best iteration ***
# - *** STAGE 4: Evaluates final model on Test Set (NOW INCLUDES PLOTS) ***
# - *** NEW: STAGE 5: SHAP Analysis for Model Interpretability ***
# - Saves model (.bin and .json), metrics, predictions, confusion matrix, feature importance, plots, and SHAP values.
# - *** NEW: Trains LightGBM, CatBoost, and Random Forest for soft-voting ensemble (averaged probabilities) ***

# ---
# This version is optimized for I/O and processing performance using the 'polars' library.
#
# You must install it first: pip install polars xgboost scikit-learn pandas optuna scipy matplotlib shap
# For ensembling extras: pip install lightgbm catboost
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
matplotlib.use('Agg')  # Use non-interactive backend for server-side plotting
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
    roc_curve,
    precision_recall_curve
)
# --- REVISION: Import GroupKFold instead of StratifiedKFold ---
from sklearn.model_selection import GroupKFold
# --- END REVISION ---

from sklearn.feature_selection import RFE, SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier  # <-- NEW: Random Forest
import xgboost as xgb
import optuna  # kept for compatibility though not used for tuning here

# --- NEW: Import SHAP for model interpretability ---
try:
    import shap
except ImportError:
    shap = None
    print("WARNING: SHAP not installed. Install with: pip install shap")
# --- END NEW ---

# --- NEW: Optional imports for diverse learners (handled gracefully if missing) ---
try:
    import lightgbm as lgb
except Exception:
    lgb = None
try:
    from catboost import CatBoostClassifier
except Exception:
    CatBoostClassifier = None
# --- END NEW ---

# ----------------------------
# Configuration and constants
# ----------------------------
print("Script started. Initializing configuration...")

# --- REVISION: Split Base Dir into Input and Output ---
INPUT_BASE_DIR = "../../public/Multiomics/Step19_ML_UpDEG_Prediction"
OUTPUT_BASE_DIR = "../../public/Multiomics/Step19_16_ML_UpDEG_Prediction_Final_Training_Multi_Models_Interpret"
print(f"  Input directory: {INPUT_BASE_DIR}")
print(f"  Output directory: {OUTPUT_BASE_DIR}")
# --- END REVISION ---


# Input files for the positive dataset splits (from INPUT_BASE_DIR)
POS_TRAIN_PATH = os.path.join(
    INPUT_BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone_Training.csv",
)
POS_VALID_PATH = os.path.join(
    INPUT_BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone_Validating.csv",
)
POS_TEST_PATH = os.path.join(
    INPUT_BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone_Testing.csv",
)

# Input files for the negative dataset splits (from INPUT_BASE_DIR)
NEG_TRAIN_PATH = os.path.join(
    INPUT_BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone_Training.csv",
)
NEG_VALID_PATH = os.path.join(
    INPUT_BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone_Validating.csv",
)
NEG_TEST_PATH = os.path.join(
    INPUT_BASE_DIR,
    "GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone_Testing.csv",
)

# --- REVISION: Path to pre-selected features (from INPUT_BASE_DIR) ---
PRESELECTED_FEATURES_JSON_PATH = os.path.join(
    INPUT_BASE_DIR, 
    "selected_top_1000_features_rfe.json"
)
# --- END REVISION ---

# --- REVISION: Path to pre-tuned parameters (from INPUT_BASE_DIR) ---
TUNED_PARAMS_JSON_PATH = os.path.join(
    INPUT_BASE_DIR, 
    "best_model_params.json"
)
# --- END REVISION ---

# --- REVISION: Model and metrics output paths (to OUTPUT_BASE_DIR) ---
MODEL_BIN_PATH = os.path.join(OUTPUT_BASE_DIR, "UpDEG_Prediction_trained_model.bin")
MODEL_JSON_PATH = os.path.join(OUTPUT_BASE_DIR, "UpDEG_Prediction_trained_model.json")
TEST_PREDICTIONS_PATH = os.path.join(OUTPUT_BASE_DIR, "test_predictions.csv")
TEST_METRICS_PATH = os.path.join(OUTPUT_BASE_DIR, "test_metrics.csv")
TEST_CONFUSION_PATH = os.path.join(OUTPUT_BASE_DIR, "test_confusion_matrix.csv")
FEATURE_IMPORTANCE_PATH = os.path.join(OUTPUT_BASE_DIR, "feature_importance.csv")

# --- MODIFIED FILE PATHS (to OUTPUT_BASE_DIR) ---
# This file is no longer loaded, but written to by this script, so it's an OUTPUT
RFE_STEP1_FEATURES_PATH = os.path.join(OUTPUT_BASE_DIR, "selected_top_1000_features_rfe.json")
SELECTED_FEATURES_PATH = os.path.join(OUTPUT_BASE_DIR, "selected_final_features.json")

# --- NEW: Additional model artifacts for ensembling ---
MODEL_LGBM_PATH = os.path.join(OUTPUT_BASE_DIR, "UpDEG_Prediction_trained_model_lgbm.txt")
MODEL_CATBOOST_PATH = os.path.join(OUTPUT_BASE_DIR, "UpDEG_Prediction_trained_model_catboost.cbm")
MODEL_RF_PATH = os.path.join(OUTPUT_BASE_DIR, "UpDEG_Prediction_trained_model_rf.pkl")
ENSEMBLE_MANIFEST_PATH = os.path.join(OUTPUT_BASE_DIR, "ensemble_manifest.json")
# --- END NEW ---

# --- NEW FILE PATHS FOR PLOTS (to OUTPUT_BASE_DIR) ---
AUC_CURVE_PATH = os.path.join(OUTPUT_BASE_DIR, "UpDEG_Prediction_AUC_curve.jpg")
AUCPR_CURVE_PATH = os.path.join(OUTPUT_BASE_DIR, "UpDEG_Prediction_AUCPR_curve.jpg")
# --- END REVISION ---

# --- NEW FILE PATHS FOR SHAP ANALYSIS ---
SHAP_VALUES_PATH = os.path.join(OUTPUT_BASE_DIR, "shap_values_test.csv")
SHAP_SUMMARY_PLOT_PATH = os.path.join(OUTPUT_BASE_DIR, "shap_summary_plot.jpg")
SHAP_BAR_PLOT_PATH = os.path.join(OUTPUT_BASE_DIR, "shap_bar_plot.jpg")
SHAP_WATERFALL_DIR = os.path.join(OUTPUT_BASE_DIR, "shap_waterfall_plots")
# --- END NEW ---

# Chromosome splits for datasets
# --- REVISION: TRAIN_CHRS and VALID_CHRS will be combined ---
TRAIN_CHRS = {f"chr{i}" for i in range(1, 14)}
VALID_CHRS = {"chr14", "chr15", "chr16", "chr17"}
# --- TEST_CHRS are still reserved exclusively for testing ---
TEST_CHRS = {"chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}

RANDOM_STATE = 42
# N_TRIALS = 2000  # <-- REVISION: Removed, Optuna is skipped
N_SPLITS_CV = 10 # Number of folds for cross-validation
warnings.filterwarnings("ignore", category=UserWarning)

# --- REVISION: COLUMNS TO REMOVE (NOW EMPTY) ---
COLUMNS_TO_REMOVE = []
# --- END REVISION ---

# --- REVISION: This variable will be defined after loading data ---
FEATURE_COLS = []
# --- REVISION: This global variable will hold chromosome info for CV ---
train_groups = None

print("Configuration loaded.")
if COLUMNS_TO_REMOVE:
    print(f"  NOTE: {len(COLUMNS_TO_REMOVE)} columns will be removed from all input files.")
else:
    print("  NOTE: No columns will be manually removed from input files.")

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
# --- REVISION: Create OUTPUT_BASE_DIR ---
print(f"Ensuring output directory exists: {OUTPUT_BASE_DIR}")
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)
# --- END REVISION ---

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
FEATURE_COLS_FULL = pos_train_pl.columns[21:-1] # Store the full list
FEATURE_COLS = FEATURE_COLS_FULL # FEATURE_COLS will be filtered later
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
print(f"  NOTE: This calculated value will be used for RFE, but overridden by the JSON for the final model.")


# -----------------------------------------------------------------
# --- REVISION: STAGE 1: LOAD PRE-SELECTED FEATURES (REPLACES RFE)
# -----------------------------------------------------------------
print("\n" + "=" * 50)
print("STAGE 1: Loading Pre-selected 1000 Features (RFE Step 1 Skipped)")
print("=" * 50)

# 1. Load the pre-selected features from the JSON file
print(f"  1. Loading features from {PRESELECTED_FEATURES_JSON_PATH}...")
try:
    with open(PRESELECTED_FEATURES_JSON_PATH, "r") as f:
        STEP1_FEATURE_COLS = json.load(f)
    print(f"  Loaded {len(STEP1_FEATURE_COLS)} features.")
except FileNotFoundError:
    print(f"\n--- ERROR: Pre-selected feature file not found! ---")
    print(f"  File not found: {PRESELECTED_FEATURES_JSON_PATH}")
    print("  This file is required to proceed.")
    exit(1)
except json.JSONDecodeError:
    print(f"\n--- ERROR: Could not parse feature file! ---")
    print(f"  File is not valid JSON: {PRESELECTED_FEATURES_JSON_PATH}")
    exit(1)

# 2. Save a copy of this feature list to the new output directory
#    (This is what the original script did in Stage 1, Step 5)
print(f"  2. Saving a copy of this feature list to {RFE_STEP1_FEATURES_PATH}...")
try:
    with open(RFE_STEP1_FEATURES_PATH, "w") as f:
        json.dump(STEP1_FEATURE_COLS, f, indent=4)
except Exception as e:
    print(f"  WARNING: Could not save feature list copy. Error: {e}")


# 3. Validate loaded features against the full dataset columns
print(f"  3. Validating loaded features against {len(FEATURE_COLS_FULL)} columns from dataset...")
missing_features = [col for col in STEP1_FEATURE_COLS if col not in FEATURE_COLS_FULL]
if missing_features:
    print(f"  WARNING: {len(missing_features)} loaded features are NOT in the main dataset.")
    print(f"    e.g.: {missing_features[:5]}")
    # Filter them out to prevent errors
    STEP1_FEATURE_COLS = [col for col in STEP1_FEATURE_COLS if col in FEATURE_COLS_FULL]
    print(f"    Proceeding with {len(STEP1_FEATURE_COLS)} valid features.")
else:
    print("  All loaded features are valid.")

# 4. Filter all datasets (X_train, X_test) to these features
print("  4. Filtering Train (Combined) and Test sets to Step 1 features...")
# CRITICAL: Overwrite the global FEATURE_COLS list
FEATURE_COLS = STEP1_FEATURE_COLS

X_train = X_train.select(FEATURE_COLS)
X_test = X_test.select(FEATURE_COLS)

print("  Datasets successfully filtered.")
print("=" * 50)
print("COMPLETED: STAGE 1: LOAD PRE-SELECTED FEATURES")
print("=" * 50 + "\n")
# -----------------------------------------------------------------
# --- END: REVISION
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# --- REVISION: STAGE 2: LOAD PRE-TUNED HYPERPARAMETERS (Replaces Optuna)
# -----------------------------------------------------------------
print("\n" + "=" * 50)
print(f"STAGE 2: Load Pre-Tuned Hyperparameters (Optuna Skipped)")
print("=" * 50)

# 1. Load the pre-tuned parameters from the JSON file
print(f"  1. Loading pre-tuned parameters from {TUNED_PARAMS_JSON_PATH}...")
try:
    with open(TUNED_PARAMS_JSON_PATH, "r") as f:
        loaded_params = json.load(f)
except FileNotFoundError:
    print(f"\n--- ERROR: Pre-tuned parameter file not found! ---")
    print(f"  File not found: {TUNED_PARAMS_JSON_PATH}")
    print("  This file is required to proceed.")
    exit(1)
except json.JSONDecodeError:
    print(f"\n--- ERROR: Could not parse parameter file! ---")
    print(f"  File is not valid JSON: {TUNED_PARAMS_JSON_PATH}")
    exit(1)

# 2. Extract parameters for RFE and Final Model
#    The 'best_iteration' will be determined by CV in Stage 3
try:
    best_n_features = loaded_params.pop("n_features_to_select", 150)
    # Remove factor, as scale_pos_weight is already provided
    loaded_params.pop("scale_pos_weight_factor", None) 
    
    final_params = loaded_params # This now contains all XGB params
    
    print(f"  Loaded n_features_to_select: {best_n_features}")
    print(f"  Loaded final XGBoost params: {final_params}")
    
    # Check for the crucial scale_pos_weight
    if "scale_pos_weight" not in final_params:
        print("  WARNING: 'scale_pos_weight' not in params file! Using calculated value.")
        final_params["scale_pos_weight"] = scale_pos_weight

except KeyError as e:
    print(f"\n--- ERROR: Parameter file is missing expected key: {e} ---")
    print("  Please ensure 'n_features_to_select' is in the JSON file.")
    exit(1)
except Exception as e:
    print(f"\n--- ERROR: Unexpected error while processing params: {e} ---")
    exit(1)

print("=" * 50)
print("COMPLETED: STAGE 2: LOAD PRE-TUNED HYPERPARAMETERS")
print("=" * 50 + "\n")
# -----------------------------------------------------------------
# --- END: REVISION
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# --- STAGE 3: RFE - Step 2 & Final Model Training
# -----------------------------------------------------------------
print("\n" + "=" * 50)
print("STAGE 3: RFE - Step 2 & Final Model Training")
print("=" * 50)

# 1. Run RFE Step 2 ONE MORE TIME using the best n_features found
print(f"  1. Running RFE - Step 2 to select best {best_n_features} features...")
# --- REVISION: This estimator MUST match the estimator from the Optuna trial ---
print("      Using lightweight RFE estimator to match tuning phase...")
rfe_estimator_final = xgb.XGBClassifier(
    n_estimators=50, # Lighter model to match original RFE
    max_depth=5,
    learning_rate=0.1,
    random_state=RANDOM_STATE,
    scale_pos_weight=scale_pos_weight, # Use *calculated* scale_pos_weight to match tuning RFE
    n_jobs=4, # Use fewer cores to match tuning RFE
)
# --- END REVISION ---

rfe_final = RFE(
    estimator=rfe_estimator_final, # Use the matched estimator
    n_features_to_select=best_n_features,
    step=0.2, # Use 0.2 step to match original RFE
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

# --- REVISION: Step 4 is now CV to find best iteration ---
print("\n  4. Determining best iteration using 10-fold Group CV...")
# Create the DMatrix with the FINAL RFE-selected features
dtrain_final = xgb.DMatrix(X_train.fill_null(0), label=y_train)
# --- FIX: Use .feature_names attribute instead of .set_feature_names() method ---
dtrain_final.feature_names = FEATURE_COLS # Good practice
# --- END FIX ---

# Create the GroupKFold generator
gkf = GroupKFold(n_splits=N_SPLITS_CV)
folds = list(gkf.split(X_train, y_train, groups=train_groups))

# Set up parameters for the CV
cv_params = final_params.copy()
cv_params["objective"] = "binary:logistic"
cv_params["eval_metric"] = "aucpr"
cv_params["n_jobs"] = 4 # Use fewer cores inside CV
cv_params["random_state"] = RANDOM_STATE

print(f"  Running {N_SPLITS_CV}-fold Group CV with {len(FEATURE_COLS)} features...")
try:
    cv_results_df = xgb.cv(
        params=cv_params,
        dtrain=dtrain_final,
        num_boost_round=2000,  # High number, will be stopped by early_stopping
        folds=folds,           # <-- Pass the chromosome-grouped folds
        metrics=["aucpr"],
        early_stopping_rounds=50,
        seed=RANDOM_STATE,
        verbose_eval=100, # Print progress every 100 rounds
    )
    
    # Get the best number of iterations (boosting rounds)
    best_iteration = len(cv_results_df)
    print(f"  CV complete. Best iteration (n_estimators): {best_iteration}")
    print(f"  Best Mean CV AUCPR: {cv_results_df['test-aucpr-mean'].iloc[-1]:.4f}")

except Exception as e:
    print(f"--- ERROR during xgb.cv: {e} ---")
    print("  Defaulting best_iteration to 100.")
    best_iteration = 100
# --- END REVISION ---

# 5. Train the FINAL XGBoost model
print("\n  5. Initializing FINAL XGBoost classifier...")
print(f"      Using best parameters from JSON: {final_params}")
print(f"      Training for {best_iteration} boosting rounds (from CV).")

# Use best_iteration, remove early_stopping
model = xgb.XGBClassifier(
    n_estimators=best_iteration,  # <-- Use best rounds from CV
    random_state=RANDOM_STATE,
    objective="binary:logistic",
    eval_metric=["logloss", "auc", "aucpr"], # Still useful for .predict_proba
    n_jobs=25, # Use all cores for final training
    early_stopping_rounds=None, # <-- No early stopping
    **final_params # Unpack the tuned XGBoost parameters
)
print("Final model initialized.")

print("Starting FINAL model training (on combined data)...")
start_time = time.time()
# Fit on combined data, no eval_set
model.fit(
    X_train.fill_null(0),  # <-- Pass FINAL filtered DataFrame (combined)
    y_train,               # <-- Pass combined labels
    verbose=False,
)
print(f"Model training complete in {time.time() - start_time:.2f} seconds.")

# 5a. Save the trained XGBoost model object
print(f"\n  5a. Saving trained XGBoost model (binary) to {MODEL_BIN_PATH}...")
model.save_model(MODEL_BIN_PATH)
print(f"      Saving trained XGBoost model (JSON) to {MODEL_JSON_PATH}...")
model.save_model(MODEL_JSON_PATH)
print("XGBoost model saved.")
print("  (Skipping train/validation metrics history export as CV was used)")

# --- NEW: 5b. Train additional diverse learners for ensembling ---
print("\n  5b. Training additional diverse learners for ensembling (LightGBM, CatBoost, Random Forest)...")

X_train_np = X_train.fill_null(0).to_numpy()
# Prepare containers for ensemble metadata
ensemble_models_info = [
    {
        "name": "xgboost",
        "path_bin": MODEL_BIN_PATH,
        "path_json": MODEL_JSON_PATH,
        "n_estimators": int(best_iteration),
        "status": "trained"
    }
]

# Train LightGBM if available
lgb_model = None
if lgb is not None:
    try:
        lgb_params = {
            "objective": "binary",
            "learning_rate": final_params.get("learning_rate", 0.1),
            "max_depth": final_params.get("max_depth", -1),
            "subsample": final_params.get("subsample", 1.0),
            "feature_fraction": final_params.get("colsample_bytree", 1.0),
            "min_child_weight": final_params.get("min_child_weight", 1.0),
            "min_split_gain": final_params.get("gamma", 0.0),
            "reg_lambda": final_params.get("reg_lambda", 1.0),
            "reg_alpha": final_params.get("reg_alpha", 0.0),
            "n_estimators": int(best_iteration),
            "random_state": RANDOM_STATE,
            "n_jobs": 25,
            "scale_pos_weight": final_params.get("scale_pos_weight", scale_pos_weight),
        }
        lgb_model = lgb.LGBMClassifier(**lgb_params)
        lgb_model.fit(X_train_np, y_train_np)
        # Save LightGBM booster if available
        try:
            booster = lgb_model.booster_
            booster.save_model(MODEL_LGBM_PATH)
        except Exception:
            # Fallback: pickle the sklearn wrapper
            with open(MODEL_LGBM_PATH, "wb") as f:
                pickle.dump(lgb_model, f)
        ensemble_models_info.append({
            "name": "lightgbm",
            "path": MODEL_LGBM_PATH,
            "n_estimators": int(best_iteration),
            "status": "trained"
        })
        print("  LightGBM trained and saved.")
    except Exception as e:
        print(f"  WARNING: LightGBM training skipped due to error: {e}")
        ensemble_models_info.append({
            "name": "lightgbm",
            "status": "error",
            "error": str(e)
        })
else:
    print("  LightGBM not installed; skipping.")
    ensemble_models_info.append({
        "name": "lightgbm",
        "status": "not_installed"
    })

# Train CatBoost if available
cat_model = None
if CatBoostClassifier is not None:
    try:
        cat_params = {
            "loss_function": "Logloss",
            "iterations": int(best_iteration),
            "learning_rate": final_params.get("learning_rate", 0.1),
            "depth": final_params.get("max_depth", 6),
            "l2_leaf_reg": final_params.get("reg_lambda", 1.0),
            "random_seed": RANDOM_STATE,
            "verbose": False,
            "scale_pos_weight": final_params.get("scale_pos_weight", scale_pos_weight),
        }
        cat_model = CatBoostClassifier(**cat_params)
        cat_model.fit(X_train_np, y_train_np)
        cat_model.save_model(MODEL_CATBOOST_PATH)
        ensemble_models_info.append({
            "name": "catboost",
            "path": MODEL_CATBOOST_PATH,
            "n_estimators": int(best_iteration),
            "status": "trained"
        })
        print("  CatBoost trained and saved.")
    except Exception as e:
        print(f"  WARNING: CatBoost training skipped due to error: {e}")
        ensemble_models_info.append({
            "name": "catboost",
            "status": "error",
            "error": str(e)
        })
else:
    print("  CatBoost not installed; skipping.")
    ensemble_models_info.append({
        "name": "catboost",
        "status": "not_installed"
    })

# --- NEW: Train Random Forest (always available via scikit-learn) ---
rf_model = None
try:
    rf_model = RandomForestClassifier(
        n_estimators=500,
        max_depth=None,
        min_samples_split=2,
        min_samples_leaf=1,
        max_features="sqrt",
        n_jobs=25,
        random_state=RANDOM_STATE,
        class_weight="balanced"
    )
    rf_model.fit(X_train_np, y_train_np)
    with open(MODEL_RF_PATH, "wb") as f:
        pickle.dump(rf_model, f)
    ensemble_models_info.append({
        "name": "random_forest",
        "path": MODEL_RF_PATH,
        "n_estimators": 500,
        "status": "trained"
    })
    print("  Random Forest trained and saved.")
except Exception as e:
    print(f"  WARNING: Random Forest training skipped due to error: {e}")
    ensemble_models_info.append({
        "name": "random_forest",
        "status": "error",
        "error": str(e)
    })

# Save ensemble manifest
try:
    manifest = {
        "models": ensemble_models_info,
        "ensemble": {
            "method": "average",
            "weights": "equal"
        },
        "features_used": FEATURE_COLS,
        "best_iteration_xgb": int(best_iteration)
    }
    with open(ENSEMBLE_MANIFEST_PATH, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"  Ensemble manifest saved to {ENSEMBLE_MANIFEST_PATH}.")
except Exception as e:
    print(f"  WARNING: Could not save ensemble manifest: {e}")

print("=" * 50)
print("COMPLETED: STAGE 3: RFE & FINAL MODEL TRAINING")
print("=" * 50 + "\n")
# -----------------------------------------------------------------
# --- END: STAGE 3
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# --- STAGE 4: Final Model Evaluation on Test Set
# -----------------------------------------------------------------
print("\n" + "=" * 50)
print("STAGE 4: Final Model Evaluation on Test Set")
print("=" * 50)
print("\nEvaluating model on the test set (using ensemble probabilities)...")

# Compute per-model probabilities
probas_list = []

# XGBoost probabilities
probas_list.append(model.predict_proba(X_test.fill_null(0))[:, 1])

# Convert test to numpy once for non-xgb models
X_test_np = X_test.fill_null(0).to_numpy()

# LightGBM probabilities if trained
if 'lgb_model' in locals() and lgb_model is not None:
    try:
        probas_list.append(lgb_model.predict_proba(X_test_np)[:, 1])
    except Exception as e:
        print(f"  WARNING: LightGBM prediction failed: {e}")

# CatBoost probabilities if trained
if 'cat_model' in locals() and cat_model is not None:
    try:
        probas_list.append(cat_model.predict_proba(X_test_np)[:, 1])
    except Exception as e:
        print(f"  WARNING: CatBoost prediction failed: {e}")

# Random Forest probabilities if trained
if 'rf_model' in locals() and rf_model is not None:
    try:
        probas_list.append(rf_model.predict_proba(X_test_np)[:, 1])
    except Exception as e:
        print(f"  WARNING: Random Forest prediction failed: {e}")

# Soft-vote ensemble (equal weights) over available models
if len(probas_list) == 0:
    print("  ERROR: No model probabilities available; aborting evaluation.")
    sys.exit(1)
elif len(probas_list) == 1:
    ensemble_prob = probas_list[0]
    print("  Ensemble reduced to single model (XGBoost).")
else:
    ensemble_prob = np.mean(np.vstack(probas_list), axis=0)
    print(f"  Ensemble combined {len(probas_list)} model probability vectors.")

# Threshold at 0.5 for labels
test_prob = ensemble_prob
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

print("=" * 50)
print("COMPLETED: STAGE 4: FINAL MODEL EVALUATION")
print("=" * 50 + "\n")
# -----------------------------------------------------------------
# --- END: STAGE 4
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# --- NEW: STAGE 5: SHAP Analysis for Model Interpretability
# -----------------------------------------------------------------
print("\n" + "=" * 50)
print("STAGE 5: SHAP Analysis for Model Interpretability")
print("=" * 50)

if shap is None:
    print("  WARNING: SHAP not installed. Skipping SHAP analysis.")
    print("  Install SHAP with: pip install shap")
else:
    try:
        print("\n  1. Computing SHAP values for XGBoost model...")
        # Convert test data to pandas for SHAP (required by TreeExplainer)
        X_test_df = X_test.fill_null(0).to_pandas()
        
        # Create SHAP explainer for XGBoost
        explainer_xgb = shap.TreeExplainer(model)
        shap_values_xgb = explainer_xgb.shap_values(X_test_df)
        
        print(f"     SHAP values computed for {X_test_df.shape[0]} test samples.")
        
        # 2. Save SHAP values to CSV
        print(f"\n  2. Saving SHAP values to {SHAP_VALUES_PATH}...")
        shap_df = pd.DataFrame(shap_values_xgb, columns=FEATURE_COLS)
        if "symbol" in test_pl.columns:
            shap_df.insert(0, 'symbol', test_pl.select("symbol").to_series().to_numpy())
        if "chr" in test_pl.columns:
            shap_df.insert(1, 'chr', test_pl.select("chr").to_series().to_numpy())
        shap_df['true_label'] = y_test_np
        shap_df['pred_label'] = test_pred
        shap_df['pred_proba'] = test_prob
        shap_df.to_csv(SHAP_VALUES_PATH, index=False)
        print("     SHAP values saved.")
        
        # 3. Generate SHAP summary plots
        print("\n  3. Generating SHAP summary plots...")
        
        # 3a. Beeswarm plot (summary plot)
        print(f"     Generating SHAP beeswarm plot and saving to {SHAP_SUMMARY_PLOT_PATH}...")
        plt.figure(figsize=(12, 8))
        shap.summary_plot(shap_values_xgb, X_test_df, show=False, max_display=20)
        plt.savefig(SHAP_SUMMARY_PLOT_PATH, dpi=300, bbox_inches='tight')
        plt.close()
        print("     SHAP beeswarm plot saved.")
        
        # 3b. Bar plot (mean absolute SHAP values)
        print(f"     Generating SHAP bar plot and saving to {SHAP_BAR_PLOT_PATH}...")
        plt.figure(figsize=(12, 8))
        shap.summary_plot(shap_values_xgb, X_test_df, plot_type="bar", 
                          show=False, max_display=20)
        plt.savefig(SHAP_BAR_PLOT_PATH, dpi=300, bbox_inches='tight')
        plt.close()
        print("     SHAP bar plot saved.")
        
        # 4. Generate waterfall plots for selected genes
        print("\n  4. Generating SHAP waterfall plots for specific predictions...")
        os.makedirs(SHAP_WATERFALL_DIR, exist_ok=True)
        
        # Select top 5 most confident positive predictions
        positive_indices = np.where(test_pred == 1)[0]
        if len(positive_indices) > 0:
            confident_positives = positive_indices[np.argsort(test_prob[positive_indices])[-5:]]
            
            print(f"     Generating waterfall plots for top 5 confident positive predictions...")
            for idx in confident_positives:
                gene_symbol = test_pl.select("symbol").to_series()[idx] if "symbol" in test_pl.columns else f"gene_{idx}"
                plt.figure(figsize=(10, 6))
                shap.waterfall_plot(
                    shap.Explanation(
                        values=shap_values_xgb[idx],
                        base_values=explainer_xgb.expected_value,
                        data=X_test_df.iloc[idx],
                        feature_names=FEATURE_COLS
                    ),
                    show=False,
                    max_display=15
                )
                waterfall_path = os.path.join(SHAP_WATERFALL_DIR, f"shap_waterfall_{gene_symbol}_pos.jpg")
                plt.savefig(waterfall_path, dpi=300, bbox_inches='tight')
                plt.close()
            print(f"     Saved {len(confident_positives)} positive waterfall plots.")
        
        # Select top 5 most confident negative predictions
        negative_indices = np.where(test_pred == 0)[0]
        if len(negative_indices) > 0:
            confident_negatives = negative_indices[np.argsort(test_prob[negative_indices])[:5]]
            
            print(f"     Generating waterfall plots for top 5 confident negative predictions...")
            for idx in confident_negatives:
                gene_symbol = test_pl.select("symbol").to_series()[idx] if "symbol" in test_pl.columns else f"gene_{idx}"
                plt.figure(figsize=(10, 6))
                shap.waterfall_plot(
                    shap.Explanation(
                        values=shap_values_xgb[idx],
                        base_values=explainer_xgb.expected_value,
                        data=X_test_df.iloc[idx],
                        feature_names=FEATURE_COLS
                    ),
                    show=False,
                    max_display=15
                )
                waterfall_path = os.path.join(SHAP_WATERFALL_DIR, f"shap_waterfall_{gene_symbol}_neg.jpg")
                plt.savefig(waterfall_path, dpi=300, bbox_inches='tight')
                plt.close()
            print(f"     Saved {len(confident_negatives)} negative waterfall plots.")
        
        # Select top 5 false positives
        false_positive_indices = np.where((test_pred == 1) & (y_test_np == 0))[0]
        if len(false_positive_indices) > 0:
            top_false_positives = false_positive_indices[np.argsort(test_prob[false_positive_indices])[-5:]]
            
            print(f"     Generating waterfall plots for top 5 false positives...")
            for idx in top_false_positives:
                gene_symbol = test_pl.select("symbol").to_series()[idx] if "symbol" in test_pl.columns else f"gene_{idx}"
                plt.figure(figsize=(10, 6))
                shap.waterfall_plot(
                    shap.Explanation(
                        values=shap_values_xgb[idx],
                        base_values=explainer_xgb.expected_value,
                        data=X_test_df.iloc[idx],
                        feature_names=FEATURE_COLS
                    ),
                    show=False,
                    max_display=15
                )
                waterfall_path = os.path.join(SHAP_WATERFALL_DIR, f"shap_waterfall_{gene_symbol}_FP.jpg")
                plt.savefig(waterfall_path, dpi=300, bbox_inches='tight')
                plt.close()
            print(f"     Saved {len(top_false_positives)} false positive waterfall plots.")
        
        # Select top 5 false negatives
        false_negative_indices = np.where((test_pred == 0) & (y_test_np == 1))[0]
        if len(false_negative_indices) > 0:
            top_false_negatives = false_negative_indices[np.argsort(test_prob[false_negative_indices])[:5]]
            
            print(f"     Generating waterfall plots for top 5 false negatives...")
            for idx in top_false_negatives:
                gene_symbol = test_pl.select("symbol").to_series()[idx] if "symbol" in test_pl.columns else f"gene_{idx}"
                plt.figure(figsize=(10, 6))
                shap.waterfall_plot(
                    shap.Explanation(
                        values=shap_values_xgb[idx],
                        base_values=explainer_xgb.expected_value,
                        data=X_test_df.iloc[idx],
                        feature_names=FEATURE_COLS
                    ),
                    show=False,
                    max_display=15
                )
                waterfall_path = os.path.join(SHAP_WATERFALL_DIR, f"shap_waterfall_{gene_symbol}_FN.jpg")
                plt.savefig(waterfall_path, dpi=300, bbox_inches='tight')
                plt.close()
            print(f"     Saved {len(top_false_negatives)} false negative waterfall plots.")
        
        print("\n  SHAP analysis complete.")
        print(f"     All SHAP outputs saved to {OUTPUT_BASE_DIR}")
        
    except Exception as e:
        print(f"\n--- WARNING: SHAP analysis failed with error: ---")
        print(f"  {e}")
        print("  Continuing without SHAP analysis.")

print("=" * 50)
print("COMPLETED: STAGE 5: SHAP ANALYSIS")
print("=" * 50 + "\n")
# -----------------------------------------------------------------
# --- END: STAGE 5
# -----------------------------------------------------------------

print("\n" + "=" * 50)
print("Script completed successfully!")
print("=" * 50)
print("\nAll outputs saved to:")
print(f"  {OUTPUT_BASE_DIR}")
print("\nKey outputs:")
print(f"  - Models: {MODEL_BIN_PATH}, {MODEL_JSON_PATH}")
print(f"  - Test predictions: {TEST_PREDICTIONS_PATH}")
print(f"  - Test metrics: {TEST_METRICS_PATH}")
print(f"  - Feature importance: {FEATURE_IMPORTANCE_PATH}")
print(f"  - AUC curves: {AUC_CURVE_PATH}, {AUCPR_CURVE_PATH}")
if shap is not None:
    print(f"  - SHAP values: {SHAP_VALUES_PATH}")
    print(f"  - SHAP plots: {SHAP_SUMMARY_PLOT_PATH}, {SHAP_BAR_PLOT_PATH}")
    print(f"  - SHAP waterfall plots: {SHAP_WATERFALL_DIR}/")
print("=" * 50)
