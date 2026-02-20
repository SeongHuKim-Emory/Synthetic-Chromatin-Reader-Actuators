# Python 3.12.3
# UpDEG vs Non-DEG classifier using XGBoost GBDT
# - Loads two CSVs (positive/negative)
# - Splits each by chromosome into train/validation/test
# - Saves model (.bin and .json), metrics, predictions, confusion matrix, and feature importance
#
# ---
# This version is optimized for I/O and processing performance using the 'polars' library.
# The expensive splitting operation is now handled entirely within polars
# before any conversion to pandas or numpy.
#
# ---
# MODIFIED:
# - Enabled CUDA for XGBoost by setting tree_method='hist' and device='cuda'.
# - Removed n_jobs=-1, as GPU parallelism is handled differently.
# - Removed sys.exit(0) to allow the ML training part of the script to run.
# ---
#
# You must install it first: pip install polars xgboost scikit-learn pandas
# (Pandas is still needed for some metric/utility functions)
# ---

import os
import json
import warnings
from typing import Dict, List, Tuple
import time
import sys

import numpy as np
import pandas as pd
import polars as pl  # <-- Import polars for fast I/O and processing
from sklearn.metrics import (
    log_loss,
    accuracy_score,
    roc_auc_score,
    average_precision_score,
    f1_score,
    confusion_matrix,
)
import xgboost as xgb

# ----------------------------
# Configuration and constants
# ----------------------------
print("Script started. Initializing configuration...")
BASE_DIR = "../../public/Multiomics/Step19_ML_UpDEG_Prediction_SRA"

POS_PATH = os.path.join(
    BASE_DIR,
    "SRA_Specific_UpDEG_with_GO.csv",
)
NEG_PATH = os.path.join(
    BASE_DIR,
    "Non_DEG_with_GO.csv",
)

# Output files for the positive dataset splits
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

# Output files for the negative dataset splits
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
TRAIN_VALID_METRICS_PATH = os.path.join(BASE_DIR, "train_validation_metrics.csv")
TEST_PREDICTIONS_PATH = os.path.join(BASE_DIR, "test_predictions.csv")
TEST_METRICS_PATH = os.path.join(BASE_DIR, "test_metrics.csv")
TEST_CONFUSION_PATH = os.path.join(BASE_DIR, "test_confusion_matrix.csv")
FEATURE_IMPORTANCE_PATH = os.path.join(BASE_DIR, "feature_importance.csv")

# Chromosome splits for datasets
TRAIN_CHRS = {f"chr{i}" for i in range(1, 14)}
VALID_CHRS = {"chr14", "chr15", "chr16", "chr17"}
TEST_CHRS = {"chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"}

RANDOM_STATE = 42
warnings.filterwarnings("ignore", category=UserWarning)

print("Configuration loaded.")

# ----------------------------
# Utility functions
# ----------------------------


def split_by_chr_pl(
    df: pl.DataFrame,
) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """
    Splits a Polars DataFrame into training, validation, and test sets
    based on chromosome using efficient Polars filtering.
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

# 1. Import files (Using polars for speed)
print(f"Loading positive dataset from {POS_PATH}...")
start_time = time.time()
# Load as Polars DataFrame and KEEP as Polars DataFrame
try:
    df_pos_pl = pl.read_csv(POS_PATH)
except Exception as e:
    print(f"ERROR: Could not load positive dataset from {POS_PATH}.")
    print(f"Details: {e}")
    sys.exit(1) # Exit if file loading fails

print(
    f"Positive dataset loaded in {time.time() - start_time:.2f} seconds. Shape: {df_pos_pl.shape}"
)

# Dynamically define feature columns
print("Defining feature columns...")
# Get columns directly from Polars DataFrame
FEATURE_COLS = df_pos_pl.columns[21:]
if not FEATURE_COLS:
    print("ERROR: No feature columns found (columns 21 and beyond).")
    print("Please check your input CSV format.")
    sys.exit(1)

print(
    f"Dynamically defined {len(FEATURE_COLS)} feature columns (from index 21 to end)."
)

print(f"Loading negative dataset from {NEG_PATH}...")
start_time = time.time()
# Load as Polars DataFrame and KEEP as Polars DataFrame
try:
    df_neg_pl = pl.read_csv(NEG_PATH)
except Exception as e:
    print(f"ERROR: Could not load negative dataset from {NEG_PATH}.")
    print(f"Details: {e}")
    sys.exit(1) # Exit if file loading fails

print(
    f"Negative dataset loaded in {time.time() - start_time:.2f} seconds. Shape: {df_neg_pl.shape}"
)


# 2 & 3. Assign labels for positive and negative sets (in Polars)
print("Assigning labels (1=positive, 0=negative)...")
# Use with_columns to add the label. This is non-mutating.
df_pos_pl = df_pos_pl.with_columns(pl.lit(1).cast(pl.Int32).alias("label"))
df_neg_pl = df_neg_pl.with_columns(pl.lit(0).cast(pl.Int32).alias("label"))

# 4 & 5. Split data and export to CSV (Using Polars for splitting and export)
print("Splitting positive and negative datasets by chromosome (using Polars)...")
start_split_time = time.time()
pos_train_pl, pos_valid_pl, pos_test_pl = split_by_chr_pl(df_pos_pl)
neg_train_pl, neg_valid_pl, neg_test_pl = split_by_chr_pl(df_neg_pl)
print(f"Splitting complete in {time.time() - start_split_time:.2f} seconds.")
print(
    f"  Positive - Train: {pos_train_pl.shape[0]}, Valid: {pos_valid_pl.shape[0]}, Test: {pos_test_pl.shape[0]}"
)
print(
    f"  Negative - Train: {neg_train_pl.shape[0]}, Valid: {neg_valid_pl.shape[0]}, Test: {neg_test_pl.shape[0]}"
)

# --- Fast CSV Export using Polars ---
print("Saving data splits to CSV...")
start_save_time = time.time()

try:
    print(f"  Saving {POS_TRAIN_PATH}...")
    pos_train_pl.write_csv(POS_TRAIN_PATH)
    print(f"  Saving {POS_VALID_PATH}...")
    pos_valid_pl.write_csv(POS_VALID_PATH)
    print(f"  Saving {POS_TEST_PATH}...")
    pos_test_pl.write_csv(POS_TEST_PATH)

    print(f"  Saving {NEG_TRAIN_PATH}...")
    neg_train_pl.write_csv(NEG_TRAIN_PATH)
    print(f"  Saving {NEG_VALID_PATH}...")
    neg_valid_pl.write_csv(NEG_VALID_PATH)
    print(f"  Saving {NEG_TEST_PATH}...")
    neg_test_pl.write_csv(NEG_TEST_PATH)
except Exception as e:
    print(f"ERROR: Could not save data splits to CSV.")
    print(f"Details: {e}")
    sys.exit(1)

print(f"All data splits saved in {time.time() - start_save_time:.2f} seconds.")

# --- XGBoost training and export sections removed as requested. ---
print("\nScript finished after saving data splits.")

