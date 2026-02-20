# Python 3.12.3
# UpDEG Prediction - Reproduction Script (With Plots + SHAP - No "Sum of Others")
# - Loads final features & pre-trained model
# - Loads Test Data
# - Exports: Metrics, Confusion Matrix, AUC Plot, AUCPR Plot
# - Exports: SHAP Bar Plot (Classic Blue Style), SHAP Beeswarm Plot (Top 20 Features ONLY)

import os
import json
import pandas as pd
import polars as pl
import xgboost as xgb
import numpy as np
import shap

# --- Plotting Imports ---
import matplotlib
matplotlib.use('Agg') # Use non-interactive backend
import matplotlib.pyplot as plt

from sklearn.metrics import (
    log_loss, accuracy_score, roc_auc_score, 
    average_precision_score, f1_score, confusion_matrix,
    roc_curve, precision_recall_curve
)

# ----------------------------
# Configuration
# ----------------------------
BASE_DIR = "../../public/Multiomics/Step19_ML_UpDEG_Prediction"

# Input Data Paths
POS_TEST_PATH = os.path.join(BASE_DIR, "GeneHancer_Genes_Elements_RNAseq_SRA_Specific_UpDEG_TSS_SRA_Max_Pool_Histone_Testing.csv")
NEG_TEST_PATH = os.path.join(BASE_DIR, "GeneHancer_Genes_Elements_RNAseq_Non_DEG_TSS_SRA_Max_Pool_Histone_Testing.csv")

# Input Artifacts
MODEL_BIN_PATH = os.path.join(BASE_DIR, "UpDEG_Prediction_trained_model.bin")
SELECTED_FEATURES_PATH = os.path.join(BASE_DIR, "selected_final_features.json")

# Output Paths for Standard Plots
AUC_CURVE_PATH = os.path.join(BASE_DIR, "reproduction_AUC_curve.jpg")
AUCPR_CURVE_PATH = os.path.join(BASE_DIR, "reproduction_AUCPR_curve.jpg")

# Output Paths for SHAP Plots
SHAP_BAR_PATH = os.path.join(BASE_DIR, "reproduction_SHAP_bar.jpg")
SHAP_BEESWARM_PATH = os.path.join(BASE_DIR, "reproduction_SHAP_beeswarm.jpg")

# Columns to remove (Alignment with training logic)
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

def compute_metrics(y_true, y_prob, y_pred):
    metrics = {}
    metrics["logloss"] = float(log_loss(y_true, y_prob))
    metrics["error"] = float(1.0 - accuracy_score(y_true, y_pred))
    metrics["auc"] = float(roc_auc_score(y_true, y_prob))
    metrics["aucpr"] = float(average_precision_score(y_true, y_prob))
    metrics["f1"] = float(f1_score(y_true, y_pred))
    return metrics

def plot_simple_shap_bar(shap_values, save_path, max_display=20):
    """
    Creates a simple bar chart of mean absolute SHAP values.
    Removes the numbers on bars and the color gradient.
    """
    # Calculate mean absolute SHAP values for each feature
    mean_shap = np.abs(shap_values.values).mean(axis=0)
    feature_names = shap_values.feature_names

    # Create a DataFrame for easier sorting and plotting
    importance_df = pd.DataFrame({
        'feature': feature_names,
        'importance': mean_shap
    }).sort_values(by='importance', ascending=True).tail(max_display)

    # Create the plot using matplotlib directly
    plt.figure(figsize=(8, len(importance_df) * 0.4 + 2)) # Adjust height based on feature count
    
    # Create horizontal bars
    bars = plt.barh(
        y=range(len(importance_df)),
        width=importance_df['importance'],
        color='#1f77b4', # Standard matplotlib blue
        height=0.7
    )

    # Add feature names as y-tick labels
    plt.yticks(range(len(importance_df)), importance_df['feature'], fontsize=10)
    
    # Add labels and title
    plt.xlabel('mean(|SHAP value|) (average impact on model output parameter)', fontsize=12)
    # plt.title('Feature Importance (SHAP)', fontsize=14) # Title is optional

    # Remove spines for cleaner look
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    # Add grid lines
    plt.grid(axis='x', linestyle='--', alpha=0.6)

    # Save plot
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight', dpi=300)
    plt.close()

# ----------------------------
# Execution
# ----------------------------
print(f"1. Loading feature list from {SELECTED_FEATURES_PATH}...")
with open(SELECTED_FEATURES_PATH, "r") as f:
    FINAL_FEATURE_COLS = json.load(f)
print(f"   Loaded {len(FINAL_FEATURE_COLS)} features.")

print("\n2. Loading Test Data...")
# Load and drop ignored columns
pos_test = pl.read_csv(POS_TEST_PATH)
neg_test = pl.read_csv(NEG_TEST_PATH)

pos_cols_drop = [c for c in COLUMNS_TO_REMOVE if c in pos_test.columns]
if pos_cols_drop: pos_test = pos_test.drop(pos_cols_drop)

neg_cols_drop = [c for c in COLUMNS_TO_REMOVE if c in neg_test.columns]
if neg_cols_drop: neg_test = neg_test.drop(neg_cols_drop)

# Combine
test_pl = pl.concat([pos_test, neg_test], how="vertical_relaxed")

# Select X (features) and y (label)
X_test = test_pl.select(FINAL_FEATURE_COLS).fill_null(0) 
y_test_series = test_pl.select("label").to_series()
y_test_np = y_test_series.to_numpy().astype(int)

print(f"   Test Data Loaded: {X_test.shape[0]} samples.")

print("\n3. Loading Pre-trained Model...")
model = xgb.XGBClassifier() 
model.load_model(MODEL_BIN_PATH)
print("   Model loaded successfully.")

print("\n4. Running Predictions...")
test_prob = model.predict_proba(X_test)[:, 1]
test_pred = (test_prob >= 0.5).astype(int)

print("\n5. Reproducing Metrics...")
metrics = compute_metrics(y_test_np, test_prob, test_pred)
metrics_df = pd.DataFrame([metrics])
print("-" * 30)
print(metrics_df.to_string(index=False))
print("-" * 30)

# ----------------------------
# PLOTTING SECTION (Metrics)
# ----------------------------
print("\n6. Generating Metric Plots...")
try:
    # --- Plot 1: ROC Curve ---
    print(f"   Saving AUC curve to {AUC_CURVE_PATH}...")
    fpr, tpr, _ = roc_curve(y_test_np, test_prob)
    auc_score = metrics.get("auc", 0.0)

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
    plt.close()

    # --- Plot 2: PR Curve ---
    print(f"   Saving AUCPR curve to {AUCPR_CURVE_PATH}...")
    precision, recall, _ = precision_recall_curve(y_test_np, test_prob)
    aucpr_score = metrics.get("aucpr", 0.0)

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
    plt.close()
except Exception as e:
    print(f"   ERROR: Could not generate metric plots. {e}")

# ----------------------------
# PLOTTING SECTION (SHAP)
# ----------------------------
print("\n7. Generating SHAP Plots...")
try:
    # Convert Polars to Pandas for SHAP compatibility
    print("   Converting Test Data to Pandas for SHAP (this may take a moment)...")
    X_test_pd = X_test.to_pandas()

    print("   Calculating SHAP values...")
    # Initialize the explainer
    explainer = shap.TreeExplainer(model)
    # Calculate SHAP values
    shap_values = explainer(X_test_pd)
    
    print("   Filtering top 20 features to remove 'Sum of other features'...")
    # 1. Compute global importance (mean absolute SHAP value)
    mean_abs_shap = np.abs(shap_values.values).mean(axis=0)
    
    # 2. Get indices of the top 20 features
    # argsort sorts ascending, so we slice from the end [::-1] and take top 20
    top_indices = np.argsort(mean_abs_shap)[::-1][:20]
    
    # 3. Slice the SHAP explanation object to strictly keep ONLY these features.
    shap_values_top = shap_values[:, top_indices]

    # --- SHAP Bar Plot (Classic Style) ---
    print(f"   Saving SHAP Bar Plot to {SHAP_BAR_PATH}...")
    # Use the custom function instead of shap.plots.bar
    plot_simple_shap_bar(shap_values_top, SHAP_BAR_PATH, max_display=20)

    # --- SHAP Beeswarm Plot ---
    print(f"   Saving SHAP Beeswarm Plot to {SHAP_BEESWARM_PATH}...")
    plt.figure()
    shap.plots.beeswarm(shap_values_top, max_display=20, show=False)
    plt.savefig(SHAP_BEESWARM_PATH, bbox_inches='tight')
    plt.close()
    
    print("   SHAP plots generated successfully.")

except Exception as e:
    print(f"   ERROR: Could not generate SHAP plots. {e}")
    print("   Note: Ensure 'shap' is installed (pip install shap).")

print("\nDone.")