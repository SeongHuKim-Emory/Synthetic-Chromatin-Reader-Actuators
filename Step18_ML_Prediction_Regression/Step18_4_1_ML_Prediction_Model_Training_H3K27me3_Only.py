import pandas as pd
import xgboost as xgb
from sklearn.metrics import (
    log_loss,
    roc_auc_score,
    average_precision_score,
    f1_score,
    accuracy_score,
    confusion_matrix,
    roc_curve,
    precision_recall_curve,
)
import os
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuration ---
# Define file paths for inputs and outputs
POSITIVE_SET_PATH = "../../public/multiomics/Step18_ML_Prediction_Regression/hg38_bin_positive_histone.csv"
NEGATIVE_SET_PATH = "../../public/multiomics/Step18_ML_Prediction_Regression/hg38_bin_negative_histone.csv"
OUTPUT_DIR = "../../public/multiomics/Step18_ML_Prediction"

# Define the suffix for all output files
FILE_SUFFIX = "_H4K20me1_only"

# Create the output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Ensured output directory exists: {OUTPUT_DIR}")

# --- Step 1 & 2: Import and Combine Data ---
print("\nStep 1 & 2: Importing and combining positive and negative sets...")
try:
    positive_df = pd.read_csv(POSITIVE_SET_PATH)
    negative_df = pd.read_csv(NEGATIVE_SET_PATH)

    hg38_bin_feature_matrix_histone = pd.concat(
        [positive_df, negative_df], ignore_index=True
    )
    print("Successfully combined dataframes.")
    print(f"Shape of the combined matrix: {hg38_bin_feature_matrix_histone.shape}")

except FileNotFoundError as e:
    print(f"Error: Could not find input file - {e}. Please check the file paths.")
    exit()

# --- Step 3: Export Combined Feature Matrix ---
print("\nStep 3: Exporting the combined feature matrix...")
combined_output_path = os.path.join(
    OUTPUT_DIR, f"hg38_bin_feature_matrix_histone{FILE_SUFFIX}.csv"
)
hg38_bin_feature_matrix_histone.to_csv(combined_output_path, index=False)
print(f"Combined feature matrix saved to: {combined_output_path}")

# --- Step 4: Split Data into Training, Validation, and Test Sets ---
print("\nStep 4: Splitting the data based on chromosomes...")
# Define chromosome groups for each dataset
train_chromosomes = [f"chr{i}" for i in range(1, 14)]
validation_chromosomes = [f"chr{i}" for i in range(14, 18)]
test_chromosomes = [f"chr{i}" for i in range(18, 23)] + ["chrX", "chrY"]

# Filter the main dataframe based on the chromosome lists
training_df = hg38_bin_feature_matrix_histone[
    hg38_bin_feature_matrix_histone["chr"].isin(train_chromosomes)
].copy()
validation_df = hg38_bin_feature_matrix_histone[
    hg38_bin_feature_matrix_histone["chr"].isin(validation_chromosomes)
].copy()
test_df = hg38_bin_feature_matrix_histone[
    hg38_bin_feature_matrix_histone["chr"].isin(test_chromosomes)
].copy()

print(f"Training set shape:   {training_df.shape}")
print(f"Validation set shape: {validation_df.shape}")
print(f"Test set shape:       {test_df.shape}")

# --- Step 5, 6, & 7: Export Data Splits ---
print("\nSteps 5, 6, & 7: Exporting data splits...")
training_df.to_csv(
    os.path.join(
        OUTPUT_DIR, f"hg38_bin_feature_matrix_histone_training{FILE_SUFFIX}.csv"
    ),
    index=False,
)
validation_df.to_csv(
    os.path.join(
        OUTPUT_DIR, f"hg38_bin_feature_matrix_histone_validation{FILE_SUFFIX}.csv"
    ),
    index=False,
)
test_df.to_csv(
    os.path.join(
        OUTPUT_DIR, f"hg38_bin_feature_matrix_histone_test{FILE_SUFFIX}.csv"
    ),
    index=False,
)
print("Training, validation, and test sets have been saved successfully.")

# --- Step 8: Model Training with XGBoost ---
print("\nStep 8: Training the XGBoost GBDT model...")

# Define features (all columns except the first four) and the label
# *** MODIFICATION: Only use H3K27me3 as the feature ***
features = [
    "H4K20me1"
]
label = "PoI"

# Prepare data for XGBoost
X_train, y_train = training_df[features], training_df[label]
X_val, y_val = validation_df[features], validation_df[label]
X_test, y_test = test_df[features], test_df[label]

# Create DMatrix objects, which are optimized for XGBoost
dtrain = xgb.DMatrix(X_train, label=y_train)
dval = xgb.DMatrix(X_val, label=y_val)
dtest = xgb.DMatrix(X_test, label=y_test)

# Set XGBoost hyperparameters
params = {
    "objective": "binary:logistic",
    "eta": 0.1,
    "max_depth": 6,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "seed": 42,
    "eval_metric": ["logloss", "error", "auc", "aucpr"],
}

# Dictionary to store evaluation results
evals_result = {}

# Train the model
model = xgb.train(
    params=params,
    dtrain=dtrain,
    num_boost_round=1000,  # Max number of rounds
    evals=[(dtrain, "train"), (dval, "validation")],
    early_stopping_rounds=50,  # Stop if validation metric doesn't improve for 50 rounds
    evals_result=evals_result,
    verbose_eval=100,  # Print evaluation results every 100 rounds
)

print("Model training completed.")

# --- Step 9 & 10: Export the Trained Model ---
print("\nSteps 9 & 10: Exporting the trained model object...")
model.save_model(
    os.path.join(OUTPUT_DIR, f"hg38_bin_histone_trained_model{FILE_SUFFIX}.bin")
)
model.save_model(
    os.path.join(OUTPUT_DIR, f"hg38_bin_histone_trained_model{FILE_SUFFIX}.json")
)
print("Model saved in both binary (.bin) and JSON (.json) formats.")


# --- Step 11: Export Training and Validation Performance Metrics ---
print("\nStep 11: Exporting training and validation performance metrics...")
# Create a dataframe from the evals_result dictionary
train_metrics_df = pd.DataFrame(evals_result["train"])
train_metrics_df["set"] = "training"
val_metrics_df = pd.DataFrame(evals_result["validation"])
val_metrics_df["set"] = "validation"

# Combine and add an 'epoch' column
performance_df = pd.concat([train_metrics_df, val_metrics_df])
performance_df["epoch"] = performance_df.index + 1

# Calculate F1 score for the validation set at the best iteration
best_iteration = model.best_iteration
val_preds_proba = model.predict(dval, iteration_range=(0, best_iteration + 1))
val_preds_class = (val_preds_proba > 0.5).astype(int)
val_f1 = f1_score(y_val, val_preds_class)

print(f"Validation F1 Score at best iteration ({best_iteration}): {val_f1:.4f}")

# Save the primary metrics
performance_df.to_csv(
    os.path.join(OUTPUT_DIR, f"training_validation_metrics{FILE_SUFFIX}.csv"),
    index=False,
)
# Save the F1 score
f1_df = pd.DataFrame([{"set": "validation", "metric": "f1_score", "value": val_f1}])
f1_df.to_csv(
    os.path.join(OUTPUT_DIR, f"validation_f1_score{FILE_SUFFIX}.csv"), index=False
)
print("Training and validation metrics saved.")


# --- Step 12: Evaluate Model on the Test Set ---
print("\nStep 12: Evaluating the model on the unseen test set...")
test_preds_proba = model.predict(dtest, iteration_range=(0, best_iteration + 1))
test_preds_class = (test_preds_proba > 0.5).astype(int)
print("Predictions generated for the test set.")


# --- Step 13: Export Test Performance and Artifacts ---
print("\nStep 13: Exporting test results...")

# 1. Predictions with class labels and probability scores
predictions_df = pd.DataFrame(
    {
        "chr": test_df["chr"],
        "start": test_df["start"],
        "end": test_df["end"],
        "true_label": y_test,
        "predicted_label": test_preds_class,
        "predicted_probability": test_preds_proba,
    }
)
predictions_df.to_csv(
    os.path.join(OUTPUT_DIR, f"test_predictions{FILE_SUFFIX}.csv"), index=False
)
print("- Test predictions saved.")

# 2. Test performance metrics
test_metrics = {
    "logloss": log_loss(y_test, test_preds_proba),
    "error": 1 - accuracy_score(y_test, test_preds_class),
    "auc": roc_auc_score(y_test, test_preds_proba),
    "aucpr": average_precision_score(y_test, test_preds_proba),
    "f1": f1_score(y_test, test_preds_class),
}
test_metrics_df = pd.DataFrame([test_metrics])
test_metrics_df.to_csv(
    os.path.join(OUTPUT_DIR, f"test_performance_metrics{FILE_SUFFIX}.csv"),
    index=False,
)
print("- Test performance metrics saved.")

# 3. Confusion matrix
tn, fp, fn, tp = confusion_matrix(y_test, test_preds_class).ravel()
cm_df = pd.DataFrame(
    {
        "Metric": ["True Negative", "False Positive", "False Negative", "True Positive"],
        "Value": [tn, fp, fn, tp],
    }
)
cm_df.to_csv(
    os.path.join(OUTPUT_DIR, f"test_confusion_matrix{FILE_SUFFIX}.csv"), index=False
)
print("- Test confusion matrix saved.")

# 4. Feature Importance Scores
importance_types = ["weight", "gain", "cover"]
feature_importance_list = []
for imp_type in importance_types:
    scores = model.get_score(importance_type=imp_type)
    temp_df = pd.DataFrame(
        {"feature": list(scores.keys()), "score": list(scores.values()), "type": imp_type}
    )
    feature_importance_list.append(temp_df)

feature_importance_df = pd.concat(feature_importance_list, ignore_index=True)
feature_importance_df.to_csv(
    os.path.join(OUTPUT_DIR, f"feature_importance_scores{FILE_SUFFIX}.csv"),
    index=False,
)
print("- Feature importance scores saved.")

print("- Generating and saving plots...")

# 5. ROC Curve
fpr, tpr, _ = roc_curve(y_test, test_preds_proba)
roc_auc = roc_auc_score(y_test, test_preds_proba)

plt.figure(figsize=(8, 6))
plt.plot(
    fpr, tpr, color="darkorange", lw=2, label=f"ROC curve (area = {roc_auc:.2f})"
)
plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver Operating Characteristic (ROC) Curve")
plt.legend(loc="lower right")
plt.grid(True)
plt.savefig(
    os.path.join(OUTPUT_DIR, f"roc_curve{FILE_SUFFIX}.jpg"), format="jpg", dpi=300
)
plt.close()
print("- ROC curve plot saved.")


# 6. Precision-Recall Curve
precision, recall, _ = precision_recall_curve(y_test, test_preds_proba)
avg_precision = average_precision_score(y_test, test_preds_proba)

plt.figure(figsize=(8, 6))
plt.step(
    recall, precision, where="post", label=f"PR curve (area = {avg_precision:.2f})"
)
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.title("Precision-Recall Curve")
plt.legend(loc="upper right")
plt.grid(True)
plt.savefig(
    os.path.join(OUTPUT_DIR, f"precision_recall_curve{FILE_SUFFIX}.jpg"),
    format="jpg",
    dpi=300,
)
plt.close()
print("- Precision-Recall curve plot saved.")


# 7. Confusion Matrix Heatmap
cm = confusion_matrix(y_test, test_preds_class)
plt.figure(figsize=(8, 6))
sns.heatmap(
    cm,
    annot=True,
    fmt="d",
    cmap="Blues",
    xticklabels=["Predicted Negative", "Predicted Positive"],
    yticklabels=["Actual Negative", "Actual Positive"],
)
plt.ylabel("Actual Label")
plt.xlabel("Predicted Label")
plt.title("Confusion Matrix Heatmap")
plt.savefig(
    os.path.join(OUTPUT_DIR, f"confusion_matrix_heatmap{FILE_SUFFIX}.jpg"),
    format="jpg",
    dpi=300,
)
plt.close()
print("- Confusion matrix heatmap saved.")

print("\n--- All tasks completed successfully! ---")
