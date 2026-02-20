#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import xgboost as xgb
import optuna

from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.model_selection import GroupKFold

import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

POSITIVE_SET_PATH = "../../public/multiomics/Step18_ML_Prediction_Regression/hg38_bin_positive_histone.csv"
NEGATIVE_SET_PATH = "../../public/multiomics/Step18_ML_Prediction_Regression/hg38_bin_negative_histone.csv"

OUTPUT_DIR = "../../public/multiomics/Step18_ML_Prediction_Regression"
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Ensured output directory exists: {OUTPUT_DIR}")

# ---------------------------------------------------------------------------
# Step 1 & 2: Import and Combine Data
# ---------------------------------------------------------------------------

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
    raise SystemExit(1)

# ---------------------------------------------------------------------------
# New: Target transform and tail‑focused sample weights
# ---------------------------------------------------------------------------

# Log‑transform target for regression
hg38_bin_feature_matrix_histone["PoI_log1p"] = np.log1p(
    hg38_bin_feature_matrix_histone["PoI"]
)

poi_all = hg38_bin_feature_matrix_histone["PoI"]

# Quantiles used to define "high‑PoI" tail
q90_all = poi_all.quantile(0.90)
q95_all = poi_all.quantile(0.95)

# Base weight 1, higher weights for the tail
weights = np.ones(len(poi_all), dtype=float)
weights[poi_all >= q90_all] = 3.0     # top 10%
weights[poi_all >= q95_all] = 8.0     # top 5% get largest weight

hg38_bin_feature_matrix_histone["poi_weight"] = weights

print(
    "Created transformed target 'PoI_log1p' and tail‑focused weights 'poi_weight'. "
    f"q90={q90_all:.3f}, q95={q95_all:.3f}, "
    f"weight range={weights.min():.1f}–{weights.max():.1f}"
)

# ---------------------------------------------------------------------------
# Step 3: Export Combined Feature Matrix
# ---------------------------------------------------------------------------

print("\nStep 3: Exporting the combined feature matrix...")
combined_output_path = os.path.join(
    OUTPUT_DIR, "hg38_bin_feature_matrix_histone_with_logtarget_tailweights.csv"
)
hg38_bin_feature_matrix_histone.to_csv(combined_output_path, index=False)
print(f"Combined feature matrix saved to: {combined_output_path}")

# ---------------------------------------------------------------------------
# Step 4: Split Data into Training, Validation, and Test Sets
# ---------------------------------------------------------------------------

print("\nStep 4: Splitting the data based on chromosomes...")

train_chromosomes = [f"chr{i}" for i in range(1, 14)]
validation_chromosomes = [f"chr{i}" for i in range(14, 18)]
test_chromosomes = [f"chr{i}" for i in range(18, 23)] + ["chrX", "chrY"]

training_df = hg38_bin_feature_matrix_histone[
    hg38_bin_feature_matrix_histone["chr"].isin(train_chromosomes)
].copy()
validation_df = hg38_bin_feature_matrix_histone[
    hg38_bin_feature_matrix_histone["chr"].isin(validation_chromosomes)
].copy()
test_df = hg38_bin_feature_matrix_histone[
    hg38_bin_feature_matrix_histone["chr"].isin(test_chromosomes)
].copy()

print(f"Training set shape: {training_df.shape}")
print(f"Validation set shape: {validation_df.shape}")
print(f"Test set shape: {test_df.shape}")

# ---------------------------------------------------------------------------
# Steps 5–7: Export Data Splits
# ---------------------------------------------------------------------------

print("\nSteps 5–7: Exporting data splits...")

training_df.to_csv(
    os.path.join(OUTPUT_DIR, "hg38_bin_feature_matrix_histone_training.csv"),
    index=False,
)
validation_df.to_csv(
    os.path.join(OUTPUT_DIR, "hg38_bin_feature_matrix_histone_validation.csv"),
    index=False,
)
test_df.to_csv(
    os.path.join(OUTPUT_DIR, "hg38_bin_feature_matrix_histone_test.csv"),
    index=False,
)

print("Training, validation, and test sets have been saved successfully.")

# ---------------------------------------------------------------------------
# Feature and label definitions
# ---------------------------------------------------------------------------

features = [
    "H3K4me1",
    "H3K4me2",
    "H3K4me3",
    "H3K9ac",
    "H3K9me3",
    "H3K27ac",
    "H3K27me3",
    "H3K36me3",
    "H4K20me1",
]

label_log = "PoI_log1p"
label_original = "PoI"
weight_col = "poi_weight"

# Base data for main model
X_train = training_df[features]
y_train_log = training_df[label_log]
w_train = training_df[weight_col]

X_val = validation_df[features]
y_val_log = validation_df[label_log]
w_val = validation_df[weight_col]

X_test = test_df[features]
y_test_log = test_df[label_log]
y_test_original = test_df[label_original]
w_test = test_df[weight_col]

dtrain = xgb.DMatrix(X_train, label=y_train_log, weight=w_train)
dval = xgb.DMatrix(X_val, label=y_val_log, weight=w_val)
dtest = xgb.DMatrix(X_test, label=y_test_log, weight=w_test)

# ---------------------------------------------------------------------------
# Step 8a: Hyperparameter Tuning with Optuna + GroupKFold
# ---------------------------------------------------------------------------

print("\nStep 8a: Starting Hyperparameter Tuning with Optuna and 10‑fold Group CV...")

train_val_df = pd.concat([training_df, validation_df], ignore_index=True)
X_train_val = train_val_df[features]
y_train_val_log = train_val_df[label_log]
w_train_val = train_val_df[weight_col]
groups = train_val_df["chr"]

dtrain_val = xgb.DMatrix(X_train_val, label=y_train_val_log, weight=w_train_val)

group_kfold = GroupKFold(n_splits=10)
folds = list(group_kfold.split(X_train_val, y_train_val_log, groups))


def objective(trial: optuna.Trial) -> float:
    params = {
        "objective": "reg:squarederror",
        "seed": 42,
        "tree_method": "hist",
        "device": "cuda",  # comment out if no GPU
        "eta": trial.suggest_float("eta", 0.01, 0.3, log=True),
        "max_depth": trial.suggest_int("max_depth", 3, 10),
        "min_child_weight": trial.suggest_float("min_child_weight", 1.0, 20.0),
        "gamma": trial.suggest_float("gamma", 0.0, 5.0),
        "subsample": trial.suggest_float("subsample", 0.5, 1.0),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
        "lambda": trial.suggest_float("lambda", 0.0, 10.0),
        "alpha": trial.suggest_float("alpha", 0.0, 10.0),
    }

    cv_results = xgb.cv(
        params=params,
        dtrain=dtrain_val,
        num_boost_round=2000,
        folds=folds,
        metrics=["rmse"],  # RMSE in log space
        maximize=False,
        early_stopping_rounds=100,
        seed=42,
        verbose_eval=False,
    )

    best_mean_rmse = cv_results["test-rmse-mean"].min()
    best_iter = cv_results["test-rmse-mean"].idxmin() + 1
    trial.set_user_attr("best_iteration", int(best_iter))
    return float(best_mean_rmse)


n_trials = 2  # increase if you have more compute
print(f"Starting Optuna optimization with {n_trials} trials...")
study = optuna.create_study(direction="minimize")
study.optimize(objective, n_trials=n_trials)

print("\nHyperparameter tuning completed.")
best_rmse = study.best_trial.value
best_params = study.best_trial.params
best_iteration = study.best_trial.user_attrs["best_iteration"]
print(f"Best Mean RMSE (log space): {best_rmse:.6f}")
print(f"Best Parameters found: {best_params}")
print(f"Optimal Iterations (from CV): {best_iteration}")

tuning_results_df = study.trials_dataframe().sort_values("value", ascending=True)
tuning_results_df.to_csv(
    os.path.join(OUTPUT_DIR, "hyperparameter_tuning_results.csv"), index=False
)
print("Tuning results saved to hyperparameter_tuning_results.csv")

# ---------------------------------------------------------------------------
# Step 8b: Final main model training
# ---------------------------------------------------------------------------

print("\nStep 8b: Training the final main XGBoost model (log target, tail weights)...")

final_params = best_params.copy()
final_params.update(
    {
        "objective": "reg:squarederror",
        "seed": 42,
        "eval_metric": ["rmse", "mae"],  # still in log space
        "tree_method": "hist",
        "device": "cuda",  # comment out if no GPU
    }
)

evals_result = {}
model = xgb.train(
    params=final_params,
    dtrain=dtrain,
    num_boost_round=2000,
    evals=[(dtrain, "train"), (dval, "validation")],
    early_stopping_rounds=100,
    evals_result=evals_result,
    verbose_eval=100,
)
print("Final main model training completed.")
print(f"Best iteration (from final train): {model.best_iteration}")

# ---------------------------------------------------------------------------
# Optional: High‑PoI expert model (trained only on tail bins)
# ---------------------------------------------------------------------------

print("\nTraining high‑PoI expert model on top 5% PoI bins...")

high_threshold = q95_all  # same global 95th percentile
train_hi_df = training_df[training_df[label_original] >= high_threshold].copy()
val_hi_df = validation_df[validation_df[label_original] >= high_threshold].copy()
test_hi_df = test_df[test_df[label_original] >= high_threshold].copy()

print(
    f"High‑PoI threshold (expert model): PoI >= {high_threshold:.3f}. "
    f"Train_hi={len(train_hi_df)}, Val_hi={len(val_hi_df)}, Test_hi={len(test_hi_df)}"
)

X_train_hi = train_hi_df[features]
y_train_hi_log = train_hi_df[label_log]
w_train_hi = train_hi_df[weight_col]

X_val_hi = val_hi_df[features]
y_val_hi_log = val_hi_df[label_log]
w_val_hi = val_hi_df[weight_col]

X_test_hi = test_hi_df[features]
y_test_hi_original = test_hi_df[label_original]
y_test_hi_log = test_hi_df[label_log]
w_test_hi = test_hi_df[weight_col]

dtrain_hi = xgb.DMatrix(X_train_hi, label=y_train_hi_log, weight=w_train_hi)
dval_hi = xgb.DMatrix(X_val_hi, label=y_val_hi_log, weight=w_val_hi)
dtest_hi = xgb.DMatrix(X_test_hi, label=y_test_hi_log, weight=w_test_hi)

expert_evals_result = {}
expert_model = xgb.train(
    params=final_params,
    dtrain=dtrain_hi,
    num_boost_round=model.best_iteration,  # reuse main model length
    evals=[(dtrain_hi, "train_hi"), (dval_hi, "validation_hi")],
    early_stopping_rounds=50,
    evals_result=expert_evals_result,
    verbose_eval=50,
)
print("High‑PoI expert model training completed.")
print(f"Expert model best iteration: {expert_model.best_iteration}")

# ---------------------------------------------------------------------------
# Steps 9 & 10: Save models
# ---------------------------------------------------------------------------

print("\nSteps 9 & 10: Exporting trained models...")

model.save_model(os.path.join(OUTPUT_DIR, "hg38_bin_histone_trained_model_main.bin"))
model.save_model(os.path.join(OUTPUT_DIR, "hg38_bin_histone_trained_model_main.json"))

expert_model.save_model(
    os.path.join(OUTPUT_DIR, "hg38_bin_histone_trained_model_expert.bin")
)
expert_model.save_model(
    os.path.join(OUTPUT_DIR, "hg38_bin_histone_trained_model_expert.json")
)

print("Both main and expert models saved.")

# ---------------------------------------------------------------------------
# Step 11: Save training / validation metrics (main model)
# ---------------------------------------------------------------------------

print("\nStep 11: Exporting training and validation performance metrics (main model)...")

train_metrics_df = pd.DataFrame(evals_result["train"])
train_metrics_df["set"] = "training"

val_metrics_df = pd.DataFrame(evals_result["validation"])
val_metrics_df["set"] = "validation"

performance_df = pd.concat([train_metrics_df, val_metrics_df])
performance_df["epoch"] = performance_df.index + 1

performance_df.to_csv(
    os.path.join(OUTPUT_DIR, "training_validation_metrics_log_space.csv"),
    index=False,
)
print("Training/validation metrics (log space) saved.")

# ---------------------------------------------------------------------------
# Step 12: Evaluate models on test set
# ---------------------------------------------------------------------------

print("\nStep 12: Evaluating models on the unseen test set...")

# Main model predictions (all test bins)
test_preds_log_main = model.predict(
    dtest, iteration_range=(0, model.best_iteration + 1)
)
test_preds_main = np.expm1(test_preds_log_main)

print("Main model predictions generated for full test set.")

# Expert model predictions (only high‑PoI test bins)
test_preds_log_expert = expert_model.predict(
    dtest_hi, iteration_range=(0, expert_model.best_iteration + 1)
)
test_preds_expert = np.expm1(test_preds_log_expert)

print("Expert model predictions generated for high‑PoI test subset.")

# ---------------------------------------------------------------------------
# Step 13: Export test results and detailed metrics
# ---------------------------------------------------------------------------

print("\nStep 13: Exporting test results and high‑PoI metrics...")

# 1) Per‑bin predictions (main model on full test set)
predictions_df = pd.DataFrame(
    {
        "chr": test_df["chr"],
        "start": test_df["start"],
        "end": test_df["end"],
        "true_label": y_test_original,
        "true_label_log1p": y_test_log,
        "predicted_value_main": test_preds_main,
        "predicted_log1p_main": test_preds_log_main,
    }
)
predictions_df.to_csv(
    os.path.join(OUTPUT_DIR, "test_predictions_log_model_main.csv"), index=False
)
print("- Test predictions (main model) saved.")

# 2) Helper to compute metrics
def compute_metrics(y_true, y_pred, subset_name):
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    return {
        "subset": subset_name,
        "count": len(y_true),
        "rmse": rmse,
        "mse": mse,
        "mae": mae,
        "r2": r2,
    }

metrics_rows = []

# Quantiles on test set
q90_test = y_test_original.quantile(0.90)
q95_test = y_test_original.quantile(0.95)

mask_q90 = y_test_original >= q90_test
mask_q95 = y_test_original >= q95_test

# Main model – all / top10% / top5%
metrics_rows.append(
    compute_metrics(y_test_original, test_preds_main, "main_all")
)
metrics_rows.append(
    compute_metrics(
        y_test_original[mask_q90], test_preds_main[mask_q90], "main_PoI>=q90"
    )
)
metrics_rows.append(
    compute_metrics(
        y_test_original[mask_q95], test_preds_main[mask_q95], "main_PoI>=q95"
    )
)

# Expert model – only defined on high‑PoI subset (>= global 95% threshold)
metrics_rows.append(
    compute_metrics(
        y_test_hi_original, test_preds_expert, "expert_PoI>=global_q95"
    )
)

metrics_df = pd.DataFrame(metrics_rows)
metrics_df.to_csv(
    os.path.join(OUTPUT_DIR, "test_performance_metrics_detailed.csv"),
    index=False,
)
print("- Detailed test metrics (overall and high‑PoI subsets) saved.")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

print("- Generating and saving plots...")

# True vs predicted (main model, original scale)
plt.figure(figsize=(8, 8))
plt.scatter(y_test_original, test_preds_main, alpha=0.3)
lims = [min(plt.xlim()[0], plt.ylim()[0]), max(plt.xlim()[1], plt.ylim()[1])]
plt.plot(lims, lims, "r-", alpha=0.75, zorder=0, label="Perfect Prediction")
plt.xlabel("True Values (PoI)")
plt.ylabel("Predicted Values (PoI)")
plt.title("True vs. Predicted Values on Test Set (log-trained, tail-weighted)")
plt.legend()
plt.grid(True)
plt.savefig(
    os.path.join(OUTPUT_DIR, "true_vs_predicted_plot_log_model_tailweights.jpg"),
    format="jpg",
    dpi=300,
)
plt.close()
print("- True vs. Predicted plot (main model) saved.")

# Residuals histogram (main model, original scale)
residuals = y_test_original - test_preds_main
plt.figure(figsize=(10, 6))
sns.histplot(residuals, bins=50, kde=True)
plt.xlabel("Residual (True - Predicted)")
plt.ylabel("Frequency")
plt.title("Histogram of Prediction Residuals (original PoI scale, main model)")
plt.axvline(x=0, color="red", linestyle="--")
plt.grid(True)
plt.savefig(
    os.path.join(OUTPUT_DIR, "residuals_histogram_log_model_tailweights.jpg"),
    format="jpg",
    dpi=300,
)
plt.close()
print("- Residuals histogram (main model) saved.")

print("\n--- All regression tasks with tail weighting, expert model, "
      "and high‑PoI metrics completed successfully! ---")
