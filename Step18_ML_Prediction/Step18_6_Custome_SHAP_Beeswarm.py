import os
import pandas as pd
import numpy as np
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
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, ListedColormap  # <--- Added ListedColormap

# --- Optional: SHAP import (graceful fallback) ---
try:
    import shap
    HAVE_SHAP = True
except Exception as e:
    print(f"SHAP not available ({e}); SHAP-based explanations will be skipped.")
    HAVE_SHAP = False

# --- Configuration ---
POSITIVE_SET_PATH = "../../public/multiomics/hg38_bin_positive_histone.csv"
NEGATIVE_SET_PATH = "../../public/multiomics/hg38_bin_negative_histone.csv"
OUTPUT_DIR = "../../public/multiomics/Step18_ML_Prediction"

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
    raise SystemExit(1)

# --- Step 3: Export Combined Feature Matrix ---
print("\nStep 3: Exporting the combined feature matrix...")
combined_output_path = os.path.join(
    OUTPUT_DIR, "hg38_bin_feature_matrix_histone.csv"
)
hg38_bin_feature_matrix_histone.to_csv(combined_output_path, index=False)
print(f"Combined feature matrix saved to: {combined_output_path}")

# --- Step 4: Split Data into Training, Validation, and Test Sets ---
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

print(f"Training set shape:    {training_df.shape}")
print(f"Validation set shape: {validation_df.shape}")
print(f"Test set shape:        {test_df.shape}")

# --- Step 5, 6, & 7: Export Data Splits ---
print("\nSteps 5, 6, & 7: Exporting data splits...")
training_df.to_csv(
    os.path.join(OUTPUT_DIR, "hg38_bin_feature_matrix_histone_training.csv"),
    index=False,
)
validation_df.to_csv(
    os.path.join(OUTPUT_DIR, "hg38_bin_feature_matrix_histone_validation.csv"),
    index=False,
)
test_df.to_csv(
    os.path.join(OUTPUT_DIR, "hg38_bin_feature_matrix_histone_test.csv"), index=False
)
print("Training, validation, and test sets have been saved successfully.")

# --- Step 8: Model Training with XGBoost ---
print("\nStep 8: Training the XGBoost GBDT model...")

# --- UPDATE: Explicitly ordered features as requested ---
features = [
    "H3K4me1", "H3K4me2", "H3K4me3", "H3K9ac", "H3K9me3",
    "H3K27ac", "H3K27me3", "H3K36me3", "H4K20me1",
]
label = "PoI"

X_train, y_train = training_df[features], training_df[label]
X_val, y_val = validation_df[features], validation_df[label]
X_test, y_test = test_df[features], test_df[label]

dtrain = xgb.DMatrix(X_train, label=y_train)
dval = xgb.DMatrix(X_val, label=y_val)
dtest = xgb.DMatrix(X_test, label=y_test)

params = {
    "objective": "binary:logistic",
    "eta": 0.1,
    "max_depth": 6,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "seed": 42,
    "eval_metric": ["logloss", "error", "auc", "aucpr"],
}

evals_result = {}

model = xgb.train(
    params=params,
    dtrain=dtrain,
    num_boost_round=1000,
    evals=[(dtrain, "train"), (dval, "validation")],
    early_stopping_rounds=50,
    evals_result=evals_result,
    verbose_eval=100,
)

print("Model training completed.")
best_iteration = model.best_iteration

# --- Step 9 & 10: Export the Trained Model ---
print("\nSteps 9 & 10: Exporting the trained model object...")
model.save_model(os.path.join(OUTPUT_DIR, "hg38_bin_histone_trained_model.bin"))
model.save_model(os.path.join(OUTPUT_DIR, "hg38_bin_histone_trained_model.json"))
print("Model saved in both binary (.bin) and JSON (.json) formats.")

# --- Step 11: Export Training and Validation Performance Metrics ---
print("\nStep 11: Exporting training and validation performance metrics...")
train_metrics_df = pd.DataFrame(evals_result["train"])
train_metrics_df["set"] = "training"
val_metrics_df = pd.DataFrame(evals_result["validation"])
val_metrics_df["set"] = "validation"

performance_df = pd.concat([train_metrics_df, val_metrics_df], axis=0)
performance_df["epoch"] = performance_df.groupby("set").cumcount() + 1

# Validation predictions at best iteration
val_preds_proba = model.predict(dval, iteration_range=(0, best_iteration + 1))

# Choose decision threshold on validation to maximize F1
prec, rec, thresh = precision_recall_curve(y_val, val_preds_proba)
f1s = (2 * prec * rec) / np.clip(prec + rec, 1e-9, None)
valid_idx = np.arange(len(thresh)) 
best_idx = valid_idx[np.argmax(f1s[:-1])] if len(thresh) > 0 else 0
best_threshold = float(thresh[best_idx]) if len(thresh) > 0 else 0.5

val_preds_class = (val_preds_proba >= best_threshold).astype(int)
val_f1 = f1_score(y_val, val_preds_class)

print(f"Validation F1 Score at best iteration ({best_iteration}) with threshold {best_threshold:.4f}: {val_f1:.4f}")

performance_df.to_csv(
    os.path.join(OUTPUT_DIR, "training_validation_metrics.csv"), index=False
)

pd.DataFrame(
    [{"set": "validation", "metric": "f1_score", "value": val_f1, "threshold": best_threshold}]
).to_csv(os.path.join(OUTPUT_DIR, "validation_f1_score.csv"), index=False)

with open(os.path.join(OUTPUT_DIR, "selected_threshold.txt"), "w") as f:
    f.write(f"{best_threshold}\n")

# --- Step 12: Evaluate Model on the Test Set ---
print("\nStep 12: Evaluating the model on the unseen test set...")
test_preds_proba = model.predict(dtest, iteration_range=(0, best_iteration + 1))
test_preds_class = (test_preds_proba >= best_threshold).astype(int)
print("Predictions generated for the test set.")

# --- Step 13: Export Test Results ---
print("\nStep 13: Exporting test results...")

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
predictions_df.to_csv(os.path.join(OUTPUT_DIR, "test_predictions.csv"), index=False)
print("- Test predictions saved.")

test_metrics = {
    "logloss": log_loss(y_test, test_preds_proba),
    "error": 1 - accuracy_score(y_test, test_preds_class),
    "auc": roc_auc_score(y_test, test_preds_proba),
    "aucpr": average_precision_score(y_test, test_preds_proba),
    "f1": f1_score(y_test, test_preds_class),
    "threshold": best_threshold,
}
test_metrics_df = pd.DataFrame([test_metrics])
test_metrics_df.to_csv(os.path.join(OUTPUT_DIR, "test_performance_metrics.csv"), index=False)
print("- Test performance metrics saved.")

tn, fp, fn, tp = confusion_matrix(y_test, test_preds_class).ravel()
cm_df = pd.DataFrame(
    {
        "Metric": ["True Negative", "False Positive", "False Negative", "True Positive"],
        "Value": [tn, fp, fn, tp],
    }
)
cm_df.to_csv(os.path.join(OUTPUT_DIR, "test_confusion_matrix.csv"), index=False)
print("- Test confusion matrix saved.")

# Feature Importance
importance_types = ["weight", "gain", "cover"]
feature_importance_list = []
for imp_type in importance_types:
    scores = model.get_score(importance_type=imp_type)
    temp_df = pd.DataFrame({
        'feature': list(scores.keys()),
        'score': list(scores.values()),
        'type': imp_type
    })
    feature_importance_list.append(temp_df)

feature_importance_df = pd.concat(feature_importance_list, ignore_index=True)
feature_importance_df.to_csv(
    os.path.join(OUTPUT_DIR, "feature_importance_scores.csv"), index=False
)
print("- Feature importance scores saved.")

# Standard Plots
print("- Generating and saving plots...")

fpr, tpr, _ = roc_curve(y_test, test_preds_proba)
roc_auc = roc_auc_score(y_test, test_preds_proba)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.grid(True)
plt.savefig(os.path.join(OUTPUT_DIR, "roc_curve.jpg"), format='jpg', dpi=300)
plt.close()

precision, recall, _ = precision_recall_curve(y_test, test_preds_proba)
avg_precision = average_precision_score(y_test, test_preds_proba)

plt.figure(figsize=(8, 6))
plt.step(recall, precision, where='post', label=f'PR curve (AP = {avg_precision:.2f})')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.title('Precision-Recall Curve')
plt.legend(loc="upper right")
plt.grid(True)
plt.savefig(os.path.join(OUTPUT_DIR, "precision_recall_curve.jpg"), format='jpg', dpi=300)
plt.close()

cm = confusion_matrix(y_test, test_preds_class)
plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
            xticklabels=['Predicted Negative', 'Predicted Positive'],
            yticklabels=['Actual Negative', 'Actual Positive'])
plt.ylabel('Actual Label')
plt.xlabel('Predicted Label')
plt.title('Confusion Matrix Heatmap')
plt.savefig(os.path.join(OUTPUT_DIR, "confusion_matrix_heatmap.jpg"), format='jpg', dpi=300)
plt.close()

# --- INTERPRETABILITY ADD-ONS ---

if HAVE_SHAP:
    print("\nSHAP: Computing SHAP values for test set...")
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)

    # Save SHAP values
    shap_df = pd.DataFrame(shap_values, columns=features)
    shap_df.insert(0, "chr", test_df["chr"].values)
    shap_df.insert(1, "start", test_df["start"].values)
    shap_df.insert(2, "end", test_df["end"].values)
    shap_df.insert(3, "true_label", y_test.values)
    shap_df.insert(4, "predicted_probability", test_preds_proba)
    shap_df.to_csv(os.path.join(OUTPUT_DIR, "shap_values_test.csv"), index=False)
    print("- SHAP values saved.")

    # -------------------------------------------------------------------------
    # MODIFICATION: Global SHAP summary with Threshold Color & Manual Ordering
    # -------------------------------------------------------------------------
    
    # 1. Prepare Coloring Data: 
    #    Define threshold: -log10(0.05)
    threshold_val = -np.log10(0.05)
    
    #    Create a binary matrix for coloring. 
    #    1 (Red) if value >= threshold_val, 0 (Blue) if value < threshold_val.
    X_color = X_test.copy()
    X_color = (X_color >= threshold_val).astype(int)

    # 2. Define Discrete Colormap:
    #    0 -> Blue, 1 -> Red
    cmap_threshold = ListedColormap(["blue", "red"])

    # 3. Plot
    shap.summary_plot(
        shap_values, 
        X_color, 
        feature_names=features, 
        show=False, 
        sort=False, 
        cmap=cmap_threshold
    )
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "shap_summary_beeswarm.jpg"), format='jpg', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"- SHAP summary beeswarm saved (Color: Red >= {threshold_val:.2f} [-log10(0.05)], Blue < {threshold_val:.2f}).")

    # -------------------------------------------------------------------------
    # ADDITIONAL PLOT: Beeswarm showing ONLY Red points (>= threshold)
    # -------------------------------------------------------------------------
    print("\nGenerating additional SHAP plot (Red points only)...")

    # 1. Create a boolean mask for feature values below the threshold.
    #    X_test is a DataFrame; .values gives the underlying numpy matrix.
    mask_below_threshold = X_test.values < threshold_val

    # 2. Create a copy of SHAP values and set indices below threshold to NaN.
    #    Matplotlib scatter plots completely ignore points with NaN coordinates, 
    #    ensuring they are not plotted at all.
    shap_values_red_only = shap_values.copy()
    shap_values_red_only[mask_below_threshold] = np.nan

    # 3. Define a colormap that is solid red.
    #    Since we've effectively removed the "blue" data points using NaN, 
    #    everything remaining should be red. We define a simple all-red map.
    cmap_solid_red = ListedColormap(["red"])

    # 4. Plot
    #    We use the masked shap_values matrix.
    #    We pass X_color just to satisfy the function's API signature for coloring data,
    #    but the cmap decides the actual color (everything remaining becomes red).
    shap.summary_plot(
        shap_values_red_only,
        X_color,
        feature_names=features,
        show=False,
        sort=False, # Keep manual order
        cmap=cmap_solid_red,
        color_bar=False, # Hide color bar as everything is one color
        alpha=1.0        # Ensure remaining red dots are solid
    )

    plt.tight_layout()
    output_path_red = os.path.join(OUTPUT_DIR, "shap_summary_beeswarm_red_only.jpg")
    
    # Explicitly set facecolor='white' during save to ensure a clean JPG background
    plt.savefig(output_path_red, format='jpg', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"- SHAP summary beeswarm (Red points only >= {threshold_val:.2f}) saved to: {output_path_red}")

    # Global bar plot of mean |SHAP| (Standard)
    # Using sort=False here as well to match the Beeswarm order
    shap.summary_plot(shap_values, X_test, feature_names=features, plot_type="bar", show=False, sort=False)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "shap_summary_bar.jpg"), format='jpg', dpi=300, bbox_inches='tight')
    plt.close()
    print("- SHAP summary bar saved.")

    # Per-sample top contributions for first N misclassified
    N = 5
    misclassified_idx = predictions_df.index[predictions_df["true_label"] != predictions_df["predicted_label"]].tolist()
    top_k = 5
    rows = []
    for idx in misclassified_idx[:N]:
        contribs = shap_values[idx]
        vals = X_test.iloc[idx].values
        order = np.argsort(np.abs(contribs))[::-1][:top_k]
        for j in order:
            rows.append({
                "row_index": int(idx),
                "chr": str(test_df.iloc[idx]["chr"]),
                "start": int(test_df.iloc[idx]["start"]),
                "end": int(test_df.iloc[idx]["end"]),
                "feature": features[j],
                "feature_value": float(vals[j]),
                "shap_contribution": float(contribs[j]),
                "predicted_probability": float(test_preds_proba[idx]),
                "true_label": int(y_test.iloc[idx]),
                "predicted_label": int(test_preds_class[idx]),
            })
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(OUTPUT_DIR, "misclassified_top_shap_contributions.csv"), index=False)
        print("- Misclassified samples' top SHAP contributions saved.")
    else:
        print("- No misclassified samples to explain with SHAP top contributions.")

    # ------------------------------------------------------------------
    # SHAP dependence plots
    # ------------------------------------------------------------------
    from shap.plots import colors as shap_colors

    feature_value_ranges = {
        feat: (
            float(hg38_bin_feature_matrix_histone[feat].min()),
            float(hg38_bin_feature_matrix_histone[feat].max()),
        )
        for feat in features
    }

    def plot_shap_dependence_with_global_colors(
        main_feat, interact_feat, shap_values_matrix, X_frame, feature_ranges_dict, cmap=None, dot_size=16, alpha=0.8
    ):
        if cmap is None:
            cmap = shap_colors.red_blue

        main_idx = features.index(main_feat)
        x_vals = X_frame[main_feat].values
        y_vals = shap_values_matrix[:, main_idx]
        c_vals = X_frame[interact_feat].values

        order = np.argsort(c_vals)
        x_vals = x_vals[order]
        y_vals = y_vals[order]
        c_vals = c_vals[order]

        vmin, vmax = feature_ranges_dict[interact_feat]

        plt.figure()
        sc = plt.scatter(
            x_vals, y_vals, c=c_vals, cmap=cmap, vmin=vmin, vmax=vmax, s=dot_size, linewidth=0.0, alpha=alpha
        )
        plt.xlabel(main_feat)
        plt.ylabel(f"SHAP value for {main_feat}")
        cb = plt.colorbar(sc)
        cb.set_label(interact_feat)

    print("\nSHAP: Generating dependence plots for all feature combinations...")
    try:
        count = 0
        for i, main_feat in enumerate(features):
            for j, interact_feat in enumerate(features):
                if i == j:
                    continue
                plot_shap_dependence_with_global_colors(
                    main_feat=main_feat,
                    interact_feat=interact_feat,
                    shap_values_matrix=shap_values,
                    X_frame=X_test,
                    feature_ranges_dict=feature_value_ranges,
                    cmap=shap_colors.red_blue,
                    dot_size=16,
                    alpha=0.8,
                )
                plt.tight_layout()
                plt.savefig(
                    os.path.join(OUTPUT_DIR, f"shap_dependence_{main_feat}__{interact_feat}.jpg"),
                    format="jpg", dpi=300, bbox_inches="tight"
                )
                plt.close()
                count += 1
        print(f"- Saved {count} SHAP dependence plots.")
    except Exception as e:
        print(f"Skipping SHAP dependence plots due to: {e}")

# Single-tree visualization
try:
    xgb.plot_tree(model, num_trees=0, rankdir='LR')
    fig = plt.gcf()
    fig.set_size_inches(20, 10)
    plt.savefig(os.path.join(OUTPUT_DIR, "tree_visualization_0.jpg"), format='jpg', dpi=150, bbox_inches='tight')
    plt.close()
    
    tree_dump = model.get_dump(dump_format='text')
    with open(os.path.join(OUTPUT_DIR, "tree_dump_first5.txt"), 'w') as f:
        for i, tree in enumerate(tree_dump[:5]):
            f.write(f"Tree {i}:\n{tree}\n\n")
    print("- Tree visualization and text dump saved.")
except Exception as e:
    print(f"Tree visualization skipped due to: {e}")

# Lightweight Partial Dependence
print("\nComputing lightweight partial dependence curves...")
pdp_rows = []
rng = np.random.default_rng(42)
sample_idx = rng.choice(len(X_test), size=min(5000, len(X_test)), replace=False)
X_ref = X_test.iloc[sample_idx].copy()

for feat in features:
    grid = np.quantile(X_test[feat].values, np.linspace(0.01, 0.99, 50))
    for val in grid:
        X_tmp = X_ref.copy()
        X_tmp[feat] = val
        preds = model.predict(xgb.DMatrix(X_tmp), iteration_range=(0, best_iteration + 1))
        pdp_rows.append({
            "feature": feat,
            "grid_value": float(val),
            "avg_pred": float(np.mean(preds))
        })

pdp_df = pd.DataFrame(pdp_rows)
pdp_df.to_csv(os.path.join(OUTPUT_DIR, "pdp_curves.csv"), index=False)

ncols = 3
nrows = int(np.ceil(len(features) / ncols))
plt.figure(figsize=(5*ncols, 3.5*nrows))
for i, feat in enumerate(features, start=1):
    plt.subplot(nrows, ncols, i)
    sub = pdp_df[pdp_df["feature"] == feat].sort_values("grid_value")
    plt.plot(sub["grid_value"], sub["avg_pred"], lw=2)
    plt.title(feat)
    plt.xlabel("Value")
    plt.ylabel("Avg predicted P(positive)")
    plt.grid(True)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "pdp_grid.jpg"), format='jpg', dpi=300)
plt.close()
print("- PDP curves saved.")

print("\n--- All tasks completed successfully! ---")