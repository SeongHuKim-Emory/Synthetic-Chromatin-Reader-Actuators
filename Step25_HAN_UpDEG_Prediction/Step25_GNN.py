#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import random
from collections import defaultdict

import numpy as np
import pandas as pd

# Threading
os.environ["OMP_NUM_THREADS"] = "36"
os.environ["OPENBLAS_NUM_THREADS"] = "36"
os.environ["MKL_NUM_THREADS"] = "36"
os.environ["VECLIB_MAXIMUM_THREADS"] = "36"
os.environ["NUMEXPR_NUM_THREADS"] = "36"

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import AdamW
from torch.optim.lr_scheduler import ReduceLROnPlateau

from torch_geometric.data import HeteroData
from torch_geometric.nn import HANConv
from torch_geometric.loader import NeighborLoader

from sklearn.metrics import (
    precision_recall_curve, roc_curve, average_precision_score,
    roc_auc_score, f1_score
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Reproducibility
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)

# -------------------------------
# Paths and I/O
# -------------------------------
BASE_DIR = "../../public/Multiomics/Step25_HAN_UpDEG_Prediction"
NODES_CSV = os.path.join(BASE_DIR, "Nodes_with_gene_enhancer_promoter_pro_en.csv")
LABELS_CSV = os.path.join(BASE_DIR, "Labels.csv")
EDGES_CSV = os.path.join(BASE_DIR, "Edges.csv")
FEATS_CSV = os.path.join(BASE_DIR, "Features.csv")

OUT_DIR = BASE_DIR
os.makedirs(OUT_DIR, exist_ok=True)

# -------------------------------
# Device
# -------------------------------
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}, CUDA available: {torch.cuda.is_available()}")

# -------------------------------
# Load data
# -------------------------------
nodes_df = pd.read_csv(NODES_CSV)
labels_df = pd.read_csv(LABELS_CSV)
edges_df = pd.read_csv(EDGES_CSV)
feats_df = pd.read_csv(FEATS_CSV)

def sanitize_node_type(nt: str) -> str:
    return nt.replace("/", "_")

nodes_df["node_type_raw"] = nodes_df["node_type"]
nodes_df["node_type"] = nodes_df["node_type_raw"].apply(sanitize_node_type)

assert nodes_df["node_id"].is_monotonic_increasing and nodes_df["node_id"].iloc[0] == 0
num_nodes_total = nodes_df.shape[0]

node_types = sorted(nodes_df["node_type"].unique().tolist())
type_to_ids = {t: nodes_df.loc[nodes_df["node_type"] == t, "node_id"].to_numpy() for t in node_types}
type_id_maps = {t: {nid: i for i, nid in enumerate(type_to_ids[t])} for t in node_types}

feat_cols = [c for c in feats_df.columns if c != "node_id"]
feature_matrix_all = feats_df.set_index("node_id").reindex(range(num_nodes_total))[feat_cols].to_numpy(dtype=np.float32)

feat_mean = feature_matrix_all.mean(axis=0, keepdims=True)
feat_std = feature_matrix_all.std(axis=0, keepdims=True) + 1e-6
feature_matrix_all_std = (feature_matrix_all - feat_mean) / feat_std

features_by_type = {t: torch.tensor(feature_matrix_all_std[type_to_ids[t], :], dtype=torch.float32) for t in node_types}

labels_sparse = labels_df.set_index("node_id")["is_upregulated"]
labels_series = labels_sparse.reindex(range(num_nodes_total), fill_value=0)

def chr_to_ord(chr_str: str) -> int:
    s = str(chr_str).lower().replace("chr", "")
    if s == "x": return 23
    if s == "y": return 24
    if s in ("m", "mt"): return 25
    try:
        return int(s)
    except (ValueError, TypeError):
        return 26

nodes_df["chr_ord"] = nodes_df["chr"].apply(chr_to_ord)

gene_type = "gene"
gene_ids = type_to_ids.get(gene_type, np.array([], dtype=int))
gene_chr_ord = nodes_df.loc[gene_ids, "chr_ord"].to_numpy()
gene_labels = labels_series.loc[gene_ids].to_numpy(dtype=np.int64)

# Chromosome-based split
train_mask_gene = torch.tensor((gene_chr_ord >= 1) & (gene_chr_ord <= 13), dtype=torch.bool)
val_mask_gene   = torch.tensor((gene_chr_ord >= 14) & (gene_chr_ord <= 17), dtype=torch.bool)
test_mask_gene  = torch.tensor((gene_chr_ord >= 18) & (gene_chr_ord <= 24), dtype=torch.bool)

# -------------------------------
# Build HeteroData graph
# -------------------------------
data = HeteroData()
for t in node_types:
    data[t].x = features_by_type[t]
    data[t].orig_id = torch.tensor(type_to_ids[t], dtype=torch.long)

data[gene_type].y = torch.tensor(gene_labels, dtype=torch.long)
data[gene_type].train_mask = train_mask_gene
data[gene_type].val_mask = val_mask_gene
data[gene_type].test_mask = test_mask_gene

edge_lists = defaultdict(list)
node_id_to_type = nodes_df.set_index("node_id")["node_type"].to_dict()

for _, row in edges_df.iterrows():
    g_id, r_id = int(row["node_id_1"]), int(row["node_id_2"])
    g_type, r_type = node_id_to_type.get(g_id), node_id_to_type.get(r_id)
    if g_type != gene_type or r_type not in {"Enhancer", "Promoter", "Promoter_Enhancer"}:
        continue
    edge_type_map = {
        "is_gene_Enhancer": (gene_type, "gene_to_Enhancer", "Enhancer"),
        "is_gene_Promoter": (gene_type, "gene_to_Promoter", "Promoter"),
        "is_gene_Promoter_Enhancer": (gene_type, "gene_to_Promoter_Enhancer", "Promoter_Enhancer"),
    }
    for bool_col, et in edge_type_map.items():
        if int(row[bool_col]) == 1 and et[2] == r_type:
            g_local, r_local = type_id_maps[g_type][g_id], type_id_maps[r_type][r_id]
            edge_lists[et].append((g_local, r_local))
            rev_et = (et[2], f"{et[2]}_to_{et[0]}", et[0])
            edge_lists[rev_et].append((r_local, g_local))

for et, pairs in edge_lists.items():
    if pairs:
        src, dst = zip(*pairs)
        data[et].edge_index = torch.tensor([src, dst], dtype=torch.long)
    else:
        data[et].edge_index = torch.empty((2, 0), dtype=torch.long)

data = data.to(device)
edge_index_dict_full = data.edge_index_dict

# -------------------------------
# Report class counts
# -------------------------------
y_train_all = data[gene_type].y[data[gene_type].train_mask].float()
num_pos_train = int(y_train_all.sum().item())
num_neg_train = int(y_train_all.numel() - y_train_all.sum().item())
print(f"[Train split] Positives: {num_pos_train}, Negatives: {num_neg_train}")

# -------------------------------
# Model
# -------------------------------
class HAN2LayerResidual(nn.Module):
    def __init__(self, metadata, in_channels_dict, hidden_dim=256, heads=2, dropout=0.5):
        super().__init__()
        # HANConv supports returning semantic attention weights. [docs]
        self.han1 = HANConv(in_channels=in_channels_dict, out_channels=hidden_dim,
                            metadata=metadata, heads=heads, dropout=dropout)
        self.han2 = HANConv(in_channels=hidden_dim, out_channels=hidden_dim,
                            metadata=metadata, heads=heads, dropout=dropout)
        self.ln1 = nn.LayerNorm(hidden_dim)
        self.ln2 = nn.LayerNorm(hidden_dim)
        gene_in = in_channels_dict["gene"]
        self.gene_proj = nn.Sequential(
            nn.Linear(gene_in, hidden_dim),
            nn.ReLU(inplace=True),
            nn.Dropout(p=dropout)
        )
        self.classifier = nn.Linear(hidden_dim, 1)
        self.dropout = nn.Dropout(p=dropout)

    def forward(self, x_dict, edge_index_dict, return_attn=False):
        h1_dict, attn1 = self.han1(x_dict, edge_index_dict, return_semantic_attention_weights=True)
        h1_gene = self.ln1(F.elu(h1_dict[gene_type]))
        h2_dict, attn2 = self.han2(h1_dict, edge_index_dict, return_semantic_attention_weights=True)
        h2_gene = self.ln2(F.elu(h2_dict[gene_type]))
        res_gene = self.gene_proj(x_dict[gene_type])
        fused = self.dropout(h2_gene + res_gene)
        out = self.classifier(fused).squeeze(-1)
        if return_attn:
            return out, attn2
        return out

metadata = data.metadata()
in_channels_dict = {nt: int(data[nt].x.size(1)) for nt in node_types}
model = HAN2LayerResidual(metadata, in_channels_dict, hidden_dim=256, heads=2, dropout=0.5).to(device)  # HANConv supports this usage. [docs]

optimizer = AdamW(model.parameters(), lr=3e-4, weight_decay=1e-4)
scheduler = ReduceLROnPlateau(optimizer, mode="max", patience=10, factor=0.5)

scaler = torch.amp.GradScaler(enabled=(device.type == "cuda"))

# -------------------------------
# Losses
# -------------------------------
LOSS_TYPE = "focal"  # "focal", "bce_cb", or "bce"
FOCAL_ALPHA = 0.75   # α applied to positives; 1-α to negatives. [focal]
FOCAL_GAMMA = 2.0    # γ controls focus on hard examples. [focal]

class FocalLoss(nn.Module):
    # Correct binary focal loss with per-class alpha_t. [focal]
    def __init__(self, alpha=0.25, gamma=2.0, eps=1e-8):
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.eps = eps
    def forward(self, logits, targets):
        targets = targets.float()
        bce = F.binary_cross_entropy_with_logits(logits, targets, reduction='none')
        # p_t = prob of true class
        p_t = torch.exp(-bce)
        # alpha_t per sample: α for positives, (1-α) for negatives
        alpha_t = self.alpha * targets + (1.0 - self.alpha) * (1.0 - targets)
        loss = alpha_t * (1.0 - p_t).pow(self.gamma) * bce
        return loss.mean()

def get_cb_pos_weight(neg_count, pos_count, beta=0.9999):
    # Effective number weighting for BCE pos_weight. [CB Loss]
    eff_neg = (1.0 - beta) / (1.0 - beta ** max(1, neg_count))
    eff_pos = (1.0 - beta) / (1.0 - beta ** max(1, pos_count))
    return torch.tensor([eff_pos / eff_neg], dtype=torch.float32, device=device)

if LOSS_TYPE == "focal":
    criterion = FocalLoss(alpha=FOCAL_ALPHA, gamma=FOCAL_GAMMA)  # Focus on hard positives. [focal]
elif LOSS_TYPE == "bce_cb":
    pos_weight_cb = get_cb_pos_weight(num_neg_train, num_pos_train, beta=0.9999)
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight_cb)  # BCEWithLogits supports pos_weight. [docs]
else:
    criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([1.0], device=device))  # Baseline BCE. [docs]

# -------------------------------
# Balanced seed sampler
# -------------------------------
train_gene_local_idx = torch.arange(data[gene_type].x.size(0), device=device)[data[gene_type].train_mask]
train_gene_labels = data[gene_type].y[data[gene_type].train_mask].cpu().numpy()

pos_mask_np = (train_gene_labels == 1)
neg_mask_np = (train_gene_labels == 0)

train_pos_local = train_gene_local_idx[pos_mask_np]
train_neg_local = train_gene_local_idx[neg_mask_np]

print(f"Train class sizes (local): pos={len(train_pos_local)}, neg={len(train_neg_local)}")

# Use 64 positives and 1:2 ratio to reduce over-repetition of the small positive pool.
SEEDS_POS = 64
POS_NEG_RATIO = 2
SEEDS_NEG = SEEDS_POS * POS_NEG_RATIO
SEEDS_PER_BATCH = SEEDS_POS + SEEDS_NEG

# Two hops for 2-layer HAN; modest fanout to keep context.
NUM_NEIGHBORS = [15, 10]
BATCHES_PER_EPOCH = max(1, int(np.ceil((len(train_pos_local) + len(train_neg_local)) / SEEDS_PER_BATCH)))

def draw_balanced_train_gene_seeds():
    pos = train_pos_local[torch.randint(low=0, high=len(train_pos_local), size=(SEEDS_POS,), device=device)]
    neg = train_neg_local[torch.randint(low=0, high=len(train_neg_local), size=(SEEDS_NEG,), device=device)]
    seeds = torch.cat([pos, neg], dim=0)
    perm = torch.randperm(seeds.size(0), device=device)
    seeds = seeds[perm]
    return seeds

def make_train_loader_for_gene_seeds(seeds_gene_local: torch.Tensor):
    # For hetero graphs, batch[gene_type].batch_size equals number of seeds and seeds are first. [NeighborLoader docs]
    loader = NeighborLoader(
        data,
        num_neighbors=NUM_NEIGHBORS,
        input_nodes=(gene_type, seeds_gene_local),
        batch_size=seeds_gene_local.numel(),
        shuffle=False,
        directed=True
    )
    return loader

# -------------------------------
# Metrics and thresholding
# -------------------------------
def compute_metrics_from_logits(logits, y_true, mask):
    with torch.no_grad():
        logits_m, y_m = logits[mask], y_true[mask].float()
        if logits_m.numel() == 0:
            return {"loss": 0, "auroc": 0.5, "auprc": 0.0, "acc": 0.0}
        probs = torch.sigmoid(logits_m)
        metric_loss = F.binary_cross_entropy_with_logits(logits_m, y_m).item()
        y_cpu = y_m.cpu().numpy()
        p_cpu = probs.cpu().numpy()
        auroc = roc_auc_score(y_cpu, p_cpu) if len(np.unique(y_cpu)) > 1 else 0.5
        auprc = average_precision_score(y_cpu, p_cpu) if y_cpu.sum() > 0 else 0.0
        acc = (((probs > 0.5).long() == y_m.long()).float().mean().item())
        return {"loss": metric_loss, "auroc": auroc, "auprc": auprc, "acc": acc}

def pick_threshold_for_pr_recall_bias(logits):
    # Favor recall for AUCPR improvements via F2 selection. [PR guidance]
    y_val = data[gene_type].y[data[gene_type].val_mask].float().cpu().numpy()
    p_val = torch.sigmoid(logits[data[gene_type].val_mask]).cpu().numpy()
    if y_val.sum() == 0:
        return 0.5, {"val_f2": 0.0, "val_prec": 0.0, "val_rec": 0.0, "best_threshold": 0.5}
    prec, rec, thr = precision_recall_curve(y_val, p_val)
    f2 = (5 * prec * rec) / np.maximum(4 * prec + rec, 1e-9)
    best_idx = int(np.nanargmax(f2[:-1])) if len(thr) > 0 else 0
    best_thr = thr[best_idx] if len(thr) > 0 else 0.5
    info = {
        "val_f2": float(f2[best_idx]),
        "val_prec": float(prec[best_idx]),
        "val_rec": float(rec[best_idx]),
        "best_threshold": float(best_thr),
    }
    return best_thr, info

# -------------------------------
# Training loop
# -------------------------------
EPOCHS = 80
best_val_auprc = -1.0
best_state = None
best_epoch = -1
epochs_no_improve = 0
early_stop_patience = 20
metrics_history = []

print("Starting training with corrected focal loss and balanced seed sampling (1:2)...")
for epoch in range(1, EPOCHS + 1):
    model.train()
    epoch_losses = []
    for _ in range(BATCHES_PER_EPOCH):
        seeds = draw_balanced_train_gene_seeds()
        loader = make_train_loader_for_gene_seeds(seeds)
        batch = next(iter(loader)).to(device)

        seed_count = seeds.numel()
        # NeighborLoader guarantees seed nodes are first and exposes batch_size on the node-type store. [NeighborLoader docs]
        assert getattr(batch[gene_type], "batch_size", seed_count) == seed_count
        seed_mask_local = torch.zeros(batch[gene_type].x.size(0), dtype=torch.bool, device=device)
        seed_mask_local[:seed_count] = True

        optimizer.zero_grad(set_to_none=True)
        with torch.amp.autocast(device_type=device.type, enabled=(device.type == "cuda")):
            logits_batch = model(batch.x_dict, batch.edge_index_dict)
            y_seed = batch[gene_type].y[seed_mask_local].float()
            logits_seed = logits_batch[seed_mask_local]
            loss = criterion(logits_seed, y_seed)

        scaler.scale(loss).backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=2.0)
        scaler.step(optimizer)
        scaler.update()
        epoch_losses.append(loss.item())

    # Full-graph evaluation
    model.eval()
    with torch.no_grad(), torch.amp.autocast(device_type=device.type, enabled=(device.type == "cuda")):
        full_logits = model(data.x_dict, edge_index_dict_full)

    train_m = compute_metrics_from_logits(full_logits, data[gene_type].y, data[gene_type].train_mask)
    val_m   = compute_metrics_from_logits(full_logits, data[gene_type].y, data[gene_type].val_mask)

    scheduler.step(val_m["auprc"])

    row = {f"train_{k}": v for k, v in train_m.items()}
    row.update({f"val_{k}": v for k, v in val_m.items()})
    row["epoch"] = epoch
    row["lr"] = optimizer.param_groups[0]['lr']
    row["train_epoch_loss_mean"] = float(np.mean(epoch_losses)) if epoch_losses else 0.0
    metrics_history.append(row)

    if epoch % 10 == 0:
        print(f"[Epoch {epoch:03d}] "
              f"Train AUROC: {train_m['auroc']:.4f}, Val AUROC: {val_m['auroc']:.4f}, "
              f"Train AUCPR: {train_m['auprc']:.4f}, Val AUCPR: {val_m['auprc']:.4f}, "
              f"LR: {row['lr']:.6f}, Train Loss: {row['train_epoch_loss_mean']:.4f}")

    improved = val_m["auprc"] > best_val_auprc
    if improved:
        best_val_auprc = val_m["auprc"]
        best_state = model.state_dict()
        best_epoch = epoch
        epochs_no_improve = 0
        torch.save(best_state, os.path.join(OUT_DIR, "han_best.ckpt"))
    else:
        epochs_no_improve += 1
        if epochs_no_improve >= early_stop_patience:
            print(f"Early stopping at epoch {epoch}, best val AUCPR at epoch {best_epoch}: {best_val_auprc:.4f}")
            break

pd.DataFrame(metrics_history).to_csv(os.path.join(OUT_DIR, "train_val_metrics.csv"), index=False)
print("Training complete. Best validation AUCPR:", best_val_auprc)

# -------------------------------
# Final evaluation
# -------------------------------
if best_state is not None:
    model.load_state_dict(best_state)
model.eval()

with torch.no_grad(), torch.amp.autocast(device_type=device.type, enabled=(device.type == "cuda")):
    test_logits, attn_dict = model(data.x_dict, edge_index_dict_full, return_attn=True)

print("Selecting threshold on validation PR with recall bias (F2)...")
best_thr, thr_info = pick_threshold_for_pr_recall_bias(test_logits)

def export_curves_and_metrics(split_name, mask, logits, threshold):
    y_true = data[gene_type].y[mask].cpu().numpy()
    probs = torch.sigmoid(logits[mask]).cpu().numpy()

    auroc = roc_auc_score(y_true, probs) if len(np.unique(y_true)) > 1 else 0.5
    auprc = average_precision_score(y_true, probs) if y_true.sum() > 0 else 0.0
    preds = (probs >= threshold).astype(np.int32)
    f1 = f1_score(y_true, preds, zero_division=0)
    prec_curve, rec_curve, _ = precision_recall_curve(y_true, probs)
    if len(np.unique(y_true)) > 1:
        fpr, tpr, _ = roc_curve(y_true, probs)
    else:
        fpr, tpr = np.array([0, 1]), np.array([0, 1])

    # Plots
    plt.figure()
    plt.plot(fpr, tpr, label=f"{split_name} ROC AUC = {auroc:.3f}")
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
    plt.xlabel("False Positive Rate"); plt.ylabel("True Positive Rate"); plt.title(f"{split_name} ROC Curve")
    plt.legend(); plt.savefig(os.path.join(OUT_DIR, f"{split_name.lower()}_roc_curve.jpg")); plt.close()

    plt.figure()
    plt.plot(rec_curve, prec_curve, label=f"{split_name} AUCPR = {auprc:.3f}")
    plt.xlabel("Recall"); plt.ylabel("Precision"); plt.title(f"{split_name} Precision-Recall Curve")
    plt.legend(); plt.savefig(os.path.join(OUT_DIR, f"{split_name.lower()}_pr_curve.jpg")); plt.close()

    metrics = {
        "auroc": float(auroc),
        "auprc": float(auprc),
        "f1_at_thr": float(f1),
        "threshold": float(threshold),
        "positives": int(y_true.sum()),
        "total": int(len(y_true)),
    }
    with open(os.path.join(OUT_DIR, f"{split_name.lower()}_metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)

    # Export predictions
    test_gene_ids = data[gene_type].orig_id[mask].cpu().numpy()
    nodes_gene_df = nodes_df.set_index("node_id").loc[test_gene_ids].reset_index()
    pred_df = pd.DataFrame({"node_id": test_gene_ids, "y_true": y_true, "y_prob": probs, "y_pred": preds})
    pred_df = nodes_gene_df.merge(pred_df, on="node_id")
    pred_df.to_csv(os.path.join(OUT_DIR, f"{split_name.lower()}_predictions.csv"), index=False)

    return metrics

with open(os.path.join(OUT_DIR, "threshold_selection.json"), "w") as f:
    json.dump(thr_info, f, indent=2)

print("Evaluating on validation and test with the selected threshold...")
val_metrics = export_curves_and_metrics("Val", data[gene_type].val_mask, test_logits, best_thr)
print("Validation metrics:", val_metrics)
test_metrics = export_curves_and_metrics("Test", data[gene_type].test_mask, test_logits, best_thr)
print("Test metrics:", test_metrics)

# Export semantic attention
with open(os.path.join(OUT_DIR, "semantic_attention.json"), "w") as f:
    attn_export = {}
    if attn_dict:
        for node_type, weights in attn_dict.items():
            if weights is not None:
                attn_export[node_type] = weights.detach().cpu().numpy().tolist()
    json.dump(attn_export, f, indent=2)

print("\nScript finished. Outputs saved to:", OUT_DIR)

