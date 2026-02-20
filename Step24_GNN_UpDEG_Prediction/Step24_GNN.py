"""
GNN for Upregulated DEG Prediction
[REFACTOR] 3-Hop Version with Fixed Train/Val/Test Split
[REFACTOR] Using [3, 3, 3] sampling for 3-edge gene nodes

[HAN IMPLEMENTATION]
This script is refactored to use a Heterogeneous Graph Attention Network (HAN).
1.  Sets up the environment.
2.  Loads all data and constructs a `HeteroData` object, splitting nodes
    by type (gene, enhancer, etc.) and edges by relation (is_loop, etc.).
3.  Creates fixed 'train_mask', 'val_mask', 'test_mask' on the 'gene' node type
    based on chromosome splits.
4.  Implements a single training run with a `HANModel`.
    - Trains on the train set (chr1-13).
    - Uses the validation set (chr14-17) for early stopping.
5.  Evaluates the best model on the held-out test set (chr18-Y).
6.  Exports test metrics, plots, and node-level predictions for 'gene' nodes.
"""

import os
import logging
import polars as pl
import numpy as np
import torch
import torch.nn.functional as F
import matplotlib
matplotlib.use('Agg') # Use a non-interactive backend
import matplotlib.pyplot as plt
from torch_geometric.data import HeteroData
from torch_geometric.loader import NeighborLoader
# [HAN] Import HANConv and new modules
from torch_geometric.nn import HANConv, Linear
from torch.nn import CrossEntropyLoss
from torch.optim import AdamW
import torch.multiprocessing  # <-- [FIX] Import multiprocessing
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    average_precision_score,
    roc_curve,
    precision_recall_curve
)
from tqdm import tqdm

# --- 1. Environment Setup ---

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
log = logging.getLogger(__name__)

# Set CPU threads
num_cpus = 25
torch.set_num_threads(num_cpus)
os.environ['OMP_NUM_THREADS'] = str(num_cpus)
os.environ['MKL_NUM_THREADS'] = str(num_cpus)
log.info(f"Set PyTorch/OMP/MKL threads to {num_cpus}")

# --- [FIX] Set multiprocessing sharing strategy ---
# This is to avoid the "Too many open files" error (errno 24) when using
# num_workers > 0 with a very complex HeteroData object (many edge types).
# 'file_system' is slower than the default 'file_descriptor' but
# avoids the OS limit on open file descriptors.
try:
    torch.multiprocessing.set_sharing_strategy('file_system')
    log.info("Set torch multiprocessing sharing strategy to 'file_system'")
except RuntimeError as e:
    log.warning(f"Could not set sharing strategy: {e}. This might be normal on some systems (e.g., Windows).")
# -------------------------------------------------

# Check for CUDA
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
log.info(f"Using device: {device}")
if device.type == 'cuda':
    log.info(f"CUDA Device Name: {torch.cuda.get_device_name(0)}")

# --- 2. File Paths and Parameters ---

# Input file paths
# [HAN] Update base path if needed. Assuming relative path from script location.
BASE_PATH = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/' # Assuming files are in the same directory as the script
NODE_FILE = os.path.join(BASE_PATH, 'Nodes.csv')
LABEL_FILE = os.path.join(BASE_PATH, 'Labels.csv')
EDGE_FILE = os.path.join(BASE_PATH, 'Edges.csv')
FEATURE_FILE = os.path.join(BASE_PATH, 'Features.csv')

# [REFACTOR] Changed output dir for fixed split 3-hop experiment
OUTPUT_DIR = os.path.join(BASE_PATH, 'results_han_fixed_split_3layer_3hop_3_3_3') # Updated dir name
os.makedirs(OUTPUT_DIR, exist_ok=True)
PREDICTION_FILE = os.path.join(OUTPUT_DIR, 'final_test_predictions.csv')
TEST_METRICS_FILE = os.path.join(OUTPUT_DIR, 'final_test_metrics_summary.csv')
EPOCH_METRICS_FILE = os.path.join(OUTPUT_DIR, 'epoch_metrics.csv')
MODEL_FILE = os.path.join(OUTPUT_DIR, 'best_model.pth')
ROC_PLOT_FILE = os.path.join(OUTPUT_DIR, 'final_test_roc_curve.png')
PR_PLOT_FILE = os.path.join(OUTPUT_DIR, 'final_test_pr_curve.png')

# Model and Training Hyperparameters
HIDDEN_DIM = 128
NUM_HEADS = 4
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 5e-4
EPOCHS = 300
PATIENCE = 30

# Mini-batching Hyperparameters
BATCH_SIZE = 1024 # <-- CHANGED: Increased from 512 for better GPU utilization
# [REFACTOR] Changed to 3-hop sampling [15, 10, 5]
# <-- CHANGED: Switched to [3, 3, 3] to match file-level comments/path and drastically reduce sampling overhead
NUM_NEIGHBORS = [-1, -1, -1] # Hops: 1-hop, 2-hop, 3-hop
NUM_WORKERS = 20 # <-- CHANGED: Increased from 4 to utilize more of the 25 available CPUs
log.info(f"Using mini-batch size: {BATCH_SIZE}")
log.info(f"Using 3-HOP neighbor sampling: {NUM_NEIGHBORS}")
log.info(f"Using {NUM_WORKERS} data loading workers") # <-- ADDED: Log new worker count

# --- 3. Data Loading and Graph Construction ([HAN] Refactored) ---

def load_graph_data():
    """
    [HAN] Loads all CSV files and constructs a PyTorch Geometric HeteroData object.
    - Nodes are split by 'node_type'.
    - Edges are split by their type (is_loop, is_overlap, etc.).
    - Creates train/val/test masks on the 'gene' node type.
    """
    try:
        log.info("Loading data from CSV files...")
        nodes_df = pl.read_csv(NODE_FILE)
        labels_df = pl.read_csv(LABEL_FILE)
        edges_df = pl.read_csv(EDGE_FILE)
        features_df = pl.read_csv(FEATURE_FILE)
        log.info(f"Nodes.csv: {nodes_df.shape}")
        log.info(f"Labels.csv: {labels_df.shape}")
        log.info(f"Edges.csv: {edges_df.shape}")
        log.info(f"Features.csv: {features_df.shape}")
    except Exception as e:
        log.error(f"An error occurred during file loading: {e}")
        return None, None, None

    data = HeteroData()
    log.info("Constructing HeteroData object...")

    # --- 3.1 Process Nodes and Features ---
    # Map global node IDs to their type
    global_id_to_type = dict(zip(nodes_df['node_id'], nodes_df['node_type']))
    node_types = nodes_df['node_type'].unique().to_list()
    
    # Store feature names (from TPM onwards)
    feature_cols = features_df.columns[6:]
    log.info(f"Using {len(feature_cols)} node features: {feature_cols}")
    
    global_to_local_maps = {}
    gene_global_ids = None # To store original IDs for prediction mapping

    for node_type in tqdm(node_types, desc="Processing node types"):
        type_nodes_df = nodes_df.filter(pl.col('node_type') == node_type)
        type_node_ids_global = type_nodes_df['node_id'].to_numpy()
        
        # Create mapping from global ID to local type-specific ID (0 to N-1)
        global_to_local_maps[node_type] = {global_id: local_id for local_id, global_id in enumerate(type_node_ids_global)}
        
        # Get features
        type_features_df = features_df.filter(pl.col('node_id').is_in(type_node_ids_global)).sort('node_id')
        data[node_type].x = torch.tensor(type_features_df.select(feature_cols).to_numpy(), dtype=torch.float)
        
        data[node_type].num_nodes = len(type_node_ids_global)

        # --- 3.2 Process Labels and Masks (only for 'gene') ---
        if node_type == 'gene':
            gene_global_ids = type_node_ids_global # Save for later
            labels_map = dict(zip(labels_df['node_id'], labels_df['is_upregulated']))
            y = [labels_map.get(gid, 0) for gid in type_node_ids_global]
            data['gene'].y = torch.tensor(y, dtype=torch.long)
            
            # Create train/val/test masks based on chromosomes
            log.info("Creating train/val/test masks for 'gene' nodes...")
            train_chroms = [f'chr{i}' for i in range(1, 14)] # chr1-13
            val_chroms = [f'chr{i}' for i in range(14, 18)] # chr14-17
            test_chroms = [f'chr{i}' for i in range(18, 23)] + ['chrX', 'chrY'] # chr18-22, X, Y
            
            chrom_map = dict(zip(type_nodes_df['node_id'], type_nodes_df['chr']))
            
            train_mask = [chrom_map.get(gid) in train_chroms for gid in type_node_ids_global]
            val_mask = [chrom_map.get(gid) in val_chroms for gid in type_node_ids_global]
            test_mask = [chrom_map.get(gid) in test_chroms for gid in type_node_ids_global]
            
            data['gene'].train_mask = torch.tensor(train_mask, dtype=torch.bool)
            data['gene'].val_mask = torch.tensor(val_mask, dtype=torch.bool)
            data['gene'].test_mask = torch.tensor(test_mask, dtype=torch.bool)
        else:
            # Add placeholder masks for other types
            data[node_type].y = torch.zeros(len(type_node_ids_global), dtype=torch.long)
            data[node_type].train_mask = torch.zeros(len(type_node_ids_global), dtype=torch.bool)
            data[node_type].val_mask = torch.zeros(len(type_node_ids_global), dtype=torch.bool)
            data[node_type].test_mask = torch.zeros(len(type_node_ids_global), dtype=torch.bool)

    # --- 3.3 Process Edges ---
    log.info("Processing edges into relations...")
    # Get edge type feature columns
    edge_type_cols = edges_df.columns[3:]
    
    # Use numpy for faster iteration
    edges_array = edges_df.select(['node_id_source', 'node_id_target'] + edge_type_cols).to_numpy()
    
    edge_index_dict = {}
    
    for row in tqdm(edges_array, desc="Building edge relations"):
        src_global, tgt_global = int(row[0]), int(row[1])
        edge_features = row[2:]
        
        # Find node types and local IDs
        src_type = global_id_to_type[src_global]
        tgt_type = global_id_to_type[tgt_global]
        src_local = global_to_local_maps[src_type][src_global]
        tgt_local = global_to_local_maps[tgt_type][tgt_global]
        
        # Iterate over edge types (is_loop, is_overlap, etc.)
        for i, edge_type_name in enumerate(edge_type_cols):
            if edge_features[i] == 1:
                relation_tuple = (src_type, edge_type_name, tgt_type)
                
                if relation_tuple not in edge_index_dict:
                    edge_index_dict[relation_tuple] = [[], []]
                
                edge_index_dict[relation_tuple][0].append(src_local)
                edge_index_dict[relation_tuple][1].append(tgt_local)

    # Add edge indices to HeteroData object
    for relation_tuple, (src_list, tgt_list) in edge_index_dict.items():
        data[relation_tuple].edge_index = torch.tensor([src_list, tgt_list], dtype=torch.long)

    log.info(f"Graph construction complete:")
    log.info(f"\n{data}")
    log.info(f"  Num 'gene' nodes: {data['gene'].num_nodes}")
    log.info(f"  Training nodes (genes): {data['gene'].train_mask.sum().item()} (Chr 1-13)")
    log.info(f"  Validation nodes (genes): {data['gene'].val_mask.sum().item()} (Chr 14-17)")
    log.info(f"  Test nodes (genes): {data['gene'].test_mask.sum().item()} (Chr 18-Y)")
    
    # 'is_TAD' might be empty, which is fine
    data = data.coalesce()
    log.info(f"\nCoalesced graph:\n{data}")

    return data, nodes_df, gene_global_ids

# --- 4. GNN Model Definition ([HAN] New Model) ---

class HANModel(torch.nn.Module):
    """
    Heterogeneous Graph Attention Network (HAN)
    - Uses three HANConv layers.
    - Applies ELU activation and dropout.
    - Final linear layer for 2-class prediction (only on 'gene' nodes).
    
    [FIX] HANConv output shape is [N, out_channels], (heads are aggregated internally).
    """
    def __init__(self, in_channels, hidden_channels, out_channels, heads, metadata):
        super(HANModel, self).__init__()
        # HANConv requires in_channels to be an int (if all node types have
        # same features) or a dict. Here it's an int.
        
        # conv1: 11 -> 128
        self.conv1 = HANConv(
            in_channels, hidden_channels, heads=heads,
            dropout=0.6, metadata=metadata
        )
        # [FIX] conv2 in_channels is hidden_channels (128), not (hidden * heads)
        # conv2: 128 -> 128
        self.conv2 = HANConv(
            hidden_channels, hidden_channels, heads=heads,
            dropout=0.6, metadata=metadata
        )
        # --- [NEW 3rd LAYER] ---
        # conv3: 128 -> 128
        self.conv3 = HANConv(
            hidden_channels, hidden_channels, heads=heads,
            dropout=0.6, metadata=metadata
        )
        # -------------------------
        
        # [FIX] Final linear layer in_channels is hidden_channels (128)
        self.lin = Linear(hidden_channels, out_channels)

    def forward(self, data):
        # data is a HeteroDataBatch object
        x_dict, edge_index_dict = data.x_dict, data.edge_index_dict
        
        # First HAN layer
        # Input x_dict has tensors of shape [N, 11]
        # [FIX] Output x_dict has tensors of shape [N, 128]
        x_dict = self.conv1(x_dict, edge_index_dict)
        
        # [FIX] Apply ELU and Dropout directly (no view/reshape)
        x_dict = {key: F.elu(x) for key, x in x_dict.items()}
        x_dict = {key: F.dropout(x, p=0.6, training=self.training) for key, x in x_dict.items()}
        
        # Second HAN layer
        # [FIX] Input x_dict has tensors of shape [N, 128]
        # [FIX] Output x_dict has tensors of shape [N, 128]
        x_dict = self.conv2(x_dict, edge_index_dict)

        # [FIX] Apply ELU and Dropout directly
        x_dict = {key: F.elu(x) for key, x in x_dict.items()}
        x_dict = {key: F.dropout(x, p=0.6, training=self.training) for key, x in x_dict.items()}
        
        # --- [NEW 3rd LAYER BLOCK] ---
        # Third HAN layer
        # Input x_dict has tensors of shape [N, 128]
        # Output x_dict has tensors of shape [N, 128]
        x_dict = self.conv3(x_dict, edge_index_dict)

        # Apply ELU and Dropout directly
        x_dict = {key: F.elu(x) for key, x in x_dict.items()}
        x_dict = {key: F.dropout(x, p=0.6, training=self.training) for key, x in x_dict.items()}
        # -----------------------------
        
        # We only care about predictions for 'gene' nodes
        # [FIX] self.lin in_features is 128.
        # This computes logits for ALL 'gene' nodes in the subgraph (seed + neighbors)
        x_gene_all = self.lin(x_dict['gene'])
        
        # [FIX for ValueError] We must only return the logits for the
        # 'gene' seed nodes, which are the first 'batch_size' nodes.
        x_gene_seed = x_gene_all[:data['gene'].batch_size]
        
        return x_gene_seed

# --- 5. Training and Evaluation Functions ---
def get_class_weights(data):
    """
    [FIX] Calculates 'balanced' class weights to handle severe imbalance.
    Formula: TotalSamples / (NumClasses * ClassCount)
    This function now automatically finds the 'gene' train_mask inside the data object.
    """
    try:
        # [FIX] Get the mask and labels from the 'gene' node type
        mask = data['gene'].train_mask
        train_labels = data['gene'].y[mask]
        
        class_counts = torch.bincount(train_labels)
        
        if len(class_counts) < 2:
            log.warning("Only one class found in labels. Using equal weights.")
            return torch.tensor([0.5, 0.5])
            
        num_classes = len(class_counts)
        total_samples = class_counts.sum()
        
        # [FIX] Use 'balanced' weighting, do NOT normalize by sum
        weights = total_samples / (num_classes * class_counts.float())
        
        log.info(f"Calculated 'balanced' class weights (for 0, 1): {weights}")
        return weights
    except Exception as e:
        log.error(f"Error in get_class_weights: {e}")
        return torch.tensor([0.5, 0.5]) # Fallback

@torch.no_grad()
def evaluate(model, loader, loss_fn, device, desc="Evaluating"):
    model.eval()
    all_probs_list = []
    all_labels_list = []
    total_loss = 0
    total_nodes = 0

    for batch in tqdm(loader, desc=desc, leave=False):
        batch = batch.to(device)
        # [HAN] Model output is already 'gene' seed node embeddings
        out = model(batch)
        # [HAN] Get labels for 'gene' seed nodes
        labels_seed = batch['gene'].y[:batch['gene'].batch_size]
        batch_size = batch['gene'].batch_size
        
        if batch_size == 0: continue
            
        loss = loss_fn(out, labels_seed)
        total_loss += loss.item() * batch_size
        total_nodes += batch_size
        probs = F.softmax(out, dim=1)
        all_probs_list.append(probs.cpu())
        all_labels_list.append(labels_seed.cpu())

    if total_nodes == 0:
        log.warning("Evaluation set is empty. Returning zero metrics.")
        metrics_empty = {'loss': 0.0, 'accuracy': 0.0, 'f1': 0.0, 'precision': 0.0, 'recall': 0.0, 'roc_auc': 0.0, 'aucpr': 0.0}
        return metrics_empty, np.array([]), np.array([])

    all_probs = torch.cat(all_probs_list, dim=0)
    all_labels = torch.cat(all_labels_list, dim=0)
    avg_loss = total_loss / total_nodes
    preds_np = all_probs.argmax(dim=1).numpy()
    labels_np = all_labels.numpy()
    probs_np = all_probs[:, 1].numpy()

    try: roc_auc = roc_auc_score(labels_np, probs_np)
    except ValueError: roc_auc = 0.5
    try: aucpr = average_precision_score(labels_np, probs_np)
    except ValueError: aucpr = 0.0

    metrics = {
        'loss': avg_loss,
        'accuracy': accuracy_score(labels_np, preds_np),
        'f1': f1_score(labels_np, preds_np, zero_division=0),
        'precision': precision_score(labels_np, preds_np, zero_division=0),
        'recall': recall_score(labels_np, preds_np, zero_division=0),
        'roc_auc': roc_auc,
        'aucpr': aucpr
    }
    return metrics, labels_np, probs_np

def train(model, loader, optimizer, loss_fn, device, desc="Training"):
    model.train()
    total_loss = 0
    total_nodes = 0
    
    for batch in tqdm(loader, desc=desc, leave=False):
        batch = batch.to(device)
        optimizer.zero_grad()
        
        # [HAN] Model output is already 'gene' seed node embeddings
        out = model(batch)
        # [HAN] Get labels for 'gene' seed nodes
        labels_seed = batch['gene'].y[:batch['gene'].batch_size]
        batch_size = batch['gene'].batch_size
        
        if batch_size == 0: continue

        loss = loss_fn(out, labels_seed)
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * batch_size
        total_nodes += batch_size

    return total_loss / total_nodes if total_nodes > 0 else 0.0

@torch.no_grad()
def predict_all(model, data, nodes_df, gene_global_ids, device):
    """
    [HAN] Run inference on ALL 'gene' nodes and map predictions back to nodes_df.
    """
    log.info("Starting inference to get predictions for all GENE nodes...")
    model.eval()
    
    # [HAN] Create a loader for all 'gene' nodes
    all_gene_nodes = torch.arange(data['gene'].num_nodes)
    all_nodes_loader = NeighborLoader(
        data, num_neighbors=NUM_NEIGHBORS, batch_size=BATCH_SIZE,
        input_nodes=('gene', all_gene_nodes), # Predict on all genes
        shuffle=False, num_workers=NUM_WORKERS
    )
    
    all_probs_list = []
    for batch in tqdm(all_nodes_loader, desc="Predicting", leave=False):
        batch = batch.to(device)
        # [HAN] Model output is for 'gene' seed nodes
        out = model(batch)
        all_probs_list.append(F.softmax(out, dim=1).cpu())
        
    all_probs = torch.cat(all_probs_list, dim=0).numpy()
    final_preds = all_probs.argmax(axis=1)
    
    log.info("Inference complete. Mapping predictions back to nodes_df.")
    
    # Create a DataFrame with predictions, indexed by the original global node_id
    gene_preds_df = pl.DataFrame({
        'node_id': gene_global_ids,
        'predicted_class': final_preds,
        'prob_upregulated': all_probs[:, 1],
        'true_label': data['gene'].y.cpu().numpy(),
        'is_train_gene': data['gene'].train_mask.cpu().numpy(),
        'is_val_gene': data['gene'].val_mask.cpu().numpy(),
        'is_test_gene': data['gene'].test_mask.cpu().numpy()
    })
    
    # Join with the main nodes_df
    # Drop old/conflicting columns if they exist, then join
    cols_to_drop = ['predicted_class', 'prob_upregulated', 'true_label', 
                    'is_train_gene', 'is_val_gene', 'is_test_gene']
    existing_cols_to_drop = [c for c in cols_to_drop if c in nodes_df.columns]
    
    if existing_cols_to_drop:
        nodes_df = nodes_df.drop(existing_cols_to_drop)

    nodes_df = nodes_df.join(gene_preds_df, on='node_id', how='left')
    
    nodes_df.write_csv(PREDICTION_FILE)
    log.info(f"All node predictions saved to {PREDICTION_FILE}")

# --- 6. Main Execution ---
def main():
    log.info("--- Starting GNN Pipeline (HAN, 3-Layer, 3-Hop Fixed Split) ---")
    
    # 1. Load Data
    # [HAN] Returns HeteroData object
    data, nodes_df, gene_global_ids = load_graph_data()
    if data is None: return
    
    data = data.to('cpu') # Keep full graph data on CPU
    persist = (NUM_WORKERS > 0)
    
    # 2. Initialize Model, Loss, Optimizer
    log.info("--- Initializing Model and Optimzer ---")
    
    # Get class weights for the *training set*
    # [HAN] Pass the HeteroData object
    # [FIX] This call is now correct, as get_class_weights no longer needs 'mask'
    class_weights = get_class_weights(data).to(device)
    loss_fn = CrossEntropyLoss(weight=class_weights)
    
    # [HAN] Get in_channels from 'gene' node features
    in_channels = data['gene'].x.shape[1]
    
    # [HAN] Initialize HANModel with metadata
    model = HANModel(
        in_channels, HIDDEN_DIM, 2, NUM_HEADS, data.metadata()
    ).to(device)
    optimizer = AdamW(model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
    
    log.info("Model architecture:")
    log.info(str(model))

    # 3. Create Loaders
    log.info("--- Creating Data Loaders ---")
    # [HAN] Specify input_nodes as a tuple ('node_type', mask)
    train_loader = NeighborLoader(
        data, num_neighbors=NUM_NEIGHBORS, batch_size=BATCH_SIZE,
        input_nodes=('gene', data['gene'].train_mask), shuffle=True,
        num_workers=NUM_WORKERS, persistent_workers=persist
    )
    val_loader = NeighborLoader(
        data, num_neighbors=NUM_NEIGHBORS, batch_size=BATCH_SIZE,
        input_nodes=('gene', data['gene'].val_mask), shuffle=False,
        num_workers=NUM_WORKERS, persistent_workers=persist
    )
    test_loader = NeighborLoader(
        data, num_neighbors=NUM_NEIGHBORS, batch_size=BATCH_SIZE,
        input_nodes=('gene', data['gene'].test_mask), shuffle=False,
        num_workers=NUM_WORKERS, persistent_workers=persist
    )
    
    # 4. Main Training Loop
    log.info("--- Starting Model Training ---")
    
    # [FIX] Track best_val_aucpr for early stopping, as loss is misleading
    best_val_aucpr = -1.0 
    epochs_no_improve = 0
    all_epoch_metrics = [] # For saving epoch-level data

    for epoch in range(1, EPOCHS + 1):
        train_loss = train(model, train_loader, optimizer, loss_fn, device, desc="Training")
        val_metrics, _, _ = evaluate(model, val_loader, loss_fn, device, desc="Validating")
        
        log.info(f"Epoch {epoch:03d} | Train Loss: {train_loss:.4f} | "
                 f"Val Loss: {val_metrics['loss']:.4f} | Val F1: {val_metrics['f1']:.4f} | "
                 f"Val AUCPR: {val_metrics['aucpr']:.4f}")
        
        epoch_data = {'epoch': epoch, 'train_loss': train_loss, **{f'val_{k}': v for k, v in val_metrics.items()}}
        all_epoch_metrics.append(epoch_data)

        # [FIX] Early stopping based on AUCPR (higher is better)
        if val_metrics['aucpr'] > best_val_aucpr:
            best_val_aucpr = val_metrics['aucpr']
            epochs_no_improve = 0
            # Save the best model
            torch.save(model.state_dict(), MODEL_FILE)
            log.info(f"New best model saved (Val AUCPR: {best_val_aucpr:.4f})")
        else:
            epochs_no_improve += 1

        if epochs_no_improve == PATIENCE:
            log.info(f"Early stopping triggered at epoch {epoch}.")
            break
            
    # Save epoch metrics
    try:
        epoch_metrics_df = pl.DataFrame(all_epoch_metrics)
        epoch_metrics_df.write_csv(EPOCH_METRICS_FILE)
        log.info(f"Epoch-level metrics saved to {EPOCH_METRICS_FILE}")
    except Exception as e:
        log.error(f"Failed to save epoch metrics: {e}")

    # 5. Final Evaluation on Held-Out Test Set (Chr18-Y)
    log.info("--- Training Complete. Loading best model for final evaluation. ---")
    
    try:
        model.load_state_dict(torch.load(MODEL_FILE))
        log.info(f"Successfully loaded best model from {MODEL_FILE}")
    except Exception as e:
        log.error(f"Could not load best model state: {e}. Using last model state.")
        
    model.to(device) # Ensure model is on device
    
    test_metrics, test_labels_np, test_probs_np = evaluate(model, test_loader, loss_fn, device, desc="Final Test")

    log.info("--- Final Test Set Performance (Held-out: Chr 18-Y) ---")
    for k, v in test_metrics.items():
        log.info(f"  Test {k}: {v:.4f}")
        
    try:
        pl.DataFrame([test_metrics]).write_csv(TEST_METRICS_FILE)
        log.info(f"Final test metrics saved to {TEST_METRICS_FILE}")
    except Exception as e:
        log.error(f"Failed to save final test metrics: {e}")

    # 6. Generate and Save Plots for the *final test set*
    log.info("Generating and saving final test set ROC and PR curves...")
    if len(test_labels_np) > 0 and len(test_probs_np) > 0:
        try:
            fpr, tpr, _ = roc_curve(test_labels_np, test_probs_np)
            plt.figure(figsize=(8, 6))
            plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'HAN Model (AUC = {test_metrics["roc_auc"]:.4f})')
            plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Baseline (Random)')
            plt.xlabel('False Positive Rate'); plt.ylabel('True Positive Rate')
            plt.title('Receiver Operating Characteristic (ROC) Curve - Final Test Set')
            plt.legend(loc="lower right"); plt.savefig(ROC_PLOT_FILE); plt.close()
            log.info(f"ROC curve saved to {ROC_PLOT_FILE}")
        except Exception as e:
            log.error(f"Could not generate ROC curve plot: {e}")

        try:
            precision, recall, _ = precision_recall_curve(test_labels_np, test_probs_np)
            baseline = (test_labels_np == 1).sum() / len(test_labels_np)
            plt.figure(figsize=(8, 6))
            plt.plot(recall, precision, color='b', lw=2, label=f'HAN Model (AUCPR = {test_metrics["aucpr"]:.4f})')
            plt.plot([0, 1], [baseline, baseline], color='red', lw=2, linestyle='--', label=f'Baseline (Prevalence = {baseline:.4f})')
            plt.xlabel('Recall'); plt.ylabel('Precision')
            plt.title('Precision-Recall Curve - Final Test Set')
            plt.legend(loc="upper right"); plt.ylim([0.0, 1.05]); plt.savefig(PR_PLOT_FILE); plt.close()
            log.info(f"PR curve saved to {PR_PLOT_FILE}")
        except Exception as e:
            log.error(f"Could not generate PR curve plot: {e}")
    else:
        log.warning("Test set was empty, skipping plot generation.")

    # 7. Run inference on all 'gene' nodes
    predict_all(model, data, nodes_df, gene_global_ids, device)
    
    log.info("--- GNN Pipeline Finished ---")

if __name__ == "__main__":
    # Ensure file paths are correct relative to where you run the script.
    # If your CSVs are in '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/'
    # relative to your CWD, change BASE_PATH back.
    # For this example, I'm assuming CSVs are in the same dir as the script.
    # BASE_PATH = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/'
    # NODE_FILE = os.path.join(BASE_PATH, 'Nodes.csv')
    # ... etc. ...
    main()