import polars as pl
import os
from collections import deque, defaultdict, Counter
from tqdm import tqdm # For a progress bar, as this will take a moment

# --- 1. File Paths ---
BASE_PATH = '../../public/Multiomics/Step24_GNN_UpDEG_Prediction/' 
NODE_FILE = os.path.join(BASE_PATH, 'Nodes.csv')
LABEL_FILE = os.path.join(BASE_PATH, 'Labels.csv')
EDGE_FILE = os.path.join(BASE_PATH, 'Edges.csv')

# --- NEW: Output files for the hop count reports ---
ALL_GENES_HOPS_FILE = os.path.join(BASE_PATH, 'all_genes_hop_report.csv')
UPREG_GENES_HOPS_FILE = os.path.join(BASE_PATH, 'upregulated_genes_hop_report.csv')

def build_graph(edges_df):
    """
    Builds an undirected adjacency list (dict of sets) from the edge DataFrame.
    This is the most efficient structure for BFS.
    """
    print("Building graph adjacency list...")
    graph = defaultdict(set)
    
    # Use .to_numpy() for a fast iteration over rows
    for source, target in tqdm(edges_df.select(['node_id_source', 'node_id_target']).to_numpy()):
        graph[source].add(target)
        graph[target].add(source)
        
    print(f"Graph built: {len(graph)} nodes with connections.")
    return graph

def find_shortest_path_bfs(graph, start_node, target_node_set):
    """
    Performs a Breadth-First Search (BFS) to find the shortest hop count
    from a start_node to *any* node in the target_node_set.
    
    Returns:
    - Hop count (int) if a path is found.
    - -1 if no path exists.
    """
    
    # Check if the start node is already a target (e.g., if a gene was also an enhancer)
    if start_node in target_node_set:
        return 0
        
    queue = deque([(start_node, 0)]) # (node_id, hop_count)
    visited = {start_node}
    
    while queue:
        current_node, hops = queue.popleft()
        
        # Check all neighbors of the current node
        for neighbor in graph[current_node]:
            if neighbor not in visited:
                visited.add(neighbor)
                
                # Check if this neighbor is our target
                if neighbor in target_node_set:
                    return hops + 1 # Found the shortest path!
                
                # Not a target, so add to queue to explore later
                queue.append((neighbor, hops + 1))
                
    # If the queue becomes empty, we've explored everything
    # and found no path to the target set.
    return -1

def check_gene_hop_distance():
    print("--- Starting Gene -> Enhancer Hop Distance Analysis ---")
    
    try:
        # --- 2. Load Data ---
        print(f"Loading data from {BASE_PATH}...")
        nodes_df = pl.read_csv(NODE_FILE)
        edges_df = pl.read_csv(EDGE_FILE)
        labels_df = pl.read_csv(LABEL_FILE)
    except Exception as e:
        print(f"Error loading files: {e}")
        return

    # --- 3. Get Node ID Sets ---
    print("Identifying node types...")
    gene_node_ids = set(nodes_df.filter(
        pl.col('node_type') == 'gene'
    )['node_id'])
    
    enhancer_node_ids = set(nodes_df.filter(
        pl.col('node_type') == 'enhancer'
    )['node_id'])
    
    upregulated_gene_ids = set(labels_df.filter(
        pl.col('is_upregulated') == 1
    )['node_id'])

    # --- 4. Build Graph ---
    graph = build_graph(edges_df)

    # --- 5. Run BFS for all genes ---
    print(f"Running BFS for all {len(gene_node_ids)} gene nodes...")
    all_gene_results = []
    upregulated_gene_results = []
    
    # --- NEW: Store gene_id with hop count ---
    all_gene_hop_data = [] # List of {'gene_id': int, 'min_hops_to_enhancer': int}
    upreg_gene_hop_data = []

    for gene_id in tqdm(gene_node_ids, desc="Analyzing gene paths"):
        # Find shortest path from this gene to *any* enhancer
        hops = find_shortest_path_bfs(graph, gene_id, enhancer_node_ids)
        
        all_gene_results.append(hops) # For Counter
        all_gene_hop_data.append({
            'gene_id': gene_id,
            'min_hops_to_enhancer': hops
        })
        
        if gene_id in upregulated_gene_ids:
            upregulated_gene_results.append(hops) # For Counter
            upreg_gene_hop_data.append({
                'gene_id': gene_id,
                'min_hops_to_enhancer': hops
            })

    # --- 6. Report Statistics (to console) ---
    all_gene_counts = Counter(all_gene_results)
    upregulated_gene_counts = Counter(upregulated_gene_results)

    print("\n--- Hop Count Report for ALL Genes ---")
    print(f"Total Genes: {len(gene_node_ids)}")
    for hops, count in sorted(all_gene_counts.items()):
        if hops == -1:
            print(f"  - Isolated (No path to enhancer): {count} genes")
        else:
            print(f"  - Min {hops} hops to enhancer: {count} genes")
    
    print("\n--- Hop Count Report for UPREGULATED Genes ---")
    print(f"Total Upregulated Genes: {len(upregulated_gene_ids)}")
    for hops, count in sorted(upregulated_gene_counts.items()):
        if hops == -1:
            print(f"  - Isolated (No path to enhancer): {count} genes")
        else:
            print(f"  - Min {hops} hops to enhancer: {count} genes")
    print("--------------------------------------------------")

    # --- 7. NEW: Export results to CSV ---
    try:
        print(f"Exporting 'All Genes' hop report to {ALL_GENES_HOPS_FILE}...")
        all_genes_df = pl.DataFrame(all_gene_hop_data)
        all_genes_df = all_genes_df.join(
            nodes_df.select(['node_id', 'identifier', 'chr']),
            left_on='gene_id',
            right_on='node_id',
            how='left'
        )
        all_genes_df.write_csv(ALL_GENES_HOPS_FILE)
        
        print(f"Exporting 'Upregulated Genes' hop report to {UPREG_GENES_HOPS_FILE}...")
        upreg_genes_df = pl.DataFrame(upreg_gene_hop_data)
        upreg_genes_df = upreg_genes_df.join(
            nodes_df.select(['node_id', 'identifier', 'chr']),
            left_on='gene_id',
            right_on='node_id',
            how='left'
        )
        upreg_genes_df.write_csv(UPREG_GENES_HOPS_FILE)
        
        print("Export complete.")
    
    except Exception as e:
        print(f"Error exporting hop count CSV files: {e}")

if __name__ == "__main__":
    check_gene_hop_distance()