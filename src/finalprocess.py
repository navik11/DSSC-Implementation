import scanpy as sc
import pandas as pd
import os
from glob import glob

name = '_osmFISH'
adata = sc.read('../data/Dataset2_osmFISH.h5ad')
sc.pp.neighbors(adata) 

csv_dir = '/home/sachida/Documents/courses/cs690/ass2/DSSC/src/out_processed_osmFISH/'

# Get all CSV files in the directory
csv_files = glob(os.path.join(csv_dir, '*.csv'))

# Loop through each CSV file in the directory
for csv_file in csv_files:
    print(f"Processing {csv_file}...")

    # Read the predicted clusters from the current CSV file
    pred_df = pd.read_csv(csv_file)
    cluster_labels = pred_df['Label'].values

    # Check if the number of labels matches the number of cells in adata
    if len(cluster_labels) == adata.n_obs:
        adata.obs['cluster_labels'] = cluster_labels.astype(str)
    else:
        raise ValueError(f"The number of labels in {csv_file} does not match the number of cells in adata")

    # If UMAP hasn't been calculated, compute it
    if 'X_umap' not in adata.obsm:
        sc.pp.neighbors(adata)  # Compute the neighborhood graph
        sc.tl.umap(adata)

    file_name = os.path.splitext(os.path.basename(csv_file))[0]
    plot_file_name = f'{name}_{file_name}.png'
    sc.pl.umap(adata, color='cluster_labels', title=f"UMAP Clustered Labels - {file_name}", save=plot_file_name)

    print(f"Saved UMAP plot as {plot_file_name}")


# pred_df = pd.read_csv('out_processed_MERFISH/pred.csv_1.csv')
# cluster_labels = pred_df['Label'].values

# if len(cluster_labels) == adata.n_obs:
#     adata.obs['cluster_labels'] = cluster_labels.astype(str)
# else:
#     raise ValueError("The number of labels in pred.csv does not match the number of cells in adata")

# if 'X_umap' not in adata.obsm:
#     sc.tl.umap(adata)
    
# sc.pl.umap(adata, color='cluster_labels', title="UMAP Clustered Labels", save='umap_clusters.png')