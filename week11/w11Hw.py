#!/usr/bin/env python

import scanpy as sc
import matplotlib.pyplot as plt

# Read the 10x dataset filtered down to just the highly-variable genes
adata = sc.read_h5ad("variable_data.h5")
adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy 

sc.pp.neighbors(adata, n_neighbors = 10, n_pcs = 40)

# 1.2
sc.tl.leiden(adata)

# 1.3
fig, axes = plt.subplots(ncols = 2, figsize = (12, 5))

# Plot UMAP
sc.tl.umap(adata, maxiter = 900)
sc.pl.umap(adata, color = 'leiden', ax = axes[0], title = 'UMAP Clusters', show = False)

# Plot t-SNE
sc.tl.tsne(adata)
sc.pl.tsne(adata, color = 'leiden', ax = axes[1], title = 't-SNE Clusters', show = False)

# Show the figure
plt.savefig('umap_tsne_plots.png', dpi = 300, bbox_inches = 'tight')
plt.show()

#2.1

wilcoxon_adata = sc.tl.rank_genes_groups(adata, groupby = 'leiden', method = 'wilcoxon', use_raw = True, copy = True)

logreg_adata = sc.tl.rank_genes_groups(adata, groupby = 'leiden', method = 'logreg', use_raw = True, copy = True, max_iter = 900)

#2.2

fig, ax = plt.subplots(figsize = (12, 8), ncols = 5, nrows = 5, sharex = True, sharey = False)

sc.pl.rank_genes_groups(wilcoxon_adata, n_genes = 25, ax = ax, title = 'Wilcoxon Rank-Sum', use_raw = True, show = False)

fig.suptitle('Top 25 Genes for Each Cluster (Wilcoxon Rank-Sum)', fontsize=16)
plt.subplots_adjust(top = 0.9)  
plt.savefig('wilcoxon_ranked_adata.png', dpi = 300, bbox_inches = 'tight')
plt.close()

fig1, ax1 = plt.subplots(figsize = (12, 8), ncols = 5, nrows = 5, sharex = True, sharey = False)

sc.pl.rank_genes_groups(logreg_adata, n_genes = 25, ax = ax, title = 'Logistic Regression', use_raw = True, show = False)

fig1.suptitle('Top 25 Genes for Each Cluster (Logistic Regression)', fontsize=16)
plt.subplots_adjust(top = 0.9)  

plt.savefig('logreg_rankged_data.png', dpi = 300)
plt.close()


# 3.1
leiden = adata.obs['leiden']
umap = adata.obsm['X_umap']
tsne = adata.obsm['X_tsne']
adata = sc.read_h5ad('filtered_data.h5')
adata.obs['leiden'] = leiden
adata.obsm['X_umap'] = umap
adata.obsm['X_tsne'] = tsne

adata.write('filtered_clustered_data.h5')

adata = sc.read_h5ad("filtered_clustered_data.h5")
adata.uns['log1p']['base'] = None # This is needed due to a bug in scanpy

sc.tl.rank_genes_groups(adata, 'leiden')

# 3.2
marker_genes = ['LDHB','S100A8','CD79A','CCL5','FCGR3A','GNLY','FCER1A','PF4']
adata.rename_categories('leiden', marker_genes)

sc.pl.umap(adata, color = 'leiden', legend_loc = 'on data', show = False)

plt.suptitle('UMAP Clusters', fontsize = 16)
plt.subplots_adjust(top = 0.9)

plt.savefig('top_gene_labelled_umap.png', dpi = 300)

# 3.3
top_marker_genes = ['Unknown','Monocytes','B Cells','CD8+ T Cells','Leukocytes','NK Cells','APCs','HSCs']
adata.rename_categories('leiden', top_marker_genes)

sc.pl.umap(adata, color = 'leiden', legend_loc = 'on data', show = False)

plt.suptitle('UMAP Clusters', fontsize = 16)
plt.subplots_adjust(top = 0.9)

plt.savefig('labelled_umap.png', dpi = 300)
plt.show()

