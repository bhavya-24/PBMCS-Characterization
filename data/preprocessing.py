import scanpy as sc
import matplotlib.pyplot as plt
# Load h5ad file
adata = sc.read_h5ad("D:\immune\project\data\pbmc_preprocessed.h5ad")

# Basic overview
print(adata)
adata.raw = adata  # Store original data

# Filter genes and cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # mitochondria genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Visualize QC
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# Filter out low-quality cells
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# Scale the data
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')  # example gene

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)  # clustering

# Plot clusters
sc.pl.umap(adata, color=['leiden'])

sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# Look at marker genes manually
marker_genes = {
    'T cells': ['CD3D'],
    'B cells': ['MS4A1'],
    'NK cells': ['GNLY'],
    'Monocytes': ['CD14', 'LYZ'],
    'Dendritic cells': ['FCER1A'],
    'Platelets': ['PPBP']
}

sc.pl.dotplot(adata, marker_genes, groupby='leiden')

# Based on what dotplot tells you, update this accordingly
cluster_to_celltype = {
    '0': 'CD4+ T cells',
    '1': 'B cells',
    '2': 'CD8+ T cells',
    '3': 'Monocytes',
    '4': 'NK cells',
    '5': 'Platelets',
    '6': 'Dendritic cells',
    # continue as needed
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_celltype)

# Count number of cells per type
print(adata.obs['cell_type'].value_counts())

# Visualize
adata.obs['cell_type'].value_counts().plot(kind='bar', title='Cell Counts by Type')

plt.xlabel("Cell Type")
plt.ylabel("Number of Cells")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()