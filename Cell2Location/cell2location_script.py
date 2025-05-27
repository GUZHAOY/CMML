import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scvi.model import SCVI
from cell2location.models import Cell2location

# Set working directory
os.chdir("/exports/eddie/scratch/s2469905/cmml/benchmark2")

# ──────────────────────────────────────────────────────────────────────────────
# 1. Reference Data Processing and Model Training
# ──────────────────────────────────────────────────────────────────────────────

# 1.1 Read 10x h5 data and annotations
adata_ref = sc.read_10x_h5("5705STDY8058285_filtered_feature_bc_matrix.h5")
adata_ref.var_names_make_unique()
meta_ref = pd.read_csv("cell_annotation.csv", index_col=0)
# Remove index prefix, keep only barcode, and remove duplicates
meta_ref.index = meta_ref.index.str.replace(r'^.*?_', '', regex=True)
meta_ref = meta_ref[~meta_ref.index.duplicated(keep='first')]
# Rename cell type column
meta_ref = meta_ref.rename(columns={'annotation_1': 'cluster'})
# Reindex by obs_names and merge
meta_sel = meta_ref.reindex(adata_ref.obs_names)
adata_ref.obs['cluster'] = meta_sel['cluster'].values

# 1.2 QC filtering and normalization
adata_ref.var['mt'] = adata_ref.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_ref, qc_vars=['mt'], inplace=True)
adata_ref = adata_ref[
  (adata_ref.obs.n_genes_by_counts > 200) &
    (adata_ref.obs.n_genes_by_counts < 5000) &
    (adata_ref.obs.pct_counts_mt < 5),
  :
].copy()
sc.pp.normalize_total(adata_ref, target_sum=1e4)
sc.pp.log1p(adata_ref)

# 1.3 Filter by cell type and create counts layer
min_cells = 25
valid_ct = adata_ref.obs['cluster'].value_counts().loc[lambda x: x >= min_cells].index
adata_ref = adata_ref[adata_ref.obs['cluster'].isin(valid_ct), :].copy()
raw_counts = adata_ref.X.copy().astype(int)
adata_ref.layers['counts'] = raw_counts

# 1.4 Train scVI model
SCVI.setup_anndata(adata_ref, labels_key="cluster", batch_key=None)
vae = SCVI(adata_ref)
vae.train(max_epochs=200)

# 1.5 Build cell state matrix
counts_mat = adata_ref.layers['counts']
counts_mat = counts_mat.toarray() if hasattr(counts_mat, "toarray") else counts_mat
df_counts = pd.DataFrame(
  counts_mat,
  index=adata_ref.obs_names,
  columns=adata_ref.var_names
)
cell_state_df = df_counts.groupby(adata_ref.obs['cluster']).mean().T

# 1.6 Train Cell2location reference model
Cell2location.setup_anndata(adata_ref, layer="counts", labels_key=None, batch_key=None)
cell2loc_ref = Cell2location(
  adata_ref,
  cell_state_df=cell_state_df,
  N_cells_per_location=1
)
cell2loc_ref.train(max_epochs=300)
cell2loc_ref.save("cell2loc_reference_card-V3")

# ──────────────────────────────────────────────────────────────────────────────
# 2. Spatial Data Processing and Spatial Model
# ──────────────────────────────────────────────────────────────────────────────

# 2.1 Read spatial data
adata_sp = sc.read_10x_h5("ST8059048_filtered_feature_bc_matrix.h5")
adata_sp.var_names_make_unique()

# 2.2 Read and filter coordinate information
pos = pd.read_csv(
  "spatial/tissue_positions_list.csv", header=None,
  names=["barcode", "in_tissue", "array_row", "array_col", "imagerow", "imagecol"]
)
pos = pos[pos.in_tissue == 1].set_index("barcode")

# 2.3 Ensure barcode format consistency (only append "-1" once)
if not pos.index.str.endswith("-1").any() and adata_sp.obs_names.str.endswith("-1").any():
  pos.index = pos.index + "-1"

# 2.4 Slice spatial AnnData
common = adata_sp.obs_names.intersection(pos.index)
adata_sp = adata_sp[common, :].copy()

# 2.5 Merge coordinates into obs
adata_sp.obs[['array_row', 'array_col']] = pos.loc[common, ['array_row', 'array_col']]

# 2.6 QC filtering for low UMI spots
sc.pp.calculate_qc_metrics(adata_sp, inplace=True)
adata_sp = adata_sp[adata_sp.obs.total_counts > 100, :].copy()

# 2.7 Create counts layer
adata_sp.layers['counts'] = adata_sp.X.copy().astype(int)

# 2.8 Gene alignment
shared_genes = adata_sp.var_names.intersection(cell_state_df.index)
adata_sp = adata_sp[:, shared_genes].copy()
cell_state_df = cell_state_df.loc[shared_genes, :].copy()

# 2.9 Train spatial model
Cell2location.setup_anndata(adata_sp, layer="counts", labels_key=None, batch_key=None)
spatial_model = Cell2location(
  adata_sp,
  cell_state_df=cell_state_df,
  N_cells_per_location=1,
  detection_alpha=200
)
spatial_model.train(max_epochs=300)

# ──────────────────────────────────────────────────────────────────────────────
# 3. Export Results and Visualize in Bulk
# ──────────────────────────────────────────────────────────────────────────────

# 3.1 Export posterior means
res = spatial_model.export_posterior(
  adata_sp,
  sample_kwargs={'num_samples': 1000, 'batch_size': 1024}
)

# 3.2 Merge into adata_sp.obs
adata_sp.obs = pd.concat(
  [adata_sp.obs, res.obsm['means_cell_abundance_w_sf']],
  axis=1
)
