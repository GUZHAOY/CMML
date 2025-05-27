import os
import scvi
import scanpy as sc
import cell2location
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

# ---------- Path Setup ----------
base_dir = "/exports/eddie/scratch/s2469905/cmml/benchmark2"
plot_dir = os.path.join(base_dir, "cell2location")
os.makedirs(plot_dir, exist_ok=True)

# ---------- Read Reference Data & Perform QC ----------
adata_ref = sc.read_10x_h5(
  os.path.join(base_dir, "5705STDY8058285_filtered_feature_bc_matrix.h5")
)
adata_ref.var_names_make_unique()

meta = pd.read_csv(os.path.join(base_dir, "cell_annotation.csv"), index_col=0)
meta.index = meta.index.str.replace(r'^.*?_', '', regex=True)
meta = meta[~meta.index.duplicated(keep="first")]
adata_ref.obs = meta.reindex(adata_ref.obs_names)

adata_ref.var['mt'] = adata_ref.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata_ref, qc_vars=['mt'], inplace=True)

# ---------- QC Plots ----------
# Violin plot for QC metrics
sc.pl.violin(
  adata_ref,
  ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
  jitter=0.4,
  show=False
)
plt.savefig(os.path.join(plot_dir, "scRNA_QC.png"), dpi=150)
plt.close()

# Scatter plot for total counts vs percentage of mitochondrial genes
sc.pl.scatter(
  adata_ref,
  x='total_counts',
  y='pct_counts_mt',
  show=False
)
plt.savefig(os.path.join(plot_dir, "scRNA_counts_vs_mt.png"), dpi=150)
plt.close()

# ---------- Clean Metadata: Remove Cells with Missing Annotations ----------
# Remove cells with missing annotations in 'annotation_1_print'
adata_ref = adata_ref[~adata_ref.obs['annotation_1_print'].isnull()].copy()

# Convert annotations to categorical for SCVI
adata_ref.obs['annotation_1_print'] = adata_ref.obs['annotation_1_print'].astype('category')

# ---------- SCVI Model Training ----------
scvi.model.SCVI.setup_anndata(
  adata_ref,
  labels_key='annotation_1_print'
)
vae = scvi.model.SCVI(adata_ref)

# Train the model (max_epochs=200 recommended)
vae.train(max_epochs=200)

# ---------- Get Latent Representations and Save to obsm ----------
latent = vae.get_latent_representation()
print("latent shape:", latent.shape, "adata_ref n_obs:", adata_ref.n_obs)

adata_ref.obsm["X_scVI"] = latent

# ---------- UMAP Visualization ----------
sc.pp.neighbors(adata_ref, use_rep="X_scVI")
sc.tl.umap(adata_ref)

sc.pl.umap(
  adata_ref,
  color="annotation_1_print",
  legend_loc='on data',
  show=False
)
plt.savefig(os.path.join(plot_dir, "reference_umap.png"), dpi=150)
plt.close()

# Load spatial transcriptomics data
adata_st = sc.read_10x_h5(
  os.path.join(base_dir, "ST8059048_filtered_feature_bc_matrix.h5")
)
adata_st.var_names_make_unique()

# Load spatial coordinate information
pos = pd.read_csv(
  os.path.join(base_dir, "spatial/tissue_positions_list.csv"),
  header=None,
  names=["barcode", "in_tissue", "array_row", "array_col", "imagerow", "imagecol"]
)
pos = pos[pos.in_tissue == 1].set_index("barcode")

# Align coordinates with data
common = adata_st.obs_names.intersection(pos.index)
adata_st = adata_st[common, :].copy()
adata_st.obs[['array_row', 'array_col']] = pos.loc[common, ['array_row', 'array_col']]
adata_st.obsm["spatial"] = adata_st.obs[['array_row', 'array_col']].to_numpy()

# Spatial data QC
sc.pp.calculate_qc_metrics(adata_st, inplace=True)
sc.pl.spatial(
  adata_st,
  color="total_counts",
  spot_size=30,
  show=False
)
plt.savefig(os.path.join(plot_dir, "UMI_heat.png"), dpi=150)
plt.close()

# ---------- Load Visium Data & Perform QC ----------
adata_st = sc.read_10x_h5(
  os.path.join(base_dir, "ST8059048_filtered_feature_bc_matrix.h5")
)
adata_st.var_names_make_unique()

# Read coordinate information again
pos = pd.read_csv(
  os.path.join(base_dir, "spatial/tissue_positions_list.csv"),
  header=None,
  names=["barcode", "in_tissue", "array_row", "array_col", "imagerow", "imagecol"]
)
pos = pos[pos.in_tissue == 1].set_index("barcode")

# Align spatial data with coordinates
common = adata_st.obs_names.intersection(pos.index)
adata_st = adata_st[common, :].copy()
adata_st.obs[["array_row", "array_col"]] = pos.loc[common, ["array_row", "array_col"]]

# Calculate QC metrics for spatial data
sc.pp.calculate_qc_metrics(adata_st, inplace=True)
adata_st.obsm["spatial"] = adata_st.obs[["array_row", "array_col"]].to_numpy()

# ---------- Spatial Plot: UMI Abundance ----------
# Plot the UMI abundance for spatial data
sc.pl.spatial(
  adata_st,
  color="total_counts",
  spot_size=30,
  show=False
)
plt.savefig(os.path.join(plot_dir, "UMI_heat.png"), dpi=150)
plt.close()

# ---------- (Assuming spatial_model has already been trained) ----------
# res = spatial_model.export_posterior(...)

# ---------- Merge Results ----------
means = res.obsm["means_cell_abundance_w_sf"]
adata_st.obs = pd.concat([adata_st.obs, means], axis=1)

# ---------- ❸ Spatial Heatmap: Top 4 Most Abundant Cell Types ----------
top_ct = means.sum(axis=0).sort_values(ascending=False).head(4).index
for ct in top_ct:
  sc.pl.spatial(
    adata_st,
    color=ct,
    cmap="magma",
    spot_size=30,
    show=False
  )
  plt.savefig(os.path.join(plot_dir, f"{ct}.png"), dpi=150)
  plt.close()

# ---------- ❹ Spot × Cell Type Composition Stacked Bar Plot ----------
comp = means[top_ct].copy()

# Ensure spatial data and composition index align
spatial_subset = adata_st.obsm["spatial"][adata_st.obs_names.get_indexer(comp.index)]

# Calculate PCA components as pseudo-regions and add to composition data
comp["region"] = PCA(n_components=1).fit_transform(spatial_subset)[:, 0].round(1)

# Reshape to long format
comp_melt = comp.melt(id_vars="region", var_name="cell_type", value_name="prop")

# Plot stacked bar plot for cell-type composition by pseudo-region
plt.figure(figsize=(6, 4))
sns.barplot(
  data=comp_melt,
  x="region",
  y="prop",
  hue="cell_type",
  estimator=np.mean,
  errorbar=None
)
plt.ylabel("Mean proportion")
plt.title("Composition by pseudo-region")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "Stack_bar_region.png"), dpi=150)
plt.close()

# ---------- ❺–❻ Ground Truth Accuracy Evaluation & Error Distribution ----------
from sklearn.neighbors import NearestNeighbors

# Load ground truth coordinates
truth = pd.read_csv(
  os.path.join(base_dir, "144600.csv"),
  names=["x","y","size"],
  header=0
)

# Spot coordinates (already loaded in pos2)
spot_coords = pos2[["imagecol","imagerow"]].to_numpy()
cell_coords = truth[["x","y"]].to_numpy()

# Nearest-neighbor matching
nn = NearestNeighbors(n_neighbors=1).fit(spot_coords)
_, idx2 = nn.kneighbors(cell_coords)
assigned = pos2.index.to_numpy()[idx2.flatten()]

# Calculate true counts
true_counts = (
  pd.Series(assigned)
  .value_counts()
  .rename_axis("barcode")
  .reset_index(name="true_count")
)

# Calculate predicted counts
pred = means.sum(axis=1).rename("pred_count").reset_index()
pred.columns = ["barcode","pred_count"]

# Clean up barcode format
true_counts["barcode"] = true_counts["barcode"].str.replace("-1$","",regex=True)
pred["barcode"] = pred["barcode"].astype(str).str.replace("-1$","",regex=True)

# Merge true and predicted counts
merged = pd.merge(true_counts, pred, on="barcode", how="inner")
print(f"✅ Merged spot count: {len(merged)}")

# Compute Pearson correlation and RMSE
r, _ = pearsonr(merged["true_count"], merged["pred_count"])
mse = mean_squared_error(merged["true_count"], merged["pred_count"])
rmse = np.sqrt(mse)

# Plot scatter plot and fit line
plt.figure(figsize=(5,5))
plt.scatter(merged["true_count"], merged["pred_count"], s=10, alpha=0.5)
m, b = np.polyfit(merged["true_count"], merged["pred_count"], 1)
plt.plot(merged["true_count"], m*merged["true_count"] + b, 'r-', linewidth=1)
plt.xlabel("Ground truth (nuclei count)")
plt.ylabel("Predicted total cells")
plt.title(f"Cell2location vs Truth\nPCC = {r:.3f}, RMSE = {rmse:.1f}")
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "Eval_scatter.png"), dpi=150)
plt.close()

# Plot error distribution histogram
err = merged.assign(error = merged["true_count"] - merged["pred_count"])
plt.figure(figsize=(5,4))
sns.histplot(err.error, kde=True, bins=30)
plt.axvline(0, color='red', linestyle='--')
plt.title("Distribution of prediction error per spot")
plt.xlabel("True - Predicted")
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "Error_hist.png"), dpi=150)
plt.close()
