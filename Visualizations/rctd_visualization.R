# RCTD Visualization
base.dir   <- "C:/Users/17733/Desktop/cmml/ICA2/data/520-1"
out.dir    <- file.path(base.dir, "plots_rctd")
dir.create(out.dir, showWarnings = FALSE)

h5_ref     <- file.path(base.dir, "5705STDY8058285_filtered_feature_bc_matrix.h5")
anno_ref   <- file.path(base.dir, "cell_annotation.csv")
h5_st      <- file.path(base.dir, "ST8059048_filtered_feature_bc_matrix.h5")
spatial_tg <- file.path(base.dir, "ST8059048_spatial.tar.gz")

# ========== Load Required R Packages ==========
library(Seurat)
library(ggplot2)
library(patchwork)

# ========== 1.1 Read Data ==========
sc <- Read10X_h5(h5_ref) |> CreateSeuratObject(project = "scRef")
meta <- read.csv(anno_ref, row.names = 1, stringsAsFactors = FALSE)

# Remove duplicate barcodes
meta$barcode <- sub("^[^_]+_", "", rownames(meta))
meta_unique <- meta[!duplicated(meta$barcode), ]
rownames(meta_unique) <- meta_unique$barcode

# Keep only cells with annotations
valid_barcodes <- intersect(colnames(sc), meta_unique$barcode)
sc <- subset(sc, cells = valid_barcodes)
sc$cluster <- meta_unique[colnames(sc), "annotation_1"]

# ========== 1.2 Quality Control (QC) ==========
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-", assay = "RNA")

# QC plots
p.vln <- VlnPlot(sc, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p.sc  <- FeatureScatter(sc, "nCount_RNA", "percent.mt") + labs(title = "nCount vs %MT")
ggsave(file.path(out.dir, "QC_scRNA_vln.png"), p.vln, width = 9, height = 3)
ggsave(file.path(out.dir, "QC_scRNA_scatter.png"), p.sc,  width = 4, height = 4)

# ========== 1.3 Dimensionality Reduction & Visualization ==========
sc <- NormalizeData(sc) |> FindVariableFeatures() |> ScaleData() |> RunPCA()
sc <- RunUMAP(sc, dims = 1:20)

p.umap <- DimPlot(sc, group.by = "cluster", label = TRUE, repel = TRUE) +
  ggtitle("Reference snRNA-seq UMAP")
ggsave(file.path(out.dir, "UMAP_reference.png"), p.umap, width = 10, height = 8)

library(spacexr)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)

# ===== 1. Set File Paths =====
base.dir <- "C:/Users/17733/Desktop/cmml/ICA2/data/520-1"
out.dir  <- file.path(base.dir, "plots_rctd")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

# ===== 2. Load RCTD Object and Extract Data =====
myRCTD  <- readRDS(file.path(base.dir, "myRCTD.rds"))
weights <- as.matrix(assay(myRCTD, "weights"))
coords  <- spatialCoords(myRCTD)

# ===== 3. Reshape to Long Format Table =====
w.long <- weights |>
  as.data.frame() |>
  rownames_to_column("cell_type") |>
  pivot_longer(-cell_type, names_to = "barcode", values_to = "prop")

coords_df <- as.data.frame(coords) |>
  rownames_to_column("barcode")

# ===== 3.1 Spatial Heatmap (Top 4 Cell Types) =====
top4 <- names(sort(rowMeans(weights), decreasing = TRUE))[1:4]

for (ct in top4) {
  sub <- w.long[w.long$cell_type == ct, ]
  sub <- merge(sub, coords_df, by = "barcode")
  
  p <- ggplot(sub, aes(x = x, y = y, color = prop)) +
    geom_point(size = 0.8) +
    scale_color_viridis_c() +
    coord_equal() +
    ggtitle(paste("RCTD:", ct)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
  
  ggsave(file.path(out.dir, paste0("RCTD_spatial_", ct, ".png")),
         plot = p, width = 4, height = 4)
}

# ===== 3.2 Region Stacked Barplot (K-means Clustering) =====
clusters <- kmeans(coords, centers = 6)$cluster
cluster_df <- data.frame(barcode = rownames(coords), cluster = factor(clusters))

spot_comp <- w.long |>
  inner_join(cluster_df, by = "barcode") |>
  group_by(cluster, cell_type) |>
  summarise(mean_prop = mean(prop), .groups = "drop")

p.bar <- ggplot(spot_comp, aes(cluster, mean_prop, fill = cell_type)) +
  geom_col(position = "fill") +
  ylab("Proportion") +
  ggtitle("Cell-type composition by region (RCTD)") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(out.dir, "RCTD_bar_cluster.png"), p.bar, width = 5, height = 4)

# ===== 3.3 Heatmap (Top 10 Cell Types × Top 50 Spots) =====
sel.ct   <- head(order(rowMeans(weights), decreasing = TRUE), 10)
sel.spot <- head(order(colSums(weights), decreasing = TRUE), 50)

heatmat <- weights[sel.ct, sel.spot]

# Remove rows or columns with zero variance or containing NA
heatmat <- heatmat[rowSums(is.na(heatmat)) == 0, ]
heatmat <- heatmat[apply(heatmat, 1, sd) > 0, ]
heatmat <- heatmat[, apply(heatmat, 2, sd) > 0]

if (nrow(heatmat) >= 2 && ncol(heatmat) >= 2) {
  pheatmap(
    heatmat,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    filename = file.path(out.dir, "RCTD_heatmap.png"),
    width = 6,
    height = 4
  )
} else {
  warning("⚠️ Visualization data for heatmap is too sparse, skipping plot")
}

library(Seurat)
library(data.table)
library(RANN)
library(dplyr)
library(ggplot2)
library(spacexr)

# ========= 1. Read Segmentation Info and Match to Spot =========
seg <- fread(seg_path, col.names = c("x","y","size"))
pos <- read.csv(spatial_pos, header = FALSE,
                col.names = c("barcode","in_tissue","array_row","array_col","imagerow","imagecol"))
pos <- pos[pos$in_tissue == 1, ]
rownames(pos) <- pos$barcode

# Nearest-neighbor match nuclei to closest spot
spot_mat <- as.matrix(pos[, c("imagecol", "imagerow")])
cell_mat <- as.matrix(seg[, .(x, y)])
nn <- nn2(data = spot_mat, query = cell_mat, k = 1)
seg$assigned_spot <- rownames(pos)[nn$nn.idx[,1]]

# True cell count table
true_counts <- as.data.frame(table(seg$assigned_spot))
colnames(true_counts) <- c("barcode", "true_cell_count")
true_counts$barcode <- sub("-1$", "", true_counts$barcode)

# ========= 2. Load Spatial Expression Data =========
counts_sp <- Read10X_h5(h5_spatial)
nUMI_sp <- colSums(counts_sp)
names(nUMI_sp) <- sub("-1$", "", names(nUMI_sp))

# ========= 3. Read RCTD Weights =========
myRCTD <- readRDS(myRCTD_path)
weights <- assay(myRCTD, "weights")
colnames(weights) <- sub("-1$", "", colnames(weights))

# ========= 4. Build Cell Size =========
counts_ref <- Read10X_h5(h5_ref)
meta_ref   <- read.csv(anno_ref, row.names = 1, stringsAsFactors = FALSE)

meta_ref$barcode <- sub("^[^_]+_", "", rownames(meta_ref))
meta_unique <- meta_ref[!duplicated(meta_ref$barcode), ]
rownames(meta_unique) <- meta_unique$barcode

shared_barcodes <- intersect(colnames(counts_ref), rownames(meta_unique))
counts_ref <- counts_ref[, shared_barcodes]
nUMI_ref   <- colSums(counts_ref)
cell_types <- meta_unique[shared_barcodes, "annotation_1"]

cell_size <- tapply(nUMI_ref, cell_types, mean)
cell_size <- cell_size[rownames(weights)]  # Match order

# ========= 5. Convert Proportions to Predicted Cell Numbers =========
est_cells <- sweep(weights, 1, cell_size, FUN = "/")
est_cells <- sweep(est_cells, 2, nUMI_sp[colnames(est_cells)], FUN = "*")
total_est_cells <- colSums(est_cells)

# ========= 6. Evaluate Predicted vs True Nuclei Counts =========
merged <- inner_join(
  true_counts,
  data.frame(barcode = names(total_est_cells), pred_raw = total_est_cells),
  by = "barcode"
)
scale_factor <- sum(merged$true_cell_count) / sum(merged$pred_raw)
merged$pred_scaled <- merged$pred_raw * scale_factor

# Pearson & RMSE
pcc  <- cor(merged$true_cell_count, merged$pred_scaled)
rmse <- sqrt(mean((merged$true_cell_count - merged$pred_scaled)^2))
cat("📊 Pearson correlation (normalized):", round(pcc, 3), "\n")
cat("📉 RMSE (normalized):", round(rmse, 2), "\n")

# ========= 7. Visualization: Prediction vs Ground Truth =========
plot_path <- file.path(base_dir, "RCTD_vs_truth_scaled.pdf")

p <- ggplot(merged, aes(x = true_cell_count, y = pred_scaled)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "RCTD Predicted vs. True Nuclei Counts",
    subtitle = paste0("Pearson = ", round(pcc, 3), ", RMSE = ", round(rmse, 2)),
    x = "True nuclei count",
    y = "Predicted cell count (scaled)"
  )

ggsave(plot_path, p, width = 6, height = 6, dpi = 150)
cat("✅ Plot saved to:", plot_path, "\n")

# Calculate error
merged$error <- merged$true_cell_count - merged$pred_scaled

# Plot error distribution histogram + density curve
library(ggplot2)

# Plot: Error histogram + density + white background
p_error <- ggplot(merged, aes(x = error)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Distribution of Prediction Error per Spot",
    x = "True nuclei count - Predicted cell count",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# Save image
ggsave(
  filename = file.path(out.dir, "RCTD_Error_Distribution.png"),
  plot     = p_error,
  width    = 6,
  height   = 4,
  dpi      = 150
)
