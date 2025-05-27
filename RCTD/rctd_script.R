#——01. Environment Setup——
library(Matrix)
library(Seurat)        # For reading 10x h5 and performing basic QC
library(spacexr)       # RCTD package
library(SpatialExperiment)
library(SummarizedExperiment)
library(dplyr)

#——02. Directory and File Paths——
## Please update the paths according to your actual directory
base.dir      <- "C:/Users/17733/Desktop/cmml/ICA2/data/520-1"
h5_ref        <- file.path(base.dir, "5705STDY8058285_filtered_feature_bc_matrix.h5")
anno_ref      <- file.path(base.dir, "cell_annotation.csv")
h5_spatial    <- file.path(base.dir, "ST8059048_filtered_feature_bc_matrix.h5")
spatial_tar   <- file.path(base.dir, "ST8059048_spatial.tar.gz")

#——03. Read and Quality Control snRNA-seq Reference Data——
# 3.1 Load counts → Create Seurat object
sc <- Read10X_h5(h5_ref) %>% 
  CreateSeuratObject(project = "scRef")

# 3.2 Calculate mitochondrial percentage
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")

# 3.3 Visualize QC (optional)
VlnPlot(sc, c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3)

# 3.4 Filter cells: for example, keep cells with 200–5000 genes and mt% < 5%
sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

# 3.5 Extract filtered counts and metadata
counts_ref <- GetAssayData(sc, slot = "counts")
meta_ref <- read.csv(anno_ref, row.names = 1, stringsAsFactors = FALSE)
## Assuming meta_ref has a column called "cluster" representing cell types
cell_types <- meta_ref[colnames(counts_ref), "cluster"]
nUMI_ref   <- colSums(counts_ref)
cell_types <- factor(meta_ref[colnames(counts_ref), "cluster"])

#——04. Construct Reference Object for RCTD——
reference <- createReference(
  counts    = counts_ref,
  cell_types= cell_types,
  nUMI      = nUMI_ref
)

#——05. Read and Construct Spatial Object——
# 5.1 Extract spatial coordinates and images
untar(spatial_tar, exdir = file.path(base.dir, "spatial"))
# 5.2 Read the coordinates
pos <- read.csv(file.path(base.dir, "spatial", "tissue_positions_list.csv"),
                header = FALSE, stringsAsFactors = FALSE)
colnames(pos) <- c("barcode","in_tissue","array_row","array_col","imagerow","imagecol")
pos <- pos[pos$in_tissue==1, ]
coords <- pos[, c("array_row","array_col")]
rownames(coords) <- pos$barcode

# 5.3 Read spatial counts and align coordinates
counts_sp  <- Read10X_h5(h5_spatial)
# Only keep spots with tissue
counts_sp  <- counts_sp[, intersect(colnames(counts_sp), rownames(coords))]
coords     <- coords[ colnames(counts_sp), ]

nUMI_sp    <- colSums(counts_sp)

# 5.4 Construct SpatialRNA object
puck <- createSpatialRNA(
  coords = coords,
  counts = counts_sp,
  nUMI   = nUMI_sp
)

#——06. Wrap Data with SpatialExperiment/SummarizedExperiment——
# 06. Wrap the data using SpatialExperiment/SummarizedExperiment
# Calculate cell counts for each cell type
cell_type_counts <- table(reference@cell_types)

# Filter out cell types with fewer than 25 cells
valid_cell_types <- names(cell_type_counts[cell_type_counts >= 25])

# Create a filtered cell_types vector
filtered_cell_types <- reference@cell_types[reference@cell_types %in% valid_cell_types]

# Filter reference data: select based on valid_cell_types
# Note: reference@counts is also a matrix (Matrix class), so we need to keep row-column consistency
filtered_counts <- reference@counts[, reference@cell_types %in% valid_cell_types]

# Create filtered reference data object
reference_se <- SummarizedExperiment(
  assays = list(counts = filtered_counts),
  colData = DataFrame(
    cell_type = filtered_cell_types,
    nUMI = reference@nUMI[reference@cell_types %in% valid_cell_types]
  )
)

# Create SpatialRNA object
spatial_spe <- SpatialExperiment(
  assays = list(counts = as.matrix(counts_sp)),
  colData = DataFrame(nUMI = nUMI_sp),
  spatialCoords = as.matrix(coords)
)

# Create RCTD object
myRCTD <- createRctd(spatial_spe, reference_se)

# Run RCTD (set doublet_mode here)
myRCTD <- runRctd(
  myRCTD,
  rctd_mode = "doublet"  # Choose different modes as needed (e.g., "full", "multi", or "doublet")
)

#——08. Inspect and Visualize Results (Example)——
weights <- assay(myRCTD, "weights")  # Cell type weights for each spot
print(head(weights))

#——09. Save RCTD and Weights to RDS——
saveRDS(myRCTD, file = "myRCTD.rds")  # Save the RCTD object
saveRDS(weights, file = "weights.rds")  # Save the weights
