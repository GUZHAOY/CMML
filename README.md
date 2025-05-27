# Cell2Location vs RCTD Benchmark

This repository contains scripts for performing deconvolution analysis using **RCTD** and **Cell2Location** methods on spatial transcriptomics data. It includes the full analysis pipeline, from data preprocessing and model training to visualization and evaluation.

## Directory Structure

- `RCTD/`: Scripts for running RCTD and generating results.
- `Cell2Location/`: Scripts for running Cell2Location and generating results.
- `Visualizations/`: Scripts for visualizing the results from both RCTD and Cell2Location, including comparison to ground truth.

## How to Run

1. **Running RCTD**:
   - Install required R packages: `Seurat`, `spacexr`, etc.
   - Execute `rctd_script.R` after setting up your data paths.

2. **Running Cell2Location**:
   - Install required Python packages: `scanpy`, `scvi`, `cell2location`, etc.
   - Execute `cell2location_script.py` after setting up your data paths.

3. **Visualizing Results**:
   - Visualize RCTD results using `rctd_visualization.R`.
   - Visualize Cell2Location results using `cell2location_visualization.py`.

## Dependencies

- Python 3.x
- R (for RCTD analysis)
- Required Python libraries: `scanpy`, `scvi`, `cell2location`, etc.
- Required R libraries: `Seurat`, `spacexr`, etc.
