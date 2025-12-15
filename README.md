# Rare Cell Detection Benchmark

This repository provides the code used in the manuscript:

**Revealing the Best Strategies for Rare Cell Type Detection in Multi-Sample Single-Cell Datasets**

The purpose of this repository is to support **methodological transparency and reproducibility** of the benchmarking analyses.

---

## 1. Overview

Rare cell type detection is a critical yet challenging task in single-cell RNA sequencing (scRNA-seq) analysis, particularly in multi-sample settings.  
In this study, we systematically benchmark representative rare cell detection methods under different analytical strategies.

### Evaluated methods
- CellSIUS  
- GapClust  
- GiniClust3  
- scCAD  
- SCISSORS  

### Detection strategies
1. Individual-sample detection  
2. Pooled detection across samples  
3. Pooled detection after batch correction  

---

## 2. Repository Structure

```text
rare-cell-detection-benchmark/
├── data/               # Example data and data description
├── scripts/            # Scripts for preprocessing, detection, and evaluation
├── results/            # Output results (generated)
├── figures/            # Figures (generated)
├── utils/              # Helper functions
├── environment.yml     # Conda environment for reproducibility
└── README.md
