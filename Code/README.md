# Code

## Structure

| Folder | Description |
|--------|-------------|
| `src/` | Core pipeline modules — distance computation, tree structures, survival analysis, clustering, and evaluation utilities |
| `notebooks/` | Jupyter notebooks for the full AML analysis and step-by-step survival pipeline |
| `converters/` | Tree format converters (DOT, MLTD, Newick) |
| `optimization/` | Differential evolution weight optimization and Group LASSO survival models (R) |
| `metrics/` | Distance metric implementations (BD, CASet/DISC, GraPhyC, MLTD, MP3) |

## Quick Start

1. Install dependencies: `pip install -r ../requirements.txt`
2. Open `notebooks/AML_final_results.ipynb` for the complete AML analysis
3. Open `notebooks/Individual_survival_pipeline_aml_gen.ipynb` for the step-by-step survival pipeline

## Pipeline Overview

```
Patient Mutation Trees
        │
        ▼
  src/distance_matrix.py  ──►  Results/distance_matrices/  (21 CSVs)
        │
        ▼
  notebooks/  ──►  Clustering, survival analysis, metric evaluation
        │
        ▼
  optimization/  ──►  Differential evolution (optimal metric weights)
```
