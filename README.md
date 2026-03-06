## Predictive Power of Distance Metrics for Cancer Patient Survival

> **Presented at ETH-DBSSE** (Computational Biology Group)
> and **EMBL Heidelberg** (OLISSIPO Annual Meeting)

### Overview

This project explored the predictive potential of distance metrics applied to mutation trees in cancer genomics. By analyzing the evolutionary pathways of tumors, we investigated whether these metrics—individually or in combination—can improve prognostic predictions beyond traditional genotype and clinical data.

### Motivation

Cancer is characterized by uncontrolled cell growth and genetic mutations, leading to intra-tumor heterogeneity that complicates treatment responses. Conventional therapies often target dominant cancer subclones, but post-remission phases may give rise to resistant subclones, challenging treatment efficacy. Mutation trees, representing the evolutionary relationships among tumor subclones, provide critical insights into drug resistance and disease progression. The application of distance metrics to analyze mutation trees has gained interest for its potential to reveal hidden patterns in tumor evolution. However, no standard metric adequately captures the complexity of all cancer types, motivating our research into optimizing metric combinations for improved survival prediction.

### Objectives

This study assesses the ability of distance metrics to predict patient survival in Acute Myeloid Leukemia and Breast Cancer. The primary goal was to determine whether these metrics—alone or combined—enhance prognostic predictions and validate their association with clinical outcomes. Using real-world datasets, we optimized metric combinations through survival data analysis, comparing stratification models based on tree-derived information with those relying solely on genotype and clinical data. The findings aimed to contribute to advancing individualized cancer treatment strategies.

### Distance Metrics

The pipeline computes pairwise distance matrices between patient mutation trees using 16 metric variants across 5 families:

| Family | Metrics | Source |
|--------|---------|--------|
| **Bourque Distance (BD)** | BD, 1BD, 2BD | Bourque et al. |
| **MP3** | MP3 (union, intersection, geometric, sigmoid) | Ciccolella et al. |
| **CASet / DISC** | CASet (intersection, union), DISC (intersection, union) | DiNardo et al. |
| **MLTD** | MLTD | Karpov et al. |
| **GraPhyC** | PD, PCD, AD, CD | Govek et al. |

Metric combinations are optimized using **differential evolution** and weighted distance matrices are evaluated via hierarchical clustering (Ward linkage), silhouette analysis, and survival modeling (Kaplan-Meier, Cox PH, log-rank tests, group LASSO).

### Datasets

- **Acute Myeloid Leukemia (AML)** — mutation trees reconstructed from patient tumor samples
- **Breast Cancer** — mutation trees from breast cancer patient cohorts

### Project Structure

```
CombinedDistances/
├── README.md
├── requirements.txt
├── .gitignore
├── Code/
│   ├── src/                        # Core pipeline modules
│   │   ├── distance_wrapper.py     # Unified interface to all 16 distance metrics
│   │   ├── distance_matrix.py      # Pairwise distance matrix computation
│   │   ├── tree.py                 # Tree data structures
│   │   ├── read_json.py            # JSON tree reader
│   │   ├── utils_survival.py       # Survival analysis utilities
│   │   ├── utils_evaluation.py     # Evaluation and plotting functions
│   │   └── utils_clustering.py     # Clustering and visualization utilities
│   ├── notebooks/                  # Jupyter analysis notebooks
│   │   ├── AML_final_results.ipynb                   # AML results notebook
│   │   └── Individual_survival_pipeline_aml_gen.ipynb # AML survival pipeline
│   ├── converters/                 # Tree format converters
│   │   ├── dot_convert.py          # DOT format converter
│   │   ├── mltd_convert.py         # MLTD format converter
│   │   └── newick_convert.py       # Newick format converter
│   ├── optimization/               # Differential evolution & R scripts
│   │   ├── diff_ev_trial_gen.py    # Cluster visualization and export
│   │   ├── diff_final.sh           # SLURM job submission script
│   │   └── Survival_R.R            # Group LASSO survival models (R)
│   └── metrics/                    # Distance metric implementations
│       ├── BD/                     # Bourque Distance
│       ├── CASetDISC/              # CASet and DISC
│       ├── GraPhyC/                # PD, PCD, AD, CD
│       ├── MLTD/                   # Multi-Labeled Tree Distance
│       └── MP3/                    # MP3 Tree Similarity
├── Results/
│   ├── distance_matrices/          # Pairwise distance CSVs (21 metrics)
│   ├── clustermap/                 # Clustermap heatmaps
│   ├── clustermap_clinical/        # Clustermaps with clinical annotation
│   ├── histograms/                 # Distance distribution histograms
│   ├── silhouette_plot/            # Silhouette score plots
│   ├── wss_plot/                   # Within-cluster sum of squares plots
│   ├── diff_evolution/             # Differential evolution outputs
│   └── survival/                   # Survival analysis results
│       ├── t_0.01/                 # Threshold = 0.01
│       ├── t_0.02/                 # Threshold = 0.02 (main results)
│       ├── t_0.03/                 # Threshold = 0.03
│       └── t_0.04/                 # Threshold = 0.04
└── Documents/
    ├── MScThesisLauraQuintas.pdf
    └── Extended_Abstract_LauraQuintas.pdf
```

### How to Run

```bash
pip install -r requirements.txt
```

The R survival analysis (`Survival_R.R`) additionally requires: `survival`, `survminer`, `glmnet`, `grpreg`, `penalized`, `gglasso`.

The main analysis workflows are in the Jupyter notebooks under `Code/notebooks/`. Start with `AML_final_results.ipynb` for the full AML analysis, or `Individual_survival_pipeline_aml_gen.ipynb` for the step-by-step survival pipeline.

### Documents

- [**MSc Thesis**](Documents/MScThesisLauraQuintas.pdf) — full thesis document
- [**Extended Abstract**](Documents/Extended_Abstract_LauraQuintas.pdf) — thesis abstract
