"""Visualization utilities for hierarchical clustering results.

Provides functions for generating histograms, cluster maps, and annotated
heatmaps with clinical data overlays using seaborn and matplotlib.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, cut_tree
from sklearn.metrics import silhouette_score
import scipy.spatial as sp


# ---------------------------------------------------------------------------
# Distance distribution
# ---------------------------------------------------------------------------

def histograms(distance_matrix: np.ndarray, output_path: str, method: str) -> None:
    """Plot a histogram of pairwise distances (upper triangle only).

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        output_path: File path to save the plot.
        method: Name of the distance metric (used in the plot title).
    """
    bin_edges = np.arange(0, 1.1, 0.05)

    upper_triangle = []
    for i in range(len(distance_matrix)):
        for j in range(i + 1, len(distance_matrix)):
            upper_triangle.append(distance_matrix[i, j])

    plt.figure(figsize=(10, 10))
    plt.hist(upper_triangle, bins=bin_edges)
    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of {method}')
    plt.savefig(output_path, dpi=200, transparent=True)
    plt.clf()


# ---------------------------------------------------------------------------
# Cluster maps
# ---------------------------------------------------------------------------

CBAR_POSITION = (1.1, 0.2, 0.03, 0.4)


def cluster_map(method: str, linkage_method: str, distance_matrix: np.ndarray,
                output_path: str) -> None:
    """Generate and save a seaborn clustermap heatmap.

    Args:
        method: Name of the distance metric (used in the plot title).
        linkage_method: Hierarchical clustering linkage method (e.g., 'ward').
        distance_matrix: Symmetric pairwise distance matrix.
        output_path: File path to save the plot.
    """
    condensed = sp.distance.squareform(distance_matrix)
    linkage_result = linkage(condensed, method=linkage_method)
    cluster = sns.clustermap(
        distance_matrix, cmap="mako_r", cbar_pos=CBAR_POSITION,
        row_linkage=linkage_result, col_linkage=linkage_result,
    )
    cluster.fig.suptitle(f'Cluster Map of {method} with {linkage_method} linkage')
    cluster.savefig(output_path, dpi=300, transparent=True)


def _build_clinical_color_df(clinical_data: pd.DataFrame,
                             clinical_cols: List[str],
                             palette_types: List[str]) -> Tuple[pd.DataFrame, Dict]:
    """Build a color-coded DataFrame for clinical annotation.

    Args:
        clinical_data: DataFrame with clinical columns.
        clinical_cols: Column names to annotate.
        palette_types: Seaborn palette names for each column.

    Returns:
        Tuple of (color DataFrame for row_colors, {col: {value: color}} lookup).
    """
    color_dict = {}
    lookup = {}
    for col_name, palette in zip(clinical_cols, palette_types):
        column = clinical_data[col_name]
        palette_colors = sns.color_palette(palette, n_colors=len(column.unique()))
        value_to_color = dict(zip(column.unique(), palette_colors))
        color_dict[col_name] = column.map(value_to_color)
        lookup[col_name] = value_to_color

    return pd.DataFrame(color_dict), lookup


def clustermap_with_clinical(distance_matrix: np.ndarray, method: str,
                             clinical_data: pd.DataFrame, linkage_method: str,
                             output_path: str,
                             clinical_cols: Optional[List[str]] = None,
                             palette_types: Optional[List[str]] = None) -> None:
    """Generate a clustermap annotated with clinical data color bars.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        method: Name of the distance metric (used in the plot title).
        clinical_data: DataFrame with clinical columns.
        linkage_method: Hierarchical clustering linkage method.
        output_path: File path to save the plot.
        clinical_cols: Column names to use for annotation. Defaults to breast cancer columns.
        palette_types: Seaborn palette names. Must match length of clinical_cols.
    """
    if clinical_cols is None:
        clinical_cols = ['Overall_Tumor_Grade', 'Vital_Status', 'Receptor_Status_Patient',
                         'Metastatic_Dz', 'Stage']
    if palette_types is None:
        palette_types = ['Spectral', 'Set2', 'hls', 'muted', 'Set2']

    color_df, _ = _build_clinical_color_df(clinical_data, clinical_cols, palette_types)
    linkage_result = linkage(sp.distance.squareform(distance_matrix), method=linkage_method)

    cluster = sns.clustermap(
        distance_matrix, cmap='mako_r', cbar_pos=CBAR_POSITION,
        row_colors=color_df, row_linkage=linkage_result, col_linkage=linkage_result,
    )
    cluster.fig.suptitle(f'Cluster Map with Clinical Data of {method} with {linkage_method} linkage')
    cluster.savefig(output_path, dpi=300, transparent=True)


def clustermap_with_clinical_legend(distance_matrix: np.ndarray, method: str,
                                    clinical_data: pd.DataFrame, linkage_method: str,
                                    output_path: str,
                                    clinical_cols: Optional[List[str]] = None,
                                    palette_types: Optional[List[str]] = None) -> None:
    """Generate a clustermap with clinical annotation and separate legends.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        method: Name of the distance metric (used in the plot title).
        clinical_data: DataFrame with clinical columns.
        linkage_method: Hierarchical clustering linkage method.
        output_path: File path to save the plot.
        clinical_cols: Column names for annotation. Defaults to AML columns.
        palette_types: Seaborn palette names.
    """
    if clinical_cols is None:
        clinical_cols = ['PriorMalig', 'Diagnosis', 'Gender', 'Tx_group', 'VitalStatus']
    if palette_types is None:
        palette_types = ['Spectral', 'Set2', 'hls', 'muted', 'Set2']

    color_df, lookup = _build_clinical_color_df(clinical_data, clinical_cols, palette_types)
    linkage_result = linkage(sp.distance.squareform(distance_matrix), method=linkage_method)

    cluster = sns.clustermap(
        distance_matrix, cmap='mako_r', cbar_pos=CBAR_POSITION,
        row_colors=color_df, row_linkage=linkage_result, col_linkage=linkage_result,
    )

    # Add a separate legend below the figure for each clinical column
    legend_y = -0.05
    for col_name in clinical_cols:
        value_to_color = lookup[col_name]
        patches = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                       markersize=10, label=label)
            for label, color in value_to_color.items()
        ]
        cluster.fig.legend(
            handles=patches, loc='lower center', bbox_to_anchor=(0.6, legend_y),
            frameon=False, title=col_name, ncol=len(value_to_color),
        )
        legend_y -= 0.05

    cluster.fig.suptitle(f'Cluster Map with Clinical Data of {method} with {linkage_method} linkage')
    cluster.savefig(output_path, dpi=300, transparent=True)


def clustermap_with_mutations(distance_matrix: np.ndarray, method: str,
                              clinical_data: pd.DataFrame, linkage_method: str,
                              output_path: str,
                              mutation_cols: Optional[List[str]] = None,
                              palette_types: Optional[List[str]] = None) -> None:
    """Generate a clustermap annotated with mutation status and legends.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        method: Name of the distance metric (used in the plot title).
        clinical_data: DataFrame with mutation columns.
        linkage_method: Hierarchical clustering linkage method.
        output_path: File path to save the plot.
        mutation_cols: Mutation column names. Defaults to common AML genes.
        palette_types: Seaborn palette names.
    """
    if mutation_cols is None:
        mutation_cols = ['TP53', 'FLT3', 'NPM1', 'NRAS', 'IDH2', 'NRAS', 'KRAS']
    if palette_types is None:
        palette_types = ['Spectral', 'Set2', 'hls', 'muted', 'Set2', 'Spectral', 'Set2']

    color_df, lookup = _build_clinical_color_df(clinical_data, mutation_cols, palette_types)
    linkage_result = linkage(sp.distance.squareform(distance_matrix), method=linkage_method)

    cluster = sns.clustermap(
        distance_matrix, cmap='mako_r', cbar_pos=CBAR_POSITION,
        row_colors=color_df, row_linkage=linkage_result, col_linkage=linkage_result,
    )

    legend_y = -0.05
    for col_name in mutation_cols:
        value_to_color = lookup[col_name]
        patches = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
                       markersize=10, label=label)
            for label, color in value_to_color.items()
        ]
        cluster.fig.legend(
            handles=patches, loc='lower center', bbox_to_anchor=(0.6, legend_y),
            frameon=False, title=col_name, ncol=len(value_to_color),
        )
        legend_y -= 0.05

    cluster.fig.suptitle(f'Cluster Map with Clinical Data of {method} with {linkage_method} linkage')
    cluster.savefig(output_path, dpi=300, transparent=True)


# ---------------------------------------------------------------------------
# Cluster label extraction
# ---------------------------------------------------------------------------

def get_cluster_labels(distance_matrix: np.ndarray, n_clusters: int) -> np.ndarray:
    """Assign cluster labels using Ward hierarchical clustering.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        n_clusters: Number of clusters.

    Returns:
        Array of integer cluster labels for each sample.
    """
    linkage_result = linkage(sp.distance.squareform(distance_matrix), method='ward')
    return cut_tree(linkage_result, n_clusters=n_clusters).flatten()
