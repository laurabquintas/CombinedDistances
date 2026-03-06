"""Clustering and survival analysis utilities.

Provides functions for:
- Silhouette-based cluster quality evaluation
- Optimal cluster number selection (LOESS smoothing, weighted silhouette)
- Kaplan-Meier and Cox PH survival model helpers
- Visualization of clustering metrics
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import scipy.spatial as sp
from scipy.cluster.hierarchy import linkage, cut_tree
from sklearn.metrics import silhouette_score, silhouette_samples, calinski_harabasz_score
from lifelines.statistics import multivariate_logrank_test, pairwise_logrank_test, logrank_test
from lifelines import KaplanMeierFitter, CoxPHFitter
import statsmodels.stats.multitest as multitest
from statsmodels.nonparametric.smoothers_lowess import lowess
from itertools import combinations
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Silhouette analysis
# ---------------------------------------------------------------------------

def silhouette_plot(method: str, linkage_method: str, distance_matrix: np.ndarray,
                    output_path: str, max_clusters: int) -> None:
    """Plot silhouette scores across a range of cluster counts.

    Args:
        method: Name of the distance metric (used in the plot title).
        linkage_method: Hierarchical clustering linkage method (e.g., 'ward').
        distance_matrix: Symmetric pairwise distance matrix.
        output_path: File path to save the plot.
        max_clusters: Maximum number of clusters to evaluate.
    """
    linkage_result = linkage(sp.distance.squareform(distance_matrix), method=linkage_method)
    scores = []

    for n_clusters in range(2, max_clusters + 1):
        labels = cut_tree(linkage_result, n_clusters=n_clusters).flatten()
        score = silhouette_score(distance_matrix, labels, metric='precomputed')
        scores.append(score)

    sns.scatterplot(x=range(2, max_clusters + 1), y=scores)
    plt.title(f'Silhouette Score of {method} with {linkage_method} linkage')
    plt.xlabel('Number of Clusters')
    plt.ylabel('Silhouette Score')
    plt.xticks(range(2, max_clusters + 1))
    plt.savefig(output_path, dpi=200, transparent=True)
    plt.clf()


def silhouette_plot_detailed(method: str, n_clusters: int, distance_matrix: np.ndarray,
                             output_path: str) -> None:
    """Plot per-cluster silhouette coefficients as a detailed bar chart.

    Args:
        method: Name of the distance metric (used in the plot title).
        n_clusters: Number of clusters to evaluate.
        distance_matrix: Symmetric pairwise distance matrix.
        output_path: File path to save the plot.
    """
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(7, 7)
    ax.set_xlim([-0.1, 1])
    ax.set_ylim([0, len(distance_matrix) + (n_clusters + 1) * 10])

    linkage_result = linkage(sp.distance.squareform(distance_matrix), method='ward')
    cluster_labels = cut_tree(linkage_result, n_clusters=n_clusters).flatten()
    avg_score = silhouette_score(distance_matrix, cluster_labels, metric='precomputed')
    sample_scores = silhouette_samples(distance_matrix, cluster_labels, metric='precomputed')

    y_lower = 9
    for cluster_idx in range(n_clusters):
        cluster_scores = sample_scores[cluster_labels == cluster_idx]
        cluster_scores.sort()
        cluster_size = cluster_scores.shape[0]
        y_upper = y_lower + cluster_size

        color = cm.nipy_spectral(float(cluster_idx) / n_clusters)
        ax.fill_betweenx(np.arange(y_lower, y_upper), 0, cluster_scores,
                         facecolor=color, edgecolor=color, alpha=0.7)
        ax.text(-0.05, y_lower + 0.5 * cluster_size, str(cluster_idx))
        y_lower = y_upper + 10

    ax.set_title(f"Silhouette plot of {method} for {n_clusters} clusters")
    ax.set_xlabel("Silhouette coefficient")
    ax.set_ylabel("Cluster label")
    ax.axvline(x=avg_score, color="red", linestyle="--")
    ax.set_yticks([])
    ax.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.clf()


def weighted_silhouette(n_clusters: int, distance_matrix: np.ndarray) -> float:
    """Compute weighted silhouette score using cluster-size-weighted medians.

    Each cluster's median silhouette is weighted by its proportion of the
    total number of samples.

    Args:
        n_clusters: Number of clusters.
        distance_matrix: Symmetric pairwise distance matrix.

    Returns:
        Weighted silhouette score.
    """
    linkage_result = linkage(sp.distance.squareform(distance_matrix), method='ward')
    cluster_labels = cut_tree(linkage_result, n_clusters=n_clusters).flatten()
    sample_scores = silhouette_samples(distance_matrix, cluster_labels, metric='precomputed')
    total_samples = len(distance_matrix)

    weighted_score = 0.0
    for cluster_idx in range(n_clusters):
        cluster_scores = sample_scores[cluster_labels == cluster_idx]
        cluster_median = np.median(cluster_scores)
        weighted_score += cluster_median * len(cluster_scores) / total_samples

    return weighted_score


# ---------------------------------------------------------------------------
# Survival analysis metadata
# ---------------------------------------------------------------------------

def get_survival_dict(methods: List[str]) -> Dict:
    """Create an empty results dictionary for survival analysis metrics.

    Args:
        methods: List of distance metric names.

    Returns:
        Dictionary mapping each method to a template of survival metrics.
    """
    return {
        method: {
            'Clusters': int,
            'Weighted Silhouette Score': float,
            'Calinski Harabasz Score': float,
            'PH Assumptions': str,
            'Log Likelihood Ratio': float,
        }
        for method in methods
    }


# ---------------------------------------------------------------------------
# Plotting utilities
# ---------------------------------------------------------------------------

def line_plot(values: List[float], max_x: int, output_path: str,
              title: str, x_label: str, y_label: str, color: str) -> None:
    """Draw and save a line plot.

    Args:
        values: Y-axis values.
        max_x: Upper bound of x-axis range (starting from 2).
        output_path: File path to save the plot.
        title: Plot title.
        x_label: X-axis label.
        y_label: Y-axis label.
        color: Line color.
    """
    sns.lineplot(x=range(2, max_x), y=values, marker='o', color=color)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xticks(range(2, max_x))
    plt.savefig(output_path, dpi=300)
    plt.clf()


def show_plot(values: List[float], min_x: int, max_x: int,
              title: str, x_label: str, y_label: str, color: str) -> None:
    """Display an interactive line plot (not saved to file).

    Args:
        values: Y-axis values.
        min_x: Lower bound of x-axis range.
        max_x: Upper bound of x-axis range.
        title: Plot title.
        x_label: X-axis label.
        y_label: Y-axis label.
        color: Line color.
    """
    sns.lineplot(x=range(min_x, max_x), y=values, marker='o', color=color)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xticks(range(min_x, max_x))
    plt.show()
    plt.clf()


# ---------------------------------------------------------------------------
# Cluster quality evolution
# ---------------------------------------------------------------------------

def prev_cluster_increase(n_clusters: int, distance_matrix: np.ndarray,
                          linkage_result: np.ndarray,
                          patients_silhouette: Optional[Dict] = None,
                          patients_cluster: Optional[Dict] = None):
    """Compute the relative silhouette change when adding one more cluster.

    For n_clusters=2, returns initial silhouette and cluster size dicts.
    For n_clusters>2, returns the relative change in silhouette proportion
    for the newly split cluster.

    Args:
        n_clusters: Number of clusters.
        distance_matrix: Symmetric pairwise distance matrix.
        linkage_result: Precomputed linkage result.
        patients_silhouette: Previous per-cluster weighted silhouette dict.
        patients_cluster: Previous per-cluster size dict.

    Returns:
        For n=2: (silhouette_dict, cluster_size_dict)
        For n>2: (relative_change, updated_silhouette_dict, updated_cluster_size_dict)
    """
    cluster_labels = cut_tree(linkage_result, n_clusters=n_clusters).flatten()
    sample_scores = silhouette_samples(distance_matrix, cluster_labels, metric='precomputed')
    total_samples = len(distance_matrix)

    if n_clusters == 2:
        sil_dict = {}
        size_dict = {}
        for cluster_idx in range(2):
            cluster_scores = sample_scores[cluster_labels == cluster_idx]
            representative_idx = cluster_labels.tolist().index(cluster_idx)
            cluster_median = np.median(cluster_scores)
            sil_dict[representative_idx] = len(cluster_scores) * cluster_median / total_samples
            size_dict[representative_idx] = len(cluster_scores)
        return sil_dict, size_dict

    new_sil_dict = {}
    new_size_dict = {}
    new_cluster_total = 0.0
    prev_cluster_sil = 0.0
    unchanged_total = 0.0
    prev_unchanged_total = 0.0

    for cluster_idx in range(n_clusters):
        cluster_scores = sample_scores[cluster_labels == cluster_idx]
        representative_idx = cluster_labels.tolist().index(cluster_idx)
        cluster_size = len(cluster_scores)
        new_size_dict[representative_idx] = cluster_size
        cluster_weighted_sil = len(cluster_scores) * np.median(cluster_scores) / total_samples
        new_sil_dict[representative_idx] = cluster_weighted_sil

        # Detect if this cluster is new or changed (different size)
        if representative_idx not in patients_cluster or cluster_size != patients_cluster[representative_idx]:
            new_cluster_total += cluster_weighted_sil
            if representative_idx in patients_silhouette:
                prev_cluster_sil = patients_silhouette[representative_idx]
                for key in patients_silhouette:
                    if key != representative_idx:
                        prev_unchanged_total += patients_silhouette[key]
        else:
            unchanged_total += cluster_weighted_sil

    new_ratio = new_cluster_total / (unchanged_total + new_cluster_total)
    prev_ratio = prev_cluster_sil / (prev_unchanged_total + prev_cluster_sil)

    if new_ratio > prev_ratio:
        change = abs(new_ratio - prev_ratio)
    else:
        change = -abs(prev_ratio - new_ratio)

    return change, new_sil_dict, new_size_dict


def min_cluster_increase(n_clusters: int, distance_matrix: np.ndarray,
                         patients_cluster: Dict) -> Tuple[float, Dict]:
    """Compute the minimum silhouette ratio when splitting into n clusters.

    Args:
        n_clusters: Target number of clusters.
        distance_matrix: Symmetric pairwise distance matrix.
        patients_cluster: Previous per-cluster size dict.

    Returns:
        Tuple of (minimum_ratio, updated_cluster_size_dict).
    """
    linkage_result = linkage(sp.distance.squareform(distance_matrix), method='ward')
    cluster_labels = cut_tree(linkage_result, n_clusters=n_clusters).flatten()
    sample_scores = silhouette_samples(distance_matrix, cluster_labels, metric='precomputed')
    total_samples = len(distance_matrix)

    if n_clusters == 2:
        scores_by_cluster = []
        new_size_dict = {}
        for cluster_idx in range(2):
            cluster_scores = sample_scores[cluster_labels == cluster_idx]
            representative_idx = cluster_labels.tolist().index(cluster_idx)
            cluster_median = np.median(cluster_scores)
            weighted = len(cluster_scores) * cluster_median / total_samples
            scores_by_cluster.append(weighted)
            new_size_dict[representative_idx] = len(cluster_scores)
        ratio = min(scores_by_cluster[0] / scores_by_cluster[1],
                    scores_by_cluster[1] / scores_by_cluster[0])
        return ratio, new_size_dict

    new_size_dict = {}
    new_cluster_scores = []
    unchanged_total = 0.0

    for cluster_idx in range(n_clusters):
        cluster_scores = sample_scores[cluster_labels == cluster_idx]
        representative_idx = cluster_labels.tolist().index(cluster_idx)
        cluster_size = len(cluster_scores)
        new_size_dict[representative_idx] = cluster_size
        cluster_weighted = len(cluster_scores) * np.median(cluster_scores) / total_samples

        if representative_idx not in patients_cluster or cluster_size != patients_cluster[representative_idx]:
            new_cluster_scores.append(cluster_weighted)
        else:
            unchanged_total += cluster_weighted

    assert len(new_cluster_scores) == 2, "Expected exactly 2 new clusters from the split"

    score_0 = new_cluster_scores[0]
    ratio_0 = new_cluster_scores[0] / (unchanged_total + new_cluster_scores[1])
    ratio_1 = new_cluster_scores[1] / (unchanged_total + score_0)
    ratio = min(ratio_0, ratio_1)

    return ratio, new_size_dict


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def cumulative_values(values: List[float]) -> List[float]:
    """Compute a cumulative sum of a list of values.

    Args:
        values: Input values.

    Returns:
        List where element i is the sum of values[0..i].
    """
    cumulative = [values[0]]
    for i in range(1, len(values)):
        cumulative.append(values[i] + cumulative[i - 1])
    return cumulative


def difference_dict(total_silhouettes: Dict[str, List[float]]) -> Dict[str, List[float]]:
    """Compute signed differences between consecutive silhouette values per method.

    Args:
        total_silhouettes: Dict mapping method names to lists of silhouette scores.

    Returns:
        Dict of signed first-order differences for each method.
    """
    diffs = {}
    for method_name, scores in total_silhouettes.items():
        method_diffs = []
        for j in range(len(scores) - 1):
            if scores[j] > scores[j + 1]:
                method_diffs.append(abs(scores[j] - scores[j + 1]))
            else:
                method_diffs.append(-abs(scores[j] - scores[j + 1]))
        diffs[method_name] = method_diffs
    return diffs


def difference(silhouettes: List[float]) -> List[float]:
    """Compute signed first-order differences of a silhouette series.

    Args:
        silhouettes: List of silhouette scores.

    Returns:
        List of signed differences (positive = improvement, negative = decline).
    """
    diffs = []
    for i in range(len(silhouettes) - 1):
        if silhouettes[i + 1] > silhouettes[i]:
            diffs.append(abs(silhouettes[i] - silhouettes[i + 1]))
        else:
            diffs.append(-abs(silhouettes[i] - silhouettes[i + 1]))
    return diffs


def check_tree_clusters(tree_map: Dict, distance_matrix: np.ndarray,
                        n_clusters: int, n_trees: int) -> None:
    """Print the tree contents of each cluster for inspection.

    Args:
        tree_map: Mapping of patient IDs to trees.
        distance_matrix: Symmetric pairwise distance matrix.
        n_clusters: Number of clusters.
        n_trees: Number of trees to print per cluster.
    """
    linkage_result = linkage(sp.distance.squareform(distance_matrix), method='ward')
    cluster_labels = cut_tree(linkage_result, n_clusters=n_clusters).flatten()

    cluster_members = {}
    for i in range(len(cluster_labels)):
        label = cluster_labels[i]
        if label not in cluster_members:
            cluster_members[label] = [i + 1]
        else:
            cluster_members[label].append(i + 1)

    for cluster_idx in range(n_clusters):
        print(f'Cluster {cluster_idx + 1}:\n')
        for i in range(n_trees):
            for value in tree_map[cluster_members[cluster_idx][i]].values():
                print(value)


def silhouette_except_worst(distance_matrix: np.ndarray, n_clusters: int) -> float:
    """Compute weighted silhouette excluding the worst-performing cluster.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        n_clusters: Number of clusters.

    Returns:
        Sum of weighted median silhouettes after removing the worst cluster.
    """
    linkage_result = linkage(sp.distance.squareform(distance_matrix), method='ward')
    cluster_labels = cut_tree(linkage_result, n_clusters=n_clusters).flatten()
    sample_scores = silhouette_samples(distance_matrix, cluster_labels, metric='precomputed')
    total_samples = len(distance_matrix)

    cluster_weighted_scores = []
    for cluster_idx in range(n_clusters):
        cluster_scores = sample_scores[cluster_labels == cluster_idx]
        weighted = np.median(cluster_scores) * len(cluster_scores) / total_samples
        cluster_weighted_scores.append(weighted)

    cluster_weighted_scores.sort()
    cluster_weighted_scores.pop(0)  # Remove the worst cluster
    return sum(cluster_weighted_scores)


def number_patients_new_cluster(distance_matrix: np.ndarray, n_clusters: int) -> int:
    """Count patients in the smaller of the two newly formed clusters.

    When going from (n_clusters - 1) to n_clusters, one cluster splits
    into two. This returns the size of the smaller new cluster.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        n_clusters: Target number of clusters.

    Returns:
        Number of patients in the smaller new cluster.
    """
    linkage_result = linkage(sp.distance.squareform(distance_matrix), method='ward')
    prev_labels = cut_tree(linkage_result, n_clusters=n_clusters - 1).flatten()
    curr_labels = cut_tree(linkage_result, n_clusters=n_clusters).flatten()

    prev_sizes = {}
    for cluster_idx in range(n_clusters - 1):
        representative = prev_labels.tolist().index(cluster_idx)
        prev_sizes[representative] = int(np.sum(prev_labels == cluster_idx))

    new_cluster_sizes = []
    for cluster_idx in range(n_clusters):
        representative = curr_labels.tolist().index(cluster_idx)
        cluster_size = int(np.sum(curr_labels == cluster_idx))
        if representative not in prev_sizes or cluster_size != prev_sizes[representative]:
            new_cluster_sizes.append(cluster_size)

    assert len(new_cluster_sizes) == 2, "Expected exactly 2 new clusters from the split"
    return min(new_cluster_sizes)


# ---------------------------------------------------------------------------
# Optimal cluster number selection (LOESS-based)
# ---------------------------------------------------------------------------

def loess_plot(distance_matrix: np.ndarray, method: str,
               min_clusters: int, max_clusters: int,
               output_dir: Optional[str] = None) -> None:
    """Plot LOESS-smoothed cluster quality across cluster counts.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        method: Name of the distance metric (used in the plot title).
        min_clusters: Minimum cluster count.
        max_clusters: Maximum cluster count.
        output_dir: Directory to save the plot. If None, displays interactively.
    """
    x = np.array(range(min_clusters, max_clusters + 1))
    y = get_prev_inc_values(distance_matrix, min_clusters, max_clusters)
    smoothed = lowess(y, x)

    plt.figure()
    plt.scatter(x, y, label='Data')
    plt.plot(smoothed[:, 0], smoothed[:, 1], 'r', label='LOESS', color='darkorange')
    plt.legend()
    plt.xlabel('Number of clusters')
    plt.ylabel('Next Cluster Relative Increase (%)')
    plt.xticks(range(min_clusters, max_clusters + 1))
    plt.title(f'LOESS Smoothing for {method}')
    if output_dir is None:
        plt.show()
    else:
        plt.savefig(f'{output_dir}/{method}_loess.png', dpi=300, bbox_inches='tight')


def loess_simple_plot(distance_matrix: np.ndarray, method: str,
                      min_clusters: int, max_clusters: int,
                      output_dir: Optional[str] = None) -> None:
    """Plot LOESS-smoothed WSS differences across cluster counts.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        method: Name of the distance metric (used in the plot title).
        min_clusters: Minimum cluster count.
        max_clusters: Maximum cluster count.
        output_dir: Directory to save the plot. If None, displays interactively.
    """
    x = np.array(range(min_clusters, max_clusters + 1))
    y = diff_weighted_silhouette(distance_matrix, max_clusters)
    smoothed = lowess(y, x)

    plt.figure()
    plt.scatter(x, y, label='Data')
    plt.plot(smoothed[:, 0], smoothed[:, 1], 'r', label='LOESS', color='darkorange')
    plt.legend()
    plt.xlabel('Number of clusters')
    plt.ylabel('WSS difference')
    plt.xticks(range(min_clusters, max_clusters + 1))
    plt.title(f'LOESS Smoothing for {method}')
    if output_dir is None:
        plt.show()
    else:
        plt.savefig(f'{output_dir}/{method}_loess.png', dpi=300, bbox_inches='tight')


def n_clusters_loess(values: List[float], min_clusters: int, max_clusters: int,
                     threshold: float = 0.05) -> Optional[float]:
    """Find the first cluster count where LOESS-smoothed values stabilize.

    Args:
        values: Quality metric values across cluster counts.
        min_clusters: Minimum cluster count.
        max_clusters: Maximum cluster count.
        threshold: Values below this are considered stable.

    Returns:
        Cluster count where stabilization occurs, or None if not found.
    """
    x = np.array(range(min_clusters, max_clusters + 1))
    smoothed = lowess(values, x)
    stable_indices = np.where(np.abs(smoothed[:, 1]) < threshold)[0]
    return smoothed[stable_indices[0], 0] if stable_indices.size > 0 else None


def n_clusters_max_silhouette(distance_matrix: np.ndarray, min_clusters: int,
                              max_clusters: int) -> int:
    """Find the cluster count that maximizes weighted silhouette score.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        min_clusters: Minimum cluster count.
        max_clusters: Maximum cluster count.

    Returns:
        Optimal number of clusters.
    """
    scores = []
    for n in range(min_clusters, max_clusters + 1):
        scores.append(weighted_silhouette(n, distance_matrix))
    return scores.index(max(scores)) + min_clusters


def get_n_clusters(distance_matrix: np.ndarray, min_clusters: int = 3,
                   max_clusters: int = 20, threshold: float = 0.05) -> float:
    """Select optimal cluster count using LOESS stabilization or max silhouette.

    Prefers the LOESS-based estimate if it exists and is smaller than
    the silhouette-based estimate.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        min_clusters: Minimum cluster count.
        max_clusters: Maximum cluster count.
        threshold: LOESS stabilization threshold.

    Returns:
        Recommended number of clusters.
    """
    values = get_prev_inc_values(distance_matrix, max_x=max_clusters)
    loess_n = n_clusters_loess(values, min_clusters, max_clusters, threshold)
    sil_n = n_clusters_max_silhouette(distance_matrix, min_clusters, max_clusters)
    return loess_n if loess_n is not None and loess_n < sil_n else sil_n


def get_prev_inc_values(distance_matrix: np.ndarray, min_x: int = 3,
                        max_x: int = 20) -> List[float]:
    """Compute cluster-increase values for a range of cluster counts.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        min_x: Minimum number of clusters.
        max_x: Maximum number of clusters.

    Returns:
        List of relative silhouette change values.
    """
    linkage_result = linkage(sp.distance.squareform(distance_matrix), method='ward')
    sil_dict, size_dict = prev_cluster_increase(2, distance_matrix, linkage_result)
    values = []
    for n in range(min_x, max_x + 1):
        change, sil_dict, size_dict = prev_cluster_increase(
            n, distance_matrix, linkage_result, sil_dict, size_dict
        )
        values.append(change)
    return values


def diff_weighted_silhouette(distance_matrix: np.ndarray, max_clusters: int) -> List[float]:
    """Compute first-order differences of weighted silhouette scores.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        max_clusters: Maximum number of clusters.

    Returns:
        List of signed differences.
    """
    scores = [weighted_silhouette(n, distance_matrix) for n in range(3, max_clusters + 2)]
    return difference(scores)


def get_n_clusters_simple(distance_matrix: np.ndarray, min_clusters: int = 3,
                          max_clusters: int = 20, threshold: float = 0.02) -> float:
    """Select optimal cluster count using simplified WSS-based LOESS method.

    Args:
        distance_matrix: Symmetric pairwise distance matrix.
        min_clusters: Minimum cluster count.
        max_clusters: Maximum cluster count.
        threshold: LOESS stabilization threshold.

    Returns:
        Recommended number of clusters.
    """
    values = diff_weighted_silhouette(distance_matrix, max_clusters)
    loess_n = n_clusters_loess(values, min_clusters, max_clusters, threshold)
    sil_n = n_clusters_max_silhouette(distance_matrix, min_clusters, max_clusters)
    return loess_n if loess_n is not None and loess_n < sil_n else sil_n
