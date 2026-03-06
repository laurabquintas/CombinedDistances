"""Compute pairwise distance matrices between patient mutation trees.

For each distance metric, builds a symmetric matrix of distances between
all pairs of patients and saves the result as a CSV file.
"""

import os
import pickle
from typing import Dict, List

import pandas as pd

import distance_wrapper as dist

# All available distance metrics to compute
ALL_METHODS = [
    'DISC_union', 'DISC_intersection',
    'CASet_intersection', 'CASet_union',
    'BD', '1BD', '2BD',
    'AD', 'CD', 'PD', 'PCD',
    'MLTD',
    'MP3_union', 'MP3_sigmoid', 'MP3_intersection', 'MP3_geometric',
]

DECIMAL_PLACES = 3


def load_tree_map(pickle_path: str) -> Dict:
    """Load the patient-to-tree mapping from a pickle file.

    Args:
        pickle_path: Path to the pickle file containing the tree map.

    Returns:
        Dictionary mapping patient IDs to their mutation trees.
    """
    with open(pickle_path, "rb") as f:
        return pickle.load(f)


def compute_distance_matrix(tree_map: Dict, method: str) -> List[List[float]]:
    """Compute a symmetric pairwise distance matrix for a given metric.

    For patients with multiple trees, the distance is the average across
    all tree pair comparisons.

    Args:
        tree_map: Dictionary mapping patient IDs to dicts of mutation trees.
        method: Name of the distance metric (see ``ALL_METHODS``).

    Returns:
        Symmetric matrix of pairwise distances (list of lists).
    """
    num_patients = len(tree_map)
    matrix = [[None for _ in range(num_patients)] for _ in range(num_patients)]

    patient_items = list(tree_map.items())

    for i, (patient_id_1, trees_1) in enumerate(patient_items):
        print(f'{method}: patient {i + 1}/{num_patients}')
        for j, (patient_id_2, trees_2) in enumerate(patient_items):
            if matrix[i][j] is not None:
                continue

            if i == j:
                matrix[i][j] = 0.0
            else:
                distances = set()
                for tree1 in trees_1.values():
                    for tree2 in trees_2.values():
                        d = dist.distance(tree1, tree2, method)
                        distances.add(d)
                matrix[i][j] = round(sum(distances) / len(distances), DECIMAL_PLACES)

            matrix[j][i] = matrix[i][j]

    return matrix


def save_distance_matrix(matrix: List[List[float]], patient_ids: List[str], output_path: str) -> None:
    """Save a distance matrix as a CSV file.

    Args:
        matrix: Symmetric distance matrix.
        patient_ids: List of patient identifiers (used as column headers).
        output_path: Path to the output CSV file.
    """
    df = pd.DataFrame(matrix, columns=patient_ids)
    df.to_csv(output_path, index=False)


if __name__ == '__main__':
    current_directory = os.getcwd()
    tree_map = load_tree_map(os.path.join(current_directory, "AML", "tree_map.pickle"))

    for method in ALL_METHODS:
        matrix = compute_distance_matrix(tree_map, method)
        patient_ids = list(tree_map.keys())
        output_path = os.path.join(current_directory, "AML", "Matrix Output", f"AML_{method}.csv")
        save_distance_matrix(matrix, patient_ids, output_path)
