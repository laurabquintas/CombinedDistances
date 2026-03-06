"""DISC distance metric for comparing mutation trees.

DISC (Distinctly Inherited Set Comparison) measures the distance between
two trees by comparing the distinct ancestor sets of all ordered mutation
pairs. For each ordered pair (i, j), it computes the Jaccard distance
between the set of ancestors of i that are NOT ancestors of j.

Supports intersection mode (only shared mutations) and union mode (all mutations).
"""

import argparse
import pickle
from typing import Dict, List, Set

import numpy as np

from CASetDISC.utils import ancestor_sets, jaccard, get_newicks, get_matrix_string, get_min_max_string


def disc(mutations: List[str], t1_ancestors: Dict[str, Set],
         t2_ancestors: Dict[str, Set]) -> float:
    """Compute the DISC distance between two trees stored as ancestor sets.

    Args:
        mutations: List of mutations over which to compute pairwise distances.
        t1_ancestors: Dict mapping each mutation to its ancestor set in tree 1.
        t2_ancestors: Dict mapping each mutation to its ancestor set in tree 2.

    Returns:
        Average Jaccard distance across all ordered mutation pairs.
        Returns 1.0 if there is only one mutation.
    """
    num_mutations = len(mutations)

    if num_mutations == 1:
        return 1.0

    total = 0.0
    for i in range(num_mutations):
        for j in range(num_mutations):
            if i != j:
                t1_distinct = t1_ancestors[mutations[i]] - t1_ancestors[mutations[j]]
                t2_distinct = t2_ancestors[mutations[i]] - t2_ancestors[mutations[j]]
                total += jaccard(t1_distinct, t2_distinct)

    return total / (num_mutations * (num_mutations - 1))


def disc_intersection(t1: str, t2: str) -> float:
    """Compute DISC distance using only mutations present in both trees.

    Args:
        t1: Newick string representation of tree 1.
        t2: Newick string representation of tree 2.

    Returns:
        DISC intersection distance.
    """
    t1_ancestors = ancestor_sets(t1)
    t2_ancestors = ancestor_sets(t2)
    shared_mutations = list(set(t1_ancestors.keys()) & set(t2_ancestors.keys()))
    return disc(shared_mutations, t1_ancestors, t2_ancestors)


def disc_union(t1: str, t2: str) -> float:
    """Compute DISC distance using all mutations from both trees.

    Mutations absent from a tree are assigned an empty ancestor set.

    Args:
        t1: Newick string representation of tree 1.
        t2: Newick string representation of tree 2.

    Returns:
        DISC union distance.
    """
    t1_ancestors = ancestor_sets(t1)
    t2_ancestors = ancestor_sets(t2)
    all_mutations = list(set(t1_ancestors.keys()) | set(t2_ancestors.keys()))

    for mutation in all_mutations:
        if mutation not in t1_ancestors:
            t1_ancestors[mutation] = set()
        if mutation not in t2_ancestors:
            t2_ancestors[mutation] = set()

    return disc(all_mutations, t1_ancestors, t2_ancestors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find pairwise DISC distances of a set of trees.')
    parser.add_argument("inputFile")
    parser.add_argument('-o', '--outputFile')
    parser.add_argument('-u', '--union', action='store_true')
    parser.add_argument('-p', '--pickle')
    parser.add_argument('-t', '--treePrint', action='store_true')
    parser.add_argument('-m', '--minmax', action='store_true')
    args = parser.parse_args()

    trees = get_newicks(args.inputFile)

    if args.treePrint:
        for tree in trees:
            print(tree)

    distance_measure = disc_union if args.union else disc_intersection
    distance_matrix = [[None for _ in range(len(trees))] for _ in range(len(trees))]
    for i, t1 in enumerate(trees):
        for j, t2 in enumerate(trees):
            distance_matrix[i][j] = distance_measure(t1, t2)
            distance_matrix[j][i] = distance_matrix[i][j]

    output = get_matrix_string(distance_matrix)

    if args.minmax:
        output += get_min_max_string(distance_matrix)

    if args.pickle:
        with open(args.pickle, 'wb') as out:
            pickle.dump(distance_matrix, out, pickle.HIGHEST_PROTOCOL)

    if args.outputFile:
        with open(args.outputFile, 'w') as out:
            out.write(output)
    else:
        print(output)
