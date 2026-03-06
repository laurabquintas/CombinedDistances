"""CASet distance metric for comparing mutation trees.

CASet (Common Ancestor Set) measures the distance between two trees by
comparing the common ancestor sets of all mutation pairs. For each pair
of mutations, it computes the Jaccard distance between their common
ancestor sets in tree1 and tree2, then averages over all pairs.

Supports intersection mode (only shared mutations) and union mode (all mutations).
"""

import argparse
import pickle
from typing import Dict, List, Set

import numpy as np

from CASetDISC.utils import ancestor_sets, jaccard, get_newicks, get_matrix_string, get_min_max_string


def caset(mutations: List[str], t1_ancestors: Dict[str, Set],
          t2_ancestors: Dict[str, Set]) -> float:
    """Compute the CASet distance between two trees stored as ancestor sets.

    Args:
        mutations: List of mutations over which to compute pairwise distances.
        t1_ancestors: Dict mapping each mutation to its ancestor set in tree 1.
        t2_ancestors: Dict mapping each mutation to its ancestor set in tree 2.

    Returns:
        Average Jaccard distance across all mutation pairs. Returns 1.0 if
        there is only one mutation (undefined comparison).
    """
    num_mutations = len(mutations)

    if num_mutations == 1:
        return 1.0

    total = 0.0
    for i in range(num_mutations):
        for j in range(i + 1, num_mutations):
            t1_ca_set = t1_ancestors[mutations[i]] & t1_ancestors[mutations[j]]
            t2_ca_set = t2_ancestors[mutations[i]] & t2_ancestors[mutations[j]]
            total += jaccard(t1_ca_set, t2_ca_set)

    return total / (num_mutations * (num_mutations - 1) / 2)


def caset_intersection(t1: str, t2: str) -> float:
    """Compute CASet distance using only mutations present in both trees.

    Args:
        t1: Newick string representation of tree 1.
        t2: Newick string representation of tree 2.

    Returns:
        CASet intersection distance.
    """
    t1_ancestors = ancestor_sets(t1)
    t2_ancestors = ancestor_sets(t2)
    shared_mutations = list(set(t1_ancestors.keys()) & set(t2_ancestors.keys()))
    return caset(shared_mutations, t1_ancestors, t2_ancestors)


def caset_union(t1: str, t2: str) -> float:
    """Compute CASet distance using all mutations from both trees.

    Mutations absent from a tree are assigned an empty ancestor set.

    Args:
        t1: Newick string representation of tree 1.
        t2: Newick string representation of tree 2.

    Returns:
        CASet union distance.
    """
    t1_ancestors = ancestor_sets(t1)
    t2_ancestors = ancestor_sets(t2)
    all_mutations = list(set(t1_ancestors.keys()) | set(t2_ancestors.keys()))

    for mutation in all_mutations:
        if mutation not in t1_ancestors:
            t1_ancestors[mutation] = set()
        if mutation not in t2_ancestors:
            t2_ancestors[mutation] = set()

    return caset(all_mutations, t1_ancestors, t2_ancestors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find pairwise CASet distances of a set of trees.')
    parser.add_argument('inputFile')
    parser.add_argument('-o', '--outputFile')
    parser.add_argument('-u', '--union', action='store_true')
    parser.add_argument('-p', '--pickle')
    parser.add_argument('-t', '--treePrint', action='store_true')
    parser.add_argument('-m', '--minmax', action='store_true')
    args = parser.parse_args()

    trees = get_newicks(args.inputFile)

    if args.treePrint:
        for i, tree in enumerate(trees):
            print('Tree {}:\n{}\n'.format(i, tree))

    distance_measure = caset_union if args.union else caset_intersection
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
