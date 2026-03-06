"""Bourque Distance (BD) metric for comparing mutation trees.

Computes the Bourque Distance and k-Bourque Distance between two trees.
BD measures structural similarity by comparing edge bipartitions:
for each edge in one tree, it checks whether an equivalent edge
(same descendant set and complement set) exists in the other tree.

When trees have different label sets, the intersection of labels is used
and the distance is computed via minimum-weight bipartite matching.
"""

import copy
from typing import Optional, Set, Tuple, List

import numpy as np
from anytree import PreOrderIter, ZigZagGroupIter
from scipy.optimize import linear_sum_assignment


def get_node_names(root, intersection: Optional[Set] = None,
                   include_node_list: bool = False):
    """Collect all node names from a tree, optionally filtered by an intersection set.

    Args:
        root: Root node of the tree.
        intersection: If provided, only include names present in this set.
        include_node_list: If True, also return a list of non-root nodes.

    Returns:
        If include_node_list is False: set of node names.
        If include_node_list is True: (set of names, list of non-root nodes).
    """
    names = set()
    non_root_nodes = []

    if intersection is None:
        for node in PreOrderIter(root):
            if not node.is_root:
                non_root_nodes.append(node)
            if isinstance(node.name, list):
                names.update(node.name)
            else:
                names.add(node.name)
    else:
        for node in PreOrderIter(root):
            if isinstance(node.name, list):
                for label in node.name:
                    if label in intersection:
                        names.add(label)
            else:
                if node.name in intersection:
                    names.add(node.name)

    if not include_node_list:
        return names
    return names, non_root_nodes


def BourqueDistance(root1, root2) -> float:
    """Compute the Bourque Distance between two mutation trees.

    If both trees share the same label set, uses Robinson-Foulds-style
    edge comparison. Otherwise, restricts to shared labels and uses
    bipartite matching to find the minimum distance.

    Args:
        root1: Root node of the first tree (will be deep-copied).
        root2: Root node of the second tree (will be deep-copied).

    Returns:
        Normalized distance in [0, 1].
    """
    root1 = copy.deepcopy(root1)
    root2 = copy.deepcopy(root2)

    names1, nodes1 = get_node_names(root1, include_node_list=True)
    names2, nodes2 = get_node_names(root2, include_node_list=True)
    total_edges = len(nodes1) + len(nodes2)

    if names1 == names2:
        # Same label set: count differing edges (Robinson-Foulds style)
        diff_count = _count_differing_edges(root1, root2, nodes1, nodes2)
        diff_count += _count_differing_edges(root2, root1, nodes2, nodes1)
        return diff_count / total_edges
    else:
        # Different label sets: use intersection-based matching
        shared_labels = names1.intersection(names2)
        matching1 = _count_matching_edges(root1, root2, nodes1, nodes2, shared_labels)
        matching2 = _count_matching_edges(root2, root1, nodes2, nodes1, shared_labels)
        return (total_edges - min(matching1, matching2)) / total_edges


def _count_differing_edges(root1, root2, nodes1: List, nodes2: List) -> int:
    """Count edges in tree1 that have no equivalent in tree2.

    An edge is equivalent if removing it produces the same descendant set
    and complement set in both trees.
    """
    diff_count = 0
    for node1 in nodes1:
        original_parent = node1.parent
        node1.parent = None
        found_match = False

        for node2 in nodes2:
            original_parent2 = node2.parent
            node2.parent = None

            descendants1 = get_node_names(node1)
            complement1 = get_node_names(root1)
            descendants2 = get_node_names(node2)
            complement2 = get_node_names(root2)

            if len(descendants1) > 0 and descendants1 == descendants2 and complement1 == complement2:
                node2.parent = original_parent2
                found_match = True
                break

            node2.parent = original_parent2

        if not found_match:
            diff_count += 1
        node1.parent = original_parent

    return diff_count


def _count_matching_edges(root1, root2, nodes1: List, nodes2: List,
                          shared_labels: Set) -> int:
    """Count edges in tree1 with an equivalent in tree2 (restricted to shared labels)."""
    match_count = 0
    for node1 in nodes1:
        original_parent = node1.parent
        node1.parent = None

        for node2 in nodes2:
            original_parent2 = node2.parent
            node2.parent = None

            descendants1 = get_node_names(node1, intersection=shared_labels)
            complement1 = get_node_names(root1, intersection=shared_labels)
            descendants2 = get_node_names(node2, intersection=shared_labels)
            complement2 = get_node_names(root2, intersection=shared_labels)

            if len(descendants1) > 0 and descendants1 == descendants2 and complement1 == complement2:
                node2.parent = original_parent2
                match_count += 1
                break

            node2.parent = original_parent2

        node1.parent = original_parent

    return match_count


def kBourqueDistance(root1, root2, k: int) -> float:
    """Compute the k-Bourque Distance using subtree-level comparisons.

    Extracts all subtrees of depth k from both trees, computes BD between
    all pairs, and uses the Hungarian algorithm to find the minimum-cost
    assignment.

    Args:
        root1: Root node of the first tree (will be deep-copied).
        root2: Root node of the second tree (will be deep-copied).
        k: Maximum subtree depth for comparison.

    Returns:
        Normalized distance in [0, 1].
    """
    root1 = copy.deepcopy(root1)
    root2 = copy.deepcopy(root2)

    internal_nodes1 = [node for node in PreOrderIter(root1) if not node.is_leaf]
    internal_nodes2 = [node for node in PreOrderIter(root2) if not node.is_leaf]

    # Ensure list1 is the longer list for the cost matrix
    if len(internal_nodes1) < len(internal_nodes2):
        internal_nodes1, internal_nodes2 = internal_nodes2, internal_nodes1

    cost_matrix = np.full((len(internal_nodes1), len(internal_nodes1)), None)

    for node1 in internal_nodes1:
        original_parent1 = node1.parent
        node1.parent = None

        # Truncate to depth k
        saved_children1 = []
        levels1 = [[n for n in level] for level in ZigZagGroupIter(node1, maxlevel=k + 1)]
        if len(levels1) == k + 1:
            for leaf in levels1[-1]:
                saved_children1.append(leaf.children)
                leaf.children = []

        for node2 in internal_nodes2:
            original_parent2 = node2.parent
            node2.parent = None

            saved_children2 = []
            levels2 = [[n for n in level] for level in ZigZagGroupIter(node2, maxlevel=k + 1)]
            if len(levels2) == k + 1:
                for leaf in levels2[-1]:
                    saved_children2.append(leaf.children)
                    leaf.children = []

            bd = BourqueDistance(node1, node2)
            idx1 = internal_nodes1.index(node1)
            idx2 = internal_nodes2.index(node2)
            cost_matrix[idx1][idx2] = bd

            if saved_children2:
                for i, leaf in enumerate(levels2[-1]):
                    leaf.children = saved_children2[i]
            node2.parent = original_parent2

        if saved_children1:
            for i, leaf in enumerate(levels1[-1]):
                leaf.children = saved_children1[i]
        node1.parent = original_parent1

    cost_matrix = np.where(cost_matrix is None, 1.0, cost_matrix)
    cost_matrix = np.array(cost_matrix, dtype=float)

    row_indices, col_indices = linear_sum_assignment(cost_matrix)
    return cost_matrix[row_indices, col_indices].sum() / len(internal_nodes1)
