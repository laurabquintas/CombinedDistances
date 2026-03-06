"""Ancestor-Descendant (AD) distance metric.

Measures the Jaccard distance between the sets of ancestor-descendant
pairs in two mutation trees. A pair (a, d) is included if node a is an
ancestor of node d (at any depth).
"""

from typing import Set, Tuple

from anytree import PreOrderIter


def AD(tree1, tree2) -> float:
    """Compute the Ancestor-Descendant distance between two trees.

    Args:
        tree1: Root node of the first tree.
        tree2: Root node of the second tree.

    Returns:
        Jaccard distance in [0, 1] between the ancestor-descendant pair sets.
    """
    pairs1 = _create_ancestor_descendant_pairs(tree1)
    pairs2 = _create_ancestor_descendant_pairs(tree2)
    return 1 - (len(pairs1 & pairs2) / len(pairs1 | pairs2))


def _create_ancestor_descendant_pairs(tree) -> Set[Tuple]:
    """Extract all (ancestor, descendant) name pairs from a tree.

    Only includes pairs where the ancestor is an internal node.

    Args:
        tree: Root node of the tree.

    Returns:
        Set of (ancestor_name, descendant_name) tuples.
    """
    pairs = set()
    for node in PreOrderIter(tree):
        if not node.is_leaf:
            for descendant in node.descendants:
                pairs.add((node.name, descendant.name))
    return pairs
