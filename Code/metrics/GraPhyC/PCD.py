"""Parent-Child Distance (PCD) metric.

Measures the Jaccard distance between the sets of direct parent-child
pairs in two mutation trees (unlike AD which includes all ancestor-descendant
pairs at any depth).
"""

from typing import Set, Tuple

from anytree import PreOrderIter


def PCD(tree1, tree2) -> float:
    """Compute the Parent-Child distance between two trees.

    Args:
        tree1: Root node of the first tree.
        tree2: Root node of the second tree.

    Returns:
        Jaccard distance in [0, 1] between the parent-child pair sets.
    """
    pairs1 = _create_parent_child_pairs(tree1)
    pairs2 = _create_parent_child_pairs(tree2)
    intersection = pairs1 & pairs2
    union = pairs1 | pairs2
    return 1 - (len(intersection) / len(union))


def _create_parent_child_pairs(tree) -> Set[Tuple]:
    """Extract all direct (parent, child) name pairs from a tree.

    Args:
        tree: Root node of the tree.

    Returns:
        Set of (parent_name, child_name) tuples.
    """
    pairs = set()
    for node in PreOrderIter(tree):
        if not node.is_leaf:
            for child in node.children:
                pairs.add((node.name, child.name))
    return pairs
