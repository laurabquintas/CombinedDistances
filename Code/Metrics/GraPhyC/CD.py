"""Root-to-node path (CD) distance metric.

Measures the Jaccard distance between the sets of root-to-node paths
in two mutation trees. Each path is represented as a sorted tuple of
node names from the root to the node.

Note: Not fully prepared for multi-labeled nodes.
"""

from typing import Set, Tuple

from anytree import PreOrderIter


def CD(tree1, tree2) -> float:
    """Compute the path-based (CD) distance between two trees.

    Args:
        tree1: Root node of the first tree.
        tree2: Root node of the second tree.

    Returns:
        Jaccard distance in [0, 1] between the root-to-node path sets.
    """
    paths1 = _create_root_to_node_paths(tree1)
    paths2 = _create_root_to_node_paths(tree2)
    return 1 - (len(paths1 & paths2) / len(paths1 | paths2))


def _create_root_to_node_paths(tree) -> Set[Tuple]:
    """Extract all root-to-node paths from a tree.

    Each path is a sorted tuple of node names from the root to a given node.
    Multi-label ancestor nodes have their individual labels added to the path.

    Args:
        tree: Root node of the tree.

    Returns:
        Set of sorted tuples representing root-to-node paths.
    """
    paths = set()
    for node in PreOrderIter(tree):
        path = [node.name]
        for ancestor in node.ancestors:
            if isinstance(ancestor.name, list):
                path.extend(ancestor.name)
            else:
                path.append(ancestor.name)
        paths.add(_sorted_tuple(path))
    return paths


def _sorted_tuple(path: list) -> Tuple:
    """Convert a path to a canonical sorted tuple representation.

    Args:
        path: List of node names.

    Returns:
        Sorted tuple (or single-element tuple if path has one element).
    """
    if len(path) == 1:
        return tuple(path)
    return tuple(sorted(path))
