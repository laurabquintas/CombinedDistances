"""Convert mutation trees to Newick format.

Newick format represents trees as nested parentheses. Multi-label nodes
are encoded as ``{label1,label2}``. Results are cached in a pickle file
to avoid redundant conversions.
"""

import os
import pickle
from typing import List

current_directory = os.getcwd()

# Default path for the Newick tree cache
NEWICK_CACHE_PATH = os.path.join(current_directory, "AML", "aml converts", "newick_trees.pickle")


def tree_to_newick(nodes: List) -> str:
    """Recursively convert a list of tree nodes to a Newick string.

    Args:
        nodes: List of anytree nodes (typically the children of a parent).

    Returns:
        Newick-formatted string representation.
    """
    items = []
    for node in nodes:
        subtree_str = ''
        if len(node.children) > 0:
            sub_tree = tree_to_newick(node.children)
            if sub_tree != '':
                subtree_str += '(' + sub_tree + ')'

        # Multi-label nodes use {label1,label2} notation
        if isinstance(node.name, list):
            subtree_str += '{' + ','.join(str(label) for label in node.name) + '}'
        else:
            subtree_str += str(node.name)

        items.append(subtree_str)

    return ','.join(items)


def newick_converter(tree, cache_path: str = NEWICK_CACHE_PATH) -> str:
    """Convert an anytree mutation tree to Newick format with disk caching.

    Args:
        tree: An anytree root node. Must have a ``data`` attribute used as
              a unique cache key.
        cache_path: Path to the pickle file used for caching conversions.

    Returns:
        Newick-formatted string for the tree.
    """
    # Load existing cache if available
    newick_cache = {}
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            newick_cache = pickle.load(f)

    if tree.data in newick_cache:
        return newick_cache[tree.data]

    # Convert and update cache
    newick_cache[tree.data] = tree_to_newick([tree])

    with open(cache_path, "wb") as f:
        pickle.dump(newick_cache, f)

    return newick_cache[tree.data]
