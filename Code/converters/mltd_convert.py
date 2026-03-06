"""Convert mutation trees to the MLTD text format.

The MLTD algorithm expects trees in a custom text format where:
- Lines with ``=`` map a node to its mutation labels.
- Lines with ``:`` map a node to its children (adjacency list).
- The root node is always labelled ``Root``.
- The first ``:`` line defines the root and its children.

Example::

    Root=Label1,Label2
    N1=Label3,Label4
    Root:N1

Files are cached on disk to avoid redundant conversions.
"""

import os
from typing import List

from anytree import PreOrderIter

current_directory = os.getcwd()

# Default output directory for MLTD tree files
MLTD_OUTPUT_DIR = os.path.join(current_directory, "AML", "aml converts", "MLTD trees")


def tree_to_mltd_string(tree) -> str:
    """Convert an anytree mutation tree to the MLTD text format.

    Args:
        tree: An anytree root node.

    Returns:
        Multi-line string in MLTD format.
    """
    node_lines: List[str] = []
    edge_lines: List[str] = []

    for node in PreOrderIter(tree):
        # Format mutation labels (handle single or multi-label nodes)
        if isinstance(node.name, list):
            node_labels = ','.join(map(str, node.name))
        else:
            node_labels = str(node.name)

        # Node-to-label mapping
        if node.is_root:
            node_lines.append(f'Root={node_labels}')
        else:
            node_lines.append(f'N{node.data}={node_labels}')

        # Node-to-children mapping (adjacency list)
        if not node.is_leaf:
            children_str = ','.join(f'N{child.data}' for child in node.children)
            if node.is_root:
                edge_lines.append(f'Root:{children_str}')
            else:
                edge_lines.append(f'N{node.data}:{children_str}')

    return '\n'.join(node_lines + edge_lines)


def mltd_converter(tree, output_dir: str = MLTD_OUTPUT_DIR) -> str:
    """Convert an anytree mutation tree to MLTD format and cache the result.

    Args:
        tree: An anytree root node. Must have a ``data`` attribute used as
              a unique identifier for the output filename.
        output_dir: Directory where text files are stored.

    Returns:
        Path to the generated text file.
    """
    path = os.path.join(output_dir, f"str_tree{tree.data}.txt")

    if not os.path.exists(path):
        with open(path, 'w') as f:
            f.write(tree_to_mltd_string(tree))

    return path
