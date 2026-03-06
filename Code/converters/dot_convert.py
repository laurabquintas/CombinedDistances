"""Convert mutation trees to DOT (GraphViz) format.

Exports anytree nodes to .gv files used by the MP3 and PD distance metrics.
Files are cached on disk to avoid redundant conversions.
"""

import os

from anytree.exporter import DotExporter

current_directory = os.getcwd()

# Default output directory for DOT tree files
DOT_OUTPUT_DIR = os.path.join(current_directory, "AML", "aml converts", "MP3 trees")


def dot_converter(tree, output_dir: str = DOT_OUTPUT_DIR) -> str:
    """Convert an anytree mutation tree to DOT format and cache the result.

    Args:
        tree: An anytree root node. Must have a ``data`` attribute used as
              a unique identifier for the output filename.
        output_dir: Directory where .gv files are stored.

    Returns:
        Path to the generated .gv file.
    """
    path = os.path.join(output_dir, f"dot_tree{tree.data}.gv")

    if not os.path.exists(path):
        DotExporter(tree).to_dotfile(path)

    return path
