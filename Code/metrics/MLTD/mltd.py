"""Multi-Labeled Tree Distance (MLTD) metric.

Computes similarity between two mutation trees using the MLTD algorithm
via a compiled C++ binary. The binary reads two tree files in MLTD text
format and outputs the similarity score.
"""

import os
import subprocess

current_directory = os.getcwd()

# Path to the compiled MLTD binary
MLTD_BINARY = os.path.join(current_directory, "MLTD", "main")


def calculate_mltd(path_tree1: str, path_tree2: str,
                   binary_path: str = MLTD_BINARY) -> float:
    """Compute the MLTD similarity between two trees.

    Calls the MLTD C++ binary as a subprocess and parses the similarity
    score from its stdout output.

    Args:
        path_tree1: File path to the first tree in MLTD text format.
        path_tree2: File path to the second tree in MLTD text format.
        binary_path: Path to the compiled MLTD binary.

    Returns:
        Similarity score in [0, 1], where 1 means identical trees.
    """
    args = [binary_path, path_tree1, path_tree2]
    result = subprocess.run(args, capture_output=True, text=True)

    # The last line of output contains the similarity score as the last word
    output_lines = result.stdout.strip().split("\n")
    similarity = float(output_lines[-1].split()[-1])

    return similarity
