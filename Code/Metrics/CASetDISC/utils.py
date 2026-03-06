"""Utility functions for the CASet and DISC distance metrics.

Provides Newick format parsing, ancestor set computation, Jaccard distance,
and matrix display helpers.
"""

from typing import Dict, Generator, List, Set, Tuple


def jaccard(set_a: Set, set_b: Set) -> float:
    """Compute the Jaccard distance between two sets.

    Note: Jaccard({}, {}) = 0 by convention.

    Args:
        set_a: First set.
        set_b: Second set.

    Returns:
        Jaccard distance in [0, 1].
    """
    union_size = len(set_a | set_b)
    intersection_size = len(set_a & set_b)
    return (1 - intersection_size / union_size) if union_size > 0 else 0


def ancestor_sets(newick: str) -> Dict[str, Set[str]]:
    """Compute the ancestor set for every mutation in a Newick tree.

    Args:
        newick: Newick string representation of a tree.

    Returns:
        Dict mapping each mutation label to its set of ancestors
        (including itself and co-located labels).
    """
    ancestors = {}
    _ancestor_sets_recursive(newick.replace(' ', ''), ancestors, set())
    return ancestors


def _ancestor_sets_recursive(newick: str, ancestors: Dict[str, Set[str]],
                             current_ancestors: Set[str]) -> None:
    """Recursive helper for ancestor set computation.

    Parses the Newick string to extract labels and build ancestor sets.

    Args:
        newick: Newick substring representing a subtree.
        ancestors: Running dict being populated (modified in place).
        current_ancestors: Mutations ancestral to everything in this subtree.
    """
    last_paren_index = newick.rfind(')')
    label_string = newick[last_paren_index + 1:]

    # Multi-label nodes use {label1,label2} notation
    if label_string[0] == '{':
        labels = set(label_string[1:-1].split(','))
    else:
        labels = {label_string}

    for label in labels:
        ancestors[label] = labels | current_ancestors

    if last_paren_index != -1:
        for subtree_newick in outer_comma_split(newick[1:last_paren_index]):
            _ancestor_sets_recursive(subtree_newick, ancestors, labels | current_ancestors)


def outer_comma_split(newick: str) -> Generator[str, None, None]:
    """Split a Newick substring on commas, ignoring nested parentheses and brackets.

    Yields substrings between top-level commas. Parenthesized and
    bracketed sections are treated as atomic units.

    Args:
        newick: The Newick string to split.

    Yields:
        Substrings between top-level commas.
    """
    chunk_start = 0
    paren_depth = 0
    bracket_depth = 0

    for i, char in enumerate(newick):
        if char == '(':
            paren_depth += 1
        elif char == ')':
            paren_depth -= 1
        elif char == '{':
            bracket_depth += 1
        elif char == '}':
            bracket_depth -= 1
        elif char == ',' and paren_depth == 0 and bracket_depth == 0:
            yield newick[chunk_start:i]
            chunk_start = i + 1

    yield newick[chunk_start:]


def get_matrix_string(matrix: List[List[float]]) -> str:
    """Format a distance matrix as a tab-separated string.

    Args:
        matrix: 2D distance matrix.

    Returns:
        Human-readable string with row/column indices.
    """
    header = 'tree\t' + '\t'.join(map(str, range(len(matrix)))) + '\n\n'
    rows = []
    for index, row in enumerate(matrix):
        rows.append('{}\t'.format(index) + '\t'.join(str(round(entry, 4)) for entry in row))
    return header + '\n'.join(rows) + '\n'


def get_min_max_string(matrix: List[List[float]]) -> str:
    """Find the tree pairs with minimum and maximum distances.

    Args:
        matrix: 2D distance matrix.

    Returns:
        String describing the min- and max-distance tree pairs.
    """
    n = len(matrix)
    min_value = float('inf')
    min_indices: Set[Tuple[int, int]] = set()
    max_value = -float('inf')
    max_indices: Set[Tuple[int, int]] = set()

    for row in range(n):
        for col in range(row + 1, n):
            entry = matrix[row][col]
            if entry == min_value:
                min_indices.add((row, col))
            elif entry < min_value:
                min_indices = {(row, col)}
                min_value = entry

            if entry == max_value:
                max_indices.add((row, col))
            elif entry > max_value:
                max_indices = {(row, col)}
                max_value = entry

    return (f'Max distance: {round(max_value, 6)} between tree pairs {max_indices}\n'
            f'Min distance: {round(min_value, 6)} between tree pairs {min_indices}')


def get_newicks(input_file: str) -> List[str]:
    """Read Newick strings from a file (one tree per line).

    Lines after ``#`` are treated as comments. Blank lines are skipped.

    Args:
        input_file: Path to the input file.

    Returns:
        List of Newick strings.
    """
    trees = []
    with open(input_file) as f:
        for line in f.readlines():
            line = line.strip().strip(';').strip()
            if '#' in line:
                line = line[:line.index('#')]
            if line and not line.isspace():
                trees.append(line)
    return trees
