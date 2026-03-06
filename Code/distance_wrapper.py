"""Unified interface for computing distances between mutation trees.

Dispatches to the appropriate metric implementation (BD, MP3, CASet, DISC,
MLTD, GraPhyC) based on the requested method name.
"""

from typing import Set

from anytree import PreOrderIter

from dot_convert import dot_converter
from mltd_convert import mltd_converter
from newick_convert import newick_converter
import BD.BDpython as BD
import MP3.mp3treesim.mp3treesim as MP3
import CASetDISC.CASet as CASet
import CASetDISC.DISC as DISC
import MLTD.mltd as MLTD
from GraPhyC.PD import PD
from GraPhyC.PCD import PCD
from GraPhyC.AD import AD
from GraPhyC.CD import CD


SUPPORTED_METHODS = [
    'BD', '1BD', '2BD',
    'MP3_union', 'MP3_sigmoid', 'MP3_geometric', 'MP3_intersection',
    'CASet_intersection', 'CASet_union',
    'DISC_intersection', 'DISC_union',
    'MLTD',
    'PD', 'PCD', 'AD', 'CD',
]

# Minimum node count thresholds required by certain metrics
_MIN_NODES_MP3 = 3
_MIN_NODES_PD = 2


def _get_node_names(tree) -> Set:
    """Extract all node names from a tree, flattening multi-label nodes."""
    names = set()
    for node in PreOrderIter(tree):
        if isinstance(node.name, list):
            names.update(node.name)
        else:
            names.add(node.name)
    return names


def _has_fewer_than_n_nodes(tree, n: int) -> bool:
    """Check if a tree has fewer than ``n`` nodes."""
    return len(set(PreOrderIter(tree))) < n


def _intersection_smaller_than(tree1, tree2, n: int) -> bool:
    """Check if the intersection of node names has fewer than ``n`` elements."""
    names1 = _get_node_names(tree1)
    names2 = _get_node_names(tree2)
    return len(names1.intersection(names2)) < n


def _compute_mp3_distance(tree1, tree2, mode: str) -> float:
    """Compute MP3-based distance for a given aggregation mode.

    Args:
        tree1: First mutation tree (anytree root node).
        tree2: Second mutation tree (anytree root node).
        mode: MP3 similarity mode ('union', 'sigmoid', 'geometric', 'intersection').

    Returns:
        Distance value in [0, 1], where 0 means identical and 1 means maximally different.
    """
    dot_path1 = dot_converter(tree1)
    dot_path2 = dot_converter(tree2)
    similarity = MP3.similarity(
        MP3.read_dotfile(dot_path1),
        MP3.read_dotfile(dot_path2),
        mode=mode,
    )
    return 1 - similarity


def distance(tree1, tree2, method: str) -> float:
    """Compute the distance between two mutation trees using the specified metric.

    Args:
        tree1: First mutation tree (anytree root node).
        tree2: Second mutation tree (anytree root node).
        method: Name of the distance metric. Must be one of ``SUPPORTED_METHODS``.

    Returns:
        Distance value in [0, 1].

    Raises:
        ValueError: If ``method`` is not supported.
    """
    # --- MP3 variants (require >=3 nodes and >=3 shared nodes) ---
    if method == 'MP3_union':
        if _has_fewer_than_n_nodes(tree1, _MIN_NODES_MP3) and _has_fewer_than_n_nodes(tree2, _MIN_NODES_MP3):
            return 1.0
        return _compute_mp3_distance(tree1, tree2, mode='union')

    if method.startswith('MP3'):
        too_few = (_has_fewer_than_n_nodes(tree1, _MIN_NODES_MP3)
                   or _has_fewer_than_n_nodes(tree2, _MIN_NODES_MP3)
                   or _intersection_smaller_than(tree1, tree2, _MIN_NODES_MP3))
        if too_few:
            return 1.0
        mp3_mode = method.split('_', 1)[1]  # e.g. 'sigmoid', 'geometric', 'intersection'
        return _compute_mp3_distance(tree1, tree2, mode=mp3_mode)

    # --- Bourque Distance variants ---
    if method == 'BD':
        return BD.BourqueDistance(tree1, tree2)
    if method == '1BD':
        return BD.kBourqueDistance(tree1, tree2, k=1)
    if method == '2BD':
        return BD.kBourqueDistance(tree1, tree2, k=2)

    # --- CASet variants ---
    if method == 'CASet_intersection':
        newick1 = newick_converter(tree1)
        newick2 = newick_converter(tree2)
        return CASet.caset_intersection(newick1, newick2)
    if method == 'CASet_union':
        newick1 = newick_converter(tree1)
        newick2 = newick_converter(tree2)
        return CASet.caset_union(newick1, newick2)

    # --- DISC variants ---
    if method == 'DISC_intersection':
        newick1 = newick_converter(tree1)
        newick2 = newick_converter(tree2)
        return DISC.disc_intersection(newick1, newick2)
    if method == 'DISC_union':
        newick1 = newick_converter(tree1)
        newick2 = newick_converter(tree2)
        return DISC.disc_union(newick1, newick2)

    # --- MLTD ---
    if method == 'MLTD':
        mltd_path1 = mltd_converter(tree1)
        mltd_path2 = mltd_converter(tree2)
        return 1 - MLTD.calculate_mltd(mltd_path1, mltd_path2)

    # --- GraPhyC variants ---
    if method == 'PD':
        if _intersection_smaller_than(tree1, tree2, _MIN_NODES_PD):
            return 1.0
        dot_path1 = dot_converter(tree1)
        dot_path2 = dot_converter(tree2)
        return PD(dot_path1, dot_path2)
    if method == 'PCD':
        return PCD(tree1, tree2)
    if method == 'AD':
        return AD(tree1, tree2)
    if method == 'CD':
        return CD(tree1, tree2)

    raise ValueError(f"Method '{method}' is not supported. Choose from: {SUPPORTED_METHODS}")
