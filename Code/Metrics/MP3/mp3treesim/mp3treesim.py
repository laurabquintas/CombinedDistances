#!/usr/bin/env python3
# coding: utf-8
"""MP3 tree similarity metric (Ciccolella et al.).

Computes similarity between two labeled trees by comparing the topology
of all label triples. Supports multiple aggregation modes:
- 'sigmoid': Weighted combination of intersection and union scores.
- 'intersection': Only considers triples of shared labels.
- 'union': Considers all triples including unshared labels.
- 'geometric': Geometric mean of intersection and union scores.
"""

from functools import partial
from itertools import combinations
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx
import numpy as np
from pygraphviz import AGraph

from multiprocessing import Pool, set_start_method


class Tree:
    """MP3-treesim tree representation.

    Wraps a NetworkX directed graph with label-to-node and node-to-label
    mappings, and precomputes all-pairs LCA (Lowest Common Ancestor) data.

    Attributes:
        T: The underlying NetworkX DiGraph.
        label_to_nodes: Dict mapping each label to the set of node IDs containing it.
        node_to_labels: Dict mapping each node ID to the set of labels it contains.
        label_set: Sorted list of all unique labels.
        LCA: Precomputed LCA lookup object.
    """

    def __init__(self, T: nx.DiGraph, label_to_nodes: Dict[str, Set],
                 node_to_labels: Dict[str, Set], label_set: List[str]):
        self.T = T
        self.label_to_nodes = label_to_nodes
        self.node_to_labels = node_to_labels
        self.label_set = label_set
        self.LCA = LCA(self.T, self.label_to_nodes, self.node_to_labels)


class LCA:
    """Precomputed Lowest Common Ancestor (LCA) lookup.

    Computes all-pairs LCA at construction time for efficient triple comparison.

    Attributes:
        LCA_dict: Dict mapping (node1, node2) pairs to their LCA node.
        LCA_label_dict: Dict mapping labels to the set of nodes containing them.
        LCA_node2lbl: Dict mapping nodes to the set of labels they contain.
    """

    def __init__(self, T: nx.DiGraph, label_to_node: Dict[str, Set],
                 node_to_labels: Dict[str, Set]):
        self.LCA_dict: Dict[Tuple, str] = {}
        self.LCA_label_dict = label_to_node
        self.LCA_node2lbl = node_to_labels
        for lca in nx.tree_all_pairs_lowest_common_ancestor(T):
            self.LCA_dict[(lca[0][0], lca[0][1])] = lca[1]

    def lca_nodes(self, node1: str, node2: str):
        """Find the LCA of two nodes.

        Args:
            node1: First node ID.
            node2: Second node ID.

        Returns:
            The LCA node ID.
        """
        try:
            return self.LCA_dict[node1, node2]
        except KeyError:
            return self.LCA_dict[node2, node1]

    def lca_labels(self, label1: str, label2: str) -> Counter:
        """Compute the multiset of LCA labels for two labels.

        Args:
            label1: First label.
            label2: Second label.

        Returns:
            Counter of LCA node labels across all node pair combinations.
        """
        nodes1 = self.LCA_label_dict[label1]
        nodes2 = self.LCA_label_dict[label2]

        lca_mset = Counter()
        for n1 in nodes1:
            for n2 in nodes2:
                lca_mset.update(self.lca_nodes(n1, n2))

        return lca_mset

    def label_to_node(self, label: str) -> Set:
        """Get the set of nodes containing a given label."""
        return self.LCA_label_dict[label]

    def node_to_labels(self, node: str) -> Set:
        """Get the set of labels at a given node."""
        return self.LCA_node2lbl[node]

    def __str__(self) -> str:
        return str(self.LCA_dict)


class ExtValue:
    """Generates unique external labels for unlabeled LCA nodes.

    Used when the LCA of a triple falls on a node with no relevant labels.
    """

    def __init__(self):
        self.counter = 1
        self.nodes: Dict[str, Set[str]] = defaultdict(set)

    def get(self, node: str) -> Set[str]:
        """Get or create an external label for a node.

        Args:
            node: Node ID.

        Returns:
            Set containing the external label for this node.
        """
        if node not in self.nodes:
            self.nodes[node].add(f'EXT{self.counter}')
            self.counter += 1
        return self.nodes[node]


def sigmoid(x: float, mult: float = 10.0) -> float:
    """Sigmoid function for weighting intersection vs union similarity.

    Maps x in [0, 1] to [0, 1] with a steep transition around 0.5.
    Returns exact 0 or 1 at the boundaries.

    Args:
        x: Input value.
        mult: Steepness multiplier (higher = steeper transition).

    Returns:
        Sigmoid-transformed value.
    """
    if x == 0:
        return 0
    if x == 1:
        return 1
    return 1 / (1 + np.exp(-mult * (x - 0.5)))


def intersect_mset_card(list_lca1: List, list_lca2: List) -> int:
    """Compute the intersection cardinality of two multisets of LCA vectors.

    Args:
        list_lca1: List of LCA comparison vectors from tree 1.
        list_lca2: List of LCA comparison vectors from tree 2.

    Returns:
        Number of matching elements (minimum count for each unique vector).
    """
    mset1: Dict[str, int] = defaultdict(int)
    mset2: Dict[str, int] = defaultdict(int)

    for lca in list_lca1:
        mset1[','.join(str(x) for x in lca)] += 1
    for lca in list_lca2:
        mset2[','.join(str(x) for x in lca)] += 1

    card = 0
    for key in mset1:
        if key in mset2:
            card += min(mset1[key], mset2[key])

    return card


def is_equal_struct(triple: Tuple[str, str, str], LCA1: LCA, LCA2: LCA) -> Tuple[bool, int, int]:
    """Compare the topological structure of a label triple between two trees.

    For each pair within the triple, computes the LCA in both trees and
    builds comparison vectors. The similarity is the multiset intersection
    cardinality of these vectors.

    Args:
        triple: A tuple of three labels to compare.
        LCA1: Precomputed LCA data for tree 1.
        LCA2: Precomputed LCA data for tree 2.

    Returns:
        Tuple of (is_missing, numerator, denominator) where:
        - is_missing: True if one tree has no node combinations for this triple.
        - numerator: Number of matching LCA structures.
        - denominator: Maximum number of node combinations across both trees.
    """
    sorted_triple = sorted(triple)
    triple_set = set(triple)

    # Generate all node combinations for this triple in both trees
    triples_nodes_t1 = [
        [n1, n2, n3]
        for n1 in LCA1.label_to_node(sorted_triple[0])
        for n2 in LCA1.label_to_node(sorted_triple[1])
        for n3 in LCA1.label_to_node(sorted_triple[2])
    ]

    triples_nodes_t2 = [
        [n1, n2, n3]
        for n1 in LCA2.label_to_node(sorted_triple[0])
        for n2 in LCA2.label_to_node(sorted_triple[1])
        for n3 in LCA2.label_to_node(sorted_triple[2])
    ]

    # Build comparison vectors for each node combination
    cmp_vecs_t1 = [[None, None, None] for _ in range(len(triples_nodes_t1))]
    ext_t1 = [ExtValue() for _ in range(len(triples_nodes_t1))]

    cmp_vecs_t2 = [[None, None, None] for _ in range(len(triples_nodes_t2))]
    ext_t2 = [ExtValue() for _ in range(len(triples_nodes_t2))]

    # For each pair of labels in the triple, compute the LCA comparison
    for pair_idx, couple in enumerate(combinations(sorted_triple, 2)):
        label_idx1 = sorted_triple.index(couple[0])
        label_idx2 = sorted_triple.index(couple[1])
        other_idx = list({0, 1, 2} - {label_idx1, label_idx2})[0]

        for node_combo_idx, node_combo in enumerate(triples_nodes_t1):
            lca_node = LCA1.lca_nodes(node_combo[label_idx1], node_combo[label_idx2])
            lca_labels = LCA1.node_to_labels(lca_node)
            if len(lca_labels & triple_set) > 0:
                cmp_vecs_t1[node_combo_idx][pair_idx] = lca_labels & triple_set
            else:
                lca_with_other = LCA1.lca_nodes(lca_node, node_combo[other_idx])
                other_labels = LCA1.node_to_labels(lca_with_other)
                if len(other_labels & triple_set) > 0:
                    cmp_vecs_t1[node_combo_idx][pair_idx] = other_labels & triple_set
                else:
                    cmp_vecs_t1[node_combo_idx][pair_idx] = ext_t1[node_combo_idx].get(lca_node)

        for node_combo_idx, node_combo in enumerate(triples_nodes_t2):
            lca_node = LCA2.lca_nodes(node_combo[label_idx1], node_combo[label_idx2])
            lca_labels = LCA2.node_to_labels(lca_node)
            if len(lca_labels & triple_set) > 0:
                cmp_vecs_t2[node_combo_idx][pair_idx] = lca_labels & triple_set
            else:
                lca_with_other = LCA2.lca_nodes(lca_node, node_combo[other_idx])
                other_labels = LCA2.node_to_labels(lca_with_other)
                if len(other_labels & triple_set) > 0:
                    cmp_vecs_t2[node_combo_idx][pair_idx] = other_labels & triple_set
                else:
                    cmp_vecs_t2[node_combo_idx][pair_idx] = ext_t2[node_combo_idx].get(lca_node)

    is_missing = len(cmp_vecs_t1) == 0 or len(cmp_vecs_t2) == 0
    numerator = intersect_mset_card(cmp_vecs_t1, cmp_vecs_t2)
    denominator = max(len(cmp_vecs_t1), len(cmp_vecs_t2))

    return is_missing, numerator, denominator


def similarity(tree1: Tree, tree2: Tree, mode: str = 'sigmoid',
               sigmoid_mult: float = 10.0, cores: int = 1) -> float:
    """Compute the MP3 similarity score between two trees.

    Args:
        tree1: MP3-treesim tree representation.
        tree2: MP3-treesim tree representation.
        mode: Similarity mode — 'sigmoid', 'intersection', 'union', or 'geometric'.
        sigmoid_mult: Steepness multiplier for the sigmoid weighting.
        cores: Number of CPU cores for parallel computation. Use 1 for serial.

    Returns:
        Similarity score in [0, 1], where 1 means identical topology.
    """
    if mode not in ['sigmoid', 'intersection', 'union', 'geometric']:
        raise AttributeError(f"Invalid mode '{mode}'. Choose from: sigmoid, intersection, union, geometric.")

    # No shared labels means zero similarity
    if len(set(tree1.label_set) & set(tree2.label_set)) == 0:
        return 0.0

    # Select label set based on mode
    if mode == 'intersection':
        labels = set(tree1.label_set) & set(tree2.label_set)
    else:
        labels = set(tree1.label_set) | set(tree2.label_set)

    label_triples = combinations(labels, 3)

    # Accumulate numerator and denominators across all triples
    numerator = 0
    denominator_intersection = 0
    denominator_union = 0

    if cores is None or cores <= 0:
        cores = None

    if cores == 1:
        # Serial computation
        for triple in label_triples:
            is_missing, num, den = is_equal_struct(triple, tree1.LCA, tree2.LCA)
            numerator += num
            if is_missing:
                denominator_union += den
            else:
                denominator_intersection += den
    else:
        # Parallel computation
        try:
            set_start_method("spawn")
        except RuntimeError:
            pass  # Already set in this process

        with Pool(processes=cores) as pool:
            func = partial(is_equal_struct, LCA1=tree1.LCA, LCA2=tree2.LCA)
            results = pool.map(func, label_triples)

            for is_missing, num, den in results:
                numerator += num
                if is_missing:
                    denominator_union += den
                else:
                    denominator_intersection += den

    # Compute final similarity based on mode
    if mode == 'intersection':
        return float(numerator) / denominator_intersection
    elif mode == 'union':
        return float(numerator) / (denominator_intersection + denominator_union)
    else:
        score_intersection = float(numerator) / denominator_intersection
        score_union = float(numerator) / (denominator_intersection + denominator_union)

        if mode == 'sigmoid':
            return (score_union
                    + sigmoid(score_intersection, mult=sigmoid_mult)
                    * min(score_intersection - score_union, score_union))
        elif mode == 'geometric':
            return np.sqrt(score_intersection * score_union)

    return 0.0


def build_tree(T: nx.DiGraph, labeled_only: bool = False,
               exclude: Optional[List[str]] = None) -> Tree:
    """Build an MP3-treesim Tree from a NetworkX directed graph.

    Each node should have a ``label`` attribute with comma-separated labels.
    Nodes without labels are assigned their node ID as the label (unless
    ``labeled_only=True``, in which case they are skipped).

    Args:
        T: NetworkX directed graph representing a tree.
        labeled_only: If True, skip nodes without a ``label`` attribute.
        exclude: Labels to exclude from the tree representation.

    Returns:
        MP3-treesim Tree object.

    Raises:
        ValueError: If T is not a valid tree.
    """
    if not nx.is_tree(T):
        raise ValueError("Not a valid tree.")

    label_to_nodes: Dict[str, Set] = defaultdict(set)
    label_set: Set[str] = set()
    node_to_labels: Dict[str, Set] = defaultdict(set)

    for node in T.nodes(data=True):
        node_id = node[0]

        if 'label' not in node[1] and not labeled_only:
            node[1]['label'] = str(node[0])
        if 'label' not in node[1] and labeled_only:
            node[1]['label'] = ''
            continue

        labels = node[1]['label'].split(',')
        for label in labels:
            if exclude and label in exclude:
                continue

            label_set.add(label)
            label_to_nodes[label].add(node_id)
            node_to_labels[node_id].add(label)

    sorted_labels = sorted(list(label_set))
    return Tree(T, label_to_nodes, node_to_labels, sorted_labels)


def read_dotfile(path: str, labeled_only: bool = False,
                 exclude: Optional[List[str]] = None) -> Tree:
    """Read a DOT file and build an MP3-treesim tree.

    The tree should have a ``label`` attribute for each node. If not present,
    the node's ID is used as its label.

    Args:
        path: Path to the DOT (.gv) file.
        labeled_only: If True, skip unlabeled nodes.
        exclude: Labels to exclude from computation.

    Returns:
        MP3-treesim Tree object.
    """
    T = nx.DiGraph(nx.drawing.nx_agraph.read_dot(path))
    return build_tree(T, labeled_only=labeled_only, exclude=exclude)


def read_dotstring(string: str, labeled_only: bool = False,
                   exclude: Optional[List[str]] = None) -> Tree:
    """Parse a DOT-format string and build an MP3-treesim tree.

    Args:
        string: DOT graph as a string.
        labeled_only: If True, skip unlabeled nodes.
        exclude: Labels to exclude from computation.

    Returns:
        MP3-treesim Tree object.
    """
    T = nx.DiGraph(nx.drawing.nx_agraph.from_agraph(AGraph(string=string)))
    return build_tree(T, labeled_only=labeled_only, exclude=exclude)


def draw_tree(tree: Tree) -> None:
    """Draw the tree using NetworkX's matplotlib visualization.

    In Jupyter/Colab, the tree displays automatically. From command line,
    call ``plt.show()`` after this function.

    Args:
        tree: MP3-treesim tree representation.
    """
    drawtree = nx.convert_node_labels_to_integers(tree.T)
    labels = nx.get_node_attributes(drawtree, 'label')

    for node in drawtree.nodes(data=True):
        del node[1]['label']

    try:
        pos = nx.nx_pydot.pydot_layout(drawtree, prog='dot')
    except Exception:
        pos = None

    nx.draw_networkx(drawtree, pos=pos, labels=labels)
