"""Path Distance (PD) metric for comparing mutation trees.

Computes the normalized difference in shortest-path distances between
all pairs of shared nodes in two trees. Trees are read from DOT files
and converted to NetworkX graphs.

Each pairwise distance difference is normalized by the maximum distance
observed for that pair minus 1 (to account for the minimum possible
distance of 1 in a connected tree).
"""

from typing import Dict, Set, Tuple

import networkx as nx


def dot_to_graph(dotfile: str) -> nx.Graph:
    """Read a DOT file and return an undirected NetworkX graph.

    Args:
        dotfile: Path to the .gv (GraphViz DOT) file.

    Returns:
        Undirected graph parsed from the DOT file.
    """
    return nx.Graph(nx.nx_pydot.read_dot(dotfile))


def all_pairs_distances(graph: nx.Graph) -> Dict:
    """Compute shortest path lengths between all pairs of nodes.

    Args:
        graph: An undirected NetworkX graph.

    Returns:
        Nested dict: distances[node1][node2] = shortest_path_length.
    """
    return dict(nx.all_pairs_shortest_path_length(graph))


def PD(file1: str, file2: str) -> float:
    """Compute the Path Distance between two trees given as DOT files.

    For each pair of nodes present in both trees, computes the normalized
    absolute difference in shortest-path distances.

    Args:
        file1: Path to the first DOT file.
        file2: Path to the second DOT file.

    Returns:
        Average normalized path distance across all compared node pairs.
    """
    graph1 = dot_to_graph(file1)
    graph2 = dot_to_graph(file2)

    distances1 = all_pairs_distances(graph1)
    distances2 = all_pairs_distances(graph2)

    total_difference = 0.0
    num_compared = 0
    seen_pairs: Set[Tuple] = set()

    for node_a in distances1:
        for node_b in distances1[node_a]:
            if node_a == node_b:
                continue

            pair = (min(node_a, node_b), max(node_a, node_b))
            if pair in seen_pairs:
                continue

            # Find the matching distance in graph2 (check both orderings)
            dist1 = distances1[node_a][node_b]
            dist2 = None

            if node_a in distances2 and node_b in distances2.get(node_a, {}):
                dist2 = distances2[node_a][node_b]
            elif node_b in distances2 and node_a in distances2.get(node_b, {}):
                dist2 = distances2[node_b][node_a]

            if dist2 is not None:
                # Normalize by max distance - 1 (minimum connected distance is 1)
                max_dist = max(dist1, dist2)
                normalizer = max_dist - 1 if max_dist > 1 else 1
                total_difference += abs(dist1 - dist2) / normalizer
                num_compared += 1

            seen_pairs.add(pair)

    return total_difference / num_compared
