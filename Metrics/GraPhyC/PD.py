
import networkx as nx

# in the wrapper file we need to verify if they have at least 2 muations in common, if not distance equals 1 - done


# came up with the strategy to normalize by dividing by the max distance between two nodes for each comparison - confirm with Xiang and Monica


#transform dot file to graph with networkx
def dot2graph(dotfile):
    """Reads a dot file and returns a graph"""
    import networkx as nx
    G = nx.Graph(nx.nx_pydot.read_dot(dotfile))
    return G


#calculate the shortest path between two nodes
def all_pairs(G):
    """Calculates the shortest path between two nodes"""
    return dict(nx.all_pairs_shortest_path_length(G))


#calculate the absolute difference between two graphs
def PD(file1, file2):
    """Calculates the absolute difference between two graphs"""
    
    G1 = dot2graph(file1)
    G2 = dot2graph(file2)

    graph1 = all_pairs(G1)
    graph2 = all_pairs(G2)
    difference = 0
    seen_pairs = set()
    compared = 0

    distances = []
    
    # order the nodes in a tuple
    def order(node1, node2):
        if node1 < node2:
            return (node1, node2)
        else:
            return (node2, node1)

    for key in graph1.keys():
        for value in graph1[key].keys():
            
            if key!=value and order(key,value) not in seen_pairs:

            
                
                if key in graph2 and value in graph2[key]:
                    div = ( max(graph1[key][value], graph2[key][value]) - 1 )
                    if div == 0:
                        div = 1
                    difference += abs(graph1[key][value] - graph2[key][value]) / div
                    compared +=1


                elif value in graph2 and key in graph2[value]:
                    div = ( max(graph1[key][value], graph2[value][key]) - 1 )
                    if div == 0:
                        div = 1
                    difference += abs(graph1[key][value] - graph2[value][key]) / div
                    compared +=1

            seen_pairs.add(order(key,value))
    
    return difference / compared




