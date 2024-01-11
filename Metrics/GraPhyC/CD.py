from anytree import PreOrderIter

#not prepared for multilabeled nodes

def CD(tree1, tree2):
    paths1 = create_paths(tree1)
    paths2 = create_paths(tree2)
    return 1 - (len(paths1.intersection(paths2)) / len(paths1.union(paths2)))


def create_paths(tree):
    paths = set()
    for node in PreOrderIter(tree):
        path = []
        path.append(node.name)
        for nodes in node.ancestors:
            if isinstance(nodes.name,list):
                for i in nodes.name:
                    path.append(i)
            path.append(nodes.name)
        paths.add(order(path))
    return paths


#ooutput of the path as a ordered tuple
def order(path):
    if len(path) == 1:
        return tuple(path)
    else:
        return tuple(sorted(path))
    


