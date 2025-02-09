from anytree import Node, RenderTree, PreOrderIter


def PCD(tree1, tree2):

    pairs1 = create_pairs(tree1)
    pairs2 = create_pairs(tree2)
    int = pairs1.intersection(pairs2)
    union = pairs1.union(pairs2)
    return 1 - (len(int) / len(union))


def create_pairs(tree):
    pairs = set()
    for node in PreOrderIter(tree):
        if not node.is_leaf:
            for child in node.children:
                pairs.add((node.name, child.name))
    return pairs


