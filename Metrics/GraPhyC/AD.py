from anytree import PreOrderIter


def AD(tree1, tree2):

    pairs1 = create_pairs(tree1)
    pairs2 = create_pairs(tree2)

    return 1 - (len(pairs1.intersection(pairs2)) / len(pairs1.union(pairs2)))



def create_pairs(tree):
    pairs = set()
    for node in PreOrderIter(tree):
        if not node.is_leaf:
            for child in node.descendants:
                pairs.add((node.name, child.name))
    return pairs