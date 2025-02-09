import pandas as pd
import os
import pickle
from tree import TreeNode
current_directory = os.getcwd()

# do I need to include locks? if the Euler works paralelly, it might be a problem with the load and dump


def tree_to_newick(nodes):
    items = []
    for node in nodes:
        s = ''
        if len(node.children) > 0:
            sub_tree = tree_to_newick(node.children)
            if sub_tree != '':
                s += '(' + sub_tree + ')'
        if isinstance(node.name, list):  # Check if node has multiple labels, in this case {label1,label2}
            s += '{' + ','.join(str(label) for label in node.name) + '}'
        else:
            s += str(node.name)
        items.append(s)
    return ','.join(items)

#create locks

def newick_converter(tree):

    newick_trees = {}

    if os.path.exists(current_directory + "/AML/aml converts/newick_trees.pickle"):
        with open(current_directory + "/AML/aml converts/newick_trees.pickle", "rb") as myfile:
            newick_trees = pickle.load(myfile)
    
    if tree.data in newick_trees.keys():
        return newick_trees[tree.data]
    
    else:
        newick_trees[tree.data] = tree_to_newick([tree])

        with open(current_directory + "/AML/aml converts/newick_trees.pickle", "wb") as myfile:
            pickle.dump(newick_trees, myfile)

        return newick_trees[tree.data]
    

