#Bourque Distance python
import numpy as np
from anytree import Node, RenderTree, PreOrderIter, LevelOrderIter, ZigZagGroupIter
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching #try this as well for k-BD
import tree
from scipy.optimize import linear_sum_assignment
import copy




def setsofchildren(root, intersection = None, lista = False):

    f = set()
    list1 = []

    if intersection == None:
        
        nodes = [node for node in PreOrderIter(root)]
        for node in nodes:
            if not node.is_root:
                list1.append(node)
            if type(node.name)!=list:
                f.add(node.name)
            else:
                for i in node.name:
                    f.add(i)
    else:
        for node in PreOrderIter(root):
            if type(node.name)!=list:
                if node.name in intersection:
                    f.add(node.name)
            else:
                for i in node.name:
                    if i in intersection:
                        f.add(i)
    if not lista:
        return f
    else:
        return f, list1


def BourqueDistance(root1, root2): 
    
    root1 = copy.deepcopy(root1)
    root2 = copy.deepcopy(root2)

    nodes1, list1 = setsofchildren(root1, lista = True)
    nodes2, list2 = setsofchildren(root2, lista = True)
    total_edges = len(list1) + len(list2) 
    
    if nodes1 == nodes2: 
        RF = comparison(root1, root2, list1, list2)
        RF += comparison(root2, root1, list2, list1)
        return RF / total_edges
    
    else:
        intersection = nodes1.intersection(nodes2) #maybe check if intersection is empty
        BD = (total_edges - min(comparison_bd(root1, root2, list1 , list2, intersection), comparison_bd(root2, root1, list2, list1, intersection))) / total_edges
        return BD


def comparison(root1, root2, list1, list2):
    diff = 0
    for mainnode1 in list1: 
        original_parent = mainnode1.parent
        mainnode1.parent = None
        edge_diff = 1
        for mainnode2 in list2:
            original_parent2 = mainnode2.parent
            mainnode2.parent = None

            desc1  = setsofchildren(mainnode1)
            asc1 = setsofchildren(root1)
            desc2  = setsofchildren(mainnode2)
            asc2 = setsofchildren(root2)

            if not len(desc1) == 0 and desc1 == desc2 and asc1 == asc2:
                mainnode2.parent = original_parent2
                edge_diff = 0
                break

            mainnode2.parent = original_parent2

        diff += edge_diff
        mainnode1.parent = original_parent
    
    return diff


def comparison_bd(root1, root2, list1, list2, intersection):
    eq = 0

    for mainnode1 in list1: 
        original_parent = mainnode1.parent
        mainnode1.parent = None
        edge_eq = 0
        for mainnode2 in list2:
            original_parent2 = mainnode2.parent
            mainnode2.parent = None

            desc1  = setsofchildren(mainnode1, intersection=intersection)
            
            asc1 = setsofchildren(root1, intersection=intersection)
            desc2  = setsofchildren(mainnode2, intersection=intersection)
            
            asc2 = setsofchildren(root2, intersection=intersection)

            if not len(desc1) == 0 and desc1 == desc2 and asc1 == asc2:
                mainnode2.parent = original_parent2
                edge_eq += 1
                break
                
            mainnode2.parent = original_parent2

        eq += edge_eq
        mainnode1.parent = original_parent
    
    return eq




def kBourqueDistance(root1, root2, k):

    root1 = copy.deepcopy(root1)
    root2 = copy.deepcopy(root2)
   
    # do the Bourque distance for every subtree with k descendants and the find the minimum weight path
    
    list1 = [node for node in PreOrderIter(root1) if not node.is_leaf] # if not node.is_leaf
    list2 = [node for node in PreOrderIter(root2) if not node.is_leaf] # if not node.is_leaf

    if len(list1)<len(list2):

        list1, list2 = list2, list1

    #matrix with numpy, filled with None values
    matrix = np.full((len(list1), len(list1)), None)

    
    for mainnode1 in list1: 

        original_parent = mainnode1.parent
        mainnode1.parent = None
        children1 = []
        tree1 = [[node for node in children] for children in ZigZagGroupIter(mainnode1, maxlevel=k+1)] 

        if len(tree1)==k+1:
            for i in tree1[-1]:
                children1.append(i.children)
                i.children = []


        

        for mainnode2 in list2:


            original_parent2 = mainnode2.parent
            mainnode2.parent = None

            children2 = []        
            tree2 = [[node for node in children] for children in ZigZagGroupIter(mainnode2, maxlevel=k+1)]

            if len(tree2)==k+1:
                for i in tree2[-1]:
                    children2.append(i.children)
                    i.children = []


            # do the Bourque distance for the two trees
            BD = BourqueDistance(mainnode1, mainnode2)
            
            matrix[list1.index(mainnode1)][list2.index(mainnode2)] = BD

            if children2 != []:
                for i in range(len(children2)):
                    tree2[-1][i].children = children2[i]

            mainnode2.parent = original_parent2
 

        if children1 != []:
            for i in range(len(children1)):
                tree1[-1][i].children = children1[i]
        mainnode1.parent = original_parent


    matrix = np.where(matrix==None, 1, matrix)

    matrix = np.array(matrix, dtype=float)
    row_ind, col_ind = linear_sum_assignment(matrix)

    return matrix[row_ind, col_ind].sum() / (len(list1))





