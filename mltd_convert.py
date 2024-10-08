import os
from anytree import Node, RenderTree, PreOrderIter
import pickle
import pandas as pd

"""
## Tree Structure
We represent trees in the following way:

```
Root=Label1,Label2, ...
Node1=Label3,Label4, ...
Node2=Label5,Label6, ...
Root:Node1,Node2
```
Each line of input contains either '=' or ':'. If a line contains '=', it means that left side of this sign represents a node and right side of it represents set of assigned labels to this node. If a line contains ':', it means that left side of it contains a node of the tree and right side of it contains children of that node (similar to adjacent list). First line that contains ':', represents the root node and its corresponding children. 

Check t1.txt file for an example. t1.txt, t2.txt, t3.txt and t4.txt are example trees used in the paper which represent 'true tree', 'inferred tree 1', 'inferred tree 2' and 'inferred tree 3', respectively. **N.B. Do not put a blank line in the tree file.**
"""

current_directory = os.getcwd()


def path_mltd(tree):
    node_lines = []
    edge_lines = []

    for node in PreOrderIter(tree):
        if not isinstance(node.name, list):
            node_labels = str(node.name)
        else:
            node_labels = ','.join(map(str, node.name))
        
        if node.is_root:
            node_lines.append(f'Root={node_labels}')
        else:
            node_lines.append(f'N{node.data}={node_labels}')

        if not node.is_leaf:
            if len(node.children)==1:
                node_children = f'N{node.children[0].data}'
            else:
                node_children = ','.join(f'N{child.data}' for child in node.children)
            if node.is_root:
                edge_lines.append(f'Root:{node_children}')
            else:
                edge_lines.append(f'N{node.data}:{node_children}')
   
    return '\n'.join(node_lines + edge_lines)


def mltd_converter(tree):
    path = current_directory + f"/AML/aml converts/MLTD trees/str_tree{tree.data}.txt"

    if not os.path.exists(path):
        with open(path, 'w') as f:
            f.write(path_mltd(tree))

    return path