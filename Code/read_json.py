
import json
from anytree.importer import JsonImporter
from anytree import Node, RenderTree, PreOrderIter
import pickle

def check_duplicate_labels(root):
    """
    Checks for duplicate labels in an anytree tree structure.

    Parameters:
    - root: The root node of the tree.

    Returns:
    - A list of duplicate node labels. If the list is empty, there are no duplicates.
    """
    labels = [node.label for node in PreOrderIter(root)]
    duplicates = set([label for label in labels if labels.count(label) > 1])
    return list(duplicates)

def read_json_trees(path):
    with open(path) as json_file:
        return json.load(json_file)

import re

def extract_aml_number(s):
    match = re.match(r'(AML-\d+)-\d+', s)
    if match:
        return match.group(1)
    return None

sample_map = read_json_trees("AML/aml_mutation_trees.js")
tree_map = dict()
importer = JsonImporter()
a = 200
m = 0
mut_to_number = dict()
sample_to_number = dict()
for sample_name, sample in sample_map.items():
    root = importer.import_(json.dumps(sample["event_tree"]))
    for node in PreOrderIter(root):
        
        name_child = dict()
        for child in node.children:
            if child.label not in name_child:
                name_child[child.label] = child
            else:
                if child.children != ():
                    prev_node = name_child[child.label]
                    for c in child.children:
                        prev_node.children = prev_node.children + (c,)
                child.parent = None
    
        for ascendant in node.ancestors:
            if ascendant.label == node.label:
                if ascendant.children != ():
                    prev_node = node.parent
                    for c in node.children:
                        prev_node.children = prev_node.children + (c,)     
                node.parent = None  
    

    sample_to_number[sample_name] = [a]
    extract = extract_aml_number(sample_name)
    new_sample = extract + '-001'
    if new_sample not in tree_map:
        tree_map[new_sample] = dict()
    tree_map[new_sample][a] = root
    a += 1

parallel_trees = dict()
for tree in tree_map:
    for k in tree_map[tree]:
        l = check_duplicate_labels(tree_map[tree][k])
        if len(l) > 0:
            if tree not in parallel_trees:
                parallel_trees[tree] = dict()
            parallel_trees[tree][k] = l


for tree in tree_map:
    for k in tree_map[tree]:
        trees = tree_map[tree][k]
        for node in PreOrderIter(trees):
            if tree in parallel_trees and k in parallel_trees[tree] and node.label in parallel_trees[tree][k]:
                label = node.label
                node.label = [label]
                for asc in node.ancestors:
                    if asc.label != 'Root':
                        if isinstance(asc.label, list):
                            node.label += asc.label
                        else:
                            node.label += [asc.label]
                label = node.label
                if len(label)==1:
                    node.label = label[0]
                

            if str(node.label) not in mut_to_number:
                if isinstance(node.label, list):
                    for i in mut_to_number:
                        #jaccard distance between labels in the dicitonary and the current label
                        set1 = set()
                        label = node.label
                        imp_label1 = label[0]
                        for j in label:
                            set1.add(j)
                        
                        set2 = set()

                        if i.startswith('['):
                            i = i.strip('[]')
                            lbl = list(i.split(','))
                            imp_label = lbl[0]
                            for j in lbl:
                                set2.add(j)

                        else:
                            imp_label = i
                            set2.add(i)

                        jaccard = len(set1.intersection(set2)) / len(set1.union(set2))
                        if jaccard > 0.5 and imp_label1 == imp_label:
                            mut_to_number[str(node.label)] = mut_to_number[i]
                            break
                if str(node.label) not in mut_to_number:        
                    mut_to_number[str(node.label)] = m
                    m += 1
                
        
            node.name = mut_to_number[str(node.label)]
            
            if node.label == 'Root':
                node.data = str(k)
            else:
                node.data = node.node_id
        



# dictionary to datafram and save it as csv
import pandas as pd
df = pd.DataFrame.from_dict(mut_to_number, orient='index', columns=['Index'])
df.to_csv('AML/aml_mut_map.csv')

# dictionary to datafram and save it as csv

df = pd.DataFrame.from_dict(sample_to_number, orient='index', columns=['Index'])
df.to_csv('AML/aml_sample_map.csv')

sorted_dict = dict(sorted(tree_map.items(), key=lambda kv: int(kv[0].split('-')[1])))

#save the tree map in pickle file
with open('AML/tree_map.pickle', 'wb') as handle:
    pickle.dump(sorted_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

