
import distance_wrapper as dist
import os
import pickle
import pandas as pd
import sys
from anytree import RenderTree, Node

current_directory = os.getcwd()

with open(current_directory + "/AML/tree_map.pickle", "rb") as f:
    tree_map = pickle.load(f)

methods = ['DISC_union', 'DISC_intersection', 'CASet_intersection', 'CASet_union', 'BD','1BD', '2BD', 'AD', 'CD', 'PD', 'PCD', 'MLTD', 'MP3_union', 'MP3_sigmoid', 'MP3_intersection', 'MP3_geometric']


for method in methods:
    matrix = [[None for _ in range(len(tree_map))] for _ in range(len(tree_map))]
    for i, (patient_id_1, k) in enumerate(tree_map.items()):
        print(f'{method}: i = {i}')
        for j, (patient_id_2, l) in enumerate(tree_map.items()):
            
            if matrix[i][j] is None:
                result = set()
                
                if i == j:
                    matrix[i][j] = float(0)
                
                else:
                    for t1 in k.values():

                        for t2 in l.values():   

                            distance = dist.distance(t1, t2, method)
                            result.add(distance)
                
                #if i == j and list(result).count(0) == 1:
                    #matrix[i][j] = float(0)
                
                    d = round(sum(result) / len(result), 3)
                    matrix[i][j] = d

                matrix[j][i] = matrix[i][j]

    df = pd.DataFrame(matrix, columns=list(tree_map.keys())) 
    df.to_csv(f'AML/Matrix Output/AML_{method}.csv', index=False)

