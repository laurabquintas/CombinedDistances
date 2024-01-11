# o ficheiro output é um csv, ver código do DiNardo

'''distance_measure = disc_union if args.union else disc_intersection
distance_matrix = [[None for x in range(len(trees))] for y in range(len(trees))]
for i, t1 in enumerate(trees):
    for j, t2 in enumerate(trees):
        distance_matrix[i][j] = distance_measure(t1, t2)
        distance_matrix[j][i] = distance_matrix[i][j]

output = get_matrix_string(distance_matrix)

if args.minmax:
    output += get_min_max_string(distance_matrix)

if args.pickle:
    with open(args.pickle, 'wb') as out:
        pickle.dump(distance_matrix, out, pickle.HIGHEST_PROTOCOL)

if args.outputFile:
    with open(args.outputFile, 'w') as out:
        out.write(output)'''
import os
import panda as pd

def get_matrix_string(matrix):
    """
    Construct a human-readable string displaying a distance matrix.
    :param matrix: the matrix to display
    :return: a string displaying the matrix
    """
    string = 'tree\t' + '\t'.join(map(str, range(len(matrix)))) + '\n\n'
    for index, row in enumerate(matrix):
        string += '{}\t'.format(index) + '\t'.join(str(round(entry, 4)) for entry in row) + '\n'

    return string

current_directory = os.getcwd()
file = current_directory + '/Karpov/main.cpp'

df = pd.read_csv(current_directory + "/Cancer data/brca_trees.csv", sep = ',')
tree_ids = df["Tree_ID"].unique()

!g++ -std=c++11 file -o main

distance_matrix = [[None for x in range(len(tree_ids))] for y in range(len(tree_ids))]

for i in range(len(tree_ids)):
    t1= current_directory + f'/Karpov/Trees/tree_{i+1}.txt'
    for j in range(len(tree_ids)):
        t2= current_directory + f'/Karpov/Trees/tree_{j+1}.txt'
        distance_matrix[i][j] = !./main t1.txt t2.txt
        distance_matrix[j][i] = distance_matrix[i][j]

output = get_matrix_string(distance_matrix)
with open(current_directory+'/MLTD_distances.csv', 'w') as out:
        out.write(output)
        
    
