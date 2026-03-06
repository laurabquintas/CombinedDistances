'''

For the visualization:

Add the mutations to the node attributes:
{"node_id": 2, "matching_label": 19, "mutation_label": "BAP1;VHL;12p11.21;3p25.3_del;8p23.2_del;9p21.3_del", "name": 2, "node_depth": 1 (root)}

Ouptput the clusters in this format:
{"cluster1": {"trees":[list of anytrees], "metadata":{pid1:{...}, pid2:{...}}

Dict example for metadata:
"metadata": {"K375": {"Overall Stage": "#F9E0C1", "Necrosis": "white", "Size of primary tumour (mm)": "#4FC3F7", "Outcome": "
#F9CB9C", "Cancer related death": "white", "Progression group": "white", "ITH Index": "#CE93D8"}, "K110": {....}}


iterator = LevelOrderIter(root) node.depth
'''

import os
import pickle
import pandas as pd
import sys
from anytree.exporter import JsonExporter, DictExporter
from anytree import Node, RenderTree, PreOrderIter

current_directory = os.getcwd()

with open(current_directory + "/AML/tree_map.pickle", "rb") as f:
    tree_map = pickle.load(f)

clinical_data = pd.read_csv(current_directory + "/AML/GEN/Results/survival/t = 0.02/clinicaldata_wclusterlabels.csv", sep = ',')

metric = 'DISC_intersection'

clinical = ['sample','PriorMalig', 'Diagnosis', 'Gender','VitalStatus'] + [metric]

color_code_Vital_Status = {'Alive NOS' : '#228B22' , 'Dead NOS' : '#B22222'}
color_code_Gender = {'Male' : '#F1C40F',  'Female' : '#D35400'}
color_code_Dia = {'AML' : '#98af6a', 'AMML' : '#7d8eb0', 'AEL' : '#b0675a', 'AMOL': '#633f38', 'AUL': 'grey'}
color_code_PriorMalig = {'No' : '#94c9f3', 'Yes' : '#faa7a2'}

# Replace values in each column with corresponding colors
clinical_data['VitalStatus'].replace(color_code_Vital_Status, inplace=True)
clinical_data['Gender'].replace(color_code_Gender, inplace=True)
clinical_data['Diagnosis'].replace(color_code_Dia, inplace=True)
clinical_data['PriorMalig'].replace(color_code_PriorMalig, inplace=True)


df = clinical_data[clinical]  


output = dict()

for cluster in list(df[metric].unique()):
    output[f'cluster{cluster}'] = {"trees": [], "metadata": {}}
    df_cluster = df[df[metric] == cluster]
    df_cluster = df_cluster.drop(columns = [metric], axis = 1)
    

    #turn df_cluster into a dictionary and index as the patient id
    new_df_cluster = df_cluster.copy()
    new_df_cluster = new_df_cluster.set_index('sample')

    df_cluster_json = new_df_cluster.to_dict('index')

    output[f'cluster{cluster}']['metadata'] = df_cluster_json

    for patient in list(df_cluster['sample']):
        #get the tree
        if isinstance(tree_map[patient],dict):
            for i in tree_map[patient]:
                tree = tree_map[patient][i]
                a = 0
                for node in PreOrderIter(tree):
                    #save_label = node.label
                    #node.label = node.matching_label
                    #node.matching_label = save_label
                    node.node_depth = node.depth
                    node.node_id = a
                    a += 1
                    node.label = str(node.label)
                    if node.label != 'Root':
                        node.size_percent = 0.2
                    node.name = str(node.name)
                    del node.data 
                    
                exporter = DictExporter()
                dict_tree = exporter.export(tree)
                output[f'cluster{cluster}']['trees'].append(dict_tree)
        else:
            tree = tree_map[patient]
            a = 0
            for node in PreOrderIter(tree):
                #save_label = node.label
                #node.label = node.matching_label
                #node.matching_label = save_label
                node.node_depth = node.depth
                node.node_id = a
                a += 1
                node.label = str(node.label)
                if node.label != 'Root':
                    node.size_percent = 0.2
                node.name = str(node.name)
                del node.data 
                
            exporter = DictExporter()
            dict_tree = exporter.export(tree)
            output[f'cluster{cluster}']['trees'].append(dict_tree)
        

#export the output to text .txt file
with open(current_directory + "/AML/Results/survival/t = 0.02/cluster_analysis_AML_DISCint.json", "w") as f:
    f.write('sample_map=' + str(output))







