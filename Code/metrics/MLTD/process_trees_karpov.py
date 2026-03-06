import pandas as pd
import os

current_directory = os.getcwd()
# Read the CSV file into a DataFrame
df = pd.read_csv(current_directory + "/Cancer data/brca_trees.csv", sep = ',')

# Create a separate text file for each tree
tree_ids = df['Tree_ID'].unique()
for tree_id in tree_ids:
    # Filter the DataFrame to include only the rows for this tree
    tree_df = df[df['Tree_ID'] == tree_id]

    # Initialize the dictionary for this tree
    labels = {}
    children = {}

    # Iterate through each row in the DataFrame
    for index, row in tree_df.iterrows():
        # Parse the row
        
        node_id = 'N'+str(row['Node_ID'])
        
        label = str(row['Mutation_ID'])
        
        parent_id = 'N'+str(row['Parent_ID'])
        
        
        labels[node_id] = []
        # If this is the first node in the tree, set it as the root
        if node_id == 'N1':
            labels[node_id].append(label)
            children[node_id] = []
        else:
            # Add this node to the dictionary
            labels[node_id].append(label)

            # Add this node as a child of its parent
            
            if parent_id not in children:
                children[parent_id] = []
            children[parent_id].append(node_id)

    # Create the output string for this tree
    output_string = ''
    for Node_ID in labels:
        output_string += Node_ID + '=' + ','.join(labels[Node_ID]) + '\n'
    for parent_id in children:
        output_string += parent_id + ':' + ','.join(children[parent_id]) + '\n'


    # Write the output string to a text file
    with open(current_directory + '/Karpov/' + f'tree_{tree_id}.txt', 'w') as outfile:
        outfile.write(output_string)