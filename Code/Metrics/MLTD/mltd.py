import os
import subprocess

current_directory = os.getcwd()

def calculate_mltd(path_tree1, path_tree2):
    """
    This function calculates the similarity between two trees using the MLTD algorithm.

    Parameters:
    path_tree1 (str): The file path to the first tree.
    path_tree2 (str): The file path to the second tree.

    Returns:
    float: The similarity score between the two trees.

    """
    # Construct the arguments for the command to be run in the subprocess
    args = [current_directory + "/MLTD/main", path_tree1, path_tree2]
    
    # Run the command, capturing the output and converting it to text
    result = subprocess.run(args, capture_output=True, text=True)
    
    # Process the output, which is expected to be multiple lines of text
    output_lines = result.stdout.strip().split("\n")
    
    # The last line contains the similarity score, which is the last word on the line
    similarity = float(output_lines[-1].split()[-1])
    
    return similarity