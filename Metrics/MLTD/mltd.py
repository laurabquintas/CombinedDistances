import os
import subprocess

current_directory = os.getcwd()

def calculate_mltd(path_tree1,path_tree2):
    args = [current_directory + "/MLTD/main", path_tree1, path_tree2]
    result = subprocess.run(args, capture_output=True, text=True)
    output_lines = result.stdout.strip().split("\n")
    similarity = float(output_lines[-1].split()[-1])
    return similarity

