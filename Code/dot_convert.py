from anytree.exporter import DotExporter
import os

#should I already use the read_dot_file as the return value?

current_directory = os.getcwd()

def dot_converter(tree):
    path = current_directory + f"/AML/aml converts/MP3 trees/dot_tree{tree.data}.gv"

    if not os.path.exists(path):
        DotExporter(tree).to_dotfile(path)

    return path
