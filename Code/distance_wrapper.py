from dot_convert import dot_converter
from mltd_convert import mltd_converter
from newick_convert import newick_converter
import BD.BDpython as BD
import MP3.mp3treesim.mp3treesim as MP3
import CASetDISC.CASet as CASet
import CASetDISC.DISC as DISC
import MLTD.mltd as MLTD
from anytree import PreOrderIter
from GraPhyC.PD import PD
from GraPhyC.PCD import PCD
from GraPhyC.AD import AD
from GraPhyC.CD import CD

def check_minus3nodes(tree) -> bool:
    if len({node for node in PreOrderIter(tree)})<3:
        return True
    else:
        return False

def check_minus3nodesintersection(tree1, tree2) -> bool:
    nodes1 = {node.name for node in PreOrderIter(tree1)}
    nodes2 = {node.name for node in PreOrderIter(tree2)}
    if len(nodes1.intersection(nodes2)) <3:
        return True
    else:
        return False
    
def check_minus2nodesintersection(tree1, tree2) -> bool:
    nodes1 = {node.name for node in PreOrderIter(tree1)}
    nodes2 = {node.name for node in PreOrderIter(tree2)}
    if len(nodes1.intersection(nodes2)) < 2:
        return True
    else:
        return False

def distance(tree1, tree2, method) -> float:

    if method=='MP3_union' and check_minus3nodes(tree1) and check_minus3nodes(tree2):
        return float(1)
    
    elif method =='MP3_union':
        dot1 = dot_converter(tree1)
        dot2 = dot_converter(tree2)
        return 1 - MP3.similarity(MP3.read_dotfile(dot1), MP3.read_dotfile(dot2), mode='union')
    
    elif method.startswith('MP3') and (check_minus3nodes(tree1) or check_minus3nodes(tree2) or check_minus3nodesintersection(tree1,tree2)):
        return float(1)

    elif method == 'BD':
        return BD.BourqueDistance(tree1, tree2)
    
    elif method == 'MP3_sigmoid':
        dot1 = dot_converter(tree1)
        dot2 = dot_converter(tree2)
        return 1 - MP3.similarity(MP3.read_dotfile(dot1), MP3.read_dotfile(dot2)) 
    
    elif method =='MP3_geometric':
        dot1 = dot_converter(tree1)
        dot2 = dot_converter(tree2)
        return 1 - MP3.similarity(MP3.read_dotfile(dot1), MP3.read_dotfile(dot2), mode='geometric')
    
    elif method =='MP3_intersection':
        dot1 = dot_converter(tree1)
        dot2 = dot_converter(tree2)
        return 1 - MP3.similarity(MP3.read_dotfile(dot1), MP3.read_dotfile(dot2), mode='intersection')
    
    elif method == 'CASet_intersection':
        newick1 = newick_converter(tree1)
        newick2 = newick_converter(tree2)
        return CASet.caset_intersection(newick1, newick2)
    
    elif method == 'CASet_union':
        newick1 = newick_converter(tree1)
        newick2 = newick_converter(tree2)
        return CASet.caset_union(newick1, newick2)
    
    elif method == 'DISC_intersection':
        newick1 = newick_converter(tree1)
        newick2 = newick_converter(tree2)
        return DISC.disc_intersection(newick1, newick2)
    
    elif method == 'DISC_union':
        newick1 = newick_converter(tree1)
        newick2 = newick_converter(tree2)
        return DISC.disc_union(newick1, newick2)
    
    elif method == 'MLTD':
        str_t1 = mltd_converter(tree1)
        str_t2 = mltd_converter(tree2)
        return 1 - MLTD.calculate_mltd(str_t1, str_t2) 
    
    elif method == 'PD':
        if check_minus2nodesintersection(tree1, tree2):
            return float(1)
        
        dot1 = dot_converter(tree1)
        dot2 = dot_converter(tree2)
        
        return PD(dot1, dot2)
    
    elif method == 'PCD':
        return PCD(tree1, tree2)
    
    elif method == 'AD':
        return AD(tree1, tree2)
    
    elif method == 'CD':
        return CD(tree1, tree2)
    
    elif method == '1BD':
        return BD.kBourqueDistance(tree1, tree2, k=1)
    
    elif method == '2BD':
        return BD.kBourqueDistance(tree1, tree2, k=2)

    else:
        raise ValueError('Method not supported')