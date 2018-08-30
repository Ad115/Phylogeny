from phylogeny.reconstruction import infer_clocklike_tree1, infer_clocklike_tree2
from phylogeny import Tree
 
# A representation of the following tree:
#         /-A
#      /-|
#   /-|   \-E
#  |  |
#--|   \-D
#  |
#  |   /-B
#   \-|
#      \-C
real = Tree('(((A,E),D),(B,C));')
    
# An ultrametric matrix representing the tree
matrix = [ [0, 8, 8, 5, 3],
           [8, 0, 3, 8, 8],
           [8, 3, 0, 8, 8],
           [5, 8, 8, 0, 5],
           [3, 8, 8, 5, 0] ]
nodes = ['A', 'B', 'C', 'D', 'E'] 


def test_infer_clocklike_tree1():
    # Get the reconstruction
    rec = infer_clocklike_tree1(matrix, nodes)
    
    assert real.compare(rec)['rf'] == 0
# ---

def test_infer_clocklike_tree2():
    # Get the reconstruction
    rec = infer_clocklike_tree2(matrix, nodes)
    
    assert real.compare(rec)['rf'] == 0
# ---
