from phylogeny.ultrametric import ultrametric_tree 
from ete3 import Tree
 
def test_ultrametric_tree():
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
    
    # Get the reconstruction
    rec = ultrametric_tree(matrix, nodes)
    
    assert real.compare(rec)['rf'] == 0
# ---
