# %load phylogeny/core/tree.py

import ete3
import numpy as np
import itertools as itr
from .distance import DistanceMatrix

class Tree(ete3.Tree):
    """Wrapper for the ETE Tree class adapted for the 
    purposes of phylogeny reconstruction.
    
    For documentation for the ETE Tree, read the `ETE3 documentation`_.
    
    .. _ETE3 documentation:
       http://etetoolkit.org/docs/latest/tutorial/index.html
    """
    
    def __init__(self, *args, leaves=None, **kwargs):
        """
        Args:
            leaves (int, optional): Populate randomly with this number of leaves.
        """
        super().__init__(*args, **kwargs)
        
        # Populate until you have as many leaves
        if leaves:
            self.populate(leaves)
    # ---
    
    def __repr__(self):
        return (  self.__class__.__name__
                + "('" + self.write(features=[]) + "')" )
    # ---
    
    @classmethod
    def from_tree(cls, tree, *args, **kwargs):
        "Create a new instance based on an existing tree."
        return cls(newick=tree.write(features=[]))
    # ---
    
    @classmethod
    def from_newick(cls, newick, *args, **kwargs):
        "Read from the newick representation."
        return cls(newick, *args, **kwargs)
    # ---
    
    @classmethod
    def from_quartet(cls, quartet):
        "Transform the quartet structure to a tree."
        # Disassemble the quartet structure
        ((a,b),(c,d)) = quartet
        # Assemble the pairs of siblings
        ab_cherry = cls.make_cherry_of(a,b)
        cd_cherry = cls.make_cherry_of(c,d)
        # Assemble the subtrees
        quartet_tree = cls.join_trees(ab_cherry, cd_cherry)

        return quartet_tree
    # ---
    
    @classmethod
    def join_trees(cls, a, b):
        "Make a 'cherry' of the two trees a and b."
        tree = cls()
        tree.add_child(a)
        tree.add_child(b)
        return tree
    # ---
    
    @classmethod
    def make_cherry_of(cls, a,b):
        "Get a cherry tree out of both items."
        cherry = cls()
        cherry.add_child(name=a)
        cherry.add_child(name=b)
        return cherry
    # ---
    
    @staticmethod
    def replace_node(old, new):
        old_parent = old.up
        old.detach()
        old_parent.add_child(new)
        return new
    # ---
    
    def show(self, mode=None, inline=False, styling=None, **kwargs):
        "Display the tree."
        
        # Check if a tree styling dict was specified
        if styling:
            # Assemble the treestyle object
            ts = ete3.TreeStyle()
            for key,value in styling.items():
                setattr(ts, key, value)
                
            # Update the arguments list
            kwargs['tree_style'] = ts
                
        # Inline Jupyter output
        if inline or (mode == 'inline'):
            return self.render('%%inline', **kwargs)
        else:
            # Tree GUI rendering
            return super().show(**kwargs)
    # ---
    
    def total_nodes(self):
        return len(list(self.traverse()))
    # ---

    def add_as_sibling(self, a, b):
        "Add leaf a as sibling of b in the tree."
        # Find the node corresponding to b and add a as sibling
        b_node = self.search_nodes(name=b)[0]
        # Make a new cherry out of a and b
        # and attach it in place of b
        cherry = self.make_cherry_of(a,b)
        self.replace_node(b_node, cherry)
    # ---

    def prune_leaves(self, to_stay):
        'Prune tree branches to leave only the leaves in `to_stay`.'
        t = self.tree
        
        # Fetch leaf nodes with those names
        stay_nodes = set()
        for leaf in t:
            if leaf.name in to_stay:
                stay_nodes.add(leaf)
        
        t.prune(stay_nodes, preserve_branch_length=True)
    # ---
    
    def distance_matrix(self):
        "Get the matrix of distances between each pair of leaves."
        # Fetch the tree's leaves
        leaves = self.get_leaves()
        n = len(leaves)
        # Get all pairs of leaves
        pairs = itr.combinations(range(n), 2)
        # Create an empty matrix
        distances = DistanceMatrix.zeros(n, names=[l.name for l in leaves])
        # Fill the distance matrix
        for i,j in pairs:
            d = self.get_distance(leaves[i],
                                  leaves[j])
            distances[i,j] = distances[j,i] = d 
            
        return distances
    # ---
# --- Tree