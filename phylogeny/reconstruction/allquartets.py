"""
All quartets method.

Given an n×n additive matrix M with n ≥ 5 associated to a 
binary tree T with positive branch lengths, we can construct
T using a two-step technique that we now describe. 

In Step 1, we compute a quartet tree on every four leaves by
applying the Four Point Method to each 4×4 submatrix of M. 

In Step 2, we assemble the quartet trees into a tree on the
full set of leaves. Step 1 is straightforward. The technique 
we use in Step 2 is called the “All Quartets Method”.

"""

import itertools as itr
from ..core import Tree
from ..core.fpc import fpc_sums


def induced_quartet(dist_matrix, idx_quartet=None):
    """Get the induced quartet ordering of 4 items."""
    if idx_quartet is None:
        idx_quartet = range(4) # The first 4 elements
    q = tuple(idx_quartet)
    
    # Calculate the relevant pairwise sums
    sums = fpc_sums(dist_matrix, q)
    # Get the quartet with smallest sum
    quartet = min(sums, key=lambda x:sums[x])
    
    return quartet
# ---

def map_names_to_quartet(quartet, names=None):
    "Map the names to the quartet's indices."
    if names:
        ((a,b),(c,d)) = quartet
        return ((names[a],names[b]), (names[c],names[d]))
    else:
        return quartet
# ---

def four_point_method(additive, names=None):
    """Method for inferring a tree from a 4x4 additive matrix.
    
    If we are given a 4×4 additive matrix 'D' that 
    corresponds to a tree 'T' with positive branch weights, 
    then we can easily compute 'T' from 'D': We calculate 
    the three pairwise sums from the four point condition, 
    we determine which of the three pairwise sums is the 
    smallest, and use that one to define the split for the 
    four leaves into two sets of two leaves each.
    """
    if names is None:
        try:
            names = additive.names
        except AttributeError:
            pass
    
    # Calculate the quartet inferred by the distances
    quartet = induced_quartet(additive)
    
    # Map the names to the quartet
    quartet = map_names_to_quartet(quartet, names)
    
    # Assemble the quartet into a tree structure
    tree = Tree.from_quartet(quartet)
    
    return tree
# ---

def all_quartets(dist_matrix, names=None):
    "Get all inferred quartet subtrees."
    if names is None:
        try:
            names = dist_matrix.names
        except AttributeError:
            pass
    
    n = len(dist_matrix)
    quartets = itr.combinations(range(n), 4)
    
    return [map_names_to_quartet(induced_quartet(dist_matrix,q), 
                                 names)
             for q in quartets]
# ---

def infer_siblings(quartets):
    """From the tree quartets, infer pairs of sibling leafs.
    
    We search for a pair x,y of leaves that is always 
    together in any quartet that contains both x and y. 
    (In other words, for all a,b, any quartet on {x,y,a,b} 
    is ((x,y),(a,b))). Any pair of leaves that are siblings 
    in the quartets tree T will satisfy this property.
    """
    together = set()
    separated = set()
    
    for q in quartets:
        quartet = { frozenset(pair) for pair in q }

        together |= quartet

        ((a,b), (c,d)) = q
        separated |= { frozenset(i) 
                          for i in [(a,c), (a,d), 
                                    (b,c), (b,d)] }

    return {frozenset(pair) for pair in together - separated}
# ---

def tree_from_quartets(quartets):
    "From the given quartets, assemble the tree."
    if len(quartets) == 1:
        q = quartets[0]
        return Tree.from_quartet(q)
    else:
        # Fetch a pair of sibling leafs
        a,b = list(infer_siblings(quartets).pop())
            
        # Recourse in quartets \ {a}
        new_quartets = [q for q in quartets if (a not in q[0]) and (a not in q[-1])]
        tree = tree_from_quartets(new_quartets)
        
        # Add a as sibling of b
        tree.add_as_sibling(a,b)
        return tree
# ---

def all_quartets_method(dist_matrix, names=None):
    "Reconstruct the tree from the dist. matrix using the all quartets method."
    if names is None:
        try:
            names = dist_matrix.names
        except AttributeError:
            pass
    quartets = all_quartets(dist_matrix, names)
    return tree_from_quartets(quartets)
# ---