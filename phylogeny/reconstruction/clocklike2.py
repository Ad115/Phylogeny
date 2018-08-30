"""
Module implementing phylogeny estimation assuming clocklike evolution.

"An assumption that is sometimes made is that sequence 
evolution is clocklike (also referred to as obeying 
the strict molecular clock), which means that the 
expected number of changes is proportional to time. 
If we assume that the leaves represent extant (i.e., 
living) species, then under the assumption of a 
strict molecular clock, the total expected number of 
changes from the root to any leaf is the same. Under 
the assumption of a strict molecular clock, the matrix 
of expected distances between the leaves in the tree 
has properties that make it “ultrametric”."

    -- From the book: "Computational Phylogenetics. An introduction 
       to designing methods for phylogeny estimation" by Tandy Warnow
"""
from ..core import DistanceMatrix, Tree


def infer_clocklike_tree2(distances):
    """Assumming the sequences evolved in a clocklike process, 
    infer the tree.
    
    "One very natural approach to estimating the tree would 
    be to select as siblings the two sequences that are the 
    most similar to each other from the three possible pairs. 
    Because the sequence evolution model is clocklike, this 
    technique will correctly construct rooted three-leaf trees 
    with high probability. Furthermore, the method can even be 
    extended to work  on  trees  with  more  than  three  leaves,  
    using  recursion...
    
    Hence, to estimate this tree, we would first compare all 
    pairs of sequences to find which pair is the most similar, 
    and we’d select 'a' and 'b' as this pair. We’d then correctly 
    infer that species 'a' and 'b' are siblings. We could then 
    remove one of these two sequences (say 'a'), and reconstruct 
    the tree on what remains. Finally, we would add 'a' into the 
    tree we construct with the other vertices, by making it a 
    sibling to 'b'."
    
        -- From the book.
    """
    if len(distances) == 2:
        # Return a cherry tree
        cherry = Tree.make_cherry_of(*distances.names)
        return cherry
    else:
        # Find closest taxa a,b in S
        (a,b), _ = min(distances.name_all(), key=lambda item: item[-1])
        # Recurse on (sequences \ a)
        chopped = distances.remove(a)
        tree = infer_clocklike_tree2(chopped)
        # Find the node corresponding to b and add a as sibling
        tree.add_as_sibling(a, b)
        return tree
# ---