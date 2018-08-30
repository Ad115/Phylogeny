from phylogeny.models import CFN_Tree


def test_creation():
    n_leaves = 50
    t = CFN_Tree(n_leaves)
    
    # Tree must have the right number of leaves
    assert len(t.get_leaves()) == n_leaves
    
    # The node probabilities must be strictly in the interval [0, 0.5)
    for node in t.traverse():
        assert node.probability < 0.5
# ---

def test_evolve_traits():
    n_leaves = 50
    n = 1_000
    
    t = CFN_Tree(n_leaves)
    sequences = t.evolve_traits([1]*n)
    
    # There must be as many resulting sequences as
    # the number of leaves in the tree
    assert len(sequences) == n_leaves
    
    # The resulting sequences must have the same 
    # length as the initial sequence
    for seq in sequences.values():
        assert len(seq) == n
    
    # The sequences must be binary
    for seq in sequences.values():
        assert set(seq) ^ {1,0} == set()
