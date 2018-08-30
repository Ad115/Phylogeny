from phylogeny.core.distance import simpledistance

def test_simpledistance():
    assert simpledistance([0]*10, [1]*4 + [0]*6) == 4
