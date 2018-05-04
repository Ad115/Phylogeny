from phylogeny.clocklike_reconstruction import simpledistance

def test_simpledistance():
    assert simpledistance([0]*10, [1]*4 + [0]*6) == 4