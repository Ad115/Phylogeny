import numpy as np
import itertools as itr
from .fpc import four_point_condition

def simple_distance(seq_1, seq_2):
    "From two binary sequences, compute their distance."
    differences = (1 for a,b in zip(seq_1, seq_2) if a != b)
    return sum(differences)
# ---

class DistanceMatrix(np.ndarray):
    """Wrapper for the Numpy array class with methods proper of a 
    distance matrix.
    
    For documentation for the Numpy array, read the `Numpy documentation`_.
    
    .. _Numpy documentation:
       http://www.numpy.org/
    """
    
    def __new__(cls, data, names=None):
        """        
        Args:
            data (iterable): Existing array.
            names (int|str, optional): Names of the nodes.
        """
        matrix = np.asarray(data).view(cls)
        
        if not names:
            n = len(matrix) 
            names = range(n)
            
        matrix.names = tuple(names)
        matrix.idx = {name:i for i,name in enumerate(matrix.names)}
        
        return matrix
    # ---
    
    def __array_finalize__(self, obj):
        if obj is None: 
            # (we're in the middle of the __new__
            # constructor, and self.names, self.idx 
            # will be set when we return to
            # __new__)
            return
        self.names = None
        self.idx = None
    # ---
    
    @classmethod
    def zeros(cls, n, names=None):
        """Return a zeroed n by n matrix"""
        return cls(data=np.zeros((n,n)), 
                   names=names)
    # ---
    
    @classmethod
    def from_sequences(cls, sequences, distance_fn=simple_distance):
        "From the given sequences, compute pairwise edit distances."
        # Get all the pairs
        pairs = itr.combinations(sequences, 2)
        # Compute distances
        distances = cls.zeros(len(sequences), names=sequences.keys())
        for i,j in pairs:
            d_ij = distance_fn(sequences[i], sequences[j])
            distances.set((i,j), d_ij)
        return distances
    # ---
    
    def __repr__(self):
        return (  f"{super().__repr__()[:-1]}, names={self.names})" )
    # ---
    
    def remove(self, name):
        """Return a new matrix with the column and row with that name deleted."""
        i = self.idx[name]
        
        # Remove row
        m = np.delete(self, (i), axis=0)
        # Remove column
        m = np.delete(m, (i), axis=1)
        
        return DistanceMatrix(m, 
                              names=(n for n in self.names if n != name))
    # ---        
    
    def is_additive(self, tolerance=1e-2):
        """Is the distances matrix additive?

        Check the four point condition on each quartet of
        indices of the matrix.
        """
        n = len(self)
        quartets = list(itr.combinations(range(n), 4))

        return all(four_point_condition(self, q, tolerance) 
                       for q in quartets)
    # ---
    
    def distances_to(self, name):
        "Get all the distances to the named sequence."
        i = self.idx[name]
        return {n:dist for n,dist in zip(self.names,self[i]) if n != name}
    # ---
    
    def get(self, item):
        "Get item by name."
        i,j = item
        idx = self.idx
        return self[idx[i], idx[j]]
    # ---
    
    def set(self, item, value):
        "Set item by name."
        i,j = item
        idx = self.idx
        self[idx[i], idx[j]] = value
        self[idx[j], idx[i]] = value
    # ---
    
    def name_all(self):
        names = self.names
        n = len(names)
        for i in range(n):
            for j in range(i+1, n):
                yield (names[i],names[j]), self[i,j]
    # ---
        
# --- DistanceMatrix