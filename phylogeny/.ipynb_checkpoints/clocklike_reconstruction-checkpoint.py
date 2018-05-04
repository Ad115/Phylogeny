def simpledistance(seq_1, seq_2):
    "From two binary sequences, compute their distance."
    differences = (1 for a,b in zip(seq_1, seq_2) if a != b)
    return sum(differences)
# ---