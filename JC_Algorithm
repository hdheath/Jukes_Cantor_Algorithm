import os
import math
from scipy.cluster.hierarchy import dendrogram, to_tree
from scipy.spatial.distance import squareform
import numpy as np

def jukes_cantor_distance(seq1, seq2):
    """
    Calculates the Jukes-Cantor distance between two DNA sequences.
    """
    n_diff = sum(1 for a, b in zip(seq1, seq2) if a != b)
    p = float(n_diff) / len(seq1)
    d = -3/4 * math.log(1 - 4/3 * p)
    return d

def neighbor_joining(dist_matrix):
    """
    Performs the neighbor joining algorithm on a distance matrix and returns the resulting tree.
    """
    n = len(dist_matrix)
    D = np.copy(dist_matrix)
    tree = {}
    for i in range(n):
        tree[i] = []
    node_names = list(range(n))
    for _ in range(n - 2):
        r = np.sum(D, axis=1)
        q = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                q[i][j] = (n - 2) * D[i][j] - r[i] - r[j]
        min_q = np.inf
        min_i = None
        min_j = None
        for i in range(n):
            for j in range(i + 1, n):
                if q[i][j] < min_q:
                    min_q = q[i][j]
                    min_i = i
                    min_j = j
        delta = (r[min_i] - r[min_j]) / (n - 2)
        limb_i = (D[min_i][min_j] + delta) / 2
        limb_j = (D[min_i][min_j] - delta) / 2
        node_names.append(node_names[-1] + 1)
        tree[node_names[-1]] = [node_names[min_i], limb_i - limb_j]
        tree[node_names[min_i]].append(limb_j)
        tree[node_names[-1]].append(node_names[min_j])
        tree[node_names[min_j]].append(limb_j)
        D_new = np.zeros((n - 1, n - 1))
        node_indices = [i for i in range(n) if i != min_i and i != min_j]
        D_new[:-1, :-1] = D[np.ix_(node_indices, node_indices)]
        D_new[-1, :-1] = (D[min_i][node_indices] + D[min_j][node_indices] - D[min_i][min_j]) / 2
        D_new[:-1, -1] = D_new[-1, :-1]
        D = D_new
        n = len(D)
    tree[node_names[0]] = [node_names[1], D[0][1] / 2]
    tree[node_names[1]] = [D[0][1] / 2]
    return tree

directory = '/Users/HMans_MacBook_Pro/Desktop/ECSLAB4'
sequences = []

# Read in sequences from fasta files
for filename in os.listdir(directory):
    with open(os.path.join(directory, filename), 'r', encoding='ISO-8859-1') as f:
        lines = [line.strip() for line in f.readlines()]
        sequence = ''.join(lines[1:])
        sequences.append(sequence)

# Calculate Jukes-Cantor distance between all pairs of sequences
distances = []
for i, seq1 in enumerate(sequences):
    for j, seq2 in enumerate(sequences):
        if j <= i:
            continue
        if len(seq1) == 0 or len(seq2) == 0:
            continue
        distance = jukes_cantor_distance(seq1, seq2)
        distances.append(distance)
        
# Convert distances to a distance matrix
dist_matrix = squareform(distances)

tree = neighbor_joining(dist_matrix)
print(tree)
