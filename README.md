# Jukes-Cantor Distance and Neighbor-Joining Phylogenetic Tree

## Overview
This Python module is designed to calculate the Jukes-Cantor genetic distance between sequences and to construct a phylogenetic tree using the neighbor-joining method. The Jukes-Cantor model is a formula used to estimate the number of nucleotide substitutions that have occurred in the DNA since two species have diverged from a common ancestor. The neighbor-joining method is a bottom-up clustering approach for the creation of phylogenetic trees, usually based on genetic distance.

https://en.wikipedia.org/wiki/Substitution_model

## Requirements
- Python 3.x
- NumPy
- SciPy

These dependencies can be installed using pip:
```bash
pip install numpy scipy
```

## Usage 

The module contains two main functions:

jukes_cantor_distance(seq1, seq2): Calculates the Jukes-Cantor distance between two DNA sequences.
neighbor_joining(dist_matrix): Performs the neighbor-joining algorithm on a distance matrix and returns the resulting tree.
To use these functions, sequences must be provided in a directory with FASTA files. The script reads these sequences, calculates the pairwise distances, and then constructs the phylogenetic tree.

## Steps 

1. Place your FASTA files in a directory.
2. Update the directory variable in the script with the path to your FASTA files. (directory = '/path/to/your/fasta/files')
3. Run the script.

## Output 

The output of the script is a phylogenetic tree printed in the console in the format of a dictionary where each key represents a node in the tree, and the value is a list containing the connected nodes and their respective distances.

## Notables 

- Ensure that all the sequences are of the same length and properly aligned for accurate distance calculations.
- The distance matrix is computed using the Jukes-Cantor model, and it assumes that all base substitutions are equally likely.

## Author 
Harrison Heath 
