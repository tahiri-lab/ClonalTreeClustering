# ClonalTreeClustering

## Overview
ClonalTreeClustering is part of a project to cluster B cell lineage trees. The objective of ClonalTreeClustering is to read and store B cell line trees in a C++ data structure. It also aims to construct a matrix of distances between all pairs of nodes. This README provides instructions on how to compile and run the program.

### 1st version
Files with "V1" in their name correspond to the first version of ClonalTreeClustering. Its purpose was to read and store B cell line trees in a C++ data structure before constructing a matrix of distances between all pairs of nodes. This version is no longer in use, see version 2.

### 2nd version
The aim of the second version of ClonalTree clustering is to read B-cell line trees and construct the distance matrix between all the nodes in the tree, as well as building a dictionary containing the abundance of each node. This is the version used and compiled by the Makefile.

This version takes as input a Newick file containing a single B cell line tree, a fasta file corresponding to the tree and an integer corresponding to the number of sequences in the fasta file. The fasta file and the integer are needed to have the abundance information to build the map.

This is not the ideal way of proceeding, as the ultimate aim is to enter a single Newick file containing several lineage trees. We would then need to provide a fasta file and the number of sequences for each of the trees contained in the Newick file. Eventually, it would be better to have the abundance information directly in the sequence identifiers in the Newick file, as in the file below (abundance is indicated after the @):
```
((((((seq18@0.66:1)seq11@6.92:1)seq4@5.27:1,(seq12@0.66:1,((seq23@0.49:1)seq19@0.49:1)seq13@0.66:1)seq5@1.81:1,seq6@1.32:1,((seq20@0.66:1)seq14@4.78:1,(seq21@1.81:1)seq15@0.66:1)seq7@0.99:1,seq8@0.82:1,((seq22@0.99:1)seq16@0.49:1,seq17@7.41:1)seq9@0.66:1,seq10@0.49:1)seq3@:1)seq1@5.44:1,seq2@25.21:1)naive@31.3:1);
```

With version 2, we are still proceeding with the fasta file and the number of sequences, in order to develop the distance matrix construction method for a single case, and while waiting for the final data supplied by the collaborators, which will contain the abundance information directly in the sequence identifiers in the Newick file.

## Prerequisites
Before compiling and running ClonalTreeClustering, ensure that you have the following:

- C++ compiler (e.g., GCC)
- Make utility

## Compilation
To compile the program, follow these steps :

1. Open a terminal and navigate to the "src" directory of the ClonalTreeClustering repository.
2. Run the following command to compile the program:
   ```
   make
   ```

## Usage
Once the program is compiled, you can run it using the following steps:

1. Make sure you are still in the "src" directory of the ClonalTreeClustering repository.
2. Run the following command to execute the program:
   ```
   ./ClonalTreeClustering <newick_file> <fasta_file> <#sequences>
   ```

   Replace `<newick_file>` with the path to the Newick file to be processed. For example:
   ```
   ../data/simulated-data/casX/casX.nk
   ```

   Replace `<fasta_file>` with the path to the FASTA file containing the sequences. For example:
   ```
   ../data/simulated-data/casX/casX-N.fasta
   ```

   Note: Replace `X` in the file paths with the specific case number, and replace `N` with the number of sequences (`<#sequences>`).

## Clean Up
To clean up the compiled files, you can run the following command in the "src" directory:
```
make clean
```

This will remove the object files and the executable generated during the compilation process.
