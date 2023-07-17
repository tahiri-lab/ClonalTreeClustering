# ClonalTreeClustering

## Overview
ClonalTreeClustering is part of a project to cluster B cell lineage trees. The objective of ClonalTreeClustering is to read and store B cell line trees in a C++ data structure. It also aims to construct a matrix of distances between all pairs of nodes. This README provides instructions on how to compile and run the program.

### 1st version
Files with "V1" in their name correspond to the first version of ClonalTreeClustering. Its purpose was to read and store B cell line trees in a C++ data structure before constructing a matrix of distances between all pairs of nodes. This version is no longer in use, see version 2.

### 2nd version
The aim of the second version of ClonalTree clustering is to read B-cell line trees and construct the distance matrix between all the nodes in the tree, as well as building a dictionary containing the abundance of each node. This is the version used and compiled by the Makefile.



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
