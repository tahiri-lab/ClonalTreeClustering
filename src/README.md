# ClonalTreeClustering

## Overview
ClonalTreeClustering is part of a project to cluster B cell lineage trees. The aim of ClonalTreeClustering is to read a Newick file containing lineage trees and construct a distance matrix between all pairs of nodes for each of the trees contained in the Newick file provided as input to the program. The program also builds a map containing abundance information for all nodes in each tree.

This is then used to build clusters of B cell line trees using the [KMeansSuperTreeClustering](https://github.com/tahiri-lab/KMeansSuperTreeClustering) program.

This README provides instructions on how to compile and run the program.

### 1st version (V1_ClonalTreeClustering.cpp)
Files with "V1" in their name correspond to the first version of ClonalTreeClustering. Its purpose was to read and store B cell line trees in a C++ data structure before constructing a matrix of distances between all pairs of nodes. This version is no longer in use, see version 2.

### 2nd version (ClonalTreeClustering.cpp)
The aim of the second version of ClonalTree clustering is to read B-cell line trees and construct the distance matrix between all the nodes in the tree, as well as building a map containing the abundance of each node. This is the version used and compiled by the Makefile.

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





Here’s a polished and professional version of your README file with a more visually appealing structure and design elements included:

---

# 🧬 ClonalTreeClustering  

**ClonalTreeClustering** is a bioinformatics tool designed for clustering B-cell lineage trees. It processes Newick files containing lineage trees, constructs a distance matrix for all node pairs, and builds a map with abundance information for each node. This data is used to perform clustering of B-cell lineage trees using the [KMeansSuperTreeClustering](https://github.com/tahiri-lab/KMeansSuperTreeClustering) program.  

---

## 📖 Overview  
- **Purpose**: To construct distance matrices and abundance maps for nodes in lineage trees.  
- **Input**:  
  - A Newick file with lineage trees.  
  - A FASTA file corresponding to the tree.  
  - An integer specifying the number of sequences in the FASTA file.  
- **Output**:  
  - Distance matrices between nodes.  
  - Abundance maps for nodes in each tree.  

---

## 🛠️ Versions  

### **1st Version: V1_ClonalTreeClustering.cpp**  
- Reads and stores B-cell lineage trees in a C++ data structure.  
- Constructs a distance matrix for node pairs.  
- **Status**: Deprecated.  

### **2nd Version: ClonalTreeClustering.cpp**  
- Reads lineage trees and constructs distance matrices and abundance maps.  
- Accepts a single Newick file with one tree, a FASTA file, and the number of sequences as input.  
- **Future Improvement**: The aim is to directly include abundance information in sequence identifiers within the Newick file.  

**Example Newick Format with Abundance**:  
```
((((((seq18@0.66:1)seq11@6.92:1)seq4@5.27:1,...))naive@31.3:1);
```  

---

## ✅ Prerequisites  
Make sure you have the following tools installed before running **ClonalTreeClustering**:  
- **C++ Compiler** (e.g., GCC)  
- **Make Utility**  

---

## 🚀 How to Compile  

1. Navigate to the `src` directory of the repository:  
   ```bash
   cd src
   ```  

2. Compile the program using the `make` command:  
   ```bash
   make
   ```  

---

## 🖥️ Usage  

Once compiled, run the program with:  
```bash
./ClonalTreeClustering <newick_file> <fasta_file> <#sequences>
```  

### Example:  
```bash
./ClonalTreeClustering ../data/simulated-data/casX/casX.nk ../data/simulated-data/casX/casX-50.fasta 50
```  
- **`<newick_file>`**: Path to the Newick file (e.g., `casX.nk`).  
- **`<fasta_file>`**: Path to the FASTA file (e.g., `casX-50.fasta`).  
- **`<#sequences>`**: Number of sequences (e.g., `50`).  

---

## 🧹 Clean-Up  

Remove all compiled files and executables with:  
```bash
make clean
```  

---

## 📄 Notes and Future Plans  

- Currently, abundance information is provided separately via a FASTA file and sequence count.  
- Future improvements aim to streamline the process by including abundance information directly in sequence identifiers within the Newick file.  

---  

## 🎯 Example Directory Structure  

```
ClonalTreeClustering/  
│  
├── src/  
│   ├── ClonalTreeClustering.cpp  
│   ├── V1_ClonalTreeClustering.cpp  
│   ├── Makefile  
│  
├── data/  
│   ├── simulated-data/  
│       ├── casX/  
│           ├── casX.nk  
│           ├── casX-50.fasta  
│  
└── README.md  
```  

---  

## 🤝 Acknowledgments  
This project is developed as part of ongoing research in B-cell lineage tree clustering. Contributions from collaborators in the bioinformatics community are greatly appreciated.  

---

Feel free to customize further with tools like [Shields.io](https://shields.io/) for badges (e.g., build status, license) and additional markdown features for tables or icons. This enhanced README makes it more professional and user-friendly!
