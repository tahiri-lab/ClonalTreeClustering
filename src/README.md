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
