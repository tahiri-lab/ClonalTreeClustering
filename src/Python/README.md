# 🌿 GBLD Program: Generalized Branch Length Distance Metric

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)  
[![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-green)](#)  
[![License](https://img.shields.io/badge/License-MIT-brightgreen)](#)

## Overview  
The **GBLD Program** implements the *Generalized Branch Length Distance (GBLD)* metric to compare lineage trees by considering overlapping nodes, tree structure, and branch weights. This tool processes lineage tree data (in Newick format) and computes key metrics like distances, penalties, and normalized differences to derive the GBLD score.

### Key Features  
- Parse **Newick** format to extract node weights.  
- Normalize branch weights and compute pairwise differences.  
- Calculate the **GBLD metric** components:
  - **W**: Weight differences between trees.
  - **D**: Structural distance between tree branches.
  - **P**: Penalty based on overlapping nodes.
- Support for large-scale lineage tree comparisons.  

---

## 📜 Workflow  
The program operates through several structured steps:

1. **Parse Newick**  
   Extracts node weights and abundance values from the Newick-formatted strings.  

2. **Weight Normalization**  
   Scales weights to a 0–1 range for uniform comparison.  

3. **Calculate Distances**  
   Computes pairwise distances between sequences in trees using the Newick structure.  

4. **Matrix Normalization**  
   Normalizes distance matrices to a 0–1 range.  

5. **Calculate GBLD Components**  
   - **W**: Pairwise weight differences.  
   - **D**: Pairwise branch distance differences.  
   - **P**: Penalty based on shared vs. total nodes.  

6. **Compute GBLD Metric**  
   Final GBLD values combine the components \( W + D + P \).  

---

## 🔧 Installation and Prerequisites  

### Prerequisites  
Ensure the following are installed on your system:  

- Python 3.8+  
- [Biopython](https://biopython.org/)  

### Installation  
Clone the repository and install dependencies:  
```bash
git clone <repository-url>
cd <repository-folder>
pip install -r requirements.txt
```

---

## 🚀 Usage  

### Input Requirements  
The program expects:  
1. **Input File**: Text file where each line contains a tree label and a Newick-formatted string.  
2. **Output Path**: File path to save results.  

### Run the Program  
Use the following command to run the program:  
```bash
python gbld_program.py
```

### Example Input  
A text file (`example.txt`) containing trees:  
```text
Tree1: (((seq1@2.0:1,seq2@3.0:2):3,seq3@1.5:1):1);
Tree2: (((seq1@2.0:1,seq2@3.0:2):3,seq4@1.0:2):1);
```

### Output  
Results are saved in a structured file containing:  
1. Original and normalized weight matrices.  
2. Pairwise GBLD components (\( W, D, P \)).  
3. Final GBLD scores.  

---

## 📂 Example Results  

### Input  
Example Newick-formatted trees:  
```text
Tree1: (((seq1@2.0:1,seq2@3.0:2):3,seq3@1.5:1):1);
Tree2: (((seq1@2.0:1,seq2@3.0:2):3,seq4@1.0:2):1);
```

### Output Snippet  
```
Original Weights:
Node, Tree1, Tree2
seq1, 2.0, 2.0
seq2, 3.0, 3.0
seq3, 1.5, 0.0
seq4, 0.0, 1.0

Normalized Weights:
Node, Tree1, Tree2, Diff Tree1-Tree2
seq1, 1.0, 1.0, 0.0
seq2, 1.0, 1.0, 0.0
seq3, 1.0, 0.0, 1.0
seq4, 0.0, 1.0, 1.0

GBLD(T1, T2):
W = 0.3333
D = 0.5000
P = 0.6667
Final GBLD = 1.5000
```

---

## 📚 Methodology  

### 1. Parsing Newick Format  
Extracts node weights, branch lengths, and abundance values.  

### 2. Normalization  
- **Weights**: Scaled to a 0–1 range for uniformity.  
- **Matrices**: Distance matrices normalized to eliminate bias from differing scales.  

### 3. Metrics Calculation  
- \( W \): Mean difference in normalized weights.  
- \( D \): Root-mean-square of pairwise distances.  
- \( P \): Overlapping node penalty.  

### 4. GBLD Computation  
The GBLD score is the sum of the three metrics:  
\[ \text{GBLD} = W + D + P \]  

---

## 📄 License  
This project is licensed under the MIT License. See `LICENSE` for details.  

---

## 🤝 Contributing  

Contributions, issues, and feature requests are welcome! Feel free to submit a pull request or raise an issue in the repository.  

---

## 💡 Acknowledgements  
This tool was developed to support bioinformatics research in lineage tree analysis and comparison. Special thanks to the Biopython community for providing powerful libraries for sequence analysis.

