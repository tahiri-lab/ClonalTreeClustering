# WMFD Calculation and Clustering

This repository contains Python scripts to calculate the Weighted Minimum Feature Difference (WMFD) between tree pairs based on various metrics and to perform clustering on the calculated distances using DBSCAN.

## Features

- **WMFD Calculation**: Computes the weighted distance between tree pairs based on normalized metrics (Branch Length, Weight, Degree, Height, and Hamming Distance).
- **Distance Matrix**: Generates a symmetric distance matrix for the tree pairs.
- **DBSCAN Clustering**: Uses DBSCAN for clustering based on the computed distance matrix.

## Requirements

- pandas
- numpy
- scikit-learn

## Usage

1. **Input CSV**: The script reads a CSV file containing the following columns:
    - `Tree_Pair`: A string representing the pair of trees.
    - `Normalized_*`: Normalized values for various tree metrics (e.g., `Normalized_Common_BL`, `Normalized_Uncommon_Weight`).
    - `Penalty`: A penalty value for each pair.
  
2. **Parameters**: 
    - Enter the lambda values for each metric.
    - Set the parameters for DBSCAN clustering: `epsilon` (ε) and `minPoints`.

3. **Output**:
    - A symmetric distance matrix.
    - DBSCAN clustering results saved in a text file.

## Example

Run the script, follow the prompts to enter the lambda values and DBSCAN parameters, and the results will be saved in the output file.

## License

This project is licensed under the MIT License.



# WMFD Calculation and Clustering

This repository provides Python scripts for calculating the **Weighted Minimum Feature Difference (WMFD)** between tree pairs based on multiple metrics, followed by clustering the results using **DBSCAN**. The goal is to analyze similarities between phylogenetic trees based on different tree features, calculate distances between them, and apply clustering algorithms to group similar tree pairs.

## Features

- **WMFD Calculation**: 
  - The script computes the WMFD value between two trees by considering several metrics such as:
    - **Branch Length**
    - **Weight**
    - **Degree**
    - **Height**
    - **Hamming Distance**
  - The final WMFD is calculated based on these metrics and a penalty factor applied to the uncommon part.
  
- **Distance Matrix**:
  - The program creates a symmetric distance matrix from all tree pair comparisons.
  - Each entry represents the computed distance (WMFD) between two trees.

- **DBSCAN Clustering**:
  - **DBSCAN (Density-Based Spatial Clustering of Applications with Noise)** is used to cluster the trees based on the distance matrix.
  - Clusters are formed by evaluating tree distances, and any tree pair that exceeds a given threshold (epsilon) is clustered together.
  
- **User Inputs**:
  - The user can define the lambda values for each metric and specify DBSCAN parameters like `epsilon` (ε) and `minPoints` to control the clustering behavior.

- **Output**:
  - A text file containing:
    - The **distance matrix** between trees.
    - **Clustering results**, listing trees in each cluster.
  - This output provides insight into the relationships between the trees and groups them based on calculated distances.

## Requirements

Make sure to install the required Python libraries:

```bash
pip install pandas numpy scikit-learn
```

## Input Format

The script reads an input CSV file containing the following columns:

- `Tree_Pair`: A string representing the tree pair (e.g., `("Tree_1", "Tree_2")`).
- `Normalized_Common_*`: The normalized values for metrics (e.g., `Normalized_Common_BL` for Branch Length).
- `Normalized_Uncommon_*`: The normalized values for metrics in the uncommon part.
- `Penalty`: A numerical value indicating the penalty applied to the uncommon part in the calculation.

Example CSV:

| Tree_Pair           | Normalized_Common_BL | Normalized_Common_Weight | Normalized_Uncommon_BL | Normalized_Uncommon_Weight | Penalty | ... |
|---------------------|----------------------|--------------------------|------------------------|----------------------------|---------|-----|
| ("Tree_1", "Tree_2")| 0.75                 | 0.30                     | 0.50                   | 0.20                       | 1.5     | ... |
| ("Tree_1", "Tree_3")| 0.85                 | 0.40                     | 0.60                   | 0.35                       | 1.2     | ... |

## How to Run

1. Clone the repository:

```bash
git clone https://github.com/yourusername/WMFD-Clustering.git
cd WMFD-Clustering
```

2. Ensure your input CSV is ready (or modify the `input_path` in the script to point to your file).

3. Run the script:

```bash
python wmfd_clustering.py
```

4. Follow the prompts to enter the **lambda values** for each metric and the **DBSCAN parameters** (`epsilon` and `minPoints`).

5. Check the generated **output text file** for the distance matrix and clustering results.

## Output

- **Distance Matrix**: A symmetric matrix of WMFD values between trees.
- **Clustering Results**: A listing of clusters formed using DBSCAN, with each cluster's tree members.
- **File Format**: Results will be saved in a `.txt` file with the following sections:
  - Parameters used (lambda values, DBSCAN parameters)
  - Distance matrix
  - Clustering results

## Example Results

The output will include:

```
Parameters:
λ₁ (Branch Length) = 1.0
λ₂ (Weight) = 0.5
λ₃ (Degree) = 0.3
λ₄ (Height) = 0.2
λ₅ (Hamming Distance) = 0.7
epsilon (ε) = 0.1
minPoints = 3

Distance Matrix:
Tree    Tree_1  Tree_2  Tree_3
Tree_1    0.0     0.5     0.3
Tree_2    0.5     0.0     0.7
Tree_3    0.3     0.7     0.0

Clustering Results:
Cluster 0: Tree_1, Tree_2
Cluster 1: Tree_3
```

## License

This project is licensed under the MIT License.
