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
