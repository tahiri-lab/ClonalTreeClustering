One small set of the newicks of 3 lineage trees is provided to show the steps to calculate the GBLD score between two lineage trees.

The output file contains the Normalized weights, differences of the normalized weights, averaged weight score based on the maximum number of nodes between two lineage trees under comparison, normalized branch length distances, penality, and finally GBLD score between each pair of lineage trees.




# GBLD Score Calculation for Lineage Trees 🌳🔬

## Overview 📚

This repository provides a small set of Newick formatted lineage trees to demonstrate the steps involved in calculating the **GBLD** (Generalized Branch Length Distance) score between two lineage trees. The GBLD score is a metric used to evaluate the similarity between lineage trees based on branch lengths and node structure.

## Input 🔄

A set of **Newick formatted** lineage trees is included. These trees serve as input for the GBLD score calculation.

## Output 📊

The output file contains the following information:

- **Normalized Weights** 📏: Normalized values of the branch lengths for each tree.
- **Differences of the Normalized Weights** ⚖️: The differences between the normalized weights of the two trees.
- **Averaged Weight Score** 🧮: The averaged weight score, calculated based on the maximum number of nodes between the two trees under comparison.
- **Normalized Branch Length Distances** 🌿: The distances between each pair of nodes of the trees, normalized.
- **Penalty** ❌: A penalty value applied to account for uncommon nodes.
- **GBLD Score** 🔢: The final **Generalized Branch Length Distance** score between each pair of lineage trees.

## Usage 🚀

To calculate the GBLD score between two lineage trees, follow these steps:

1. Prepare the input **Newick** trees.
2. Run the script to compute the GBLD score.
3. Retrieve the output containing the normalized weights, differences, averaged weight score, normalized branch length distances, penalty, and the final GBLD score.

## Example 🖥️

```bash
python calculate_gbld.py input_tree1.nwk input_tree2.nwk
```

## Dependencies 📦

- Python 3.x
- Necessary libraries (e.g., `BioPython`, `NumPy`)
