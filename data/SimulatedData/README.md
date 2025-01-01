# simulated-data

## B-cell lineage trees inferred by GCtree

In this folder, we used folder 50 of the "simulated-data" file.

folder "simulated-data" contains several files with a total of 92 simulated BCR lines obtained with the GCtree [[WS DeWitt III et al., 2018]](https://academic.oup.com/mbe/article/35/5/1253/4893244) simulator available here: https://github.com/matsengrp/gctree

Folder 50 contains several fasta files corresponding to the BCR lineage tree, containing at most the same number of sequences as the folder name.

For each given lineage file "tree.fasta", there is:
- the real tree "tree.GT.naive.nk"
- the tree deduced by GCtree "tree.GT.nk"

The GBLD scores of these 10 lineage trees are saved in the output file.


  
# Simulated Data: B-cell Lineage Trees Inferred by GCtree 🧬🌳

## Overview 🔍

This folder contains data from **Folder 50** of the "simulated-data" collection. The folder includes a total of **92 simulated B-cell receptor (BCR) lineages**, generated using the **GCtree simulator** [[WS DeWitt III et al., 2018]](https://academic.oup.com/mbe/article/35/5/1253/4893244). The simulator is available for download here: [GCtree GitHub Repository](https://github.com/matsengrp/gctree).

## Folder Structure 📁

- **Simulated-Data**: Main directory containing all simulated BCR lineage data.
- **Folder 50**: Contains several **FASTA files** for BCR lineage trees. The number of sequences in each file corresponds to the folder name (50 sequences in this case).

For each lineage file in **Folder 50**, you will find:
1. **tree.fasta**: FASTA file containing BCR lineage sequences.
2. **tree.GT.naive.nk**: The real, **ground truth tree** in Newick format.
3. **tree.GT.nk**: The tree **deduced by GCtree** in Newick format.

## Output 📝

For each lineage tree, the **GBLD (Generalized Branch Length Distance)** scores are calculated and saved in an output file. This includes the **GBLD scores** for the **10 lineage trees** in Folder 50, as well as the **intermediate steps** leading to the final score.

## Citation 📚

If you use this data or methodology in your research, please cite the following article:

> DeWitt III, W.S., et al. (2018). "GCtree: A simulator for B-cell lineage tree construction." *Molecular Biology and Evolution*, 35(5), 1253–1264. [Link to paper](https://academic.oup.com/mbe/article/35/5/1253/4893244).
