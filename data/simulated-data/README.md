# 🌳 **simulated-data**

## 🧬 **B-cell Lineage Trees Inferred by GCtree**

This folder contains 92 simulated BCR lineage trees generated using the **GCtree** simulator by [WS DeWitt III et al., 2018](https://academic.oup.com/mbe/article/35/5/1253/4893244). You can access the GCtree simulator [here](https://github.com/matsengrp/gctree).

Each folder represents a lineage tree and contains several **FASTA** files with sequences. The number of sequences corresponds to the folder name. For each lineage file **tree.fasta**, you will find:
- **tree.GT.naive.nk**: The real tree.
- **tree.GT.nk**: The tree deduced by GCtree.

---

## 🧑‍🔬 **Simple Cases (Case 1, Case 2, Case 3)**

For each of the manually created cases, the following files are provided:
- **casX-N.fasta**: The alignment file (X is the case number and N is the number of sequences).
- **casX.nk**: The corresponding tree.
