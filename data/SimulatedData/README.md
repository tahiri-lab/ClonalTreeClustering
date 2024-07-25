# simulated-data

## B-cell lineage trees inferred by GCtree

In this folder, we used folder 50 of the "simulated-data" file.

folder "simulated-data" contains several files with a total of 92 simulated BCR lines obtained with the GCtree [[WS DeWitt III et al., 2018]](https://academic.oup.com/mbe/article/35/5/1253/4893244) simulator available here: https://github.com/matsengrp/gctree

Folder 50 contains several fasta files corresponding to the BCR lineage tree, containing at most the same number of sequences as the folder name.

For each given lineage file "tree.fasta", there is:
- the real tree "tree.GT.naive.nk"
- the tree deduced by GCtree "tree.GT.nk"

The GBLD scores of these 10 lineage trees are saved in the output file.

