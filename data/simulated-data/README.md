# simulated-data

## B-cell lineage trees infered by GCtree

This folder contains several files with a total of 92 simulated BCR lines obtained with the GCtree [[WS DeWitt III et al., 2018]](https://academic.oup.com/mbe/article/35/5/1253/4893244) simulator available here: https://github.com/matsengrp/gctree

Each folder contains several fasta files corresponding to the BCR lineage tree, containing at most the same number of sequences as the folder name.

For each given lineage file "tree.fasta", there is:
- the real tree "tree.GT.naive.nk"
- the tree deduced by GCtree "tree.GT.nk"

## Simple cases (case 1, case 2 and case 3)

For each given case, created manually, there is:
- the alignment file "casX-N.fasta", where X is the number of the case and N the number of sequences
- the tree "casX.nk"
