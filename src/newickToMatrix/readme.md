# n2m_Lineage README
Author : Anna ARTIGES.

### Description :
This program computes the matrix of distance between lineage trees. The distance is calculated from the GLBD metric.

### Requirements :
 - Linux
 - c++ (any c++ compiler)

## How to create executable file :

With Linux/UNIX, you need to compile the n2m_Lineage.cpp file, to do it carry the command :
 > c++ n2m_Lineage.cpp -o n2m

## How to run the n2m program :
The program works if you have one file containing multiple Newick sequences with a name (name of the tree they represent). If a sequence does not have a name, it will be assigned a name by the program. The program also works with multiple Newick file as input, each file must contain only one Newick sequence.

### ⚠ Nodes cannot have the word "node" in their names or they won't be compared between trees. Be sure to check the Newick sequences before running the program.

#### Run the following :
Single file with X sequences in it :
> ./n2m -s X FILE output_file.txt

Multiple files :
> ./n2m -m FILE_1 FILE_2 output_file.txt

## Output file

The output_file.txt contains the distance matrix of each lineage tree and the matrix of distances between trees.
Before each distance matrix the name of the tree and number of nodes is written as follows :
``` tree_1 	 number of nodes :12 ```

At the end of the file is the comparison matrix between all trees. A value of -1.0 means that the trees have less than 3 nodes in common and are therefore not compared.
```
Comparison Matrix
tree_1	 0.000000	  13.228757	 7.681146  -1.0
tree_2	 13.228757	 0.000000	  14.966630 -1.0
tree_3	 7.681146	  14.966630	 0.000000  -1.0
tree_4  -1.0       -1.0       -1.0      0.000000
```



