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

#### Run the following :
Single file with X sequences in it :
> ./n2m -s X FILE output_file.txt

Multiple files :
> ./n2m -m FILE_1 FILE_2 output_file.txt

The output_file.txt contains the distance matrix of each lineage tree and the matrix of distances between trees.



