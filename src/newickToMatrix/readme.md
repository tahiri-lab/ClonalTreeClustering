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


Here’s the revised, professional version of your `README.md` with stickers:

```markdown
# n2m_Lineage README
**Author**: Anna ARTIGES

### 📜 Description:
The `n2m_Lineage` program computes the distance matrix between lineage trees using the **GLBD metric**. It processes Newick sequences, calculating pairwise distances to compare lineage trees effectively.

### ⚙️ Requirements:
- 🐧 **Linux** operating system
- 🖥️ **C++ compiler** (any compatible C++ compiler)

## 🔨 How to Create the Executable File:

To compile the `n2m_Lineage.cpp` file on Linux/UNIX, run the following command in the terminal:
```bash
c++ n2m_Lineage.cpp -o n2m
```

## 🏃‍♂️ How to Run the `n2m` Program:

The program can process files containing multiple Newick sequences, where each sequence represents a lineage tree. If a sequence does not have a name, the program will automatically assign one.

You can run the program with a single file containing multiple Newick sequences, or with multiple files, each containing a single Newick sequence.

### ⚠️ Important Notes:
- **Nodes cannot contain the word "node"** in their names, or they will not be compared between trees. Ensure that your Newick sequences are correctly formatted before running the program.

### 📝 Command Syntax:

#### Single File with X Sequences:
```bash
./n2m -s X FILE output_file.txt
```

#### Multiple Files:
```bash
./n2m -m FILE_1 FILE_2 output_file.txt
```

## 📂 Output File:

The `output_file.txt` will contain:
1. The distance matrix for each lineage tree.
2. The matrix of distances between trees.

### Example Format for Distance Matrix:
Before each distance matrix, the name of the tree and the number of nodes will be listed, for example:
```
tree_1   number of nodes: 12
```

At the end of the file, you'll find the comparison matrix between all trees. A value of `-1.0` indicates that the trees have fewer than 3 nodes in common and were not compared.

### Example Comparison Matrix:
```
Comparison Matrix
tree_1    0.000000    13.228757   7.681146   -1.0
tree_2    13.228757   0.000000    14.966630  -1.0
tree_3    7.681146    14.966630   0.000000   -1.0
tree_4    -1.0        -1.0        -1.0        0.000000
```
