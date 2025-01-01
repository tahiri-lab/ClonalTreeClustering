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


```markdown
# NewickToMatrix Conversion v0.1

**Authors**: UQAM Lab  
**Date**: June 2009

### 📜 Description:
The `NewickToMatrix` program converts lineage trees represented in Newick format into distance matrices. The program computes tree distances based on the **GLBD metric**. This utility is ideal for comparing lineage trees and analyzing the relationships between them in a structured manner.

### 🧰 Requirements:
- 🐧**Linux** or UNIX-based system
- 🖥️**C++ Compiler** (any standard C++ compiler)

### 🔨 How to Compile the Program:

To compile the `NewickToMatrix` program, use the following command on a Linux/UNIX system:

```bash
c++ NewickToMatrix.cpp -o n2m
```

This will create the executable `n2m` to run the program.

### 🏃‍♂️ How to Run the Program:

The program can process multiple input files containing Newick sequences or a single file with multiple sequences. The format of the input is determined by the command-line argument passed.

#### Command Syntax:

##### 1. Single File with Multiple Newick Sequences:
```bash
./n2m -s X input_file output_file.txt
```
Where `X` is the number of sequences in `input_file`, and `output_file.txt` is the destination for the resulting distance matrix.

##### 2. Multiple Files with One Newick Sequence Each:
```bash
./n2m -m file1 file2 output_file.txt
```
This command allows you to provide multiple files, each containing one Newick sequence.

### ⚠️ Important Notes:
- **Node names cannot contain the word "node"** as this will prevent comparisons between trees.
- Ensure that your Newick sequences are correctly formatted before running the program.

### 📝 Output File:
The `output_file.txt` will contain:
1. The distance matrix for each lineage tree.
2. A comparison matrix between all trees.

The format of the output will be as follows:

#### Example Output for Individual Distance Matrices:
```
tree_1    number of nodes: 12
```

#### Comparison Matrix:
The comparison matrix lists the distance between each tree. A value of `-1.0` indicates that the trees have fewer than 3 nodes in common and were not compared.

Example:
```
Comparison Matrix
tree_1    0.000000    13.228757    7.681146    -1.0
tree_2    13.228757   0.000000     14.966630   -1.0
tree_3    7.681146    14.966630    0.000000    -1.0
tree_4    -1.0        -1.0         -1.0         0.000000
```

### 🛠️ How it Works:
The program first reads Newick sequences, processes the lineage trees, and then computes the pairwise distances between them using the GLBD metric. It stores the resulting distances in a matrix and outputs this information into the specified output file.

1. **Input Handling**: The program can handle both single and multiple files containing Newick sequences.
2. **Distance Calculation**: For each pair of trees, the distance is computed and stored.
3. **Output**: The results are saved in a matrix format, making it easy to analyze the distances and relationships between the trees.

### 🔄 Future Improvements:
- Add more distance metrics for comparison.
- Improve error handling for malformed Newick sequences.

### 👨‍💻 Contributors:
- **UQAM Lab** for developing and maintaining the code.
```

