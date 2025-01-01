# A Public Database of Memory and Naive B-Cell Receptor Sequences

One public database of more than 37 million unique BCR sequences from three healthy adult donors.
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0160853

-150 nucleotides of this data set are sampled randomly and converted to the weighted Newick format to use in GBLD metric.

The output file contains GBLD scores along with the process until reaching these scores.



# A Public Database of Memory and Naive B-Cell Receptor Sequences

This project leverages a public database containing over **37 million unique B-cell receptor (BCR) sequences** from three healthy adult donors. The data can be accessed via the following publication:

[PLoS ONE: A Public Database of Memory and Naive B-Cell Receptor Sequences](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0160853)

## Dataset Description 📊

The dataset consists of **B-cell receptor sequences** from memory and naive B-cells. For the purposes of this project, **150 nucleotides** from this dataset are randomly sampled and converted into a **weighted Newick format**. This format is used for further analysis, particularly for calculating the **GBLD (Genetic Branch Length Distance)** metric.

## Process Overview 🔄

1. **Data Sampling**: 150 nucleotides are randomly selected from the available BCR sequences.
2. **Newick Conversion**: The sampled data is converted into a **weighted Newick format**.
3. **GBLD Calculation**: The Genetic Branch Length Distance (GBLD) score is calculated using the converted Newick sequences.

## Output 🖥️

The output file generated from this process contains:

- **GBLD Scores**: The final calculated GBLD score for each pair of BCR sequences.
- **Intermediate Process**: The steps leading to the final GBLD score, including the normalization of branch lengths, differences, penalties, and other intermediary calculations.

## Usage 🔧

To replicate the process and compute the GBLD scores for the BCR sequences:

1. Download the dataset or prepare your own set of BCR sequences in Newick format.
2. Follow the instructions in the repository to convert sequences and calculate the GBLD score.

## Citation 📚

If you use this dataset or methodology in your work, please cite the following article:

> [A Public Database of Memory and Naive B-Cell Receptor Sequences](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0160853)  
> *PLoS ONE, 2016*
