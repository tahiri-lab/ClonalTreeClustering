CXX = g++
CXXFLAGS = -std=c++11

TARGET = ClonalTreeClustering

$(TARGET): ClonalTreeClustering.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
	@echo ""
	@echo "//==========================================================================================="
	@echo "Program   : ClonalTreeClustering - 2023"
	@echo "Authors   : Amel Chadli and Nadia Tahiri (Universite de Sherbrooke)"
	@echo "> ClonalTreeClustering is part of a project to cluster B cell lineage trees."
	@echo "> The aim of ClonalTreeClustering is to read a Newick file containing lineage trees and"
	@echo "  construct a distance matrix between all pairs of nodes for each of the trees contained"
	@echo "  in the Newick file provided as input to the program."
	@echo "> The program also builds a map containing abundance information for all nodes in each tree."
	@echo "//==========================================================================================="
	@echo "==> Run the following command to execute the program:"
	@echo "./ClonalTreeClustering <newick_file> <fasta_file> <#sequences>"
	@echo ""
	@echo "==> Command line execution example:"
	@echo "./ClonalTreeClustering ../data/simulated-data/cas1/cas1.nk ../data/simulated-data/cas1/cas1-14.fasta 14"

.PHONY: clean
clean:
	rm -f $(TARGET)
