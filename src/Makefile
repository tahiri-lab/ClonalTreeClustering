CXX = g++
CXXFLAGS = -std=c++11

TARGET = ClonalTreeClustering

$(TARGET): ClonalTreeClustering.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
	@echo ""
	@echo "//=================================================================================="
	@echo "Program   : ClonalTreeClustering - 2023"
	@echo "Authors   : Amel Chadli and Nadia Tahiri (Universite de Sherbrooke)"
	@echo "> ClonalTreeClustering is part of a project to cluster B cell lineage trees."
	@echo "> The objective of ClonalTreeClustering is to read and store B cell line trees in a"
	@echo "  C++ data structure."
	@echo "> It also aims to construct a matrix of distances between all pairs of nodes."
	@echo "//=================================================================================="
	@echo "==> Run the following command to execute the program:"
	@echo "./ClonalTreeClustering <newick_file> <fasta_file> <#sequences>"
	@echo ""
	@echo "==> Command line execution example:"
	@echo "./ClonalTreeClustering ../data/simulated-data/cas1/cas1.nk ../data/simulated-data/cas1/cas1-14.fasta 14"

.PHONY: clean
clean:
	rm -f $(TARGET)
