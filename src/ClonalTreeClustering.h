#ifndef HEADER_H
#define HEADER_H

#include <fstream>
#include <map>
#include <vector>

using namespace std;

vector<vector<char>> readNewick(ifstream& filename); // we open and close the file only once

void initInputLineage(struct InputLineage *aLineage);

//int nbNodes(string newickLineage); // for the penalty

void readNewickLineage(vector<char> newickLineageLine, struct InputLineage * aLineage);

map<string, int> abundance(ifstream& filename2, int seq);

void abundanceToInputLineage(struct InputLineage *aLineage, map<string, int> abundanceMap);

void degree(vector<char> newickLineageLine, struct InputLineage * aLineage);
  
#endif
