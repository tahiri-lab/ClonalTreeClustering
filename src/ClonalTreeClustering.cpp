#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <cctype>
#include <string>
#include <map>
#include <sstream>

#include "ClonalTreeClustering.h"
#include "struct-lineage.h"

using namespace std;

//================================================================
// FUNCTIONS
//================================================================

vector<vector<char>> readNewick(ifstream& file1) {
  vector<vector<char>> newickLineageWhole = {};
  vector<char> newickLineageLine = {};
  char c;
  
  while(file1.get(c)) {
    newickLineageLine.push_back(c);
    if(c == ';'){
      newickLineageWhole.push_back(newickLineageLine);
      newickLineageLine.clear();
    }
  }
  
  if (!newickLineageLine.empty()) {
      newickLineageWhole.push_back(newickLineageLine);
  }

  cout<<"---------------------------------------------------"<<endl;
  cout<<"lineage file :"<<endl;
  for(const vector<char>& line : newickLineageWhole) {
    for(char c : line) {
      cout<<c;
    }
    cout<<endl;
  }
  cout<<"---------------------------------------------------"<<endl;
  
  return newickLineageWhole;
}

void initInputLineage(struct InputLineage *aLineage) {
  aLineage->idSeq.push_back(vector<char>());
  //aLineage->root = 0;
  aLineage->size = 0;
  aLineage->abundance = {};
  aLineage->degree = {};
  aLineage->depth = {};
  aLineage->edge = {};
  aLineage->length.push_back(vector<double>());
  aLineage->ADJACENCE = NULL;  
}

//nbNodes -> uses "newickLineage", this function will be used to calculate de penalty

void readNewickLineage(vector<char> newickLineageLine, struct InputLineage *aLineage) {
  int cpt = 0;
  int s = 0; // the number of sequences = size of the lineage
  char c;
  vector<char> nameSeq;
  vector<double> length;
  regex reg("[seqnaiv]");
  char previousChar = '\0';
  char previousChar2;

  for(int i=0; newickLineageLine[i]!='\0'; i++) {
    c = newickLineageLine[i];
    if(c == '(')
      cpt++;
    else if(c == ')')
      cpt--;
    if((regex_search(string(1, c), reg))||(isdigit(c) && nameSeq.back() == 'q')||(isdigit(c) && previousChar2 == 'q')) {
      nameSeq.push_back(c);
    }
    if(c == ':') {
      aLineage->idSeq.push_back(nameSeq);
      nameSeq.clear();
      s++;
    }
    if((isdigit(c) && previousChar == ':')||(isdigit(c) && previousChar2 == ':')) {
      double doubleValue = stod(string(1, c));
      length.push_back(doubleValue);
    }
    if((c == ',')||(c == ')')) {
      aLineage->length.push_back(length);
      length.clear();
    }

    previousChar2 = previousChar;
    previousChar = c;
  }
  aLineage->size = s;

  if(cpt==0){
    // display cpt to ensure it is equal to 0
    cout<<"0. Number of parentheses = "<<cpt<<" (it's OK!)"<<endl;
  }else {
    cerr<<"Error: There's a problem with the number of parentheses."<<endl;
    //return -1;
  }

  // display the vector of vectors idSeq
  cout<<"1. idSeq :";
  for(const vector<char>& seq : aLineage->idSeq) {
    for(char c : seq) {
      cout<<c;
    }
    cout<<endl;
  }
  
  cout<<"2. Size of the lineage (nb of nodes) = "<<aLineage->size<<endl;
  
  // display the vector of vectors length
  cout<<"5. length :";
  for(const vector<double>& lgth : aLineage->length) {
    for(double l : lgth) {
      cout<<l;
    }
    cout<<endl;
  }
  cout<<"---------------------------------------------------"<<endl;
}

map<string, int> abundance(ifstream& file2, int seq) {
  map<string, int> dictionary;
  string line;
  
  while(getline(file2, line)) { // we go through the file line by line
      if(line.find(">naive") != string::npos) // count the occurrences of ">naive"
          dictionary["naive"]++;
      else { // we count the occurrence for each sequence from ">seq1"
          for(int i=1; i<=(seq+1); i++) {
              string key = "seq" + to_string(i);
              if(line.find(key) != string::npos)
                  dictionary[key]++;
          }
      }
  }

  // case of the tens adding to the units
  for(int i=1; i<=static_cast<int>(seq/10); i++) {
      string keyUnit = "seq" + to_string(i);
      for(int j=i*10; j<i*10+10; j++) {
          string key10 = "seq" + to_string(j);
          dictionary[keyUnit] -= dictionary[key10];
      }
  }
  
  auto it = dictionary.begin();
  while(it != dictionary.end()) {
    if(it->second == 0) {
      it = dictionary.erase(it);
    }else {
      it++;
    }
  }
  
  return dictionary;
}

void abundanceToInputLineage(struct InputLineage *aLineage, map<string, int> abundanceMap) {
  vector<string> convertedData;
  for(const auto& vec : aLineage->idSeq) {
    string str(vec.begin(), vec.end());
    convertedData.push_back(str);
  }

  for(const auto& str : convertedData) {
    for(const auto& pair : abundanceMap) {
      if(pair.first == str)
        aLineage->abundance.push_back(pair.second);
    }
  }

  cout<<"---------------------------------------------------"<<endl;
  cout<<"3. abundance :"<<endl;
  for(const auto& integer : aLineage->abundance) {
    cout<<integer<<endl;
  }
}

void degree(vector<char> newickLineageLine, struct InputLineage * aLineage) {
  char c;
  int cptChildren = 0;
  vector<char> idSeqDegree;
  map<string, int> degrees;
  //regex reg("[seqnaiv]");
  //char previousChar = '\0';
  //char previousChar2;

  int i = 0;
  int j = 0;
  int k = 0;
  while(newickLineageLine[i]!='\0'){
    c = newickLineageLine[i];
    if(c == 's') {
      while(c != ')'){
        k = i;
        while(c != ':') {
          c = newickLineageLine[k];
          if( ((isdigit(c))||(isalpha(c))) && (c != ':') )
          idSeqDegree.push_back(c);
          k++;
        }
        i++;
        c = newickLineageLine[i];
        if(c == ':')
          cptChildren++;
      }
      j = i+1;
      while(c != ':') {
        c = newickLineageLine[j];
        if( ((isdigit(c))||(isalpha(c))) && (c != ':') )
        idSeqDegree.push_back(c);
        j++;
      }
      stringstream ss;
      for(char letter : idSeqDegree)
        ss<<static_cast<char>(letter);
      idSeqDegree.clear();
      string key = ss.str();
      degrees[key] = cptChildren;
      for(const auto& pair : degrees)
        cout<<"Clé : "<<pair.first<<", Valeur : "<<pair.second<<endl;
      /*cout<<"idSeqDegree = ";
      for(char letter : idSeqDegree)
        cout<<letter;
      cout<<endl;
      cout<<"Chaine de caracteres : "<<key<<endl;*/
      i = j;
    }
    i++;
  }

  cout<<"Nb de ':' avant la 1ere ')' : "<<cptChildren<<endl;
}


/*void AdjacenceMatrix( double **Adjacence, long int *ARETE, double *LONGUEUR,int size,int kt){
	
	int i,j;
	
	for(i=1;i<=2*size-2;i++) //(n+1)
		for(j=1;j<=2*size-2;j++){
			Adjacence[i][j] = Adjacence[j][i] = INFINI;
}
	for(i=1;i<=2*size-3-kt;i++){
		Adjacence[ARETE[2*i-2]][ARETE[2*i-1]] = LONGUEUR[i-1];//(LONGUEUR[i-1]>5*epsilon)?LONGUEUR[i-1]:5*epsilon;
		Adjacence[ARETE[2*i-1]][ARETE[2*i-2]] = LONGUEUR[i-1]; //(LONGUEUR[i-1]>5*epsilon)?LONGUEUR[i-1]:5*epsilon;
	}
}
*/

//================================================================
// PROGRAMME
//================================================================

int main(int argc, char* argv[]) {

  if(argc != 4) {
    cerr<<"Usage: ./ClonalTreeClustering <newick_file> <fasta_file> <#sequences>"<<endl;
    return 1;
  }
  
  string filename1 = argv[1] ;
  string filename2 = argv[2] ;
  int nbSeq = stoi(argv[3]);

  // open the newick file in read mode
  ifstream file1(filename1);
  // check that the file has been opened correctly
  if(!file1) {
      cerr<<"Error: Unable to open the file."<<endl;
      return -1;
  }
  vector<vector<char>> newickLineageWhole = readNewick(file1);
  file1.close();

  // we read the "newick file" line by line
  for(const vector<char>& newickLineageLine : newickLineageWhole) {
    
    struct InputLineage aLineage;
    
    initInputLineage(&aLineage);
    readNewickLineage(newickLineageLine, &aLineage);

    // collect the info of the occurrences (different files)
    cout<<"---------------------------------------------------"<<endl;

    // open the fasta file in read mode
    ifstream file2(filename2);
    if(!file2) {
      cerr<<"Error: Unable to open the file."<<endl;
      return -1;
    }
    map<string, int> abundanceDico = abundance(file2, nbSeq);
    file2.close();
    
    // display the counts
    int total_count=0;
    cout << "Occurrences (abundances) :" << endl;
    for (const auto& pair : abundanceDico) {
      cout<<pair.first<<" : "<<pair.second<<endl;
      total_count += pair.second;
    }
    cout<<"Total count of sequences : "<<total_count<<endl;
    // add the info of the occurrences to the structure InputLineage
    abundanceToInputLineage(&aLineage, abundanceDico);

    degree(newickLineageLine, &aLineage);

    cout<<"---------------------------------------------------"<<endl;
    cout<<"**************** END OF 1 LINEAGE *****************"<<endl;
    cout<<"---------------------------------------------------"<<endl;
  }

  return 0;
}
