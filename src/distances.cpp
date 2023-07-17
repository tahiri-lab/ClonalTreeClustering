#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <cctype>
#include <string>
#include <map>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>

using namespace std;

//================================================================
// FUNCTIONS
//================================================================

//= convert newick file to vector<vector<chat>>
//= each line corresponds to a lineage tree
vector<vector<char>> readNewick(ifstream& file1) {
  vector<vector<char>> newickLineageWhole = {};
  vector<char> newickLineageLine = {};
  string line;
  
  while(getline(file1, line)) {
    if(!line.empty()) {
      for(char c : line) {
        newickLineageLine.push_back(c);
        if(c == ';') {
          newickLineageWhole.push_back(newickLineageLine);
          newickLineageLine.clear();
        }
      }
    }
  }

  cout<<"---------------------------------------------------"<<endl;
  cout<<"lineage file: ";
  for(const vector<char>& line : newickLineageWhole) {
    for(char c : line) {
      cout<<c;
    }
    cout<<endl;
  }
  cout<<"---------------------------------------------------"<<endl;
  
  return newickLineageWhole;
}

//==================== RELATIONSHIPS BETWEEN NODES ====================//

typedef struct node_t{
  const char *name;
  const char *parentName;
  //node_t *parent;//add this for the relationship map
  int parent;//add this for the relationship vector
  double distance;//add this for the relationship vector
  double blen;//branch length
  node_t *up;
  int nbranches;
  int mbranches;
  node_t **branches;
  int serial;
}node_t;

int serial =0;

node_t *alloc(){

  node_t *nd = new node_t;
  nd->name =NULL;
  nd->blen=-1;
  nd->mbranches = 2;
  nd->nbranches = 0;
  nd->up =NULL;
  nd->branches =(node_t **) calloc(nd->mbranches,sizeof(node_t*));
  nd->serial = serial++;
  return nd;
}

vector<string> mysplit2(char *str){
  vector<string> vec;
  int i=0;
  char buf[strlen(str)];//THIS IS ENOUGH
  memset(buf,0,strlen(str));
  int at=0;
  for(;i<strlen(str);i++){
    if(str[i]=='\t'||str[i]==' '||str[i]==';')
      continue;
      
    if(str[i]=='('||str[i]==')'||str[i]==':'||str[i]==','){
      if(strlen(buf)>0){
	vec.push_back(buf);
	memset(buf,0,strlen(str));
	at=0;
      }
      vec.push_back(string()+str[i]);
      continue;
    }else{
      buf[at++]=str[i];
    }
  }
  if(strlen(buf)>0){
    vec.push_back(buf);
    memset(buf,0,strlen(str));
    at=0;
  }
  
  return vec;
}

int atter=0;//dag

//production rule: prod := (prod,prod)nam:len
node_t *parse(vector<string> &vec){

  node_t *nd = alloc();
  int seri = nd->serial;
  if(atter>=vec.size())
    return nd;
  //catch leaflike case
  if(vec[atter]!="("){
    if(vec[atter]!=":"){
      nd->name = strdup(vec[atter].c_str());
      atter++;
    }
    if(vec[atter]==":"){
      nd->blen = atof(vec[atter+1].c_str());
      atter+=2;
    }
    return nd;
  }

  //catch recursive cases;
  if(vec[atter]=="("){
    while(1){
      atter++;
      node_t *anode =parse(vec);
      anode->up = nd;
      if(nd->nbranches==nd->mbranches){
	nd->mbranches = 2*nd->mbranches;
	nd->branches =(node_t**) realloc(nd->branches,sizeof(node_t*)*nd->mbranches);
	fprintf(stderr,"reallocing\n");
      }
      nd->branches[nd->nbranches++] = anode;
      if(vec[atter]==",")
	continue;
      else
	break;
    }
  }
  if(vec[atter]==")"){
    atter++;
    if(vec[atter]==")")//catch the case when there is no name or branch length
      return nd;
  }

  if(vec[atter]!=":"){
    nd->name = strdup(vec[atter].c_str());
    atter++;
  }
  if(vec[atter]==":"){
    nd->blen = atof(vec[atter+1].c_str());
    atter+=2;
  }
  return nd;
}

void print_node(FILE *fp,node_t *nd){
  fprintf(fp,"----------\nnd->serial: %d this: %p\n",nd->serial,nd);
  //  fprintf(fp,"nd->left: %p nd->right: %p nd->up: %p\n",nd->left,nd->right,nd->up);
  fprintf(fp,"nd->name: %s nd->blen: %f nd->up:%p\n----------\n",nd->name,nd->blen,nd->up);
  
}

void serialize(node_t *nd,node_t **lst){
  assert(nd);
  lst[nd->serial] = nd;
  for(int i=0;i<nd->nbranches;i++)
    if(nd->branches[i])
      serialize(nd->branches[i],lst);
  
}

// add a function to store relationships between nodes in a dictionary
vector<node_t> createParentTable(node_t* root) {
  vector<node_t> parentTable;

  map<node_t*, string> nodeParentMap; //map to store the parent name for each node

  //tree prefix browsing to fill parent table
  vector<node_t*> stack;
  stack.push_back(root);
  while(!stack.empty()) {
    node_t* currentNode = stack.back();
    stack.pop_back();

    //add the node to the parent table with its parent's name
    node_t node;
    node.name = currentNode->name;
    node.parentName = (currentNode->up != NULL) ? currentNode->up->name : "none"; //stores the name of the parent, or "none" if there is no parent
    node.blen = currentNode->blen;
    parentTable.push_back(node);

    //adds node branches to the stack
    for(int i=0; i<currentNode->nbranches; ++i)
      stack.push_back(currentNode->branches[i]);
  }
  return parentTable;
}

void printParentTable(const vector<node_t>& parentTable) {
  for(int i=0; i<parentTable.size(); ++i) {
    cout<<"Node "<<parentTable[i].name<<": Parent="<<parentTable[i].parentName<<", Distance="<<parentTable[i].blen<<endl;
  }
}

//========================= DISTANCE MATRIX =========================//
/*
vector<vector<double>> buildDistanceMatrix(const vector<node_t>& parentTable) {
  int n = parentTable.size();
  vector<vector<double>> distanceMatrix(n, vector<double>(n, 0.0));

  // Parcourir chaque nœud de la parentTable
  for (int i = 0; i < n; ++i) {
    int parent = parentTable[i].parent;
    double distance = parentTable[i].distance;

    // Remplir la distance entre le nœud actuel et son parent
    if (parent != -1) {
      distanceMatrix[i][parent] = distance;
      distanceMatrix[parent][i] = distance;
    }
  }

  // Parcourir la distanceMatrix pour mettre à jour les distances entre les nœuds
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (distanceMatrix[i][k] != 0 && distanceMatrix[k][j] != 0) {
          double distance = distanceMatrix[i][k] + distanceMatrix[k][j];
          if (distanceMatrix[i][j] == 0 || distance < distanceMatrix[i][j]) {
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
          }
        }
      }
    }
  }

  return distanceMatrix;
}*/
/*
double** buildDistanceMatrix(const vector<node_t>& parentTable) {
  int n = parentTable.size();
  double** distanceMatrix = new double*[n];
  for (int i = 0; i < n; ++i) {
    distanceMatrix[i] = new double[n];
    memset(distanceMatrix[i], 0, sizeof(double) * n);
  }

  // Remplir la distance entre chaque nœud et son parent
  for (int i = 0; i < n; ++i) {
    int parent = parentTable[i].parent;
    double distance = parentTable[i].distance;

    if (parent != -1) {
      distanceMatrix[i][parent] = distance;
      distanceMatrix[parent][i] = distance;
    }
  }

  // Calculer les distances entre les nœuds en utilisant l'algorithme de Floyd-Warshall
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (distanceMatrix[i][k] != 0 && distanceMatrix[k][j] != 0) {
          double distance = distanceMatrix[i][k] + distanceMatrix[k][j];
          if (distanceMatrix[i][j] == 0 || distance < distanceMatrix[i][j]) {
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
          }
        }
      }
    }
  }

  return distanceMatrix;
}*/
double** convertToDistanceMatrix(const vector<node_t>& parentTable, int size) {
  // Créer une matrice symétrique carrée de taille 'size'
  double** distanceMatrix = new double*[size];
  for (int i = 0; i < size; ++i) {
    distanceMatrix[i] = new double[size];
    for (int j = 0; j < size; ++j) {
      distanceMatrix[i][j] = 0.0; // Initialiser toutes les valeurs à 0
    }
  }

  // Remplir la matrice avec les distances à partir du tableau parentTable
  for (const node_t& node : parentTable) {
    int parent = node.parent;
    double distance = node.distance;
    if (parent >= 0 && parent < size) {
      distanceMatrix[parent][node.serial] = distance;
      distanceMatrix[node.serial][parent] = distance;
    }
  }

  return distanceMatrix;
}
/*
void printDistanceMatrix(const vector<vector<double>>& distanceMatrix, const vector<node_t>& parentTable) {
  int numNodes = parentTable.size();

  // Print the distance matrix
  for (int i = 0; i < numNodes; ++i) {
    for (int j = 0; j < numNodes; ++j) {
      cout << distanceMatrix[i][j] << " ";
    }
    cout << endl;
  }
}*/
/*
void printDistanceMatrix(double** distanceMatrix, int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cout << distanceMatrix[i][j] << " ";
    }
    cout << endl;
  }
}

void freeDistanceMatrix(double** distanceMatrix, int n) {
  for (int i = 0; i < n; ++i) {
    delete[] distanceMatrix[i];
  }
  delete[] distanceMatrix;
}*/

void printDistanceMatrix(double** distanceMatrix, int size) {
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      cout << distanceMatrix[i][j] << " ";
    }
    cout << endl;
  }
}

void freeDistanceMatrix(double** distanceMatrix, int size) {
  for (int i = 0; i < size; ++i) {
    delete[] distanceMatrix[i];
  }
  delete[] distanceMatrix;
}



//================================================================
// PROGRAMME
//================================================================

int main(int argc, char* argv[]) {

  if(argc != 2) {
    cerr<<"Usage: ./ClonalTreeClustering <newick_file>"<<endl;
    return 1;
  }
  
  string filename = argv[1] ;

  // open the newick file in read mode
  ifstream file(filename);
  // check that the file has been opened correctly
  if(!file) {
      cerr<<"Error: Unable to open the file."<<endl;
      return -1;
  }
  vector<vector<char>> newickLineageWhole = readNewick(file);
  file.close();

  //browse each lineage
  for(const vector<char>& newickLineageLine : newickLineageWhole) {
    //modify the vector to work on the relationships
    vector<char> modifiedLineageLine = newickLineageLine;
    //check that the vector contains at least two elements
    if(modifiedLineageLine.size() >= 2) {
      //delete the first item by moving the following items to the left
      modifiedLineageLine.erase(modifiedLineageLine.begin());
      //delete the penultimate element
      modifiedLineageLine.erase(modifiedLineageLine.end() - 2);
    }
    // display the vector content after deletion
    cout<<"vector: ";
    for(char c : modifiedLineageLine) {
      cout<<c;
    }
    cout<<endl;

    //convert the vector to a char ptr to work on the relationships
    modifiedLineageLine.push_back('\0');
    char * newick = modifiedLineageLine.data();
    cout<<"char ptr: "<<newick<<endl;

    //==================== RELATIONSHIPS BETWEEN NODES ====================//
    fprintf(stderr,"newick: %s\n",newick);
    vector<string> vec=mysplit2(strdup(newick));
    for(int i=0;0&&i<vec.size();i++)
      fprintf(stderr,"%d) %s\n",i,vec[i].c_str());

    node_t *nd = parse(vec);
    node_t **lst = new node_t*[serial];
    serialize(nd,lst);
      
    // store the relationships in a map
    vector<node_t> parentTable = createParentTable(nd);
    printParentTable(parentTable);

    //========================= DISTANCE MATRIX =========================//
/*
    // Construire la matrice de distances
    vector<vector<double>> distanceMatrix = buildDistanceMatrix(parentTable);
      
    // Afficher la matrice de distances
    printDistanceMatrix(distanceMatrix, parentTable);*/
      
    double** distanceMatrix = convertToDistanceMatrix(parentTable, parentTable.size());
    printDistanceMatrix(distanceMatrix, parentTable.size());
    freeDistanceMatrix(distanceMatrix, parentTable.size());
  }

  return 0;
}
