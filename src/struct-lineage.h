#include <vector>
using namespace std;

//================================================================
// STRUCTURES
//================================================================

struct LNode {	
  	int NoNode;
  	struct TNode **child;
  	int nbchildren;
  };

struct InputLineage{
  vector<vector<char>> idSeq;
  //int root;
  int size;
  vector<int> abundance;
  vector<int> degree;
  vector<int> depth;
  vector<long int> edge;
  vector<vector<double>> length;
  double ** ADJACENCE;
};

//================================================================
// CONSTANTES
//================================================================
#define SPECIES 1
#define GENE 2
#define SPECIES_NAME_LENGTH 50