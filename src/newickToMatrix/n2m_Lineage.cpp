//================================================================================================
//=  NewickToMatrix Conversion v0.1 
//=  Authors : UQAM lab
//=  Date : June 2009 
//=
//=  Description : 
//=
//= AnnARTIGES =» Temporary script for lineage trees.
//================================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <typeinfo>

#include "hgt_3.4_interactive/structures.h"
#include "hgt_3.4_interactive/utils_tree.cpp"
#include "hgt_3.4_interactive/fonctions.cpp"

#define FAIL -1
#define TRUE 1
#define FALSE 0
#define INFINI 999999.99

void newickToMatrixLineage(const char *newick,FILE *out);

// Function to check the Newick Format and detect any errors before parsing the sequence for the distance matrix
void checkFormat(const char *newickLineageLine){

  int i = 0, j = 0, k = 0, a = 0;
  char symbol;
  int cpt = 0;


  do
	{
		symbol = newickLineageLine[cpt++];

		if (symbol == ':') i++;

    if (symbol == '(') j++;
    else if (symbol == ')') j--;

    if (symbol == ';') k++;

    if (symbol == '%') a++; 

	}  while(symbol != '\0');


  if (newickLineageLine[0] != '(') { printf("Incorrect Newick file format. Newick string must begin with a '(' character."); exit(FAIL);}

  if (i == 0) { printf("Incorrect Newick file format. Edge lengths must be indicated after a ':' characters."); exit(FAIL);}

  if (j != 0) { printf("Incorrect Newick file format. Number of right parentheses must be equal to number of left parentheses."); exit(FAIL); }

  if (k == 0) { printf("Incorrect Newick file format. Newick string must be followed by a ';' character."); exit(FAIL);}
	else if (k > 1) { printf("Incorrect Newick file format. Newick string must contain (in the end) only one ';' character."); exit(FAIL);}

  if (a>0) { printf("Incorrect Newick file format. Newick string cannot contain \'%%\' character."); exit(FAIL);}

}

int nbNodesNewick(const char * newick, int& nodes){
	int i=0;
	int n = 0;
	char symbol,parent=' ';
	char symbolOld =' ';  	
	int temoin =0;

	do{
		symbol = newick[i];
		i++;
		if (symbol == ':'){
			if(symbolOld != ')' && symbolOld != ',') n++;
			else nodes++;
		}
		if(symbol >= '0' && symbol <= '9' && symbolOld==')') temoin=1;
		if(symbol==':' && temoin==1) temoin=0;
		symbolOld = symbol;
	}while(symbol != ';');

	return n;
}

int getNamesNaive(const char * newick, char ** lesNoms, char * newStr, int size){

	int i = 0, j, x = 0, y = 0, idStart = 0, idStop = 0, num = 0, root = 1;
	char symbol, *numero;
	char ** namesTemp;
	namesTemp = (char**)malloc(20*sizeof(char*));
	numero = (char*)malloc((100) * sizeof(char));

	for(i = 0; i <= size; i++)
		namesTemp[i] = (char*)malloc(50);

	i = 0;

	do{

		symbol = newick[i];

		if((symbol == '(') || (symbol == ',') || (symbol == ')'))
		{
			idStart = i;
		}

		if(symbol == ':')
		{
			// Write in the new string everything between the ":" and the last character before the name of the node.

			for(j = idStop; j <= idStart; j++)
			{
				newStr[x++] = newick[j];
			}


			// Only retrieve the name of the node and assign a number if the node already has a name.
			if(i != idStart+1)
			{
				// If the node is "naive" (root), its name goes in the first position of the vector and the number assigned is 0.
				// Also change root variable to 0 for later use. 
				if(newick[idStart+1] == 'n')
				{
					newStr[x++] = '0';
					for(j = 1; j <= 5; j++)
					{
						namesTemp[0][y++] = newick[idStart+j];
					}
					root = 0;
				}
				
				else
				{
					num++;

					for(j = idStart+1; j < i; j++)
					{
						namesTemp[num][y++] = newick[j];
					}

					namesTemp[num][y++] = '\0';
					itoa_(num, numero, 10);
					for(y=0; y<strlen(numero); y++)
					{
						newStr[x++] = numero[y];
					}			
				}
				y = 0;

				
			}
			idStop = i;
		}
		i++;

	}while(symbol != ';');


	// Write the last information of the newick line in the new string.
	for(j = idStop; j <= idStart+1; j++)
	{
		newStr[x++] = newick[j];
	}

	int names = 0;

	// If the naive cell is present root = 0 (otherwise root = 1)
	// This allows us to know at what line to start parsing the namesTemp tab.
	for(i = root; i <= num; i++)
	{
		for(j = 0; j <= strlen(namesTemp[i]); j++)
		{
			lesNoms[names][j] = namesTemp[i][j];
		}
		names++;
	}

	return names;
}

//========================================================================================================
//============================================ MAIN ======================================================
//========================================================================================================
int main(int nargc,char **argv){
	
	if (nargc != 3){
		printf("\nNombre de parametres incorrects !!");
		exit(1);
	}
	
	char * newick = (char*)malloc(100000*sizeof(char));
	int size;
	FILE * in = fopen(argv[1],"r");
	FILE * out = fopen(argv[2],"w");
	newick = readNewick(in);
	checkFormat(newick);
	printf("\nLe format est bon");
	newickToMatrixLineage(newick,out);

	fclose(in);
	fclose(out);
	
	
	
	return 0;
}
//==================================================================
//=
//==================================================================
void newickToMatrixLineage(const char *newick,FILE *out){
	int i,j;
	int pos_racine = -1;
	long int ** ARETEBcell;
	double *  LONGUEUR;
	char   ** NAMES;
	char * newString;
	double ** ADJACENCE;
	double ** ADD;
	int kt;
	int nbNodes = 0;
	int size = nbNodesNewick(newick, nbNodes);
	int fullSize = size + nbNodes;
	printf("\nNb nodes avec noms %d et nodes sans noms %d", size, nbNodes);

	newString = (char*)malloc(100000*sizeof(char));
	
	ADD = (double**)malloc(2*size*sizeof(double*));
	ADJACENCE = (double**)malloc(2*size*sizeof(double*));
	for(i=0;i<2*size;i++){
		ADD[i] = (double*)malloc(2*size*sizeof(double));
		ADJACENCE[i] = (double*)malloc(2*size*sizeof(double));
	}

	LONGUEUR = (double*)malloc((4*(size))*sizeof(double));


	ARETEBcell = (long int**)malloc(fullSize*sizeof(long int));
	for(i = 0; i <= fullSize; i++)
	{
		ARETEBcell[i] = (long int*)malloc(2); 
	}


	NAMES=(char**)malloc(2*size*sizeof(char*));
	for(i=0;i<=size;i++)
		NAMES[i] = (char*)malloc(50);
	//Method to retrieve the nodes' names and store them in a list
	int test = getNamesNaive(newick, NAMES, newString, size);
	
	for(i=0; i<size; i++){ printf("\n %s\t %d", NAMES[i], i); }

	
	printf("\nici");	
	pos_racine = lectureNewickBcell(newString,ARETEBcell,LONGUEUR,NAMES,&kt);

	printf("\n allo");
    loadAdjacenceMatrixLineage(ADJACENCE,ARETEBcell,LONGUEUR,fullSize,kt);
	printf("le \n");

    Floyd(ADJACENCE,ADD,fullSize,kt); 
    printf("\n monde \n");


	int max_taille=0;
	for(j=1;j<=size;j++){
		max_taille = (strlen(NAMES[j])>max_taille)?strlen(NAMES[j]):max_taille;
	}


	fprintf(out,"\n%d",size);
	for(i = 0; i < fullSize; i++){
		fprintf(out,"\n%s",NAMES[i]);
		if(strlen(NAMES[i])<max_taille){
			for(j = strlen(NAMES[i]); j <= max_taille; j++){
				fprintf(out," ");
			}
		}
		for(j = 0; j < fullSize; j++){
			fprintf(out,"  %3.5lf",ADD[i][j]);
		}
	}

}

//===================================================================
//= 
//===================================================================
/*void MatrixToNewick(char ** NAMES,double **ADD,int size,char *newick){
	
	char * chaine;
	int * ARETE=(long int*)malloc(4*(2*(size))*sizeof(long int));
	double * LONGUEUR=(double*)malloc((4*(size))*sizeof(double));
//	kt = Tree_edges (aTree->ADD,aTree->ARETE,aTree->LONGUEUR,n,binaire);
	SAVEASNewick(chaine,LONGUEUR,ARETE,NAMES,null, int size); 
}
*/