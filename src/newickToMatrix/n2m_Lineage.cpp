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

void newickToMatrix(const char *newick,FILE *out);

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
		//if(symbolOld == ')' && symbol == ':') nodes++;

		//if(symbol == ':' && symbolOld !=')' && temoin != 1) n++;
		if(symbol >= '0' && symbol <= '9' && symbolOld==')') temoin=1;
		if(symbol==':' && temoin==1) temoin=0;
		symbolOld = symbol;
	}while(symbol != ';');

	return n;
}

int getNamesNaive(const char * newick, char ** lesNoms, char * newStr){

	int i = 0, j, k, x = 0, y = 0, idStart = 0, idStop = 0, num = 0, root = 1;
	char symbol, *numero;
	char ** namesTemp;
	namesTemp = (char**)malloc(20*sizeof(char*));
	numero = (char*)malloc((100) * sizeof(char));

	for(k = 0; k <= 10; k++)
		namesTemp[k] = (char*)malloc(50);


	//printf("\n nouvelle str de base %s", newStr);

	do{

		symbol = newick[i];
		//printf("\n je peux au moins rentrer dans la boucle ou pas ? %d", i);

		if((symbol == '(') || (symbol == ',') || (symbol == ')'))
		{
			idStart = i;
		}

		if(symbol == ':')
		{
			// Write in the new string everything between the ":" and the last character before the name of the node.

			//printf("\n début de chaîne %d", idStop);
			//printf("\n fin de chaîne %d", idStart);
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
					newStr[x++] = numero[0];				
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

	
	//printf("\n SORTIE BOUCLE VERIFICATION %s", newStr);

	// If the naive cell is present root = 0 (otherwise root = 1)
	// This allows us to know at what line to start parsing the namesTemp tab.
	for(i = root; i <= num; i++)
	{
		for(j = 0; j <= strlen(namesTemp[i]); j++)
		{
			lesNoms[names][j] = namesTemp[i][j];
		}
		names++;
		//printf("\nverification qu'on a au moins les noms, %s", namesTemp[i]);
	}
	//printf("\n SORTIE BOUCLE VERIFICATION %s", newStr);

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
	newickToMatrix(newick,out);

	fclose(in);
	fclose(out);
	
	
	
	return 0;
}
//==================================================================
//=
//==================================================================
void newickToMatrix(const char *newick,FILE *out){
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
	printf("Nb nodes avec noms %d et nodes sans noms %d", size, nbNodes);

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

	//Temporary method before finding better solution to store names somewhere
	NAMES=(char**)malloc(2*size*sizeof(char*));
	for(i=0;i<=size;i++)
		NAMES[i] = (char*)malloc(50);
	int test = getNamesNaive(newick, NAMES, newString);
	
	for(i=0; i<size; i++){ printf("\n %s\t %d", NAMES[i], i); }

	
	printf("\nici");	
	pos_racine = lectureNewickBcell(newString,ARETEBcell,LONGUEUR,NAMES,&kt);

	printf("\n allo");
    loadAdjacenceMatrixLineage(ADJACENCE,ARETEBcell,LONGUEUR,size,kt);
	printf("le \n");

    Floyd(ADJACENCE,ADD,size,kt); 
    printf("\n monde \n");
	
	for(i=0; i<=size; i++){
		for(j=0; j<=size; j++){
			printf("%lf \t", ADD[i][j]);
		}
		printf("\n");
	}


	int max_taille=0;
	for(j=1;j<=size;j++){
		max_taille = (strlen(NAMES[j])>max_taille)?strlen(NAMES[j]):max_taille;
	}


	fprintf(out,"\n%d",size);
	for(i = 0; i <= size; i++){
		fprintf(out,"\n%s",NAMES[i]);
		if(strlen(NAMES[i])<max_taille){
			for(j = strlen(NAMES[i]); j <= max_taille; j++){
				fprintf(out," ");
			}
		}
		for(j = 0; j <= fullSize; j++){
			fprintf(out,"  %3.5lf",ADD[i][j]);
		}
	}
	//}
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
