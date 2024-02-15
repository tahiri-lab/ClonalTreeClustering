//================================================================================================
//=  NewickToMatrix Conversion v0.1 
//=  Authors : UQAM lab
//=  Date : June 2009 
//=
//=  Description : 
//=
//================================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

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
		//if (symbol == ':' && symbolOld != ')') n++;
		if(symbol == '(') nodes++;

		if(symbol == ':' && symbolOld !=')' && temoin != 1) n++;
		if(symbol >= '0' && symbol <= '9' && symbolOld==')') temoin=1;
		if(symbol==':' && temoin==1) temoin=0;
		symbolOld = symbol;
	}while(symbol != ';');

	return n;
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
	int pos_racine=-1;
	long int *  ARETE;
	double *  LONGUEUR;
	char   ** NAMES;
	double ** ADJACENCE;
	double ** ADD;
	int kt;
	int nbNodes = 0;
	int size = nbNodesNewick(newick, nbNodes);
	int fullSize = size + nbNodes;
	
	ADD = (double**)malloc(2*size*sizeof(double*));
	ADJACENCE = (double**)malloc(2*size*sizeof(double*));
	for(i=0;i<2*size;i++){
		ADD[i] = (double*)malloc(2*size*sizeof(double));
		ADJACENCE[i] = (double*)malloc(2*size*sizeof(double));
	}
	ARETE=(long int*)malloc(4*(2*(size))*sizeof(long int));
	LONGUEUR=(double*)malloc((4*(size))*sizeof(double));
	
	NAMES=(char**)malloc(2*size*sizeof(char*));
	for(i=0;i<=size;i++)
		NAMES[i] = (char*)malloc(50);
printf("ici");	
	pos_racine = lectureNewick(newick,ARETE,LONGUEUR,NAMES,&kt);

	//for(i=0; i<=size; i++){ printf("\n %lf\t %d", LONGUEUR[i], i); }
	
	printf("\n allo");
    loadAdjacenceMatrix(ADJACENCE,ARETE,LONGUEUR,size,kt);
    printf("le \n");

	//printf("\nListe toutes les branches :");
	//for(i=0; i<=3*size; i++) {printf("\n ARETE %d\t %ld", i, ARETE[i]);}
	//printf("\nListe toutes les distances :");
	//for(i = 0; i <= size+nbNodes; i++) {printf("\n LONGUEUR %d\t %f", i, LONGUEUR[i]);}


    Floyd(ADJACENCE,ADD,size,kt); 
    	printf("\n monde \n");
	/*
	for(i=1; i<=2*size-1; i++){
		for(j=1; j<=2*size-1; j++){
			printf("%lf \t", ADD[i][j]);
		}
		printf("\n");
	}*/


	//if(pos_racine != -1){
		//=recherche de la plus longue taille de nom d'especes
		int max_taille=0;
		for(j=1;j<=size;j++){
			max_taille = (strlen(NAMES[j])>max_taille)?strlen(NAMES[j]):max_taille;
		}
		
		/*for(i=size+1; i<=2*size-2; i++){
			//char * c = i;
			char c = static_cast<char>(i);
			NAMES[i] = c;
		}
		for(i=1; i<=2*size; i++){printf("\n le nom : %s", NAMES[i]);}*/


		fprintf(out,"\n%d",size);
		for(i=1;i<=size;i++){
			fprintf(out,"\n%s",NAMES[i]);
			if(strlen(NAMES[i])<max_taille){
				for(j=strlen(NAMES[i]);j<=max_taille;j++){
					fprintf(out," ");
				}
			}
			for(j=1;j<=fullSize;j++){
				fprintf(out,"  %3.5lf",ADD[i][j]);
			}
		}
		for(i=size+1; i<=fullSize; i++){
			fprintf(out, "\n%d", i);
			for(j=1; j<=fullSize; j++){
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
