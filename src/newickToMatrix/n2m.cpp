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
	int size = nbSpeciesNewick(newick);
	
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
	printf("allo");
    loadAdjacenceMatrix(ADJACENCE,ARETE,LONGUEUR,size,kt);
    printf("le");
    Floyd(ADJACENCE,ADD,size,kt); 
    	printf("monde");
	//if(pos_racine != -1){
		//=recherche de la plus longue taille de nom d'especes
		int max_taille=0;
		for(j=1;j<=size;j++){
			max_taille = (strlen(NAMES[j])>max_taille)?strlen(NAMES[j]):max_taille;
		}
		fprintf(out,"\n%d",size);
		for(i=1;i<=size;i++){
			fprintf(out,"\n%s",NAMES[i]);
			if(strlen(NAMES[i])<max_taille){
				for(j=strlen(NAMES[i]);j<=max_taille;j++){
					fprintf(out," ");
				}
			}
			for(j=1;j<=size;j++){
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
