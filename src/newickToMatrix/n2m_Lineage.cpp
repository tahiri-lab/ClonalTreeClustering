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
#include <string>
#include <map>
#include <iostream>

#include "hgt_3.4_interactive/structures.h"
#include "hgt_3.4_interactive/utils_tree.cpp"
#include "hgt_3.4_interactive/fonctions.cpp"

using namespace std;

#define FAIL -1
#define TRUE 1
#define FALSE 0
#define INFINI 999999.99


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
	int i = 0;
	int n = 0;
	char symbol,parent = ' ';
	char symbolOld = ' ';  	
	int temoin =0;

	do{
		symbol = newick[i];
		i++;
		if (symbol == ':'){
			//printf("\n symbol vaut : %d", newick[i-2]);
			if(symbolOld != ')' && symbolOld != ',' && symbolOld != '(' && temoin != 2) { n++; }
			else {
				nodes++;
				temoin = 0;
			}
		}
		if(symbol >= '0' && symbol <= '9' && symbolOld==')') {
			temoin=1;
			//printf("\nrentre dans la boucle ou on assigne valeur 1 à témoin");
			}
		if(symbol==':' && temoin==1) {
			temoin=0;
			//printf("\nrentre dans la boucle ou on assigne valeur 0 à témoin");
			}

		if(symbol == ')' || symbol == ',' || symbol == '(')
		{
			if(newick[i] == '@' || newick[i] == ':') { temoin = 2; }
		}
		symbolOld = symbol;
	}while(symbol != ';');

	return n;
}

void getNamesNaive(const char * newick, char ** lesNoms, char * newStr, int size, std::map <std::string, int> &abondMap, std::map <std::string, int> &namesMap){

	int i = 0, j, x = 0, y = 0, idStart = 0, idStop = 0, num = 0, endId;
	int root = 1, temoin_ab = 0, stop = 0, a, b, abondance;
	char symbol, *numero, * key_names, *abond;
	char ** namesTemp;
	namesTemp = (char**)malloc(20*sizeof(char*));
	numero = (char*)malloc((100) * sizeof(char));
	abond = (char*)malloc((20) * sizeof(char));

	for(i = 0; i <= size; i++)
		{namesTemp[i] = (char*)malloc(50);}

	i = 0;

	do{

		symbol = newick[i];

		if((symbol == '(') || (symbol == ',') || (symbol == ')'))
		{
			idStart = i;
		}

		if(symbol == '@')
		{
			temoin_ab = 1;
			stop = i;
			a = i+1;
			b = 0;

			while(newick[a] != ':')
			{
				abond[b++] = newick[a];
				a++;
			}
			abondance = atoi(abond);

		}

		if(symbol == ':')
		{
			// Write in the new string everything between the ":" and the last character before the name of the node.
			for(j = idStop; j <= idStart; j++)
			{
				newStr[x++] = newick[j];
			}
			// Only retrieve the name of the node and assign a number if the node already has a name.
			if(i != idStart+1 && newick[idStart+1] != '@')
			{
				if(temoin_ab == 1) {
					endId = stop;
				}
				else {
					endId = i;
					abondance = 1;
				}
				// If the node is "naive" (root), its name goes in the first position of the vector and the number assigned is 0.
				// Also change root variable to 0 for later use. 
				if(newick[idStart+1] == 'n')
				{
					newStr[x++] = '0';
					/*for(j = 1; j <= 5; j++)
					{
						namesTemp[0][y] = newick[idStart+j];
						y++;
					}*/
					namesTemp[0] = "naive";
					key_names = namesTemp[0];
					namesMap[key_names] = 0;
					abondMap[key_names] = abondance;
					
					root = 0;
					y=0;
				}
				
				else
				{
					num += 1;

					for(j = idStart+1; j < endId; j++)
					{
						namesTemp[num][y] = newick[j];
						y++;
					}

					namesTemp[num][y++] = '\0';
					key_names = namesTemp[num];
					namesMap[key_names] = num;
					abondMap[key_names] = abondance;
					itoa_(num, numero, 10);
					for(y=0; y<strlen(numero); y++)
					{
						newStr[x++] = numero[y];
					}
					temoin_ab = 0;		
				}
				y = 0;
				idStop = i;
			}
			else if (i != idStart+1 && newick[idStart+1] == '@') {
				idStop = stop;} 
			else if (i == idStart+1) {
				newStr[x++] = '@';
				newStr[x++] = '1';
				idStop = i;
				temoin_ab = 0;
				}
		temoin_ab = 0;
			
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

	auto it = namesMap.begin();
	while(it != namesMap.end())
	{
		if(it->second == 0 && it->first != "naive")
		{
			it = namesMap.erase(it);
		}
		else {it++;}
	}

	//return namesMap;
}

void newickToMatrixLineage(const char *newick,FILE *out, std::map <std::string, int>& dicNames, std::map <std::string, int>& dicAbond, double **&ADD){
	int i,j;
	int pos_racine = -1;
	int dist_naive;
	long int ** ARETEBcell;
	double *  LONGUEUR;
	char ** NAMES;
	char * newString;
	double ** ADJACENCE;
	int kt;
	int nbNodes = 0;
	int size = nbNodesNewick(newick, nbNodes);
	int fullSize = size + nbNodes;
	printf("\n\nNb nodes avec noms %d et nodes sans noms %d", size, nbNodes);


	newString = (char*)malloc(100000*sizeof(char));
	
	ADD = (double**)malloc(2*size*sizeof(double*));
	ADJACENCE = (double**)malloc(2*size*sizeof(double*));
	for(i=0;i<2*size;i++){
		ADD[i] = (double*)malloc(2*size*sizeof(double));
		ADJACENCE[i] = (double*)malloc(2*size*sizeof(double));
	}

	LONGUEUR = (double*)malloc((4*(size))*sizeof(double));


	// Modification to the initialization of ARETEBcell
	// Only way I found to initialize it without having any problems for the following trees
	ARETEBcell = new long int*[fullSize];
	for(i = 0; i < fullSize; i++)
	{
		ARETEBcell[i] = new long int[2];
		ARETEBcell[i][0] = 0;
		ARETEBcell[i][1] = 0;
	}

	NAMES=(char**)malloc(2*size*sizeof(char*));
	for(i=0;i<=size;i++)
		NAMES[i] = (char*)malloc(50);
	//Method to retrieve the nodes' names and store them in a list
	getNamesNaive(newick, NAMES, newString, size, dicAbond, dicNames);
	//printf("\nLa récupération des noms est faite !");
	//printf("\n Nouvelle séquence Newick: %s", newString);
	

	
	//printf("\nici");	
	dist_naive = lectureNewickBcell(newString,ARETEBcell,LONGUEUR,NAMES,&kt, size, dicNames, dicAbond);
	//for(i=0; i < fullSize; i++){printf("\n okay let's see %s   %d", NAMES[i], i);}
	/*for(i = 0; i < fullSize; i++){
		char *key = NAMES[i];
		printf("\nAbundancy %s\t %d", key, dicAbond[key]);
	}*/

	printf("\n allo ");
    loadAdjacenceMatrixLineage(ADJACENCE, ARETEBcell, LONGUEUR, fullSize-1, kt);
	printf("le ");

    FloydLineage(ADJACENCE, ADD, fullSize-1, kt, dist_naive); 
    printf("monde \n");


	int max_taille=0;
	for(j = 1; j <= size;j++){
		max_taille = (strlen(NAMES[j])>max_taille)?strlen(NAMES[j]):max_taille;
	}


	fprintf(out,"\n\n\n%d",size);
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

	free(ADJACENCE);
	delete[] ARETEBcell;
	free(LONGUEUR);
	free(newString);

}

float calculMetric(double ** ADDT1, double ** ADDT2, std::map <std::string, int> namesT1, std::map <std::string, int> namesT2,
	std::map <std::string, int> abondT1, std::map <std::string, int> abondT2){

	int j, i = 1, sizeT1 = namesT1.size(), sizeT2 = namesT2.size(), a = 0, b = 0;
	float Penality, comNod = 0, totNod;
	float Weight = 0.0, wt1, wt2;
	float Dist = 0.0, dt1, dt2;
	float metric;
	char ** TN;
	TN = (char**)malloc(20*sizeof(char*));

	/*=======================
			PENALITY
			Retrieve number of common nodes, total nodes and the names of all nodes.
	=======================*/
	// The total number of nodes (without repetitions) is the sum of nodes in each trees minus the number of common nodes.
	totNod = sizeT1 + sizeT2;


	auto it = namesT1.begin(), it2 = namesT2.begin();
	if(namesT1.count("naive") == 1 && namesT2.count("naive") == 1) {
		comNod++;
		totNod--;
		TN[0] = "naive";
		it++;
		it2++;
	}

	for(it; it != namesT1.end(); it++)
	{
		const std::string& key = it->first;

		char* node = new char[key.length() + 1];
    	strcpy(node, key.c_str());
		TN[i] = node;

		if(namesT2[key] != 0 && node != "naive"){
			comNod++;
			totNod--;

		}
		else {
			namesT2[key] = sizeT2 + b;
			b++;
		}
		i++;
	}

	for(it2; it2 != namesT2.end(); it2++)
	{
		const std::string& key = it2->first;
		char* node = new char[key.length() + 1];
    	strcpy(node, key.c_str());

		if((namesT1[key] == 0) && node[0]!='n'){
			TN[i] = node;
			namesT1[key] = sizeT1 + a;
			a++;
			i++;
		}
	}
	Penality = (comNod/totNod);
	printf("\n Penalité vaut : %f", Penality);


	/*=======================
			ABUNDANCE
			calculate the abundance difference in the trees
	=======================*/
	for(i = 0; i < totNod; i++)
	{
		wt1 = abondT1[TN[i]];
		wt2 = abondT2[TN[i]];
		//printf("\nChecking le node %s à la pos %d et les poids dans T1 %f et dans T2 %f", TN[i], i, wt1, wt2);
		Weight += abs(wt1 - wt2);
	}

	/*=======================
			DISTANCE
			calculate the distances difference between the trees
	=======================*/
	for(i = 0; i < totNod-1; i++)
	{
		char* nodeI = TN[i];

		for(j = i+1; j < totNod; j++)
		{
			char* nodeJ = TN[j];
			dt1 = 0;
			dt2 = 0;

			if(namesT1.count(nodeI) == 1 && namesT1.count(nodeJ) == 1)
			{
				dt1 = ADDT1[namesT1[nodeI]][namesT1[nodeJ]];
			}
			
			if(namesT2.count(nodeI) == 1 && namesT2.count(nodeJ) == 1)
			{
				dt2 = ADDT2[namesT2[nodeI]][namesT2[nodeJ]];
			}
			Dist += pow((dt1 - dt2), 2);
		}
		
		//printf("\n Distance %f au noeud num %d \n", Dist, i);
	}
	Dist = sqrt(Dist);

	metric = Penality * (Weight + Dist);

	printf("\n okay, donc on a Weight %f et Distance %f", Weight, Dist);
	printf("\n la métrique vaut : %f", metric);

	free(TN);
	return metric;

}

//========================================================================================================
//============================================ MAIN ======================================================
//========================================================================================================
int main(int nargc,char **argv){
	
	if (nargc < 3){
		printf("\nPas assez de paramètres !");
		exit(1);
	}
	//printf("\n the hell ? %d", nargc);
	
	char * newick1 = (char*)malloc(100000*sizeof(char));
	char * newick2 = (char*)malloc(100000*sizeof(char));
	int i, j;
	double ** forest;
	double **ADD = nullptr;
	double **ADD2 = nullptr;

	forest = (double**)malloc(2*nargc*sizeof(double*));
	for(i=0;i<2*nargc;i++){
		forest[i] = (double*)malloc(2*nargc*sizeof(double));
	}

	FILE * out = fopen(argv[nargc-1],"w");

	for(i = 1; i < nargc-2; i++)
	{

		class map <string, int> dicAbond;
		class map <string, int> dicNames;

		FILE * in1 = fopen(argv[i], "r");
		newick1 = readNewick(in1);
		checkFormat(newick1);
		//printf("\nT1");
		newickToMatrixLineage(newick1, out, dicNames, dicAbond, ADD);

		for(j = i+1; j <= nargc-2; j++)
		{

			class map <string, int> dicAbond2;
			class map <string, int> dicNames2;
			
			FILE * in2 = fopen(argv[j], "r");
			newick2 = readNewick(in2);
			checkFormat(newick2);
			newickToMatrixLineage(newick2, out, dicNames2, dicAbond2, ADD2);
			double trying = calculMetric(ADD, ADD2, dicNames, dicNames2, dicAbond, dicAbond2);
			forest[i][j] = trying;
			forest[j][i] = trying;
			free(newick2);
			free(ADD2);
			fclose(in2);
		}

		free(newick1);
		fclose(in1);
	}
	
	for(i = 1; i <= nargc-2; i++)
	{
		printf("\n %d", i);
		for(j = 1; j <= nargc-2; j++)
		{printf("\t %f", forest[i][j]);}
	}
	printf("\n");

	//fclose(in1);
	//fclose(in2);
	fclose(out);
	//fclose(out2);
	
	
	
	return 0;
}
//==================================================================
//=
//==================================================================


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
