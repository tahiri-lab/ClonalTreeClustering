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
#include <sstream>
#include <fstream>
#include <vector>

#include "hgt_3.4_interactive/structures.h"
#include "fonctionsLineage.cpp"

using namespace std;

#define FAIL -1
#define TRUE 1
#define FALSE 0
#define INFINI 999999.99


//========================================================================================================
//============================================ MAIN ======================================================
//========================================================================================================
int main(int nargc,char **argv){
	
	if (nargc < 3){
		printf("\nPas assez de paramètres !");
		exit(1);
	}
	char * format = argv[1];
	if (format[0] != '-') {
		printf("Incorrect argument for the format of the input \nCommand format : ./n2m -s X FILE output_file.txt \n\t\t ./n2m -m FILE_1 FILE_2 output_file.txt");
		printf("\n Argument must be \n -s X \t if one file contains all the Newick sequences with X the number of sequences \n -m \t if multiple files contain one Newick sequence each \n");
		exit(FAIL);
	}

	FILE * out = fopen(argv[nargc-1],"w");

	int nbTrees;
	char * newick = (char*)malloc(100000*sizeof(char));
	char ** treeNames;
	int i, j;
	double ** forest;
	string delimiter = ".";
	//double **ADD = nullptr;
	double ***Connect;
	double *** lineageForest;
	std::vector<std::map<std::string, int>> weightForest;
	std::vector<std::map<std::string, int>> namesForest;

	if(format[1] == 'm')
	{

		nbTrees = nargc - 3;
		printf("\t nombre d'arbres %d", nbTrees);
		lineageForest = (double***)malloc(2*nbTrees*sizeof(double**));
		Connect = (double***)malloc(2*nbTrees*sizeof(double**));
		treeNames = (char**)malloc(nbTrees*sizeof(char*));
		for(i = 0; i < 2*nbTrees; i++)
		{
			lineageForest[i] = nullptr;
			Connect[i] = nullptr;
		}
		
		for(i = 2; i < nargc - 1; i++)
		{
			FILE * in = fopen(argv[i], "r");
			char * file = argv[i];
			char * treeId = (char*)malloc(20*sizeof(char));
			
			newick = readNewick(in);
			class map <string, int> dicAbond;
			class map <string, int> dicNames;
			j = 0;
			while(file[j] != '.'){
				treeId[j] = file[j];
				j++;
			}
			treeNames[i-2] = treeId;

			checkFormat(newick);
			newickToMatrixLineage(newick, out, dicNames, dicAbond, lineageForest[i-2], Connect[i-2], treeId);
			weightForest.push_back(dicAbond);
			namesForest.push_back(dicNames);
			
			free(newick);
			fclose(in);
		}
	}

	else if (format[1] == 's')
	{
		nbTrees = atoi(argv[2]);
		lineageForest = (double***)malloc(2*nbTrees*sizeof(double**));
		Connect = (double***)malloc(2*nbTrees*sizeof(double**));
		treeNames = (char**)malloc(2*nbTrees*sizeof(char*));
		for(i = 0; i < 2*nbTrees; i++)
		{
			lineageForest[i] = nullptr;
			Connect[i] = nullptr;
		}
		string line;
		string in = argv[3];
		ifstream file(in);
		i = 0;
		char * treeId;
		char * oldTreeId;
		while(getline(file, line))
		{
			if(line.find("((") != string::npos)
			{
				if(treeId == oldTreeId || treeId[0] == 0) { treeId = "unnamed"; }
				treeNames[i] = treeId;
				char * newick = new char[line.length() + 1];
				strcpy(newick, line.c_str());
				class map <string, int> dicAbond;
				class map <string, int> dicNames;

				checkFormat(newick);
				newickToMatrixLineage(newick, out, dicNames, dicAbond, lineageForest[i], Connect[i], treeId);
				weightForest.push_back(dicAbond);
				namesForest.push_back(dicNames);
				i++;
				oldTreeId = treeId;
			}
			else
			{
				treeId = new char[line.length() + 1];
    			strcpy(treeId, line.c_str());
			}
		}
		file.close();
	}


	forest = (double**)malloc(2*nbTrees*sizeof(double*));
	for(i=0;i<2*nbTrees;i++){
		forest[i] = (double*)malloc(2*nbTrees*sizeof(double));
		forest[i][i] = 0.00;
	}

	for(i = 0; i < nbTrees-1; i++)
	{
		class map <string, int> dicAbond = weightForest[i];
		class map <string, int> dicNames = namesForest[i];
		
		for(j = i+1; j <= nbTrees-1; j++)
		{

			class map <string, int> dicAbond2 = weightForest[j];
			class map <string, int> dicNames2 = namesForest[j];
			
			double trying = calculMetric(lineageForest[i], lineageForest[j], Connect[i], Connect[j], dicNames, dicNames2, dicAbond, dicAbond2);
			forest[i][j] = trying;
			forest[j][i] = trying;
		}
	}
	
	fprintf(out, "\n\n Comparison Matrix");
	for(i = 0; i < nbTrees; i++)
	{
		fprintf(out, "\n%s", treeNames[i]);
		for(j = 0; j < nbTrees; j++)
		{ fprintf(out, "\t %f", forest[i][j]); }
	}
	printf("\n");

	fclose(out);
	
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
