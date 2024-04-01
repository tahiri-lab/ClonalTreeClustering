// 11/03/2024
// ARTIGES Anna
/*
File containing all the functions needed to compare lineage trees. Some of the functions are the same as the ones in the hgt_3.4_interactive/utils_tree.cpp.

Only xto, itoa_ and readNewick have not been modified from the original repository hgt_3.4_interactive.
FloydLineage, loadAdjacenceMatrixLineage and lectureNewickBcell are copies from the Floyd, loadAdjacenceMatrix and lectureNewick functions with modifications
in order to work for lineage trees.
checkFormat, nbNodesNewick, getNamesNaive and calculMetric are functions written/developped specifically for this program.
*/

static void xtoa (unsigned long val,char *buf,unsigned radix,int is_neg){
	char *p;                /* pointer to traverse string */
	char *firstdig;         /* pointer to first digit */
	char temp;              /* temp char */
	unsigned digval;        /* value of digit */

	p = buf;

	if (is_neg) {
		/* negative, so output '-' and negate */
		*p++ = '-';
		val = (unsigned long)(-(long)val);
	}

	firstdig = p;           /* save pointer to first digit */

	do {
		digval = (unsigned) (val % radix);
		val /= radix;       /* get next digit */

		/* convert to ascii and store */
		if (digval > 9)
			*p++ = (char) (digval - 10 + 'a');  /* a letter */
		else
			*p++ = (char) (digval + '0');       /* a digit */
	} while (val > 0);

	/* We now have the digit of the number in the buffer, but in reverse
	order.  Thus we reverse them now. */

	*p-- = '\0';            /* terminate string; p points to last digit */

	do {
		temp = *p;
		*p = *firstdig;
		*firstdig = temp;   /* swap *p and *firstdig */
		--p;
		++firstdig;         /* advance to next two digits */
	} while (firstdig < p); /* repeat until halfway */
}

char * itoa_(int val,char *buf,int radix){
	if (radix == 10 && val < 0)
		xtoa((unsigned long)val, buf, radix, 1);
	else
		xtoa((unsigned long)(unsigned int)val, buf, radix, 0);
	return buf;
}

/*/=================================================================================================
// Storing the newick sequence from a file in a character string
//=================================================================================================*/
char * readNewick(FILE *in){
	
	int cpt=0;
	char c;
	char * newick = (char*)malloc(100000);

	do{
		c=(char)fgetc(in);
		newick[cpt] = c;
		cpt++;
	}while(c!=';');
	
	newick[cpt] = '\0';

	return newick;
}

/*/=================================================================================================
// Checking the format of the Newick sequence
// Some checks to add for lineage sequences (if @Weight, node must have a name)
//=================================================================================================*/
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

/*/=================================================================================================
// Retrieving the number of nodes in a Newick sequence (distinction between nodes with or without a name)
//=================================================================================================*/
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
			// n = number of nodes with a name and a weight (there shouldn't be a node with a name and distance but without weight)
			// nodes = number of nodes without a name or weight (not really a node)
			if(symbolOld != ')' && symbolOld != ',' && symbolOld != '(' && temoin != 2) { n++; }
			
			else {
				nodes++;
				temoin = 0;
			}
		}
		if(symbol >= '0' && symbol <= '9' && symbolOld==')') {
			temoin=1;
			}
		if(symbol==':' && temoin==1) {
			temoin=0;
			}

		if(symbol == ')' || symbol == ',' || symbol == '(')
		{
			if(newick[i] == '@' || newick[i] == ':') { temoin = 2; }
		}
		symbolOld = symbol;
	}while(symbol != ';');

	return n;
}

/*/=================================================================================================
// Storing the nodes' names in a dictionary and assigning keys to the names, also stores the weights of nodes with name
//=================================================================================================*/
void getNamesNaive(const char * newick, char ** &lesNoms, char * newStr, int size, std::map <std::string, int> &abondMap, std::map <std::string, int> &namesMap){

	int i = 0, j, x = 0, y = 0, idStart = 0, idStop = 0, num = 1, endId, pos;
	int root = 1, temoin_ab = 0, stop = 0, a, b, abondance;
	char symbol, *numero, * key_names, *abond;
	char ** namesTemp;

	namesTemp = (char**)malloc(20*sizeof(char*));
	key_names = (char*)malloc(20*sizeof(char*));
	numero = (char*)malloc((100) * sizeof(char));
	abond = (char*)malloc((20) * sizeof(char));

	for (i = 0; i <= size; i++)
		{ namesTemp[i] = (char*)malloc(50); }

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
				a = 0;
				if(temoin_ab == 1) { endId = stop; }
				else { endId = i; }

				for(j = idStart+1; j < endId; j++)
				{
					key_names[y] = newick[j];
					y++;
				}
				key_names[y] = '\0';

				for(j = stop+1; j < i; j++)
				{ 	abond[a++] = newick[j];  }
				abondance = atoi(abond);
				
				// If the node is "naive" (root), its name goes in the first position of the vector and the number assigned is 0.
				// Also change root variable to 0 for later use, we assign it to 0 so that it is in the same position for all matrices / trees
				if(strcmp( key_names, "naive" ) == 0)
				{
					pos = 0;
					root = 0;
				}
				else {
					pos = num;
					num++;
				}

				itoa_(pos, numero, 10);
				for(y = 0; y < strlen(numero); y++)
				{
					newStr[x++] = numero[y];
				}

				y = 0;
				for(j = 0; j <= strlen(key_names); j++)
				{
					lesNoms[pos][y++] = key_names[j];
				}
				namesMap[key_names] = pos;
				abondMap[key_names] = abondance;
				y = 0;
				idStop = i;
			}

			else if (i != idStart+1 && newick[idStart+1] == '@') {
				// This case should not happen but this is kept in case of a mistake
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

	//int names = 0;

	// If the naive cell is present root = 0 (otherwise root = 1)
	// This allows us to know at what line to start parsing the namesTemp tab.
	/*for(i = root; i < num; i++)
	{
		for(j = 0; j < strlen(namesTemp[i]); j++)
		{
			lesNoms[names][j] = namesTemp[i][j];
		}
		printf("\n Checking %s", namesTemp[i]);
		names++;
	}*/

	// Remove any key with the value 0 (except naive) from the names' dictionary
	// A 0 means that the key was created in the dictionary but not assigned any value
	auto it = namesMap.begin();
	while(it != namesMap.end())
	{
		if(it->second == 0 && it->first != "naive")
		{
			it = namesMap.erase(it);
		}
		else {it++;}
	}

}


/*/=================================================================================================
// Calculating the comparison metric between two lineage trees
// Returns the metric if the trees have at least 3 common nodes, returns -1 otherwise
//=================================================================================================*/
float calculMetric(double ** ADDT1, double ** ADDT2, double ** AdjT1, double ** AdjT2,
	std::map <std::string, int> namesT1, std::map <std::string, int> namesT2, std::map <std::string, int> abondT1, std::map <std::string, int> abondT2){

	int j, i = 1, sizeT1 = namesT1.size(), sizeT2 = namesT2.size(), a = 1, b = 1;
	float Penality, comNod = 0, totNod;
	float Weight = 0.0, wt1, wt2;
	float Dist = 0.0, Dij_T1, Dij_T2;
	float Connect = 0.0, con_T1, con_T2;
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
		std::size_t test = key.find("node");
		if( test >= 0)
		{
			char* node = new char[key.length() + 1];
    		strcpy(node, key.c_str());
			TN[i] = node;
			if(namesT2[key] != 0 && node != "naive"){

				comNod++;
				totNod--;

			}
			else {
				namesT2[key] = sizeT2 + 2*b;
				b++;
			}
		}
		i++;
	}
	if(comNod >= 3)
	{
		for(it2; it2 != namesT2.end(); it2++)
		{
			const std::string& key = it2->first;
			char* node = new char[key.length() + 1];
    		strcpy(node, key.c_str());

			if((namesT1[key] == 0) && node !="naive"){
				TN[i] = node;
				namesT1[key] = sizeT1 + a;
				a++;
				i++;
			}
		}
		Penality = (comNod/totNod);
		//printf("\n Penality is : %f", Penality);


		/*=======================
			ABUNDANCE
			calculate the abundance difference in the trees
		=======================*/
		for(i = 0; i < totNod; i++)
		{
			wt1 = abondT1[TN[i]];
			wt2 = abondT2[TN[i]];
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

				Dij_T1 = ADDT1[namesT1[nodeI]][namesT1[nodeJ]];
				Dij_T2 = ADDT2[namesT2[nodeI]][namesT2[nodeJ]];
				if(isinf(Dij_T1)) { Dij_T1 = 0; }
				if(isinf(Dij_T2)) { Dij_T2 = 0; }


				con_T1 = AdjT1[namesT1[nodeI]][namesT1[nodeJ]];
				con_T2 = AdjT2[namesT2[nodeI]][namesT2[nodeJ]];
				//if(isinf(con_T1)) { con_T1 = 0; }
				//if(isinf(con_T2)) { con_T2 = 0; }

				Connect += abs(con_T1 - con_T2);
				Dist += pow((Dij_T1 - Dij_T2), 2);
				//printf("\n On regarde les nodes %s et %s dans arbre 1 %f et arbre 2 %f", nodeI, nodeJ, con_T1, con_T2);
			}
		}
		Dist = sqrt(Dist);

		metric = Penality * (Weight + Dist + Connect);

		//printf("\n We have Distance %f and Adjacence %f", Dist, Connect);
		//printf("\t The metric is : %f", metric);

		free(TN);
		return metric;
	}
	else {return -1;}

}

/*/=================================================================================================
// Initializing the distance matrix
//=================================================================================================*/
void loadAdjacenceMatrixLineage( double **tempDist, double **Adjacence, long int **ARETEB, double *LONGUEUR, int size, int kt){
	
	int i,j;
	
	for(i = 0; i <= size; i++) /*/(n+1)*/
	{	for(j = 0; j <= size; j++){
			tempDist[i][j] = tempDist[j][i] = INFINI;
			Adjacence[i][j] = Adjacence[j][i] = 0;

		}
	}
	
	for(i = 0; i <= size; i++){
		tempDist[ARETEB[i][0]][ARETEB[i][1]] = LONGUEUR[i];
		tempDist[ARETEB[i][1]][ARETEB[i][0]] = LONGUEUR[i];
		Adjacence[ARETEB[i][0]][ARETEB[i][1]] = 1;
		Adjacence[ARETEB[i][1]][ARETEB[i][0]] = 1;
	}
	
}

/*/=================================================================================================
// Calculating the full distance matrix for one tree
//=================================================================================================*/
void FloydLineage(double ** tempDist , double ** DIST,int n,int kt, int dist_naive)
{
	int i, j, k, root;

	for(i = 0; i <= n; i++)
		for(j = 0; j <= n; j++)
		{
			if(i == j)
				DIST[i][j] = 0;
			else
				DIST[i][j] = tempDist[i][j];
		}


		for(i = 0; i <= n; i++)
			for(j = 0; j <= n; j++)
				for(k = 0; k <= n; k++)
				{
					if((DIST[j][i] + DIST[i][k]) < DIST[j][k])
					{
						if(i == 0) { root = dist_naive; }
						else { root = 0; }

						DIST[j][k] = DIST[j][i] + DIST[i][k] - 2*root;
					}
						
				}
				
}

/*/=================================================================================================
// Reading the Newick sequence and retrieving the branches' distance
//=================================================================================================*/
int lectureNewickBcell(const char * newick, long int ** ARETEB, double * LONGUEUR, char ** &lesNoms, int *kt, int size,
						std::map <std::string, int>& namesMap, std::map <std::string, int>& abondMap)
{
	/*	19-02-2024 Modifications by AnnArtiges and comments in English
	This function is supposed to read a lineage tree in the Newick format and return the distances between all nodes in the tree.
	It is based on the previous "lectureNewick" function (see above).
	*/

	// TODO: Add your command handler code here
	int n;                                     
	int cpt_x;
	int i, j, k, a, a1, a2, a3, VertexNumber, numero, idSeq, ancetre, suiv, no_name = 1;
	int abondance;
	char symbol, *string, *string1, *string2, *string4, *anc, *int_node, *key_names;
	char *abond;
	char symbolOld =' ';
	int zz, xx, jj, ii;
	double longueur, dist_root = 0;		// Consider the naive cell as the root of the lineage tree
	char * tempString;
	int cpt=0;
	string4 = (char*)malloc((100000) * sizeof(char));
	abond = (char*)malloc((20) * sizeof(char));

	int temoin = 0;
	int cpt_parenthese = 0;

	i=0;
	n = 0;

	do
	{
		symbol = newick[cpt++];
		if (symbol==':') i++;
		if(symbol == ':' && symbolOld !=')' && symbolOld !=',') n++;
		if(symbol >= '0' && symbol <= '9' && symbolOld==')') temoin=1;
		if(symbol==':' && temoin==1) temoin=0;
		symbolOld = symbol;
	}  while(symbol != '\0');

	cpt=0;

	if(i<=2*n-3)(*kt)=i;
	else (*kt)=2*n-3;

	k=0;
	a=0;
	cpt=0;
	VertexNumber = size;
	int pos_racine = -1;

	string = (char*)malloc((100000) * sizeof(char));
	string2 = (char*)malloc((100000) * sizeof(char));
	string1 = (char*)malloc((100000) * sizeof(char));
	anc = (char*)malloc((100) * sizeof(char));


	do{		
		symbol = newick[cpt++];
		if ((symbol != ' ')&&(symbol != '\n')&&(symbol != '\t'))
		{
			string[a++] = symbol;
		}
	}while(symbol !='\0');

	numero = 0;
	while (string[0] == '(')   // traiter toute la chaine
	{
		a1 = 0;
		a2 = 0;
		temoin = 0;
		while( string[a2] != ')')  // traiter la paire () la plus profonde
		{
			if(string[a2] == '('){a1 = a2;}  // retrouver la parenthèse ouvrante
			a2++;
		}

		suiv = int(a2) + 1;
		int next_node = a2 + 1;
		i = 0;
		zz = a1 +1;
		
		/*Retrieve the name of the ancestor-node.
		 Specifically to get whether or not it is the naive cell, in which case we retrieve the distance associated and store it in the dist_root variable.*/
		if(string[suiv] == '@' || string[suiv] == ':')
		{
			if(string[suiv] == '@')
			{
				int ab = suiv + 1;
				int ba = 0;
				while(string[ab] != ':')
				{
					abond[ba++] = string[ab];
					ab++;
				}
				next_node = ab;
				abondance = atoi(abond);
			}
			else{
				next_node = suiv;
				abondance = 1;
			}
			temoin = 1;
			ancetre = VertexNumber;

			std::string str = "node" + std::to_string(no_name);
			char* int_node = new char[str.length() + 1];
    		strcpy(int_node, str.c_str());
			
			lesNoms[ancetre] = int_node;
			namesMap[int_node] = VertexNumber;
			abondMap[int_node] = abondance;

			VertexNumber++;
			no_name++;

		}

		else
		{
			while(string[suiv] != ':')
			{
				anc[i++] = string[suiv];
				suiv++;
			}
			anc[i++] = '\0';

			if(anc[0] == '0')
			{
				suiv++;
				i = 0;
				temoin = 2;
				while(string[suiv] != ')')
				{
					string1[i++] = string[suiv];
					suiv++;
				}
				string1[i++] = '\0';
				dist_root = atof(string1);
				ancetre = 0;
			}
			else {
				ancetre = atoi(anc);
			}
		}


		for ( ii = a1+1; ii <= a2; ii++)
		{
			if (string[ii] == ':')
			{
				xx = 0;
				a3 = ii+1;

				for(jj = zz; jj < ii; jj++)
				{
					string1[xx] = string[jj];
					xx++;
				}
				string1[xx++] = '\0';
				idSeq = atof(string1);

				if(string1[0] == NULL) {
					idSeq = VertexNumber;
					itoa_(VertexNumber, int_node, 10);
					lesNoms[VertexNumber] = int_node;
					VertexNumber++;
					
				}

			}

			else if((string[ii] == ',') || (string[ii] == ')'))
			{
				xx = 0;
				zz = ii +1;
				for ( jj = a3; jj < ii; jj++)
				{
					string1[xx++] = string[jj]; 
				}
				string1[xx++] = '\0';
				longueur = atof(string1);

				/*
				The ARETEB contains the nodes pair (child-ancestor) in order of processing. This organisation will make it easier to build the distance matrix later.
				It also requires less juggling than the previous function.
				*/
				ARETEB[numero][0] = idSeq;
				ARETEB[numero][1] = ancetre;
				LONGUEUR[numero] = longueur + dist_root;
				numero++;

			}
		}

		// fin for pour traiter noeud
		//transcrire la nouvelle chaine

		// Case in which we reached the root node
		if(temoin == 2)
		{
			string[0] = '0';
		}

		xx = 0;
		for ( jj = 0; jj < (int)a1; jj++)
		{ string2[xx++] = string[jj]; }

		/*
		Since many internal nodes are already named in lineage trees there is not always a need to add a name for the ancestor node.
		We only add it if the node was not named before, otherwise the ancestor node will be named as in the sequence.
		*/
		if(temoin == 1)
		{
			itoa_(ancetre,string1,10);   // 
			for( jj = 0; jj < (int) strlen(string1); jj++)
			{ string2[xx++] = string1[jj]; }
		}
		
		// transcrire la fin
		for( jj = next_node; jj <= a; jj++)
		{
			string2[xx++] = string[jj];
		}		
		tempString = string;
		string = string2;
		string2 = tempString;
		tempString = 0;

		a = xx;  // mettre la longueur à jour

	}


	/*printf("\nBranches list :");
	for(i = 0; i <= numero-1; i++){
		printf("\n%ld-%ld --> %lf",ARETEB[i][0],ARETEB[i][1],LONGUEUR[i]);
	}
    printf("\nEnd of list\n");*/

	free(string);
	free(string1);
	free(string2);
	free(tempString);

	(*kt) = 2*n-3 - (*kt);
    
	return dist_root;

}

/*/=================================================================================================
// Function using all the previous functions to build a distance matrix and write it in the output file
//=================================================================================================*/
void newickToMatrixLineage(const char *newick,FILE *out, std::map <std::string, int>& dicNames, std::map <std::string, int>& dicAbond, double **&ADD, double **&Adjacence, char* treeID){
	int i,j;
	int pos_racine = -1;
	int dist_naive;
	long int ** ARETEBcell;
	double *  LONGUEUR;
	char ** NAMES;
	double ** tempDIST;
	int kt;
	int nbNodes = 0;
	int size = nbNodesNewick(newick, nbNodes);
	int fullSize = size + nbNodes;

	char * newString = (char*)malloc(100000*sizeof(char));
	
	ADD = (double**)malloc(2*fullSize*sizeof(double*));
	Adjacence = (double**)malloc(2*fullSize*sizeof(double*));
	tempDIST = (double**)malloc(2*fullSize*sizeof(double*));
	for(i=0;i<=2*fullSize;i++){
		ADD[i] = (double*)malloc(2*fullSize*sizeof(double));
		Adjacence[i] = (double*)malloc(2*fullSize*sizeof(double));
		tempDIST[i] = (double*)malloc(2*fullSize*sizeof(double));
	}

	LONGUEUR = (double*)malloc((4*(size))*sizeof(double));
	ARETEBcell = new long int*[fullSize];
	for(i = 0; i < fullSize; i++)
	{
		ARETEBcell[i] = new long int[2];
		ARETEBcell[i][0] = 0;
		ARETEBcell[i][1] = 0;
	}

	NAMES=(char**)malloc(2*size*sizeof(char*));
	for(i=0;i<=size;i++)
	{	NAMES[i] = (char*)malloc(50); }

	//Method to retrieve the nodes' names and store them in a list
	getNamesNaive(newick, NAMES, newString, size, dicAbond, dicNames);
	//printf("\n Nouvelle séquence Newick: %s \n", newString);
	
	// Retrieve the distances between nodes then build the distance matrix
	dist_naive = lectureNewickBcell(newString,ARETEBcell,LONGUEUR,NAMES,&kt, size, dicNames, dicAbond);
    loadAdjacenceMatrixLineage(tempDIST, Adjacence, ARETEBcell, LONGUEUR, fullSize-1, kt);
    FloydLineage(tempDIST, ADD, fullSize-1, kt, dist_naive); 



	int max_taille=0;
	for(j = 1; j <= size;j++){
		max_taille = (strlen(NAMES[j])>max_taille)?strlen(NAMES[j]):max_taille;
	}

	// Write the name of the tree, number of nodes and distance matrix in the output file
	fprintf(out,"%s \t number of nodes :%d",treeID,size);
	for(i = 0; i < fullSize; i++){
		fprintf(out,"\n%s",NAMES[i]);
		//printf("\n %s", NAMES[i]);
		if(strlen(NAMES[i])<max_taille){
			for(j = strlen(NAMES[i]); j <= max_taille; j++){
				fprintf(out," ");
			}
		}
		for(j = 0; j < fullSize; j++){
			fprintf(out,"  %3.5lf",ADD[i][j]);
			//printf("\t %f", Adjacence[i][j]);
		}
	}
	fprintf(out, "\n\n\n");
	//printf("\n\n");

	free(tempDIST);
	delete[] ARETEBcell;
	free(LONGUEUR);
	free(newString);

}