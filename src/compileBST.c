/*! \file compileBST.c
 *  \brief	   This implements the applyPatch program.
 *  \author    Lucie Pansart
 *  \author    Jean-Louis Roch
 *  \version   1.0
 *  \date      30/9/2016
 *  \warning   Usage: compileBST n originalFile 
 *  \copyright GNU Public License.
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include <assert.h>
#include <string.h>

// Déclaration du nombre d'élèment dans le dictionnaire en variable globale
// pour gérer la fin de l'écriture de l'arbre : il ne faut pas mettre de , 
// sur le dernier couple

int valeur_max; 

// Lecture du fichier et calcul des probabilités d'apparition de chaque élément
// Le résultat est stocké dans un tableau de flotant


float * lecture_fichier(int n, FILE *F){
	int i;
	int S = 0;
	int * res = malloc(n*sizeof(int));
	float * res2 = malloc(n*sizeof(float));
	for (i=0;i<n;i++){
		fscanf(F, "%d", res+i);
		S=S+res[i];
	}

	for (i=0;i<n;i++){
		*(res2+i)=(float)*(res+i)/S;
	}
	free(res);
	return res2;
}

// Calcul des sommes partielles de probabilités
// On stocke ces sommes dans une matrice dont la case (i,j)
// contient la somme partielle des probabilités des éléments i à j.


float ** somme(int n, float *f){
	int i;  
	int j;  
	float ** res = (float**)malloc(n*sizeof(float*));
	for (i = 0; i < n; i++) 
	{
		res[i] = (float*)calloc(n, sizeof(float));
	}
	for (i = 0; i<n; i++){
		res[i][i] = f[i];
		for (j= i+1; j<n; j++){ 
			res[i][j]= res[i][j-1] + f[j];
		}
	}
	free(f);
	return res;
}


// Fonctions de libération des allocations de matrices

void libere_matrice_flottant(  float** M, int taille ) {
	int i ;
	for ( i =0; i < taille; i++) {
		free(M[i]);
		M[i] = NULL;
	}
	free(M);
	M = NULL;
} 
/*
void libere_matrice_entier(  int** M, int taille ) {
	int i ;
	for ( i =0; i < taille; i++) {
		free(M[i]);
		M[i] = NULL;
	}
	free(M);
} 

*/

// Calcul de l'ABR optimal en créant deux matrices cout et racine que l'on instancie
// sur les valeurs possibles ( racine de l'ABRopt ne possédant qu'un seul élèment, ..)
// puis calcul progressif en utilisant la formule de Bellman obtenue dans la première 
// partie du projet
// En sortie la matrice cout contient les couts optimaux pour les ABR contenant les éléments
// entre i et j et la matrice racine contient les racines de ces ABR
// Renvoit la matrice racine



int ** ABRopt(float ** S, int n){
	int i,j,k,s;
	float t;
	float ** cout = (float**)malloc((n+1)*sizeof(float*));
	int ** r= (int**)malloc(n*sizeof(int*));
	for (i = 0; i < n; i++)
	{
		r[i]=(int*)calloc(n,sizeof(int));
		r[i][i]=i;
	}
	for (i=0;i<n+1;i++){
		cout[i] = (float*)calloc((n+1), sizeof(float));
	}
	for(i=0; i<n; i++){
		cout[i][i] = 0;
		cout[i][i+1] = S[i][i];
	}
	for( i = 0; i < n; i++){
		for(j=i+2;j<n+1;j++){
			cout[i][j]=(float)n;
		}
	}

	for (k=1;k<n;k++){
		for (i=0;i<n-k;i++){
			j=i+k+1;
			for (s=i+1;s< j+1;s++){
				t=cout[i][s-1]+cout[s][j]+S[i][j-1];
				if (t<cout[i][j]){
					cout[i][j]=t;
					r[i][j-1]=s-1;
				}
			}
		}
	}
	libere_matrice_flottant( cout, n + 1);
	return r;
}


// Affichage de l'ABR par récursivité

void ABR(int ** racine, int k1, int k2){
	// Si l'ABR ne contient en fait qu'un élèment
	if (k1==k2){
		if ( racine[k1][k2] == valeur_max-1 ) { // alors dernier élèment à écrire donc pas de ","
			printf(" {-1, -1}");
		} else {
			printf(" {-1, -1},");
		}

		// si la racine est l'élèment tout à gauche
	} else if(racine[k1][k2]==k1){
		printf(" {-1, %d},",racine[k1+1][k2]);
		ABR(racine,k1+1,k2);

		//Si la racine est l'élément tout à droite
	} else if ( racine[k1][k2] == k2 ){
		if ( racine[k1][k2] == valeur_max-1 ) { // dernier élément : pas de ","
			(ABR(racine,k1,k2-1));
			printf(" {%d, -1}",racine[k1][k2-1]);
		} else {
			(ABR(racine,k1,k2-1));
			printf(" {%d, -1},",racine[k1][k2-1]);
		}

		// la racine a des éléments à gauche et à droite
	} else {
		ABR(racine,k1,racine[k1][k2]-1);
		printf(" {%d, %d},",racine[k1][racine[k1][k2]-1],racine[racine[k1+1][k2]][k2]);              
		ABR(racine,racine[k1][k2]+1,k2);
	}
}


/**
 * Main function
 * \brief Main function
 * \param argc  A count of the number of command-line arguments
 * \param argv  An argument vector of the command-line arguments.
 * \warning Must be called with a positive long integer, n,  and a filename, freqFile, as commandline parameters and in the given order.
 * \returns { 0 if succeeds and prints C code implementing an optimal ABR on stdout; exit code otherwise}
 */
int main (int argc, char *argv[]) {
	char *useless_ptr;
	// conversion de l'argument en entrée en entier
	valeur_max = (int) strtol(argv[1], &useless_ptr, 10);
	long n = 0 ; // Number of elements in the dictionary
	FILE *freqFile = NULL ; // File that contains n positive integers defining the relative frequence of dictinary elements 

	if(argc != 3){
		fprintf(stderr, "!!!!! Usage: ./compileBST n  originalFile !!!!!\n");
		exit(EXIT_FAILURE); /* indicate failure.*/
	}

	{ // Conversion of parameter n in a long 
		int codeRetour = EXIT_SUCCESS;
		char *posError;
		long resuLong;
		n = atol(argv[1] ) ;

		assert(argc >= 2);
		// Conversion of argv[1] en long
		resuLong = strtol(argv[1], &posError, 10);
		// Traitement des erreurs
		switch (errno)
		{
			case EXIT_SUCCESS :
				// Conversion du long en int
				if (resuLong > 0)
				{
					n = (long)resuLong;
					fprintf(stderr, "Number of elements in the dictionary: %ld\n", n);         
				}
				else
				{
					(void)fprintf(stderr, "%s cannot be converted into a positive integer matching the number of elements in the dicionary.\n", argv[1]) ; 
					codeRetour = EXIT_FAILURE;
				}
				break;

			case EINVAL :
				perror(__func__);
				(void)fprintf(stderr, "%s does not match a long integer value. \n", argv[1]);
				codeRetour = EXIT_FAILURE;
				break;

			case ERANGE :
				perror(__func__);
				(void)fprintf(stderr, "%s does not fit in a long int :\n" "out of bound [%ld;%ld]\n",
						argv[1], LONG_MIN, LONG_MAX);
				codeRetour = EXIT_FAILURE;
				break;
			default :
				perror(__func__);
				codeRetour = EXIT_FAILURE;
		} // switch (errno)
		if  (codeRetour != EXIT_SUCCESS) return codeRetour ;
	}

	freqFile = fopen(argv[2] , "r" );
	if (freqFile==NULL) {fprintf (stderr, "!!!!! Error opening originalFile !!!!!\n"); exit(EXIT_FAILURE);}
	float ** S=somme(valeur_max, lecture_fichier(valeur_max, freqFile));
	int ** P=ABRopt(S, valeur_max);

	printf("static int BSTroot = %d;\n",*(P[0]+valeur_max - 1));
	printf("static int BSTtree[%d][2] = {",valeur_max);
	ABR(P, 0, valeur_max - 1);

	printf(" };\n");
	fclose(freqFile);

//	libere_matrice_flottant(  S, n);
//	libere_matrice_entier( P, n);


	return 0;
}
