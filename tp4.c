//TP4 : jeu de la vie

/* Bibliothèques */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */


//-------------------------------------------------------------------//
// Code de Makoto Mastumoto pour le Mersenne Twister                 //
//-------------------------------------------------------------------//

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */





//-------------------------------------------------------------------//
// Code du jeu de la vie                                             //
//-------------------------------------------------------------------//


//-------------------------------------------------------------------//
//creationGrille : création grille taille n*n et initialisation      //
//                 vide (avec des points)                            //
//                                                                   //
//En entrée : n : entier égal à la taille du tableau                 //
//                                                                   //
//En  sortie : tableau de caractères créé                            // 
//-------------------------------------------------------------------//

char ** creationGrille (int n)
{
	int     i, j, test = 1;
	char ** grille  = NULL ;
	
	grille = malloc ( n * sizeof (char *) );
	
	// test réussite allocation mémoire
	if ( grille != NULL )
	{
		for ( i = 0; i < n; i ++)
		{
			grille [i] = malloc ( n * sizeof (char) ) ;
			
			// test réussite allocation mémoire
			if ( grille [i] == NULL )
			{
				test = 0;
			}
		}
	}
	
	// en cas d'erreur d'allocation mémoire
	if ( grille == NULL || test == 0 )
	{
		printf ("Erreur allocation mémoire\n");
		//quitte le programme
		exit (0);
	}  
	
	// initialisation de la grille avec des points
	for (i = 0; i < n; i ++)
	{
		for (j = 0; j < n; j ++)
		{
			grille [i][j] = '.';
		}
	}
	
	return grille;
}


//-------------------------------------------------------------------//
//afficheGrille : affiche la grille de taille n*n                    //
//                                                                   //
//En entrée : grille : tableau 2D à afficher                         //
//            n : entier égal à la taille du tableau                 //
//                                                                   //
//En  sortie : void                                                  // 
//-------------------------------------------------------------------//

void afficheGrille (char ** grille, int n)
{
	int i, j;
	
	for (i = 0; i < n ; i ++)
	{
		for (j = 0; j < n; j ++)
		{
			printf ("%2c", grille [i][j]);
		}
		printf ("\n");
	}
	printf ("\n");
}


//-------------------------------------------------------------------//
//initGlider : initialise la grille de taille n*n avec un glider     //
//                                                                   //
//En entrée : grille : tableau 2D contenant la grille                //
//            n : entier égal à la taille du tableau                 //
//                                                                   //
//En  sortie : grille : tableau modifié                              // 
//-------------------------------------------------------------------//

char ** initGlider (char ** grille, int n)
{
	grille [1][1] = 'X';
	grille [2][2] = 'X';
	grille [3][0] = 'X';
	grille [3][1] = 'X';
	grille [3][2] = 'X';

	return grille;
}


//-------------------------------------------------------------------//
//initGalaxie : initialise la grille de taille 50² avec une galaxie  //
//                                                                   //
//En entrée : grille : tableau 2D contenant la grille                //
//            n : entier égal à la taille du tableau                 //
//                                                                   //
//En  sortie : grille : tableau modifié                              // 
//-------------------------------------------------------------------//

char ** initGalaxie (char ** grille, int n)
{
	if (n != 50)
	{
		printf ("Grille trop petite pour contenir une galaxie !\n");
		//quitte le programme
		exit (0);
	}
	
	int i, j;
	
	// initialisation de la galaxie
	for (i = 20; i < 22; i ++)
	{
		for (j = 20; j < 26; j++)
		{
			grille [i][j] = 'X';
		}
	}
	
	for (i = 27; i < 29; i ++)
	{
		for (j = 23; j < 29; j++)
		{
			grille [i][j] = 'X';
		}
	}
	
	for (i = 23; i < 29; i ++)
	{
		for (j = 20; j < 22; j++)
		{
			grille [i][j] = 'X';
		}
	}
	
	for (i = 20; i < 26; i ++)
	{
		for (j = 27; j < 29; j++)
		{
			grille [i][j] = 'X';
		}
	}
	
	return grille;
}


//-------------------------------------------------------------------//
//initPulsar : initialise la grille de taille 50² avec un pulsar     //
//                                                                   //
//En entrée : grille : tableau 2D contenant la grille                //
//            n : entier égal à la taille du tableau                 //
//                                                                   //
//En  sortie : grille : tableau modifié                              // 
//-------------------------------------------------------------------//

char ** initPulsar (char ** grille, int n)
{
	if (n != 50)
	{
		printf("Grille trop petite pour contenir un pulsar !\n");
		// quitte le programme
		exit (0);
	}
	
	int i , j = 20;
	
	// initialisation du pulsar
	for (i = 20; i < 25; i ++)
	{
		grille [i][j] = 'X';
	}
	
	j = 24;
	
	for (i = 20; i < 25; i ++)
	{
		grille [i][j] = 'X';
	}
	
	grille [20][22] = 'X';
	grille [24][22] = 'X';
	
	return grille;
}


//-------------------------------------------------------------------//
//initRandom : initialise la grille de taille n*n aléatoirement      //
//             avec une proportion choisie de cellules vivantes      //
//                                                                   //
//En entrée : grille : tableau 2D contenant la grille                //
//            n : entier égal à la taille du tableau                 //
//            prop : réel égal à la prop choisie de cellules         //
//            vivantes au début du jeu                               //
//                                                                   //
//En  sortie : grille : tableau modifié                              // 
//-------------------------------------------------------------------//

char ** initRandom (char ** grille, int n, float prop)
{
	double rdm;
	int    i, j;
	
	//initialisation aléatoire
	for (i = 0; i < n ; i ++)
	{
		for (j = 0; j < n; j ++)
		{
			rdm = genrand_real2();
			
			if (rdm < prop)
			{
				grille [i][j] = 'X';
			}
		}
	}
	
	return grille;
}


//-------------------------------------------------------------------//
//libereGrille : libère la mémoire de la grille de taille n*n        //
//                                                                   //
//En entrée : grille : tableau 2D contenant la grille                //
//            n : entier égal à la taille du tableau                 //
//                                                                   //
//En  sortie : void                                                  // 
//-------------------------------------------------------------------//

void libereGrille (char ** grille, int n)
{
	int i;
	
	if (grille != NULL)
	{
		for (i = 0; i < n; i ++)
		{
			free (grille[i]);
		}
		free (grille);
	}    
}


//-------------------------------------------------------------------//
//indiceValide : renvoie 1 si l'indice [k][l] de la cellule voisine  //
//               de la case considérée [i][j] est valide, 0 sinon    //
//               L'indice est valide si la case voisine de [i][j]    //
//               n'est pas [i][j] et si elle ne se trouve pas en     //
//               dehors de la grille                                 //
//               Cette fonction sert uniquement pour l'univers non   //
//               torique                                             //
//                                                                   //
//En entrée : grille : tableau 2D contenant la grille                //
//            n : entier égal à la taille du tableau                 //
//            i : ligne de la case considérée                        //
//            j : colonne de la case considérée                      //
//            k : ligne de la case voisine                           //
//            l : colonne de la case voisine                         //
//                                                                   //
//En  sortie : entier 1 ou 0                                         // 
//-------------------------------------------------------------------//

int indiceValide (char ** grille, int n, int i, int j, int k, int l)
{
	// cellule étudiée
	if ( k == i && l == j )
	{
		return 0;
	}
	// cellule en dehors de la grille
	else if ( k < 0 || l < 0 || k >= n || l >= n )
	{
		return 0;
	}
	else
	{
		return 1;
	}
}


//-------------------------------------------------------------------//
//nbVoisins : renvoie le nombre de voisins vivants de la case        //
//            d'indice [i][j] de la grille dans un univers           //
//            non torique                                            //
//                                                                   //
//En entrée : grille : tableau 2D contenant la grille                //
//            n : entier égal à la taille du tableau                 //
//            i : ligne de la case considérée                        //
//            j : colonne de la case considérée                      //
//                                                                   //
//En  sortie : nombre de voisins (entier)                            // 
//-------------------------------------------------------------------//

int nbVoisins (char ** grille, int n, int i, int j)
{
	int nv = 0, k, l, test = 0 ;
	
	for (k = (i - 1); k <= (i + 1); k ++)
	{ 
		for (l = (j - 1); l <= (j + 1); l ++)
		{
			test = indiceValide (grille, n, i, j, k , l);
			
			// si l'indice du voisin est valide et si le voisin est vivant
			if ( test  == 1 && grille [k][l] == 'X' )
			{
				nv ++;
			}
		}
	}
	
	return nv;
}


//-------------------------------------------------------------------//
//nbVoisinsTore  : renvoie le nombre de voisins vivants de la case   //
//                 d'indice [i][j] de la grille dans un univers      //
//                 torique                                           //
//                                                                   //
//En entrée : grille : tableau 2D contenant la grille                //
//            n : entier égal à la taille du tableau                 //
//            i : ligne de la case considérée                        //
//            j : colonne de la case considérée                      //
//                                                                   //
//En  sortie : nombre de voisins (entier)                            // 
//-------------------------------------------------------------------//

int nbVoisinsTore (char ** grille, int n, int i, int j)
{
	int nv = 0, k, l, kv, lv;
	
	for (k = (i - 1); k <= (i + 1); k ++)
	{
		for (l = (j - 1); l <= (j + 1); l ++)
		{
			kv = k;
			lv = l;
			
			// adaptation de l'indice dans l'univers torique
			if ( k < 0 )
			{
				kv = k + n;
			}
			
			if ( l < 0 )
			{
				lv = l + n;
			}
			
			if ( k >= n )
			{
				kv = k % n;
			}
			
			if ( l >= n )
			{
				lv = l % n;
			}
			
			// si l'indice du voisin n'est pas la cellule actuelle et si le voisin est vivant
			if ( ( kv != i || lv != j ) && grille [kv][lv] == 'X' )
			{
				nv ++;
			}
		}
	}
	
	return nv;
}


//-------------------------------------------------------------------//
//jeu : modifie la grille2 à partir de grille1 en suivant les règles //
//      du jeu de la vie                                             //
//                                                                   //
//En entrée : grille1 : tableau 2D contenant la grille à l'instant t //
//            grille2 :tableau 2D contenant la grille à l'instant t+1//
//            taille : entier égal à la taille du tableau            //
//            tore : entier valant 1 si on est dans un univers       //
//            torique et 0 dans le cas contraire                     //
//                                                                   //
//En  sortie : void (passage par référence)                          // 
//-------------------------------------------------------------------//

void jeu (char ** grille1, char ** grille2, int taille, int tore)
{
	int i, j, nv = 0;
	
	for (i = 0; i < taille ; i ++)
	{
		for (j = 0; j < taille; j ++)
		{
			// univers non torique
			if (tore == 0)
			{
				nv = nbVoisins (grille1, taille, i, j);
			}
			// univers torique
			else
			{
				nv = nbVoisinsTore (grille1, taille, i, j);
			}
			
			// la cellule vivante survit si elle a 2 ou 3 voisins
			if ( grille1 [i][j] == 'X' && (nv == 2 || nv == 3) )
			{
				grille2 [i][j] = 'X';
			}
			// la cellule morte nait si elle a 3 voisins
			else if ( grille1 [i][j] == '.' &&  nv == 3 )
			{
				grille2 [i][j] = 'X';
			}
			// la cellule meurt ou reste morte dans les autre cas
			else
			{
				grille2 [i][j] = '.';
			}
		}
	}
}



//-------------------------------------------------------------------//
// programme principal                                               //
//-------------------------------------------------------------------//

int main (int argc, char *argv[])
{
	// initialisation du Mersenne Twister
	unsigned long init[4] = {0x123, 0x234, 0x345, 0x456}, length = 4;
	init_by_array(init, length);
	
	// déclaration des variables
	int     tailleGrille, choix, i, tore, debug = 1, iter;
	float   prop;
	char ** grille1;
	char ** grille2;
	
	
	printf ("Bienvenue dans le jeu de la vie !\n\n");
	
	printf ("Veuillez choisir une option pour le type d'univers :\n\n");
	printf ("1 : univers torique \n");
	printf ("2 : univers non torique \n");
	printf ("\nVotre choix : ");
	
	scanf ("%d", &choix);
	
	// choix du type d'univers
	switch (choix) 
	{
		case 1 :
			tore = 1;
		break;
		
		case 2 :
			tore = 0;
		break;

		default :
			printf ("Erreur de saisie\n");
			//quitte le programme
			exit (0);
		break;
	}
	
	printf ("\n");
	
	printf ("Veuillez choisir la taille de la grille de jeu :\n\n");
	printf ("1 : grille de taille 10*10 \n");
	printf ("2 : grille de taille 50*50 \n");
	printf ("\nVotre choix : ");
	
	scanf ("%d", &choix);
	
	// choix de la taille de la grille
	switch (choix) 
	{
		case 1 :
			tailleGrille = 10;
		break;
		
		case 2 :
			tailleGrille = 50;
		break;

		default :
			printf ("Erreur de saisie\n");
			//quitte le programme
			exit (0);
		break;
	}

	printf ("\n");
	
	// création de 2 grilles de taille tailleGrille
	grille1 = creationGrille (tailleGrille);
	grille2 = creationGrille (tailleGrille);
	
	
	
	printf ("Veuillez choisir l'initialisation de la grille de jeu :\n\n");
	printf ("1 : glider \n");
	printf ("2 : galaxie (seulement pour une grille 50*50)\n");
	printf ("3 : pulsar (seulement pour une grille 50*50)\n");
	printf ("4 : aléatoire \n");
	printf ("\nVotre choix : ");
	
	scanf ("%d", &choix);
	
	printf ("\n");
	
	// choix de l'initialisation
	switch (choix)
	{
		case 1 :
			initGlider (grille1, tailleGrille);
		break;
		
		case 2 :
			initGalaxie (grille1, tailleGrille);
		break;
		
		case 3 :
			initPulsar (grille1, tailleGrille);
		break;
		
		case 4 :
			printf ("Veuillez entrer la proportion de cellules vivantes au départ : ");
	
			scanf ("%f", &prop);
			
			if ( prop < 0 || prop > 1 )
			{
				printf ("Le nombre saisi doit etre compris entre 0 et 1 \n");
				//quitte le programme
				exit (0);
			}
			else
			{
				printf ("\n");
				initRandom (grille1, tailleGrille, prop);
			}
		break;

		default :
			printf ("Erreur de saisie\n");
			//quitte le programme
			exit (0);
		break;
	}
		
	// choix du nombre d'itérations
	printf ("Veuillez choisir le nombre d'itérations : ");
	
	scanf ("%d", &iter);
	
	
	printf ("\nDébut du jeu :\n");
	
	afficheGrille (grille1, tailleGrille);
	
	// pause de 1 seconde pour voir l'état initial
	usleep (1000000);
	
	for (i = 0; i < iter; i ++)
	{
		// alternance des grilles à t et t+1
		if ( i % 2 == 0)
		{
			jeu (grille1, grille2, tailleGrille, tore);
			afficheGrille (grille2, tailleGrille);
		}
		else
		{
			jeu (grille2, grille1, tailleGrille, tore);
			afficheGrille (grille1, tailleGrille);
		}
		
		// pause de 0.3 secondes en mode debug
		if (debug == 1)
		{
			usleep (300000);
		}
	}
	
	printf ("Fin du jeu après %d itérations\n", iter);
	

	// libération mémoire
	libereGrille (grille1, tailleGrille);
	libereGrille (grille2, tailleGrille);
	
	return 0;
}
