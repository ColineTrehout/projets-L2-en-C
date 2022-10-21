//TP 3 

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */
#define NB_CLASSE 6


//valeur de t de la loi de Student pour alpha = 0.05 et n=30
#define T 2.042


//-------------------------------------------------------------------//
// Code de Makoto Mastumoto pour le Mersenne Twister                 //
//-------------------------------------------------------------------//

static unsigned long mt[N]; /* the array for the state vector  */
static int mti= N + 1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti < N; mti ++) {
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
        if (i>=N) { mt[0] = mt[N-1]; i = 1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i >= N) { mt[0] = mt[N-1]; i=1; }
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

        for (kk=0; kk < N-M; kk ++) {
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


/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0 / 4294967295.0); 
    /* divided by 2^32-1 */ 
}



//-------------------------------------------------------------------//
//simPi : calcule une valeur approchée de Pi par une méthode de      //
//        Monte Carlo                                                //
//                                                                   //
//En entrée : nbpoints : nombre de points aléatoires à tirer         //
//                                                                   //
//En  sortie : valeur de Pi calculée (double)                        // 
//-------------------------------------------------------------------//

double simPi (long nbpoints)
{
	long inDisk = 0;
	double tiragex;
	double tiragey;
	long i;
	
	for ( i = 0; i < nbpoints; i ++)
	{
		tiragex = genrand_real1();
		tiragey = genrand_real1();
		
		if ( (tiragex * tiragex + tiragey * tiragey) <= 1) 
		{
			inDisk ++;
		}
	}
	
	return ( (double ) inDisk / (double) nbpoints * 4);
}



//-------------------------------------------------------------------//
//calculMoyenne : remplit un tableau avec les valeurs de Pi calculées//
//                par la fonction simPi et calcule la moyenne de     //
//                ces valeurs                                        //
//                                                                   //
//En entrée : nbexp : nombre d'expériences à réaliser                //
//            nbpoints : nombre de points aléatoires à chaque        //
//            expérience                                             //
//            valpi : tableau qui contiendra les valeurs de Pi       //
//            moypi: tableau qui contiendra la moyenne calculée de Pi//
//            indiceMoy : indice où stocker la moyenne pour la       //
//            simulation en cours                                    //
//                                                                   //
//En  sortie : void (passage par référence)                          // 
//-------------------------------------------------------------------//

void calculMoyenne (int nbexp, long nbpoints, double * valpi, double * moypi, int indiceMoy)
{
	long i;
	double pi;
	
	for ( i = 0; i < nbexp; i ++)
	{
		pi = simPi( nbpoints );
		valpi[i] = pi;
		moypi[indiceMoy] = moypi[indiceMoy] + pi;
	}
	moypi[indiceMoy] = moypi[indiceMoy] / nbexp;
}



//-------------------------------------------------------------------//
//calculCarresEcarts : remplit un tableau avec les valeurs des       //
//                     carrés des écarts entre les valeurs calculées //
//                     de Pi et la moyenne calculée                  //
//                                                                   //
//En entrée : nbexp : nombre d'expériences à réaliser                //
//            sommeEcartsCarres : tableau qui contiendra les valeurs //
//            des carrés des écarts à la moyenne calculée            //
//            valpi : tableau contenant les valeurs de pi calculées  //
//            moypi : tableau contenant la moyenne calculée          //
//            indiceMoy : indice où est stocké la moyenne pour la    //
//            simulation en cours                                    //
//                                                                   //
//En  sortie : void (passage par référence)                          // 
//-------------------------------------------------------------------//

void calculCarresEcarts (int nbexp, double * sommeEcartsCarres, double * valpi, double * moypi, int indiceMoy)
{
	int i;
	for ( i = 0; i < 30; i ++)
	{
		sommeEcartsCarres[i] = (valpi[i] - moypi[indiceMoy]) * (valpi[i] - moypi[indiceMoy]);
	}
}



//-------------------------------------------------------------------//
//calculVarianceEstimee : calcul de l'estimation de la variance      //
//                        stockée dans une case d'un tableau         //
//                                                                   //
//En entrée : nbexp : nombre d'expériences à réaliser                //
//            S2 : tableau qui contiendra les valeurs de             //
//            l'estimation sans biais de la variance S²              //
//            indiceS2 : indice où stocker la valeur de              //
//            S² pour la simulation en cours                         //
//            sommeEcartsCarres : tableau contenant les valeurs des  //
//            carrés des écarts à la moyenne calculée                //
//                                                                   //
//En  sortie : void (passage par référence)                          // 
//-------------------------------------------------------------------//

void calculVarianceEstimee (int nbexp, double * S2, int indiceS2, double * sommeEcartsCarres)
{
	int i;
	for ( i = 0; i < 30; i ++)
	{
		S2[indiceS2] =  S2[indiceS2] + sommeEcartsCarres[i];
	}
	S2[indiceS2] = S2[indiceS2] / (nbexp - 1);
}



//-------------------------------------------------------------------//
//calculR : calcul du rayon d'erreur au seuil de 95%, stocké dans    //
//          une case d'un tableau                                    //
//                                                                   //
//En entrée : nbexp : nombre d'expériences à réaliser                //
//            S2 : tableau contenant les valeurs de S²               //
//            indiceS2 : indice où est stocké la valeur              //
//            de S² pour la simulation en cours                      //
//            rayon : tableau qui contiendra le rayon d'erreur R     //
//            indicerayon : indice où doit etre stocké R             //
//                                                                   //
//En  sortie : void (passage par référence)                          // 
//-------------------------------------------------------------------//

void calculR (int nbexp,  double * S2, int indiceS2, double * rayon, int indicerayon)
{
	rayon[indicerayon] = T * sqrt ( S2[indiceS2] / nbexp );
}



//-------------------------------------------------------------------//
// programme principal                                               //
//-------------------------------------------------------------------//

int main(void)
{
	/*Initialisation du Mersenne Twister*/
	
    unsigned long	init[4] = {0x123, 0x234, 0x345, 0x456}, length=4;
    init_by_array(init, length);
	

	//Question 1
	
	double pi;
	
	/*Calcul et affichage valeur approchée de Pi et écart à la valeur théorique*/
	
	pi = simPi(1000); 
	printf("mille tirages aléatoires\npi = %.8lf, écart : %e\n\n", pi, pi-M_PI);
	
	pi = simPi(1000000); 
	printf("1 million tirages aléatoires\npi = %.8lf, écart : %e\n\n", pi, pi-M_PI); 
	
	pi = simPi(100000000); 
	printf("100 millions tirages aléatoires\npi = %.8lf, écart : %e\n\n", pi, pi-M_PI);
	
	pi = simPi(1000000000);
	printf("1 milliard tirages de aléatoires\npi = %.8lf, écart : %e\n\n", pi, pi-M_PI);
	
	/*Affichage valeur théorique de Pi*/
	
	printf("Valeur de pi contenue dans la librairie math.h\n%.8lf\n\n", M_PI);
	
	
	
	//Question 2
	
	/*Calcul des valeurs de Pi, moyenne et écart à la valeur théorique*/
	
	double * valpi0 = calloc ( 30, sizeof (double) );
	double * valpi1 = calloc ( 30, sizeof (double) );
	double * valpi2 = calloc ( 30, sizeof (double) );
	
	double * moypi = calloc ( 3, sizeof (double) );
	
	calculMoyenne( 30, 1000, valpi0, moypi, 0);
	printf("Moyenne pi sur mille tirages = %.8lf, écart : %e\n", moypi[0], moypi[0]-M_PI);
	
	calculMoyenne( 30, 1000000, valpi1, moypi, 1);
	printf("Moyenne pi sur 1 million tirages = %.8lf, écart : %e\n", moypi[1],moypi[1]-M_PI);
	
	//calculMoyenne( 30, 10000000000, valpi2, moypi, 2);
	printf("Moyenne pi sur 10 milliards tirages = %.8lf, écart : %e\n", moypi[2], moypi[2]-M_PI);
	
	//Question 3
	
	printf("\n");
	
	double * sommeEcartsCarres0 = calloc ( 30, sizeof (double) );
	double * sommeEcartsCarres1 = calloc ( 30, sizeof (double) );
	double * sommeEcartsCarres2 = calloc ( 30, sizeof (double) );
	
	double * S2 = calloc ( 3, sizeof (double) );
	
	double * rayon = calloc (3, sizeof (double) );
	

	/*Calcul des carrés des écarts à la moyenne*/
	
	//mille  tirages aléatoires
	calculCarresEcarts (30, sommeEcartsCarres0, valpi0, moypi, 0);
	
	//1 million tirages aléatoires
	calculCarresEcarts (30, sommeEcartsCarres1, valpi1, moypi, 1);
	
	//10 milliards tirages aléatoires
	calculCarresEcarts (30, sommeEcartsCarres2, valpi2, moypi, 2);
	
	/*Calcul et affichage de l'estimation sans biais de la variance S²*/
	
	calculVarianceEstimee ( 30, S2, 0, sommeEcartsCarres0);
	printf("S²(n) pour mille tirages = %e\n", S2[0]);
	
	calculVarianceEstimee ( 30, S2, 1, sommeEcartsCarres1);
	printf("S²(n) pour 1 million tirages = %e\n", S2[1]);
	
	calculVarianceEstimee ( 30, S2, 2, sommeEcartsCarres2);
	printf("S²(n) pour 10 milliards tirages = %e\n", S2[2]);
	
	printf("\n");
	
	/*Calcul et affichage du rayon d'erreur au seuil de 95%*/
	
	calculR ( 30,  S2, 0, rayon, 0);
	printf("Rayon d'erreur pour mille tirages R = %e\n", rayon[0]);
	
	calculR ( 30,  S2, 1, rayon, 1);
	printf("Rayon d'erreur pour 1 million tirages R = %e\n", rayon[1]);
	
	calculR ( 30,  S2, 2, rayon, 2);
	printf("Rayon d'erreur pour 10 milliards tirages R = %e\n", rayon[2]);
	
	printf("\n");
	
	/*Affichage des intervalles de confiance à 95%*/
	
	printf("Intervalle de confiance pour mille tirages [%.8lf; %.8lf]\n", moypi[0] - rayon[0], moypi[0] + rayon[0]);
	
	printf("Intervalle de confiance pour 1 million tirages [%.8lf; %.8lf]\n", moypi[1] - rayon[1], moypi[1] + rayon[1]);
	
	printf("Intervalle de confiance pour 10 milliards tirages [%.8lf; %.8lf]\n", moypi[2] - rayon[2], moypi[2] + rayon[2]);
	
	
	/*Libération mémoire*/
	free (valpi0);
	free (valpi1);
	free (valpi2);
	free (moypi);
	free (sommeEcartsCarres0);
	free (sommeEcartsCarres1);
	free (sommeEcartsCarres2);
	free (S2);
	free (rayon);
	
	return 1;
}
