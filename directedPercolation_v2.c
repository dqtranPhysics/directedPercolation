//This program uses the Domany-Kinzel algorithm to simulate bond directed percolation on lattice for 1+1d and 2+1d dimensions
//by Dat Tran, written 11/2022
//Input
//L - size of lattice
//tmax - number of total iterations
//d - dimension
//p - percolation
//seed - a large negative integer to generate random number
//Output
//Text files containing measurement of number of active sites, rms spreading, and survivability probability

#include <stdio.h>        // standard io library
#include <string.h>   // standard string library
#include <stdlib.h>       // standard library with lots of functions
#include <math.h>           // standard math library
#define NRANSI              // needed for NR
#include "my_nrutil.h"    // include NR header files
int *A;		//the lattice array
int *N_a;    //number of active sites
int *L_old;		//the integer list with coordinates of the wet sites
int * L_new;		// the integer list with coordinates of the neigboring sites to the wet sites
long L, seed;
double p;

int dkModel(int s_left, int s_right);
float ran2(long *idum);
void initA();
void plotA();
int countActiveSites();
//define the constructor for dynamic array
typedef struct {
  int *array;
  size_t used;
  size_t size;
} Array;
void initArray(Array *a, size_t initialSize);
void insertArray(Array *a, int element);
void freeArray(Array *a);

int main(int argc, char *argv[]) {
	int seed_pos;	
	long t, tmax;
    
    FILE *numberActiveSites;

        if (argc == 5){
                fprintf(stderr,"\n  ** Good Initialization **\n");
                L = atol(argv[1]);
                tmax = atol(argv[2]);
                p = atof(argv[3]);
                seed = atol(argv[4]);
                fprintf(stderr,"Lattice size: %ld\n",L);
                fprintf(stderr,"Number of iterations: %ld\n",tmax);
                fprintf(stderr,"Percolation probability: %f\n",p);
                fprintf(stderr,"Random seed: %ld\n\n",seed);
        }
        else {
                fprintf(stderr,"\n ** Initialization error **\n");
                fprintf(stderr,"Usage: ./directedPercolation.x L tmax p seed\n");  // correct input syntax
                return 1;
        }
	
	seed_pos = L/2;
    L_old = ivector(0,L);
    L_old[0] = seed_pos;
    //Array L_old;
    //initArray(&L_old,100);
    //insertArray(&L_old, seed_pos);
    N_a = ivector(0,tmax); //number of active sites

	//////////////////////////////////
	//Perform MC move over time tmax//
	//////////////////////////////////

    //Initialize A with zeros except for the first iteration
	for (int t = 1; t <= tmax; t++) {
        A = ivector(0,L-1);
        if (t == 1) {
            for (int i = 0; i < L; i++) {
                if (i == seed_pos) {
                    A[i] = 1;
                }
                else A[i] = 0;
            }
        }
        else {
            for (int i = 0; i < L; i++) {
                A[i] = 0;
            }
        }
        L_new = ivector(0,L);
        for (int i = 0; i < L; i++) L_new[i] = 0;
        //Array L_new;
        //initArray(&L_new,100);
        
		 //Update lattice A and L_new with L_old
        int j = 0;
        for (int i = 0; i < L; i++) {
		    int temp_val = L_old[i];
            if (temp_val == 0 || temp_val > L || temp_val < 0) break;
            //printf("%d\n",temp_val);
		    
		    //Updating array A with L_old, making sure to make the previous wet sites unwet.
            A[temp_val] = 1;
		    //Updating L_new
		    if (temp_val > 0 && temp_val < L - 1 && j > 1) {
                //L_new[j] = temp_val - 1;
                L_new[j] = temp_val + 1;
                j = j+1;
		    }
            else if (temp_val > 0 && temp_val < L - 1 && j == 0) {
                 L_new[j] = temp_val - 1;
                 L_new[j+1] = temp_val + 1;
                 j = j+2;
             }
		    //else if (temp_val == 0)  insertArray(&L_new, temp_val + 1);
            else if (temp_val == 0) {
                L_new[j] = temp_val + 1;
                j = j + 1;
            }
		    //else if (temp_val == L - 1)  insertArray(&L_new, temp_val - 1);
            else if (temp_val == L - 1)  {
                L_new[j] = temp_val - 1;
                j = j + 1;
            }
		 }
        //fprintf(stderr,"MC Move number %d \n", t);

		//free L_old
        //freeArray(&L_old);
        //free_ivector(L_old,0,L);
        //L_old = ivector(0,L);
        for (int i = 0; i < L; i++) L_old[i] = 0;

        //print L_new for debugging
        /*fprintf(stderr,"L_new: ");
        for (int i = 0; i < L; i++) {fprintf(stderr,"%d ",L_new[i]); };
        printf("\n");*/
        
		//Use DK model to check whether L_new sites are wetted?
        //Here we will iterate over the sites in L_new (not wet) to determine whether to wet them
        
        j = 0;
        for (int i = 0; i < L; i++) {
            int temp_val = L_new[i];
            int is_wet;
            if (temp_val == 0) break;
            is_wet = dkModel(A[temp_val - 1], A[temp_val + 1]);
            //Write wet sites to L_old
            if (is_wet == 1) {
                //numSites = numSites + 1;
                //insertArray(&L_old,temp_val);
                L_old[j] = temp_val;
                j = j + 1;
            }
        }

        //prtint L_old for debugging
        /*fprintf(stderr,"L_old: ");
        for (int i = 0; i < L; i++) {fprintf(stderr,"%d ",L_old[i]); };
        printf("\n");*/
        
        //Draw the lattice
        plotA();
        
        //reset array A
        free_ivector(A,0,L-1);
        free_ivector(L_new,0,L);
        
        //record number of active sites
        N_a[t-1] = countActiveSites();

	}
    free_ivector(L_old,0,L);
    
    //Write observables to file
    char head[150] = "numberOfActiveSites";
    char tail[50];
    sprintf(tail,"-L%ld-p%.2f.txt",L,p);
    numberActiveSites = fopen(strcat(head,tail),"w");
    for (int t = 0; t < tmax; t++) {
        fprintf(numberActiveSites,"%d\n",N_a[t]);
    }
    
    fclose(numberActiveSites);
    free_ivector(N_a,0,tmax-1);
    
}

int dkModel(int s_left,int s_right) {
	int s_now;
    double z = ran2(&seed);
    double p1 = p;
    double p2 = 1 - (1-p)*(1-p);
    if (s_left != s_right && z < p1) s_now = 1;
	else if (s_left == 1 && s_right == 1 && z < p2) s_now = 1;
	else s_now = 0;

	return s_now;
}

void plotA() {
	for (int i = 0; i < L; i++) {
        if (A[i] == 1) {
            fprintf(stderr,"\u2022 ");
        }
		else fprintf(stderr,"- ");
	}
	fprintf(stderr,"\n");

	return;
}

int countActiveSites() {
    int num = 0;
    for (int i = 0; i < L; i++) {
        if(A[i] == 1)
            num++;
    }
    return num;
}
//below is the functions that enable dyanmical arrays since such things don't exist in any libraries in C.

void initArray(Array *a, size_t initialSize) {
  a->array = malloc(initialSize * sizeof(int));
  a->used = 0;
  a->size = initialSize;
}

void insertArray(Array *a, int element) {
  // a->used is the number of used entries, because a->array[a->used++] updates a->used only *after* the array has been accessed.
  // Therefore a->used can go up to a->size 
  if (a->used == a->size) {
    a->size *= 2;
    a->array = realloc(a->array, a->size * sizeof(int));
  }
  a->array[a->used++] = element;
}

void freeArray(Array *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}

//Example of using them
/*
Array a;
int i;

initArray(&a, 5);  // initially 5 elements
for (i = 0; i < 100; i++)
  insertArray(&a, i);  // automatically resizes as necessary
printf("%d\n", a.array[9]);  // print 10th element
printf("%d\n", a.used);  // print number of elements
freeArray(&a);
*/

// below is a NR random number generator. It generated float numbers evenly over range (0,1)
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        float temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
