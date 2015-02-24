/***
	The Ising Model
	- 	T Mathialakan
		PHY 905
		ECE, MSU.
***/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>	

/*
	The energy between a pair of spin
*/
int pair(bool a, bool b){ // XOR 
	return (a^b) ? 2:-2; // -J if parallel, J if anti parallel;
}

/*
	Total magnetization in a Lattice
*/
int magnetization(bool L[], int n){
	int sum=0;
	 for(int i=0; i<n-1; i++)
		 sum=(int)L[i];
	return 2*sum - n;
}

/*
	The energy by a single spin formed with neighbours
*/
int spotvalue(bool L[], int k, bool spin, int c, int n){
	int val=0;
	if(k%c > 0) val += pair(L[k-1], spin ); // Left neighbour
	if(k%c < c-1) val += pair(L[k+1], spin ); // Right neighbour
	if(k>c) val += pair(L[k-c], spin ); // Above neighbour
	if(k<n-c) val += pair(L[k+c], spin); // Below neighbour
	return val;
}

/*
	The energy by a single spin formed with neighbours in periodic boundary condition
*/
int spotvalue_PBC(bool L[], int k, bool spin, int c, int n){
	int val=0, r=c;
	if(k%c > 0) val += pair(L[k-1], spin );  else val += pair(L[k+c-1], spin ); // Left neighbour
	if(k%c < c-1) val += pair(L[k+1], spin ); else val += pair(L[k-c+1], spin ); // Right neighbour
	if(k>c) val += pair(L[k-c], spin ); else val += pair(L[k+c*(r-1)], spin ); // Above neighbour
	if(k<n-c) val += pair(L[k+c], spin); else val += pair(L[k-c*(r-1)], spin ); // Below neighbour
	return val;
}

/*
	Print the 1D matrix
*/
void print_matrix_1D(bool *a, int m, int n){
	int i, j;
	printf("------------------------\n");
	for(i = 0; i< m; i++){
		for(j = 0; j< n; j++)
			printf("%d\t", a[i*m +j]);
		printf("\n");
	}
}

int main(int argc, char** arg){
	printf("Ising Model \n");
	int r =10, c=10, k=0;
	int n = r*c; // number of spins          
	bool *L = (bool *) malloc(sizeof(bool)*n); // Lattice
	float T = 0, THigh=5.0, interval=0.1; //Temperature from T to THigh with interval interval
	int samples =(int)((THigh-T)/interval)+1 ; //Sample size
	float R[3*samples]; // Result array set of temperature, average energy, and magnetization
	float delE=0.0; //newE - oldE
	float delR=0.0; //the random number between 0 and 1
	float  delMag = 0.0; //newM - oldM 
	int j=0;
	bool spin;
	const clock_t begin_time = clock();
	while(T<=THigh){
		// initialize all to  up
		for(int i=0; i<n; i++)
			L[i]=true;
		//print_matrix_1D(L, r, c);
		float beta = 1.0/T;
		float newE =0.0; // New energy
		float newMag =0.0; // New magnetization
		float curMag = n*1.0;
		float avgMag = 1.0; // initialize by the value when all are up
		float curE = -2*n + (r+c);//2*r*c -(r+c) - initialize by the value when all are up, - 2*n if Periodic BC
		float avgE = curE/n; 
		int i=0,
		it=10000; //Number of iterations
		while(i<it){
			k = rand() % (n-1); //index of a random spin
			delE = spotvalue_PBC(L, k, !L[k], r, n);//  Change in Energy around the neighbours 
			if (delE<=0){
				L[k] = !L[k]; //Accept the change
				newE = curE+delE;
				curE = newE;
				delMag = (L[k]) ? 2:-2;
				newMag = curMag+ delMag;
				curMag = newMag;
				}
			else {
					delR = (float)rand()/RAND_MAX;
					if(exp(-delE/T)>=delR){
						L[k] = !L[k]; //Accept the change
						newE = curE+delE;
						curE = newE;
						delMag = (L[k]) ? 2:-2;
						newMag = curMag+ delMag;
						curMag = newMag;
					}else
					{
						newE = curE;
						newMag = curMag;
					}
				} 
			avgE += (float) newE/n;
			avgMag += (float) newMag/n;
			i++;
		}
		R[j++]=T;
		R[j++]= avgE/(it+1);
		R[j++]=avgMag/(it+1);
		T+=interval;
	}
	printf( "time taken   %f ", (double)( clock () - begin_time ) / CLOCKS_PER_SEC );
	int i=0;
	printf("Temperature\t Energy\t Magnetization\n");
	while(i<3*samples){
		for(j = 0; j< 3; j++)
			printf("%f\t", R[i+j]);
		printf("\n");
		i+=j;
		}
	return 0;
}
