#include <stdio.h>
#include <stdlib.h> 
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "omp.h"
using namespace std;
int main() {
	//declarations
	const double pi = 4 * atan(1.0);
	int M; //size of M by M grid
	int j;
	int nthreads;
	double h = 1 / double(M + 1); //step size, set so that the actual size of box is 1 by 1
	vector<vector<double> > phi; // potential to be found, a vector of vectors makes a 2d array
	vector<vector<double> > phiNew; //new potential after each step
	int steps;//number of iterations
	double sqrdiff;//difference of squares between iteration
	double norm;//
	//inputs
	cout << "Enter size for M x M grid to be used in evaluation, M=";
	cin >> M;
	//intialization
	double tstart = omp_get_wtime();//begin timer
	steps = 0;
	norm = 1;//intialize to some value larger then that which is acceptable
	phi = phiNew = vector<vector<double> >(M, vector<double>(M, 1.0));//intialize arrays with ones
	//set boundary conditions
	for (int k = 0; k<M; k++){
		phi[k][0] = phiNew[k][0] = sin(pi*k*h);
		phi[k][M] = phiNew[k][M] = sin(pi*k*h)*exp(-pi);
		phi[0][k] = phiNew[0][k] = 0;
		phi[M - 1][k] = phiNew[M - 1][k] = 0;//don't know why the point [M][M] is segfaulting here, but it should be defined two lines above anyways
	}
	while (norm > 0.0001){//main loop
		//jacobi iteration
		sqrdiff = 0;//clear the total for each jacobi step
		for (int i = 1; i < M - 1; i++){//loop through x's, we only go to M-1 because M-1 is our boundary(fixed)
			#pragma omp parallel reduction(+:sqrdiff)//Dividing tasks up into columns, so value of i is shared (they can talk left and right)
			// j stays in one thread(this is ok since y is completely contained within the "forked" region).
			// sqrdiff will be reduced by summation, so I no longer need to define sqrdifftot as in the serial case
			{//start of parallel region
				#pragma omp for
				for (int j = 1; j < M - 1; j++){//loop through y's 
					phiNew[i][j] = 0.25 * (phi[i - 1][j] + phi[i + 1][j] + phi[i][j - 1] + phi[i][j + 1]);
					sqrdiff+= (phiNew[i][j]-phi[i][j])*(phiNew[i][j]-phi[i][j]);
				}//end of y loop
			}//end of parallel region 
		}//end of x loop 
		//loop to reset the old phi as the current phi
		#pragma omp parallel
		{
			#pragma omp for
			for (int i = 1; i < M - 1; i++){//loop through x's, we only go to M-1 because M-1 is our boundary(fixed)
				for (int j = 1; j < M - 1; j++){//loop through y's
					phi[i][j] = phiNew[i][j];
				}//end of x loop
			}//end of y loop 
		}
		norm = sqrt(sqrdiff); 
		steps++;
	}//end of while loop
	double tend = omp_get_wtime();
	printf("Calculation took %i jacobi steps to converge.\n", steps);
	printf("This took %f seconds \n ", (tend - tstart)); 
	return 0;
}












