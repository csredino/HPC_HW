#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
using namespace std;
int main() {
	//declarations 
	const double pi = 4 * atan(1.0);
	int M; //size of M by M grid
	double h = 1 / double(M + 1); //step size, set so that the actual size of box is 1 by 1
	vector<vector<double> > phi; // potential to be found, a vector of vectors makes a 2d array
	vector<vector<double> > phiNew; //new potential after each step
	int steps;//number of iterations
	double norm;//rms difference between iterations of the solution
	double normtot;//running sum of norms over all elements in the M x M array
	double avgnorm;//
	//inputs
	cout << "Enter size for M x M grid to be used in evaluation, M=";
	cin >> M;
	//intialization
	steps = 0;
	normtot = 0;//intialize normtot to be zero
	avgnorm = 1;//intialize to some value larger than that which is acceptable
	phi = phiNew = vector<vector<double> >(M, vector<double>(M, 1.0));//intialize arrays with ones
	//set boundary conditions
	for (int k = 0; k < M; k++){
		phi[k][0] = phiNew[k][0] = sin(pi*k*h);
		phi[k][M] = phiNew[k][M] = sin(pi*k*h)*exp(-pi);
		phi[0][k] = phiNew[0][k] = 0;
		phi[M - 1][k] = phiNew[M - 1][k] = 0;//don't know why the point [M][M] is segfaulting here, but it should be defined two lines above anyways
	}
	while (avgnorm > 0.00005){//main loop
		//jacobi iteration
		normtot = 0;//clear the total for each jacobi step
		for (int i = 1; i < M - 1; i++){//loop through x's, we only go to M-1 because M-1 is our boundary(fixed)
			for (int j = 1; j < M - 1; j++){//loop through y's
				phiNew[i][j] = 0.25 * (phi[i - 1][j] + phi[i + 1][j] + phi[i][j - 1] + phi[i][j + 1]);
				norm = sqrt((phiNew[i][j] - phi[i][j])*(phiNew[i][j] - phi[i][j]));
				normtot = normtot + norm;
			}//end of y loop
		}//end of x loop
		vector< vector<double> > swap = phi;
		phi = phiNew;
		phiNew = swap;//update old phi for next step, i've been told the "swap" array will help to avoid pointer errors,don't fully understand this
		avgnorm = normtot / (M*M);
		steps++;
	}//end of while loop
	printf("Calculation took %i jacobi steps to converge.\n", steps);
	return 0;
}


