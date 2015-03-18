#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
using namespace std;
int main(int argc, char*argv[]) {
	//declarations
	double totaltime;
	double runningtotal = 0;//intialize total time
	const double pi = 4 * atan(1.0);
	int M = 600; //size of M by M grid, hardcoded this time, so I don't have to specify which process is     
	//taking the input	
	int my_id;
	int num_procs;
	int elem_name;
	int j_first;//first and last indices for loops need to be variables, since first and last processors 
	//(bottom and top) will not be updating the bottom and top boundaries
	int j_last;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	double  h = 1 / double(M + 1); //step size, set so that the actual size of box is 1 by 1                 
	int steps;//number of iterations
	double sqrdiff;//difference of squares between iteration
	double norm;//the L2 norm that is actually tested for convergence
	MPI_Status status;

	//intialization
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Get_processor_name(processor_name, &elem_name);
	double phi[int(M / num_procs) + 2][M];//2 extra "ghost points", these overlap with other arrays for 
	//communications
	double phiNew[int(M / num_procs) + 1][M];//the updated vector only needs one overlap point
	double tstart = MPI_Wtime();;//begin timer
	steps = 0;
	norm = 1;//intialize to some value larger then that which is acceptable
	j_first = 1;
	j_last = int(M / num_procs);
	if (my_id == 0){
		j_first++;//first processor doesn't update bottom boundary
	}
	if (my_id == num_procs - 1){
		j_last--;//bottom processor doesn't update top boundary
	}
	//intialize  internal points
	for (int i = 0; i < M; i++){//loop through x's
		for (int j = 0; j < int(M / num_procs) + 2; j++){//loop through y's
			phi[j][i] = 1;
		}
	}
	//set boundary conditions
	for (int k = 0; k<M; k++){//set y boundaries, only for two processors
		if (my_id == num_procs - 1){
			phi[int(M / num_procs) + 1][k] = sin(pi*k*h);//boundary for top should only be set by one processor
		}
		if (my_id == 0){
			phi[0][k] = sin(pi*k*h)*exp(-pi);//boundary for bottom should only be set by one processor
		}
	}
	for (int k = 0; k<int(M / num_procs) + 2; k++){//set x boundaries, this is for all processors, should not 
		//matter that I am also setting the overlap points because they will be written over in the first 
		//iteration of the while loop during communication and before any jacobi steps occur
		phi[k][0] = 0;//left side is zero
		phi[k][M - 1] = 0;//right side is zero
	}
	while (norm > 0.0001){//main loop
		//do all communication before the actual jacobi step
		if (my_id<num_procs - 1){
			MPI_Send(phi[int(M / num_procs)], M, MPI_DOUBLE, my_id + 1, 0, MPI_COMM_WORLD);//send top row up if this process is not on top
			//printf("Process %d is pitching up. \n",my_id) ;//debug line
		}
		if (my_id>0){
			MPI_Recv(phi[0], M, MPI_DOUBLE, my_id - 1, 0, MPI_COMM_WORLD, &status);//recieve if this process is not on the bottom
			//printf("Process %d is catching. \n",my_id) ;//debug line
		}
		if (my_id>0){
			MPI_Send(phi[1], M, MPI_DOUBLE, my_id - 1, 1, MPI_COMM_WORLD);//send down if not on the bottom
			//printf("Process %d is pitching down. \n",my_id) ;//debug line
		}
		if (my_id<num_procs - 1){
			MPI_Recv(phi[int(M / num_procs) + 1], M, MPI_DOUBLE, my_id + 1, 1, MPI_COMM_WORLD, &status);
			//printf("Process %d is catching. \n",my_id) ;//debug line
		}
		//jacobi iteration
		sqrdiff = 0;//clear the total for each jacobi step
		for (int i = 1; i < M - 1; i++){//loop through x's, we only go to M-1 because M-1 is our boundary(fixed) 
			for (int j = j_first; j <= j_last; j++){//loop through y's
				phiNew[j][i] = 0.25 * (phi[j][i - 1] + phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]);
				sqrdiff += (phiNew[j][i] - phi[j][i])*(phiNew[j][i] - phi[j][i]);
			}//end of y loop
		}//end of x loop
		//loop to reset the old phi as the current phi
		for (int i = 1; i < M - 1; i++){//loop through x's, we only go to M-1 because M-1 is our boundary(fixed) 
			for (int j = j_first; j <= j_last; j++){//loop through y's
				phi[j][i] = phiNew[j][i];
			}//end of x loop
		}//end of y loop
		MPI_Allreduce(&sqrdiff, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		norm = sqrt(norm);
		steps++;
		if (my_id == 0){
			printf("At iteration %d, diff is %e\n", steps, norm);//debug line
		}
	}//end of while loop
	double tend = MPI_Wtime();
	runningtotal = runningtotal + (tend - tstart);//this total will have to be reduced
	MPI_Reduce(&runningtotal, &totaltime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);//arguements are (partial sums,total sum, number of elements, data type, 
	//operation to be used in reduciton, rank, and method of communication)
	if (0 == my_id) {
		printf("Calculation took %i jacobi steps to converge.\n", steps);
		printf("This took %f seconds \n ", totaltime);
	}
	MPI_Finalize();
	return 0;
}
