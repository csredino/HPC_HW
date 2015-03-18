#include "mpi.h"
#include <time.h> 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char *argv[]) {
	double totaltime;
	double runningtotal = 0;//intialize total time
	int N; //number of terms in the summation
	int max;
	int iterationmax;
	printf("Enter number of iterations for summation, N=");
	scanf("%d", &N);
	int my_id;
	double pi;
	int num_procs;
	int elem_name;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Get_processor_name(processor_name, &elem_name);
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);//arguements are (what i broadcast, how many elements it is, what its type is, what its rank will be, and how I communicate it)
	int x = int(N / num_procs);
	int low = my_id * x;
	int high = low + x;
	if (N<10000){//only perform average if N is small enough that time granularity would be an issue
		iterationmax = 100000;
	}
	else{
		iterationmax = 1;
	}
	for (int i = 0; i<iterationmax; i++){//loop to average the time
		double sumpi = 0;//intialize sum
		double tstart = MPI_Wtime();
		for (int n = low; n<high; n++){
			sumpi = sumpi + (pow(-1, n)) / (2 * n + 1);//rounding errors will be different depending on which threads go faster
		};
		MPI_Reduce(&sumpi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);//arguements are (partial sums,total sum, number of elements, data type, operation to be used in reduciton, rank, and method of communication)
		pi = 4 * pi;
		double tend = MPI_Wtime();
		runningtotal = runningtotal + (tend - tstart);//this total will have to be reduced
	}//end of time averaging loop
	MPI_Reduce(&runningtotal, &totaltime, 1, MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);//arguements are (partial sums,total sum, number of elements, data
	type, operation to be used in reduciton, rank, and method of communication)
		if (0 == my_id) {
		printf("This took %f seconds \n ", (totaltime / iterationmax));
		printf("The value of pi for N = %i is %.5f.\n", N, pi);
		}
	MPI_Finalize(); 
	return 0;
}