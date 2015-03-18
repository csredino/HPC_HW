#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include "omp.h" 
int num_threads = 2;//don't need if I set threads in bash 
//int my_id; //didn't end up using, let threads figure out the four loop 
double pi;
//void omp_set_num_threads(num_threads ); //could do this in bash instead using
int main() { 
	double totaltime=0;//intialize total time 
	int N; //number of terms in the summation 
	int max; 
	printf("Enter number of iterations for summation, N=");
	scanf("%d",&N); 
	if(N<10000){//only perform average if N is small enough that time granularity would be an issue 
	max=100000;
	} 
	else{ 
		max=1; 
	}
	for(int i=0;i<max;i++){//loop to average the time 
		double tstart = omp_get_wtime() ;//begin timer
		double sumpi = 0;//intialize sum
		#pragma omp parallel private(n) shared(N) reduction(+:sumpi)//make n private, make N public, make sumpi a reduction
		{//start of parallel region
			//my_id = omp_get_thread_num ();//didn't end up using this 
			#pragma omp for 
			for(int n =0;n<N;n++){ 
				sumpi = sumpi + (pow(-1,n))/(2*n+1);//rounding errors will be different depending on which threads go faster
		};
		}//end of parallel region
		pi=4*sumpi; 
		double tend = omp_get_wtime ( ) ;
		totaltime=totaltime+(tend-tstart);//no need to perform a reduction on total time since time calipers are outside the parallel region 
	}//end of time averaging loop 
printf ( "This took %f seconds \n " , (totaltime/max)); 
printf("The value of pi for N = %i is %.5f.\n",N,pi);
return 0;
}

