#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
extern "C" double ddot_(int *N, double *dx, int *incx, double *dy, int *incy);//prototypes blas ddot function
int main() {
	for (int n = 1; n<6; n++){//master loop to change the size of vectors to test
		int N = n * 10000;//number of elements in each vector
		int inc1 = 1;//stride of the vectors to be multiplied,
		int inc2 = inc1;
		clock_t totalt1 = 0; //intialize total time to zero
		for (int j = 0; j<10000; j++){ //loop to perform dot product many times (10000)
			double vector1[N];
			double vector2[N];
			for (int i = 0; i<N; i++){//loop to build vectors
				vector1[i] = 1;//filling vectors with 1s just for simplicity
				vector2[i] = 1;
			};
			clock_t t1 = clock();//start timer
			double product = ddot_(&N, vector1, &inc1, vector2, &inc2);
			clock_t t2 = clock();//stop timer
			clock_t t = t2 - t1;//calc time elapsed
			totalt1 = totalt1 + t;//update running total time
		}
		double averaget1 = ((double)totalt1 / (double)CLOCKS_PER_SEC) / 10000;//calculates average and casts as a double, also, in this step the clock time is converted to seconds
		printf("Vectors with %i elements took an average time of %.5f seconds to multiply.\n", N, averaget1);
	}
}