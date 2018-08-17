#include <stdio.h>
#include <stdlib.h>

void addVectors(int N, float *a, float *b, float *c) {

	int n;
	
	for(n=0;n<N;++n) {
	 c[n] = a[n] + b[n];
	}

}


int main(int argc, char **argv) {

	int N = 100;
	float *a = (float*) malloc(N*sizeof(float));
	float *b = (float*) malloc(N*sizeof(float));
	float *c = (float*) malloc(N*sizeof(float));

	int n;

	for(n=0;n<N;n++) {
 	  a[n] = 1+n;
	  b[n] = 1-n;
	}

	addVectors(N,a,b,c);
	
	for(n=0;n<5;++n) {
	 printf("c[%d] = %g\n",n,c[n]);
	}

	free(a);
	free(b);
	free(b);

	return 0;

}
