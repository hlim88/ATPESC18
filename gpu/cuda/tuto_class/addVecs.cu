#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>

__global__ void addVectors(int N, float *a, float *b, float *c) {

	int n = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(n<N) {
	 c[n] = a[n] + b[n];
	}

}


int main(int argc, char **argv) {

	int N = 100;

	//Host memory allocation
	float *h_a = (float*) malloc(N*sizeof(float));
	float *h_b = (float*) malloc(N*sizeof(float));
	float *h_c = (float*) malloc(N*sizeof(float));

	int n;

	for(n=0;n<N;n++) {
 	  h_a[n] = 1+n;
	  h_b[n] = 1-n;
	}

	// Device memory allocation
	float *d_a, *d_b, *d_c;

	cudaMalloc(&d_a, N*sizeof(float));
	cudaMalloc(&d_b, N*sizeof(float));
	cudaMalloc(&d_c, N*sizeof(float));

	// Copy data from host to device
	cudaMemcpy(d_a, h_a, N*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, h_b, N*sizeof(float), cudaMemcpyHostToDevice);
	
	//save this for later
	int NthreadsPerBlock = 10;
	int NthreadBlocks = (N+NthreadsPerBlock-1)/NthreadsPerBlock ;
	addVectors<<<NthreadBlocks, NthreadsPerBlock>>>(N,d_a,d_b,d_c);

	//copy result from device to host
	cudaMemcpy(h_c, d_c, N*sizeof(float), cudaMemcpyDeviceToHost);


	for(n=0;n<5;++n) {
	 printf("h_c[%d] = %g\n",n,h_c[n]);
	}

	free(h_a);
	free(h_b);
	free(h_b);

	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	return 0;

}
