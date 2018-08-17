#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>

__global__ void revArray(int N, float *a, float *b) {

	int n = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(n<N) {
	 b[N-1-n] = a[n];
	}

}


int main(int argc, char **argv) {

	int N = 100;

	//Host memory allocation
	float *h_a = (float*) malloc(N*sizeof(float));
	float *h_b = (float*) malloc(N*sizeof(float));
	int n;

	for(n=0;n<N;n++) {
 	  h_a[n] = 1+n;
	}

	// Device memory allocation
	float *d_a, *d_b;

	cudaMalloc(&d_a, N*sizeof(float));
	cudaMalloc(&d_b, N*sizeof(float));

	// Copy data from host to device
	cudaMemcpy(d_a, h_a, N*sizeof(float), cudaMemcpyHostToDevice);
	
	//save this for later
	int NthreadsPerBlock = 10;
	int NthreadBlocks = (N+NthreadsPerBlock-1)/NthreadsPerBlock ;
	revArray<<<NthreadBlocks, NthreadsPerBlock>>>(N,d_a,d_b);

	//copy result from device to host
	cudaMemcpy(h_b, d_b, N*sizeof(float), cudaMemcpyDeviceToHost);


	for(n=0;n<N;++n) {
	 printf("h_b[%d] = %g\n",n,h_b[n]);
	}

	free(h_a);
	free(h_b);

	cudaFree(d_a);
	cudaFree(d_b);

	return 0;

}
