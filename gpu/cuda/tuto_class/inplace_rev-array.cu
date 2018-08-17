#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>

__global__ void revArray(int N, float *a) {

	int n = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(n<N/2) {
	  float a1 = a[n];
	  a[n] = a[N-1-n];
	  a[N-1-n] = a1;
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
	float *d_a;

	cudaMalloc(&d_a, N*sizeof(float));

	// Copy data from host to device
	cudaMemcpy(d_a, h_a, N*sizeof(float), cudaMemcpyHostToDevice);
	
	//save this for later
	int NthreadsPerBlock = 10;
	int NthreadBlocks = ((N/2)+NthreadsPerBlock-1)/NthreadsPerBlock ;
	revArray<<<NthreadBlocks, NthreadsPerBlock>>>(N,d_a);

	//copy result from device to host
	cudaMemcpy(h_a, d_a, N*sizeof(float), cudaMemcpyDeviceToHost);


	for(n=0;n<N;++n) {
	 printf("h_a[%d] = %g\n",n,h_a[n]);
	}

	free(h_a);

	cudaFree(d_a);

	return 0;

}
