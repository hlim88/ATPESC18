#include<stdio.h>
#include<omp.h>

static long num_steps = 100000000;
double step;
#define MIN_BLK 100000000
#define NUM_THREADS 2

void main() {

	int i, nthreads;
	double pi, sum[NUM_THREADS]; 
	
	step = 1.0/(double) num_steps;
	omp_set_num_threads(NUM_THREADS);	
	double tdata = omp_get_wtime();

	#pragma omp parallel
	{
		int i, id, nthrds;
		double x;
		id = omp_get_thread_num();
		nthrds = omp_get_num_threads();
		if(id == 0) nthreads = nthrds;	
	
		for(i = id, sum[id]=0.0; i<num_steps; i=i+nthrds) 
		{
		  x = (i+0.5)*step;
		  sum[id] += 4.0/(1.0+x*x);
		}
	}
	
	for(i=0,pi=0.0;i<nthreads;i++) pi += sum[i]*step;
	tdata = omp_get_wtime() - tdata;
	printf("pi = %f in %f secs\n", pi, tdata);

}
