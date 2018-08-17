//OpenCL AOC Kernel
__kernel 
void SimpleKernel (__global float * restrict a, 
				  __global float * restrict b, 
				  __global float * restrict x, 
				  uint array_size) {
int i;
for (i = 0; i < array_size; i++)
       x[i] = a[i] * b[i];
}

