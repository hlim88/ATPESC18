//OpenCL AOC Kernel
__kernel void SimpleKernel (__global float * restrict a, 
				  __global float * restrict b, 
				  __global float * restrict x) {

       size_t i = get_global_id(0);
       x[i] = a[i] * b[i];
}

