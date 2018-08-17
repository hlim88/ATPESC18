#include <math.h>
#include <fstream>
#include <assert.h>
#include <stdio.h>
#include <string>

#include "CL/cl.hpp"
#include "utility.h"

static const cl_uint vectorSize = 4096; //must be evenly disible by workSize
static const cl_uint workSize = 256;

int main(void)
{
  cl_int err;

  //Setup Platform

  //Get Platform ID
  std::vector<cl::Platform> PlatformList;
  //Create the OpenCL program object from the aocx binary
  err = cl::Platform::get(&PlatformList);
  assert(err==CL_SUCCESS);
  checkErr(PlatformList.size()==1 ? CL_SUCCESS : -1, "cl::Platform::get");
  print_platform_info(&PlatformList);
  
  //Setup Device
  //Get Device ID
  std::vector<cl::Device> DeviceList;
  //Create the OpenCL kernel object from the program
  err = PlatformList[0].getDevices(CL_DEVICE_TYPE_ALL, &DeviceList);
  assert(err==CL_SUCCESS);
  print_device_info(&DeviceList);
  
  //Create Context
  //Set up the kernel argument for each buffer in the DeviceList 
  cl::Context mycontext(DeviceList);
  assert(err==CL_SUCCESS);
  
  //Create Command queue
  //Launch the kernel by enqueueing the tasks
  cl::CommandQueue myqueue(mycontext, DeviceList[0]);
  assert(err==CL_SUCCESS);

  //Create Buffers for input and output
  cl::Buffer Buffer_In(mycontext, CL_MEM_READ_ONLY, sizeof(cl_float)*vectorSize);
  cl::Buffer Buffer_In2(mycontext, CL_MEM_READ_ONLY, sizeof(cl_float)*vectorSize);
  cl::Buffer Buffer_Out(mycontext, CL_MEM_WRITE_ONLY, sizeof(cl_float)*vectorSize);

  //Inputs and Outputs to Kernel, X and Y are inputs, Z is output
  //The aligned attribute is used to ensure alignment
  //so that DMA could be used if we were working with a real FPGA board
  cl_float X[vectorSize]  __attribute__ ((aligned (64)));
  cl_float Y[vectorSize]  __attribute__ ((aligned (64)));
  cl_float Z[vectorSize]  __attribute__ ((aligned (64)));

  //Allocates memory with value from 0 to 1000
  cl_float LO= 0;   cl_float HI=1000;
  fill_generate(X, Y, Z, LO, HI, vectorSize);

  //Write data to device
  err = myqueue.enqueueWriteBuffer(Buffer_In, CL_FALSE, 0, sizeof(cl_float)*vectorSize, X);
  err = myqueue.enqueueWriteBuffer(Buffer_In2, CL_FALSE, 0, sizeof(cl_float)*vectorSize, Y);
  assert(err==CL_SUCCESS);
  myqueue.finish();

  // create the kernel
  const char *kernel_name = "SimpleKernel";

  //Read in binaries from file
  std::ifstream aocx_stream("../device/SimpleKernel.aocx", std::ios::in|std::ios::binary);
  checkErr(aocx_stream.is_open() ? CL_SUCCESS:-1, "SimpleKernel.aocx");
  std::string prog(std::istreambuf_iterator<char>(aocx_stream), (std::istreambuf_iterator<char>()));
  cl::Program::Binaries mybinaries (1, std::make_pair(prog.c_str(), prog.length()+1));

  // Create the Program from the AOCX file.
  cl::Program program(mycontext, DeviceList, mybinaries);

  // build the program
  //////////////      Compile the Kernel.... For Intel FPGA, nothing is done here, but this comforms to the standard
  err=program.build(DeviceList);
  assert(err==CL_SUCCESS);

  // create the kernel
  //////////////       Find Kernel in Program
  cl::Kernel kernel(program, kernel_name, &err);
  assert(err==CL_SUCCESS);
  
  //////////////     Set Arguments to the Kernels
  err = kernel.setArg(0, Buffer_In);
  assert(err==CL_SUCCESS);
  err = kernel.setArg(1, Buffer_In2);
  assert(err==CL_SUCCESS);
  err = kernel.setArg(2, Buffer_Out);
  assert(err==CL_SUCCESS);
  //err = kernel.setArg(3, vectorSize);
  //assert(err==CL_SUCCESS);

  printf("\nLaunching the kernel...\n");
  
  // Launch Kernel
  err=myqueue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(vectorSize), cl::NDRange(workSize));// It needs to uncomment

  // Below two commented one for single case
  //err=myqueue.enqueueTask(kernel); 
  //assert(err==CL_SUCCESS);

  // read the output
  err=myqueue.enqueueReadBuffer(Buffer_Out, CL_TRUE, 0, sizeof(cl_float)*vectorSize, Z);
  assert(err==CL_SUCCESS);

  err=myqueue.finish();
  assert(err==CL_SUCCESS);
  
  float CalcZ[vectorSize];

  for (int i=0; i<vectorSize; i++)
    {
      //////////////  Equivalent Code running on CPUs
      CalcZ[i] = X[i] * Y[i]; 
      
    }

  //Print Performance Results
  verification (X, Y, Z, CalcZ, vectorSize);

  return 1;
}

