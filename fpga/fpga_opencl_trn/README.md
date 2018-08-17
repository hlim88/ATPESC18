This is example on FPGA of ATPESC 2018
(running on JLSE)

Part A: set up the env and copy the lab files
1. Login to a shell on JLSE:
$ ssh login.jlse.anl.gov
2. Please submit for an interacive queue so that we don’t tie up resources on the jlse login node. E.g.:
$ qsub -q gomez (or it or skylake_8180 -- some queue you have access to) -I -t 120 -n 1
3. If you do not wish to work out of your home directory, cd into another directory. E.g.
$ mkdir atpesc_lab
$ cd atpesc_lab
4. Obtain the lab files from my home directory and cd into that directory:
$ cp -r /home/jmoawad/fpga_opencl_trn/ .
$ cd fpga_opencl_trn
5. Source the environment setup script:
$ source aoc-env.sh
Note that this is a slightly modified version of /soft/fpga/altera/pro/18.0.0.219/aocenv.sh
which sets the environment to use the OpenCL emulator rather than an FPGA board

Part B: compile a single work-item kernel
1. Change directories into the OpenCL_single directory. Here we have a simple OpenCL example
showing a single work-item kernel:
$ cd OpenCL_single
2. Change directories into the device directory where our OpenCL kernel code resides:
$ cd device
3. Compile the kernel SimpleKernel.cl by executing the following command:
$ aoc -v -march=emulator SimpleKernel.cl
4. Double check the kernel compiled successfully. You should now have a *.aocx file which is the
kernel executable FPGA image:
$ ls
SimpleKernel/ SimpleKernel.aoco SimpleKernel.aocr SimpleKernel.aocx
SimpleKernel.cl

Part C: compile and execute the host application
1. Change directory into the folder containing our main() c++ source project for the host:
$ cd ../host/
2. List the directory files and you should see 3 source files and a Makefile:
$ ls
main.cpp Makefile utility.cpp utility.h
3. Compile the software using a make command:
$ make
4. Run the executable that was just created:
$ ./SimpleOpenCLApp

Part D: compile a NDRange kernel
1. Change directories into the OpenCL_ndrange directory. Here we have a simple OpenCL example
showing a single work-item kernel:
$ cd ../../OpenCL_ndrange
2. Change directories into the device directory where our OpenCL kernel code resides:
$ cd device
3. Compile the kernel SimpleKernel.cl by executing the following command:
$ aoc -v -march=emulator SimpleKernel.cl
4. Double check the kernel compiled successfully. You should now have a *.aocx file which is the
kernel executable FPGA image:
$ ls
SimpleKernel/ SimpleKernel.aoco SimpleKernel.aocr SimpleKernel.aocx
SimpleKernel.cl

Part E: compile and execute the host application
1. Change directory into the folder containing our main() c++ source project for the host:
$ cd ../host/
2. List the directory files and you should see 3 source files and a Makefile:
$ ls
main.cpp Makefile utility.cpp utility.h
3. Compile the software using a make command:
$ make
4. Run the executable that was just created:
$ ./SimpleOpenCLApp

Today we used targeted our OpenCL kernel to the emulator. To target an actual
FPGA card, there is a Nallatech 385A with Intel Arria 10 on “Ruth”. The queue
name is fpga_385a
$ qsub -q fpga_385a -I -t 120 -n 1
Launch aoc without the emulator flag. Adding the -report flag will generate the
html report.
$ aoc -v –report SimpleKernel.cl

Design examples available:
• 1D FFT https://www.intel.com/content/www/us/en/programmable/support/supportresources/design-examples/design-software/opencl/fft-1d.html
• 2D FFT https://www.intel.com/content/www/us/en/programmable/support/supportresources/design-examples/design-software/opencl/fft-2d.html
• Finite Difference Computation
https://www.intel.com/content/www/us/en/programmable/support/support-resources/designexamples/design-software/opencl/fdtd-3d.html
• Many others at https://www.intel.com/content/www/us/en/programmable/products/designsoftware/embedded-software-developers/opencl/developer-zone.html#design-examples
