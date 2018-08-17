#!/bin/bash
# Intel Quartus and FPGA OpenCL env
#
TOPDIR=/soft/fpga/altera


export QUARTUS_ROOTDIR="$TOPDIR/pro/18.0.0.219"
#export AOCL_BOARD_PACKAGE_ROOT=${QUARTUS_ROOTDIR}/hld/board/nalla_pcie
export AOCL_BOARD_PACKAGE_ROOT=${QUARTUS_ROOTDIR}/hld/board/a10_ref
export CL_CONTEXT_EMULATOR_DEVICE_INTELFPGA=a10gx

#
# COMMON
#

export INTELFPGAOCLSDKROOT=${QUARTUS_ROOTDIR}"/hld"

export PATH="$QUARTUS_ROOTDIR"/quartus/linux64/jre64/bin:$PATH
export PATH=$PATH:"$QUARTUS_ROOTDIR"/bin:"$INTELFPGAOCLSDKROOT"/linux64/bin
export PATH=$PATH:"$INTELFPGAOCLSDKROOT"/bin
export PATH=$PATH:"$QUARTUS_ROOTDIR"/quartus/bin
# for HLS
export QUARTUS_ROOTDIR_OVERRIDE=$QUARTUS_ROOTDIR/quartus
export PATH=$PATH:"$QUARTUS_ROOTDIR"/quartus/sopc_builder/bin/ # for qsys-script
export PATH=$PATH:"$QUARTUS_ROOTDIR"/modelsim_ase/linuxaloem/ # for vsim

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"$INTELFPGAOCLSDKROOT"/linux64/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"$INTELFPGAOCLSDKROOT"/host/linux64/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${AOCL_BOARD_PACKAGE_ROOT}/linux64/lib

export QUARTUS_64BIT=1

export QSYS_ROOTDIR="$QUARTUS_ROOTDIR"/quartus/sopc_builder/bin

#
#
export LM_LICENSE_FILE=28518@ftsn3.ftm.alcf.anl.gov

#

echo
echo LM_LICENSE_FILE: $LM_LICENSE_FILE
echo QUARTUS_ROOTDIR: $QUARTUS_ROOTDIR
echo AOCL_BOARD_PACKAGE_ROOT: $AOCL_BOARD_PACKAGE_ROOT
echo INTELFPGAOCLSDKROOT: $INTELFPGAOCLSDKROOT                        
echo
