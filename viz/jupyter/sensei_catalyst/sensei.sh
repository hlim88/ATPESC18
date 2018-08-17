#!/bin/bash

export DISPLAY=:0.0
mpiexec -n 8  /projects/ATPESC2018/jupyter/sensei-catalyst/bin/oscillator  -f /home/hyunlim/atpesc/viz/jupyter-ATPESC2018/sensei_catalyst/oscillator.xml  -b 16  -t 0.2  -s "128,128,128"  /home/hyunlim/atpesc/viz/jupyter-ATPESC2018/sensei_catalyst/sample.osc
 
