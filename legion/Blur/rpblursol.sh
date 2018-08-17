#!/bin/bash -l
#
#PBS -l nodes=1
#PBS -l walltime=00:05:00
#PBS -d .

regent blur_solution.rg  -logfile $PBS_O_WORKDIR/prof0 -lg:prof 1 -ll:cpu 4 -i images/earth.png -p 8 -lg:serializer ascii


