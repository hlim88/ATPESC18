#!/bin/bash -l
#
#PBS -l nodes=1
#PBS -l walltime=00:05:00
#PBS -d .

regent 2.rg -logfile $PBS_O_WORKDIR/prof0 -lg:prof 1 -ll:cpu 4 -lg:serializer ascii


