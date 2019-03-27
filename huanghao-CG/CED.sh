#!/bin/bash
#PBS -o CED.out
#PBS -e CED.err
#PBS -N CED
#PBS -q fast
#PBS -l nodes=1:ppn=3
#PBS -l walltime=1:00:00

cd $PBS_O_WORKDIR
mpirun -np 3 lmp-stable -i CED.in > CED.log
