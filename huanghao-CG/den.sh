#!/bin/bash
#SBATCH -J 5144147
#SBATCH --partition=cpu
#SBATCH --time=333:0:0
#SBATCH --ntasks=8

mpirun -np 8 lmp-stable -i den.in > den.log
