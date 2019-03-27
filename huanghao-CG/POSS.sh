#!/bin/bash
#PBS -N POSS
#PBS -q gtx
#PBS -l nodes=1:ppn=2
#PBS -l walltime=48:00:00

cd $PBS_O_WORKDIR

gmx grompp -f mini.mdp -c POSS-40.gro -p POSS-40.top -o mini.tpr
mpirun -np 4 gmx_gpu mdrun -deffnm mini -ntomp 8 -gpu_id 0101
gmx grompp -f anneal.mdp -c mini.gro -p POSS-40.top -o anneal.tpr
mpirun -np 4 gmx_gpu mdrun -deffnm anneal -ntomp 8 -gpu_id 0101
gmx grompp -f equ.mdp -c anneal.gro -p POSS-40.top -o equ.tpr
mpirun -np 4 gmx_gpu mdrun -deffnm equ -ntomp 8 -gpu_id 0101
gmx grompp -f md.mdp -c equ.gro -p POSS-40.top -o md.tpr
mpirun -np 4 gmx_gpu mdrun -deffnm md -ntomp 8 -gpu_id 0101



