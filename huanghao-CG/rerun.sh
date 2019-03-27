#!/bin/bash
#PBS -N rerun
#PBS -q fast
#PBS -l nodes=1:ppn=3
#PBS -l walltime=1:00:00

cd $PBS_O_WORKDIR

gmx grompp -f md.mdp -c equ.gro -p rerun.top -o rerun.tpr
gmx_fast mdrun -ntomp 3 -deffnm rerun -rerun md.xtc

