#!/usr/bin/env python 
#coding=utf-8

import os, sys

e_list = [0.488,0.489,0.490,0.491,0.492,0.493,0.494,0.495,0.496]

def write_in(x1,y1):
    x1.write("""
units          real
atom_style     full
boundary       p p p

pair_style     lj/cut     20.0000
#pair_modify    mix arithmetic
pair_modify    tail no
#kspace_style   pppm 1.0e-4
#dielectric     1.0
special_bonds  lj 0.0 0.0 1.0
bond_style     harmonic
angle_style    harmonic
dihedral_style  opls
improper_style none

read_data      den.data

pair_coeff     1   1   0.5230    5.673
pair_coeff     1   2   %.4f      5.621
pair_coeff     2   2   %.4f      5.569

reset_timestep 0
compute        Ei all inter all
variable       T equal temp
variable       P equal press
variable       E equal epair
variable       hov equal -c_Ei[1]/250+8.314*300/4184
variable       den equal density

timestep       10
fix            1 all npt temp 300 300 1000.0 iso 1.0 1.0 10000.0
velocity       all create 300 1166140691 mom yes rot yes dist gaussian
fix            2 all ave/time 5000 20 100000 v_T v_P v_E v_den v_hov file hov
thermo         1000
thermo_style   custom step temp press pe v_T v_E v_P
dump           1 all custom 1000 hov.atom id type xu yu zu
dump_modify    1 sort id
run            1000000
write_data     hov.data
""" %((0.523*y1)**0.5,y1))

def write_sh(x2,y2,z2):
    path = os.getcwd()
    x2.write("""
#!/bin/bash
#PBS -o PDMS.out
#PBS -e PDMS.err
#PBS -N %s
#PBS -q batch
#PBS -l nodes=1:ppn=%i

cd  %s
mpirun -np %i lmp-mod -i hov.in > hov.log
""" %(y2,z2,path,z2))

for y1 in e_list:
    dirname = "%.3f" %y1
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    os.system("cp ../den.data ./")
    write_in(open("hov.in","w"),y1)
    write_sh(open("hov.sh","w"),dirname,8)
    print os.getcwd()
    os.system("qsub hov.sh")
    os.chdir("..")
