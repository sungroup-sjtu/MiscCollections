#!/usr/bin/env python 
#coding=utf-8

import os, sys

e_list = [0.529]
d_list = [4.185]

def write_in(x1,y1,y2):
    x1.write("""
units          real
atom_style     full
boundary       p p p
pair_style     lj/cut     20.0000
pair_modify    tail no
special_bonds  lj 0.0 0.0 1.0
pair_modify    mix arithmetic
bond_style     class2
angle_style    class2
dihedral_style opls
improper_style none

read_data      PEO-st.data

pair_coeff     1   1   0.444   4.329 
pair_coeff     2   2    %.3f   %.3f 

thermo_style   multi
thermo         100
min_style      cg
minimize       1.0e-4 1.0e-6 5000 1000

variable       T equal temp
variable       P equal press
variable       E equal epair
variable       st equal (lz/2.0)*(pzz-(pxx+pyy)/2)*0.01

fix            1 all nvt temp 423 423 500.0 
fix            2 all ave/time 4000 50 200000 v_T v_E v_st file equ
velocity       all create 423 1166140691 mom yes rot yes dist gaussian
thermo         1000
thermo_style   custom step temp press pe v_E v_P v_T 
reset_timestep 0
timestep       5
run            200000
unfix          1
unfix          2
write_restart  equ.rst

reset_timestep 0
timestep       5
fix            1 all nvt temp 423 423 500.0 
fix            2 all ave/time 5000 40 200000 v_T v_E v_st file 423K
thermo         1000
thermo_style   custom step temp press pe v_E v_P v_T
dump           1 all custom 1000 423K.atom id type xu yu zu
dump_modify    1 sort id
run            2000000
write_data     423K.data
""" %(y1,y2))

def write_sh(x2,y2,z2):
    path = os.getcwd()
    x2.write("""
#!/bin/bash
#PBS -o PEO-st.out
#PBS -e PEO-st.err
#PBS -N %s
#PBS -q fast
#PBS -l nodes=1:ppn=%i

cd  %s
mpirun -np %i lmp-stable -i 423K.in > 423K.log
""" %(y2,z2,path,z2))

for y1 in e_list:
	for y2 in d_list:
		dirname = "%3d%3d" %(y1*1000,y2*1000)
		if not os.path.exists(dirname):
			os.mkdir(dirname)
			os.chdir(dirname)
			os.system("cp ../PEO-st.data ./")
			write_in(open("423K.in","w"),y1,y2)
			write_sh(open("423K.sh","w"),dirname,3)
			print os.getcwd()
			os.system("qsub 423K.sh")
			os.chdir("..")
