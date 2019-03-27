#coding=utf-8
import os, sys

e_list = [0.509,0.504,0.494,0.489]
d_list = [4.180]

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

read_data      PEO.data

pair_coeff     1   1    %.3f   %.3f  
pair_coeff     2   2    %.3f   %.3f 

thermo_style   multi
thermo         100
min_style      cg
minimize       1.0e-4 1.0e-6 5000 1000

compute        Ei all inter all
variable       T equal temp
variable       P equal press
variable       E equal epair
variable       den equal density

reset_timestep 0
timestep       5
fix            1 all npt temp 363 363 500.0 iso 1 1 5000
fix            2 all ave/time 5000 40 200000 v_T v_P v_E v_den file equ
velocity       all create 363 1166140691 mom yes rot yes dist gaussian
thermo         1000
thermo_style   custom step temp press pe v_E v_P v_T 
run            2000000
unfix          1
unfix          2
write_data     equ.data

reset_timestep 0
timestep       5
fix            1 all npt temp 363 363 500.0 iso 1 1 5000
fix            2 all ave/time 5000 40 200000 v_T v_P v_E v_den file den
thermo         1000
thermo_style   custom step temp press pe v_E v_P v_T
dump           1 all custom 1000 den.atom id type xu yu zu
dump_modify    1 sort id
run            5000000
write_data     den.data
""" %((0.48328-9.1741*363*10**-5),(4.29887+7.0922*363*10**-5),y1,y2))

def write_sh(x2,y2,z2):
    path = os.getcwd()
    x2.write("""#!/bin/bash
#SBATCH -J %s
#SBATCH --partition=fast
#SBATCH --time=333:0:0
#SBATCH --ntasks=%i

mpirun -np %i lmp-stable -i den.in > den.log
""" %(y2,z2,z2))

for y1 in e_list:
	for y2 in d_list:
		dirname = "%3d%3d" %(y1*1000,y2*1000)
		if not os.path.exists(dirname):
			os.mkdir(dirname)
			os.chdir(dirname)
			os.system("cp ../PEO.data ./")
			write_in(open("den.in","w"),y1,y2)
			write_sh(open("den.sh","w"),dirname,6)
			print os.getcwd()
			os.system("sbatch den.sh")
			os.chdir("..")
		
