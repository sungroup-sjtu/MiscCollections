#coding=utf-8
import os, sys

T_list=range(100,420,20)

def write_in(x1,y1):
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

read_data      PEO-15.data

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
timestep       1
fix            1 all npt temp %d 800 100.0 iso 1 1 1000
fix            2 all ave/time 5000 40 200000 v_T v_P v_E v_den file increaseT
velocity       all create %d 1166140691 mom yes rot yes dist gaussian
thermo         1000
thermo_style   custom step temp press pe v_E v_P v_T 
run            400000
unfix          1
unfix          2
write_data     increaseT.data

reset_timestep 0
timestep       5
fix            1 all nvt temp 800 800 500.0
fix            2 all ave/time 5000 40 200000 v_T v_P v_E v_den file 800K
thermo         1000
thermo_style   custom step temp press pe v_E v_P v_T 
run            1000000
unfix          1
unfix          2
write_data     800K.data

reset_timestep 0
timestep       5
fix            1 all npt temp 800 %d 500.0 iso 1 1 5000
fix            2 all ave/time 5000 40 200000 v_T v_P v_E v_den file decreaseT
thermo         1000
thermo_style   custom step temp press pe v_E v_P v_T 
run            2000000
unfix          1
unfix          2
write_data     800K.data

reset_timestep 0
timestep       5
fix            1 all npt temp %d %d 500.0 iso 1 1 5000
fix            2 all ave/time 5000 40 200000 v_T v_P v_E v_den file equ
thermo         1000
thermo_style   custom step temp press pe v_E v_P v_T 
run            2000000
unfix          1
unfix          2
write_data     equ.data

reset_timestep 0
timestep       5
fix            1 all npt temp %d %d 500.0 iso 1 1 5000
fix            2 all ave/time 5000 40 200000 v_T v_P v_E v_den file den
thermo         1000
thermo_style   custom step temp press pe v_E v_P v_T
dump           1 all custom 1000 den.atom id type xu yu zu
dump_modify    1 sort id
run            2000000
write_data     den.data
""" %((0.48328-9.1741*y1*10**-5),(4.29887+7.0922*y1*10**-5),(0.5713-1*y1*10**-4),(4.147+9.0909*y1*10**-5),y1,y1,y1,y1,y1,y1,y1))

def write_sh(x2,y2,z2):
    path = os.getcwd()
    x2.write("""#!/bin/bash
#SBATCH -J %s
#SBATCH --partition=fast
#SBATCH --time=333:0:0
#SBATCH --ntasks=%i

mpirun -np %i lmp-stable -i den.in > den.log
""" %(y2,z2,z2))

for y1 in T_list:
	dirname = "%d" %y1
	if not os.path.exists(dirname):
		os.mkdir(dirname)
	os.chdir(dirname)
	os.system("cp ../PEO-15.data ./")
	write_in(open("den.in","w"),y1)
	write_sh(open("den.sh","w"),dirname,6)
	print os.getcwd()
	os.system("sbatch den.sh")
	os.chdir("..")
