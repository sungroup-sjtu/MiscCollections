#!/usr/bin/env python
"""This script generate and submit the simulation job for Materials Simulation Data.
usage creat_job.py -b f f f -t f f f -p den hov ... -m water methane ... -c isobar -j s -computer cluster 
    -b: the pressure range used in simulation
    -t: the temperature range used in simulation
	These two keywords can follow three or one real number. 
	Three real number represent start, end, and interval value.
	One means only one thermodynamic point is considered.
    -p: the properties that would be calculated.
	available properties incluing:
	den: density	               hov: heat of vaporization
	cp:  heat of capacity          ap: coefficient of thermal expansion 
	kt: isothermal compressibility er: permittivity
	rdf:                           rg: radius of gyration                           
        dc: diffusion coefficient      vis: viscosity                          
        td: thermal conductivity       st: surface tension            
        vle: vapor liquid equilibrium  svp: saturated vapor pressure
	ct: critical temperature       cd: critical density
	cb: critical pressure          nbp: normal boiling point
    -m: the molecules that would be simulated
	Each molecule have to exist a data file named *.data
    -c: "isobar" or "isothermo" curves
    -j: job running type
	s: a curve is simulated using a job through step by step changing thermodynamic state.
	i: each thermodynamic state is simulated using individual job.
    -computer: which computer will be used
	cluster: computer in our lab 
	HPC: HPC of SJTU

    -b, -t, and -j terms have default values:
    -b: 1.0
    -t: 300.0
    -j: i
"""


import os
import sys
import time
JOB_HOME = os.getcwd()
NPT_properties = ["den","hov","cp","kt","ap","dc","rdf","rg","vis","td","er"]
NVT_properties = ["st","vle","ct","cd","cb","svp","nbp"]
all_properties = NPT_properties + NVT_properties


quence = "power"
nproc = "8"
sh_cluster_script = """
#!/bin/bash
 
#PBS -l nodes=1:ppn=%NPROC%
#PBS -q %QUENCE%  
echo PBS work @ $PBS_NODEFILE
JOB_HOME=%PATH%
cd /state/partition1/jobs
mkdir ${USER}
cd ${USER}
mkdir %NAME%
cd %NAME%
cp -rf $JOB_HOME/%NAME%.in .
cp -rf $JOB_HOME/%MOLE%_%MODEL%.data .
date >date_initial.txt
mpirun -np %NPROC% lmp-stable < %NAME%.in > %NAME%.log
date >date_final.txt
cp ./* $JOB_HOME/
cd ..
rm -r %NAME%
"""

sh_hpc_script = """
#!/bin/bash

#BSUB -J %NAME%.sh
#BSUB -q cpu
#BSUB -o %NAME%.out
#BSUB -e %NAME%.err
#BSUB -n 16
#BSUB -R "span[hosts=1]"

############################################

echo ""
echo "----------------------- INTIALIZATIONS -----------------------------"
echo ""

EXE=lammps_cao

if test ! -x `which $EXE` ; then
    echo
    echo "ERROR: `which $EXE` not existent or not executable"
    echo "Aborting"
    exit 1
fi

CURDIR=$PWD
cd $CURDIR
cd %PATH%

rm -f nodelist >& /dev/null
touch nodelist

for host in `echo $LSB_HOSTS`
do
echo $host >> nodelist
done

NP=`cat nodelist |wc -l`
NN=`cat nodelist |sort |uniq|tee nodes |wc -l`

echo
echo "Executable : `which $EXE`"
echo "Working directory is $CURDIR"
echo "Running on host `hostname`"
echo "Directory is `pwd`"
echo "This jobs runs on $NN nodes"
echo "This job has allocated $NP core(s)"
echo

ls -al

echo ""
echo "----------------------- RUN -----------------------------"
echo ""

date '+RUN STARTED ON %m/%d/%y AT %H:%M:%S'

#/lustre/utility/mpich2-1.4.1p1-gcc/bin/mpirun -np $NP -machinefile nodelist lammps_cao < %NAME%.in > %NAME%.log

#/lustre/utility/intel/impi/4.1.1.036/intel64/bin/mpirun -np $NP -machinefile nodelist lmp_cao_Nov29 -i %NAME%.in > %NAME%.log
mpirun -np $NP -machinefile nodelist lmp-sz -i %NAME%.in > %NAME%.log
#%EXE% %NAME%.in > %NAME%.log
date '+RUN ENDED ON %m/%d/%y AT %H:%M:%S'

echo ""
echo "----------------------- DONE ----------------------------"
echo ""

ls -al

"""

inscript_head="""
units              real
atom_style         full
boundary           p p p 

pair_style         lj/cut/coul/long 12.0
pair_modify        mix arithmetic
pair_modify        tail yes
kspace_style       pppm 0.0001
bond_style         class2
angle_style        class2
dihedral_style     opls
improper_style     class2

read_data          $NAME$_$MODEL$.data

group              $NAME$ molecule <> 1 $MOL_NUMBER$
minimize           1.0e-5 1.0e-7 1000 10000

timestep           1.0	
thermo             1000	
thermo_style       custom step temp press vol etotal pe	
"""

def script_add(en,model):
	inscript_tmp  = "dump              1 all atom 1000 $NAME$_$TEMP$_$PRESS$_%s_%s.lammpstrj\n" %(en,model)
	inscript_tmp += "dump_modify       1 sort id\n" 
	inscript_tmp += "run               1000000\n" 
	inscript_tmp += "write_restart     $NAME$_$TEMP$_$PRESS$_%s_%s.rst\n" %(en,model)
	inscript_tmp += "undump            1\n"  
	return(inscript_tmp)

def creat_inscript(properties):
	sur_flag = 0
	inscript_body = [["","","","","",""],["","","","","",""]]
	inscript_body[0][1] += "variable          T equal temp\n" 
	inscript_body[0][1] += "variable          P equal press\n" 
	inscript_body[0][1] += "variable          V equal vol\n" 
	inscript_body[0][1] += "variable          H equal enthalpy\n" 
	inscript_body[0][1] += "variable          nmole equal $MOL_NUMBER$\n" 
	inscript_body[0][1] += "variable          Tref equal $TEMP$\n" 
	inscript_body[0][1] += "variable          Pref equal $PRESS$\n" 
	inscript_body[0][1] += "variable          kb equal 1.38\n" 
	inscript_body[0][1] += "variable          Na equal 6.02\n" 
	inscript_body[0][2] += "fix               npt all npt temp $TEMP$ $TEMP$ 100.0 iso $PRESS$ $PRESS$ 1000.0 \n" 
	inscript_body[0][2] += "run               1000000\n" 
	inscript_body[1][2] += "fix               nvt all nvt temp $TEMP$ $TEMP$ 100.0\n" 
	inscript_body[1][2] += "run               1000000\n" 
	for p in properties:
		if p == "den":
			inscript_body[0][1] += "variable          den equal density\n" 
			inscript_body[0][2] += "fix               denout all ave/time 100 1000 100000 v_den file den_$NAME$_$TEMP$_$PRESS$.log\n" 
			inscript_body[0][3] += "unfix             denout\n" 
		elif p == "hov":
			inscript_body[0][0] += "#compute           inter all inter $NAME$\n"
			inscript_body[0][0] += "compute            inter $NAME$ group/group $NAME$\n"
			inscript_body[0][1] += "variable          hov equal -c_inter/$MOL_NUMBER$+8.314*v_Tref/4184\n"
			inscript_body[0][2] += "fix               hovout all ave/time 100 1000 100000 v_hov file hov_$NAME$_$TEMP$_$PRESS$.log\n"
			inscript_body[0][3] += "unfix             hovout\n" 
		elif p == "cp":
			inscript_body[0][1] += "variable          HH equal v_H*v_H\n" 
			inscript_body[0][1] += "variable          cpfactor equal 503.2166*4184/v_Tref/v_Tref\n" 
			inscript_body[0][2] += "fix               cpout all ave/time 100 1000 100000 v_cpfactor v_nmole v_HH v_H v_H file cp_$NAME$_$TEMP$_$PRESS$.log\n"
			inscript_body[0][3] += "unfix             cpout\n" 
		elif p == "kt":
			inscript_body[0][1] += "variable          VV equal v_V*v_V\n" 
			inscript_body[0][1] += "variable          ktfactor equal 0.07246/v_Tref\n" 
			inscript_body[0][2] += "fix               ktout all ave/time 100 1000 100000 v_ktfactor v_V v_VV v_V v_V file kt_$NAME$_$TEMP$_$PRESS$.log\n"
			inscript_body[0][3] += "unfix             ktout\n" 
		elif p == "ap":
			inscript_body[0][1] += "variable          VH equal v_V*v_H\n" 
			inscript_body[0][1] += "variable          apfactor equal 503.2166/v_Tref/v_Tref\n" 
			inscript_body[0][2] += "fix               apout all ave/time 100 1000 100000 v_apfactor v_V v_VH v_V v_H file ap_$NAME$_$TEMP$_$PRESS$.log\n"
			inscript_body[0][3] += "unfix             apout\n"
		elif p == "rdf":
			inscript_body[0][0] += "compute           11rdf  all rdf 100 1 1\n" 
			inscript_body[0][2] += "fix               rdfout all ave/time 1000 1000 1000000 c_11rdf file rdf_$NAME$_$TEMP$_$PRESS$.log mode vector\n" 
			inscript_body[0][3] += "unfix             rdfout\n"
		elif p == "rg":
			inscript_body[0][0] += "compute           chunkrg all chunk/atom molecule nchunk once ids once\n" 
			inscript_body[0][0] += "compute           rg all gyration/chunk chunkrg\n"
			inscript_body[0][2] += "fix               rgout all ave/time 1000 1000 1000000 c_rg file rg_$NAME$_$TEMP$_$PRESS$.log mode vector\n"
			inscript_body[0][3] += "unfix             rgout\n"
		elif p == "dc":
			inscript_body[0][0] += "compute           chunkmsd all chunk/atom molecule nchunk once ids once\n" 
			inscript_body[0][0] += "compute           msd all msd/chunk chunkmsd\n"
			inscript_body[0][0] += "#compute           vacf all vacf\n"
			inscript_body[0][4] += "fix               msdout all ave/time 10 10 100 c_msd[4] file dc_$NAME$_$TEMP$_$PRESS$.log mode vector\n"
			inscript_body[0][4] += "#fix               msdout all ave/chunk 100 1 100 chunk c_msd[4] file dc_$NAME$_$TEMP$_$PRESS$.log\n"
			inscript_body[0][4] += "#fix               vacfout all ave/time 100 1 100 c_vacf[4] file vacf_$NAME$_$TEMP$_$PRESS$.log\n"
			inscript_body[0][5] += "#unfix             vacfout\n"
			inscript_body[0][5] += "unfix             msdout\n"
		elif p == "vis":
			inscript_body[0][1] += "variable          visfactor equal 0.000000001*0.74397*100*v_V/v_Tref\n"
			inscript_body[0][1] += "variable          pxy equal pxy\n"
			inscript_body[0][1] += "variable          pyz equal pyz\n"
			inscript_body[0][1] += "variable          pxz equal pxz\n"
			inscript_body[0][4] += "fix               nonPout all ave/time 100 1 100 v_pxy v_pyz v_pxz v_visfactor file nonp_$NAME$_$TEMP$_$PRESS$.log\n"
			inscript_body[0][4] += "fix               visacfout1 all ave/correlate 100 1000 1000000 v_pxy v_pxy type upper file visacf1_$NAME$_$TEMP$_$PRESS$.log\n" 
			inscript_body[0][4] += "fix               visacfout2 all ave/correlate 100 1000 1000000 v_pyz v_pyz type upper file visacf2_$NAME$_$TEMP$_$PRESS$.log\n" 
			inscript_body[0][4] += "fix               visacfout3 all ave/correlate 100 1000 1000000 v_pxz v_pxz type upper file visacf3_$NAME$_$TEMP$_$PRESS$.log\n" 
			inscript_body[0][5] += "unfix             nonPout\n" 
			inscript_body[0][5] += "unfix             visacfout1\n"
			inscript_body[0][5] += "unfix             visacfout2\n"
			inscript_body[0][5] += "unfix             visacfout3\n"
		elif p == "td":
			inscript_body[0][1] += "variable          tdfactor equal 10000000*3.50035*100*v_V/v_Tref/v_Tref\n"
			inscript_body[0][0] += "compute           ke all ke/atom\n" 
			inscript_body[0][0] += "compute           pe all pe/atom\n" 
			inscript_body[0][0] += "compute           stress all stress/atom NULL virial\n" 
			inscript_body[0][0] += "compute           hflux all heat/flux ke pe stress\n" 
			inscript_body[0][4] += "fix               hfluxout all ave/time 100 1 100 c_hflux[1] c_hflux[2] c_hflux[3] v_tdfactor file fhlux_$NAME$_$TEMP$_$PRESS$.log\n"
			inscript_body[0][4] += "fix               tdacfout1 all ave/correlate 100 1000 1000000 c_hflux[1] c_hflux[1] type upper file tdacf1_$NAME$_$TEMP$_$PRESS$.log\n" 
			inscript_body[0][4] += "fix               tdacfout2 all ave/correlate 100 1000 1000000 c_hflux[2] c_hflux[2] type upper file tdacf2_$NAME$_$TEMP$_$PRESS$.log\n" 
			inscript_body[0][4] += "fix               tdacfout3 all ave/correlate 100 1000 1000000 c_hflux[3] c_hflux[3] type upper file tdacf3_$NAME$_$TEMP$_$PRESS$.log\n" 
			inscript_body[0][5] += "unfix             hfluxout\n" 
			inscript_body[0][5] += "unfix             tdacfout1\n"
			inscript_body[0][5] += "unfix             tdacfout2\n"
			inscript_body[0][5] += "unfix             tdacfout3\n"
		elif p == "st":
			inscript_body[1][1] += "variable          st equal 0.01*(lz/2.0)*(pzz-(pxx+pyy)/2.0)*1.01325\n"
			inscript_body[1][2] += "fix               stout all ave/time 100 1000 100000 v_st file st_$NAME$_$TEMP$.log\n"
			inscript_body[1][3] += "unfix             stout\n"
		elif p == "vle" or p == "ct" or p == "cd" or p == "cb" or p == "svp" or p == "nbp":
			if sur_flag == 0:
				inscript_body[1][1] += "variable          svp equal pzz\n"
				inscript_body[1][2] += "fix               svpout all ave/time 100 1000 100000 v_svp file svp_$NAME$_$TEMP$.log\n"
				inscript_body[1][2] += "fix               vleout all ave/spatial 1000 1000 1000000 z lower 1.0 density/mass ave running file vle_$NAME$_$TEMP$.log\n"
				inscript_body[1][3] += "unfix             svpout\n"
				inscript_body[1][3] += "unfix             vleout\n"
			else:
				pass
			sur_flag += 1
	if inscript_body[0][3] != "":
		inscript_body[0][2] += script_add("npt","bulk")
	inscript_body[0][3] += "unfix             npt\n"
	
	if inscript_body[0][5] != "":
		inscript_body[0][4] += "fix               nvt all nvt temp $TEMP$ $TEMP$ 100.0\n" 
		inscript_body[0][4] += script_add("nvt","bulk")
		inscript_body[0][5] += "unfix             nvt\n"
	
	if inscript_body[1][3] != "":
		inscript_body[1][2] += script_add("nvt","sur")
	inscript_body[1][3] += "unfix             nvt\n"

	return(inscript_body)

def creat_inscript_body(temp,bar,script1,script2,script3,script4):
	this_inscript_tmp = script1+script2+script3+script4 
	this_inscript_tmp = this_inscript_tmp.replace("$TEMP$",str(temp)) 
	this_inscript_tmp = this_inscript_tmp.replace("$PRESS$",str(bar))
	return(this_inscript_tmp)

def submit_job (computer,path1,shname,mole,model):
	if computer == "cluster":
		shscript = sh_cluster_script
	elif computer == "hpc":
		shscript = sh_hpc_script
	shscript = shscript.replace("%QUENCE%",quence)
	shscript = shscript.replace("%NPROC%",nproc)
	shscript = shscript.replace("%PATH%",path1)
	name = shname.strip(".sh")
	shscript = shscript.replace("%NAME%",name)
	shscript = shscript.replace("%MOLE%",mole)
	shscript = shscript.replace("%MODEL%",model)
	shf = open (shname,'w')
	shf.write(shscript)
	shf.close()
	if computer == "cluster":
		os.system("qsub %s"%shname)
		time.sleep(1)		
	elif computer == "hpc":
		os.system("bsub < ./%s"%shname)		
		time.sleep(1)		

def creat_mole_curve(iso,mole,mole_n,type,model,arr1,arr2,path1,script,ins,computer):
	if iso == "isobar":
		name1 = "bar"
		name2 = "K"
	elif iso == "isothermo":
		name1 = "K"
		name2 = "bar"
	for i in arr1: 
		if type == "s":
			this_inscript = script
			if model == "bulk":
				inname = mole + "_" + str(i) + name1 + "_" + model + ".in"
				shname = mole + "_" + str(i) + name1 + "_" + model + ".sh"
			elif model == "sur":
				inname = mole + "_" + model + ".in"
				shname = mole + "_" + model + ".sh"
			inoutf = open("%s/%s" %(path1,inname),'w')
			for j in arr2: 
				if iso == "isobar":
					this_inscript =this_inscript.replace("$TEMP$",str(j)) 
					this_inscript =this_inscript.replace("$PRESS$",str(i)) 
					this_inscript += creat_inscript_body(j,i,ins[2],ins[3],ins[4],ins[5])
				elif iso == "isothermo":
					this_inscript =this_inscript.replace("$TEMP$",str(i)) 
					this_inscript =this_inscript.replace("$PRESS$",str(j)) 
					this_inscript += creat_inscript_body(i,j,ins[2],ins[3],ins[4],ins[5])
			this_inscript = this_inscript.replace("$NAME$",mole)
			this_inscript = this_inscript.replace("$MOL_NUMBER$",mole_n)
			this_inscript = this_inscript.replace("$MODEL$",model)
			inoutf.write(this_inscript)
			inoutf.close()
			submit_job (computer,path1,shname,mole,model)
		elif type == "i":
			for j in arr2:
				this_inscript = script
				if iso == "isobar":
					if model == "bulk":
						inname = mole + "_" + str(j) + name2 + "_" + str(i) + name1 + "_" + model + ".in"
						shname = mole + "_" + str(j) + name2 + "_" + str(i) + name1 + "_" + model + ".sh"
					elif model == "sur":
						inname = mole + "_" + str(j) + name2 + "_" + model + ".in" 
						shname = mole + "_" + str(j) + name2 + "_" + model + ".sh" 
					this_inscript =this_inscript.replace("$TEMP$",str(j)) 
					this_inscript =this_inscript.replace("$PRESS$",str(i)) 
					this_inscript += creat_inscript_body(j,i,ins[2],ins[3],ins[4],ins[5])
				elif iso == "isothermo":
					inname = mole + "_" + str(i) + name1 + "_" + str(V) + name2 + "_" + model + ".in" 
					shname = mole + "_" + str(i) + name1 + "_" + str(V) + name2 + "_" + model + ".sh" 
					this_inscript =this_inscript.replace("$TEMP$",str(i)) 
					this_inscript =this_inscript.replace("$PRESS$",str(j)) 
					this_inscript += creat_inscript_body(i,j,ins[2],ins[3],ins[4],ins[5])
				this_inscript = this_inscript.replace("$NAME$",mole)
				this_inscript = this_inscript.replace("$MOL_NUMBER$",mole_n)
				this_inscript = this_inscript.replace("$MODEL$",model)
				inoutf = open("%s/%s" %(path1,inname),'w')
				inoutf.write(this_inscript)
				inoutf.close()
				submit_job (computer,path1,shname,mole,model)

def generate_job(type,temp,press,properties,moles,curve,computer):
	ins_tmp = creat_inscript(properties)	
	for mole in moles:
		os.system("mkdir %s"%mole)
		mole_path = JOB_HOME + "/" + mole

		inf = open("%s.data"%mole,'r')
		countline = 0
		for line in inf:
			countline += 1
			if line.find("Bonds") != -1:
				max_mole = countline - 2
				break
		inf.close()
		inf = open("%s.data"%mole,'r')
		countline = 0
		for line in inf:
			countline += 1
			if countline == max_mole:
				string = line.strip().split()
				mole_n = string[1]
				break
		inf.close()
	
		inscript = inscript_head
		if mole == "water":
			inscript += "fix                fshakespce all shake 0.0001 10 100 b 1 a 1\n"

		tem_arr = [val for val in NPT_properties if val in properties]
		if len(tem_arr) != 0:
			os.system("cp %s.data %s/%s_bulk.data"%(mole,mole_path,mole))
			bulk_inscript = inscript
			bulk_inscript += ins_tmp[0][0]	
			bulk_inscript += ins_tmp[0][1]	
			if curve == "isobar":
				creat_mole_curve(curve,mole,mole_n,type,"bulk",press,temp,mole_path,bulk_inscript,ins_tmp[0],computer)
			if curve == "isothermo":
				creat_mole_curve(curve,mole,mole_n,type,"bulk",temp,press,mole_path,bulk_inscript,ins_tmp[0],computer)

		tem_arr = [val for val in NVT_properties if val in properties]
		if len(tem_arr) != 0:
			outf = open("%s/%s_sur.data"%(mole_path,mole),'w')
			inf = open("%s.data"%mole,'r')
			for line in inf:
				if line.find("zlo") != -1:
					string = line.strip().split()
					z = float(string[1])*5.0
					outf.write("     %s    %.4f   zlo zhi\n" %(string[0],z))
				else:
					outf.write(line)
			inf.close()
			outf.close()
			sur_inscript = inscript
			sur_inscript += ins_tmp[1][0]
			sur_inscript += ins_tmp[1][1]
			creat_mole_curve("isobar",mole,mole_n,type,"sur",press,temp,mole_path,sur_inscript,ins_tmp[1],computer)
					

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print __doc__
    else:
	computer = "cluster"
	TYPE = "i"
	CURVE = "isobar"

	properties = []
	molecule = []
	temperatures = []
	pressures = []

	temp_start = 300.0
	temp_end = temp_start + 1.0
	temp_var = 10.0
	bar_start = 1.0
	bar_end = bar_start + 0.1
	bar_var = 1.0
	t_n = 0
	b_n = 0
	start = 1
	for para in sys.argv[start:]:
		if para == "-b":
			for bar_para in sys.argv[start+1:]:
				if bar_para.find("-") == -1:
					b_n = b_n + 1
				else:
					break
			if b_n == 1:
				CURVE = "isobar"
				bar_start = float(sys.argv[start+1])
				bar_end = bar_start + 0.1
				bar_var = 1.0
			elif b_n == 3:
				bar_start = float(sys.argv[start+1])
				bar_end = float(sys.argv[start+2])
				bar_var = float(sys.argv[start+3])
			else:
				sys.exit("parameters for pressure error: neithor one or three parameters")

		elif para == "-t":
			for temp_para in sys.argv[start+1:]:
				if temp_para.find("-")==-1:
					t_n =t_n + 1
				else:
					break
			if t_n == 1:
				CURVE = "isothermo"
				temp_start = float(sys.argv[start+1])
				temp_end = temp_start + 1.0
				temp_var = 10.0
			elif t_n == 3:
				temp_start = float(sys.argv[start+1])
				temp_end = float(sys.argv[start+2])
				temp_var = float(sys.argv[start+3])
			else:
				sys.exit("parameters for temperature error: neithor one or three parameters")
		elif para == "-p":
			for proper_para in sys.argv[start+1:]:
				if proper_para.find("-") == -1:
					properties.append(proper_para)
				else:
					break
		elif para == "-m":
			for mole_para in sys.argv[start+1:]:
				if mole_para.find("-") == -1:
					molecule.append(mole_para)
				else:
					break
		elif para == "-c":
			CURVE = sys.argv[start+1]
		
		elif para == "-j":
			TYPE = sys.argv[start+1]
			
		elif para == "-computer":
			computer = sys.argv[start+1]
		start += 1

	if properties[0] == "all":
		properties = all_properties
	
	if (t_n == 0) and (b_n == 0):
		CURVE == "isobar" 
	elif (t_n == 1) and (b_n == 1):
		CURVE == "isobar" 
	elif (t_n == 3) and (b_n == 3):
		if CURVE == "":
			sys.exit("you should assign the type of curve: isobar or isothermo") 

	if len(properties) == 0:
		sys.exit("no property would be calculated")
	if len(molecule) == 0:
		sys.exit("no molecule would be calculated")
	tmp_arr = list(set(properties).difference(set(all_properties)))
	if len(tmp_arr) != 0 and tmp_arr[0] != "all":
		sys.exit("some properties set error")

	
	temp_N = int((temp_end-temp_start)/temp_var)
	if temp_N < 0:
		sys.exit("temperature range error: end temperature must larger than start temperature or using minus interval value")
	for i in range(0,temp_N+1):
		V1 = temp_start + i*temp_var
		temperatures.append(V1)
	
	press_N = int((bar_end-bar_start)/bar_var)
	if press_N < 0:
		sys.exit("pressure range error: end press must larger than start press or using minus interval value")
	for i in range(0,press_N+1):
		V1 = bar_start + i*bar_var
		pressures.append(V1)
	
	outf = open("job-parameters.log",'w')
	outf.write("molecules ")
	for para in molecule:
		outf.write("%s "%para)
	outf.write("\n")

	outf.write("properties ")
	for para in properties:
		outf.write("%s "%para)
	outf.write("\n")

	outf.write("temperatures ")
	for para in temperatures:
		outf.write("%s "%para)
	outf.write("\n")
	
	outf.write("pressures ")
	for para in pressures:
		outf.write("%s "%para)
	outf.write("\n")
	
	outf.write("curve %s\n" % CURVE)
	outf.write("job_type %s\n" %TYPE)
	outf.write("computer %s\n" %computer)

	generate_job(TYPE,temperatures,pressures,properties,molecule,CURVE,computer)
