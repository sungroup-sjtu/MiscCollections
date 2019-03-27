#!/bin/env python

from sys import argv

import os

cgrepeat = 30
REPEAT = 302
iatom  = 10
midatom = range(iatom+2,REPEAT-iatom)
RULES = [tuple(range(1,iatom+2))] + [tuple(midatom[i:i+iatom]) for i in range(0,len(midatom),iatom)] + [tuple(range(REPEAT-iatom,REPEAT+1))]

mass = (12,12,1,1,12,16,1,1,1,1)
MASS = [(12,12,1,1,12,16,1,1,1,1,1)] + [mass[:] for i in range(1,cgrepeat-1)] + [(1,12,12,1,1,12,16,1,1,1,1)]

FORMAT = "%.4f"

def read_frame(f):
	if len(ATOMS) > 0:
		del ATOMS[:len(ATOMS)]
	if len(BOXLINES) > 0:
		del BOXLINES[:len(BOXLINES)]

	nAtoms = 0

	line = f.readline()
	while line:
		if line.startswith("Generated"):
			line = f.readline()
			nAtoms = int(line.strip())
		else:
			count = 1
			while count <= nAtoms:
				strs = line.strip().split()
				length = len(strs)
				if length == 6:
					ATOMS.append((float(strs[3])*10, float(strs[4])*10, float(strs[5])*10))
				if length == 5:
					ATOMS.append((float(strs[2])*10, float(strs[3])*10, float(strs[4])*10))
				count = count + 1
				line = f.readline()

			strs = line.strip().split()
			BOXLINES.append(float(strs[0])*10)
			BOXLINES.append(float(strs[1])*10)
			BOXLINES.append(float(strs[2])*10)

			return True

		line = f.readline()

	return False

def convert(out):
	for i in range(len(ATOMS) / REPEAT):
		for j in range(len(RULES)):
			t_x = 0.0
			t_y = 0.0
			t_z = 0.0
			t_mass = 0.0
			for k in range(len(RULES[j])):
				t_x = t_x + ATOMS[RULES[j][k]-1+i*REPEAT][0] * MASS[j][k]
				t_y = t_y + ATOMS[RULES[j][k]-1+i*REPEAT][1] * MASS[j][k]
				t_z = t_z + ATOMS[RULES[j][k]-1+i*REPEAT][2] * MASS[j][k]
				t_mass = t_mass + MASS[j][k]

			t_x = float(t_x) / t_mass
			t_y = float(t_y) / t_mass
			t_z = float(t_z) / t_mass
			f_x = FORMAT % t_x
			f_y = FORMAT % t_y
			f_z = FORMAT % t_z

			out.write(str(len(RULES)*i+j+1) + "	  " + str(j+1) + "	" + f_x + "  " + f_y + "  " + f_z + "\n")


time = 0
f = open(argv[1], 'r')
o = open("PDMS54-CG.lammpstrj", 'w')

ATOMS = []
BOXLINES = []
while read_frame(f):
	o.write("ITEM: TIMESTEP\n")
	o.write(str(time) + "\n")
	o.write("ITEM: NUMBER OF ATOMS\n")
	o.write(str(len(ATOMS) / REPEAT * len(RULES)) + "\n")
	o.write("ITEM: BOX BOUNDS pp pp pp\n")
	o.write("0	" + str(BOXLINES[0]) + "\n")
	o.write("0	" + str(BOXLINES[1]) + "\n")
	o.write("0	" + str(BOXLINES[2]) + "\n")
	o.write("ITEM: ATOMS id type xu yu zu\n")

	convert(o)

	time = time + 1000

f.close()
