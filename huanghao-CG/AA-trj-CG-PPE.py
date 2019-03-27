#!/bin/env python

from sys import argv
import os

REPEAT = 865
RULES = [(18,19,20,21,22,23,24,25),(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),(26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,1)] + \
[tuple(range(i,i+17)+[i-17]) for i in range(43,825,17)] + \
[(825,826,827,828,829,830,831,832,833,834,835,836,837,838,839,840,841,808), \
(842,843,844,845,846,847,848,849,850,851,852,853,854,855,856,857,858,825),(842,859,860,861,862,863,864,865)]

mass = (8,12,12,12,12,12,12,12,12,1,1,1,1,1,1,1,1,8)
MASS = [(8,12,12,1,1,1,1,1)] + [mass[:] for i in range(1,51)] + [(8,12,12,1,1,1,1,1)]


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

			if int(j+1)%52 == 1 or int(j+1)%52 == 0:
				l=1
			else:
				l=2
			out.write(str(len(RULES)*i+j+1) + "	  " + str(l) + "	" + f_x + "  " + f_y + "  " + f_z + "\n")

time = 0
f = open(argv[1], 'r')
o = open("PPE-CG.lammpstrj", 'w')

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
