#!/bin/env python

from sys import argv
import os

REPEAT = 220
RULES = [(16,15,17,20),(14,13,15,19),(12,18,13,11),(7,6,8,18),\
(1,2,8,9),(3,2,4,20),(10,9,11,17),(5,4,6,19),\
(21,54,55,60,85,86,155,156,157,167,168,169),(22,25,26,31,87,88,94,95,105,106,107),\
(23,32,33,38,89,90,108,109,110,120,121,122),(24,39,40,45,91,92,124,125,135,136,137),\
(69,70,71,76,187,188,189,190,191,201,202,203),(46,47,48,53,138,139,140,141,142,152,153,154),\
(61,62,63,68,170,171,172,173,174,184,185,186),(77,78,79,84,204,205,206,207,208,218,219,220),\
(56,57,58,59,158,159,160,161,162,163,164,165,166),(27,28,29,30,96,97,98,99,100,101,102,103,104),\
(34,35,36,37,111,112,113,114,115,116,117,118,119),(41,42,43,44,126,127,128,129,130,131,132,133,134),\
(72,73,74,75,192,193,194,195,196,197,198,199,200),(49,50,51,52,143,144,145,146,147,148,149,150,151),\
(64,65,66,67,175,176,177,178,179,180,181,182,183),(79,80,81,82,207,208,209,210,211,212,213,214,215)]

MASS = [(28,8,8,8)]*8 + [(12,12,12,12,1,1,1,1,1,1,1,1)]*8 + [(12,12,12,12,1,1,1,1,1,1,1,1,1)]*8

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
	l=0
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
			
			if int(j+1)%24 >= 1 and int(j+1)%24 <= 8:
				l=1
			elif int(j+1)%24 >= 9 and int(j+1)%24 <= 16:
				l=2
			else:
				l=3
			out.write(str(len(RULES)*i+j+1) + "	  " + str(l) + "	" + f_x + "  " + f_y + "  " + f_z + "\n")


time = 0
f = open(argv[1], 'r')
o = open("POSS-CG.lammpstrj", 'w')

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