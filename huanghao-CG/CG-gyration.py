#!/bin/env python

from sys import argv
import os

MASS=[2,1,2]
RULES=range(3)

def read_frame(f):
	line = f.readline()
	while line:
		if line.startswith(ITEM: TIMESTEP):
			line = f.readline()
			line = f.readline()
			line = f.readline()
			line = f.readline()
			line = f.readline()
			line = f.readline()
			line = f.readline()
			line = f.readline()
		else:
			while Ture:
				line = f.readline()
				strs = line.strip().split()
				length = len(strs)
				if length == 5:
				ATOMS.append((float(strs[2]), float(strs[3]), float(strs[4])))
				print ATOMS
			return True
		line = f.readline()
	return False

def convert(out):
	for i in range(len(ATOMS)/len(RULES)):
		t_x = 0.0
		t_y = 0.0
		t_z = 0.0
		g_x = 0.0
		g_y = 0.0
		g_z = 0.0
		g = 0.0 
		t_mass = 0.0
		for j in range(len(RULES)):
			t_x = t_x + ATOMS[RULES[j]+i*len(RULES)][0] * MASS[j]
			t_y = t_y + ATOMS[RULES[j]+i*len(RULES)][1] * MASS[j]
			t_z = t_z + ATOMS[RULES[j]+i*len(RULES)][2] * MASS[j]
			t_mass = t_mass + MASS[j]
		t_x = float(t_x) / t_mass
		t_y = float(t_y) / t_mass
		t_z = float(t_z) / t_mass
		for k in range(len(RULES)):
			g_x = g_x + MASS[k]*(ATOMS[RULES[k]+i*len(RULES)][0]-t_x)**2
			g_y = g_y + MASS[k]*(ATOMS[RULES[k]+i*len(RULES)][1]-t_y)**2
			g_z = g_z + MASS[k]*(ATOMS[RULES[k]+i*len(RULES)][2]-t_z)**2
		g = ((g_x + g_y + g_z) / t_mass)**0.5
		f_g = '%.4f' %g
		out.write(f_g + "\n")

f = open(argv[1], 'r')
o = open("gyration", 'w')

ATOMS = []
read_frame(f)
convert(o)
f.close()
