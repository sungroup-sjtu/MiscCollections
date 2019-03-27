#!/usr/bin/env python
#encoding=utf-8

"""
[Usage] python PolymerCG_Bonded.py [MD trajectory file] [CG_INPUT] [nbins] [BOND/ANGLE/DIHEDRAL] [CG_id1] [CG_id2] [CG_id3] [CG_id4] [AA results]
[CG_INPUT] CG rule definition
"""

import os, sys, math
import numpy as np
import matplotlib.pyplot as plt

if sys.argv[1] == "-h" or sys.argv[1] == "h" or sys.argv[1] == "help" or sys.argv[1] == "-help":
	print __doc__
	sys.exit()
	
FORMAT = "%.4f"	

try:
	N_BINS = int(sys.argv[3])+1
	MODULE = sys.argv[4].strip().upper()
	ATYPE = int(sys.argv[5])
	BTYPE = int(sys.argv[6])
	CTYPE = 0
	DTYPE = 0
	AARESULT = ""

	if MODULE == "DIHEDRAL":
		CTYPE = int(sys.argv[7])
		DTYPE = int(sys.argv[8])
		AARESULT = sys.argv[9]
	elif MODULE == "OOPA":
		CTYPE = int(sys.argv[7])
		DTYPE = int(sys.argv[8])
		AARESULT = sys.argv[9]
	elif MODULE == "ANGLE":
		CTYPE = int(sys.argv[7])
		AARESULT = sys.argv[8]
	elif MODULE == "BOND":
		AARESULT = sys.argv[7]
	else:
		print "UNKOWN MODULE TYPE"
		sys.exit()
	print "[N_BINS]", N_BINS, "[MODULE]",MODULE 
	print "input parameters are correct."
except:
	print __doc__
	sys.exit()

nATOMS = 0
	
def CG_Read(CGIn):
	print "starting reading CG rule input"
	global nMono, nMole
	global a1, a2, a3, a4

	line = CGIn.readline()
	while line:
		if line.startswith("Atoms_per_monomer"):
			tps = line.strip().split()
			nMono = int(tps[1])  #Atoms_per_monomer
		if line.startswith("Atoms_per_molecule"):
			tps = line.strip().split()
			nMole = int(tps[1])  #Atoms_per_molecule
		if line.startswith("CG_Bead_ID"):
			count = 0
			cgRule = []
			while count < nMono:
				line = CGIn.readline()
				tps = line.strip().split()
				cgRule.append(((int(tps[0])),(int(tps[1])),(str(tps[2])),(int(tps[3])),(int(tps[4])),(str(tps[5]))))
				count +=1;
		line = CGIn.readline()

	a1=[]; a2=[]; a3=[]; a4=[]
	i = 0; id = cgRule[0][1]
	while (i < nMono) and (id < nMole-nMono-1):
		if MODULE == "BOND":
			a1.append((cgRule[i][0],cgRule[i][1],cgRule[i][2],cgRule[i][3],cgRule[i][4],cgRule[i][5]))
			a2.append((cgRule[i][0],cgRule[i][1]+nMono,cgRule[i][2],cgRule[i][3],cgRule[i][4],cgRule[i][5]))
			id=cgRule[i][1]+nMono 
			i+=1 
		elif MODULE == "ANGLE":
			a1.append((cgRule[i][0],cgRule[i][1],cgRule[i][2],cgRule[i][3],cgRule[i][4],cgRule[i][5]))
			a2.append((cgRule[i][0],cgRule[i][1]+nMono,cgRule[i][2],cgRule[i][3],cgRule[i][4],cgRule[i][5]))
			a3.append((cgRule[i][0],cgRule[i][1]+nMono*2,cgRule[i][2],cgRule[i][3],cgRule[i][4],cgRule[i][5]))
			id=cgRule[i][1]+nMono*2 
			i+=1
		elif MODULE == "DIHEDRAL":
			a1.append((cgRule[i][0],cgRule[i][1],cgRule[i][2],cgRule[i][3],cgRule[i][4],cgRule[i][5]))
			a2.append((cgRule[i][0],cgRule[i][1]+nMono,cgRule[i][2],cgRule[i][3],cgRule[i][4],cgRule[i][5]))
			a3.append((cgRule[i][0],cgRule[i][1]+nMono*2,cgRule[i][2],cgRule[i][3],cgRule[i][4],cgRule[i][5]))
			a4.append((cgRule[i][0],cgRule[i][1]+nMono*3,cgRule[i][2],cgRule[i][3],cgRule[i][4],cgRule[i][5]))
			id=cgRule[i][1]+nMono*3
			i+=1

	print "CG rule input file reading finished"

				
def CoM(AA_Coor): #calculating Center of Mass of a CG bead from the coordinates of atoms grouped
	tx=0.0; ty=0.0; tz=0.0; tmass=0.0
	CG_Coor=[]
	for aa in AA_Coor:
		tx = tx + aa[0]*aa[3]
		ty = ty + aa[1]*aa[3]
	   	tz = tz + aa[2]*aa[3]
		tmass = tmass + aa[3]
	tx = FORMAT % (float(tx)/tmass)
	ty = FORMAT % (float(ty)/tmass)
	tz = FORMAT % (float(tz)/tmass)
	CG_Coor.append(tx); CG_Coor.append(ty);CG_Coor.append(tz)
	return CG_Coor
	

def read_frame(f):
	global nATOMS
	global nMono, nMole
	global a1, a2, a3, a4

	line = f.readline()
	while line:
		if line.startswith("ITEM: NUMBER"):
			line = f.readline()
			nATOMS = int(line.strip())
		elif line.startswith("ITEM: ATOMS"):
			count = 1; iMole = 0
			while count <= nATOMS:
				repeat = 1
				aid_vec = [];bid_vec = [];cid_vec = [];did_vec = []
				AA_Coor1 = [];AA_Coor2 = [];AA_Coor3 = [];AA_Coor4 = []
				CG_Coor1 = [];CG_Coor2 = [];CG_Coor3 = [];CG_Coor4 = []
				temp1 = []; temp2 = []; temp3 = []
				iMono = 0 
				while repeat <= nMole:
					
					line = f.readline()
					strs = line.strip().split()
					ii=0	
					if MODULE == "BOND":
						while ii < nMono:
							if int(strs[0])-iMole*nMole == int(a1[ii][1])+iMono*nMono:
								AA_Coor1.append((float(strs[2]), float(strs[3]), float(strs[4]), float(a1[0][4])))
								break
							if int(strs[0])-iMole*nMole == int(a2[ii][1])+iMono*nMono:
								AA_Coor2.append((float(strs[2]), float(strs[3]), float(strs[4]), float(a2[0][4])))
								break
							ii+=1	
						
						if (len(AA_Coor1)==nMono and len(AA_Coor2)==nMono):
							if not temp1: 				 # if temp1 is empty	
								CG_Coor1 = CoM(AA_Coor1) # calculating the centre of mass of a group of atoms
								CG_Coor2 = CoM(AA_Coor2)
							else:	
								CG_Coor1 = temp1
								CG_Coor2 = CoM(AA_Coor1)
							aid_vec.append((float(CG_Coor1[0]), float(CG_Coor1[1]), float(CG_Coor1[2])))
							bid_vec.append((float(CG_Coor2[0]), float(CG_Coor2[1]), float(CG_Coor2[2])))
							temp1 = CG_Coor2
							CG_Coor1 = [];CG_Coor2 = []
							iMono = iMono + 1
						

							
					if MODULE == "ANGLE":
						while ii < nMono:
						#	print iMole, iMono, ii, a1[ii][1],strs[0]
							if int(strs[0])-iMole*nMole == int(a1[ii][1]+iMono*nMono):
								AA_Coor1.append((float(strs[2]), float(strs[3]), float(strs[4]), float(a1[0][4]))); break
							if int(strs[0])-iMole*nMole == int(a2[ii][1]+iMono*nMono):
								AA_Coor2.append((float(strs[2]), float(strs[3]), float(strs[4]), float(a2[0][4]))); break
							if int(strs[0])-iMole*nMole == int(a3[ii][1]+iMono*nMono):
								AA_Coor3.append((float(strs[2]), float(strs[3]), float(strs[4]), float(a3[0][4]))); break
							ii+=1
							
						if 	(len(AA_Coor1)==nMono and len(AA_Coor2)==nMono and len(AA_Coor3)==nMono):
							if (not temp1) and (not temp2):
								CG_Coor1 = CoM(AA_Coor1)
								CG_Coor2 = CoM(AA_Coor2)
								CG_Coor3 = CoM(AA_Coor3)
							else:
								CG_Coor1 = temp1
								CG_Coor2 = temp2
								CG_Coor3 = CoM(AA_Coor1)
							aid_vec.append((float(CG_Coor1[0]), float(CG_Coor1[1]), float(CG_Coor1[2])))
							bid_vec.append((float(CG_Coor2[0]), float(CG_Coor2[1]), float(CG_Coor2[2])))
							cid_vec.append((float(CG_Coor3[0]), float(CG_Coor3[1]), float(CG_Coor3[2])))
							temp1 = CG_Coor2; temp2 = CG_Coor3
							CG_Coor1 = [];CG_Coor2 = [];CG_Coor3 = []
							iMono = iMono + 1
							
					if MODULE == "DIHEDRAL":
						while ii < nMono:
							if int(strs[0])-iMole*nMole == int(a1[ii][1]+iMono*nMono):
								AA_Coor1.append((float(strs[2]), float(strs[3]), float(strs[4]), float(a1[0][4])));break
							if int(strs[0])-iMole*nMole == int(a2[ii][1]+iMono*nMono):
								AA_Coor2.append((float(strs[2]), float(strs[3]), float(strs[4]), float(a2[0][4]))); break
							if int(strs[0])-iMole*nMole == int(a3[ii][1]+iMono*nMono):
								AA_Coor3.append((float(strs[2]), float(strs[3]), float(strs[4]), float(a3[0][4]))); break
							if int(strs[0])-iMole*nMole == int(a4[ii][1]+iMono*nMono):
								AA_Coor4.append((float(strs[2]), float(strs[3]), float(strs[4]), float(a4[0][4]))); break
							ii+=1
						
						if 	(len(AA_Coor1)==nMono and len(AA_Coor2)==nMono and len(AA_Coor3)==nMono and len(AA_Coor4)==nMono):						
							if (not temp1) and (not temp2) and (not temp3):
								CG_Coor1 = CoM(AA_Coor1)
								CG_Coor2 = CoM(AA_Coor2)
								CG_Coor3 = CoM(AA_Coor3)
								CG_Coor4 = CoM(AA_Coor4)
							else:
								CG_Coor1 = temp1
								CG_Coor2 = temp2
								CG_Coor3 = temp3
								CG_Coor4 = CoM(AA_Coor1)
							aid_vec.append((float(CG_Coor1[0]), float(CG_Coor1[1]), float(CG_Coor1[2])))
							bid_vec.append((float(CG_Coor2[0]), float(CG_Coor2[1]), float(CG_Coor2[2])))
							cid_vec.append((float(CG_Coor3[0]), float(CG_Coor3[1]), float(CG_Coor3[2])))
							did_vec.append((float(CG_Coor4[0]), float(CG_Coor4[1]), float(CG_Coor4[2])))
							temp1 = CG_Coor2; temp2 = CG_Coor3; temp3 = CG_Coor4
							CG_Coor1 = [];CG_Coor2 = [];CG_Coor3 = [];CG_Coor4 = []
							iMono = iMono + 1
					

					repeat = repeat + 1
					count = count + 1
				
				if MODULE == "BOND":
					calcBond(aid_vec, bid_vec)
				elif MODULE == "ANGLE":
					calcAngle(aid_vec, bid_vec, cid_vec)
				elif MODULE == "DIHEDRAL":
					calcDihedral(aid_vec, bid_vec, cid_vec, did_vec)
				elif MODULE == "OOPA":
					calcOOPA(aid_vec, bid_vec, cid_vec, did_vec)
				
				aid_vec = []
				bid_vec = []
				cid_vec = []
				did_vec = []
				iMole = iMole + 1
			return True

		line = f.readline()

	return False


DIS_VEC = []
def calcBond(avec, bvec):
	for va in avec:
		for vb in bvec:
			dis = (va[0] - vb[0])**2 + (va[1] - vb[1])**2 + (va[2] - vb[2])**2
			dis = math.sqrt(dis)
			DIS_VEC.append(dis)

def calcAngle(avec, bvec, cvec):
	for va in avec:
		for vb in bvec:
			for vc in cvec:
				w1 = (va[0]-vb[0], va[1]-vb[1], va[2]-vb[2])
				w2 = (vc[0]-vb[0], vc[1]-vb[1], vc[2]-vb[2])
				w1len = math.sqrt(w1[0]**2 + w1[1]**2 + w1[2]**2)
				w2len = math.sqrt(w2[0]**2 + w2[1]**2 + w2[2]**2)
				dot = w1[0]*w2[0] + w1[1]*w2[1] + w1[2]*w2[2]
				cosphi = dot / w1len / w2len
				if cosphi > 1.0:
					cosphi = 1.0
				elif cosphi < -1.0:
					cosphi = -1.0

				DIS_VEC.append(math.acos(cosphi))

def calcDihedral(avec, bvec, cvec, dvec):
	for va in avec:
		for vb in bvec:
			for vc in cvec:
				for vd in dvec:
					w1 = cross3p(va, vb, vc)
					w2 = cross3p(vb, vc, vd)
					w1len = math.sqrt(dotp(w1, w1))
					w2len = math.sqrt(dotp(w2, w2))
					w3 = triProduct(va, vb, vc)
					
					sign = dotp(w3, w2)
					isn = -1
					if sign > 0:
						isn = -1
					else:
						isn = 1

					cosphi = dotp(w1, w2) / w1len / w2len
					if cosphi > 1.0:
						cosphi = 1.0
					elif cosphi < -1.0:
						cosphi = -1.0

					DIS_VEC.append(isn*math.acos(cosphi))

def calcOOPA(avec, bvec, cvec, dvec):
	for va in avec:
		for vb in bvec:
			for vc in cvec:
				for vd in dvec:
					w1 = (va[0]-vb[0], va[1]-vb[1], va[2]-vb[2])
					w2 = cross3p(vc, vb, vd)
					w1len = math.sqrt(dotp(w1, w1))
					w2len = math.sqrt(dotp(w2, w2))
					sinphi = dotp(w1, w2) / w1len / w2len
					if sinphi > 1.0:
						sinphi = 1.0
					elif sinphi < -1.0:
						sinphi = -1.0

					DIS_VEC.append(math.asin(sinphi))

def dotp(va, vb):
	d = va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2]
	return d

def cross3p(va, vb, vc):
	v1 = (va[0]-vb[0], va[1]-vb[1], va[2]-vb[2])
	v2 = (vc[0]-vb[0], vc[1]-vb[1], vc[2]-vb[2])
	w = (v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0])
	return w

def triProduct(r1, r2, r3):
	x = (r1[0]-r2[0], r1[1]-r2[1], r1[2]-r2[2])
	y = (r3[0]-r2[0], r3[1]-r2[1], r3[2]-r2[2])
	a = x[0]*y[0] + x[1]*y[1] + x[2]*y[2]
	b = y[0]**2 + y[1]**2 + y[2]**2
	w = (a*y[0]-b*x[0], a*y[1]-b*x[1], a*y[2]-b*x[2])
	return w

def getDistribution(frames):
	DIS_VEC.sort()
	
	output = []
	distribution = [0 for i in range(N_BINS)]
	if MODULE == "BOND":
		min = DIS_VEC[0] - 0.05
		max = DIS_VEC[len(DIS_VEC)-1] + 0.05
		interval = float(max - min) / N_BINS

		for dis in DIS_VEC:
			location = int(math.ceil((dis - min) / interval))
			distribution[location - 1] = distribution[location - 1] + 1
	
		x = []
		y = []
		for i in range(len(distribution)):
			x.append(min + i * interval)
			y.append(float(distribution[i]))

		integrate = np.trapz(y, x)

		for i in range(len(distribution)):
			a = x[i]
			b = y[i] / integrate
			output.append((a, b))
		
		return output
	elif MODULE == "ANGLE":
		min = DIS_VEC[0] - 0.01
		max = DIS_VEC[len(DIS_VEC)-1] + 0.01
		interval = float(max - min) / N_BINS

		for dis in DIS_VEC:
			location = int(math.ceil((dis - min) / interval))
			distribution[location - 1] = distribution[location - 1] + 1

		x = []
		y = []
		for i in range(len(distribution)):
			x.append(min + i * interval)
			y.append(float(distribution[i]))

		integrate = np.trapz(y, x)

		for i in range(len(distribution)):
			a = x[i] * 180 / math.pi
			b = y[i] / integrate
			output.append((a, b))
		
		return output
	elif MODULE == "DIHEDRAL":
		min = DIS_VEC[0] - 0.01
		max = DIS_VEC[len(DIS_VEC)-1] + 0.01
		interval = float(max - min) / N_BINS

		for dis in DIS_VEC:
			location = int(math.ceil((dis - min) / interval))
			distribution[location - 1] = distribution[location - 1] + 1

		x = []
		y = []
		for i in range(len(distribution)):
			x.append(min + i * interval)
			y.append(float(distribution[i]))

		integrate = np.trapz(y, x)

		for i in range(len(distribution)):
			a = x[i] * 180 / math.pi
			b = y[i] / integrate
			output.append((a, b))
		
		return output
	elif MODULE == "OOPA":
		min = DIS_VEC[0] - 0.01
		max = DIS_VEC[len(DIS_VEC)-1] + 0.01
		interval = float(max - min) / N_BINS

		for dis in DIS_VEC:
			location = int(math.ceil((dis - min) / interval))
			distribution[location - 1] = distribution[location - 1] + 1

		x = []
		y = []
		for i in range(len(distribution)):
			x.append(min + i * interval)
			y.append(float(distribution[i]))

		integrate = np.trapz(y, x)

		for i in range(len(distribution)):
			a = x[i] * 180 / math.pi
			b = y[i] / integrate
			output.append((a, b))
		
		return output

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit()
	else:
		CGIn = open(sys.argv[2],'r')
		CG_Read(CGIn); CGIn.close()

		trjIn = open(sys.argv[1],'r')
		frames = 0
		while read_frame(trjIn):
			frames += 1
		trjIn.close()
		print frames, "Frames is read from the trajectory file."
		
		output = getDistribution(frames)
		
		outfile = ""
		if MODULE == "BOND":
			outfile = str(ATYPE) + "-" + str(BTYPE) + "." + MODULE
		elif MODULE == "ANGLE":
			outfile = str(ATYPE) + "-" + str(BTYPE) + "-" + str(CTYPE) + "." + MODULE
		elif MODULE == "DIHEDRAL":
			outfile = str(ATYPE) + "-" + str(BTYPE) + "-" + str(CTYPE) + "-" + str(DTYPE) + "." + MODULE
		elif MODULE == "OOPA":
			outfile = str(ATYPE) + "-" + str(BTYPE) + "-" + str(CTYPE) + "-" + str(DTYPE) + "." + MODULE

		site_cg = 0
		site_aa = 0
		max_pro = 0

		writeDis = open(outfile, 'w')
		x_cg = []
		y_cg = []
		for (i, j) in output:
			x_cg.append(i)
			y_cg.append(j)
			writeDis.write(str(i) + "   " + str(j) + "\n")
			if j > max_pro:
				site_cg = i
				max_pro = j

		writeDis.close()

		max_pro = 0
		aafile = open(AARESULT, 'r')
		x_aa = []
		y_aa = []
		for line in aafile:
			i = 0
			j = 0
			strs = line.strip().split()
			if MODULE == "BOND":
				i = float(strs[0])
				j = float(strs[1])
			elif MODULE == "ANGLE":
				i = float(strs[0]) 
				j = float(strs[1])
			elif MODULE == "DIHEDRAL":
				i = float(strs[0]) 
				j = float(strs[1])
			elif MODULE == "OOPA":
				i = float(strs[0])
				j = float(strs[1])
				
			x_aa.append(i)
			y_aa.append(j)

			if j > max_pro:
				site_aa = i
				max_pro = j

		aafile.close()

		print "Max AA site = ", site_aa
		print "Max CG site = ", site_cg

		## plot
		plt.figure()
		plt.plot(x_aa, y_aa, label="AA")
		plt.plot(x_cg, y_cg, label="CG")
		plt.legend()
		plt.show()
