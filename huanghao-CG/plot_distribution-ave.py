#!/usr/bin/env python
#encoding=utf-8

"""
[USAGE]
plot_intra_distribution_ave.py [trjfile] [molecule_particles] [nbins] [MODULE:BOND|ANGLE|DIHEDRAL|OOPA] [OUTPUT]

This script is used to calculate histogram of CG trajectory by assigning atom id, and compared with AA distribution.
"""

import os, sys, math
import numpy as np
import matplotlib.pyplot as plt


if sys.argv[1] == "-h" or sys.argv[1] == "h" or sys.argv[1] == "help" or sys.argv[1] == "-help":
    print __doc__
    sys.exit()

try:
    N_PARTICLES = int(sys.argv[2])
    N_BINS = int(sys.argv[3]) + 1
    MODULE = sys.argv[4].strip().upper()
    AARESULT = sys.argv[5]
    ATYPE = None
    BTYPE = None
    CTYPE = None
    DTYPE = None
#    ATYPE = int(sys.argv[5])
#    BTYPE = int(sys.argv[6])
#    CTYPE = 0
#    DTYPE = 0

except:
    print __doc__
    sys.exit()

nATOMS = 0

def read_frame(f):
    global nATOMS

    line = f.readline()
    while line:
        if line.startswith("ITEM: NUMBER"):
            line = f.readline()
            nATOMS = int(line.strip())
        elif line.startswith("ITEM: ATOMS"):
            count = 1
            while count <= nATOMS:
                repeat = 1
                aid_vec = []
                bid_vec = []
                cid_vec = []
                did_vec = []
                while repeat <= N_PARTICLES:
                    line = f.readline()
                    strs = line.strip().split()
                    if int(strs[0])%N_PARTICLES == ATYPE:
                        aid_vec.append((float(strs[2]), float(strs[3]), float(strs[4])))
                    elif int(strs[0])%N_PARTICLES == BTYPE:
                        bid_vec.append((float(strs[2]), float(strs[3]), float(strs[4])))
                    elif MODULE == "ANGLE" and int(strs[0])%N_PARTICLES == CTYPE:
                        cid_vec.append((float(strs[2]), float(strs[3]), float(strs[4])))
                    elif MODULE == "DIHEDRAL" or MODULE == "OOPA":
                        if int(strs[0])%N_PARTICLES == CTYPE:
                            cid_vec.append((float(strs[2]), float(strs[3]), float(strs[4])))
                        elif int(strs[0])%N_PARTICLES == DTYPE:
                            did_vec.append((float(strs[2]), float(strs[3]), float(strs[4])))

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
        output_raw = []
        outfile = ""
        has_input_before = False
        while True:
            input_str = raw_input()
            input_split = map(int, input_str.strip().split())

            if input_split[0] == 0:
                break

            DIS_VEC = []
            ATYPE = input_split[0]
            BTYPE = input_split[1]

            if MODULE in ["ANGLE", "DIHEDRAL", "OOPA"]:
                CTYPE = input_split[2]
            if MODULE in ["DIHEDRAL", "OOPA"]:
                DTYPE = input_split[3]

            trjfile = open(sys.argv[1], 'r')
            frames = 0
            while read_frame(trjfile):
                frames = frames + 1

            trjfile.close()

            if not has_input_before:
                if MODULE == "BOND":
                    outfile = str(ATYPE) + "-" + str(BTYPE) + "." + MODULE
                elif MODULE == "ANGLE":
                    outfile = str(ATYPE) + "-" + str(BTYPE) + "-" + str(CTYPE) + "." + MODULE
                elif MODULE == "DIHEDRAL":
                    outfile = str(ATYPE) + "-" + str(BTYPE) + "-" + str(CTYPE) + "-" + str(DTYPE) + "." + MODULE
                elif MODULE == "OOPA":
                    outfile = str(ATYPE) + "-" + str(BTYPE) + "-" + str(CTYPE) + "-" + str(DTYPE) + "." + MODULE

                has_input_before = True

            output_raw += getDistribution(frames)

        output_raw = np.array(output_raw)

        # resample
        hist, edges = np.histogram(output_raw[:,0], bins=N_BINS)
        edge_begin = edges[0]
        edge_step = edges[1] - edges[0]
        
        output = np.zeros((N_BINS, 2))
        output[:,0] = (edges + edge_step/2.0)[:-1]
        for row in output_raw:
            idx = np.floor((row[0] - edge_begin)/edge_step)
            try:
                output[idx,1] += row[1] / hist[idx]
            except IndexError:
                print idx

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
