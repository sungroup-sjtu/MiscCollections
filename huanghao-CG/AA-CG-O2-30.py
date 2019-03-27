#!/usr/bin/env python

import os, sys

FORMAT   = "%.4f"

def write_numbers(dataFile, Natom, Nbond, Nangle, Ndihedral, Nimproper):
    dataFile.write("LAMMPS data file\n\n")
    dataFile.write(str(Natom) + " atoms\n")
    dataFile.write(str(Nbond) + " bonds\n")
    dataFile.write(str(Nangle) + " angles\n")
    dataFile.write(str(Ndihedral) + " dihedrals\n")
    dataFile.write(str(Nimproper) + " impropers\n\n")

def write_types(dataFile, Natom, Nbond, Nangle, Ndihedral, Nimproper):
    dataFile.write(str(Natom) + " atom types\n")
    dataFile.write(str(Nbond) + " bond types\n")
    dataFile.write(str(Nangle) + " angle types\n")
    dataFile.write(str(Ndihedral) + " dihedral types\n")
    dataFile.write(str(Nimproper) + " improper types\n\n")

def write_pbc(dataFile, PBC_x, PBC_y, PBC_z):
#    pbc_x = FORMAT % float(str(x))
#    pbc_y = FORMAT % float(str(y))
#    pbc_z = FORMAT % float(str(z))
#    dataFile.write("     0.0000     " + pbc_x + "   xlo xhii\n")
#    dataFile.write("     0.0000     " + pbc_y + "   ylo yhi\n")
#    dataFile.write("     0.0000     " + pbc_z + "   zlo zhi\n\n")
    dataFile.write(PBC_x)
#    dataFile.write("\n")
    dataFile.write(PBC_y)
#    dataFile.write("\n")
    dataFile.write(PBC_z)
    dataFile.write("\n")

def write_mass(dataFile):
    dataFile.write("Masses\n\n")
    dataFile.write("          1     44   #   CCO\n")
#    dataFile.write("          2     74.15420   #   SiM\n")

def write_atoms(dataFile, NUM):
    dataFile.write("Atoms\n\n")
    mid = 1
    i   = 1
    ii = 1
    tid = 1
    type = "    #  CCO"
    while i <= len(NUM):
#        if ii == 1:
#            tid = 1
#            type = "    #  SiT"
#        elif ii > 1 and ii <= 19:
#            tid = 2
#            type = "    #  SiM"
#        else:
#            tid = 1
#            type = "    #  SiT"
        type = "    #  CCO"
        x = FORMAT % float(str(NUM[i-1][0]))
        y = FORMAT % float(str(NUM[i-1][1]))
        z = FORMAT % float(str(NUM[i-1][2]))
#        if ii == 5:
#        dataFile.write("          " + str(i) + "        " + str(mid) + "        " + str(tid) + "    -1.0000     " + x + "     " + y + "     " + z + type + "\n")
#        else:
        dataFile.write("          " + str(i) + "        " + str(mid) + "        " + str(tid) + "     0.0000     " + x + "     " + y + "     " + z + type + "\n")
        
        i = i + 1
        ii = ii + 1
        if ii > 30:
            ii = 1
            mid = mid + 1

#    tid = 4
#    type = "    #  SOD"
#    for n in range(len(ION)):
#            x = FORMAT % float(str(ION[n][0]))
#            y = FORMAT % float(str(ION[n][1]))
#            z = FORMAT % float(str(ION[n][2]))
#            dataFile.write("          " + str(i) + "        " + str(mid) + "        " + str(tid) + "     1.0000     " + x + "     " + y + "     " + z + type + "\n")
#        
#            i = i + 1
#            mid = mid + 1
    dataFile.write("\n")


def write_bonds(dataFile, NUM):
    dataFile.write("Bonds\n\n")
    total_NUM = len(NUM) * 29 / 30
    i = 1
    y = 1
    ii = 1
    tid = 1
    type = "    #  CCO-CCO"
    while i <= total_NUM:
#        if ii == 1:
#            tid = 1
#            type = "    #  SiT-SiM"
#        elif ii > 1 and ii <= 18:
#            tid = 2
#            type = "    #  SiM-SiM"
#        elif ii == 19:
#            tid = 1
#            type = "    #  SiT-SiM"
        type = "    #  CCO-CCO"
        dataFile.write("          " + str(i) + "        " + str(tid) + "        " + str(y) + "        " + str(y+1) + type + "\n")

        i = i + 1
        y = y + 1
        ii = ii + 1
        if ii > 29:
            ii = 1
            y = y + 1

    dataFile.write("\n")

def write_angles(dataFile, NUM):
    dataFile.write("Angles\n\n")
    total_NUM = len(NUM) * 28 / 30
    i = 1
    y = 1
    ii = 1
    tid = 1
    type = "    #  CCO-CCO-CCO"
    while i <= total_NUM:
#        if ii == 1:
#            tid = 1
#            type = "    #  SiT-SiM-SiM"
#        elif ii > 1 and ii <= 17:
#            tid = 2
#            type = "    #  SiM-SiM-SiM"
#        elif ii == 18:
#            tid = 1
#            type = "    #  SiT-SiM-SiM"
        type = "    # CCO-CCO-CCO" 
        dataFile.write("          " + str(i) + "        " + str(tid) + "        " + str(y) + "        " + str(y+1) + "        " + str(y+2) + type + "\n")

        i = i + 1
        y = y + 1
        ii = ii + 1
        if ii > 28:
            ii = 1
            y = y + 2

    dataFile.write("\n")

def write_dihedrals(dataFile, NUM):
    dataFile.write("Dihedrals\n\n")
    total_NUM = len(NUM) * 27 / 30
    i = 1
    y = 1
    ii = 1
    tid = 1
    type = "    #  CCO-CCO-CCO-CCO"
    while i <= total_NUM:
#        if ii == 1:
#            tid = 1
#            type = "    #  SiT-SiM-SiM-SiM"
#        elif ii > 1 and ii <= 16:
#            tid = 2
#            type = "    #  SiM-SiM-SiM-SiM"
#        elif ii == 17:
#            tid = 1
#            type = "    #  SiT-SiM-SiM-SiM"
        type = "    #  CCO-CCO-CCO-CCO"
        dataFile.write("          " + str(i) + "        " + str(tid) + "        " + str(y) + "        " + str(y+1) + "        " + str(y+2)+ "        " + str(y+3) + type + "\n")

        i = i + 1
        y = y + 1
        ii = ii + 1
        if ii > 27:
            ii = 1
            y = y + 3

    dataFile.write("\n")

def write_paras(dataFile):
#    dataFile.write("Pair Coeffs\n\n")
#    dataFile.write("          1       0.469       4.5850   #  SiT-SiT\n")
#    dataFile.write("          2       0.420       4.5060   #  SiM-SiM \n\n")
    dataFile.write("Bond Coeffs\n\n")
    dataFile.write("          1      2.43  60.0  60.0  60.0   #  CCO-CCO\n\n")
#    dataFile.write("          2      33.4600       3.3500   #  SiM-SiM \n\n")

    dataFile.write("Angle Coeffs\n\n")
    dataFile.write("          1      100.0  4.2  4.2  4.2     #  CCO-CCO-CCO\n\n")
#    dataFile.write("          2       6.3300     104.200    #  SiM-SiM-SiM \n")

    dataFile.write("BondBond Coeffs\n\n")
    dataFile.write("          1      0.0  2.43  2.43   #  CCO-CCO\n\n")

    dataFile.write("BondAngle Coeffs\n\n")
    dataFile.write("          1      0.0  0.0  2.43  2.43   #  CCO-CCO\n\n")

    dataFile.write("Dihedral Coeffs\n\n")
    dataFile.write("          1       -0.140  0.03  0.58  0    #  CCO-CCO-CCO-CCO\n\n")
#    dataFile.write("          2       -0.140  0.03  0.58  0    #  SiM-SiM-SiM-SiM \n")

if __name__ == "__main__":
    input = sys.argv[1]
    output = sys.argv[2]
    NUM = []
    inf = open(input, 'r')
    for line in inf:
        if line.strip().startswith("Atoms"):
            break
        elif line.find("xlo") != -1:
            PBC_x = line 
        elif line.find("ylo") != -1:
            PBC_y = line 
        elif line.find("zlo") != -1:
            PBC_z = line 
    for line in inf:
        strs = line.strip().split()
        if len(strs) > 1:
            if  strs[2] == "1" or strs[2] == "2":
                x = float(strs[4])
                y = float(strs[5])
                z = float(strs[6])
                NUM.append((x, y, z))
        elif line.strip().startswith("Bonds"):
            break
    inf.close()

    dataFile = open(output, 'w')
    write_numbers(dataFile, len(NUM), (len(NUM))*29/30, (len(NUM))*28/30, (len(NUM))*27/30, 0)
    write_types(dataFile, 1, 1, 1, 1, 0)
    write_pbc(dataFile, PBC_x, PBC_y, PBC_z)
    write_mass(dataFile)
    write_atoms(dataFile, NUM)
    write_bonds(dataFile, NUM)
    write_angles(dataFile, NUM)
    write_dihedrals(dataFile, NUM)
    write_paras(dataFile)
    dataFile.close()

