#!/usr/bin/env pypy
# coding=utf-8

import sys
import math
import argparse

parser = argparse.ArgumentParser(description='Analyse the cluster of surfactants in saline solution. '
                                             'Modify lines start with -_- . '
                                             'Be aware that polar bead should be at first or last position.')
parser.add_argument('trj_file', help='input trajectory')
parser.add_argument('ion_types', nargs='*', type=str, help='types of ion beads')
parser.add_argument('-b', dest='BEGIN', type=int, default=0, help='begin from this frame')
parser.add_argument('-e', dest='END', type=int, default=0, help='end at this frame')
parser.add_argument('-s', dest='SKIP', type=int, default=1, help='only read every nr-th frame')
parser.add_argument('-g', dest='GYRATION', default=False, action='store_true',
                    help='calculate the gyration tensor of micelles')
parser.add_argument('-t', dest='WRITE_TRJ', default=False, action='store_true',
                    help='write the trajectory of calculated atoms')
parser.add_argument('-r', dest='RDF', default=False, action='store_true',
                    help='calculate RDF between beads to COM of the largest micelle')
parser.add_argument('--br', dest='RDF_BEGIN', type=int, default=0, help='calculate RDF from this frame')
parser.add_argument('--mr', dest='MAX_R', type=int, default=40, help='max distance when calculate RDF')
parser.add_argument('--dr', dest='D_R', type=float, default=0.1, help='distance interval when calculate RDF')
parser.add_argument('--cut', dest='CUT_SPH', default=False, action='store_true', help='cut sphere when calculate RDF')
parser.add_argument('--head', dest='HEAD', default=False, action='store_true', help='polar group is at first')
parser.add_argument('--gmx', dest='GMX', default=False, action='store_true', help='GROMACS gro file')
args = parser.parse_args()

BEGIN = args.BEGIN
END = args.END
SKIP = args.SKIP
CALC_GYRATION = args.GYRATION
WRITE_TRJ = args.WRITE_TRJ
CALC_RDF = args.RDF
RDF_BEGIN = args.RDF_BEGIN if args.RDF_BEGIN != 0 else args.BEGIN
MAX_R = args.MAX_R
D_R = args.D_R
CUT_SPH = args.CUT_SPH
HEAD = args.HEAD
GMX = 1 if args.GMX else 0
SCALE = 0 if GMX else 1

# -_- modify this according to your need
nsSteps = 1e5
DIS_TAIL = 0.0
R_CG = 1.0
box = [0, 0, 0]

# -_- modify this according to your need
SUR_STRUC = {
    'PDMS': ['1','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2',
    '2','2','2','2','2','2','2','2','2','2','2','2','2','2','2',
    '2','2','2','2','2','2','2','2','2','2','2','2','2','2','2',
    '2','2','2','2','2','2','2','1'],
}

SUR_NUMBER = {i: 0 for i in SUR_STRUC.keys()}
SUR_TYPES = []
for sur in SUR_STRUC.values():
    for i in sur:
        if not i in SUR_TYPES:
            SUR_TYPES.append(i)
ION_TYPES = args.ion_types
ATOM_TYPES = SUR_TYPES + ION_TYPES
print ATOM_TYPES


# make the rdf list
r_segment = int(MAX_R / D_R)
r_list = range(r_segment + 1)
r = [round((r_list[i] + r_list[i + 1]) / 2. * D_R, 3) for i in range(r_segment)]
rdf_count = [[0 for i in r] for j in ATOM_TYPES]
rdf = [[0 for i in r] for j in ATOM_TYPES]


def find_composition(all, components):
    composition = []
    sort_components = sorted(components.values(), reverse=1, key=len)
    checked = 0
    while checked < len(all):
        found = 0
        for i in sort_components:
            if checked + len(i) > len(all):
                continue
            if all[checked:checked + len(i)] == i:
                checked += len(i)
                component = [k for k, v in components.items() if v == i][0]  # dict reverse lookup
                SUR_NUMBER[component] += 1
                composition.append(component)
                found = 1
                break
        if found == 0:
            print 'surfactants not found'
            sys.exit()
    return composition


def move_mic_x(xMIC, x):
    new_x = []
    xMIC.sort()
    dx = [xMIC[i + 1] - xMIC[i] for i in range(len(xMIC) - 1)]
    max_dx = max(dx)
    dx_left = min(xMIC)
    dx_right = 1 - max(xMIC)
    if max_dx > dx_left and max_dx > dx_right:
        hehe = dx.index(max_dx)
        move_left = xMIC[hehe] + max_dx / 2.
    else:
        move_left = (dx_left - dx_right) / 2.
    for i in x:
        tmp_x = i - move_left
        if tmp_x < 0:
            tmp_x += 1
        if tmp_x > 1:
            tmp_x -= 1
        new_x.append(tmp_x)
    return new_x


def move_mic(xyzMIC, xyz):
    xMIC = [i[0] for i in xyzMIC]
    x = [i[0] for i in xyz]
    new_x = move_mic_x(xMIC, x)
    yMIC = [i[1] for i in xyzMIC]
    y = [i[1] for i in xyz]
    new_y = move_mic_x(yMIC, y)
    zMIC = [i[2] for i in xyzMIC]
    z = [i[2] for i in xyz]
    new_z = move_mic_x(zMIC, z)
    return [[new_x[i], new_y[i], new_z[i]] for i in range(len(new_x))]


def calc_center(xyz):
    center = [0, 0, 0]
    for i in xyz:
        for j in range(3):
            center[j] += i[j]
    n = len(xyz)
    return [center[i] / n for i in range(3)]


def calc_dis(xyz1, xyz2):
    dx = abs(xyz2[0] - xyz1[0]) * box[0]
    dy = abs(xyz2[1] - xyz1[1]) * box[1]
    dz = abs(xyz2[2] - xyz1[2]) * box[2]
    if dx > box[0] / 2.:
        dx = box[0] - dx
    if dy > box[1] / 2.:
        dy = box[1] - dy
    if dz > box[2] / 2.:
        dz = box[2] - dz
    if dx < MAX_R and dy < MAX_R and dz < MAX_R:
        return math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    else:
        return 0


def dis_SUR_tail(xyzSUR1, xyzSUR2):
    if not HEAD:
        for i in xyzSUR1[:]:
            for j in xyzSUR2[:]:
                dis = calc_dis(i, j)
                if 0 < dis < DIS_TAIL:
                    return 1
        return 0
    else:
        for i in xyzSUR1[1:]:
            for j in xyzSUR2[1:]:
                dis = calc_dis(i, j)
                if 0 < dis < DIS_TAIL:
                    return 1
        return 0


def find_cluster(xyzSUR):
    Cluster = []
    Monomer = []
    flag = [0 for i in range(nSUR)]
    Ncluster = 0
    for i in range(nSUR):
        if flag[i] != 0:
            continue
        path = [i]
        expand = []
        Ncluster += 1
        flag[i] = Ncluster
        for j in range(i + 1, nSUR):
            if flag[j] != 0:
                continue
            if dis_SUR_tail(xyzSUR[i], xyzSUR[j]) > 0:
                path.append(j)
                expand.append(j)
                flag[j] = Ncluster
        while expand != []:
            m = expand[0]
            for n in range(i + 1, nSUR):
                if flag[n] != 0:
                    continue
                if dis_SUR_tail(xyzSUR[m], xyzSUR[n]) > 0:
                    path.append(n)
                    expand.append(n)
                    flag[n] = Ncluster
            expand.remove(m)
        if len(path) >= 1:
            Cluster.append(path)
        else:
            Monomer.append(path)
        Cluster.sort(key=len, reverse=True)
    return Cluster, len(Monomer)


def cut_sph(R, r, d):
    if R >= d + r:
        Vpercent = 1
    elif R <= r - d:
        Vpercent = (R / r) ** 3
    elif R <= d - r:
        Vpercent = 0
    else:
        x = (R * R - r * r + d * d) / (2 * d)
        H = R - x
        h = r - d + x
        V1percent = H * H * (3 * R - H) / 4. / r ** 3
        V2percent = h * h * (3 * r - h) / 4. / r ** 3
        Vpercent = V1percent + V2percent
    return Vpercent


def calc_rdf(com, xyz, rdf_count):
    for i in xyz:
        dis = calc_dis(com, i)
        if CUT_SPH:
            if dis > 0:
                rdf_count[0] += cut_sph(r[0], R_CG, dis)
                for j in range(r_segment - 1):
                    rdf_count[j + 1] += cut_sph(r[j + 1], R_CG, dis) - cut_sph(r[j], R_CG, dis)
        else:
            if dis < MAX_R and dis > 0:
                rdf_count[int(dis / D_R)] += 1


def cubic_root(a, b, c, d):
    alpha = -b ** 3 / 27 / a ** 3 - d / 2 / a + b * c / 6 / a ** 2
    beta = c / 3 / a - b ** 2 / 9 / a ** 2
    theta = alpha / pow(-beta, 1.5)
    x1 = -b / 3 / a + 2 * math.sqrt(-beta) * math.cos(math.acos(theta) / 3)
    x2 = -b / 3 / a + 2 * math.sqrt(-beta) * math.cos(math.acos(theta) / 3 + 2 * math.pi / 3)
    x3 = -b / 3 / a + 2 * math.sqrt(-beta) * math.cos(math.acos(theta) / 3 - 2 * math.pi / 3)
    root = [x1, x2, x3]
    root.sort(reverse=True)
    return root


def calc_eig3(matrix):
    a00 = matrix[0][0]
    a01 = matrix[0][1]
    a02 = matrix[0][2]
    a10 = matrix[1][0]
    a11 = matrix[1][1]
    a12 = matrix[1][2]
    a20 = matrix[2][0]
    a21 = matrix[2][1]
    a22 = matrix[2][2]
    a = -1
    b = a00 + a11 + a22
    c = a01 * a10 - a00 * a11 - a00 * a22 + a02 * a20 - a11 * a22 + a12 * a21
    d = a00 * a11 * a22 - a00 * a12 * a21 - a01 * a10 * a22 + a01 * a12 * a20 + a02 * a10 * a21 - a02 * a11 * a20
    return cubic_root(a, b, c, d)


def calc_rg(xyz):
    center = calc_center(xyz)
    G = [[0 for i in range(3)] for j in range(3)]
    for m in range(3):
        for n in range(3):
            for i in xyz:
                G[m][n] += (i[m] - center[m]) * (i[n] - center[n])
            G[m][n] = G[m][n] / len(xyz)
    l1, l2, l3 = calc_eig3(G)
    Rg = (l1 + l2 + l3) ** 0.5  # radius of gyration
    b = l1 - (l2 + l3) / 2  # asphercity
    c = l2 - l3  # acylindricity
    k2 = (b ** 2 + 0.75 * c ** 2) / Rg ** 4  # anisotropy
    return [round(l1, 1), round(l2, 1), round(l3, 1), round(Rg, 1), round(k2, 4)]


def write_trj(dumpfile, timestep, box, xyz, typeBead):
    nAtoms = 0
    for i in xyz:
        nAtoms += len(i)
    if not GMX:
        dumpfile.write(('ITEM: TIMESTEP\n'
                        '%i\n'
                        'ITEM: NUMBER OF ATOMS\n'
                        '%i\n'
                        'ITEM: BOX BOUNDS pp pp pp\n'
                        '0 %f\n'
                        '0 %f\n'
                        '0 %f\n'
                        'ITEM: ATOMS id type xs ys zs\n'
                        ) % (timestep, nAtoms, box[0], box[1], box[2]))
    else:
        dumpfile.write('Generated by python : t= %f\n%i\n' % (timestep * 1000 / nsSteps, nAtoms))

    nAtom = 0
    for i in range(1 + len(ION_TYPES)):
        for j in xyz[i]:
            nAtom += 1
            if i == 0:
                if not GMX:
                    dumpfile.write('%i %s %f %f %f\n' % (nAtom, typeBead[nAtom - 1], j[0], j[1], j[2]))
                else:
                    dumpfile.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' % (
                        1, 'SUR', typeBead[nAtom - 1], nAtom, j[0] * box[0] / 10, j[1] * box[1] / 10,
                        j[2] * box[2] / 10))  # convert to nanometer
            else:
                if not GMX:
                    dumpfile.write('%i %s %f %f %f\n' % (nAtom, ION_TYPES[i - 1], j[0], j[1], j[2]))
                else:
                    dumpfile.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' % (
                        2, 'ION', ION_TYPES[i - 1], nAtom, j[0] * box[0] / 10, j[1] * box[1] / 10,
                        j[2] * box[2] / 10))  # convert to nanometer
    if GMX:
        dumpfile.write('%f %f %f\n' % (box[0] / 10, box[1] / 10, box[2] / 10))  # convert to nanometer


class AnalyseTrajectory():
    def __init__(self):
        self.inf = open(args.trj_file)
        self.cluster_file = open('CLUSTER.txt', 'w')
        # self.cluster_file.write(
            # '#Frame\tTimestep\tTime(ns)\tNcluster\tNmonomer\tCluster_size\tGyration[l1^2,l2^2,l3^2,Rg,k^2]\n')
        if WRITE_TRJ or CALC_RDF:
            if not GMX:
                self.dump_big_file = open('BIG.ltrj', 'w')
                self.dump_all_file = open('ALL.ltrj', 'w')
            else:
                self.dump_big_file = open('BIG.gro', 'w')
                self.dump_all_file = open('ALL.gro', 'w')

        self.frame = 0
        self.total_frame = 0
        self.start = 0
        self.nAtoms = 0


    def read_trj(self):
        if not GMX:
            for line in self.inf:
                if line.find('TIMESTEP') != -1:
                    n = 1
                    self.xyz_sur_bead = []
                    self.type_sur_bead = []
                    self.type_all_bead = []
                    self.xyz_ion = [[] for i in ION_TYPES]
                    self.frame += 1
                    if self.frame >= BEGIN:
                        self.start = 1
                    if END != 0 and self.frame > END:
                        break
                else:
                    if self.start == 0 or (self.frame - 1) % SKIP != 0:
                        continue
                    n += 1
                    word = line.strip().split()
                    if n == 2:
                        self.timestep = int(word[0])
                    elif n == 4:
                        self.nAtoms = int(word[0])
                    elif n == 6:
                        box[0] = float(word[1]) - float(word[0])
                        xlow = float(word[0])
                    elif n == 7:
                        box[1] = float(word[1]) - float(word[0])
                        ylow = float(word[0])
                    elif n == 8:
                        box[2] = float(word[1]) - float(word[0])
                        zlow = float(word[0])
                    elif n == 9 and word[4] != 'xs':
                        global SCALE
                        SCALE = 0
                    elif n > 9 and n <= 9 + self.nAtoms:
                        type = word[1]
                        if type in SUR_TYPES:
                            if SCALE:
                                xyz_tmp = [float(word[2]), float(word[3]), float(word[4])]
                            else:
                                xyz_tmp = [(float(word[2]) - xlow) / box[0], (float(word[3]) - ylow) / box[1],
                                           (float(word[4]) - zlow) / box[2]]
                            self.xyz_sur_bead.append(xyz_tmp)
                            self.type_sur_bead.append(type)
                            self.type_all_bead.append(type)
                        elif type in ION_TYPES:
                            self.type_all_bead.append(type)
                            i = ION_TYPES.index(type)
                            if SCALE:
                                xyz_tmp = [float(word[2]), float(word[3]), float(word[4])]
                            else:
                                xyz_tmp = [(float(word[2]) - xlow) / box[0], (float(word[3]) - ylow) / box[1],
                                           (float(word[4]) - zlow) / box[2]]
                            self.xyz_ion[i].append(xyz_tmp)

                    if n == 9 + self.nAtoms:
                        self.ana_trj()
        else:
            for line in self.inf:
                if line.find('Generated') != -1:
                    word = line.strip().split()
                    time = float(word[-1])
                    self.timestep = int(time / 1000 * nsSteps)
                    n = 1
                    self.xyz_sur_bead = []
                    self.type_sur_bead = []
                    self.type_all_bead = []
                    self.xyz_ion = [[] for i in ION_TYPES]
                    self.frame += 1
                    if self.frame >= BEGIN:
                        self.start = 1
                    if END != 0 and self.frame > END:
                        break
                else:
                    if self.start == 0 or (self.frame - 1) % SKIP != 0:
                        continue
                    n += 1
                    if n == 2:
                        self.nAtoms = int(line.strip())
                    elif n > 2 and n <= 2 + self.nAtoms:
                        type = line[10:15].strip()
                        if type in SUR_TYPES:
                            xyz_tmp = [float(line[20:28]), float(line[28:36]), float(line[36:44])]
                            self.xyz_sur_bead.append(xyz_tmp)
                            self.type_sur_bead.append(type)
                            self.type_all_bead.append(type)
                        elif type in ION_TYPES:
                            self.type_all_bead.append(type)
                            i = ION_TYPES.index(type)
                            xyz_tmp = [float(line[20:28]), float(line[28:36]), float(line[36:44])]
                            self.xyz_ion[i].append(xyz_tmp)

                    if n == 3 + self.nAtoms:
                        word = line.strip().split()
                        box[0] = float(word[0]) * 10  # convert the length to Angstrom
                        box[1] = float(word[1]) * 10  # convert the length to Angstrom
                        box[2] = float(word[2]) * 10  # convert the length to Angstrom
                        for i in range(len(self.xyz_sur_bead)):
                            for j in range(3):
                                self.xyz_sur_bead[i][j] /= (box[j] / 10)  # convert the length to Angstrom
                        for i in range(len(self.xyz_ion)):
                            for j in range(len(self.xyz_ion[i])):
                                for k in range(3):
                                    self.xyz_ion[i][j][k] /= (box[k] / 10)  # convert the length to Angstrom
                        self.ana_trj()
        if CALC_RDF:
            self.write_rdf()


    def ana_trj(self):
        self.total_frame += 1
        if self.total_frame == 1:
            self.composition = find_composition(self.type_sur_bead, SUR_STRUC)
            print SUR_NUMBER
        sys.stdout.write('\r\t%i' % self.total_frame)
        sys.stdout.flush()

        global nSUR
        nSUR = sum(SUR_NUMBER.values())
        xyzSUR = [[] for i in range(nSUR)]
        checked = 0
        for i in range(nSUR):
            xyzSUR[i] = self.xyz_sur_bead[checked:checked + len(SUR_STRUC[self.composition[i]])]
            checked += len(SUR_STRUC[self.composition[i]])

        cluster, Nmonomer = find_cluster(xyzSUR)
        Ncluster = len(cluster)
        if Ncluster == 0:
            cluster_size = []
            cluster_big_sur_id = [1]  # fake a cluster, for the sake of bug
        else:
            cluster_size = [len(i) for i in cluster]
            cluster_big_sur_id = cluster[0]
        # self.cluster_file.write('%i\t%i\t%.1f\t%i\t%i\t%s' % (self.frame, self.timestep, 1. * self.timestep / nsSteps,
                                                              # Ncluster, Nmonomer, cluster_size))

        if CALC_GYRATION:
            for m in cluster:
                xyz_tmp = []
                for n in m:
                    for p in xyzSUR[n]:
                        xyz_tmp.append(p)
                xyzTMP_move = move_mic(xyz_tmp, xyz_tmp)
                xyz_scale = [[i[0] * box[0], i[1] * box[1], i[2] * box[2]] for i in xyzTMP_move]
                gyration = calc_rg(xyz_scale)
                self.cluster_file.write('\t%s' % gyration[-2])

        self.cluster_file.write('\n')

        if WRITE_TRJ or (CALC_RDF and self.frame >= RDF_BEGIN):
            xyz_big_sur = []
            type_big_sur = []
            for j in cluster_big_sur_id:
                xyz_big_sur.append(xyzSUR[j])
                type_big_sur.append(self.composition[j])

            xyz_big_bead = []
            type_big_bead = []
            xyz_big_ion = [[] for j in ION_TYPES] + [[]]
            xyz_all = [[] for j in ION_TYPES] + [[]]

            for sur in xyz_big_sur:
                for bead in sur:
                    xyz_big_bead.append(bead)
            for sur in type_big_sur:
                type_big_bead += SUR_STRUC[sur]

            xyz_big_ion[0] = xyz_big_bead
            xyz_all[0] = self.xyz_sur_bead
            for j in range(len(ION_TYPES)):
                xyz_big_ion[j + 1] = self.xyz_ion[j]
                xyz_all[j + 1] = self.xyz_ion[j]

            xyz_big_ion_move = [move_mic(xyz_big_bead, xyz_big_ion[j]) for j in range(1 + len(ION_TYPES))]
            write_trj(self.dump_big_file, self.timestep, box, xyz_big_ion_move, type_big_bead)
            xyz_all_move = [move_mic(xyz_big_bead, xyz_all[j]) for j in range(1 + len(ION_TYPES))]
            write_trj(self.dump_all_file, self.timestep, box, xyz_all_move, self.type_all_bead)

        if CALC_RDF and self.frame >= RDF_BEGIN:
            COM = calc_center(xyz_big_ion_move[0])
            xyz_sur_bead_separated = [[] for j in SUR_TYPES]
            for j in range(len(type_big_bead)):
                jj = SUR_TYPES.index(type_big_bead[j])
                xyz_sur_bead_separated[jj].append(xyz_big_ion_move[0][j])
            for j in range(len(SUR_TYPES)):
                calc_rdf(COM, xyz_sur_bead_separated[j], rdf_count[j])
            for j in range(len(ION_TYPES)):
                calc_rdf(COM, xyz_big_ion_move[j + 1], rdf_count[j + len(SUR_TYPES)])

    def write_rdf(self):
        for i in range(len(ATOM_TYPES)):
            for j in range(r_segment):
                rdf[i][j] = rdf_count[i][j] / 4 / math.pi / r[j] ** 2 / D_R / self.total_frame
        rdf_file = open('RDF-COM.txt', 'w')
        rdf_file.write('r')
        for i in SUR_TYPES:
            rdf_file.write('\tBEAD%s' % i)
        for i in range(len(ION_TYPES)):
            rdf_file.write('\tION%s' % ION_TYPES[i])
        rdf_file.write('\n')
        for j in range(r_segment):
            rdf_file.write('%f' % r[j])
            for i in range(len(rdf)):
                rdf_file.write('\t%f' % rdf[i][j])
            rdf_file.write('\n')
        rdf_file.close()


if __name__ == '__main__':
    process = AnalyseTrajectory()
    process.read_trj()
    print ''

