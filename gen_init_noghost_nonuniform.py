# generate init lattice for lammps
# no arguments needed

import numpy as np

L = 13 # number of cells (length)
# using cubic lattice
H = L
D = L

# number of lattice points (fcc)
N = L * H * D * 4
# number of vacancies
Nvacancy = 8
# number of atoms
Natom = N - Nvacancy

# lattice constant
a0 = 1.6 # np.sqrt(2./np.sqrt(3.)/rho)

# for lattice coordinates 
xh = L * a0   
yh = H * a0
zh = D * a0

# output file name
filename = '8void13cell.LJ.init'

# generate Vkl using g(V)   V[-1,-0.25]
G0 = 0.3
def gV(G0):
    if np.random.rand() < G0:
        # uniform component
        return np.random.uniform(0.25, 1.0)
    else:
        # delta component
        return 0.25

with open(filename, 'w') as fw:
    fw.write('LAMMPS data file via genConf.py, version 7 Jul 2017, timestep = 0\n\n')
    fw.write('%d atoms\n%d atom types\n\n' % (Natom, Natom))
    fw.write('0.0e+00 %.16g xlo xhi\n 0.0e+00 %.16g ylo yhi\n 0.0e+00 %.16g zlo zhi\n ' % (xh, yh, zh))

    # mass
    fw.write('\nMasses\n\n')
    for i in range(Natom):
        fw.write('%d 1\n' % (i + 1))


    # PairIJ Coeffs
    fw.write('\nPairIJ Coeffs # lj/cut\n\n')
    # uniform interaction energy
    #Vkl = np.random.uniform(low=0.25, high=1.0, size=int(Natom * (Natom + 1) / 2))
    # g(V)
    Vkl = [gV(G0) for _ in range(int(Natom * (Natom + 1) / 2))]

    k = 0
    for i in range(Natom):
        for j in range(i, Natom):
            fw.write('%d %d %f 1\n' % ((i + 1), (j + 1), Vkl[k]))
            k += 1


    # Atoms
    fw.write('\nAtoms # atomic\n\n')
    m = 0
    ID = 0
    # randomly choose vacancies
    vacancyIDs = np.random.choice(N, Nvacancy, replace=False) + 1
    
    
    x1 = 0.0
    y1 = 0.0
    z1 = 0.0
    for k in range(D):
        for j in range(H):
            for i in range(L):
                ID += 1
                if ID in vacancyIDs:
                    continue
                m += 1
                x = i * a0 + x1
                y = j * a0 + y1
                z = k * a0 + z1
                # lattice.append([x,y,z])
                fw.write('%d %d %.16g %.16g %.16g\n' % (m, m, x, y, z))

    x2 = 0.5 * a0
    y2 = 0.5 * a0
    z2 = 0.0
    for k in range(D):
        for j in range(H):
            for i in range(L):
                ID += 1
                if ID in vacancyIDs:
                    continue
                m += 1
                x = i * a0 + x2
                y = j * a0 + y2
                z = k * a0 + z2
                # lattice.append([x,y,z])
                fw.write('%d %d %.16g %.16g %.16g\n' % (m, m, x, y, z))

    x3 = 0.0
    y3 = 0.5 * a0
    z3 = 0.5 * a0
    for k in range(D):
        for j in range(H):
            for i in range(L):
                ID += 1
                if ID in vacancyIDs:
                    continue
                m += 1
                x = i * a0 + x3
                y = j * a0 + y3
                z = k * a0 + z3
                # lattice.append([x,y,z])
                fw.write('%d %d %.16g %.16g %.16g\n' % (m, m, x, y, z))

    x4 = 0.5 * a0
    y4 = 0.0
    z4 = 0.5 * a0
    for k in range(D):
        for j in range(H):
            for i in range(L):
                ID += 1
                if ID in vacancyIDs:
                    continue
                m += 1
                x = i * a0 + x4
                y = j * a0 + y4
                z = k * a0 + z4
                # lattice.append([x,y,z])
                fw.write('%d %d %.16g %.16g %.16g\n' % (m, m, x, y, z))
    print(vacancyIDs)