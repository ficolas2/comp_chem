#!/usr/bin/env python 

import numpy as np
from scipy.linalg import eigh

molecule = input("Titulo: ")
n = int(input("Numero de atomos: "))
ne = int(input("Numero de electrones: "))
bonds = input("Indices de los enlaces: ")
c = [int(s.strip()) for s in bonds.split (" ")]
nbonds = int(len(c)/2)
nloc = int(input("Numero de enlaces localizados: "))

for i in range(len(c)):
    c[i] -= 1

H = np.zeros((n, n))
for j in range(0, len(c), 2):
    H[c[j],c[j+1]] = 1.0
    H[c[j+1],c[j]] = 1.0

print()
print("Full hamiltonian for %s" %molecule)
print(H)
print()

energies, eigv = eigh(H)
sort_perm = (-energies).argsort()
e = energies[sort_perm]
orbitals = eigv[:,sort_perm]

noc = np.zeros(n) # Calculamos la ocupacion de cada OM
nocc = int(ne/2)
for i in range(nocc) :
    noc[i] = 2.0
if (ne%2) == 1:
    if (abs(e[nocc] - e[nocc-1]) < 1e-5):
        noc[nocc-1] = 1.5
        noc[nocc] = 1.5
    else :
        noc[nocc] = 1.0
    nocc += 1
if (abs(e[nocc] - e[nocc-1]) < 1e-5):
    noc[nocc] = 0.5 * noc[nocc-1]
    noc[nocc-1] = noc[nocc]
    nocc += 1

total_energy = np.dot(noc, e)
print(total_energy)

print("     OM     Energy    ocupation coefficients")
for j in range (n):
    print(" %6d " % (j+1), " %8.4f " % e[j], " %8.4f " % noc[j], end = " ")
    for i in range (n) :
        print("%8.5f " % orbitals[i, j], end = "")
    print()

# Energia de deslocalizacion
de = total_energy - 2* nloc
print("Total energy = %8.4f " % ne , " alpha + %8.4f beta " % total_energy )
print("Delocalization energy = %8.4f beta " % de )
print()

# Ordenes de enlace
bo = np.zeros(nbonds) # Vector para calcular los ordenes de enlace
print ( "\n Atoms           Bond order ")
for j in range(nbonds):
    k1 = c[2*j]
    k2 = c[2*j+1]
    for i in range (nocc):
        bo [j] += noc[i]* orbitals[k1, i]* orbitals[k2, i]
    print ( " %4d " % ( k1 +1) , "-%3d " % (k2 +1) , " %8.3f " % bo [j])

# Densidad pi sobre cada atomo
ch = np.zeros(n)
for i in range (nocc):
    for j in range (n):
         ch[j] += noc[i] * orbitals[j, i]**2

# Indice de valencia libre
free_valence = np.full(n, np.sqrt(3))
for j in range (nbonds):
    k1 = c[2*j]
    k2 = c[2*j+1]
    free_valence[k1] -= bo[j]
    free_valence[k2] -= bo[j]

print("\n Atom          Charge density   Net charge   Free valence")
for i in range (n):
    print(" %6d " % (i+1) , " %10.3f " % ch[i] , " %14.3f " % (1.0 - ch[i]) ," %10.3f " % free_valence[i])
