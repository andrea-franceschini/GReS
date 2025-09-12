import gmsh

import os

import numpy as np

import sys

gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 0)  # Disable terminal output
gmsh.option.setNumber("General.Verbosity", 0)  # Reduce verbosity to minimum

gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.1)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.1)


# read centers and radius array
R = np.loadtxt('rad.txt')
P = np.loadtxt('pts.txt')

b = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
vTag = np.zeros((1,R.shape[0]))

count = 0
for i in range(R.shape[0]):
    s = gmsh.model.occ.addSphere(P[i][0], P[i][1], P[i][2], R[i])
    gmsh.model.occ.intersect([(3, s)], [(3, b)],removeObject=True,removeTool=False)

gmsh.model.occ.remove([(3, b)], recursive=True)

gmsh.model.occ.synchronize()


# define physical groups
c = 1

for i in range(2,R.shape[0]+2):
    gmsh.model.addPhysicalGroup(3, [i], c)
    c += 1
    print(i)

gmsh.model.mesh.generate(3)

# ... and save it to disk
outfile = 'pores.vtk'

gmsh.write(outfile)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()



