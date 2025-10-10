import gmsh
import sys

from share.doc.gmsh.examples.api.plugin import numElements

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)

flagRecombine = True

gmsh.model.add('domainTop')

# set number of layers for extrusion
xMin = 0
yMin = 0
zMin = 0

X = 8
Y = 1
Z = 1

P1 = gmsh.model.occ.addPoint(0,     0,     0)
P2 = gmsh.model.occ.addPoint(X,    0,     0)
P3 = gmsh.model.occ.addPoint(X,   1,   0)
P4 = gmsh.model.occ.addPoint(0,     1,   0)

L11 = gmsh.model.occ.addLine(P1, P2)
L12 = gmsh.model.occ.addLine(P2, P3)
L13 = gmsh.model.occ.addLine(P3, P4)
L14 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L11, L12, L13, L14])
surf1 = gmsh.model.occ.addPlaneSurface([loop])

outExtrude1 = gmsh.model.occ.extrude([(2, surf1)], 0, 0, 1, numElements=[20], recombine=flagRecombine)

gmsh.model.occ.synchronize()

if flagRecombine:
    gmsh.model.mesh.setRecombine(2, surf1)


gmsh.model.mesh.set_transfinite_curve(L11,21)
gmsh.model.mesh.set_transfinite_curve(L12,13)
gmsh.model.mesh.set_transfinite_curve(L13,21)
gmsh.model.mesh.set_transfinite_curve(L14,13)

gmsh.model.mesh.set_transfinite_surface(surf1)

# outExtrude[0][1] stores the tag of the surface opposite to the original surface after extrusion
gmsh.model.addPhysicalGroup(2, [outExtrude1[3][1]], 1, name="interface_left")
gmsh.model.addPhysicalGroup(3, [outExtrude1[1][1]], 1, name="left_block")



P1 = gmsh.model.occ.addPoint(X,     0,     0)
P2 = gmsh.model.occ.addPoint(2*X,    0,     0)
P3 = gmsh.model.occ.addPoint(2*X,   1,   0)
P4 = gmsh.model.occ.addPoint(X,     1,   0)

L21 = gmsh.model.occ.addLine(P1, P2)
L22 = gmsh.model.occ.addLine(P2, P3)
L23 = gmsh.model.occ.addLine(P3, P4)
L24 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L21, L22, L23, L24])
surf2 = gmsh.model.occ.addPlaneSurface([loop])

outExtrude2 = gmsh.model.occ.extrude([(2, surf2)], 0, 0, 1, numElements=[12], recombine=flagRecombine)

gmsh.model.occ.synchronize()

if flagRecombine:
    gmsh.model.mesh.setRecombine(2, surf2)


gmsh.model.mesh.set_transfinite_curve(L21,21)
gmsh.model.mesh.set_transfinite_curve(L22,13)
gmsh.model.mesh.set_transfinite_curve(L23,21)
gmsh.model.mesh.set_transfinite_curve(L24,13)

gmsh.model.mesh.set_transfinite_surface(surf2)

# outExtrude[0][1] stores the tag of the surface opposite to the original surface after extrusion
gmsh.model.addPhysicalGroup(2, [outExtrude2[5][1]], 2, name="interface_right")
gmsh.model.addPhysicalGroup(2, [outExtrude1[5][1],outExtrude2[3][1]], 3, name="fixed")
gmsh.model.addPhysicalGroup(2, [outExtrude1[0][1], outExtrude2[0][1]], 4, name="top_load")
gmsh.model.addPhysicalGroup(3, [outExtrude2[1][1]], 2, name="right_block")

# ---------------------- MESH SETTINGS -------------------------

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.write("beam.vtk")

gmsh.finalize()