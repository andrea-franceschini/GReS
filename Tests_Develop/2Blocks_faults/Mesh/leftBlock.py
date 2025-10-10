import gmsh
import sys

from share.doc.gmsh.examples.api.plugin import numElements

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)

outfile = sys.argv[1]
NX = int(sys.argv[2])
NY = int(sys.argv[3])
NZ = int(sys.argv[4])
elem_type = sys.argv[5]

if elem_type == 'hexa':
    flagRecombine = True
elif elem_type == 'tetra':
    flagRecombine = False


if not NY % 2 == 0:
    raise ValueError("Elements on Y direction must be even for BCs imposition")






gmsh.model.add('domainTop')
# ---------------------- TOP LAYER-------------------------

# set number of layers for extrusion
xMin = 0
yMin = 0
zMin = 0

X = 2.5
Y = 10
Z = 15

P1 = gmsh.model.occ.addPoint(xMin,     yMin,     zMin)
P2 = gmsh.model.occ.addPoint(xMin+X,   yMin,     zMin)
P3 = gmsh.model.occ.addPoint(xMin+X,   yMin+Y,   zMin)
P4 = gmsh.model.occ.addPoint(xMin,     yMin+Y,   zMin)

L1 = gmsh.model.occ.addLine(P1, P2)
L2 = gmsh.model.occ.addLine(P2, P3)
L3 = gmsh.model.occ.addLine(P3, P4)
L4 = gmsh.model.occ.addLine(P4, P1)

loop1 = gmsh.model.occ.addCurveLoop([L1, L2, L3, L4])
surf1 = gmsh.model.occ.addPlaneSurface([loop1])

outExtrude = gmsh.model.occ.extrude([(2, surf1)], 0, 0, Z, numElements=[NZ], recombine=flagRecombine)

gmsh.model.occ.synchronize()

if flagRecombine:
    gmsh.model.mesh.setRecombine(2, surf1)


gmsh.model.mesh.set_transfinite_curve(L1,NX+1)
gmsh.model.mesh.set_transfinite_curve(L2,NY+1)
gmsh.model.mesh.set_transfinite_curve(L3,NX+1)
gmsh.model.mesh.set_transfinite_curve(L4,NY+1)

gmsh.model.mesh.set_transfinite_surface(surf1)

# outExtrude[0][1] stores the tag of the surface opposite to the original surface after extrusion
gmsh.model.addPhysicalGroup(2, [outExtrude[5][1]], 3, name="fixed_wall")
gmsh.model.addPhysicalGroup(2, [surf1], 2, name="bottom_fixed")
gmsh.model.addPhysicalGroup(2, [outExtrude[3][1]], 1, name="interface")
gmsh.model.addPhysicalGroup(3, [outExtrude[1][1]], 1, name="left_block")

# ---------------------- MESH SETTINGS -------------------------

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

#if '-nopopup' not in sys.argv:
#    gmsh.fltk.run()

gmsh.write(outfile+".vtk")

gmsh.finalize()