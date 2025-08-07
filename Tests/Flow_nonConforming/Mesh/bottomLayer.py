import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)

gmsh.model.add('domainTop')
# ---------------------- TOP LAYER-------------------------

# set number of layers for extrusion
xMin = 0
yMin = 0
zMin = 0

X = 20
Y = 10
Z = 5

mesh_density = 2.5

P1 = gmsh.model.occ.addPoint(xMin,     yMin,     zMin,   mesh_density)
P2 = gmsh.model.occ.addPoint(xMin+X,   yMin,     zMin,   mesh_density)
P3 = gmsh.model.occ.addPoint(xMin+X,   yMin+Y,   zMin,   mesh_density)
P4 = gmsh.model.occ.addPoint(xMin,     yMin+Y,   zMin,   mesh_density)

L1 = gmsh.model.occ.addLine(P1, P2)
L2 = gmsh.model.occ.addLine(P2, P3)
L3 = gmsh.model.occ.addLine(P3, P4)
L4 = gmsh.model.occ.addLine(P4, P1)

loop1 = gmsh.model.occ.addCurveLoop([L1, L2, L3, L4])
surf1 = gmsh.model.occ.addPlaneSurface([loop1])

outExtrude = gmsh.model.occ.extrude([(2, surf1)], 0, 0, Z, recombine=False)

gmsh.model.occ.synchronize()

# outExtrude[0][1] stores the tag of the surface opposite to the original surface after extrusion
gmsh.model.addPhysicalGroup(2, [outExtrude[0][1]], 1, name="bottom_interface")
gmsh.model.addPhysicalGroup(outExtrude[1][0], [outExtrude[1][1]], 1, name="bottom_layer")

# ---------------------- MESH SETTINGS -------------------------

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

#if '-nopopup' not in sys.argv:
#    gmsh.fltk.run()

gmsh.write('bottomLayer.vtk')

gmsh.finalize()