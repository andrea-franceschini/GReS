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

X = 8
Y = 1
Z = 1

NX = 32;
NY = 4;
NZ = 4;

mesh_density = 1

P1 = gmsh.model.occ.addPoint(xMin,   yMin,     zMin,   mesh_density)
P2 = gmsh.model.occ.addPoint(xMin,   yMin+Y,   zMin,   mesh_density)
P3 = gmsh.model.occ.addPoint(xMin,   yMin+Y,   zMin+Z,   mesh_density)
P4 = gmsh.model.occ.addPoint(xMin,   yMin,     zMin+Z,   mesh_density)

L1 = gmsh.model.occ.addLine(P1, P2)
L2 = gmsh.model.occ.addLine(P2, P3)
L3 = gmsh.model.occ.addLine(P3, P4)
L4 = gmsh.model.occ.addLine(P4, P1)

loop1 = gmsh.model.occ.addCurveLoop([L1, L2, L3, L4])
surf1 = gmsh.model.occ.addPlaneSurface([loop1])

outExtrude = gmsh.model.occ.extrude([(2, surf1)], X, 0, 0, numElements=[NX], recombine=True)



gmsh.model.occ.synchronize()

gmsh.model.mesh.set_transfinite_curve(L1,NY+1)
gmsh.model.mesh.set_transfinite_curve(L2,NZ+1)
gmsh.model.mesh.set_transfinite_curve(L3,NY+1)
gmsh.model.mesh.set_transfinite_curve(L4,NZ+1)

gmsh.model.mesh.set_transfinite_surface(1)

gmsh.model.mesh.setRecombine(2, surf1)

gmsh.model.addPhysicalGroup(2, [loop1], 1, name="fix_left")
gmsh.model.addPhysicalGroup(2, [outExtrude[0][1]], 2, name="interface_left")
gmsh.model.addPhysicalGroup(3, [outExtrude[1][1]], 1, name="left_beam")

# ---------------------- MESH SETTINGS -------------------------

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.write('leftBeam.vtk')

gmsh.finalize()