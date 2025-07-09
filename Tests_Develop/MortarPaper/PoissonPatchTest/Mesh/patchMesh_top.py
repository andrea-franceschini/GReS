import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)

gmsh.model.add('domain')
# ---------------------- TOP BLOCK (3x3) -------------------------
NX_top, NY_top, NZ_top = 3, 3, 3
xMin = 0
yMin = 0
zTop = 1
X = 1
Y = 1
Z = 1

P1t = gmsh.model.occ.addPoint(xMin,     yMin,     zTop)
P2t = gmsh.model.occ.addPoint(xMin+X,   yMin,     zTop)
P3t = gmsh.model.occ.addPoint(xMin+X,   yMin+Y,   zTop)
P4t = gmsh.model.occ.addPoint(xMin,     yMin+Y,   zTop)

L1t = gmsh.model.occ.addLine(P1t, P2t)
L2t = gmsh.model.occ.addLine(P2t, P3t)
L3t = gmsh.model.occ.addLine(P3t, P4t)
L4t = gmsh.model.occ.addLine(P4t, P1t)

loop2 = gmsh.model.occ.addCurveLoop([L1t, L2t, L3t, L4t])
surf2 = gmsh.model.occ.addPlaneSurface([loop2])

outTop = gmsh.model.occ.extrude([(2, surf2)], 0, 0, Z, numElements=[NZ_top], recombine=True)

gmsh.model.occ.synchronize()

gmsh.model.mesh.set_transfinite_curve(L1t, NX_top + 1)
gmsh.model.mesh.set_transfinite_curve(L2t, NY_top + 1)
gmsh.model.mesh.set_transfinite_curve(L3t, NX_top + 1)
gmsh.model.mesh.set_transfinite_curve(L4t, NY_top + 1)
gmsh.model.mesh.set_transfinite_surface(surf2)
gmsh.model.mesh.set_transfinite_volume(outTop[1][1])

gmsh.model.mesh.setRecombine(2, surf2)

gmsh.model.addPhysicalGroup(2, [loop2], 1, name="top_interface")
gmsh.model.addPhysicalGroup(2, [outTop[0][1]], 2, name="top_load")
gmsh.model.addPhysicalGroup(outTop[1][0], [outTop[1][1]], 1, name="top_block")

# ---------------------- MESH SETTINGS -------------------------

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

#if '-nopopup' not in sys.argv:
#    gmsh.fltk.run()
gmsh.write('top.vtk')

gmsh.finalize()
