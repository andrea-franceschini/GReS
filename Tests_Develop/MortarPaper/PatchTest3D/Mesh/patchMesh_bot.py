import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)

gmsh.model.add('domain')


# ---------------------- BOTTOM BLOCK (4x4) -------------------------
NX_bot, NY_bot, NZ_bot = 4, 4, 4
X, Y, Z = 1.0, 1.0, 1.0
xMin, yMin, zMin = 0.0, 0.0, 0.0

P1b = gmsh.model.occ.addPoint(xMin,     yMin,     zMin)
P2b = gmsh.model.occ.addPoint(xMin+X,   yMin,     zMin)
P3b = gmsh.model.occ.addPoint(xMin+X,   yMin+Y,   zMin)
P4b = gmsh.model.occ.addPoint(xMin,     yMin+Y,   zMin)

L1b = gmsh.model.occ.addLine(P1b, P2b)
L2b = gmsh.model.occ.addLine(P2b, P3b)
L3b = gmsh.model.occ.addLine(P3b, P4b)
L4b = gmsh.model.occ.addLine(P4b, P1b)

loop1 = gmsh.model.occ.addCurveLoop([L1b, L2b, L3b, L4b])
surf1 = gmsh.model.occ.addPlaneSurface([loop1])

outBot = gmsh.model.occ.extrude([(2, surf1)], 0, 0, Z, numElements=[NZ_bot], recombine=True)

gmsh.model.occ.synchronize()

gmsh.model.mesh.set_transfinite_curve(L1b, NX_bot + 1)
gmsh.model.mesh.set_transfinite_curve(L2b, NY_bot + 1)
gmsh.model.mesh.set_transfinite_curve(L3b, NX_bot + 1)
gmsh.model.mesh.set_transfinite_curve(L4b, NY_bot + 1)
gmsh.model.mesh.set_transfinite_surface(surf1)
gmsh.model.mesh.set_transfinite_volume(outBot[1][1])

gmsh.model.mesh.setRecombine(2, surf1)

gmsh.model.addPhysicalGroup(2, [surf1], 1, name="bottom")
gmsh.model.addPhysicalGroup(2, [outBot[0][1]], 2, name="bot_interface")
gmsh.model.addPhysicalGroup(outBot[1][0], [outBot[1][1]], 1, name="bot_block")


gmsh.model.mesh.generate(3)

#if '-nopopup' not in sys.argv:
#    gmsh.fltk.run()
gmsh.write('bottom.vtk')

gmsh.finalize()

