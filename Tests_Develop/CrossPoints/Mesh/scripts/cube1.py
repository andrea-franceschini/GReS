import gmsh

import os

import sys

gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 0)  # Disable terminal output
gmsh.option.setNumber("General.Verbosity", 0)  # Reduce verbosity to minimum

outfile = sys.argv[1]

N1 = 3
N2 = 4
N3 = 5
N4 = 6


gmsh.model.add('domain')


recombineFlag = True


# cube 1 (bottom left)

X = 1
Y = 1
Z = 1

xMin = 0
yMin = 0
zMin = 0

P1 = gmsh.model.occ.addPoint(xMin,yMin,zMin)
P2 = gmsh.model.occ.addPoint(xMin+X,yMin,zMin)
P3 = gmsh.model.occ.addPoint(xMin+X,yMin+Y,zMin)
P4 = gmsh.model.occ.addPoint(xMin,yMin+Y,zMin)

L1 = gmsh.model.occ.addLine(P1,P2)
L2 = gmsh.model.occ.addLine(P2,P3)
L3 = gmsh.model.occ.addLine(P3,P4)
L4 = gmsh.model.occ.addLine(P4,P1)

loop_1 = gmsh.model.occ.addCurveLoop([L1,L2,L3,L4])

surf_1 = gmsh.model.occ.addPlaneSurface([loop_1])

# dim tag is [(geo_size,tag_entity)]  numElements is a list!
outExt = gmsh.model.occ.extrude([(2, surf_1)], 0, 0, 1, numElements=[N1], recombine=recombineFlag)

gmsh.model.occ.synchronize()

gmsh.model.mesh.setRecombine(2, surf_1)

gmsh.model.mesh.set_transfinite_curve(L1,N1+1)
gmsh.model.mesh.set_transfinite_curve(L2,N1+1)
gmsh.model.mesh.set_transfinite_curve(L3,N1+1)
gmsh.model.mesh.set_transfinite_curve(L4,N1+1)

gmsh.model.mesh.set_transfinite_surface(loop_1)

gmsh.model.addPhysicalGroup(2, [outExt[2 + L2 - L1][1]], 1, name="interf_x")
gmsh.model.addPhysicalGroup(2, [outExt[2 + L3 - L1][1]], 2, name="interf_y")

surf_bc = [outExt[2+i - L1][1] for i in [L1,L4]]
surf_tags = [surf_1] + [outExt[0][1]] + surf_bc
gmsh.model.addPhysicalGroup(2, surf_tags, 9, name = "bound_1")


gmsh.model.addPhysicalGroup(outExt[1][0], [outExt[1][1]],1, name="cube1")

# cube 2 - bottom right

xMin = 1
yMin = 0
zMin = 0

P1 = gmsh.model.occ.addPoint(xMin,yMin,zMin)
P2 = gmsh.model.occ.addPoint(xMin+X,yMin,zMin)
P3 = gmsh.model.occ.addPoint(xMin+X,yMin+Y,zMin)
P4 = gmsh.model.occ.addPoint(xMin,yMin+Y,zMin)

L1 = gmsh.model.occ.addLine(P1,P2)
L2 = gmsh.model.occ.addLine(P2,P3)
L3 = gmsh.model.occ.addLine(P3,P4)
L4 = gmsh.model.occ.addLine(P4,P1)

loop_2 = gmsh.model.occ.addCurveLoop([L1,L2,L3,L4])

surf_2 = gmsh.model.occ.addPlaneSurface([loop_2])

# dim tag is [(geo_size,tag_entity)]  numElements is a list!
outExt = gmsh.model.occ.extrude([(2, surf_2)], 0, 0, 1, numElements=[N2], recombine=recombineFlag)

gmsh.model.occ.synchronize()

gmsh.model.mesh.setRecombine(2, surf_2)

gmsh.model.mesh.set_transfinite_curve(L1,N2+1)
gmsh.model.mesh.set_transfinite_curve(L2,N2+1)
gmsh.model.mesh.set_transfinite_curve(L3,N2+1)
gmsh.model.mesh.set_transfinite_curve(L4,N2+1)

gmsh.model.mesh.set_transfinite_surface(loop_2)

gmsh.model.addPhysicalGroup(2, [outExt[2 + L4 - L1][1]], 3, name="interf_x_2")
gmsh.model.addPhysicalGroup(2, [outExt[2 + L3 - L1][1]], 4, name="interf_y_2")

surf_bc = [outExt[2+i - L1][1] for i in [L1,L2]]
surf_tags = [surf_2] + [outExt[0][1]] + surf_bc
gmsh.model.addPhysicalGroup(2, surf_tags, 10, name = "bound_2")


gmsh.model.addPhysicalGroup(outExt[1][0], [outExt[1][1]],2, name="cube2")

# cube 3 - top right

xMin = 1
yMin = 1
zMin = 0

P1 = gmsh.model.occ.addPoint(xMin,yMin,zMin)
P2 = gmsh.model.occ.addPoint(xMin+X,yMin,zMin)
P3 = gmsh.model.occ.addPoint(xMin+X,yMin+Y,zMin)
P4 = gmsh.model.occ.addPoint(xMin,yMin+Y,zMin)

L1 = gmsh.model.occ.addLine(P1,P2)
L2 = gmsh.model.occ.addLine(P2,P3)
L3 = gmsh.model.occ.addLine(P3,P4)
L4 = gmsh.model.occ.addLine(P4,P1)

loop_3 = gmsh.model.occ.addCurveLoop([L1,L2,L3,L4])

surf_3 = gmsh.model.occ.addPlaneSurface([loop_3])

# dim tag is [(geo_size,tag_entity)]  numElements is a list!
outExt = gmsh.model.occ.extrude([(2, surf_3)], 0, 0, 1, numElements=[N3], recombine=recombineFlag)

gmsh.model.occ.synchronize()

gmsh.model.mesh.setRecombine(2, surf_3)

gmsh.model.mesh.set_transfinite_curve(L1,N3+1)
gmsh.model.mesh.set_transfinite_curve(L2,N3+1)
gmsh.model.mesh.set_transfinite_curve(L3,N3+1)
gmsh.model.mesh.set_transfinite_curve(L4,N3+1)

gmsh.model.mesh.set_transfinite_surface(loop_3)

gmsh.model.addPhysicalGroup(2, [outExt[2 + L4 - L1][1]], 5, name="interf_x_3")
gmsh.model.addPhysicalGroup(2, [outExt[2 + L1 - L1][1]], 6, name="interf_y_3")

surf_bc = [outExt[2+i - L1][1] for i in [L2,L3]]
surf_tags = [surf_3] + [outExt[0][1]] + surf_bc
gmsh.model.addPhysicalGroup(2, surf_tags, 11, name = "bound_3")


gmsh.model.addPhysicalGroup(outExt[1][0], [outExt[1][1]],3, name="cube3")

# cube 4 - top left

xMin = 0
yMin = 1
zMin = 0

P1 = gmsh.model.occ.addPoint(xMin,yMin,zMin)
P2 = gmsh.model.occ.addPoint(xMin+X,yMin,zMin)
P3 = gmsh.model.occ.addPoint(xMin+X,yMin+Y,zMin)
P4 = gmsh.model.occ.addPoint(xMin,yMin+Y,zMin)

L1 = gmsh.model.occ.addLine(P1,P2)
L2 = gmsh.model.occ.addLine(P2,P3)
L3 = gmsh.model.occ.addLine(P3,P4)
L4 = gmsh.model.occ.addLine(P4,P1)

loop_4 = gmsh.model.occ.addCurveLoop([L1,L2,L3,L4])

surf_4 = gmsh.model.occ.addPlaneSurface([loop_4])

# dim tag is [(geo_size,tag_entity)]  numElements is a list!
outExt = gmsh.model.occ.extrude([(2, surf_4)], 0, 0, 1, numElements=[N4], recombine=recombineFlag)

gmsh.model.occ.synchronize()

gmsh.model.mesh.setRecombine(2, surf_4)

gmsh.model.mesh.set_transfinite_curve(L1,N4+1)
gmsh.model.mesh.set_transfinite_curve(L2,N4+1)
gmsh.model.mesh.set_transfinite_curve(L3,N4+1)
gmsh.model.mesh.set_transfinite_curve(L4,N4+1)

gmsh.model.mesh.set_transfinite_surface(loop_4)

gmsh.model.addPhysicalGroup(2, [outExt[2 + L2 - L1][1]], 7, name="interf_x_4")
gmsh.model.addPhysicalGroup(2, [outExt[2 + L1 - L1][1]], 8, name="interf_y_4")

surf_bc = [outExt[2+i - L1][1] for i in [L3,L4]]
surf_tags = [surf_4] + [outExt[0][1]] + surf_bc
gmsh.model.addPhysicalGroup(2, surf_tags, 12, name = "bound_4")


gmsh.model.addPhysicalGroup(outExt[1][0], [outExt[1][1]],4, name="cube4")


# generate mesh

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

# ... and save it to disk
outfile += '.vtk'

gmsh.write(outfile)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()