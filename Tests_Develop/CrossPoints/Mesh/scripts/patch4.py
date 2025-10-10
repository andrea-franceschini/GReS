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

bot_1 = gmsh.model.occ.addPlaneSurface([loop_1])

# dim tag is [(geo_size,tag_entity)]  numElements is a list!
outExt = gmsh.model.occ.extrude([(2, bot_1)], 0, 0, 1, numElements=[N1], recombine=recombineFlag)


gmsh.model.occ.synchronize()

gmsh.model.mesh.setRecombine(2, bot_1)

gmsh.model.mesh.set_transfinite_curve(L1,1+1)
gmsh.model.mesh.set_transfinite_curve(L2,1+1)
gmsh.model.mesh.set_transfinite_curve(L3,1+1)
gmsh.model.mesh.set_transfinite_curve(L4,1+1)

gmsh.model.mesh.set_transfinite_surface(loop_1)

top_1 = outExt[0][1]

interf_1_x = outExt[2 + L2 - L1][1]
interf_1_y = outExt[2 + L3 - L1][1]

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

bot_2 = gmsh.model.occ.addPlaneSurface([loop_2])

# dim tag is [(geo_size,tag_entity)]  numElements is a list!
outExt = gmsh.model.occ.extrude([(2, bot_2)], 0, 0, 1, numElements=[N2], recombine=recombineFlag)


gmsh.model.occ.synchronize()

gmsh.model.mesh.setRecombine(2, bot_2)

gmsh.model.mesh.set_transfinite_curve(L1,1+1)
gmsh.model.mesh.set_transfinite_curve(L2,1+1)
gmsh.model.mesh.set_transfinite_curve(L3,1+1)
gmsh.model.mesh.set_transfinite_curve(L4,1+1)

gmsh.model.mesh.set_transfinite_surface(loop_2)


top_2 = outExt[0][1]

interf_2_x = outExt[2 + L4 - L1][1]
interf_2_y = outExt[2 + L3 - L1][1]

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

bot_3 = gmsh.model.occ.addPlaneSurface([loop_3])

# dim tag is [(geo_size,tag_entity)]  numElements is a list!
outExt = gmsh.model.occ.extrude([(2, bot_3)], 0, 0, 1, numElements=[N3], recombine=recombineFlag)

top_3 = outExt[0][1]

gmsh.model.occ.synchronize()

gmsh.model.mesh.setRecombine(2, bot_3)

gmsh.model.mesh.set_transfinite_curve(L1,1+1)
gmsh.model.mesh.set_transfinite_curve(L2,1+1)
gmsh.model.mesh.set_transfinite_curve(L3,1+1)
gmsh.model.mesh.set_transfinite_curve(L4,1+1)

gmsh.model.mesh.set_transfinite_surface(loop_3)

interf_3_x = outExt[2 + L4 - L1][1]
interf_3_y = outExt[2 + L1 - L1][1]

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

bot_4 = gmsh.model.occ.addPlaneSurface([loop_4])

# dim tag is [(geo_size,tag_entity)]  numElements is a list!
outExt = gmsh.model.occ.extrude([(2, bot_4)], 0, 0, 1, numElements=[N4], recombine=recombineFlag)

top_4 = outExt[0][1]

gmsh.model.occ.synchronize()

gmsh.model.mesh.setRecombine(2, bot_4)

gmsh.model.mesh.set_transfinite_curve(L1,1+1)
gmsh.model.mesh.set_transfinite_curve(L2,1+1)
gmsh.model.mesh.set_transfinite_curve(L3,1+1)
gmsh.model.mesh.set_transfinite_curve(L4,1+1)

gmsh.model.mesh.set_transfinite_surface(loop_4)

interf_4_x = outExt[2 + L2 - L1][1]
interf_4_y = outExt[2 + L1 - L1][1]



gmsh.model.addPhysicalGroup(outExt[1][0], [outExt[1][1]],4, name="cube4")

# set slave interface 1
gmsh.model.addPhysicalGroup(2, [interf_1_x,interf_1_y], 1, name="slave_1")

# set slave interface 2
gmsh.model.addPhysicalGroup(2, [interf_3_x,interf_3_y], 2, name="slave_2")

# set master interface 1
gmsh.model.addPhysicalGroup(2, [interf_2_x,interf_4_y], 3, name="master_1")

# set master interface 2
gmsh.model.addPhysicalGroup(2, [interf_4_x,interf_2_y], 4, name="master_2")

# set bottom surface (fixed)
gmsh.model.addPhysicalGroup(2, [bot_1,bot_2,bot_3,bot_4], 5, name="bottom_fixed")

# set top surface (loaded)
gmsh.model.addPhysicalGroup(2, [top_1,top_2,top_3,top_4], 6, name="top_load")

# generate mesh

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

# ... and save it to disk
outfile += '.vtk'

gmsh.write(outfile)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()