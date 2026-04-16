import gmsh

import sys

gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 0)  # Disable terminal output
gmsh.option.setNumber("General.Verbosity", 0)  # Reduce verbosity to minimum

outfile = sys.argv[1]
N1 = int(sys.argv[2])

gmsh.model.add('domain')

# Mesh.Format = 16 # msh output format
# Mesh.MshFileVersion = 2.2 # Version of the MSH file format to use

NX = N1
NY = N1
NZ = N1
# This variable can then be used in the definition of Gmsh's simplest
# `elementary entity', a `Point'. A Point is uniquely identified by a tag (
# strictly positive integer; here `1') and defined by a list of four numbers:
# three coordinates (X, Y and Z) and the target mesh size (lc) close to the
# point:

X = 1
Y = 1
Z = 1
xMin = 0
yMin = 0
zMin = 0

P1l = gmsh.model.occ.addPoint(xMin,yMin,zMin)
P2l = gmsh.model.occ.addPoint(xMin+X,yMin,zMin)
P3l = gmsh.model.occ.addPoint(xMin+X,yMin+Y,zMin)
P4l = gmsh.model.occ.addPoint(xMin,yMin+Y,zMin)

L1l = gmsh.model.occ.addLine(P1l,P2l)
L2l = gmsh.model.occ.addLine(P2l,P3l)
L3l = gmsh.model.occ.addLine(P3l,P4l)
L4l = gmsh.model.occ.addLine(P4l, P1l)

gmsh.model.occ.addCurveLoop([1,2,3,4],1)

surfLeft = gmsh.model.occ.addPlaneSurface([1])

# dim tag is [(geo_size,tag_entity)]  numElements is a list!
outL = gmsh.model.occ.extrude([(2, surfLeft)], 0, 0, 1, numElements=[NZ], recombine=True)

gmsh.model.occ.synchronize()

gmsh.model.mesh.setRecombine(2, surfLeft)

gmsh.model.mesh.set_transfinite_curve(L1l,NX+1)
gmsh.model.mesh.set_transfinite_curve(L2l,NY+1)
gmsh.model.mesh.set_transfinite_curve(L3l,NX+1)
gmsh.model.mesh.set_transfinite_curve(L4l,NY+1)

gmsh.model.mesh.set_transfinite_surface(1)

gmsh.model.addPhysicalGroup(2, [outL[2+ P2l - P1l][1]], 1, name="left_interface")

surf_tags = [outL[2+i - P1l][1] for i in [P1l, P3l, P4l]]
surf_tags = [surfLeft] + [outL[0][1]] + surf_tags
gmsh.model.addPhysicalGroup(2, surf_tags, 2, name = "left_lateralBound")


gmsh.model.addPhysicalGroup(outL[1][0], [outL[1][1]],1, name="left_mesh")

gmsh.model.mesh.generate(3)

# ... and save it to disk
outfile += '.vtk'
gmsh.write(outfile)

#if '-nopopup' not in sys.argv:
#    gmsh.fltk.run()

gmsh.finalize()
