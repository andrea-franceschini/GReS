import gmsh

import gmsh
import sys

gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 0)  # Disable terminal output
gmsh.option.setNumber("General.Verbosity", 0)  # Reduce verbosity to minimum

outfile = sys.argv[1]
N2 = int(sys.argv[2])

gmsh.model.add('domain')

# Mesh.Format = 16 # msh output format
# Mesh.MshFileVersion = 2.2 # Version of the MSH file format to use

NX = N2
NY = N2
NZ = N2
# This variable can then be used in the definition of Gmsh's simplest
# `elementary entity', a `Point'. A Point is uniquely identified by a tag (
# strictly positive integer; here `1') and defined by a list of four numbers:
# three coordinates (X, Y and Z) and the target mesh size (lc) close to the
# point:

X = 1
Y = 1
Z = 1
xMin = 1
yMin = 0
zMin = 0

P1r = gmsh.model.occ.addPoint(xMin,yMin,zMin)
P2r = gmsh.model.occ.addPoint(xMin+X,yMin,zMin)
P3r = gmsh.model.occ.addPoint(xMin+X,yMin+Y,zMin)
P4r = gmsh.model.occ.addPoint(xMin,yMin+Y,zMin)

L1r = gmsh.model.occ.addLine(P1r,P2r)
L2r = gmsh.model.occ.addLine(P2r,P3r)
L3r = gmsh.model.occ.addLine(P3r,P4r)
L4r = gmsh.model.occ.addLine(P4r, P1r)

curveRight = gmsh.model.occ.addCurveLoop([L1r,L2r,L3r,L4r])

surfRight = gmsh.model.occ.addPlaneSurface([curveRight])

# dim tag is [(geo_size,tag_entity)]  numElements is a list!
outR = gmsh.model.occ.extrude([(2,surfRight)],0,0,1,numElements=[NZ],recombine=True)

gmsh.model.occ.synchronize()

gmsh.model.mesh.setRecombine(2, surfRight)

gmsh.model.mesh.set_transfinite_curve(L1r,NX+1)
gmsh.model.mesh.set_transfinite_curve(L2r,NY+1)
gmsh.model.mesh.set_transfinite_curve(L3r,NX+1)
gmsh.model.mesh.set_transfinite_curve(L4r,NY+1)

gmsh.model.mesh.set_transfinite_surface(surfRight)

gmsh.model.addPhysicalGroup(2, [outR[2 + P4r - P1r][1]], 1, name="right_interface")

surf_tags = [outR[2+i - P1r][1] for i in [P1r, P2r, P3r]]
surf_tags = [surfRight] + [outR[0][1]] + surf_tags
gmsh.model.addPhysicalGroup(2, surf_tags, 2, name = "right_lateralBound")


gmsh.model.addPhysicalGroup(outR[1][0], [outR[1][1]],1, name="right_mesh")





gmsh.model.mesh.generate(3)

# ... and save it to disk
outfile += '.vtk'
gmsh.write(outfile)

#if '-nopopup' not in sys.argv:
#    gmsh.fltk.run()

gmsh.finalize()

