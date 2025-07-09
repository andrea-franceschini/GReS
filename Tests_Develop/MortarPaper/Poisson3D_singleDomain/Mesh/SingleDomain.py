import gmsh

import gmsh
import sys

gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 0)  # Disable terminal output
gmsh.option.setNumber("General.Verbosity", 0)  # Reduce verbosity to minimum

outfile = sys.argv[1]
N = int(sys.argv[2])
elem_type = sys.argv[3] # tetra (no recombine), hexa (element order = 1), hexa27 (element order = 2)

gmsh.model.add('domain')

recombineFlag = True

if elem_type == 'tetra':
    recombineFlag = False

# Mesh.Format = 16 # msh output format
# Mesh.MshFileVersion = 2.2 # Version of the MSH file format to use

NX = 2*N
NY = N
NZ = N
# This variable can then be used in the definition of Gmsh's simplest
# `elementary entity', a `Point'. A Point is uniquely identified by a tag (
# strictly positive integer; here `1') and defined by a list of four numbers:
# three coordinates (X, Y and Z) and the target mesh size (lc) close to the
# point:

X = 2
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
L4 = gmsh.model.occ.addLine(P4, P1)

gmsh.model.occ.addCurveLoop([1,2,3,4],1)

surfLeft = gmsh.model.occ.addPlaneSurface([1])

# dim tag is [(geo_size,tag_entity)]  numElements is a list!
outL = gmsh.model.occ.extrude([(2, surfLeft)], 0, 0, 1, numElements=[NZ], recombine=recombineFlag)

gmsh.model.occ.synchronize()

if recombineFlag:
    gmsh.model.mesh.setRecombine(2, surfLeft)

gmsh.model.mesh.set_transfinite_curve(L1,NX+1)
gmsh.model.mesh.set_transfinite_curve(L2,NY+1)
gmsh.model.mesh.set_transfinite_curve(L3,NX+1)
gmsh.model.mesh.set_transfinite_curve(L4,NY+1)

gmsh.model.mesh.set_transfinite_surface(1)

surf_tags = [surfLeft] + [outL[0][1]] +  [e[1] for e in outL[2:]]

gmsh.model.addPhysicalGroup(2, surf_tags, 1, name = "lateralBound")


gmsh.model.addPhysicalGroup(outL[1][0], [outL[1][1]],1, name="domain")

if elem_type == "tetra" or elem_type == "hexa":
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
elif elem_type == "hexa27":
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
else:
    print("Unknown element type")

gmsh.model.mesh.generate(3)

# ... and save it to disk
outfile += '.vtk'
gmsh.write(outfile)

#if '-nopopup' not in sys.argv:
#    gmsh.fltk.run()

gmsh.finalize()



"""
Line(1) = {1, 2}
Line(2) = {2, 3}
Line(3) = {3, 4}
Line(4) = {4, 1}

Transfinite Line{1} = NX+1
Transfinite Line{2} = NY+1
Transfinite Line{3} = NX+1
Transfinite Line{4} = NY+1;

Curve Loop(1) = {1, 2, 3, 4}

Plane Surface(1) = {1}

Transfinite Surface {1} # structured grid
#Recombine Surface {1}; # using hexahedra

# extruding mesh along the z direction
Extrude {0,0,Z} { Surface{1}; Layers{NZ}};



Extrude {0, 0, Z} { Surface{1}; Layers{NZ};}
Physical Volume("cube",1) = {1}
Physical Surface("surfBound",1) = {1,13,17,21,25,26}

# Mesh second-order
#Mesh.ElementOrder = 1;
#Mesh.SecondOrderIncomplete = 0;  # <== ENSURES FULL HEXA-27

#Mesh 3;
#Save "domain.vtk";
"""