import gmsh

import os

import sys

gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 0)  # Disable terminal output
gmsh.option.setNumber("General.Verbosity", 0)  # Reduce verbosity to minimum

outfile = sys.argv[1]
N1 = int(sys.argv[2])
N2 = int(sys.argv[3])
elem_type = sys.argv[4]

gmsh.model.add('domain')


recombineFlag = True

if elem_type == 'tetra':
    recombineFlag = False

# Mesh.Format = 16 # msh output format
# Mesh.MshFileVersion = 2.2 # Version of the MSH file format to use

NX = int(0.25*N1)
NY = N1
NZ = N1
# This variable can then be used in the definition of Gmsh's simplest
# `elementary entity', a `Point'. A Point is uniquely identified by a tag (
# strictly positive integer; here `1') and defined by a list of four numbers:
# three coordinates (X, Y and Z) and the target mesh size (lc) close to the
# point:

X = 0.5
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
outL = gmsh.model.occ.extrude([(2, surfLeft)], 0, 0, 1, numElements=[NZ], recombine=recombineFlag)

gmsh.model.occ.synchronize()


if recombineFlag:
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

##############################################################################################

# RIGHT MESH

NX = int(0.25*N2)
NY = N2
NZ = N2

xMin = xMin+0.5

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
outR = gmsh.model.occ.extrude([(2,surfRight)],0,0,1,numElements=[NZ],recombine=recombineFlag)

gmsh.model.occ.synchronize()

if recombineFlag:
    gmsh.model.mesh.setRecombine(2, surfRight)


gmsh.model.mesh.set_transfinite_curve(L1r,NX+1)
gmsh.model.mesh.set_transfinite_curve(L2r,NY+1)
gmsh.model.mesh.set_transfinite_curve(L3r,NX+1)
gmsh.model.mesh.set_transfinite_curve(L4r,NY+1)

gmsh.model.mesh.set_transfinite_surface(surfRight)

gmsh.model.addPhysicalGroup(2, [outR[2 + P4r - P1r][1]], 3, name="right_interface")

surf_tags = [outR[2+i - P1r][1] for i in [P1r, P2r, P3r]]
surf_tags = [surfRight] + [outR[0][1]] + surf_tags
gmsh.model.addPhysicalGroup(2, surf_tags, 4, name = "right_lateralBound")


gmsh.model.addPhysicalGroup(outR[1][0], [outR[1][1]],2, name="right_mesh")



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

# Define the output folder and file
output_folder = os.path.join(os.getcwd(), 'Mesh', 'meshes')
os.makedirs(output_folder, exist_ok=True)

# Define the output file path
output_file = os.path.join(output_folder, outfile)
gmsh.write(output_file)

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