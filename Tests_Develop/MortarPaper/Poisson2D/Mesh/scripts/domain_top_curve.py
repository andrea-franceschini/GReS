import gmsh

import os

import sys

gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 0)  # Disable terminal output
gmsh.option.setNumber("General.Verbosity", 0)  # Reduce verbosity to minimum

outfile = sys.argv[1]
NX = int(sys.argv[2])
NY = int(sys.argv[3])

elem_type = sys.argv[4]

gmsh.model.add('domain')


recombineFlag = True

if elem_type == 'tetra' or elem_type == "tetra6":
    recombineFlag = False

# Mesh.Format = 16 # msh output format
# Mesh.MshFileVersion = 2.2 # Version of the MSH file format to use


# This variable can then be used in the definition of Gmsh's simplest
# `elementary entity', a `Point'. A Point is uniquely identified by a tag (
# strictly positive integer; here `1') and defined by a list of four numbers:
# three coordinates (X, Y and Z) and the target mesh size (lc) close to the
# point:

X = 1
Y = 0.5

xMin = 0
yMin = 0.5
zMin = 0

P1l = gmsh.model.occ.addPoint(xMin,yMin,zMin)
P2l = gmsh.model.occ.addPoint(xMin+X,yMin,zMin)
P3l = gmsh.model.occ.addPoint(xMin+X,yMin+Y,zMin)
P4l = gmsh.model.occ.addPoint(xMin,yMin+Y,zMin)
P5 = gmsh.model.occ.addPoint(xMin+0.5*X,0.1,zMin)

L1l = gmsh.model.occ.addCircleArc(P1l,P5,P2l)
L2l = gmsh.model.occ.addLine(P2l,P3l)
L3l = gmsh.model.occ.addLine(P3l,P4l)
L4l = gmsh.model.occ.addLine(P4l, P1l)

gmsh.model.occ.addCurveLoop([1,2,3,4],1)

bot = gmsh.model.occ.addPlaneSurface([1])

gmsh.model.occ.synchronize()


if recombineFlag:
    gmsh.model.mesh.setRecombine(2, bot)

gmsh.model.mesh.set_transfinite_curve(L1l,NX+1)
gmsh.model.mesh.set_transfinite_curve(L2l,NY+1)
gmsh.model.mesh.set_transfinite_curve(L3l,NX+1)
gmsh.model.mesh.set_transfinite_curve(L4l,NY+1)

gmsh.model.mesh.set_transfinite_surface(1)

gmsh.model.addPhysicalGroup(1, [L1l], 1, name="interf")
gmsh.model.addPhysicalGroup(1, [L2l,L3l,L4l], 2, name="bound")


gmsh.model.addPhysicalGroup(2, [1], 1, name="topMesh")

##############################################################################################


if elem_type == "tetra" or elem_type == "hexa":
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
elif elem_type == "hexa27" or "tetra6":
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
else:
    print("Unknown element type")


gmsh.model.mesh.generate(2)

# ... and save it to disk
outfile += '.vtk'

# Define the output folder and file
output_folder = os.path.join(os.getcwd(), 'Mesh', 'meshes')
os.makedirs(output_folder, exist_ok=True)

# Define the output file path
output_file = os.path.join(output_folder, outfile)
gmsh.write(output_file)

#if '-nopopup' not in sys.argv:
#   gmsh.fltk.run()

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