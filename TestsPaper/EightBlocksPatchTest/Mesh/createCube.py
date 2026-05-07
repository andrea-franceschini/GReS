import gmsh
import sys
import argparse




gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)

numb = int(sys.argv[1])
nel = int(sys.argv[2])
x0 = float(sys.argv[3])
y0 = float(sys.argv[4])
z0 = float(sys.argv[5])

parser = argparse.ArgumentParser()
parser.add_argument("numb", type=int)
parser.add_argument("nel", type=int)
parser.add_argument("x0", type=float)
parser.add_argument("y0", type=float)
parser.add_argument("z0", type=float)
args = parser.parse_args()

#print(args.arg1, args.arg2, args.arg3, args.arg4, args.arg5)


gmsh.model.add('cube')
# ---------------------- TOP LAYER-------------------------

# set number of layers for extrusion
X = 0.5
Y = 0.5
Z = 0.5

P1 = gmsh.model.occ.addPoint(x0,     y0,     z0)
P2 = gmsh.model.occ.addPoint(x0+X,   y0,     z0)
P3 = gmsh.model.occ.addPoint(x0+X,   y0+Y,   z0)
P4 = gmsh.model.occ.addPoint(x0,     y0+Y,   z0)

L1 = gmsh.model.occ.addLine(P1, P2)
L2 = gmsh.model.occ.addLine(P2, P3)
L3 = gmsh.model.occ.addLine(P3, P4)
L4 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L1, L2, L3, L4])
surf = gmsh.model.occ.addPlaneSurface([loop])

outExtrude = gmsh.model.occ.extrude([(2, surf)], 0, 0, Z, numElements=[nel], recombine=True)

gmsh.model.occ.synchronize()

gmsh.model.mesh.setRecombine(2, surf)

gmsh.model.mesh.set_transfinite_curve(L1,nel+1)
gmsh.model.mesh.set_transfinite_curve(L2,nel+1)
gmsh.model.mesh.set_transfinite_curve(L3,nel+1)
gmsh.model.mesh.set_transfinite_curve(L4,nel+1)

gmsh.model.mesh.set_transfinite_surface(surf)

# outExtrude[0][1] stores the tag of the surface opposite to the original surface after extrusion
gmsh.model.addPhysicalGroup(2, [surf], 1, name="face_1")
gmsh.model.addPhysicalGroup(2, [outExtrude[2][1]], 2, name="face_2")
gmsh.model.addPhysicalGroup(2, [outExtrude[3][1]], 3, name="face_3")
gmsh.model.addPhysicalGroup(2, [outExtrude[4][1]], 4, name="face_4")
gmsh.model.addPhysicalGroup(2, [outExtrude[5][1]], 5, name="face_5")
gmsh.model.addPhysicalGroup(2, [outExtrude[0][1]], 6, name="face_6")
gmsh.model.addPhysicalGroup(3, [outExtrude[1][1]], 1, name="vol"+str(numb))

# ---------------------- MESH SETTINGS -------------------------

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.write("cube"+str(numb)+".vtk")

gmsh.finalize()