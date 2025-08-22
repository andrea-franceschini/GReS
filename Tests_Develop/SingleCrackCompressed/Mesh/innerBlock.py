import gmsh
import sys

from share.doc.gmsh.examples.api.plugin import numElements

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)

outfile = sys.argv[1]
NXb = int(sys.argv[2])
NYb = int(sys.argv[3])
NXt = int(sys.argv[4])
NYt = int(sys.argv[5])
nz = int(sys.argv[6])


gmsh.model.add('domainInner')
# ---------------------- TOP LAYER-------------------------

# set number of layers for extrusion

P1 = gmsh.model.occ.addPoint(-1.5,     -1,     -0.5)
P2 = gmsh.model.occ.addPoint(1.5,   -1,     -0.5)
P3 = gmsh.model.occ.addPoint(1.5,   0.5460,   -0.5)
P4 = gmsh.model.occ.addPoint(-1.5,     -0.5460,   -0.5)

L11 = gmsh.model.occ.addLine(P1, P2)
L12 = gmsh.model.occ.addLine(P2, P3)
L13 = gmsh.model.occ.addLine(P3, P4)
L14 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L11, L12, L13, L14])
r_bot = gmsh.model.occ.addPlaneSurface([loop])

e1 = gmsh.model.occ.extrude([(2, r_bot)], 0, 0, 1, numElements=[nz], recombine=True)

P1 = gmsh.model.occ.addPoint(-1.5,     -0.5460,     -0.5)
P2 = gmsh.model.occ.addPoint(1.5,   0.5460,     -0.5)
P3 = gmsh.model.occ.addPoint(1.5,   1,   -0.5)
P4 = gmsh.model.occ.addPoint(-1.5,     1,   -0.5)

L21 = gmsh.model.occ.addLine(P1, P2)
L22 = gmsh.model.occ.addLine(P2, P3)
L23 = gmsh.model.occ.addLine(P3, P4)
L24 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L21, L22, L23, L24])
r_top = gmsh.model.occ.addPlaneSurface([loop])

e2 = gmsh.model.occ.extrude([(2, r_top)], 0, 0, 1, numElements=[nz], recombine=True)

gmsh.model.occ.synchronize()

gmsh.model.mesh.set_transfinite_curve(L11,NXb)
gmsh.model.mesh.set_transfinite_curve(L12,NYb)
gmsh.model.mesh.set_transfinite_curve(L13,NXb)
gmsh.model.mesh.set_transfinite_curve(L14,NYb)

gmsh.model.mesh.set_transfinite_surface(r_bot)


gmsh.model.mesh.set_transfinite_curve(L21,NXt)
gmsh.model.mesh.set_transfinite_curve(L22,NYt)
gmsh.model.mesh.set_transfinite_curve(L23,NXt)
gmsh.model.mesh.set_transfinite_curve(L24,NYt)

gmsh.model.mesh.set_transfinite_surface(r_top)

gmsh.model.mesh.setRecombine(2, r_bot)
gmsh.model.mesh.setRecombine(2, r_top)

# define physical surfaces

outer_interf = [e1[2][1],e1[3][1],e1[5][1],e2[3][1],e2[4][1],e2[5][1]]
inner_interf_top = [e2[2][1]]
inner_interf_bot = [e1[4][1]]

z_fixed = [r_bot,r_top,e1[0][1],e2[0][1]]

vol = [e1[1][1],e2[1][1]]

gmsh.model.addPhysicalGroup(2,outer_interf,1,'outer_interf')
gmsh.model.addPhysicalGroup(2,inner_interf_top,2,'inner_interf_top')
gmsh.model.addPhysicalGroup(2,inner_interf_bot,3,'inner_interf_bot')
gmsh.model.addPhysicalGroup(2,z_fixed,4,'z_fixed')
gmsh.model.addPhysicalGroup(3,vol,1,'inner_block')


# ---------------------- MESH SETTINGS -------------------------

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.write(outfile+".vtk")

gmsh.finalize()