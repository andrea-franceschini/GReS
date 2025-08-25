import gmsh
import sys
import math

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

c = 0.7
c2 = 0.85

gmsh.model.add('domainOuter')

Nout = 11
Nout_inter = 5


# define coordinates of model
cOut = 15 # size of outer grid
xIn = 2 # X size of inner grid
yIn = 1 # Y size of inner grid
slope = 20 # crack slope in deg
rad = math.radians(slope)

# rectangle 1
P1 = gmsh.model.occ.addPoint(-cOut,     -cOut,     -0.5)
P2 = gmsh.model.occ.addPoint(cOut,   -cOut,     -0.5)
P3 = gmsh.model.occ.addPoint(cOut,   -yIn,   -0.5)
P4 = gmsh.model.occ.addPoint(-cOut,     -yIn,   -0.5)

L11 = gmsh.model.occ.addLine(P1, P2)
L12 = gmsh.model.occ.addLine(P2, P3)
L13 = gmsh.model.occ.addLine(P3, P4)
L14 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L11, L12, L13, L14])
r1 = gmsh.model.occ.addPlaneSurface([loop])

# rectangle 2
P1 = gmsh.model.occ.addPoint(-cOut,     yIn,     -0.5)
P2 = gmsh.model.occ.addPoint(cOut,   yIn,     -0.5)
P3 = gmsh.model.occ.addPoint(cOut,   cOut,   -0.5)
P4 = gmsh.model.occ.addPoint(-cOut,     cOut,   -0.5)

L21 = gmsh.model.occ.addLine(P1, P2)
L22 = gmsh.model.occ.addLine(P2, P3)
L23 = gmsh.model.occ.addLine(P3, P4)
L24 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L21, L22, L23, L24])
r2 = gmsh.model.occ.addPlaneSurface([loop])

# rectangle 3
P1 = gmsh.model.occ.addPoint(-cOut,     -yIn,     -0.5)
P2 = gmsh.model.occ.addPoint(-xIn,   -yIn,     -0.5)
P3 = gmsh.model.occ.addPoint(-xIn,   yIn,   -0.5)
P4 = gmsh.model.occ.addPoint(-cOut,     yIn,   -0.5)

L31 = gmsh.model.occ.addLine(P1, P2)
L32 = gmsh.model.occ.addLine(P2, P3)
L33 = gmsh.model.occ.addLine(P3, P4)
L34 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L31, L32, L33, L34])
r3 = gmsh.model.occ.addPlaneSurface([loop])


# rectangle 4
P1 = gmsh.model.occ.addPoint(xIn,     -yIn,     -0.5)
P2 = gmsh.model.occ.addPoint(cOut,   -yIn,     -0.5)
P3 = gmsh.model.occ.addPoint(cOut,   yIn,   -0.5)
P4 = gmsh.model.occ.addPoint(xIn,     yIn,   -0.5)

L41 = gmsh.model.occ.addLine(P1, P2)
L42 = gmsh.model.occ.addLine(P2, P3)
L43 = gmsh.model.occ.addLine(P3, P4)
L44 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L41, L42, L43, L44])
r4 = gmsh.model.occ.addPlaneSurface([loop])


outExtrude = gmsh.model.occ.extrude([(2, r1), (2, r2), (2, r3), (2, r4)], 0, 0, 1, numElements=[nz], recombine=True)




gmsh.model.occ.synchronize()


for s in [r1, r2, r3, r4]:
    gmsh.model.mesh.setRecombine(2, s)


gmsh.model.mesh.set_transfinite_curve(L11,Nout)
gmsh.model.mesh.set_transfinite_curve(L12,round(0.5*Nout))
gmsh.model.mesh.set_transfinite_curve(L13,Nout)
gmsh.model.mesh.set_transfinite_curve(L14,round(0.5*Nout))

gmsh.model.mesh.set_transfinite_curve(L21,Nout)
gmsh.model.mesh.set_transfinite_curve(L22,round(0.5*Nout))
gmsh.model.mesh.set_transfinite_curve(L23,Nout)
gmsh.model.mesh.set_transfinite_curve(L24,round(0.5*Nout))

gmsh.model.mesh.set_transfinite_curve(L31,2*Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L32,Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L33,2*Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L34,Nout_inter)

gmsh.model.mesh.set_transfinite_curve(L41,2*Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L42,Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L43,2*Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L44,Nout_inter)

for s in [r1, r2, r3, r4]:
    gmsh.model.mesh.set_transfinite_surface(s)



P1 = gmsh.model.occ.addPoint(-xIn,     -yIn,     -0.5)
P2 = gmsh.model.occ.addPoint(xIn,   -yIn,     -0.5)
P3 = gmsh.model.occ.addPoint(xIn,   xIn*math.tan(rad),   -0.5)
P4 = gmsh.model.occ.addPoint(-xIn,     -xIn*math.tan(rad),   -0.5)

L11 = gmsh.model.occ.addLine(P1, P2)
L12 = gmsh.model.occ.addLine(P2, P3)
L13 = gmsh.model.occ.addLine(P3, P4)
L14 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L11, L12, L13, L14])
r_bot = gmsh.model.occ.addPlaneSurface([loop])

e1 = gmsh.model.occ.extrude([(2, r_bot)], 0, 0, 1, numElements=[nz], recombine=True)

P1 = gmsh.model.occ.addPoint(-xIn, -xIn*math.tan(rad),     -0.5)
P2 = gmsh.model.occ.addPoint(xIn,   xIn*math.tan(rad),     -0.5)
P3 = gmsh.model.occ.addPoint(xIn,   yIn,   -0.5)
P4 = gmsh.model.occ.addPoint(-xIn,     yIn,   -0.5)

L21 = gmsh.model.occ.addLine(P1, P2)
L22 = gmsh.model.occ.addLine(P2, P3)
L23 = gmsh.model.occ.addLine(P3, P4)
L24 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L21, L22, L23, L24])
r_top = gmsh.model.occ.addPlaneSurface([loop])

e2 = gmsh.model.occ.extrude([(2, r_top)], 0, 0, 1, numElements=[nz], recombine=True)

gmsh.model.occ.synchronize()

gmsh.model.mesh.set_transfinite_curve(L11,NXb,"Bump",c)
gmsh.model.mesh.set_transfinite_curve(L12,NYb,"Progression",c2)
gmsh.model.mesh.set_transfinite_curve(L13,NXb,"Bump",c)
gmsh.model.mesh.set_transfinite_curve(L14,NYb,"Progression",1/c2)

gmsh.model.mesh.set_transfinite_surface(r_bot)


gmsh.model.mesh.set_transfinite_curve(L21,NXt,"Bump",c)
gmsh.model.mesh.set_transfinite_curve(L22,NYt,"Progression",1/c2)
gmsh.model.mesh.set_transfinite_curve(L23,NXt,"Bump",c)
gmsh.model.mesh.set_transfinite_curve(L24,NYt,"Progression",c2)

gmsh.model.mesh.set_transfinite_surface(r_top)

gmsh.model.mesh.setRecombine(2, r_bot)
gmsh.model.mesh.setRecombine(2, r_top)

# define physical surfaces

master_interf_1 = [e1[2][1],e2[4][1],15,17,20,22]
slave_interf_1 = [7,10]
master_interf_2 = [e1[3][1],e1[5][1],e2[3][1],e2[5][1]]
slave_interf_2 = [16,23]
inner_interf_top = [e2[2][1]]
inner_interf_bot = [e1[4][1]]
z_fixed = [r_bot,r_top,e1[0][1],e2[0][1],1,2,3,4,9,14,19,24]
volInner = [e1[1][1],e2[1][1]]



gmsh.model.addPhysicalGroup(2, [8,13,18], 1, name="out_left")
gmsh.model.addPhysicalGroup(2, [6,11,21], 2, name="out_right")
gmsh.model.addPhysicalGroup(2, z_fixed, 3, name="z_fixed")
gmsh.model.addPhysicalGroup(2, master_interf_1, 4, name="interface_master_1")
gmsh.model.addPhysicalGroup(2, slave_interf_1, 5, name="interface_slave_1")
gmsh.model.addPhysicalGroup(2, master_interf_2, 6, name="interface_master_2")
gmsh.model.addPhysicalGroup(2, slave_interf_2, 7, name="interface_slave_2")
gmsh.model.addPhysicalGroup(2,inner_interf_top,8,'inner_interf_top')
gmsh.model.addPhysicalGroup(2,inner_interf_bot,9,'inner_interf_bot')

gmsh.model.addPhysicalGroup(3, [1,2,3,4], 1, name="outer_block")
gmsh.model.addPhysicalGroup(3,volInner,2,'inner_block')


# ---------------------- MESH SETTINGS -------------------------

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

#if '-nopopup' not in sys.argv:
#    gmsh.fltk.run()

gmsh.write(outfile+".vtk")

gmsh.finalize()