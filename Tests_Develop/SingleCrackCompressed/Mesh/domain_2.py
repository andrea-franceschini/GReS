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

c = 0.15
c2 = 0.88

gmsh.model.add('domainOuter')

Nout = int(sys.argv[7])
Nout_inter = int(sys.argv[8])



# define coordinates of model
cOut = 100 # size of outer grid
xIn = 2 # X size of inner grid
yIn = 1 # Y size of inner grid
slope = 20 # crack slope in deg
rad = math.radians(slope)

k = (xIn-math.cos(rad))/(2*xIn)
NX_tip = round(k*NXb)
NX_center = NXb-NX_tip


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

gmsh.model.mesh.set_transfinite_curve(L31,10*Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L32,Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L33,10*Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L34,Nout_inter)

gmsh.model.mesh.set_transfinite_curve(L41,10*Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L42,Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L43,10*Nout_inter)
gmsh.model.mesh.set_transfinite_curve(L44,Nout_inter)

for s in [r1, r2, r3, r4]:
    gmsh.model.mesh.set_transfinite_surface(s)


# bottom block (divided in 3 parts to distinguish fault from stick regions )
P1 = gmsh.model.occ.addPoint(-xIn,     -yIn,     -0.5)
P2 = gmsh.model.occ.addPoint(-math.cos(rad),   -yIn,     -0.5)
P3 = gmsh.model.occ.addPoint(-math.cos(rad),   -math.sin(rad),   -0.5)
P4 = gmsh.model.occ.addPoint(-xIn,     -xIn*math.tan(rad),   -0.5)

L11 = gmsh.model.occ.addLine(P1, P2)
L12 = gmsh.model.occ.addLine(P2, P3)
L13 = gmsh.model.occ.addLine(P3, P4)
L14 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L11, L12, L13, L14])
r_bot_1 = gmsh.model.occ.addPlaneSurface([loop])

e1_bot = gmsh.model.occ.extrude([(2, r_bot_1)], 0, 0, 1, numElements=[nz], recombine=True)


P1 = gmsh.model.occ.addPoint(-math.cos(rad),     -yIn,     -0.5)
P2 = gmsh.model.occ.addPoint(math.cos(rad),   -yIn,     -0.5)
P3 = gmsh.model.occ.addPoint(math.cos(rad),   math.sin(rad),   -0.5)
P4 = gmsh.model.occ.addPoint(-math.cos(rad),     -math.sin(rad),   -0.5)

L21 = gmsh.model.occ.addLine(P1, P2)
L22 = gmsh.model.occ.addLine(P2, P3)
L23 = gmsh.model.occ.addLine(P3, P4)
L24 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L21, L22, L23, L24])
r_bot_2 = gmsh.model.occ.addPlaneSurface([loop])

e2_bot = gmsh.model.occ.extrude([(2, r_bot_2)], 0, 0, 1, numElements=[nz], recombine=True)


P1 = gmsh.model.occ.addPoint(math.cos(rad),     -yIn,     -0.5)
P2 = gmsh.model.occ.addPoint(xIn,   -yIn,     -0.5)
P3 = gmsh.model.occ.addPoint(xIn,   xIn*math.tan(rad),   -0.5)
P4 = gmsh.model.occ.addPoint(math.cos(rad),     math.sin(rad),   -0.5)

L31 = gmsh.model.occ.addLine(P1, P2)
L32 = gmsh.model.occ.addLine(P2, P3)
L33 = gmsh.model.occ.addLine(P3, P4)
L34 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L31, L32, L33, L34])
r_bot_3 = gmsh.model.occ.addPlaneSurface([loop])

e3_bot = gmsh.model.occ.extrude([(2, r_bot_3)], 0, 0, 1, numElements=[nz], recombine=True)




P1 = gmsh.model.occ.addPoint(-xIn,     -xIn*math.tan(rad),     -0.5)
P2 = gmsh.model.occ.addPoint(xIn,   xIn*math.tan(rad),     -0.5)
P3 = gmsh.model.occ.addPoint(xIn,   yIn,   -0.5)
P4 = gmsh.model.occ.addPoint(-xIn,     yIn,   -0.5)

L41 = gmsh.model.occ.addLine(P1, P2)
L42 = gmsh.model.occ.addLine(P2, P3)
L43 = gmsh.model.occ.addLine(P3, P4)
L44 = gmsh.model.occ.addLine(P4, P1)

loop = gmsh.model.occ.addCurveLoop([L41, L42, L43, L44])
r_top = gmsh.model.occ.addPlaneSurface([loop])

e_top = gmsh.model.occ.extrude([(2, r_top)], 0, 0, 1, numElements=[nz], recombine=True)

gmsh.model.occ.synchronize()

gmsh.model.mesh.set_transfinite_curve(L11,NX_tip,"Progression",c2)
gmsh.model.mesh.set_transfinite_curve(L12,NYb)
gmsh.model.mesh.set_transfinite_curve(L13,NX_tip,"Progression",1/c2)
gmsh.model.mesh.set_transfinite_curve(L14,NYb)

gmsh.model.mesh.set_transfinite_surface(r_bot_1)


gmsh.model.mesh.set_transfinite_curve(L21,NX_center,"Bump",c)
gmsh.model.mesh.set_transfinite_curve(L22,NYb)
gmsh.model.mesh.set_transfinite_curve(L23,NX_center,"Bump",c)
gmsh.model.mesh.set_transfinite_curve(L24,NYb)

gmsh.model.mesh.set_transfinite_surface(r_bot_2)

gmsh.model.mesh.set_transfinite_curve(L31,NX_tip,"Progression",1/c2)
gmsh.model.mesh.set_transfinite_curve(L32,NYb)
gmsh.model.mesh.set_transfinite_curve(L33,NX_tip,"Progression",c2)
gmsh.model.mesh.set_transfinite_curve(L34,NYb)

gmsh.model.mesh.set_transfinite_surface(r_bot_3)

gmsh.model.mesh.set_transfinite_curve(L41,NXt)
gmsh.model.mesh.set_transfinite_curve(L42,NYt)
gmsh.model.mesh.set_transfinite_curve(L43,NXt)
gmsh.model.mesh.set_transfinite_curve(L44,NYt)

gmsh.model.mesh.set_transfinite_surface(r_top)


gmsh.model.mesh.setRecombine(2, r_bot_1)
gmsh.model.mesh.setRecombine(2, r_bot_2)
gmsh.model.mesh.setRecombine(2, r_bot_3)
gmsh.model.mesh.setRecombine(2, r_top)

# define physical surfaces

master_interf_1 = [e1_bot[2][1],e2_bot[2][1],e2_bot[2][1],e_top[4][1],15,17,20,22]
slave_interf_1 = [7,10]
master_interf_2 = [e1_bot[3][1],e1_bot[5][1],e3_bot[3][1],e3_bot[5][1],e_top[3][1],e_top[5][1]]
slave_interf_2 = [16,23,e2_bot[3][1],e2_bot[5][1]]
inner_interf_top = [e_top[2][1]]
inner_interf_bot_stick = [e1_bot[4][1],e3_bot[4][1]]
inner_interf_bot_slip = [e2_bot[4][1]]
z_fixed = [r_bot_1,r_bot_2,r_bot_3,r_top,e1_bot[0][1],e2_bot[0][1],e3_bot[0][1],e_top[0][1],1,2,3,4,9,14,19,24]
volInner = [e1_bot[1][1],e2_bot[1][1],e3_bot[1][1],e_top[1][1]]



gmsh.model.addPhysicalGroup(2, [8,13,18], 1, name="out_left")
gmsh.model.addPhysicalGroup(2, [6,11,21], 2, name="out_right")
gmsh.model.addPhysicalGroup(2, z_fixed, 3, name="z_fixed")
gmsh.model.addPhysicalGroup(2, master_interf_1, 4, name="interface_master_1")
gmsh.model.addPhysicalGroup(2, slave_interf_1, 5, name="interface_slave_1")
gmsh.model.addPhysicalGroup(2, master_interf_2, 6, name="interface_master_2")
gmsh.model.addPhysicalGroup(2, slave_interf_2, 7, name="interface_slave_2")
gmsh.model.addPhysicalGroup(2,inner_interf_top,8,'inner_interf_top')
gmsh.model.addPhysicalGroup(2,inner_interf_bot_stick,10,'inner_interf_bot_stick')
gmsh.model.addPhysicalGroup(2,inner_interf_bot_slip,9,'inner_interf_bot_slip')

gmsh.model.addPhysicalGroup(3, [1,2,3,4], 1, name="outer_block")
gmsh.model.addPhysicalGroup(3,volInner,2,'inner_block')


# ---------------------- MESH SETTINGS -------------------------

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.write(outfile+".vtk")

gmsh.finalize()