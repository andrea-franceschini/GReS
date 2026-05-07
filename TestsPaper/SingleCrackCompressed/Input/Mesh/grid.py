import gmsh
import sys
import math

from share.doc.gmsh.examples.api.plugin import numElements

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)

outfile = sys.argv[1]
szBot = float(sys.argv[2])
szBotFault = float(sys.argv[3])
szTop = float(sys.argv[4])
szTopFault = float(sys.argv[5])
nz = int(sys.argv[6])

gmsh.model.add('domainOuter')


# define coordinates of model
lOut = 50 # half size of outer edge

# half size of the fault
lIn = 1


phi = 20 # crack slope in deg
phi = math.radians(phi)

# bottom
P1 = gmsh.model.occ.addPoint(-lOut,     -lOut,     -0.5,szBot)
P2 = gmsh.model.occ.addPoint(0,     -lOut,     -0.5,szBot)
P3 = gmsh.model.occ.addPoint(lOut,   -lOut,     -0.5,szBot)
P4 = gmsh.model.occ.addPoint(lOut,   0,     -0.5,szBot)
P5 = gmsh.model.occ.addPoint(lOut,   lOut*math.tan(phi),   -0.5,szBot)
P6 = gmsh.model.occ.addPoint(lIn*math.cos(phi),   lIn*math.sin(phi),   -0.5,szBotFault)
P7 = gmsh.model.occ.addPoint(-lIn*math.cos(phi),   -lIn*math.sin(phi),   -0.5, szBotFault)
P8 = gmsh.model.occ.addPoint(-lOut,   -lOut*math.tan(phi),   -0.5, szBot)



L11 = gmsh.model.occ.addLine(P1, P2)
L12 = gmsh.model.occ.addLine(P2, P3)
L13 = gmsh.model.occ.addLine(P3, P4)
L14 = gmsh.model.occ.addLine(P4, P5)
L15 = gmsh.model.occ.addLine(P5, P6)
L16 = gmsh.model.occ.addLine(P6, P7)
L17 = gmsh.model.occ.addLine(P7, P8)
L18 = gmsh.model.occ.addLine(P8, P1)

loop = gmsh.model.occ.addCurveLoop([L11, L12, L13, L14, L15, L16, L17, L18])
bottom = gmsh.model.occ.addPlaneSurface([loop])

# top
P1 = gmsh.model.occ.addPoint(-lOut,   -lOut*math.tan(phi),   -0.5,szTop)
P2 = gmsh.model.occ.addPoint(-lIn*math.cos(phi),   -lIn*math.sin(phi),   -0.5,szTopFault)
P3 = gmsh.model.occ.addPoint(lIn*math.cos(phi),   lIn*math.sin(phi),   -0.5,szTopFault)
P4 = gmsh.model.occ.addPoint(lOut,   lOut*math.tan(phi),   -0.5,szTop)
P5 = gmsh.model.occ.addPoint(lOut,   lOut,   -0.5,szTop)
P6 = gmsh.model.occ.addPoint(0,   lOut,   -0.5,szTop)
P7 = gmsh.model.occ.addPoint(-lOut,   lOut,   -0.5,szTop)
P8 = gmsh.model.occ.addPoint(-lOut,   0,   -0.5,szTop)

L21 = gmsh.model.occ.addLine(P1, P2)
L22 = gmsh.model.occ.addLine(P2, P3)
L23 = gmsh.model.occ.addLine(P3, P4)
L24 = gmsh.model.occ.addLine(P4, P5)
L25 = gmsh.model.occ.addLine(P5, P6)
L26 = gmsh.model.occ.addLine(P6, P7)
L27 = gmsh.model.occ.addLine(P7, P8)
L28 = gmsh.model.occ.addLine(P8, P1)

loop = gmsh.model.occ.addCurveLoop([L21, L22, L23, L24, L25, L26, L27, L28])
top = gmsh.model.occ.addPlaneSurface([loop])


extBot = gmsh.model.occ.extrude([(2, bottom)], 0, 0, 1, numElements=[nz], recombine=True)
extTop = gmsh.model.occ.extrude([(2, top)], 0, 0, 1, numElements=[nz], recombine=True)


gmsh.model.occ.synchronize()


for s in [top, bottom]:
    gmsh.model.mesh.setRecombine(2, s)


#gmsh.model.mesh.set_transfinite_curve(L11,Nbot)
#gmsh.model.mesh.set_transfinite_curve(L12,Nbot)
#gmsh.model.mesh.set_transfinite_curve(L13,Nbot)
#gmsh.model.mesh.set_transfinite_curve(L14,Nfault,"Bump",k)
#gmsh.model.mesh.set_transfinite_curve(L15,Nbot)
#gmsh.model.mesh.set_transfinite_curve(L16,Nbot)
#
#gmsh.model.mesh.set_transfinite_curve(L21,Ntop)
#gmsh.model.mesh.set_transfinite_curve(L22,Ntop)
#gmsh.model.mesh.set_transfinite_curve(L23,Ntop)
#gmsh.model.mesh.set_transfinite_curve(L24,Nfault_top,"Bump",k)
#gmsh.model.mesh.set_transfinite_curve(L25,Ntop)
#gmsh.model.mesh.set_transfinite_curve(L26,Ntop)


# define physical surfaces

z_fixed = [bottom,top,extBot[0][1],extTop[0][1]]


interf_stick_slave = [extBot[6][1],extBot[8][1]]
interf_stick_master = [extTop[2][1],extTop[4][1]]

fault_slave = [extBot[7][1]]
fault_master = [extTop[3][1]]

left = [extBot[9][1],extTop[8][1],extTop[9][1]]
right = [extBot[4][1],extBot[5][1],extTop[5][1]]

gmsh.model.addPhysicalGroup(2, z_fixed, 1, name="z_fixed")
gmsh.model.addPhysicalGroup(2, left, 2, name="load_left")
gmsh.model.addPhysicalGroup(2, right, 3, name="load_right")
gmsh.model.addPhysicalGroup(2, interf_stick_master, 4, name="stick_master")
gmsh.model.addPhysicalGroup(2, interf_stick_slave, 5, name="stick_slave")
gmsh.model.addPhysicalGroup(2, fault_master, 6, name="fault_master")
gmsh.model.addPhysicalGroup(2, fault_slave, 7, name="fault_slave")

gmsh.model.addPhysicalGroup(3,[1,2],1,'domain')


# ---------------------- MESH SETTINGS -------------------------

gmsh.option.setNumber("Mesh.ElementOrder", 1)

gmsh.model.mesh.generate(3)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.write(outfile+".vtk")

gmsh.finalize()