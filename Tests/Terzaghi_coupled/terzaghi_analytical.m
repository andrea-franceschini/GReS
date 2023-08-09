%clear; clc;
%Script for Terzaghi analytical solution, uses intrinsic permeability.
%Found in JCP_2018.


plotMode = false;
%true ---> outside from GReS main
%false ----> inside GReS main (error checking)

%%%%%%%%%%%%%      Fluid and poroelastic properties       %%%%%%%%%%%%%%%%
rho = 998.2; %[kg/m3] fluid density
k = 1.0e-12; %[m^2] intrinsic permeability
khydro = 1.0e-5;
porosity = 0.375; 
dyn = 1.00e-6; %[kPa*s]
cf = 4.4e-7 ; % [kPa^-1] fluid compressibility
lambda = 40000; %[kPa] costante lamè
G = 40000; %[kPa] lamè constant
biot = 1.0; %coefficiente di Biot
M = (porosity*cf)^-1; %Biot Modulus, assuming cbr=0
Ku = lambda + 2*(G/3) + biot^2*M;
B = biot*M/Ku;
nu = lambda/(2*(lambda+G));
nuU = (3*nu+biot*B*(1-2*nu))/(3-biot*B*(1-2*nu)); %undrained Poisson Coefficient
c = 2*k*G*(1-nu)*(nuU-nu)/(dyn*biot^2*(1-nuU)*(1-2*nu)^2);
gamma = B*(1+nuU)/(3*(1-nuU));
cm = (lambda+2*G)^-1;
%c_alternative = khydro/ (9.81*(M^-1+biot^2*cm));
L = 10; %[m] length of soil column
pL = 10; %[kPa] top pressure
tau = L^2/c; %[s] charateristic consolidation time
targetTime=2*tau*[0.01,0.1,0.33,1];

%number of series terms
nm = 1000;


%array for time instants read from OutTime.dat
if plotMode == false
    fileName = 'OutTime.dat';
    t = readtime(fileName);
    zmesh = zvector;
    zmesh = zmesh';
    [p0fem,u0fem] = iniSol(zmesh,M,pL,Ku,biot,G);
    [pfem,ufem] = TerzaghiSol(u0fem,zmesh,t,nm,L,c,pL,biot,gamma,G,nu);
else
    t = [30 300 600 1800];
    nt = length(t);
    nz = 50; %number of calculation points along z-axis
    z = linspace(0,L,nz);   
    [p0,u0] = iniSol(z,M,pL,Ku,biot,G);    
    [p,u] = TerzaghiSol(u0,z,t,nm,L,c,pL,biot,gamma,G,nu);
    figure(1)
    plot(p,repmat(z',[1,nt]))
    figure(2)
    plot(u,repmat(z',[1,nt]))
end





%%
if plotMode == true
%convergence profile - temporary
h = repmat([0.5;0.333;0.25;0.2],1,3);
err_press=zeros(4,2);
err_press(:,1) = [0.0704, 0.0384, 0.0315,0.0301];
err_press(:,2)=[0.07,0.0322,0.0207, 0.0168]; 
err_press(:,3)=[0.07,0.0317,0.0181, 0.0123]; 
loglog(h(:,1),err_press(:,1),'-o','color','k')
hold on
loglog(h(:,2),err_press(:,2),'-x','color','k')
loglog(h(:,2),err_press(:,3),'-*','color','k')
xlim([0.15 0.6])
ylim([0.01 0.1])
grid on
xlabel('h (m) - Dimensione caratteristica dell` elemento ')
ylabel('Norma errore  L^2 della pressione')
legend('\Delta t=0.2 s','\Delta t =0.1 s','\Delta t =0.05 s')
L = legend;
L.AutoUpdate = 'off';

%triangle indicating slope of the diagram in log-log scale
%starting point of triangle
X = 0.42;
Y = 0.048;
dimX = 0.05;
slope = 2;
dimY = slope*dimX/X;
px = [X, X+dimX, X+dimX,X];
py = [Y, Y, Y+Y*dimY, Y];
plot(px,py,'color','k')
text(X+dimX/2,Y-0.25*Y*dimY,'1.0')
text(X+1.1*dimX,Y+0.5*Y*dimY,'2.0')
end
%slope

