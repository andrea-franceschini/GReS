clear; clc;
%Script for Terzaghi analytical solution, uses intrinsic permeability.
%Found in JCP_2018.



%Fluid and poroelastic properties
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
c_alternative = khydro/ (9.81*(M^-1+biot^2*cm));
%z vector
L = 10; %[m] lunghezza colonna di terreno
nz = 21; %number of calculation points along z-axis
z = linspace(0,L,nz);


pL = 10; %[kPa] top pressure
tau = L^2/c; %[s] charateristic consolidation time

p0 = zeros(nz,1);
p0(1:end) = (biot*M*pL)/(Ku+4*G/3);
u0 = arrayfun(@(z) 1/(Ku+4*G/3)*pL*(z),z);


%number of series terms
nm = 1000;
m = [0:nm];


%array for time instants
t = [30 300 600 1800];
nt = length(t);

cellsP = arrayfun(@(m) bsxfun(@(z,t) (1/(2*m+1))*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*sin(((2*m+1)*pi*(L-z))/(2*L)),z',t),m,'UniformOutput',false);
cellsU = arrayfun(@(m) bsxfun(@(z,t) (1/(2*m+1)^2)*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*cos(((2*m+1)*pi*(L-z))/(2*L)),z',t),m,'UniformOutput',false);

matP = cell2mat(cellsP);
matU = cell2mat(cellsU);

seriesP = sum(reshape(matP,nz,nt,[]),3);
seriesU = sum(reshape(matU, nz,nt,[]),3);

p = (4*gamma*pL/pi)*seriesP; %[kPa]
u = repmat(u0',1,nt) + ((1-2*nu)*biot*gamma*pL/2/G/(1-nu))*(repmat(z',1,nt) - (8*L/pi^2)*seriesU); 

figure(1)
plot(p,repmat(z',[1,nt]),'-o')


figure(2)
plot(u,repmat(z',[1,nt]),'-o')


%calculating analytical solution on all nodes of the mesh
zmesh = load('depht_H02.dat');
zmesh = zmesh';
nzm = length(zmesh);

p0fem = zeros(nzm,1);
p0fem(1:end) = (biot*M*pL)/(Ku+4*G/3);
u0fem = arrayfun(@(zmesh) 1/(Ku+4*G/3)*pL*(zmesh),zmesh);
u0femt = u0fem';

%save pIni02.dat p0fem -ascii
%save uIni02.dat u0femt -ascii


cellsPfem = arrayfun(@(m) bsxfun(@(zmesh,t) (1/(2*m+1))*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*sin(((2*m+1)*pi*(L-zmesh))/(2*L)),zmesh',t),m,'UniformOutput',false);
cellsUfem = arrayfun(@(m) bsxfun(@(zmesh,t) (1/(2*m+1)^2)*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*cos(((2*m+1)*pi*(L-zmesh))/(2*L)),zmesh',t),m,'UniformOutput',false);

matPfem = cell2mat(cellsPfem);
matUfem = cell2mat(cellsUfem);

seriesPfem = sum(reshape(matPfem,nzm,nt,[]),3);
seriesUfem = sum(reshape(matUfem, nzm,nt,[]),3);

pfem = (4*gamma*pL/pi)*seriesPfem; %[kPa]
ufem = repmat(u0fem',1,nt) + ((1-2*nu)*biot*gamma*pL/2/G/(1-nu))*(repmat(zmesh',1,nt) - (8*L/pi^2)*seriesUfem); 

%saving p for future 
save analytical_pressure_H02.dat pfem -ascii
save analytical_displacement_H02.dat ufem -ascii


position = repmat(z',[1,nt]);
save zcoord.dat position -ascii