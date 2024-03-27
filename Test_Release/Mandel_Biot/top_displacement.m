% MANDEL PROBLEM utility to generate input files for top displacement
%Solution according to Verrujit: Theory and Problems of Poroelasticity

%% COMPUTE SOLUTION
% flag to save initial solution for FEM analysis
plotMode = false;

%Proprietà materiale (sabbia)
K = 1.0e-12; %[m^2] Intrinsic permability
porosity = 0.375;
dyn = 1.00e-6; %[kPa*s]
cf = 4.4e-7 ; % [kPa^-1] fluid compressibility
lambda = 40000; %[kPa] costante lamè
mu = 40000; %[kPa] lamè constant
alpha = 1.0; %coefficiente di Biot
nu = lambda/(2*(lambda+mu));

%domain dimensions
F = 10; %[kN/m]
nx = 50; %numb. of calc.points along axis x
nz = 50; %numb. of calc.points along axis z (vertical)
a = 1; %[m]
b = 1; %[m]
x = linspace(0,a,nx);
z = linspace(0,a,nz);

M = (porosity*cf)^-1; %Biot Modulus, assuming cbr=0
Ku = lambda + 2*(mu/3) + alpha^2*M;
cm = (lambda+2*mu)^-1; %vertical uniaxial compressibility

B = alpha*M/Ku;


nuU = (3*nu+alpha*B*(1-2*nu))/(3-alpha*B*(1-2*nu));
fac = (1-nu)/(nuU-nu);


%version of consolidation coefficient with hydraulic conductivity
%c = k/(9.81*(M^-1+alpha^2*cm));

%version of consolidation coefficient with intrinsic permability
c = (2*K*mu*(1-nu)*(nuU-nu))/(dyn*alpha^2*(1-nuU)*(1-2*nu)^2);

%calculating zeros of mandel's function tan(alpha)+fac*alpha;
syms w
f = (tan(w)-fac*w);

%number of serie terms
nm = 20;
rangeIn = 0.1;
%calculate a specified number of series terms by solving a non linear
%equation
alphan = mandelZeros(f, rangeIn, nm);

%c = 2*k*G*(1-nu)*(nuU-nu)/(dyn*biot^2*(1-nuU)*(1-2*nu)^2); %consolidation coefficient
% list of required time steps for computazion of uz;
t = load('list_top_DZ.dat');
% list of nodes belonging to top surface
nList = unique(topology.surfaces(topology.surfaceTag == 2,:));
% get coordinates
zmesh = (topology.coordinates(nList,3))';
xmesh = (topology.coordinates(nList,1))';
[p0fem,ux0fem,uz0fem] = mandelIni(B,F,xmesh,zmesh,nuU,mu);
[pfem,uxfem,uzfem,norm_p] = mandelSol(alphan,xmesh,t(2:end)',c,p0fem,F,nu,nuU,mu,zmesh);
uzfem = [uz0fem' uzfem];
% writing file with nodes list
%% PRINT FILES

%fList = fopen('topDisp\list', 'w');
% fprintf(fList, '%i \n', nList);

% writing list of files for general BCS input file
% fid = fopen('BCList', 'w');
% print list of files
% for i = 1:length(t)
%     str = strcat('topDisp/time',num2str(i-1),'.dat');
%     fprintf(fid, '%2.5f  %s \n', t(i), str);
% end


for i = 1:length(t)
    str = strcat('topDisp\time',num2str(i-1),'.dat');
    f = fopen(str, 'w');
    fprintf(f, '% TIME %2.4f \n', t(i));
    fprintf(f, '%2.6e \n', uzfem(:,i));
end 

