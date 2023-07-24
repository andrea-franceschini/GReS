clear; clc;
%Solution according to Verrujit: Theory and Problems of Poroelasticity

% flag to save initial solution for FEM analysis 
writemode = true;

%Proprietà materiale (sabbia)
K = 1.0e-12; %[m^2] Intrinsic permability
porosity = 0.375; 
dyn = 1.00e-6; %[kPa*s]
cf = 4.4e-7 ; % [kPa^-1] fluid compressibility
lambda = 40000; %[kPa] costante lamè
mu = 40000; %[kPa] lamè constant
alpha = 1.0; %coefficiente di Biot
nu = lambda/(2*(lambda+mu));

a = 1; %[m]
b = 1; %[m]
F = 10; %[kN/m]
nx = 50; %numb. of calc.points along axis x
nz = 50; %numb. of calc.points along axis z (vertical)
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

%ftry = @(q) (tan(q)+fac*q);
%number of serie terms
nm = 20;
rangeIn = 0.1;
alphan = mandelZeros(f, rangeIn, nm); 


%c = 2*k*G*(1-nu)*(nuU-nu)/(dyn*biot^2*(1-nuU)*(1-2*nu)^2); %consolidation coefficient

%calculating initial solution
p0 = (1/(3*a))*B*(1+nuU)*F;
ux0 = F*nuU*x/(2*mu*a);
uz0 = -F*(1-nuU)*z/(2*mu*a);

%array for time instants
t = [0.05 0.5 5 10];
nt = length(t);

%calculating series terms (creates cell array)
cellsP = arrayfun(@(alphan) bsxfun(@(x,t) (sin(alphan)/(alphan - ...
    sin(alphan)*cos(alphan)))*(cos(alphan*x/a)-cos(alphan))*exp(-alphan^2*c*t/a^2),x',t),alphan,'UniformOutput',false);

cellsUx1 = arrayfun(@(alphan) arrayfun(@(t) ((sin(alphan)*cos(alphan))/(alphan - ...
   sin(alphan)*cos(alphan)))*exp(-alphan^2*c*t/a^2),t),alphan,'UniformOutput',false);

cellsUx2 = arrayfun(@(alphan) bsxfun(@(x,t) ...
    ((cos(alphan)/(alphan - sin(alphan)*cos(alphan)))*sin(alphan*x/a)...
    *exp(-alphan^2*c*t/a^2)),x',t),alphan,'UniformOutput',false);

cellsUz = arrayfun(@(alphan) arrayfun(@(t) ((sin(alphan)*cos(alphan))/(alphan - ...
   sin(alphan)*cos(alphan)))*exp(-alphan^2*c*t/a^2),t),alphan,'UniformOutput',false);

%converting cell array into 2S matrices
matP = cell2mat(cellsP);
matUx1 = cell2mat(cellsUx1);
matUx2 = cell2mat(cellsUx2);
matUz = cell2mat(cellsUz);
%reshaping matrices and summation along third direction. 
%Dim 1: position
%Dim 2: time
%Dim 3: Series term

seriesP = sum(reshape(matP,nx,nt,[]),3);
seriesUx1 = sum(reshape(matUx1,1,nt,[]),3);
seriesUx2 = sum(reshape(matUx2,nz,nt,[]),3);
seriesUz = sum(reshape(matUz,1,nt,[]),3);

%calculating solutions
p = 2*p0*seriesP;
ux = x'*((F*nu)/(2*mu*a) - (F*nuU)/(mu*a)*seriesUx1)+ F/mu*seriesUx2;
uz = z'*(-(F*(1-nu)/(2*mu*a))+(F*(1-nuU)/(mu*a))*seriesUz);
norm_p = p/p0;


time_string = "t =" + t + "s";


figure(1)
plot(x/a, norm_p ,'-o')
legend(time_string)
xlabel('Normalized distance x/a')
ylabel('Normalized pressure p/p_0');

figure(2)
plot(x,ux,'-o')
legend(time_string)
xlabel('Distance (m)')
ylabel('Displacement DX (m)')

figure(3)
plot(-uz,z,'-o')
legend(time_string)
ylabel('Distance (m)')
xlabel('Displacement DZ (m)')


%%%%%%%%pressure VS time
time_factor = a^2/c;
x_time = [0, a/4, a/2];
nx_dim = length(x_time);
dimless_time = time_factor*linspace(0.001,1,1000);
nt_dim = length(dimless_time);

cellsPtime = arrayfun(@(alphan) bsxfun(@(x_time,dimless_time) (sin(alphan)/(alphan - ...
    sin(alphan)*cos(alphan)))*(cos(alphan*x_time/a)-cos(alphan))*exp(-alphan^2*c*dimless_time/a^2),x_time',dimless_time),alphan,'UniformOutput',false);

% cellsU1 = arrayfun(@(alphan) bsxfun(@(x,t) ((sin(alphan)*cos(alphan))/(alphan - ...
%     sin(alphan)*cos(alphan)))*exp(-alphan^2*c*t/a^2),x',t),alphan','UniformOutput',false);
% 
% cellsU2 = arrayfun(@(alphan) bsxfun(@(x,t) (sin(alphan)/(alphan - ...
%     sin(alphan)*cos(alphan)))*cos(alphan*x/a)*exp(-alphan^2*c*t/a^2),x',t),alphan','UniformOutput',false);


matPtime = cell2mat(cellsPtime);
% matU1 = cell2mat(cellsU1);
% matU2 = cell2mat(cellsU2);

seriesPtime = sum(reshape(matPtime,nx_dim,nt_dim,[]),3);
% seriesU1 = sum(reshape(matU1,nx,nt,[]),3);
% seriesU2 = sum(reshape(matU2,nx,nt,[]),3);


ptime = 2*p0*seriesPtime;
norm_ptime = ptime/p0;

figure(4)
semilogx(dimless_time/time_factor, norm_ptime' ,'-')
legend('x=0','x=a/4','x=a/2')


%%%%calculating initial solution to be assign to the model%%%%%%%%%%%%%%%%
%First save in the repo the vectors of coordinates for each node
if writemode
    zmesh = load('zmesh.dat');
    xmesh = load('xmesh.dat');
    nNodes = length(zmesh);
    p0fem = zeros(nNodes,1);
    ux0fem = zeros(nNodes,1);
    uz0fem = zeros(nNodes,1);
    p0fem(:) = (1/(3*a))*B*(1+nuU)*F;
    ux0fem(:) = F*nuU*xmesh/(2*mu*a);
    uz0fem(:) = -F*(1-nuU)*zmesh/(2*mu*a);

    save C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Hexa\H01\pIniH01.dat p0fem -ascii
    save C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Hexa\H01\uxIniH01.dat ux0fem -ascii
    save C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Hexa\H01\uzIniH01.dat uz0fem -ascii
    save C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Hexa\H01\pAnal.dat p -ascii
    save C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Hexa\H01\uxAnal.dat ux -ascii
    save C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Hexa\H01\uzAnal.dat uz -ascii
    save C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Hexa\H01\xAnal.dat x -ascii
    save C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Hexa\H01\zAnal.dat z -ascii
    
else 
    disp('WriteMode currently set to false')
end
