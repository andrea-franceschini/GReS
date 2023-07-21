clear; clc;
%Solution according to Verrujit: Theory and Problems of Poroelasticity


%Proprietà materiale (sabbia)
K = 1.00e-12; %[m^2] Intrinsic permability
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
x = linspace(0,a,nx);

M = (porosity*cf)^-1; %Biot Modulus, assuming cbr=0
Ku = lambda + 2*(mu/3) + alpha^2*M;
cm = (lambda+2*mu)^-1; %vertical uniaxial compressibility

B = alpha*M/Ku;

nuU = (3*nu+alpha*B*(1-2*nu))/(3-alpha*B*(1-2*nu));
fac = 2*(1-nu)/(1-2*nu);


%version of consolidation coefficient with hydraulic conductivity
%c = k/(9.81*(M^-1+alpha^2*cm));

%version of consolidation coefficient with intrinsic permability
c = (2*K*mu*(1-nu)*(nuU-nu))/(dyn*alpha^2*(1-nuU)*(1-2*nu)^2);

%calculating zeros of mandel's function tan(alpha)+fac*alpha;
syms w
f = (tan(w)-fac*w);

%ftry = @(q) (tan(q)+fac*q);
%number of serie terms
nm = 100;
rangeIn = 0.1;
alphan = mandelZeros(f, rangeIn, nm); 


%c = 2*k*G*(1-nu)*(nuU-nu)/(dyn*biot^2*(1-nuU)*(1-2*nu)^2); %consolidation coefficient

p0 = (1/(3*a))*B*(1+nuU)*F;
u0 = F*nuU*x/(2*mu*a);


%array for time instants
t = [0.05 0.5 5 10];
nt = length(t);


cellsP = arrayfun(@(alphan) bsxfun(@(x,t) (sin(alphan)/(alphan - ...
    sin(alphan)*cos(alphan)))*(cos(alphan*x/a)-cos(alphan))*exp(-alphan^2*c*t/a^2),x',t),alphan,'UniformOutput',false);

% cellsU1 = arrayfun(@(alphan) bsxfun(@(x,t) ((sin(alphan)*cos(alphan))/(alphan - ...
%     sin(alphan)*cos(alphan)))*exp(-alphan^2*c*t/a^2),x',t),alphan','UniformOutput',false);
% 
% cellsU2 = arrayfun(@(alphan) bsxfun(@(x,t) (sin(alphan)/(alphan - ...
%     sin(alphan)*cos(alphan)))*cos(alphan*x/a)*exp(-alphan^2*c*t/a^2),x',t),alphan','UniformOutput',false);


matP = cell2mat(cellsP);
% matU1 = cell2mat(cellsU1);
% matU2 = cell2mat(cellsU2);

seriesP = sum(reshape(matP,nx,nt,[]),3);
% seriesU1 = sum(reshape(matU1,nx,nt,[]),3);
% seriesU2 = sum(reshape(matU2,nx,nt,[]),3);


p = 2*p0*seriesP;
norm_p = p/p0;

figure(1)
plot(x/a, norm_p ,'-o')



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

figure(2)
semilogx(dimless_time/time_factor, norm_ptime' ,'-')
legend('x=0','x=a/4','x=a/2')


% f = @(x) (-tan(x)-fac*x);
% fder = @(x) (1/(cos(x))^2+fac);
% iter = 0;
% tol = 1.0e-5;
% for i = 1:nm
%     iter = 0;
%     x0 = 2*i - 2;
%     dx = 1;
% while dx>tol && iter<100
%     dx = -f(x0)/fder(x0);
%     alphan(i) = x0 + dx;
%     x0 =  x0 + dx;
%     iter = iter +1;
% end
% end

% x0 = 0;
% for i=1:nm
%     alphan(i) = fzero(@(a)  tan(a)+fac*a,x0);
%     x0 = x0+2;
% end




