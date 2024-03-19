function Mandel_Analytical(mesh, mat, F)
% Compute Mandel analytical solution for specified time_steps
% Save results in MAT-FILE for postprocessing

fprintf('\n Computing Analytical solution \n');

% Get model geometry
nx = 50; %numb. of calc.points along axis x
nz = 50; %numb. of calc.points along axis z (vertical)
x_min = min(mesh.coordinates(:,1));
x_max = max(mesh.coordinates(:,1));
a = abs(x_max-x_min);
y_min = min(mesh.coordinates(:,3));
y_max = max(mesh.coordinates(:,3));
b = abs(y_max-y_min);


% Get Material parameters from materials class
K = mat.getMaterial(1).PorousRock.getPermVector();
K = K(1);
porosity = mat.getMaterial(1).PorousRock.getPorosity();
dyn = mat.getFluid().getDynViscosity;
cf = mat.getFluid().getFluidCompressibility; % [kPa^-1] Fluid compressibility
E = mat.getMaterial(1).ConstLaw.E;
nu = mat.getMaterial(1).ConstLaw.nu;
alpha = mat.getMaterial(1).PorousRock.getBiotCoefficient();

% Compute derivated material parameters
lambda =(E*nu)/((1+nu)*(1-2*nu)); %[kPa] first lamè constant
mu = E/(2*(1+nu)); %[kPa] second lamè constant
M = (porosity*cf)^-1; %Biot Modulus, assuming cbr=0
Ku = lambda + 2*(mu/3) + alpha^2*M;
B = alpha*M/Ku;
nuU = (3*nu+alpha*B*(1-2*nu))/(3-alpha*B*(1-2*nu));
fac = (1-nu)/(nuU-nu);


%Consolidation coefficient with intrinsic permability
c = (2*K*mu*(1-nu)*(nuU-nu))/(dyn*alpha^2*(1-nuU)*(1-2*nu)^2);


% Compute zeros of mandel's function tan(alpha)+fac*alpha;
syms w
f = (tan(w)-fac*w);
%number of serie terms
nm = 20;
rangeIn = 0.1;
%calculate a specified number of series terms by solving a non linear equation
alphan = mandelZeros(f, rangeIn, nm);

t = readtime("outTime.dat");
x = linspace(0,a,nx);
z = linspace(0,b,nz);
[p0,~,~] = mandelIni(B,F,x,z,nuU,mu);
[p,ux,uz,~] = mandelSol(alphan,x,t,c,p0,F,nu,nuU,mu,z);
save('Mandel_Analytical.mat','p','ux','uz','x','z','t')




end