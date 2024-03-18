function Terzaghi_analytical(mesh, mat, pL)

% number of terms for analyitcal soluton series
nm = 1000;

% Get model geometry
L_min = min(mesh.coordinates(:,3));
L_max = max(mesh.coordinates(:,3));
L = abs(L_max-L_min);

% Get Material parameters from materials class
k = mat.getMaterial(1).PorousRock.getPermVector();
k = k(1);
porosity = mat.getMaterial(1).PorousRock.getPorosity();
dyn = mat.getFluid().getDynViscosity;
cf = mat.getFluid().getFluidCompressibility; % [kPa^-1] Fluid compressibility
E = mat.getMaterial(1).ConstLaw.E;
nu = mat.getMaterial(1).ConstLaw.nu;
biot = mat.getMaterial(1).PorousRock.getBiotCoefficient();

% compute material parameters for analytical formulas
lambda =(E*nu)/((1+nu)*(1-2*nu)); %[kPa] first lamè constant
G = E/(2*(1+nu)); %[kPa] second lamè constant
M = (porosity*cf)^-1; %Biot Modulus, assuming cbr=0
Ku = lambda + 2*(G/3) + biot^2*M;
B = biot*M/Ku;
nu = lambda/(2*(lambda+G));
nuU = (3*nu+biot*B*(1-2*nu))/(3-biot*B*(1-2*nu)); %undrained Poisson Coefficient
c = 2*k*G*(1-nu)*(nuU-nu)/(dyn*biot^2*(1-nuU)*(1-2*nu)^2);
gamma = B*(1+nuU)/(3*(1-nuU));


t = readtime("outTime.dat");
nz = 50; %number of calculation points along z-axis
z = linspace(L_min,L_max,nz);
[~,u0] = iniSol(z,z,M,pL,Ku,biot,G);
[p,u] = TerzaghiSol(u0,z,z,t,nm,L,c,pL,biot,gamma,G,nu);
save('Terzaghi_Analytical.mat','p','u','z','t')
end

