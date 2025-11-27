function Terzaghi_analytical(mesh, mat, pL, time)
fprintf('Computing Terzaghi Analytical solution... \n');
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

nz = 50; %number of calculation points along z-axis
z = linspace(L_min,L_max,nz);
[~,u0] = iniSol(z,z,M,pL,Ku,biot,G);
[p,u] = TerzaghiSol(u0,z,z,time,nm,L,c,pL,biot,gamma,G,nu);
save('Terzaghi_Analytical.mat','p','u','z','time')
fprintf('Done computing Terzaghi analytical solution.\n');
end

function [p,u] = TerzaghiSol(u0,zu,zp,t,nm,L,c,pL,biot,gamma,G,nu)

nzu = length(zu);
nzp = length(zp);
m = 0:nm;
nt = length(t);

cellsP = arrayfun(@(m) bsxfun(@(zp,t) (1/(2*m+1))*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*sin(((2*m+1)*pi*(L-zp))/(2*L)),zp',t),m,'UniformOutput',false);
cellsU = arrayfun(@(m) bsxfun(@(zu,t) (1/(2*m+1)^2)*exp((-(2*m+1)^2*pi^2*c*t)/(4*L^2))*cos(((2*m+1)*pi*(L-zu))/(2*L)),zu',t),m,'UniformOutput',false);

matP = cell2mat(cellsP);
matU = cell2mat(cellsU);

seriesP = sum(reshape(matP,nzp,nt,[]),3);
seriesU = sum(reshape(matU, nzu,nt,[]),3);

p = (4*gamma*pL/pi)*seriesP; %[kPa]
u = repmat(u0',1,nt) + ((1-2*nu)*biot*gamma*pL/2/G/(1-nu))*(repmat(zu',1,nt) - (8*L/pi^2)*seriesU); 

end

function [p0,u0] = iniSol(zu,zp,M,pL,Ku,biot,G)
nzp = length(zp);
p0 = zeros(nzp,1);
p0(1:end) = (biot*M*pL)/(Ku+4*G/3);
u0 = arrayfun(@(zu) 1/(Ku+4*G/3)*pL*(zu),zu);
end



