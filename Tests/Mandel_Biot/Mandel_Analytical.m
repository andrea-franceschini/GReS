function Mandel_Analytical(mesh, mat, F)
% Compute Mandel analytical solution for specified time_steps
% Save results in MAT-FILE for postprocessing

fprintf('\n Computing Mandel Analytical solution... \n');

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

t = OutState.readTime("outTime.dat");
x = linspace(0,a,nx);
z = linspace(0,b,nz);
[p0,~,~] = mandelIni(B,F,x,z,nuU,mu);
[p,ux,uz,~] = mandelSol(alphan,x,t,c,p0,F,nu,nuU,mu,z);
save('Mandel_Analytical.mat','p','ux','uz','x','z','t')

fprintf('Done Computing Mandel Analytical solution. \n');
end


function [p,ux,uz,norm_p] = mandelSol(alphan,x,t,c,p0,F,nu,nuU,mu,z)
nx = length(x);
nt = length(t);
nz = length(z);


cellsP = arrayfun(@(alphan) bsxfun(@(x,t) (sin(alphan)/(alphan - ...
   sin(alphan)*cos(alphan)))*(cos(alphan*x/max(x))-cos(alphan))*exp(-alphan^2*c*t/max(x)^2),x',t),alphan,'UniformOutput',false);

cellsUx1 = arrayfun(@(alphan) arrayfun(@(t) ((sin(alphan)*cos(alphan))/(alphan - ...
   sin(alphan)*cos(alphan)))*exp(-alphan^2*c*t/max(x)^2),t),alphan,'UniformOutput',false);

cellsUx2 = arrayfun(@(alphan) bsxfun(@(x,t) ...
   ((cos(alphan)/(alphan - sin(alphan)*cos(alphan)))*sin(alphan*x/max(x))...
   *exp(-alphan^2*c*t/max(x)^2)),x',t),alphan,'UniformOutput',false);

cellsUz = arrayfun(@(alphan) arrayfun(@(t) ((sin(alphan)*cos(alphan))/(alphan - ...
   sin(alphan)*cos(alphan)))*exp(-alphan^2*c*t/max(x)^2),t),alphan,'UniformOutput',false);

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
ux = x'*((F*nu)/(2*mu*max(x)) - (F*nuU)/(mu*max(x))*seriesUx1)+ F/mu*seriesUx2;
uz = z'*(-(F*(1-nu)/(2*mu*max(x)))+(F*(1-nuU)/(mu*max(x)))*seriesUz);
norm_p = p/p0;
end

function  alphan = mandelZeros(f, rangeIn, nm)
rangeDim = 10;
alphan = [];
count = 0;
diff = 10;
toll = 0.0001;
while diff>toll && count<nm
   alphaTmp = findzeros(f, [rangeIn, rangeIn+rangeDim]);
   alphan = [alphan double(alphaTmp)];
   count = length(alphan);
   rangeIn = rangeIn+rangeDim;
   diff = abs(alphan(end)-alphan(end-1)-pi);
end

if count<nm
   alphan = [alphan(1:end-1) alphan(end):pi:alphan(end)+(nm-count)*pi];
end
end

function [p0,ux0,uz0] = mandelIni(B,F,x,z,nuU,mu)
%calculate initial solution for mandel benchamrk
p0 = (1/(3*max(x)))*B*(1+nuU)*F;
ux0 = F*nuU*x/(2*mu*max(x));
uz0 = -F*(1-nuU)*z/(2*mu*max(x));
end

function sol = findzeros(f,range,err)
if nargin < 3
  err = 1e-3;
end
sol = vpasolve(f,range);
if(isempty(sol))
  return
else
  lowLimit = sol-err;
  highLimit = sol+err;
  temp = findzeros(f,[range(1) lowLimit],1);
  if ~isempty(temp)
    sol = sort([sol temp]);
  end
  temp = findzeros(f,[highLimit range(2)],1);
  if ~isempty(temp)
    sol = sort([sol temp]);
  end
  return
end
end

