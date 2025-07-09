function state = applyMandelIC(state,mat,msh,F)
% Get Material parameters from materials class
porosity = mat.getMaterial(1).PorousRock.getPorosity();
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

x = msh.coordinates(:,1);
z = msh.coordinates(:,3);

state.pressure = state.pressure+(1/(3*max(x)))*B*(1+nuU)*F;
state.dispConv(1:3:end) = -F*nuU*x/(2*mu*max(x));
state.dispConv(3:3:end) = F*(1-nuU)*z/(2*mu*max(x));
state.dispCurr = state.dispConv;
end

