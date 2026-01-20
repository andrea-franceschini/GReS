function applyChemoMechanicsIC(state,mat,mesh)
% Apply initial conditions for Terzaghi problem

% Get Material parameters from materials class
porosity = mat.getMaterial(1).PorousRock.getPorosity();
cf = mat.getFluid().getFluidCompressibility; % [kPa^-1] Fluid compressibility
E = mat.getMaterial(1).ConstLaw.E;
nu = mat.getMaterial(1).ConstLaw.nu;
biot = mat.getMaterial(1).PorousRock.getBiotCoefficient();

% compute material parameters for analytical formulas
lambda =(E*nu)/((1+nu)*(1-2*nu)); %[kPa] first lamè constant
G = E/(2*(1+nu)); %[kPa] second lamè constant
M = (porosity*cf)^-1; %Biot Modulus, assuming cbr=0
Ku = lambda + 2*(G/3) + biot^2*M;

c_in = 6195; %  initial concentration value

state.data.pressure = state.data.pressure + c_in;
% zu = mesh.coordinates(:,3);
% state.data.dispConv(3:3:end) = arrayfun(@(zu) 1/(Ku+4*G/3)*pL*(zu),zu);
% state.data.dispCurr = state.data.dispConv;
end