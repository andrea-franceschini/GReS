%% testing channel size computation (naive but robust technique)
% define 2 cartesian grids
msh1 = Mesh();
msh1.createCartesianGrid(2,1,[0 1],[0 1],12,12);

msh2 = Mesh();
msh2.createCartesianGrid(2,1,[0 1],[0 1],15,15);

% modify cartesian grid orientation
fac = 0.3;
msh2.coordinates(:,3) = fac*msh2.coordinates(:,1) + fac*msh2.coordinates(:,2) + 0.2; 
p = [0.2,0.3,0.1];
d = computeChannelSize(p,msh1,msh2);

[x1,y1,z1] = deal(msh1.coordinates(:,1),msh1.coordinates(:,2),msh1.coordinates(:,3));
[x2,y2,z2] = deal(msh2.coordinates(:,1),msh2.coordinates(:,2),msh2.coordinates(:,3));
scatter3(x1,y1,z1)
hold on
scatter3(x2,y2,z2)
plot3(p(1),p(2),p(3),'ro')

%% roughly predicting change in permeability given a grain compressibility and a pressure

kg = 1e-8; %[kPa-1]
p = 5e4;   % [kPa]
c = 5e-6; % initial channel size (micron)
dg = 500e-6; %[diameter of sandy grain in micron]
V = 4/3*pi*(dg/2)^3;
r = dg/2;
K = (c^3/12); % initial permeability


dV = kg*p;
Vnew = V*(1+dV); % volume change [micron^3]
rnew = (3*Vnew/(pi*4))^(1/3);
dc = 2*(rnew-r);

dK = ((c+dc)^3-c^3)/(c^3); % percentage change of permeability
delta_K = dK*K;
