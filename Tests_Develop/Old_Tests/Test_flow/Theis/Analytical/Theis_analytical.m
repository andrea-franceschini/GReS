close all
clear
%
% Theis analytical solution for fully penetrating wells in confined
% aquifers
%
% Generate Theis model input files
flOutParam = false;
%
% Compute Theis analytical solution and plot
% flTSol = false;
%
% Produce the BC files
flBCFile = true;
discrScheme = 'FVTPFA';  % of 'FEM'
%
% INPUT PARAMETERS
% Aquifer hydrological parameters
alpha = 1.e-5;      % Rock compressibility [1/kPa]
poro = 0.3;         % Porosity
beta = 4.6e-7;      % Water compressibility [1/kPa]
mu = 1.1574e-11;    % Water dynamic viscosity [kPa*d]
gamma = 9.81;       % Water specific weight [kPa/m]
k = 5.e-13;         % Rock permeability [m2]
%
b = 15;           % Aquifer thickness [m]
R = 500;          % Aquifer radius [m]
%
q = 50;           % Well rate [m3/d]
%
pIni = 392.4;     % Initial pressure [kPa]
%
% Derived properties
S = b*gamma*(alpha+poro*beta);  % Storativity
K = gamma*k/mu;                    % Hydraulic conductivity
T = b*K;                           % Transmissibility
%
if flOutParam
  % Generate MAT file with the above problem settings
  if isfile('TheisBenchSettings.mat')
    delete 'TheisBenchSettings.mat'
  end
  TModParam = struct("alpha",alpha,"poro",poro,"beta",beta,"mu",mu, ...
    "gamma",gamma,"k",k,"b",b,"q",q,"pIni",pIni,"S",S,"K",K,"T",T,"R",R);
  save("TheisBenchSettings.mat","TModParam")
end
%
%
% R = 500;     % Aquifer radial extension
% rMin = 0.5;
% rMult = 0.1;
% tMin = 0.1;
% tMax = 1000;
% tMult = 0.1;
rStep = [0.1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, ...
         24, 25, 26, 28, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 120, ...
         140, 160, 180, 200, 250, 300, 350, 400, 450, 500];
% tStep = [0.5, 1, 10, 20, 50, 100, 150, 200, 300, 400, 500, 1000, 2000, 3000];
tStep = [0.5, 1, 10, 20, 40, 60, 80, 100, 125, 150, 175, 200]; % 0 is included after to create the BC files
tStep = [0 tStep];
% rStep = [1 100 200];
% tStep = [1, 10];
%
[tVec,rVec] = meshgrid(tStep(2:end),rStep);
u = rVec(:).^2*S./(4*T*tVec(:));
W = expint(u);
dh = q/(4*pi*T)*W;
dh = reshape(dh,length(rStep),length(tStep)-1);
dh = [zeros(length(rStep),1), dh];
%
% Print the profiles
strTStep = string(tStep);
for i=1:length(tStep)
  plot(rStep,pIni-gamma*dh(:,i),'-o','DisplayName',strTStep(i));
  hold on
end
legend
hold off
%
% Print lateral boundary conditions
if flBCFile
  nItems = 336;
  fldPath = '../Hexa/dirSurfLatSurf_hexa/';
  for s = 1:length(tStep)
    fid = fopen(strcat(fldPath,'time',num2str(s-1),'.dat'),'w');
    fprintf(fid,'%s %f\n','% TIME',tStep(s));
    for j=1:nItems
      fprintf(fid,'%f\n',pIni-gamma*dh(end,s));
    end
    fclose(fid);
  end
  fName = 'dirSurfLatSurf_hexa';
  fid = fopen(strcat('../Hexa/',fName,'.dat'),'w');
  %
  if strcmp(discrScheme,'FEM')
    fprintf(fid,'%s\n','NodeBC                            % BC item');
  elseif strcmp(discrScheme,'FVTPFA')
    fprintf(fid,'%s\n','SurfBC                            % BC item');
  end
  fprintf(fid,'%s\n','Dir                               % BC type');
  fprintf(fid,'%s\n','Flow                              % Physics');
  fprintf(fid,'%s\n','latSurf                           % BC name');
  fprintf(fid,'%s\n',strcat(fName,'/list'));
  for s = 1:length(tStep)
    fprintf(fid,'%f %s\n',tStep(s),strcat(fName,'/time',num2str(s-1),'.dat'));
  end
  fprintf(fid,'%s\n','End');
end