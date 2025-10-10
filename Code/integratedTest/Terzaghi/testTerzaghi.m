%% Coupled FEM - FV

model = ModelType(["SinglePhaseFlow_FVTPFA","Poromechanics_FEM"]);
runTerzaghi;
clearvars -except domain
analyticalSol = load('Terzaghi_Analytical.mat');
zAn = analyticalSol.z;
zNodes = domain.grid.topology.coordinates(:,3);
zCells = domain.grid.topology.cellCentroid(:,3);
pNum = domain.outstate.results.expPress(:,2:end);
uzNum = domain.outstate.results.expDispl(3:3:end,2:end);


for i = 1:numel(analyticalSol.t)
  % compare pressure and displacement solution at each time step
  % interpolate analytical pressure
  pAn = interp1(zAn,analyticalSol.p(:,i),zCells); 
  % interpolate analytical pressure
  uzAn = -interp1(zAn,analyticalSol.u(:,i),zNodes); 
  %
  assert(norm((pAn-pNum(:,i))./pNum(:,i))<1e0,'Pressure error for out time %i \n',i)
  relErrU = (uzAn-uzNum(:,i))./(uzNum(:,i));
  relErrU(isinf(relErrU)) = 0;
  assert(norm(relErrU)<1e0,'Displacement error for out time %i\n',i)
end

clear domain

%% Coupled FEM - FEM

model = ModelType(["SinglePhaseFlow_FEM","Poromechanics_FEM"]);
runTerzaghi;
clearvars -except domain

analyticalSol = load('Terzaghi_Analytical.mat');
zAn = analyticalSol.z;
zNodes = domain.grid.topology.coordinates(:,3);
pNum = domain.outstate.results.expPress(:,2:end);
uzNum = domain.outstate.results.expDispl(3:3:end,2:end);


for i = 1:numel(analyticalSol.t)
  % compare pressure and displacement solution at each time step
  % interpolate analytical pressure
  pAn = interp1(zAn,analyticalSol.p(:,i),zNodes);
  % interpolate analytical pressure
  uzAn = -interp1(zAn,analyticalSol.u(:,i),zNodes);
  %
  relErrP = (pAn-pNum(:,i))./pNum(:,i);
  relErrP(isinf(relErrP)) = 0;
  relErrP(isnan(relErrP)) = 0;
  assert(norm(relErrP)<1e0,'Pressure error for out time %i \n',i)
  relErrU = (uzAn-uzNum(:,i))./(uzNum(:,i));
  relErrU(isinf(relErrU)) = 0;
  assert(norm(relErrU)<1e0,'Displacement error for out time %i\n',i)
end
