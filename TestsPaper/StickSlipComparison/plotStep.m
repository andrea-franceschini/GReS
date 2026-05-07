function plotStep(historyFile,solv,tStep)

% this function produces plot for the StickSlipCase for a requested time
% step

% to plot: sigma_n and norm tangential traction along the crack vertical
% axis
vars = load(historyFile);
mult = [vars.output.multipliers];

mult = mult(:,tStep);
sn = mult(1:3:end);
norm_tT = sqrt(mult(2:3:end).^2+mult(3:3:end).^2);

interf = solv.interfaces{1};

mshSlave = getMesh(interf,MortarSide.slave);

surfC = mshSlave.surfaceCentroid;

NX = round(10/(2*surfC(1,2)));
NY = round(15/(2*surfC(1,3)));

if interf.multiplierLocation == entityField.surface

  vertAxis1 = linspace(0.5*NX,NX*NY - 0.5*NX,NX);
  vertAxis2 = vertAxis1+1;

  snPlot = 0.5*(sn(vertAxis1)+sn(vertAxis2));
  tTPlot = 0.5*(norm_tT(vertAxis1)+norm_tT(vertAxis2));

  zPlot = surfC(vertAxis2,3);

else

  [~,i] = sort(mshSlave.coordinates(:,3));
  coord = mshSlave.coordinates(i,:);
  nodes = coord(:,2) == 5;
  sn = sn(i);
  tT = norm_tT(i);
  snPlot = sn(nodes);
  tTPlot = tT(nodes);
  zPlot = linspace(0,15,NY+1);


end

figure
plot(snPlot,zPlot,'r-s')
xlabel('sigma_n')
ylabel('z(m)')

figure
plot(tTPlot,zPlot,'b-s')
xlabel('norm tangential traction')
ylabel('z(m)')



end