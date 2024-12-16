Fx = 0; %[kPa]
Fy = -10;
% Elastic properties and constitutive tensor
E = 100000;
nu = 0.25;
D = zeros(3);
D([1 5]) = 1-nu;
D([2 4]) = nu;
D(9) = 0.5*(1-2*nu);
D = (E/((1+nu)*(1-2*nu)))*D;
Dmat = D;
%%
fprintf('Processing non conforming with dual multipliers \n')
[dual_tx,dual_ty,dual_ux] = RunNonConfPatchTest('dual',Fx,Fy,Dmat);
fprintf('Processing non conforming with standard multipliers \n')
[standard_tx,standard_ty,standard_ux] = RunNonConfPatchTest('standard',Fx,Fy,Dmat);
fprintf('Processing conforming \n')
[tx_anal,ty_anal,conf_ux] = RunConfPatchTest(Fx,Fy,Dmat);

%% plotting
xNode = linspace(0,1,numel(dual_tx));
xCell = linspace(0+0.5/numel(tx_anal),1-0.5/numel(tx_anal),numel(tx_anal));
xNodeAnal = linspace(0,1,numel(conf_ux));
tiledlayout(1,3)
nexttile
plot(xNode,dual_tx,'k.-','LineWidth',1,'MarkerSize',10)
hold on
plot(xNode,standard_tx,'r.-','LineWidth',1,'MarkerSize',10)
plot(xCell,tx_anal,'b--','LineWidth',1)
xlabel('x (m)')
ylabel('t_x (kPa)')
legend('dual','standard','anal')

nexttile
plot(xNode,dual_ty,'k.-','LineWidth',1,'MarkerSize',10)
hold on
plot(xNode,standard_ty,'r.-','LineWidth',1,'MarkerSize',10)
plot(xCell,ty_anal,'b--','LineWidth',1)
xlabel('x (m)')
ylabel('t_y (kPa)')
legend('dual','standard','anal')

nexttile
plot(xNode,dual_ux,'k.','LineWidth',1,'MarkerSize',10)
hold on
plot(xNode,standard_ux,'r.','LineWidth',1,'MarkerSize',10)
plot(xNodeAnal,conf_ux,'b--','LineWidth',1)
xlabel('x (m)')
ylabel('t_y (kPa)')
legend('dual','standard','anal')