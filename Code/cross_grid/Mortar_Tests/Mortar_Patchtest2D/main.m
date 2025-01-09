patch = 3;
switch patch
    case 1
        Fx = 0; %[kPa]
        Fy = -10;
        nu = 0;
    case {2,3}
        Fx = 10; %[kPa]
        Fy = 0;
        nu = 0.25;
end
% Elastic properties and constitutive tensor
E = 100000;
D = zeros(3);
D([1 5]) = 1-nu;
D([2 4]) = nu;
D(9) = 0.5*(1-2*nu);
D = (E/((1+nu)*(1-2*nu)))*D;
Dmat = D;

nel = 4;
nXs = nel+1;
rat = 1/2;
nYs = round(nXs);
nXm = round(nel*rat+1);
nYm = round(nXm);

stab = 'unstable';
%%
fprintf('Processing non conforming with dual multipliers \n')
[dual_tx,dual_ty,dual_ux,dual_uy,xNode] = RunNonConfPatchTest('dual',Fx,Fy,Dmat,patch,stab,nXs,nYs,nXm,nYm,nTip);
fprintf('Processing non conforming with standard multipliers \n')
[standard_tx,standard_ty,standard_ux,standard_uy,xNode] = RunNonConfPatchTest('standard',Fx,Fy,Dmat,patch,stab,nXs,nYs,nXm,nYm,nTip);
fprintf('Processing conforming \n')
[tx_anal,ty_anal,conf_ux,conf_uy,xNodeAnal] = RunConfPatch2(Fx,Fy,Dmat,patch,stab,nXs,nYs,nTip);

%% plotting
figure(1)
tiledlayout(1,4)
nexttile
plot(xNode,dual_tx,'k.-','LineWidth',1,'MarkerSize',10)
hold on
plot(xNode,standard_tx,'r.-','LineWidth',1,'MarkerSize',10)
plot(xNodeAnal,tx_anal,'b--','LineWidth',1)
xlabel('x (m)')
ylabel('t_x (kPa)')
legend('dual','standard','anal','Location','best')

nexttile
plot(xNode,dual_ty,'k.-','LineWidth',1,'MarkerSize',10)
hold on
plot(xNode,standard_ty,'r.-','LineWidth',1,'MarkerSize',10)
plot(xNodeAnal,ty_anal,'b--','LineWidth',1)
xlabel('x (m)')
ylabel('t_y (kPa)')
legend('dual','standard','anal','Location','best')

nexttile
plot(xNode,dual_ux,'k.','LineWidth',1,'MarkerSize',10)
hold on
plot(xNode,standard_ux,'r.','LineWidth',1,'MarkerSize',10)
plot(xNodeAnal,conf_ux,'b--','LineWidth',1)
xlabel('x (m)')
ylabel('u_x (kPa)')

nexttile
plot(xNode,dual_uy,'k.','LineWidth',1,'MarkerSize',10)
hold on
plot(xNode,standard_uy,'r.','LineWidth',1,'MarkerSize',10)
plot(xNodeAnal,conf_uy,'b--','LineWidth',1)
xlabel('x (m)')
ylabel('u_y (kPa)')
legend('dual','standard','anal','Location','best')

set(gcf, 'Position', [100, 100, 1000, 400])
%exportgraphics(gcf,strcat('Plots/Solution_Patch_',num2str(patch),'_',stab,'.png'),'Resolution',300)

% Plotting displacements on both interfaces (just for standard basis
% figure(2)
% tiledlayout(1,2)
% nexttile
% err_tx_dual = abs((dual_tx-tx_anal)./tx_anal);
% err_tx_stand = abs((standard_tx-tx_anal)./tx_anal);
% plot(xNode,err_tx_dual,'k.-','LineWidth',1,'MarkerSize',10)
% hold on
% plot(xNode,err_tx_stand,'r.-','LineWidth',1,'MarkerSize',10)
% xlabel('x (m)')
% ylabel('Rel. error t_x (kPa)')
% legend('dual','standard','Location','best')
% 
% nexttile
% err_ty_dual = abs((dual_ty-ty_anal)./ty_anal);
% err_ty_stand = abs((standard_ty-ty_anal)./ty_anal);
% plot(xNode,err_ty_dual,'k.-','LineWidth',1,'MarkerSize',10)
% hold on
% plot(xNode,err_ty_stand,'r.-','LineWidth',1,'MarkerSize',10)
% xlabel('x (m)')
% ylabel('Rel. error t_y (kPa)')
% legend('dual','standard','Location','best')
%exportgraphics(gcf,strcat('Plots/Rel_error_Patch_',num2str(patch),'_',stab,'.png'),'Resolution',300)
