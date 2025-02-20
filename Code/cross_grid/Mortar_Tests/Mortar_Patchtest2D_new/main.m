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

E1 = 100000;
D = zeros(3);
D([1 5]) = 1-nu;
D([2 4]) = nu;
D(9) = 0.5*(1-2*nu);
D = (E1/((1+nu)*(1-2*nu)))*D;
DmatM = D;

E2 = 100000;
D = zeros(3);
D([1 5]) = 1-nu;
D([2 4]) = nu;
D(9) = 0.5*(1-2*nu);
D = (E2/((1+nu)*(1-2*nu)))*D;
DmatS = D;

nel = 10;
nXs = nel+1;
rat = 0.4;
nYs = round(nXs);
nXm = round(nel*rat+1);
nYm = round(nXm);
h = 1/nel;
alpha = 100;
gamma = alpha*h/E1;


stab = 'unstable';
nTip = 0;
%%
fprintf('Processing non conforming with dual multipliers \n')
[dual_tx,dual_ty,dual_us,dual_um] = RunNonConfPatchTest('standard',Fx,Fy,DmatS,DmatM,patch,stab,nXs,nYs,nXm,nYm);
fprintf('Processing non conforming with standard multipliers \n')
[standard_tx,standard_ty,standard_us,standard_um,x] = RunNonConfConstant(Fx,Fy,DmatS,DmatM,patch,nXs,nYs,nXm,nYm,E1,E2);
fprintf('Processing conforming \n')
[tx_conf,ty_conf,conf_us,conf_um] = RunConfPatch(Fx,Fy,DmatS,DmatM,patch,'unstable',nXs,nYs,nYm);

xNodeM = linspace(0,1,size(dual_um,1));
xNodeS = linspace(0,1,size(dual_us,1));
xNodeConf = linspace(0,1,size(conf_us,1));



%% plotting
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultAxesFontName','Times')

ms = 12;
figure(1)
tiledlayout(1,2)
nexttile
plot(xNodeS,dual_tx,'b.-','LineWidth',1,'MarkerSize',ms)
hold on
plot(x,standard_tx,'r.-','LineWidth',1,'MarkerSize',ms)
plot(xNodeConf,tx_conf,'k-','LineWidth',1)
xlabel('x (m)','Interpreter','latex')
ylabel('$t_x$ (kPa)','Interpreter','latex')
legend('Dual','Constant','Conforming','Location','best','Interpreter','latex')
if patch==1
   ylim([-0.1 0.1])
end

nexttile
plot(xNodeS,dual_ty,'b.-','LineWidth',1,'MarkerSize',ms)
hold on
plot(x,standard_ty,'r.-','LineWidth',1,'MarkerSize',ms)
plot(xNodeConf,ty_conf,'k-','LineWidth',1)
xlabel('x (m)','Interpreter','latex')
ylabel('$t_y$ (kPa)','Interpreter','latex')
legend('Dual','Constant','Conforming','Location','best','Interpreter','latex')
if patch==1
   ylim([-11 -9])
end

% comparing displacements
% nexttile
% plot(xNodeS,dual_us(:,2),'ko','LineWidth',1,'MarkerSize',ms)
% hold on
% plot(xNodeM,dual_um(:,2),'ks','LineWidth',1,'MarkerSize',ms,'Ma0rkerFaceColor','k')
% plot(xNodeConf,conf_us(:,2),'k-','LineWidth',1,'MarkerSize',ms)
% xlabel('x (m)')
% ylabel('u_y (kPa)')
% legend('dual_slave','dual_master','conforming')
% title('dual u_y')
% 
% nexttile
% plot(xNodeS,dual_us(:,1),'ko','LineWidth',1,'MarkerSize',ms)
% hold on
% plot(xNodeM,dual_um(:,1),'ks','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','k')
% plot(xNodeConf,conf_us(:,1),'k-','LineWidth',1,'MarkerSize',ms)
% xlabel('x (m)')
% ylabel('u_x (kPa)')
% legend('standard slave','standard master','conforming')
% title('dual u_x')




% nexttile
% plot(xNode,dual_ux,'k.','LineWidth',1,'MarkerSize',10)
% hold on
% plot(xNode,standard_ux,'r.','LineWidth',1,'MarkerSize',10)
% plot(xNodeAnal,conf_ux,'b--','LineWidth',1)
% xlabel('x (m)')
% ylabel('u_x (kPa)')

set(gcf, 'Position', [100, 100, 600, 400])
%exportgraphics(gcf,strcat('Results/Solution_Patch_',num2str(patch),'.png'),'Resolution',400)


%exportgraphics(gcf,strcat('Plots/Rel_error_Patch_',num2str(patch),'_',stab,'.png'),'Resolution',300)
