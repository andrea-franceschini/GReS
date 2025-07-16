%
clear 
close all

ms = 6;


out_P0 = load('Traction_results/outP0_v1.mat');
out_dual = load('Traction_results/outDual.mat');
out_conforming = load('Traction_results/outConforming.mat');

tP0 = 14;
tConf = tP0+1;
snP0 = out_P0.res(tP0).s_n;
tauP0 = out_P0.res(tP0).tauNorm;
dz = 0.5/numel(out_P0.res(1).gap);
zP0 = linspace(dz,15-dz,numel(out_P0.res(1).gap));
snConf = out_conforming.res(tConf).s_n;
tauConf = out_conforming.res(tConf).tauNorm;
snDual = out_dual.res(tConf).s_n;
tauDual = out_dual.res(tConf).tauNorm;
zConf = linspace(0,15,numel(out_conforming.res(1).gap));

figure(1)
plot(snP0,zP0,'g-o','LineWidth',0.5,'MarkerSize',ms,'MarkerEdgeColor','g','MarkerFaceColor','g');
hold on
plot(snConf,zConf,'k-','LineWidth',1.5);
plot(snDual,zConf,'r-o','LineWidth',0.5,'MarkerSize',ms,'MarkerEdgeColor','r','MarkerFaceColor','r');
legend('P0 stabilized','Dual Conforming','Dual non conforming')
set(gca, 'FontName', 'Times', 'FontSize', 16, 'TickLabelInterpreter', 'latex');

figure(2)
plot(tauP0,zP0,'g-o','LineWidth',0.5,'MarkerSize',ms,'MarkerEdgeColor','g','MarkerFaceColor','g');
hold on
plot(tauConf,zConf,'k-s','LineWidth',1.5);
plot(tauDual,zConf,'r-o','LineWidth',0.5,'MarkerSize',ms,'MarkerEdgeColor','r','MarkerFaceColor','r');
legend('P0 stabilized','Dual Conforming','Dual non conforming')
set(gca, 'FontName', 'Times', 'FontSize', 16, 'TickLabelInterpreter', 'latex');



%figure(2)

