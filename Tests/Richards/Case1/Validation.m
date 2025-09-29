% This file is a comparison of two approaches to describing two variables,
% the saturation and the relative permeability. The first approach proposes
% an analytical expression for these variables, while the second approach
% relies on interpolated values from a table.

close all;
clear;
clc;
figures_dir = 'Figs/';
image_dir = strcat(pwd,'/',figures_dir);
load("Inputs/Solution/output1B.mat")
pressplotA=pressplot;
swplotA=swplot;

load("Inputs/Solution/output3B.mat")
pressplotC=pressplot;
swplotC=swplot;

t = [10 50 100];
t_max = t(end);
t = t/t_max;
tAnalytical = strcat("Ana - ",num2str(t',"%.1f")," T");
tTabular = strcat("Tab - ",num2str(t',"%.1f")," T");
tstr = strcat([tAnalytical; tTabular]);


figure('Position', [100, 100, 700, 700])
hold on
plot(pressplotA./(pressplotA(pos,:)),ptsZ,'k-', 'LineWidth', 2, 'MarkerSize', 14);
plot(pressplotC./(pressplotC(pos,:)),ptsZ,'k.', 'LineWidth', 2, 'MarkerSize', 14);
xlabel('p/p_{top}')
ylabel('z/H')
legend(tstr, 'Location', 'southeast', 'NumColumns', 2)
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = strcat(image_dir,'pressure_validation','.png');
exportgraphics(gcf,stmp,'Resolution',400)


figure('Position', [100, 100, 700, 700])
hold on
plot(swplotA,ptsZ,'k-', 'LineWidth', 2, 'MarkerSize', 14);
plot(swplotC,ptsZ,'k.', 'LineWidth', 2, 'MarkerSize', 14);
xlabel('S_w')
ylabel('z/H')
legend(tstr, 'Location', 'southwest', 'NumColumns', 2)
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = strcat(image_dir,'saturation_validation', '.png');
exportgraphics(gcf,stmp,'Resolution',400)
