% This file compares the analysis carried out in GReS using the analytical
% solution against that of MRST.
% The simulations are stored in .mat files and, if necessary, the analyses
% carried out by GReS can be updated by running the Main* files again.

% Note:
% - For comparison, the same conditions must be guaranteed in both
% GReS and MRST. That's why it's important to compare the same mesh
% refinement scheme and the same control for the non-linear solver.

%% CONVERGENCY IN TIME
close all;
clear;
clc;
figures_dir = 'Figs/';
load("Inputs/Solution/MRST_30.mat")
load("Inputs/Solution/GReS_30.mat")

tind = 1:length(MRST_time);
t = MRST_time(tind)/86400;
tMRST = strcat("MRST ",num2str(t',"%.1f")," Days");
tind = 1:length(GReS_time);
t = GReS_time(tind)/86400;
tGReS = strcat("GReS ",num2str(t,"%.1f")," Days");
tstr = strcat([tGReS; tMRST]);

figure('Position', [100, 100, 700, 700])
hold on
plot(GReS_pres,GReS_elev,'k-', 'LineWidth', 2, 'MarkerSize', 14);
plot(MRST_pres,MRST_elev,'k.', 'LineWidth', 2, 'MarkerSize', 18);
xlabel('Head Pressure (m)')
ylabel('Height (m)')
legend(tstr, 'Location', 'northwest', 'NumColumns', 2)
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = strcat(figures_dir,'GReSxMRST_head_pressure','.png');
exportgraphics(gcf,stmp,'Resolution',400)
    
figure('Position', [100, 100, 700, 700])
hold on
plot(GReS_satu,GReS_elev,'k-', 'LineWidth', 2, 'MarkerSize', 14);
plot(MRST_satu,MRST_elev,'k.', 'LineWidth', 2, 'MarkerSize', 18);
xlabel('Saturation')
ylabel('Height (m)')
legend(tstr, 'Location', 'northwest', 'NumColumns', 2)
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = strcat(figures_dir,'GReSxMRST_saturation','.png');
exportgraphics(gcf,stmp,'Resolution',400)

%% MESH CONVERGENCY 
load("Inputs/Solution/GReS_Mesh_Conv.mat")
load("Inputs/Solution/MRST_Mesh_Conv.mat")

tMRST = strcat("MRST ",num2str(GReS_mesh',"%i")," Cells");
tGReS = strcat("GReS ",num2str(GReS_mesh',"%i")," Cells");
tstr = strcat([tGReS; tMRST]);

figure('Position', [100, 100, 700, 700])
hold on
for i=1:length(GReS_result)
    plot(-GReS_result(i).pressure/GReS_weight,GReS_result(i).height,'k-', 'LineWidth', 2, 'MarkerSize', 14);
end
for i=1:length(MRST_result)
    plot(-MRST_result(i).pressure,MRST_result(i).height,'k.', 'LineWidth', 2, 'MarkerSize', 14);
end
xlabel('Head Pressure (m)')
ylabel('Height (m)')
legend(tstr, 'NumColumns', 2)
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = strcat(figures_dir,'GReSxMRST_mesh_pressure','.png');
exportgraphics(gcf,stmp,'Resolution',400)

%Plotting saturation
figure('Position', [100, 100, 700, 700])
hold on
for i=1:length(GReS_result)
    plot(GReS_result(i).saturation,GReS_result(i).height,'k-', 'LineWidth', 2, 'MarkerSize', 14);
end
for i=1:length(MRST_result)
    plot(MRST_result(i).saturation,MRST_result(i).height,'k.', 'LineWidth', 2, 'MarkerSize', 14);
end
xlabel('Saturation')
ylabel('Height (m)')
legend(tstr, 'Location', 'northwest', 'NumColumns', 2)
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = strcat(figures_dir, 'GReSxMRST_mesh_saturation', '.png');
exportgraphics(gcf,stmp,'Resolution',400)