clear
close all
clc

%-------------------------------
% Settings
%-------------------------------

ng = 4:400;                 % number of Gauss points
el_type = 'hexa27';         % element type
nInt = 25;                  % number of integration points
nIntSupp = 25;              % same as above
it_avg = 3;                 % average iterations for EB scheme
nN = 9;                     % number of nodes per element for hexa27

% Cost constants
c_hexa27 = 250;             % cost for EB (hexa27)
c_hexa8  = 118;             % cost for EB (hexa8)

%----------------------------------------
% Initialize figure
%----------------------------------------

figure(1); hold on; grid on;
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

%----------------------------------------
% Loop over refinement ratios
%----------------------------------------
lw=1.5;

for rh = 2
    
    % ----- EB cost -----
    if rh >= 1
        R_slave = 0.5 * (rh^2 + 1);
    else
        R_slave = 1;
    end

    cost_per_gp = R_slave * (it_avg * c_hexa27) + 16 + nN;
    eb_cost = @(x) x * cost_per_gp;
    eb = eb_cost(ng);

    if rh == 0.5
        plot(ng, eb, 'r-','LineWidth',lw, 'DisplayName', 'EB - $r=1/2$');
    else
        plot(ng, eb, 'r-','LineWidth',lw, 'DisplayName', 'EB');
    end

    % ----- RBF cost -----
    if rh >= 1
        R_master = rh^2;
        R_slave2 = 1;
    else
        R_master = 1;
        R_slave2 = 1;
    end

    c_rbf = R_slave * 2 * (2 * nInt) + R_slave2 * 9 * (2 * nInt);
    rbf_cost = @(x) R_master * ((1/3) * nInt^3 + (nN + 3) * nInt^2) + x * (c_rbf);
    rbf = rbf_cost(ng);

    if rh == 0.5
        plot(ng, rbf, 'b-', 'LineWidth',lw, 'DisplayName', 'RBF - $r=1/2$');
    else
        plot(ng, rbf, 'b-', 'LineWidth',lw, 'DisplayName', 'RBF');
    end
end

%----------------------------------------
% Labels and export
%----------------------------------------

ax = gca;
ax.FontSize = 14;               % Tick labels

xlabel('$N_G$', 'FontSize', 16)
ylabel('$\mathrm{FLOPs}$', 'FontSize', 16)
legend('Location', 'northwest', 'FontSize', 14)

exportgraphics(gcf, 'cost_comparison_eb_rbf.pdf', 'ContentType', 'vector');