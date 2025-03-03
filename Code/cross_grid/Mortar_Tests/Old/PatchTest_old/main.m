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

nel = 30;
nXs = nel+1;
rat = 0.3;
nYs = round(nXs);
nXm = round(nel*rat+1);
nYm = round(nXm);
h = 1/(rat*nel);
alpha = 10;
gamma = (alpha*h)/E1;



stab = 'unstable';
nTip = 0;
%%
fprintf('Processing non conforming with standard multipliers \n')
[standard_tx,standard_ty,standard_us,standard_um] = RunNonConfPatchTest('standard',Fx,Fy,DmatS,DmatM,patch,stab,nXs,nYs,nXm,nYm);
fprintf('Processing non conforming with dual multipliers \n')
[dual_tx,dual_ty,dual_us,dual_um] = RunNonConfPatchTest('dual',Fx,Fy,DmatS,DmatM,patch,stab,nXs,nYs,nXm,nYm);
fprintf('Processing Unbiased \n')
[const_tx,const_ty,const_us,const_um,x] = RunNonConfConstant(Fx,Fy,DmatS,DmatM,patch,nXs,nYs,nXm,nYm,E1,E2);
fprintf('Processing Unbiased \n')
[ub_tx,ub_ty,ub_us,ub_um] = RunPuso2020(Fx,Fy,DmatS,DmatM,gamma,patch,nXs,nYs,nXm,nYm);
fprintf('Processing conforming \n')
[conf_tx,conf_ty,conf_us,conf_um] = RunConfPatch(Fx,Fy,DmatS,DmatM,patch,'unstable',nXs,nYs,nYm);

xNodeM = linspace(0,1,size(dual_um,1));
xNodeS = linspace(0,1,size(dual_us,1));
xNodeConf = linspace(0,1,size(conf_us,1));



%% plotting
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultAxesFontName','Times')

ms = 12;
figure(1)
tiledlayout(1,2)

% Plot tx
nexttile
plot(xNodeS, dual_tx, 'b.-', 'LineWidth', 1, 'MarkerSize', ms)
hold on
plot(xNodeS, standard_tx, 'r.-', 'LineWidth', 1, 'MarkerSize', ms)
plot(x, const_tx, 'g.-', 'LineWidth', 1, 'MarkerSize', ms)
plot(xNodeS, ub_tx, 'm.-', 'LineWidth', 1, 'MarkerSize', ms)
plot(xNodeConf, conf_tx, 'k-', 'LineWidth', 1.5)
xlabel('x (m)', 'Interpreter', 'latex')
ylabel('$t_x$ (kPa)', 'Interpreter', 'latex')
legend('Dual', 'Standard', 'Constant', 'Unbiased', 'Conforming', 'Location', 'best', 'Interpreter', 'latex')
grid on
title('$t_x$ distribution', 'Interpreter', 'latex')
if patch == 1
   ylim([-0.1 0.1])
end

% Plot ty
nexttile
plot(xNodeS, dual_ty, 'b.-', 'LineWidth', 1, 'MarkerSize', ms)
hold on
plot(xNodeS, standard_ty, 'r.-', 'LineWidth', 1, 'MarkerSize', ms)
plot(x, const_ty, 'g.-', 'LineWidth', 1, 'MarkerSize', ms)
plot(xNodeS, ub_ty, 'm.-', 'LineWidth', 1, 'MarkerSize', ms)
plot(xNodeConf, conf_ty, 'k-', 'LineWidth', 1)
xlabel('x (m)', 'Interpreter', 'latex')
ylabel('$t_y$ (kPa)', 'Interpreter', 'latex')
legend('Dual', 'Standard', 'Constant', 'Unbiased', 'Conforming', 'Location', 'best', 'Interpreter', 'latex')
grid on
title('$t_y$ distribution', 'Interpreter', 'latex')
if patch == 1
   ylim([-11 -9])
end

exportgraphics(gcf,strcat('Patch_',num2str(patch),'.png'),'Resolution',400)

hold off

