% 1D interpolation between two 1D interfaces
% The code exploits Mortar1D class to compute the interpolation operator
% and perform error analysis
close all; clear
%% INPUT DATAS

% INPUT PARAMETERS
degree = 1;

type = 'gauss';

nSizes = 8; % number of uniform refinment for convergence analysis
nMnodes = 5; % number of nodes for the first mesh refinment
errNormRBF = zeros(nSizes,1);
errNormSB = errNormRBF;
errNormEB = errNormRBF;
errNormRBF_lump = errNormRBF;
for sizeCount = 1:nSizes  % loop trough different mesh refinments
    nMaster = nMnodes*(2^(sizeCount-1));
    nSlave = round(1.3*nMaster);
    master = zeros(nMaster,2); % coordinates of master side
    slave = zeros(nSlave,2);
    x1 = -1;
    x2 = 1;
    x3 = -1;
    x4 = 1;
    % x-coordinates
    master(:,1) = (linspace(x1,x2,nMaster))';
    slave(:,1) = (linspace(x3,x4,nSlave))';

    % y-coordinates - k controls the curvature of each grid (non
    % conformity)
    k = 0;
    curve = @(x) k*x.^2 - k;
    master(:,2) = curve(master(:,1));
    slave(:,2) = curve(slave(:,1));

    % number of RBF interpolation points for each element
    nInt = 4;

    % Number of integration points for RBF testing (GP class taken from GReS)
    nGP = 3;

    % Build a topology matrix for master/slave surfs based on nodes position
    mastertop = build_topol(master(:,1));
    slavetop = build_topol(slave(:,1));

    % Compute mortar operator
    mortar = Mortar2D(1,'set',mastertop,slavetop,master,slave);
    [D_RBF,M_RBF] = mortar.computeMortarRBF(nGP,nInt,type);
    [D_EB,M_EB, E_EB] = mortar.computeMortarElementBased(nGP);
    [D_SB,M_SB,E_SB] = mortar.computeMortarSegmentBased(nGP);
    
    E_RBF = D_RBF\M_RBF;

    D_RBF_lump = diag(sum(D_RBF,2));
    E_RBF_lump = D_RBF_lump\M_RBF;
    
    % remove useless components to make it sparse
    E_RBF(abs(E_RBF) < 1e-4) = 0;
    E_RBF = E_RBF./sum(E_RBF,2);
    spE_RBF = nnz(E_RBF)/numel(E_RBF);
    % Analytical function to interpolate
    f = @(x) sin(4*x) + x.^2;
    
    % Compute interpolation error
    errNormSB(sizeCount) = mortar.computeInterpError(E_SB,f);
    errNormEB(sizeCount) = mortar.computeInterpError(E_EB,f);
    errNormRBF(sizeCount) = mortar.computeInterpError(E_RBF,f);
    errNormRBF_lump(sizeCount) = mortar.computeInterpError(E_RBF_lump,f);
    %errNormRBF_w(sizeCount) = mortar.computeInterpError(E_RBF_w,f);
end

%% PLOT CONVERGENCE PROFILE AND INTERPOLATED FUNCTION

name = strcat('L2_',type,'_Int',num2str(nInt));
fID = fopen(strcat('Results_lin\',name,'.dat'),'w');
fprintf(fID,'%2.6e \n',errNormRBF);

fID = fopen(strcat('Results_lin\L2_eb.dat'),'w');
fprintf(fID,'%2.6e \n',errNormEB);


% convergence profile
figure(1)
h = 1./(nMnodes*2.^(0:nSizes-1));
fID = fopen(strcat('Results_lin\h.dat'),'w');
fprintf(fID,'%2.6e \n',h);
%loglog(h,errNormSB,'r-o')
loglog(h,errNormRBF,'b-s')
hold on
loglog(h,errNormEB,'g-^')
loglog(h,errNormRBF_lump,'r-o')
legend('RBF','Element-based','lumped RBF')

%% Plot results of analytical function





















