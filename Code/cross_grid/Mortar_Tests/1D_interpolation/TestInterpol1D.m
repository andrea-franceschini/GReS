% 1D interpolation between two 1D interfaces
% The code exploits Mortar1D class to compute the interpolation operator
% and perform error analysis
close all; clear
%% INPUT DATAS

% INPUT PARAMETERS
degree = 1;

nSizes = 4; % number of uniform refinment for convergence analysis
nMnodes = 8; % number of nodes for the first mesh refinment
errNormRBF = zeros(nSizes,1);
errNormSB = errNormRBF;
errNormEB = errNormRBF;
for sizeCount = 1:nSizes  % loop trough different mesh refinments
    nMaster = nMnodes*(2^(sizeCount-1));
    nSlave = round(2*nMaster);
    master = zeros(nMaster,2); % coordinates of master side
    slave = zeros(nSlave,2);
    x1 = -1;
    x2 = 1;
    % x-coordinates
    master(:,1) = (linspace(x1,x2,nMaster))';
    slave(:,1) = (linspace(x1,x2,nSlave))';

    % y-coordinates - k controls the curvature of each grid (non
    % conformity)
    k = 0;
    curve = @(x) k*x.^2 - k;
    master(:,2) = curve(master(:,1));
    slave(:,2) = curve(slave(:,1));

    % number of RBF interpolation points for each element
    nInt = 10;

    % Number of integration points for RBF testing (GP class taken from GReS)
    nGP = 15;

    % Build a topology matrix for master/slave surfs based on nodes position
    mastertop = build_topol(master(:,1));
    slavetop = build_topol(slave(:,1));

    % Compute mortar operator
    mortar = Mortar1D(1,'set',mastertop,slavetop,master,slave);
    [E_RBF,~,tRBF] = mortar.computeMortarRBF(nGP,nInt);
    [E_EB,~,tEB] = mortar.computeMortarElementBased(nGP);
    [E_SB] = mortar.computeMortarSegmentBased(nGP);

    % Analytical function to interpolate
    f = @(x) sin(3*x) + exp(2*x);

    % Compute interpolation error
    errNormSB(sizeCount) = mortar.computeInterpError(E_SB,f);
    errNormEB(sizeCount) = mortar.computeInterpError(E_EB,f);
    errNormRBF(sizeCount) = mortar.computeInterpError(E_RBF,f);

end

%% PLOT CONVERGENCE PROFILE AND INTERPOLATED FUNCTION

% convergence profile
figure(1)
h = 1./(nMnodes*2.^(0:nSizes-1));
loglog(h,errNormSB,'r-o')
hold on
loglog(h,errNormRBF,'b-s')
loglog(h,errNormEB,'g-^')
legend('Segment-based', 'Radial Basis Functions', 'Element-based')





















