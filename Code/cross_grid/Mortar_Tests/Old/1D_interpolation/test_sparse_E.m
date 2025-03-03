% 1D interpolation between two 1D interfaces
% The code exploits Mortar1D class to compute the interpolation operator
% and perform error analysis
close all; clear
%% INPUT DATAS

% INPUT PARAMETERS
degree = 1;

type = 'gauss';

%nSizes = 1; % number of uniform refinment for convergence analysis
nMnodes = 20; % number of nodes for the first mesh refinment
ratVals = 0.2:0.2:3;
vals = zeros(length(ratVals),nMnodes);
i=1;
for ratio = ratVals  % loop trough different mesh refinments
    nMaster = nMnodes;
    nSlave = round(ratio*nMaster);
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
    % [D_EB,M_EB, E_EB] = mortar.computeMortarElementBased(nGP);
    % [D_SB,M_SB,E_SB] = mortar.computeMortarSegmentBased(nGP);
    
    E_RBF = D_RBF\M_RBF;
    % E_RBF_mod = E_RBF;
    % 
    % D_RBF_lump = diag(sum(D_RBF,2));
    % E_RBF_lump = D_RBF_lump\M_RBF;
    % 
    % k = 1e-5;
      
    % remove useless components to make it sparse
    % E_RBF_mod(abs(E_RBF) < k) = 0;
    % E_RBF_mod = E_RBF_mod./sum(E_RBF_mod,2);
    % spE_RBF = nnz(E_RBF_mod)/numel(E_RBF_mod);
    nr = round(size(E_RBF,1)/2);
    vals(i,:) = E_RBF(nr,:);
    i = i+1;
end

%% PLOT VALUES OF E ON A GENERIC ROW
plot(1:nMaster,abs(vals(1,:)),'k','LineWidth',1.2)
hold on
plot(1:nMaster,abs(vals(9,:)),'r','LineWidth',1.2)





















