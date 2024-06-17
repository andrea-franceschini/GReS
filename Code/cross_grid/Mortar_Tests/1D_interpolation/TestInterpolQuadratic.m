% 1D interpolation between two 1D interfaces
% The code exploits Mortar2D class to compute the interpolation operator
% and perform error analysis
% This scripts tests quadratic basis functions: effect of Gauss points
% evene befor the effect of the interpolation?
close all; clear
warning('off','MATLAB:nearlySingularMatrix');
%% INPUT DATAS

% INPUT PARAMETERS
type = 'wendland';

nSizes = 9; % number of uniform refinment for convergence analysis
nMnodes = 8; % number of nodes for the first mesh refinment
errNormRBF = zeros(nSizes,1);
errNormSB = errNormRBF;
errNormEB = errNormRBF;
for sizeCount = 1:nSizes  % loop trough different mesh refinments
    nMaster = nMnodes*(2^(sizeCount-1));
    nSlave = round(1.5*nMaster);
    % setting odd number of nodes
    if mod(nMaster,2) == 0
        nMaster = nMaster+1;
    end
    if mod(nSlave,2) == 0
        nSlave = nSlave+1;
    end
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
    nInt = 8;

    % Number of integration points for RBF testing (GP class taken from GReS)
    nGP = 3;

    % Build a topology matrix for master/slave surfs based on nodes position
    mastertop = build_topol_quad(master(:,1));
    slavetop = build_topol_quad(slave(:,1));

    % Compute mortar operator
    mortar = Mortar2D(2,'set',mastertop,slavetop,master,slave);
    [E_RBF,~,tRBF] = mortar.computeMortarRBF(nGP,nInt,type);
    [E_EB,~,tEB] = mortar.computeMortarElementBased(nGP);
    %[E_SB] = mortar.computeMortarSegmentBased(nGP);

    % Analytical function to interpolate
    f = @(x) sin(4*x) + x.^2;

    % Compute interpolation error
    %errNormSB(sizeCount) = mortar.computeInterpError(E_SB,f);
    errNormEB(sizeCount) = mortar.computeInterpError(E_EB,f);
    errNormRBF(sizeCount) = mortar.computeInterpError(E_RBF,f);
end

%% PLOT CONVERGENCE PROFILE AND INTERPOLATED FUNCTION
name = strcat('L2_',type,'_Int',num2str(nInt));
fID = fopen(strcat('Results\',name,'.dat'),'w');
fprintf(fID,'%2.6e \n',errNormRBF);

fID = fopen(strcat('Results\L2_eb.dat'),'w');
fprintf(fID,'%2.6e \n',errNormEB);



figure(1)
h = 1./(nMnodes*2.^(0:nSizes-1));
fID = fopen(strcat('Results\h.dat'),'w');
fprintf(fID,'%2.6e \n',h);
%loglog(h,errNormSB,'r-o')
loglog(h,errNormRBF,'b-s')
hold on
loglog(h,errNormEB,'g-^')
legend('Element-based')

% convergence rates
rateRBF = errNormRBF(1:end-1)./errNormRBF(2:end);
rateEB = errNormEB(1:end-1)./errNormEB(2:end);


%%


% % convergence profile
% eb = load("EB.mat");
% ebNorm = eb.errNormEB;
% rbf4 = load("RB_int4.mat");
% rbf4Norm = rbf4.errNormRBF;
% rbf8 = load("RB_int8.mat");
% rbf8Norm = rbf8.errNormRBF;
% rbf_wend = load("RB_wendInt8.mat");
% rbf_wend = rbf_wend.errNormRBF;
% 
% figure(1)
% h = 1./(nMnodes*2.^(0:nSizes-1));
% %loglog(h,errNormSB,'r-o')
% loglog(h,ebNorm,'g-^')
% hold on
% loglog(h,rbf4Norm,'b-^')
% loglog(h,rbf8Norm,'r-^')
% loglog(h,rbf_wend,'k-^')
% legend('Element-based','Gauss nInt = 4','Gauss nInt = 8','Wendland nInt = 8')






















