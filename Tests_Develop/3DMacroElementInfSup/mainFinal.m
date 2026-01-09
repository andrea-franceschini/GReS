% simple 3D patch test to check minimal kernel in 3D

clear
close all
clc
% set mesh
NM = [2 4 6 8 10 12 14];
NS = [2 2 2 2 2 2 2];
% NM = 6;
% NS = 12;

infSup = zeros(numel(NM),1);
infSupStab1 = zeros(numel(NM),1);

for r = 1:numel(NM)

  fprintf("Processing Refinement %i \n",r)
  Nm = NM(r);
  Ns = NS(r);

  b1 = BlockStructuredMesh([0,1;0 1;1 2],[Nm Nm 1],1);
  mesh1 = processGeometry(b1);
  b2 = BlockStructuredMesh([0,1;0 1;0 1],[Ns Ns 1],1);
  mesh2 = processGeometry(b2);

  % define model

  simParam = SimulationParameters("infSupTest.xml");

  elems1 = Elements(mesh1,2);
  faces1 = Faces(mesh1);
  grid1 = struct('topology',mesh1,'cells',elems1,'faces',faces1);
  mat1 = Materials("infSupTest.xml");

  elems2 = Elements(mesh2,2);
  faces2 = Faces(mesh2);
  grid2 = struct('topology',mesh2,'cells',elems2,'faces',faces2);
  mat2 = Materials("infSupTest.xml");

  % Create and set the print utility
  printUtils1 = OutState(mesh1,"folderName","OUT/cubeTop","timeList",1);
  printUtils2 = OutState(mesh1,"folderName","OUT/cubeBot","timeList",1);

  % Create an object of the "Boundaries" class
  bcTop = Boundaries("BCs/fixTop.xml",grid1);
  bcBot = Boundaries("BCs/fixBot.xml",grid2);

  % Create object handling construction of Jacobian and rhs of the model
  domainTop = Discretizer('Boundaries',bcTop,...
    'OutState',printUtils1,...
    'Materials',mat1,...
    'Grid',grid1);

  domainBot = Discretizer('Boundaries',bcBot,...
    'OutState',printUtils2,...
    'Materials',mat2,...
    'Grid',grid2);

  domainTop.addPhysicsSolver("infSupTest.xml");
  domainBot.addPhysicsSolver("infSupTest.xml");

  domainTop.simparams = simParam;
  domainBot.simparams = simParam;

  domains = [domainTop; domainBot];

  % build the mortar interface with xml input
  interfaces = buildInterfaces('infSupTest.xml',domains);


  % get mortar matrices and remove boundary dofs
  interfaces{1}.computeConstraintMatrices();

  % get boundary dof
  [bcDof1,~] = getBC(getPhysicsSolver(domainTop,"Poromechanics"),"fixTop",0);
  [bcDof2,~] = getBC(getPhysicsSolver(domainBot,"Poromechanics"),"fixBot",0);

  % get interface dofs belonging to bc dofs
  loc2glob1 = DoFManager.dofExpand(interfaces{1}.interfMesh.local2glob{1},3);
  loc2glob2 = DoFManager.dofExpand(interfaces{1}.interfMesh.local2glob{2},3);

  isBC1 = ismember(loc2glob1,bcDof1);
  isBC2 = ismember(loc2glob2,bcDof2);
  % activeDof1 = loc2glob1(~isBC1);
  % activeDof2 = loc2glob2(~isBC2);


  M = interfaces{1}.M;
  D = interfaces{1}.D;

  M(:,bcDof1) = [];
  D(:,bcDof2) = [];


  % get kernel of B^T (same as the Schur complement!)
  % this kernel is crazy!
  B = full([D -M]');
 % [u,e,w] = svd(B);

  %get matrix A removing extra entries
  poroTop = getPhysicsSolver(domainTop,"Poromechanics");
  poroTop.domain.stateOld = copy(poroTop.domain.state);
  poroBot = getPhysicsSolver(domainBot,"Poromechanics");
  poroBot.domain.stateOld = copy(poroBot.domain.state);
  poroBot.assembleSystem([]);
  poroTop.assembleSystem([]);
  K1 = poroTop.K;
  K2 = poroBot.K;

  K1(:,bcDof1) = [];  K1(bcDof1,:) = [];
  K2(:,bcDof2) = [];  K2(bcDof2,:) = [];
  K = blkdiag(K2,K1);


  %% check kernel orthogonality

  S = B'*(K\B);

  interfaces{1}.assembleConstraint();

  H = interfaces{1}.stabilizationMat;

  %Hold = interfaces{1}.computeStabilizationMatrixOld();

  % scale H to investigate the stability effectiveness
  %H = 1e-1*H;
  % 
  % nH = null(full(H));
  % nS = null(full(S));


  % test kernel
  % checkered kernel
  % v1 = [1;0;0;1;0;0;1;0;0;1;0;0];
  % v2 = [-1;0;0;-1;0;0;-1;0;0;-1;0;0];
  % v = [v1;v2;v1;v2];
  %
  %[v,e] = eig(S);
  %v = real(v); e = real(e);

%  fprintf("norm(S*v) = %2.4e \n",norm(S*v))

  % for i = 1:size(nS,2)
  %   r1 = nS(:,i)'*nH(:,1);
  %   r2 = nS(:,i)'*nH(:,2);
  %   r3 = nS(:,i)'*nH(:,3);
  %   fprintf('%1.4e %1.4e %1.4e \n',r1,r2,r3);
  % end

  % inf-sup constant evaluation (according to ElAbbasi_2001 paper

  h_s = 1/Ns;

  % scaled mass matrix (h*Q)^-1
  invQ = diag(1./(h_s^3*ones(size(S,1),1)));

  % Q = diag(h_s^2*ones(size(S,1),1));
  % 
  % 
  % v = [1e-5 1e-2 1e-1 0.5 1 1.5 2 5 10 20 100 1e3 1e4 1e5];
  % 
  % infSup1 = sqrt(min(real(eig(invQ*S))));
  % 
  % 
  eUnstab = abs(real(eig(invQ*S)));
  eStab = eig(invQ*(S + H));
  infSup(r) = sqrt(min(eUnstab));
  infSupStab1(r) = sqrt(min(eStab));
  % % eNew = real(eig(S+H,Q));
  % % infSupStab2(r) = sqrt(min(eNew));
  % 
  % % stabilization computed as in Franc 2022
  % % Hfranc = H;
  % % Hfranc(Hfranc~=0) = -1;
  % % Hfranc = Hfranc + diag(2*ones(size(Hfranc,1),1));
  % % Hfranc = h_s^3*Hfranc;
  % 
  % 
  % % rescale H using Schur complement sum
  % 
  % infSupStabOld = zeros(numel(v),1);
  % 
  % infSupStabNew = zeros(numel(v),1);
  % 
  % 
  % figure
  % for i = 1:numel(v)
  %   scaledSchurStabOld = invQ*(S+v(i)*Hold);
  %   scaledSchurStabNew = invQ*(S+v(i)*H);
  %   infSupStabOld(i) = sqrt(min(eig(scaledSchurStabOld)));
  %   infSupStabNew(i) = sqrt(min(eig(scaledSchurStabNew)));
  % end
  % 
  % loglog(v,infSupStabOld,'k-o',"LineWidth",2)
  % grid on
  % hold on
  % loglog(v,infSupStabNew,'b-s',"LineWidth",2)
  % loglog(v(v==1),infSupStabOld(v==1.0),'ro')
  % loglog(v(v==1),infSupStabNew(v==1.0),'rs')
  % xlabel("Stabilization scaling")
  % ylabel("Inf-sup value")

  %fprintf("\n Inf-sup constant value for unstabilized: %1.4e \n",infSup)
   %fprintf("\n Inf-sup constant value for stabilized franc: %1.4e \n",infSupStabFranc)


end

%%

figure('Color','w','Position',[100 100 520 420])

plot(log(1./NM), infSupStab1, 'k-s', ...
     'LineWidth', 1.3, ...
     'MarkerSize', 7, ...
     'MarkerFaceColor','k');
hold on

plot(log(1./NM), infSup, 'r-o', ...
     'LineWidth', 1.3, ...
     'MarkerSize', 7, ...
     'MarkerFaceColor','r');

% ---- Legend (LaTeX) ----
lgd = legend({'Unstabilized', 'Stabilized'}, ...
             'Interpreter','latex', ...
             'FontSize',12, ...
             'Location','best');
lgd.Box = 'off';

xlabel('$h_1$', 'Interpreter','latex', 'FontSize',14)
ylabel('$\beta^*$', 'Interpreter','latex', 'FontSize',14)

set(gca, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',12, ...
    'LineWidth',1.1)

xlim([-3 -1.2])
ylim([-0.1 1])

grid on

% ---- Export to PDF (vector) ----
exportgraphics(gcf,'infSupScaling3.pdf','ContentType','vector')






%%
% % Eigen-decomposition
% [U,D] = eig(S);
%
% % Sort eigenvalues (optional — use 'descend' if you prefer)
% [d, idx] = sort(diag(D), 'ascend');
% D = diag(d);
% V = V(:, idx);
%
% % Enforce consistent sign convention:
% % Make each eigenvector have positive largest component
% for i = 1:size(V,2)
%     [~, jmax] = max(abs(V(:,i)));
%     if V(jmax,i) < 0
%         V(:,i) = -V(:,i);
%     end
% end
%
% % 4Reorder columns to best match the canonical basis
% % (optional but useful if A is nearly diagonal)
% perm = zeros(1, size(V,2));
% used = false(1, size(V,2));
% for i = 1:size(V,2)
%     % find which canonical direction e_j is closest to this vector
%     [~, j] = max(abs(V(:,i)));
%     % if that column j is already taken, pick next best
%     while j <= size(V,2) && used(j)
%         [~, j] = max(abs(V(:,i)) .* ~used);
%     end
%     if j > size(V,2), j = i; end
%     perm(i) = j;
%     used(j) = true;
% end
% V = V(:, perm);
% D = D(perm, perm);  % reorder eigenvalues accordingly
%
% % Display results
% disp('Eigenvalues (sorted and reordered):');
% disp(diag(D));
% disp('Eigenvectors (aligned to canonical basis):');
% disp(V);




