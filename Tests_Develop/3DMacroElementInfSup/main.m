% simple 3D patch test to check minimal kernel in 3D

% set mesh 
Nm = 2;
Ns = 4;
b1 = BlockStructuredMesh([0,1;0 1;1 1+1/Nm],[Nm Nm 1],1);
mesh1 = processGeometry(b1);
b2 = BlockStructuredMesh([0,1;0 1;1-1/Ns 1],[Ns Ns 1],1);
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
bc1 = Boundaries("BCs/fixTop.xml",grid1);
bc2 = Boundaries("BCs/fixBot.xml",grid2);

% Create object handling construction of Jacobian and rhs of the model
domain1 = Discretizer('Boundaries',bc1,...
                     'OutState',printUtils1,...
                     'Materials',mat1,...
                     'Grid',grid1);

domain2 = Discretizer('Boundaries',bc2,...
                     'OutState',printUtils2,...
                     'Materials',mat2,...
                     'Grid',grid2);

domains = [domain1; domain2];

% build the mortar interface with xml input
[interfaces,domains] = Mortar.buildInterfaces('interfaces.xml',domains);


% get mortar matrices and remove boundary dofs
interfaces{1}.computeMortarMatrices();
% get boundary dof
[bcDof1,~] = getBC(getPhysicsSolver(domain1,"displacements"),"FixTop",0);
[bcDof2,~] = getBC(getPhysicsSolver(domain2,"displacements"),"FixBot",0);

% get interface dofs belonging to bc dofs
loc2glob1 = DoFManager.dofExpand(interfaces{1}.mesh.local2glob{1},3);
loc2glob2 = DoFManager.dofExpand(interfaces{1}.mesh.local2glob{2},3);

isBC1 = ismember(loc2glob1,bcDof1);
isBC2 = ismember(loc2glob2,bcDof2);
activeDof1 = loc2glob1(~isBC1);
activeDof2 = loc2glob2(~isBC2);

M = interfaces{1}.M(:,activeDof1);
D = interfaces{1}.D(:,activeDof2);


% get kernel of B^T (same as the Schur complement!)
% this kernel is crazy!
B = full([D -M]');
[u,e,w] = svd(B);

%get matrix A removing extra entries
getPhysicsSolver(domain1,"Poromechanics").computeMat([],[]);
getPhysicsSolver(domain2,"Poromechanics").computeMat([],[]);
K1 = getPhysicsSolver(domain1,"Poromechanics").K;
K2 = getPhysicsSolver(domain2,"Poromechanics").K;

K1(:,bcDof1) = [];  K1(bcDof1,:) = [];
K2(:,bcDof2) = [];  K2(bcDof2,:) = [];
K = blkdiag(K2,K1);


%% check kernel orthogonality

S = B'*(K\B);

computeStabilizationMatrix(interfaces{1}); 
H = obj.stabilizationMat;

nH = null(full(H));
nS = null(full(S));


% test kernel 
% checkered kernel
v1 = [1;0;0;-1;0;0;1;0;0;0;0;0];
v2 = [0;0;0;-1;0;0;0;0;0;1;0;0];
v = [v1;v2;v1;v2];
% 
% [v,e] = eig(S);
% v = real(v); e = real(e);

fprintf("norm(S*v) = %2.4e \n",norm(S*v))

for i = 1:size(nS,2)
  r1 = nS(:,i)'*nH(:,1);
  r2 = nS(:,i)'*nH(:,2);
  r3 = nS(:,i)'*nH(:,3);
  fprintf('%1.4e %1.4e %1.4e \n',r1,r2,r3);
end

%% inf-sup constant evaluation

h_s = 1/ 









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




