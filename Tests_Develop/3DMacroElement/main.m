% simple 3D patch test to check minimal kernel in 3D

% build mesh
b1 = BlockStructuredMesh([0,1;0 1;1 2],[2 2 1],1);
mesh1 = processGeometry(b1);
b2 = BlockStructuredMesh([0,1;0 1;0 1],[4 4 1],1);
mesh2 = processGeometry(b2);

% define model 
model = ModelType("Poromechanics_FEM");
simParam = SimulationParameters("simParam.dat");

elems1 = Elements(mesh1,2);
faces1 = Faces(model, mesh1);
grid1 = struct('topology',mesh1,'cells',elems1,'faces',faces1);
mat1 = Materials(model,"materialsList.dat");
dofm1 = DoFManager(mesh1,model);

elems2 = Elements(mesh2,2);
faces2 = Faces(model, mesh2);
grid2 = struct('topology',mesh2,'cells',elems2,'faces',faces2);
mat2 = Materials(model,"materialsList.dat");
dofm2 = DoFManager(mesh2,model);

% Create and set the print utility
printUtils1 = OutState(model,mesh1,'outTime.dat','folderName','blockTop','flagMatFile',false);
printUtils2 = OutState(model,mesh2,'outTime.dat','folderName','blockBot','flagMatFile',false);


% Write BC files programmatically
mkdir BCs
foldName = 'BCs';
writeBCfiles(strcat(foldName,'/fixTop'),'NodeBC','Dir',{'Poromechanics','x','y','z'},'FixTop',0,0,mesh1,[2,3,4,5,6]);
writeBCfiles(strcat(foldName,'/fixBot'),'NodeBC','Dir',{'Poromechanics','x','y','z'},'FixBot',0,0,mesh2,[1,3,4,5,6]);

% Collect BC input file in a list
fileName = ["BCs/topLoad.dat","BCs/fixBot.dat"];

%
% Create an object of the "Boundaries" class 
bc1 = Boundaries("BCs/fixTop.dat",model,grid1);
bc2 = Boundaries("BCs/fixBot.dat",model,grid2);

% Create object handling construction of Jacobian and rhs of the model
domain1 = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofm1,...
                     'Boundaries',bc1,...
                     'OutState',printUtils1,...
                     'Materials',mat1,...
                     'Grid',grid1);

domain2 = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofm2,...
                     'Boundaries',bc2,...
                     'OutState',printUtils2,...
                     'Materials',mat2,...
                     'Grid',grid2);

domains = [domain1; domain2];

% build the mortar interface with xml input
[interfaces,domains] = Mortar.buildInterfaces('interfaces.xml',domains);

% get mortar matrices and remove boundary dofs
interfaces{1}.computeMortarMatrices();
% get boundary dof
[bcDof1,~] = getBC(getSolver(domain1,"Poromechanics"),"FixTop",0);
[bcDof2,~] = getBC(getSolver(domain2,"Poromechanics"),"FixBot",0);

% get interface dofs belonging to bc dofs
loc2glob1 = dofId(interfaces{1}.mesh.local2glob{1},3);
loc2glob2 = dofId(interfaces{1}.mesh.local2glob{2},3);

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
getSolver(domain1,"Poromechanics").computeMat([],[]);
getSolver(domain2,"Poromechanics").computeMat([],[]);
K1 = getSolver(domain1,"Poromechanics").K;
K2 = getSolver(domain2,"Poromechanics").K;

K1(:,bcDof1) = [];  K1(bcDof1,:) = [];
K2(:,bcDof2) = [];  K2(bcDof2,:) = [];
K = blkdiag(K2,K1);

S = B'*(K\B);

% test kernel 
v = [1;0;0;-1;zeros(size(S,1)-4,1)];
% 
% [v,e] = eig(S);
% v = real(v); e = real(e);

fprintf("norm(S*v) = %2.4e \n",norm(S*v))









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




