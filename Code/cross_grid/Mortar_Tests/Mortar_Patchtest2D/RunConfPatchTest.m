function [s_x,s_y,u_x] = RunConfPatchTest(Fx,Fy,Dmat)
% Run a conforming 2 block patch test with standard or dual lagrange
% multipliers

% plots solution to paraview
% output: multipliers on the interface

% DEFINE MODEL

% IMPORT MESHES
mesh = Mesh();

% Set the input file name
gaussQuad = Gauss(12,2,2);


% Import the mesh data into the Mesh object
mesh.importGMSHmesh('MeshConforming/Block_hexa.msh');


% Element class for further stiffness matrix computation
elem = Elements(mesh,gaussQuad);

% computing stiffness matrix 
K = stiff(mesh, elem, Dmat, gaussQuad);
f = zeros(2*mesh.nNodes,1);

% ------------------- APPLY BCS -------------------------------

%------------------- TOP LOAD BCS -----------------------------
% get Loaded dofs on top edge
nodesLoad = find(abs(mesh.coordinates(:,2)-2)<1e-3);
[xC,id] = sort(mesh.coordinates(nodesLoad,1),'ascend');
lx = xC(2)-xC(1);
n_ext = nodesLoad(id([1 end]));
n_int = nodesLoad(id(2:end-1));
f(2*n_ext-1) = Fy*lx/2;
f(2*n_int-1) = Fy*lx;

%-------------------------LATERAL LOAD BCS (top only) ---------------------
nodesLoad = find(all([abs(mesh.coordinates(:,1)-0)<1e-3, abs(mesh.coordinates(:,2))>1-1e-3],2));
[yC,id] = sort(mesh.coordinates(nodesLoad,2),'ascend');
ly1 = yC(end)-yC(end-1);
ly2 = yC(2)-yC(1);
% get top node (it has half of the entities influence)
n1 = nodesLoad(id(end));
nIn = nodesLoad(id(3:end-1));
n2 = nodesLoad(id(2));
n3 = nodesLoad(id(1));
%f(2*n1-1) = Fx*ly1/2;
%f(2*nIn-1) = Fx*ly1;
%f(2*n2-1) = Fx*(ly2/2+ly1/2);
% f(2*n3-1) = Fx*ly2/2;

%-------------------------LATERAL LOAD BCS whole lateral side -------------
nodesLoad = find(abs(mesh.coordinates(:,1)-0)<1e-3);
[yC,id] = sort(mesh.coordinates(nodesLoad,2),'ascend');
% get node on the interface
id2 = find(abs(mesh.coordinates(nodesLoad,2)-1)<1e-3);
ly1 = max(yC(2:end)-yC(1:end-1));
ly2 = min(yC(2:end)-yC(1:end-1));
% get top node (it has half of the entities influence)
nOut = nodesLoad(id([1 end]));
%nBot = nodesLoad(id(end-1));
nHalf = [nodesLoad(id2)-1; nodesLoad(id2)];
nIn = nodesLoad(~ismember(nodesLoad,[nHalf;nOut]));
f(2*nOut-1) = Fx*ly1/2;
f(2*nHalf-1) = Fx*(ly1/2+ly2/2);
f(2*nIn-1) = Fx*ly1;

%------------------- BOTTOM FIXED BCS -----------------------------
% get fixed dofs on bottom edge
% y bottom constraint
dirNod = find(abs(mesh.coordinates(:,2)-0)<1e-3);
[K,f] = applyDir(2*dirNod, zeros(length(dirNod),1), K, f);
[K,f] = applyDir(2*dirNod-1, zeros(length(dirNod),1), K, f);

% -------------------LATERAL CONSTRAINT BCS-----------------------
% get nodes on right edge of master domain
nodesFixX = find(abs(mesh.coordinates(:,1)-1)<1e-3);
%[K,f] = applyDir(2*nodesFixX-1, zeros(length(nodesFixX),1), K, f);


% solve linear system
u = K\f;

% compute stresses
stress = computeStress(mesh,elem,Dmat,u,gaussQuad);
%
plotSolution(mesh,'solution_conforming',u,stress);
%
% get stress on slave interface
cellInterf = find(mesh.surfaceTag==2);
[~,centroids] = elem.quad.findAreaAndCentroid(cellInterf');
[~,id] = sort(centroids(:,1),'ascend');

s_x = stress(cellInterf(id),3);
s_y = stress(cellInterf(id),2);

% get id of nodes on the interface
intNod = find(abs(mesh.coordinates(:,2)-1)<1e-3);
[~,id] = sort(mesh.coordinates(intNod,1),'ascend');
u_x = u(2*intNod(id));
end

