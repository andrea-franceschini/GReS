function [K,areaNod] = stiffPoisson(mesh, element)
% computing classic stiffness matrix for Poisson problem
% uses element class from GReS to compute derivatives of basis functions
areaNod = zeros(mesh.nNodes,1);
nEntries = sum(mesh.surfaceNumVerts.^2);
iiVec = zeros(nEntries,1);
jjVec = zeros(nEntries,1);
KVec = zeros(nEntries,1);
%
l = 0;
for el = 1:mesh.nSurfaces
  % Get the right material stiffness for each element
  vtkId = mesh.surfaceVTKType(el);
  nodes = mesh.surfaces(el,:);
  if vtkId == 5
    elem = element.getElement(vtkId);
    dN = getDerBasisF(elem,el);
    [A,~] = findAreaAndCentroid(element.getElement(vtkId),el);
    Kloc = dN'*dN*A;
  else
    elem = getElement(element,vtkId);
    nCoord = mesh.coordinates(nodes,1:2);
    [gradN,dJw] = getDerBasisFAndDet(elem,nCoord);
    Kloc = pagemtimes(gradN,'transpose',gradN,'none');
    Kloc = sum(Kloc.*reshape(dJw,1,1,[]),3);
  end
  s = numel(Kloc);
  % get global DoF
  % node indices
  % update vector of nodal areas
  areaNod(nodes) = areaNod(nodes) + element.findNodeArea(el);
  % dof corresponding to dofs
  [jjLoc,iiLoc] = meshgrid(nodes,nodes);
  iiVec(l+1:l+s) = iiLoc(:);
  jjVec(l+1:l+s) = jjLoc(:);
  KVec(l+1:l+s) = Kloc(:);
  l = l + s;
end
% populate stiffness matrix in sparse format
K = sparse(iiVec, jjVec, KVec);
end
