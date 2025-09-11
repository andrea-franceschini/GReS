function processBC_grains(grainMesh,interfaces,press)

nodesInterf = [];

for i = 1:numel(interfaces)
  nodesInterf = [nodesInterf;(interfaces{i}.mesh.local2glob{1})';(interfaces{i}.mesh.local2glob{2})'];
end

nodesInterf = unique(nodesInterf);

isContactSurf = ~any(~ismember(grainMesh.surfaces,nodesInterf),2);

% get list of surfaces to apply fluid pressure:
% exclude boundary surfaces and grain2grain surfaces

tol = 1e-3;
% get surfaces having only boundary nodes
isBoundSurf = any([abs(grainMesh.surfaceCentroid-0.15)<tol,abs(grainMesh.surfaceCentroid-0.85)<tol],2);

nodesTop = find(abs(grainMesh.coordinates(:,3)-0.85)<tol);
nodesBot = find(abs(grainMesh.coordinates(:,3)-0.15)<tol);
nodesNorth = find(abs(grainMesh.coordinates(:,2)-0.85)<tol);
nodesSouth = find(abs(grainMesh.coordinates(:,2)-0.15)<tol);
nodesEast = find(abs(grainMesh.coordinates(:,1)-0.15)<tol);
nodesWest = find(abs(grainMesh.coordinates(:,1)-0.85)<tol);
nodesBound = unique([nodesWest;nodesEast;nodesSouth;nodesNorth;...
  nodesBot;nodesTop]);


isFluidSurf = ~any([isBoundSurf isContactSurf],2);
% get id of nodes on the fluid surfaces
nodeFluid = unique(grainMesh.surfaces(isFluidSurf,:));
nodeFluid = nodeFluid(~ismember(nodeFluid,nodesBound));

% get nodal forces (pressure*area*averagenodenormal)
nodeNormal = zeros(grainMesh.nNodes,3);
nodeArea = zeros(grainMesh.nNodes,1);
for i = reshape(find(isFluidSurf),1,[])
  nid = grainMesh.surfaces(i,:);
  nodeArea(nid) = nodeArea(nid)+grainMesh.surfaceArea(i)/3;
  coords = grainMesh.coordinates(nid,:);
  v1 = coords(1,:)-coords(2,:);
  v2 = coords(1,:)-coords(3,:);
  n = cross(v1,v2)/norm(cross(v1,v2));
  nodeNormal(nid,:) = nodeArea(nid).*n;
end

%nodeArea = nodeArea(nodeFluid);

surfTags = [];
% spot boundary grains lying on the boundary (here we must fix all nodes)
% fix all associated dofs
% get all effective surfaceTag (included in the input of the interfaces)
for i = 1:numel(interfaces)
  surftaglist = [interfaces{i}.mesh.msh(1).surfaceTag;...
    interfaces{i}.mesh.msh(2).surfaceTag];
  surfTags = [surfTags;unique(surftaglist)];
end

surfTags = unique(surfTags);
floatTags = ~ismember(1:grainMesh.nSurfaceTag,surfTags);

nodeFloating = unique(grainMesh.surfaces(ismember(grainMesh.surfaceTag,[find(floatTags) 138]),:));

% remove fluid pressure to fixed nodes
[~,ia,~] = intersect(nodeFluid,nodeFloating);
nodeFluid(ia) = [];

nodeNormal = nodeNormal(nodeFluid,:);   % already considering the node area!

press_nodes = press(nodeFluid).*nodeNormal;
pX = press_nodes(:,1);
pY = press_nodes(:,2);
pZ = press_nodes(:,3);




% modify bcs of domain structure





writeBCfiles('BCs/fixZ','NodeBC','Dir',{'Poromechanics','z'},'FixedZ',0,0,[nodesTop;nodesBot]);
writeBCfiles('BCs/fixX','NodeBC','Dir',{'Poromechanics','x'},'FixedX',0,0,[nodesWest;nodesEast]);
writeBCfiles('BCs/fixY','NodeBC','Dir',{'Poromechanics','y'},'FixedY',0,0,[nodesSouth;nodesNorth]);
writeBCfiles('BCs/fixAll','NodeBC','Dir',{'Poromechanics','x','y','z'},'FixedAll',0,0,nodeFloating);
writeBCfiles('BCs/pressureX','NodeBC','Neu',{'Poromechanics','x'},'LoadX',0,0,nodeFluid,-pX);
writeBCfiles('BCs/pressureY','NodeBC','Neu',{'Poromechanics','y'},'LoadY',0,0,nodeFluid,-pY);
writeBCfiles('BCs/pressureZ','NodeBC','Neu',{'Poromechanics','z'},'LoadZ',0,0,nodeFluid,-pZ);

end

