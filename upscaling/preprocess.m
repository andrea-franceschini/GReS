clear
clc

mesh = Mesh();
mesh.importMesh("Input/cylinder.vtk");
topology = Mesh();

% set bc files with variable load value (assumed constant)
setLateralBC(mesh,-1);
setTopBC(mesh,-1);


function setLateralBC(mesh,loadValue)

listName = 'BCs/listLat';
valsName = 'BCs/valsLat';

% return files for x,y,z component of lateral load in normal direction
isSurfLateral = mesh.surfaceTag == 3;
nS = sum(isSurfLateral);
surfId = find(isSurfLateral);

tri = Triangle(1,mesh);

load = deal(zeros(nS,3));
k = 0;

% compute normal of each surface element
for id = surfId'
  % compute normalized normal
  n = computeNormal(tri,id);
  load(k+1,:) = loadValue*n;
  k = k+1;
end

vals = load(:);


% write to file
fList = fopen(listName,'w');
fprintf(fList,'%i ',[nS nS nS]);
fprintf(fList,'   %% Number of fixed entities \n');
list = repmat(surfId,3,1);
fprintf(fList,'%i \n',list);

% Writing BC vals for each time step
ft = fopen(valsName,'w');
fprintf(ft,'%%Time %2.4f \n',0.0);
fprintf(ft,'%1.6e \n',vals);

end



function setTopBC(mesh,loadValue)

% set bc along the x direction at the top surface

listName = 'BCs/listTop';
valsName = 'BCs/valsTop';

% return files for x,y,z component of lateral load in normal direction
isSurfTop = mesh.surfaceTag == 1;
nS = sum(isSurfTop);
surfId = find(isSurfTop);


load = loadValue*ones(nS,1);

vals = load(:);


% write to file
fList = fopen(listName,'w');
fprintf(fList,'%i ',[nS 0 0]);
fprintf(fList,'   %% Number of fixed entities \n');
fprintf(fList,'%i \n',surfId);

% Writing BC vals for each time step
ft = fopen(valsName,'w');
fprintf(ft,'%%Time %2.4f \n',0.0);
fprintf(ft,'%1.6e \n',vals);

end
