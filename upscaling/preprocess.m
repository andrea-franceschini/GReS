clear
clc

grid = Grid();
grid.importMesh("Input/cylinder.vtk");

% set bc files with variable load value (assumed constant)
setLateralBC(grid,-1);
setTopBC(grid,-1);

% ----------------------------------------------------------------
% ------ function to set confining load on lateral surface -------
% ----------------------------------------------------------------

function setLateralBC(grid,loadValue)

listName = 'BCs/listLat';
valsName = 'BCs/valsLat';

% return files for x,y,z component of lateral load in normal direction
isSurfLateral = grid.surfaces.tag == 3;
nS = sum(isSurfLateral);
surfId = find(isSurfLateral);

load = deal(zeros(nS,3));
k = 0;

% compute normal of each surface element
for id = surfId'
   % retrieve normal
   n = grid.surfaces.normal(id,:);
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

% ----------------------------------------------------------------
% ------ function to set load on top surface ---------------------
% ----------------------------------------------------------------

function setTopBC(grid,loadValue)

% set bc along the x direction at the top surface

listName = 'BCs/listTop';
valsName = 'BCs/valsTop';

% return files for x,y,z component of lateral load in normal direction
isSurfTop = grid.surfaces.tag == 1;
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
