function setBCfiles(mesh)
mkdir BCs 
% custom BCs

% fix Z - plain strain condition
writeBCfiles('BCs/z_fix','SurfBC','Dir',{'Poromechanics','z'},'z_fix_out',0,0,mesh,3); 


% neumann condition for constant compressive stress
writeBCfiles('BCs/loadLeft','SurfBC','Neu',{'Poromechanics','x'},'left_load',0,100,mesh,1); 
writeBCfiles('BCs/loadRight','SurfBC','Neu',{'Poromechanics','x'},'right_load',0,-100,mesh,2); 


% fix X and Y at centered node locations in the outer boundary
tol = 1e-3;
C = max(mesh.coordinates,[],"all");
nY = find(all([abs(mesh.coordinates(:,1)-C)<tol,...
               abs(mesh.coordinates(:,2))<tol],2));

nY = [nY; find(all([abs(mesh.coordinates(:,1)+C)<tol,...
               abs(mesh.coordinates(:,2))<tol],2))];

nX = find(all([abs(mesh.coordinates(:,1)) <tol,...
               abs(mesh.coordinates(:,2) - C) < tol],2));


nX = [nX;find(all([abs(mesh.coordinates(:,1)) <tol,...
               abs(mesh.coordinates(:,2) + C) < tol],2))];

if isempty(nX) || isempty(nY)
   error(['Check for odd number of nodes in the outer boundary to' ...
     'correctly fix X and Y directions!']);
end

writeBCfiles('BCs/x_fix','NodeBC','Dir',{'Poromechanics','x'},'x_fix_out',0,0,nX); 
writeBCfiles('BCs/y_fix','NodeBC','Dir',{'Poromechanics','y'},'y_fix_out',0,0,nY); 





end