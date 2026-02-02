function setBCfiles(topMesh,botMesh)
mkdir BCs 
% custom BCs

% TOP BLOCK
writeBCfiles('BCs/z_fixed_top','SurfBC','Dir',{'Poromechanics','z'},'z_fixed_top',0,0,topMesh,2); % left block lateral fix
writeBCfiles('BCs/topLoad','SurfBC','Dir',{'Poromechanics','y'},'left_bottom_fix',0,-0.1,topMesh,3); % left block bottom fix


% BOTTOM BLOCk

% find nodes on bottom right corner to fix in x direction
tol = 1e-3;
nList = find(all([abs(botMesh.coordinates(:,1)-2)<tol,botMesh.coordinates(:,2)<tol],2));

writeBCfiles('BCs/z_fixed_bottom','SurfBC','Dir',{'Poromechanics','z'},'z_fixed_bottom',0,0,botMesh,2); % right block bottom fix
writeBCfiles('BCs/y_fixed','SurfBC','Dir',{'Poromechanics','y'},'y_fixed_bottom',0,0,botMesh,1); % right block bottom fix
writeBCfiles('BCs/x_fixed','NodeBC','Dir',{'Poromechanics','x'},'x_fixed',0,0,nList); % right block bottom fix
end

