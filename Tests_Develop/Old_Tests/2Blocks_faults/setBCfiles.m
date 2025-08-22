function setBCfiles(leftMesh,rightMesh)
mkdir BCs 
% custom BCs
writeBCfiles('BCs/leftLat','SurfBC','Dir',{'Poromechanics','x'},'left_lateral_roller',0,0,leftMesh,1); % left block lateral fix
writeBCfiles('BCs/leftBot','SurfBC','Dir',{'Poromechanics','z'},'left_bottom_fix',0,0,leftMesh,3); % left block bottom fix
writeBCfiles('BCs/rightBot','SurfBC','Dir',{'Poromechanics','z'},'right_bottom',0,0,rightMesh,4); % right block bottom fix
writeBCfiles('BCs/rightLatLoad','SurfBC','Neu',{'Poromechanics','x'},'right_latLoad',[0 1 6 15 20],[0 -5 -5 0 0],rightMesh,2); % right block bottom fix
writeBCfiles('BCs/rightTopLoad','SurfBC','Neu',{'Poromechanics','z'},'right_topLoad',[0 1 6 15 20],[0 0 -18 -18 0],rightMesh,3); % right block bottom fix
% fix y direction in the symmetry vertical axis
tol = 1e-2;
nL = find(all([abs(leftMesh.coordinates(:,2)-5)<tol,leftMesh.coordinates(:,3)<tol],2));
nR = find(all([abs(rightMesh.coordinates(:,2)-5)<tol,rightMesh.coordinates(:,3)<tol],2));
if isempty(nL) || isempty(nR)
   error('Zero nodes having y=5 coordinates. Mesh should be symmetric!')
end
writeBCfiles('BCs/leftYFix','NodeBC','Dir',{'Poromechanics','y'},'left_yfix',0,0,nL); % right block bottom fix
writeBCfiles('BCs/rightYFix','NodeBC','Dir',{'Poromechanics','y'},'right_yfix',0,0,nR); % right block bottom fix
end

