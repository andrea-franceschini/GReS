function setBCfiles(mesh,wellsId)
mkdir BCs 
% custom BCs

% fix displacements in outer boundaries
writeBCfiles('BCs/fix_Z','SurfBC','Dir',{'Poromechanics','z'},'fix_Z',0,0,mesh,3); 
writeBCfiles('BCs/fix_Y','SurfBC','Dir',{'Poromechanics','y'},'fix_Y',0,0,mesh,[4 6]); 
writeBCfiles('BCs/fix_X','SurfBC','Dir',{'Poromechanics','x'},'fix_X',0,0,mesh,[5 7]); 
writeBCfiles('BCs/source','VolumeForce','Neu','Poromechanics','wells',[0,3,7,10],[0,0.12,0.44,0.75], wellsId); 

end

