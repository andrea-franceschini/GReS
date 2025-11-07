function setBCfiles(meshList)
mkdir BCs 

assert(numel(meshList)==8)
% custom BCs
for i = 1:4
  writeBCfiles("BCs/block_"+num2str(i)+"_fix",'SurfBC','Dir',{'Poromechanics','x','y','z'},'fix',0,0,meshList{i},1)
end

for i = 5:8
  writeBCfiles("BCs/block_"+num2str(i)+"_load",'SurfBC','Neu',{'Poromechanics','z'},'load',0,-1,meshList{i},6)
end

end

