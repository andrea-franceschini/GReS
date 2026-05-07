function bcs = setBCfiles(grid)
mkdir BCs

assert(numel(grid)==8)

bcs = cell(8,1);
% custom BCs
for i = 1:4
  bcs{i} = Boundaries('BCs/fix.xml',grid{i});
end

for i = 5:8
  bcs{i} = Boundaries('BCs/load.xml',grid{i});
end

end

