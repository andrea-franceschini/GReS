clear
clc

runtests('Simparam/testSimparam.m');
runtests(fullfile('Mesh','testMesh.m'));
runtests(fullfile('DoFManager','testDoFManager.m'));

% runtests('Material/testMaterialN.m');
runtests('GrowningDomain/testGrowning.m');