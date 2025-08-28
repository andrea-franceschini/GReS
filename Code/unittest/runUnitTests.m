clear
clc
runtests(fullfile('Mesh','testMesh.m'));
runtests(fullfile('DoFManager','testDoFManager.m'));
% runtests('Material/testMaterialN.m');
runtests('Simparam/testSimparam.m');
runtests('GrowningDomain/testGrowning.m');
