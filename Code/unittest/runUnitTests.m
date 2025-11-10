clear
clc
runtests(fullfile('Mesh','testMesh.m'));
runtests(fullfile('DoFManager','testDoFManager.m'));
runtests('Simparam/testSimparam.m');
runtests('Materials/testMaterials.m');
