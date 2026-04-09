% test new grid utilities

grid = Grid();
grid.importMesh('Column_tetra.msh');
% 
% cells = grid.cells;
% faces = grid.faces;
% 
% [fList,ptr] = getData(cells.cells2faces);
% nFPC = diff(ptr);                       % number of faces per cell
% 
% hf2Cell = repelem((1:cells.num)',nFPC,1);
% hf2Cell == faces.neighbors(fList);

%%
% clear
% clc
%
% profile off
% profile on
% grid = structuredMesh(100,100,10,[0 1],[0 1],[0 1]);
% nC = grid.cells.num;
% profile viewer