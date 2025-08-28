clear
close all
clc

%%
runFluidProblem;
press = solverFlow.getState().data.pressure;
clearvars -except press solverFlow gridFluid

%%
fluid2grains;

clearvars -except press_surf press_grain mshGrain press_grain_inner


plotFunction(mshGrain,'Output_InterpolatedPressure',press_grain_inner)



%%
runGrainsProblem;

%%
% list = load('listNan.dat');
% 
% nodes = mapDofToNodeID(list,solver);
% 
% 
% 
% function nodeID = mapDofToNodeID(dof,solver)
% % run this if we obtain a solution
% 
% nInterf = numel(solver.interfaces);
% nDofsInterf = zeros(nInterf,1);
% nodeID = zeros(3*numel(dof),1);
% % map global dof id to cell in the pore grid
% for i = 1:numel(solver.interfaces)
%   nDofsInterf(i) = solver.interfaces{i}.totMult;
% end
% 
% nDofCum = cumsum(nDofsInterf);
% nDofCum2 = [0;nDofCum(1:end-1)];
% dof = dof - solver.nDof;
% 
% for d = 1:numel(dof)
%   id = dof(d) > nDofCum;
%   interfID = sum(id) + 1;
%   multID = dof(d) - nDofCum2(interfID);
%   surfID = ceil(multID/3);
%   nodes = solver.interfaces{interfID}.mesh.msh(2).surfaces(surfID,:);
%   nodeID(3*d-2:3*d) = solver.interfaces{interfID}.mesh.local2glob{2}(nodes);
% end
% 
% nodeID = unique(nodeID);
% 
% 
% x = nodeID';   % example column
% str = sprintf('%i, ', x-1);    % makes one long string
% str = str(1:end-1);           % remove trailing comma
% disp(str)
% end
% 
