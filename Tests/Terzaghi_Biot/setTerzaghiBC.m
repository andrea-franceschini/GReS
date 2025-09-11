function setTerzaghiBC(foldName,F,msh)
% F: value of pressure on top of the column
writeBCfiles(strcat(foldName,'/dirFlowTop'),'SurfBC','Dir','SinglePhaseFlow','NoFlowTop',0,0,msh,2);
% Top load
writeBCfiles(strcat(foldName,'/neuPorotop'),'SurfBC','Neu',{'Poromechanics','z'},'TopLoad',0,F,msh,2);
% Lateral roller
writeBCfiles(strcat(foldName,'/dirPoroLatY'),'NodeBC','Dir',{'Poromechanics','y'},'LatFixedY',0,0,msh,3);
writeBCfiles(strcat(foldName,'/dirPoroLatX'),'NodeBC','Dir',{'Poromechanics','x'},'LatFixedX',0,0,msh,4);
% Bottom fixed
writeBCfiles(strcat(foldName,'/dirPoroBottom'),'SurfBC','Dir',{'Poromechanics','x','y','z'},'BotFixed',0,0,msh,1);
% =======
%     % F: value of pressure on top of the column
%     writeBCfiles(strcat(foldName,'/dirFlowTop'),'SurfBC','Dir','Flow','NoFlowTop',0,0,msh,2);
%     % Top load
%     writeBCfiles(strcat(foldName,'/neuPorotop'),'SurfBC','Neu',{'Poro','z'},'TopLoad',0,F,msh,2);
%     % Lateral roller
%     writeBCfiles(strcat(foldName,'/dirPoroLatY'),'NodeBC','Dir',{'Poro','y'},'LatFixedY',0,0,msh,3);
%     writeBCfiles(strcat(foldName,'/dirPoroLatX'),'NodeBC','Dir',{'Poro','x'},'LatFixedX',0,0,msh,4);
%     %  Bottom fixed
%     writeBCfiles(strcat(foldName,'/dirPoroBottom'),'NodeBC','Dir',{'Poro','x','y','z'},'BotFixed',0,0,msh,1);
% >>>>>>> 1dfffa00097f21a2e1d34699913ab58ea5431391
end
