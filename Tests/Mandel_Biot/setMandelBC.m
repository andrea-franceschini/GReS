function fList = setMandelBC(foldName,F,msh)

% Fixed bottom face
f1 = strcat(foldName,'/dirBotFixed');
writeBCfiles(f1,'SurfBC','Dir',{'Poro','z'},'BotFix',0,0,msh,1);
% Top load
f2 = strcat(foldName,'/topLoad');
writeBCfiles(f2,'SurfBC','Neu',{'Poro','z'},'TopLoad',0,F,msh,2);
% Lateral roller normal y
f3 = strcat(foldName,'/dirPoroLatY');
writeBCfiles(f3,'SurfBC','Dir',{'Poro','y'},'LatFixedY',0,0,msh,3);
f4 = strcat(foldName,'/dirPoroLatX');
writeBCfiles(f4,'SurfBC','Dir',{'Poro','x'},'LatFixedX',0,0,msh,4);
% Free flow lateral face
f5 = strcat(foldName,'/dirNoFlow');
writeBCfiles(f5,'SurfBC','Dir',{'Flow'},'LatNoFlow',0,0,msh,5);


fList = string({f1,f2,f3,f4,f5});
for i = 1:numel(fList)
   fList(i) = strcat(fList(i),'.dat');
end

end
