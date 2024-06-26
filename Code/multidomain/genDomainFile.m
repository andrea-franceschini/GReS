function genDomainFile(name,nD)
% This utility writes an empty domain input file in the current folder
fID = fopen(name,'w');
fprintf(fID,'%% Input for multidomain simulations \n');
for i = 1:nD
   fprintf(fID,'%% Domain %i \n',i);
   fprintf(fID,'<Domain> \n');
   fprintf(fID,'<Name> \n');
   fprintf(fID,'<ModelType> \n');
   fprintf(fID,'<SimulationParameters> \n');
   fprintf(fID,'<Geometry> \n');
   fprintf(fID,'<Gauss> \n');
   fprintf(fID,'<Material> \n');
   fprintf(fID,'<DoFManager> \n');
   fprintf(fID,'<BoundaryConditions> \n');
   fprintf(fID,'<OutState> \n');
   fprintf(fID,'<EndDomain> \n \n');
end
end
