function fileNames = setBC(foldName,grid,faces)
%SETBC - create the files necessary to aply the boundary conditions using
%the regions defined in PhysicalNames of the gmsh file.
%
% To use the Physical Surface defined in gmsh to impose a boundary
% condition, it's necessary to pass the input variable faces as a struct
% faces = struct('name',[],'type',[],'physics',[],'field',[],'values',[],'times',[]);
% - faces(n) - number of boundaries condition to be created.
% - faces.name - name associated to this boundary condition saved in the folder.
% - faces.type - type of boundary condition, 'Dir'-Dirichlet,'Neu'-Neumann.
% - faces.physics - in which physics field the boundary is applied.
% - faces.values - a vector with the values of each time.
% - faces.times - a vector with the times the values change.
% the vectors values and times, need to be the same dimension for each
% boundary condition.

fName = strcat(foldName,'Boundary');
if ~isfolder(fName)
    mkdir(fName);
end

fileNames = "";
% Writing the principal BCs file
for ii=1:length(faces)
    % Boundary condition.
    if isempty(faces(ii).name) || ...
            isempty(faces(ii).type) || ...
            isempty(faces(ii).physics) || ...
            isempty(faces(ii).field) || ...
            isempty(faces(ii).values) || ...
            isempty(faces(ii).times)
        continue
    end
    fileName = strcat(foldName,'Bound_',faces(ii).name,".dat");
    fID = fopen(fileName,'w');
    fprintf(fID,'SurfBC                            %% BC item \n');
    fprintf(fID,'%s                               %% BC type \n',faces(ii).type);
    fprintf(fID,'%s                              %% Physics \n',faces(ii).physics);
    fprintf(fID,'%s                           %% BC name \n',faces(ii).name);
    item = strcat(fName,'/List_',faces(ii).name);
    fprintf(fID,'%s \n',item);
    ntimes = length(faces(ii).times);
    for i=1:ntimes
        step = faces(ii).times(i);
        item = strcat(fName,'/Time_',faces(ii).name,'_',int2str(i-1),'.dat');
        fprintf(fID,'%s %s \n',num2str(step,'% 10.2f'),item);
    end
    fprintf(fID,'End');
    fclose(fID);
    fileNames = fileNames+fileName+' _||_ ';

    % Writing list of nodes
    elms = [];
    for i=1:length(faces(ii).field)
        temp = find(grid.topology.surfaceTag==grid.topology.surfaceRegions.(faces(ii).field(i)));
        elms = [elms ; temp];
    end
    elms=sort(elms);
    fID = fopen(strcat(fName,'/List_',faces(ii).name),'w');
    fprintf(fID,'%i                   %% #ID_constrained',length(elms));
    for i=1:length(elms)
        fprintf(fID,'\n%i',elms(i));
    end
    fclose(fID);

    % Writing each time step
    for i=1:ntimes
        item = strcat(fName,'/Time_',faces(ii).name,'_',int2str(i-1),'.dat');
        fID = fopen(item,'w');
        fprintf(fID,'%% TIME %s',num2str(faces(ii).times(i),'% 10.2f'));
        for j=1:length(elms)
            fprintf(fID,'\n%d',faces(ii).values(i));
        end
        fclose(fID);
    end
end
fileNames = strsplit(fileNames,' _||_ ');
fileNames = fileNames(1:length(fileNames)-1);
end