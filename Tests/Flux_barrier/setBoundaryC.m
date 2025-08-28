function fileNames = setBoundaryC(foldName,grid,faces)

    fName = strcat(foldName,'/BC');
    if ~isfolder(fName)
        mkdir(fName);
    end

    fileNames = "";
    % Writing the principal BCs file
    for ii=1:length(faces)
        % Boundary condition.
        if isempty(faces(ii).type) || isempty(faces(ii).values) ...
                || isempty(faces(ii).times) || isempty(faces(ii).field)
            continue
        end
        fileName = strcat(fName,'_',faces(ii).name,".dat");
        fID = fopen(fileName,'w');
        fprintf(fID,'SurfBC                            %% BC item \n');        
        fprintf(fID,'%s                               %% BC type \n',faces(ii).type);
        fprintf(fID,'SinglePhaseFlow                              %% Physics \n');
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