function writeBCfile_sigma_nn(fName, bcName, time, sigma_n_val, ...
    mesh, surfaceTag)
    % Not trying to write a very general functions, which implies:
    % I fix type == Neu, item == 'SurfBC', physic = {'Poro','x','y','z'}
    % Also treating time as 0 and sigma_n_val as a scalar (same everywhere
    % since we are modeling it for c_max=1 bc here)
    
    type = 'Neu';
    item = 'SurfBC';
    physics = {'Poromechanics', 'x', 'y', 'z'};
    physics = string(physics);
    % Extracting physics from the coordinates
    ph = physics(1);
    dir = physics(2:end);
    
    assert(length(time)==length(sigma_n_val),['BC times and values ' ...
        'sets must have equal size']);
    
    % Writing general BCs file
    fID = fopen(strcat(fName,'.dat'),'w');
    
    % Start writing the data file - space given in the middle for clarity in
    % the data file
    fprintf(fID,'%s            %% BC item \n',item);
    fprintf(fID,'%s            %% BC type \n',type);
    fprintf(fID,'%s            %% Physics \n',ph);
    fprintf(fID,'%s            %% BC name \n',bcName);
    
    % Encode the list file name
    listName = strcat(fName,'/list'); % Create a list in a folder with fName
    fprintf(fID,'%s \n',listName); % Encode it in the bc data file
    
    % Encode the time files
    for i = 0:length(time)-1
      fprintf(fID,'%2.6f %s/time%i.dat \n',time(i+1),fName,i);
    end
    % End file
    fprintf(fID,'End');
    
    % Make a folder with fName if it doesn't exist yet
    if ~isfolder(fName)
      mkdir(fName);
    end
    
    % % Enable writing in the fName/list.dat file
    fList = fopen(listName,'w');
    
    % Extract the list of all surfaces to constrain
    surf = find(ismember(mesh.surfaceTag,surfaceTag));
    
    % For 'SurfBC': list = surf for better notation to write the /list file
    list = surf;
    
    % logical for distinguishing input directions
    which_directions = ismember(["x","y","z"], dir);
    n_directions = sum(which_directions);
    fprintf(fList,'%i ', which_directions*length(list));
    % Adding a comment in the output file
    fprintf(fList,'   %% Number of fixed entities \n');
    % Repeat the list for the number of coordinates specified in the input
    list = repmat(list, n_directions, 1);
    % Print the repeated list in the output file
    fprintf(fList,'%i \n',list);
    
    % Enable writing in the fName/time0.dat file
    t_name = strcat(fName,'/time0.dat');
    ft = fopen(t_name,'w');
    fprintf(ft,'%%Time %2.4f \n', 0); % Print time (here 0) with 4 decimals
    
    % Compute the x-,y-,z-direction sigmas to impose on each surfaces to write
    % in the fName/time0.dat file
    N_surfaces = length(surf);
    sigma_n_final = zeros(n_directions*N_surfaces, 1);
    for i = 1:N_surfaces
        surfacePoints = mesh.surfaces(surf(i),:); % Extract point indices for the i-th surface
        sigma_n_x = 0;
        sigma_n_y = 0;
        sigma_n_z = 0;
        for j = 1:length(surfacePoints)
            % Extract coordinates for j-th point in the i-th surface
            coords = mesh.coordinates(surfacePoints(j),:);
            
            % Compute unit normals and the corresponding stresses to be imposed
            unit_normals = coords / norm(coords); % unit normals
            sigma_n_x = sigma_n_x + sigma_n_val*unit_normals(1);
            sigma_n_y = sigma_n_y + sigma_n_val*unit_normals(2);
            sigma_n_z = sigma_n_z + sigma_n_val*unit_normals(3);
        end

        % Average sigma_n over all points in a surface
        sigma_n_x = sigma_n_x / length(surfacePoints);
        sigma_n_y = sigma_n_y / length(surfacePoints);
        sigma_n_z = sigma_n_z / length(surfacePoints);


        % Store the stresses in sigma_n_final
        sigma_n_final(i)                = sigma_n_x;
        sigma_n_final(N_surfaces + i)   = sigma_n_y;
        sigma_n_final(2*N_surfaces + i) = sigma_n_z;
    end
    
    % Write sigma_n_final in the fName/time0.dat file
    fprintf(ft,'%1.6e \n', sigma_n_final); % print the time value list
end
