classdef Dirichlet < handle
    % DIRICHLET boundary conditions for fixing the value of the unknowns in a
    % certain region of the domain
    % Entities IDs and the values of the unknowns are defined in the INPUT file
    % This class returns two matrices, one with the constrained indeces
    % for the stiffness matrix and another one with the value of the
    % unknown for each degree of freedom, for every entity of the mesh that 
    % is part of the Dirichlet's boundary
    
    properties 
        % INPUT DATA:
        % Vector with the total number of constrained entities for each 
        % degree of freedom
        numID = [];
        % ID of constrained entity
        ID = [];
        % Maximum number of degrees of freedom
        dofmax = []; 
        % Value of the unknown (dim = ID*dofmax)
        value = [];
        
        % OUTPUT DATA:
        % Vector with indeces of restrained nodes for the system stiffness 
        % matrix (dim = numID*dofmax)
        dir_ind = [];
        % Vector with fixed value of the unknown (dim = numID*dofmax)
        dir_value = [];
    end
    
    methods (Access = public)
     % Class constructor method
     function obj = Dirichlet(inputString)
       % Calling the function setDirInput to set input properties
       obj.setDirInput(inputString);
     end
   
    
     % Function for the generation of the output vectors
     function DirBoundary(obj,varargin)
        % Memory allocation
        obj.dir_ind = zeros(sum(obj.numID),1);
        i = 1;
        % Generating dir_ind 
        while i <= sum(obj.numID)
        for col = 1 : obj.dofmax
            for row = 1 : max(obj.numID)
                if obj.ID(row,col) ~= 0
                obj.dir_ind(i) = obj.dofmax * (obj.ID(row,col)-1) + col;
                i = i + 1;
                end
            end
        end
        end

        
        obj.dir_value = zeros(sum(obj.numID),1);
        i = 1;
        % Generating dir_value
        while i <= sum(obj.numID)
        for col = 1 : obj.dofmax
            for row = 1 : max(obj.numID)
                if obj.value(row,col) ~= 0
                obj.dir_value(i) = obj.value(row,col);
                end
                i = i + 1;
            end
        end
        end

                
    end  
   end
   
          

   methods (Access = private)
      % Function that sets the parameters coming from "data"
      % (Boundaries) inside the vector "params"
      function setDirInput(obj, inputString)
      words = strsplit(inputString, ' ');
      params = zeros(length(words),1);
      k = 0;
      for j = 1 : length(words)
        if (length(words{j}) > 0)
            k = k + 1 ; 
            params(k) = sscanf(words{j}, '%e');
        end
            % Number of degrees of freedom
            obj.dofmax = params (1);
            % Total number of constrained entinties in the mesh
            obj.numID = params(2:obj.dofmax+1);
      end
      
      constraint = params(obj.dofmax+2:end);
      % Memory allocation
      obj.ID = zeros (max(obj.numID),obj.dofmax);
      obj.value = zeros (max(obj.numID),obj.dofmax);
      
      
      count = 1;
      for m = 1 : length(obj.numID)
          for n = 1 : obj.numID(m)
             if obj.numID(m) ~= 0
             % Matrix with IDs of constrained entities (zeros will be ignored)
             obj.ID(n,m) = constraint(count);
             % Matrix with the value of the unknown (zeros not related to a specif
             % ID will be ignored)
             obj.value(n,m) = constraint(count+1);
             count = count + 2;
             end
           end
      end

      
      end
    

   end

end

     

   

