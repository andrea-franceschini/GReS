classdef NodeBC < handle
  % NODAL BOUNDARY CONDITIONS class 
  % Constrained node IDs and the values of the unknowns are defined in the INPUT file.
  % This class returns two vectors, one with the constrained dof(bound_dof) for the stiffness matrix and another one with the value of the unknown for each degree of freedom (bound_value).
    
  properties (Access = public)
    % Vector with the total number of constrained nodes for each 
    % degree of freedom
    numID = [];
    % IDs of constrained nodes
    ID_x = [];
    ID_y = [];
    ID_z = [];
    % Maximum number of degrees of freedom for each node
    dofmax = 3; 
    % Value of the unknown 
    value_x = [];
    value_y = [];
    value_z = [];
        
    % OUTPUT DATA:
    % Indeces of restrained degrees of freedom for the global stiffness matrix, dim = sum(numID)
    bound_dof = [];
    % Values of the unknowns, dim = sum(numID)
    bound_value = [];
  end
    
  methods (Access = public)
     % Class constructor method
     function obj = NodeBC(inputString)
       % Calling the function to set object properties 
       obj.setBCInput(inputString);
     end   

     % Generation of the OUTPUT VECTORS 
     function NodeBoundary(obj,varargin)
       obj.bound_dof = zeros(sum(obj.numID),1);
       i = 1;
       % Generating bound_dof 
       while i <= sum(obj.numID)
          for j = 1 : length(obj.ID_x)
             obj.bound_dof(i) = obj.dofmax*obj.ID_x(j)-2;
             i = i + 1;
          end
          for j = 1 : length(obj.ID_y)
             obj.bound_dof(i) = obj.dofmax*obj.ID_y(j)-1;
             i = i + 1;
          end
          for j = 1 : length(obj.ID_z)
             obj.bound_dof(i) = obj.dofmax*obj.ID_z(j);
             i = i + 1;
          end 
       end
     
       if sum(obj.numID) ~= length(obj.bound_dof)
          error('Error in nodal boundary conditions: wrong bound_dof size');
       end

       % Generating bound_value
       obj.bound_value = zeros(sum(obj.numID),1);
       i = 1;
       while i <= sum(obj.numID)
         for j = 1 : length(obj.value_x)
            obj.bound_value(i) = obj.value_x(j);
            i = i + 1;
         end
         for j = 1 : length(obj.value_y)
            obj.bound_value(i) = obj.value_y(j);
            i = i + 1;
         end
         for j = 1 : length(obj.value_z)
            obj.bound_value(i) = obj.value_z(j);
            i = i + 1;
         end 
       end 
       
     end  
    
  end
           
  methods (Access = private)
      % Assigning material parameters (read inside the class Boundaries) to object properties
     function setBCInput(obj, inputString)
        words = strsplit(inputString, ' ');
        params = zeros(length(words),1);
        k = 0;
        for j = 1 : length(words)
          if (length(words{j}) > 0)
             k = k + 1 ; 
             params(k) = sscanf(words{j}, '%e');
          end
          % Total number of constrained nodes in the mesh
          obj.numID = params(1:obj.dofmax);
        end
      
       constraint = params(obj.dofmax+1:end);
       ID = zeros (max(obj.numID),obj.dofmax);
       value = zeros (max(obj.numID),obj.dofmax);
      
       count = 1;
       for m = 1 : length(obj.numID)
          for n = 1 : obj.numID(m)
             if obj.numID(m) ~= 0
             % Matrix with IDs of constrained nodes (zeros will be ignored)
             ID(n,m) = constraint(count);
             % Matrix with the values of the unknowns (zero not related to a specif ID will be ignored)
             if ID(n,m) == 0
               break
             end
             value(n,m) = constraint(count+1);
             count = count + 2;
             end
          end
       end
      
       % Vectors for constrained dof indeces 
       for c = 1 : obj.dofmax
          [nrow,~] = size(ID);
          for r = 1 : nrow
             if (c==1) && (ID(r,c) ~= 0)
                obj.ID_x(r) = ID(r,c);
             elseif ID(r,c) == 0
                 break
             elseif (c==2) && (ID(r,c) ~= 0)
                obj.ID_y(r) = ID(r,c);
             elseif ID(r,c) == 0
                 break
             elseif (c==3) && (ID(r,c) ~= 0)
                obj.ID_z(r) = ID(r,c);
             elseif ID(r,c) == 0
                 break
             end
          end
       end 
      
       % Value of the constrained dof
       for i = 1:length(obj.ID_x)
         obj.value_x(i) = value(i,1);
       end
       for i = 1:length(obj.ID_y)
         obj.value_y(i) = value(i,2);
       end
       for i = 1:length(obj.ID_z)
         obj.value_z(i) = value(i,3);
       end
         
     end
  end
   
end

     

   

