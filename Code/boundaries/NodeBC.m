classdef NodeBC < handle
  % NODAL BOUNDARY CONDITIONS class 
  % Constrained node IDs and values of the BCs are defined in the INPUT
  % file.

  properties (Access = public)
    % Indeces of restrained degrees of freedom for the global 
    % stiffness matrix
    boundDof
    % Values of the boundary conditions
    boundVal
    % Type of boundary condition (either Dirichlet or Neumann)
    boundType
    % Physical process for which the BC is imposed
    boundPhysics
    % File ID
    fID
    % Time interval
    timeInt
  end

  methods (Access = public)
    % Class constructor method
    function obj = NodeBC(fid,strVec)
      % Calling the function to set object properties 
      obj.setBC(fid,strVec);
    end
    
%     function updateBC(obj)
%       if obj.timeInt(1,1) < obj.timeInt(1,2)
% %         pos = 1;
%         obj.timeInt(2,1) = 2;
%         obj.timeInt(2,2) = 1;
%       else
% %         pos = 2;
%         obj.timeInt(2,1) = 1;
%         obj.timeInt(2,2) = 2;
%       end
%       readStep(obj,obj.timeInt(2,2));
%     end
  end

  methods (Access = private)
    % Assigning material parameters (read inside the class Boundaries) to object properties
    % Structure of blockSplt:
    %       blockSplt(1)     -> label defining the discrete item on which 
    %                           the BC is applied
    %       blockSplt(2)     -> type of BC (Dir = Dirichlet, Neu = Neumann)
    %       blockSplt(3)     -> physical process for which the BC is
    %                           specified (Poro = poromechanics, 1PFlow = 
    %                           single-phase flow)
    %       blockSplt(4)     -> BC name
    %       blockSplt(5:7)   -> # constrained nodes along x,y,z
    %       blockSplt(8:end) -> node IDs and BC value alternating
    function setBC(obj,fid,blockSplt)
      obj.fID = fid;
      numVec = str2double(blockSplt(5:end))';
      if any(isnan(numVec))
        error('There are nonnumeric entries in %s BC condition',blockSplt(4));
      end
      if strcmp(blockSplt(3),'poro')
        nNod = sum(numVec(1:3));
        l = 3;
      elseif strcmp(blockSplt(3),'flow')
        nNod = numVec(1);
        l = 1;
      else
        error('Unknown physics for boundary %s',blockSplt(4))
      end
      if nNod == 0
        error('No boundary conditions are prescribed for %s BC',blockSplt(4));
      end
      if (length(numVec) ~= 3 + nNod && strcmp(blockSplt(3),'poro')) || ...
         (length(numVec) ~= 1 + nNod && strcmp(blockSplt(3),'flow'))
        error('Wrong number of nodes in BC %s',blockSplt(4));
      end
      obj.boundDof = numVec((l+1):length(numVec));
%       obj.boundVal = numVec((l+2):2:length(numVec));
      obj.boundType = blockSplt(2);
      obj.boundPhysics = blockSplt(3);
      % numVec(1) = # constrained nodes along x
      % numVec(2) = # constrained nodes along y
      % numVec(3) = # constrained nodes along z
      if strcmp(blockSplt(3),'poro')
        s = [ 0; cumsum(numVec(1:3))];
        for i=1:3
          if numVec(i) > 0
          obj.boundDof((s(i)+1):s(i+1)) = ...
                                    (obj.boundDof((s(i)+1):s(i+1))-1)*3 + i;
          end
        end
      end
      obj.boundVal = zeros(length(obj.boundDof),2);
      obj.timeInt = [0 0; 1 2];
%       obj.readStep(1);
%       obj.readStep(2);
    end
    
%     function readStep(obj,pos)
%       try
%       obj.timeInt(1,pos) = fscanf(obj.fID,['TIME' '%e'],[1 1]);
%       catch
%         error('Wrong format of the TIME header in %s %s BC at time %f', ...
%           obj.boundPhysics,obj.boundType,obj.timeInt(1,pos))
%       end
%       try
%         obj.boundVal(:,pos) = fscanf(obj.fID,'%e\n',[1 length(obj.boundDof)]);
%       catch
%         error('Wrong number of values in %s %s BC at time %f', ...
%         obj.boundPhysics,obj.boundType,obj.timeInt(1,pos))
%       end
% %       if ~feof(obj.fID); tmpLine = fgetl(obj.fID); end
%     end
  end
end