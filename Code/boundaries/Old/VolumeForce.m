classdef VolumeForce < handle
  % Class dealing with external forces in the model () 
  % Constrained node IDs and values of the BCs are defined in the INPUT
  % file.

  properties (Access = public)
    % Indeces of restrained degrees of freedom for the global 
    % stiffness matrix
    boundDof
    % Values of the boundary conditions
    boundVal
    % Type of boundary condition (either Dirichlet or Neumann)
%     boundType
    % Physical process for which the BC is imposed
    boundPhysics
    % File ID
    fID
    % Time interval
    timeInt
  end

  methods (Access = public)
    % Class constructor method
    function obj = VolumeForce(fid,strVec)
      % Calling the function to set object properties 
      obj.setVF(fid,strVec);
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
    %       blockSplt(2)     -> physical process for which the BC is
    %                           specified (Poro = poromechanics, 1PFlow = 
    %                           single-phase flow)
    %       blockSplt(3)     -> BC name
    %       blockSplt(4:7)   -> # constrained nodes along x,y,z
    %       blockSplt(8:end) -> node IDs and BC value alternating
    function setVF(obj,fid,blockSplt)
      obj.fID = fid;
      numVec = str2double(blockSplt(4:end))';
      if any(isnan(numVec))
        error('There are nonnumeric entries in %s volume force',blockSplt(3));
      end
      nNod = numVec(1);
      if nNod == 0
        error('No loaded cells for %s volumetric force',blockSplt(3));
      end
      if length(numVec) ~= 1 + nNod
        error('Wrong number of loaded cells in %s volumetric force',blockSplt(3));
      end
      obj.boundDof = numVec(2:end);
%       obj.boundVal = numVec((l+2):2:length(numVec));
      obj.boundPhysics = blockSplt(2);
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