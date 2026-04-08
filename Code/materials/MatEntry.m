classdef (Abstract) MatEntry < handle
  % MATENTRY Abstract interface for creating materials in GReS
  %
  % This abstract class defines the interface for any material in GReS.
  % Any material must implement all functionality required in this
  % class.
  %
  % Properties:
  %
  %
  % Coupled Solver Notes:
  %
  %
  % Constructor:


  properties (Abstract)
    % the fields modified by the solver
    %fields
  end


  properties (GetAccess=public, SetAccess=protected)
    % handle to domain properties
    data
  end

  methods
    function obj = MatEntry(varargin)
      % inputStruct: struct with additional solver-specific parameters
      obj.data = struct();
    end
  end

  methods (Abstract)

    % mandatory methods that need to be implemented in any material

    % read the input data of the solver and assign variables to cell tags
    materialUpdate(obj,prop,varargin);

    % compute the jacobian and the rhs
    % getProperty(obj,prop,varargin);

  end

  methods (Abstract, Static)

    % get the list of variable fields affected by the solver
    getProps();
    % getMaterial();

  end


  methods

     % function stat = getState(obj,varName)
    %   % get a copy of a state variable field
    %   if nargin < 2
    %     stat = obj.domain.getState();
    %   else
    %     if ~isfield(obj.domain.getState().data,varName)
    %       error("Variable %s does not exist in the State object",varName)
    %     end
    %     stat = obj.domain.getState().data.(varName);
    %   end
    % end

    function out = getMaterial(obj)
      out = class(obj);
    end

    function registerProp(obj,prop,value)
      obj.data.(obj.getMaterial).(prop) = value;
    end

    function out = getProp(obj,prop)
      
      obj.data.(obj.getMaterial).(prop) = value;
    end



  end


  methods (Static)
    

  end

end
