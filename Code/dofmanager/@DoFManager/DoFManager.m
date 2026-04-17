classdef DoFManager < handle
  % Degree of freedom manager
  % Map each entity to a corresponding degree of freedom in the solution
  % system
  % 2 ordering can be specified in input:
  % Field-based ordering (default)
  % Domain-based ordering
  properties (Access = private)
    grid
    dofMap                        % cell array with dof map for each variable field
    numbComponents
    fields = struct("variableName",[],...
                    "range",[],...
                    "tags",[],...
                    "fieldLocation",[]);
    nVars = 0;
    totDofs
  end


  methods (Access = public)
    function obj = DoFManager(grid)

      if nargin == 0 
        return
      end
      obj.grid = grid;
      obj.totDofs = 0;
    end

    registerVariable(obj,varName,fldLoc,nComp,varargin)


    function dofs = getLocalDoF(obj,id,ents)

      % return the DoF numbering in global indexing for input entities
      % in isempty(varargin) all dofs are returned in correct order
      % the input id is the result of a call to getVariableId(varName)

      if nargin < 3
        % return all dofs of id field
        ents = find(obj.dofMap{id} > 0);
      end

      if isscalar(id)
        ncomp = obj.numbComponents(id);
        dofs = obj.dofMap{id}(ents);
        dofs =  repelem(dofs,ncomp,1) + ...
          repmat((0:ncomp-1)',numel(ents),1);
      else
        for i = 1:numel(id)
          dofs = [dofs; getLocalDoF(id(i))];
        end
      end
    end


    function ents = getLocalEnts(obj,id,ents)

      % same as getLocalDoF but WITHOUT component expansion
      ncomp = obj.numbComponents(id);
      ents = round(((ncomp-1)+obj.dofMap{id}(ents))/ncomp);
    end

    function dofs = getDoF(obj,id,ents)
      % return the DoF numbering in global indexing for input entities
      % if isempty(varargin) all dofs are returned in correct order
      % the id is the result of a call to getVariableId(varName)

      if isscalar(id)
        dofs = getLocalDoF(obj,id) + obj.fields(id).range(1) - 1;
      else
        assert(isempty(ents),"entity input in getDoF is admitted " + ...
          "only for single variable query")

        for i = 1:numel(id)
          dofs = [dofs; getDoF(obj,id(i))];
        end
      end

    end


    function activeEnts = getActiveEntities(obj,varId,flagExpand)
      
      if ~isnumeric(varId)
        varId = string(varId);
      end

      assert(isscalar(varId),"Input variable must be a scalar string or a" + ...
        "character vector (or an integer)")
      id = obj.getVariableId(varId);
      activeEnts = find(obj.dofMap{id});

      % expand to account for component
      if nargin > 2
        if flagExpand
          activeEnts = obj.dofExpand(activeEnts,obj.numbComponents(id));
        end
      end
    end

    function tags = getTargetRegions(obj,varId)
      % get regions where variable fields are present simultaneously
      % varId - array of id or a string array of variable names
      % if the variable is not present, the result is empty

      if nargin > 2 
        error("getTargetRegions:List of variables must be a single array or  single string array of variable names")
      end

      tags = [];


      if numel(varId) > 1
        tags = obj.getTargetRegions(varId(1));
        for i = 2:numel(varId)
          tags = intersect(tags, obj.getTargetRegions(varId(i)));
        end
        return
      end

      if isVariable(obj,varId)
        id = obj.getVariableId(varId);
        tags = unique([obj.fields(id).tags]);
      end

    end

    function cells = getFieldCells(obj,varId)
      tags = getTargetRegions(obj,varId);
      cells = getEntitiesFromTags(...
        entityField.cell,obj.grid,entityField.cell,tags);
    end

    function location = getFieldLocation(obj,varId)

      id = obj.getVariableId(varId);
      location = obj.fields(id).fieldLocation;


    end

    function id = getVariableId(obj,varId)
      % return the id of the requested input variable

      assert(isVariable(obj,varId),"Requested variable is not available" + ...
          "in the DoFManager")

      if isnumeric(varId)
        id = varId;
      else     
        id = find(strcmp([obj.fields.variableName],varId));
      end
    end

    function out = isVariable(obj,varId)

      if isnumeric(varId)
        out = all([varId(:) > 0; varId(:) <= numel(obj.fields)]);
      else
        out = any(contains([obj.fields.variableName],varId));
      end


    end

    function numDof = getNumbDoF(obj,varId)
      if nargin == 1
        numDof = obj.totDofs;
      else
        id = getVariableId(obj,varId);
        r = obj.fields(id).range;
        numDof = r(2)-r(1)+1;
      end
    end

    function nVars = getNumberOfVariables(obj)
      nVars = obj.nVars;
    end

    function varNames = getVariableNames(obj,varId)
      if nargin == 1
         varNames = [obj.fields.variableName];
      elseif nargin == 2
         id = getVariableId(obj,varId);
         varNames = [obj.fields.variableName];
         varNames = varNames(id);
      end
    end

    function ncomp = getNumberOfComponents(obj,varId)
      id = getVariableId(obj,varId);
      ncomp = obj.numbComponents(id);
    end


  end

  methods (Static)

    function dofOut = dofExpand(dofIn,nComp)
      % make sure input is a column vector
      dofIn = reshape(dofIn,[],1);
      dofOut = nComp*repelem(dofIn,nComp,1)-repmat((nComp-1:-1:0)',size(dofIn,1),1);
    end
  end
end
