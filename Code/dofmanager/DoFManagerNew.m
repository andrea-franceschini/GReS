classdef DoFManagerNew < handle
  % Degree of freedom manager
  % Map each entity to a corresponding degree of freedom in the solution
  % system
  % 2 ordering can be specified in input:
  % Field-based ordering (default)
  % Domain-based ordering
  properties (Access = private)
    grid
    dofMap % cell array with dof map for each variable field
    numbComponents
    fields = struct("variableName",[],...
      "range",[],...
      "tags",[]);
    nVars = 0;
    totDofs
  end


  methods (Access = public)
    function obj = DoFManagerNew(grid)
      obj.grid = grid;
      obj.totDofs = 0;
    end

    function registerVariable(obj,varName,fieldLocation,nComp,tags)
      % varName: the name of the variable field
      % fieldLocation: a enum of type entityField
      % tags: the cellTag (or surfaceTag for lower dimensional fields) where the variable is actually present
      % numbComponents: the number of dofs per entitiy

      id = obj.nVars+1;


      if any(strcmpi([obj.fields.variableName],varName))
        % variable field already exist
        id = getVariableId(obj,varName);
        assert(isempty(intersect(tags,obj.fields(id).tags)), ...
          "Variable field %i has been already defined on specified tags", ...
          tags);

        % update the tags
        obj.fields(id).tags = sort([obj.fields(id).tags tags]);

        % update the entities with only new entities
        entList = getEntities(fieldLocation,obj.grid,tags);
        isInactive = obj.dofMap{id}(entList) == 0;
        nNewEnts = sum(isInactive);
        maxEnt = max(obj.dofMap{id});
        obj.dofMap{id}(entList(isInactive)) =  nComp*(maxEnt:maxEnts+nNewEnts)+1;

        % update the ranges with new number of entities
        obj.fields(id).range(2) = obj.fields(id).range(2) + nComp*nNewEnts;
        for k = id+1:numel(obj.fields)
          obj.fields(k).range = obj.fields(k).range + nComp*nNewEnts;
        end

      else
        % new variable field 

        obj.fields(id).variableName = varName;
        obj.fields(id).tags = tags;
        obj.numbComponents(end+1) = nComp;

        % return the entity of type fieldLocation for the given mesh tags
        entList = getEntities(fieldLocation,obj.grid,tags);
        totEnts = getNumberOfEntities(fieldLocation,obj.grid);
        totActiveEnts = length(entList);

        % populate the dof map
        obj.dofMap{id} = zeros(totEnts,1);
        obj.dofMap{id}(entList) = nComp*(0:totActiveEnts-1)'+1;

        % update number of variables and dof counter
        obj.fields.range = [obj.totDofs+1,obj.totDofs+nComp*totActiveEnts];
        obj.totDofs = obj.totDofs + nComp*totActiveEnts;
        obj.nVars = obj.nVars+1;
      end

    end


    function dofs = getLocalDoF(obj,id,ents)
      % return the DoF numbering in global indexing for input entities
      % in isempty(varargin) all dofs are returned in correct order
      % the id is the result of a call to getVariableId(varName)

      if isscalar(id)
        dofs = repelem(obj.dofMap{id}(ents),ncomp) + ...
          repmat((0:ncomp(id)-1)',numel(ents),1);
      else
        assert(isempty(ents),"entity input in getDoF is admitted " + ...
          "only for single variable query")

        for i = 1:numel(id)
          dofs = [dofs; getLocalDoF(id(i))];
        end
      end
    end

    function dofs = getDoF(obj,id,ents)
      % return the DoF numbering in global indexing for input entities
      % in isempty(varargin) all dofs are returned in correct order
      % the id is the result of a call to getVariableId(varName)

      if isscalar(id)
        dofs = getLocalDoF(obj,id) + obj.fields(id).range(1);
      else
        assert(isempty(ents),"entity input in getDoF is admitted " + ...
          "only for single variable query")

        for i = 1:numel(id)
          dofs = [dofs; getDoF(id(i))];
        end
      end

    end


    function activeEnts = getActiveEntities(obj,varName,flagExpand)
      assert(isscalar(string(varName)),"The variable name must be a " + ...
        "string scalar or a character vector")
      id = obj.getVariableId(varName);
      activeEnts = find(obj.dofMap{id});

      % expand to account for component
      if flagExpand
        activeEnts = obj.dofExpand(activeEnts,obj.numbComponents(id));
      end
    end

    function tags = getTargetRegions(obj,varName)
      id = obj.getVariableId(varName);
      tags = obj.fields(id).tags;
    end

    function id = getVariableId(obj,varName)
      % return the id of the requested input variable
      assert(isVariable(obj,varName),"Requested variable is not available" + ...
        "in the DoFManager")
      id = find(any(strcmp([obj.fields.variableName],varName)));
    end

    function fl = isVariable(obj,varName)
      varId = strcmp([obj.fields.variableName],varName);
      fl = any(varId);
    end

    function numDof = getNumbDof(obj,varName)
      if isempty(varname)
        numDof = obj.totDofs;
      else
        id = getVariableId(obj,varName);
        r = obj.fields(id).range;
        numDof = r(2)-r(1)+1;
      end
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