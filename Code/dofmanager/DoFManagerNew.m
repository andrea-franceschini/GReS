classdef DoFManagerNew < handle
  % Degree of freedom manager
  % Map each entity to a corresponding degree of freedom in the solution
  % system
  % 2 ordering can be specified in input:
  % Field-based ordering (default)
  % Domain-based ordering
  properties (Access = private)
    mesh
    dofMap % cell array with dof map for each variable field
    numbComponents
    fields = struct("variableName",[],...
      "range",[],...
      "tags",[]);
    nVars = 0;
    totDofs
  end


  methods (Access = public)
    function obj = DoFManagerNew(mesh)
      obj.mesh = mesh;
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

        % This is a Phylosophical choice: we can register a variable field
        % only one time for each domain. While we can still define more
        % then one physicsSolver in each domain, those are required to act
        % on disjoint sets of variable fields. 
        % An obvious consequence is that we cannot define two istances of
        % the same physicsSolver in a single domain.

        error("Variable %s has already been registered. GReS variables can" + ...
          "only be registered once in each domain ", varName)

        % the following  code make sense if we allow registering the same
        % variable more than once
          
        % id = getVariableId(obj,varName);
        % assert(isempty(intersect(tags,obj.fields(id).tags)), ...
        %   "Variable field %i has been already defined on specified tags", ...
        %   tags);

        % % update the tags
        % obj.fields(id).tags = sort([obj.fields(id).tags tags]);
        % 
        % % update the entities with only new entities
        % entList = getEntities(fieldLocation,obj.mesh,tags);
        % isInactive = obj.dofMap{id}(entList) == 0;
        % nNewEnts = sum(isInactive);
        % maxEnt = max(obj.dofMap{id});
        % obj.dofMap{id}(entList(isInactive)) =  nComp*(maxEnt:maxEnts+nNewEnts)+1;
        % 
        % % update the ranges with new number of entities
        % obj.fields(id).range(2) = obj.fields(id).range(2) + nComp*nNewEnts;
        % for k = id+1:numel(obj.fields)
        %   obj.fields(k).range = obj.fields(k).range + nComp*nNewEnts;
        % end

      else
        % new variable field 

        obj.fields(id).variableName = varName;
        obj.fields(id).tags = tags;
        obj.numbComponents(end+1) = nComp;

        % return the entity of type fieldLocation for the given mesh tags
        entList = getEntities(fieldLocation,obj.mesh,tags);
        totEnts = getNumberOfEntities(fieldLocation,obj.mesh);
        totActiveEnts = length(entList);

        % populate the dof map
        obj.dofMap{id} = zeros(totEnts,1);
        obj.dofMap{id}(entList) = nComp*(0:totActiveEnts-1)'+1;

        % update number of variables and dof counter
        obj.fields(id).range = [obj.totDofs+1,obj.totDofs+nComp*totActiveEnts];
        obj.totDofs = obj.totDofs + nComp*totActiveEnts;
        obj.nVars = obj.nVars+1;
      end

    end


    function dofs = getLocalDoF(obj,id,ents)
      % return the DoF numbering in global indexing for input entities
      % in isempty(varargin) all dofs are returned in correct order
      % the id is the result of a call to getVariableId(varName)


      if isscalar(id)
        ncomp = obj.numbComponents(id);

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

    function ents = getLocalEnts(obj,id,ents)

      % get entities without dof expansion
      ncomp = obj.numbComponents(id);
      ents = round(((ncomp-1)+obj.dofMap{id}(ents))/ncomp);
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
      
      if ~isnumeric(varName)
        varName = string(varName);
      end

      assert(isscalar(varName),"Input variable must be a scalar string or a" + ...
        "character vector (or an integer)")
      id = obj.getVariableId(varName);
      activeEnts = find(obj.dofMap{id});

      % expand to account for component
      if nargin > 2
        if flagExpand
          activeEnts = obj.dofExpand(activeEnts,obj.numbComponents(id));
        end
      end
    end

    function tags = getTargetRegions(obj,varName)
      id = obj.getVariableId(varName);
      tags = [obj.fields(id).tags];
    end

    function cells = getFieldCells(obj,varName)
      tags = getTargetRegions(obj,varName);
      cells = getEntities(entityField.cell,obj.mesh,tags);
    end

    function id = getVariableId(obj,varName)
      % return the id of the requested input variable
      if isnumeric(varName)
        assert(all(varName > 0) && all(varName < numel(obj.fields)),"Input variable" + ...
          "ID is not included in the domain")

        id = varName;

      else
        assert(isVariable(obj,varName),"Requested variable is not available" + ...
          "in the DoFManager")
        id = find(strcmp([obj.fields.variableName],varName));
      end
    end

    function fl = isVariable(obj,varName)
      varId = strcmp([obj.fields.variableName],varName);
      fl = any(varId);
    end

    function numDof = getNumbDoF(obj,varName)
      if isempty(varName)
        numDof = obj.totDofs;
      else
        id = getVariableId(obj,varName);
        r = obj.fields(id).range;
        numDof = r(2)-r(1)+1;
      end
    end

    function nVars = getNumberOfVariables(obj)
      nVars = obj.nVars;
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