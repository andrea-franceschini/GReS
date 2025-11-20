classdef DoFManagerNew < handle
  % Degree of freedom manager
  % Map each entity to a corresponding degree of freedom in the solution
  % system
  % 2 ordering can be specified in input:
  % Field-based ordering (default)
  % Domain-based ordering
  properties (Access = private)
    mesh
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
    function obj = DoFManagerNew(mesh)
      obj.mesh = mesh;
      obj.totDofs = 0;
    end

    function registerVariable(obj,varName,fieldLocation,nComp,tags)
      % varName: the name of the variable field
      % fieldLocation: a enum of type entityField
      % tags: the cellTag (or surfaceTag for lower dimensional fields) where the variable is actually present
      % numbComponents: the number of dofs per entitiy

      % return an instance of the registeredVariable as an instance of
      % entityField()

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
        obj.fields(id).fieldLocation = fieldLocation;
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

      % create the entity field


    end


    function dofs = getLocalDoF(obj,id,ents)

      % return the DoF numbering in global indexing for input entities
      % in isempty(varargin) all dofs are returned in correct order
      % the input id is the result of a call to getVariableId(varName)


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

      % same as getLocalDoF but without component expansion
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
      id = obj.getVariableId(varId);
      tags = [obj.fields(id).tags];
    end

    function cells = getFieldCells(obj,varId)
      tags = getTargetRegions(obj,varId);
      cells = getEntities(entityField.cell,obj.mesh,tags);
    end

    function location = getFieldLocation(obj,varId)

      id = getVariableId(varId);
      location = obj.fields(id).location;


    end

    function id = getVariableId(obj,varId)
      % return the id of the requested input variable
      if isnumeric(varId)
        assert(all(varId > 0) && all(varId <= numel(obj.fields)),"Input variable" + ...
          "ID is not included in the domain")

        id = varId;

      else
        assert(isVariable(obj,varId),"Requested variable is not available" + ...
          "in the DoFManager")
        id = find(strcmp([obj.fields.variableName],varId));
      end
    end

    function fl = isVariable(obj,varId)
      varId = strcmp([obj.fields.variableName],varId);
      fl = any(varId);
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


  end

  methods (Static)

    function dofOut = dofExpand(dofIn,nComp)
      % make sure input is a column vector
      dofIn = reshape(dofIn,[],1);
      dofOut = nComp*repelem(dofIn,nComp,1)-repmat((nComp-1:-1:0)',size(dofIn,1),1);
    end
  end
end