function registerVariable(obj,varName,fldLoc,nComp,varargin)
% REGISTERVARIABLE  Register a new variable field in the simulation domain.
%
%   REGISTERVARIABLE(obj, varName, fldLoc, nComp, tags)
%   REGISTERVARIABLE(obj, varName, fldLoc, nComp, 'nEntities', N)
%
%   Registers a variable field on a set of mesh entities and builds the
%   corresponding DOF map. Each variable can only be registered once per
%   domain. Multiple physics solvers sharing the same domain must act on
%   disjoint sets of variable fields.
%
% -------------------------------------------------------------------------
% INPUT ARGUMENTS
% -------------------------------------------------------------------------
%
%   varName  - (string) Unique name of the variable field to register
%              (e.g., 'displacements', 'pressure'). Raises an error if the
%              name has already been registered in this domain.
%
%   fldLoc   - (entityField) Enum specifying the mesh entity type on which
%              the variable is defined (e.g., entityField.node,
%              entityField.surface,
%              entityField.cell).
%
%   nComp    - (integer) Number of degrees of freedom per entity
%              (e.g., 1 for a scalar, 3 for a 3D vector field).
%
%   varargin - Two calling conventions are supported:
%
%              REGISTERVARIABLE(obj, varName, fldLoc, nComp, tags)
%                tags  — (numeric array) Cell or surface tags identifying
%                        the mesh regions where the variable is active.
%                        Entities are derived from these tags via the mesh
%                        topology.
%
%              REGISTERVARIABLE(obj, varName, fldLoc, nComp, 'nEntities', N)
%                'nEntities' — keyword flag selecting the direct mode.
%                N           — (integer) Total number of active entities.
%                              Entities are assumed to be 1:N with no tag
%                              filtering.
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%
%   - A variable can only be registered once per domain. Attempting to
%     register the same name twice raises an error. This enforces that
%     physics solvers operating on the same domain register disjoint
%     variables
%
% -------------------------------------------------------------------------
% SEE ALSO
% -------------------------------------------------------------------------
%   entityField, getEntitiesList, getFieldCells

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
  % entList = getEntities(fldLoc,obj.mesh,tags);
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
  obj.fields(id).fieldLocation = fldLoc;
  obj.numbComponents(end+1) = nComp;


  if nargin < 6
    % return the entity of type fldLoc for the given mesh tags
    tags = varargin{1};
    obj.fields(id).tags = tags;
    cells = obj.getFieldCells(id);
    entList = getEntitiesList(fldLoc,obj.mesh,entityField.cell,cells);
    totActiveEnts = length(entList);
  else
    assert(strcmp(varargin{1},"nEntities"))
    totActiveEnts = varargin{2};
    entList = reshape(1:totActiveEnts,[],1);
  end

  totEnts = numel(getEntitiesList(fldLoc,obj.mesh,fldLoc));

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
