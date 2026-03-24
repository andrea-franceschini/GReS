function addBC(obj,varargin)
% ADDBC  Add a boundary condition to the Boundaries object.
%
%   ADDBC(obj, Name, Value, ...) registers a boundary condition in the
%   simulation object OBJ using name-value pair arguments.
%
% -------------------------------------------------------------------------
% INPUT ARGUMENTS (Name-Value Pairs)
% -------------------------------------------------------------------------
%
%   'name'           - (string) Unique identifier for the boundary
%   condition.
%                      Must not already exist in the boundary condition
%                      database. An error is raised if the name is already
%                      in use. Example: 'inlet_velocity'
%
%   'variable'       - (string) Name of the physical variable to which the
%                      boundary condition applies (e.g., 'u', 'pressure').
%                      The value is case-insensitive. Example: 'velocity'
%
%   'field'          - (string) Source entity field that defines the
%   spatial
%                      support of the boundary condition. Used internally
%                      to resolve the entity type (e.g., 'node', 'face').
%                      The value is case-insensitive. Example: 'face'
%
%   'entityList'     - (double array) List of entity indices (e.g., node
%   IDs
%                      or face IDs) to which the boundary condition is
%                      applied. Example: [1, 2, 5, 10]
%
%   'entityListType' - (string) Specifies how 'entityList' is interpreted.
%                      Determines the type of entities being listed.
%                      Example: 'faceTag', 'nodeID'
%
%   'type'           - (string) Type of the boundary condition. Standard
%   types
%                      include:
%                        'dirichlet' - Essential BC (automatically sets
%                                      essential = true).
%                        'neumann'   - Natural BC (essential = false).
%                        <custom>    - Any user-defined type recognized by
%                                      BCtype.isCustomBC(). See note on
%                                      'essential' below.
%                      The value is case-insensitive. NOTE: The correctness
%                      of the BC type name is NOT validated here. It is the
%                      user's responsibility to provide a valid type.
%                      Validation occurs at the physics-solver level.
%
%   'components'     - (optional) Component indices or identifiers
%   specifying
%                      which components of a vector/tensor variable the BC
%                      applies to. If omitted or missing, the BC applies to
%                      all components. Example: [1, 2]  (apply to x- and
%                      y-components only)
%
%   'essential'      - (logical, optional) Relevant only for CUSTOM BC
%   types.
%                      Explicitly sets whether the BC is essential (true)
%                      or natural (false).
%                        - If not provided for a custom BC, the BC is
%                        treated
%                          as NON-essential and a warning is issued.
%                        - Ignored for standard types ('dirichlet' is
%                        always
%                          essential; 'neumann' is always non-essential).
%                      Example: true
%
%   'BCevent'        - (optional) One or more BC event descriptors to be
%                      registered alongside the boundary condition. Each
%                      event is passed to ADDBCEVENT individually. Example:
%                      myEventStruct
%
% -------------------------------------------------------------------------
% USAGE EXAMPLES
% -------------------------------------------------------------------------
%
%   % Dirichlet BC on a set of faces obj.addBC('name',
%   'wall_noslip',    ...
%             'variable',       'velocity',       ... 'field',
%             'face',           ... 'entityList',     [10, 11, 12],     ...
%             'entityListType', 'faceTag',        ... 'type',
%             'dirichlet');
%
%   % Custom BC with explicit essential flag obj.addBC('name',
%   'inlet_custom',   ...
%             'variable',       'pressure',       ... 'field',
%             'face',           ... 'entityList',     [5],              ...
%             'entityListType', 'faceTag',        ... 'type',
%             'myCustomType',   ... 'essential',      false);
%
%   % BC applied to specific components only obj.addBC('name',
%   'symmetry_x',     ...
%             'variable',       'velocity',       ... 'field',
%             'face',           ... 'entityList',     [3, 4],           ...
%             'entityListType', 'faceTag',        ... 'type',
%             'dirichlet',      ... 'components',     1);
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%
%   - BC names must be unique within OBJ. Reusing a name raises an error. -
%   The 'type' field is case-insensitive but must be a valid name
%     recognized at the physics-solver level. No validation is performed
%     here.
%   - For custom BC types, always specify 'essential' explicitly to avoid
%     unintended behavior and suppress the automatic warning.
%
% -------------------------------------------------------------------------
% SEE ALSO
% -------------------------------------------------------------------------
%   addBCEvent, BoundaryEntities, BCtype




default = struct('name',string.empty,...
  'variable',string.empty,...
  'field',string.empty,...
  'entityList',double.empty,...
  'entityListType',string.empty,...
  'type',string.empty,...
  'components',missing,...
  'essential',missing);

params = readInput(default,varargin{:});

name = params.name;
if obj.db.isKey(name)
  error("'%s' boundary condition name already defined", name);
end

srcEnt = entityField(lower(params.field));
variable = lower(params.variable);
type = lower(params.type);

% BC type
% THE USER IS RESPONSIBLE FOR INPUTTING THE CORRECT BCTYPE NAME.
% THIS IS ONLY CHECKED AT THE PHYSICS_SOLVER LEVEL

essentialFlag = false;

if BCtype.isCustomBC(type)
  if ~ismissing(params.essential)
    essentialFlag = logical(params.essential);
  else
    gresLog().warning(1, sprintf( ...
      "Custom boundary condition '%s' has been automatically considered to be NOT essential!\n" + ...
      "To override this behavior, specify the logical field 'essential' in the Boundary condition input.", ...
      bctype));
  end
elseif BCtype(type) == BCtype.dirichlet
  essentialFlag = true;
end

bc = struct('data', [], 'essential',essentialFlag,...
  'sourceField',srcEnt,'type', type, 'variable', variable);

if ismissing(params.components)
  comp = [];
else
  comp = params.components;
end

bcEnt = BoundaryEntities(name,srcEnt);
bcEnt.setEntities(params.entityListType,...
  params.entityList,...
  comp,obj.grid.topology);

bc.data = bcEnt;

% add BC to the database
obj.db(name) = bc;

% add the bc event
if isfield(params,"BCevent")
  for i = 1:numel(params.BCevent)
    addBCEvent(obj,name,params.BCevent(i))
  end
end

% finalize the boundary condition
% obj.computeBoundaryProperties(name);

end

