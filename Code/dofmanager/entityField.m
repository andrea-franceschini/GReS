classdef entityField

  % general class for handling a generic field in the problem

  % the class is useful for both purely geometric information and also
  % physical information for variable fields

  % the location of the field
  enumeration
    node
    face
    cell
    surface
  end

  % properties needed only for variable field
  % properties
  %   mesh
  %   variableName
  %   fieldId                % the position of the field in the system matrix of the domain
  %   numbComponents
  % end

  methods

    function varargout = getIncidenceID(target,mesh,source,sourceList)
      % given entities of type 'source', return the connected entities of
      % type 'target'.
      % sourceList - optional input with id of source entities

      % output: list: list of indices connected, ptr: pointer mapping each
      % source id with the first corresponding target entity in the list

      if nargin < 4
        sourceList = [];
      end

      switch source
        case entityField.node
          [list,ptr] = getIncidenceFromNode(target,mesh,sourceList);
        case entitiyField.surface
          [list,ptr] = getSurfaceIndicence(target,mesh,sourceList);
        case entityField.cell
          [list,ptr] = getCellIncidence(target,mesh,sourceList);
      end

      varargout{1} = list;
      
      if nargout > 1
        varargout{2} = ptr;
      end
    end


    function [list,ptr] = getIncidenceFromNode(target,mesh,sourceList)

      if isempty(sourceList)
        sourceList = 1:mesh.nNodes;
      end

      switch target
        case entityField.node
          list = sourceList;
          ptr = reshape(1:numel(list),[],1);
        otherwise
          error("Incidence from node to %s is not yet available")
      end

    end

    function [list,ptr] = getIndicenceFromSurface(target,mesh,sourceList)

      if isempty(sourceList)
        sourceList = 1:mesh.nSurfaces;
      end

      switch target
        case entityField.node
          ptr = mesh.surfaceNumVerts(sourceList);
          tmp = mesh.surfaces(sourceList,:);
          list = reshape(tmp',[],1);
        case entityField.surface
          list = sourceList;
          ptr = reshape(1:numel(list),[],1);
        otherwise
          error("Incidence from surface to %s is not yet available")
      end

    end

    function [list,ptr] = getIndicenceFromCell(target,mesh,sourceList)

      if isempty(sourceList)
        sourceList = 1:mesh.nCells;
      end

      switch target
        case entityField.node
          ptr = mesh.cellNumVerts(sourceList);
          tmp = mesh.cells(sourceList,:);
          list = reshape(tmp',[],1);
        case entityField.cell
          list = sourceList;
          ptr = reshape(1:numel(list),[],1);
        otherwise
          error("Incidence from cell to %s is not yet available")
      end

    end


    function getIncidenceMap(target,grid,source,srcList)
      % given entities of type 'source', return the connected entities of
      % type 'target'.
      % sourceList - optional input with id of source entities

      % output: list: list of indices connected, ptr: pointer mapping each
      % source id with the first corresponding target entity in the list

      if nargin < 4
        sourceList = [];
      end

      switch source
        case entityField.node
          [list,ptr] = getIncidenceFromNode(target,mesh,sourceList);
        case entitiyField.surface
          [list,ptr] = getSurfaceIndicence(target,mesh,sourceList);
        case entityField.cell
          [list,ptr] = getCellIncidence(target,mesh,sourceList);
      end

      varargout{1} = list;

      if nargout > 1
        varargout{2} = ptr;
      end

      mesh = grid.topology;
      elem = grid.cell;


    end

  end

    % function getEntitiesInterpolation(target,grid,source,sourceList)
    % end

  end


  % methods
  % 
  % 
  %   function ents = getEntities(obj, mesh, tags)
  %     % GETENTITIES returns the entities of this type from the mesh
  %     % obj      : the enum instance (entityField.node, etc.)
  %     % mesh     : mesh object containing cells, faces, etc.
  %     % tags     : optional tags to filter entities
  % 
  %     if nargin < 3
  %       tags = [];
  %     end
  % 
  %     if mesh.nCells == 0
  %       assert(obj == entityField.node || obj == entityField.surface,...
  %         "Cannot retrieve entity %s from a 2D mesh object",obj);
  %       % 2D mesh as input
  %       source = mesh.surfaces;
  %       sourceTag = mesh.surfaceTag;
  %     else
  %       source = mesh.cells;
  %       sourceTag = mesh.cellTag;
  %     end
  % 
  %     switch obj
  %       case entityField.node
  %         if isempty(tags)
  %           ents = unique(source(:));
  %         else
  %           cellsID = ismember(sourceTag, tags);
  %           ents = unique(source(cellsID,:));
  %         end
  % 
  %       case entityField.face
  %         mesh = mesh.topology;
  %         faces = mesh.faces;
  %         fN = faces.faceNeighbors;
  %         if isempty(tags)
  %           ents = 1:size(fN,1);  % all faces
  %         else
  %           cellsID = ismember(sourceTag, tags);
  %           ents = find(any(ismember(fN, find(cellsID)),2));
  %         end
  % 
  %       case entityField.cell
  %         if isempty(tags)
  %           ents = 1:size(source,1);
  %         else
  %           ents = find(ismember(sourceTag, tags));
  %         end
  % 
  %       case entityField.surface
  % 
  %         if isempty(tags)
  %           ents = 1:size(source,1);
  %         else
  %           ents = find(ismember(sourceTag, tags));
  %         end
  %       otherwise
  %         error('Unknown entityField type.');
  %     end
  %   end

    % function nEnts = getNumberOfEntities(target,mesh,tags,source)
    %   % return the number of entities
    % 
    %   if nargin < 3
    %     target == source;
    %   end
    % 
    % 
    %     ents = getEntities(obj,mesh);
    %   else
    %     ents = getEntities(obj,mesh,tags);
    %   end
    %   nEnts = numel(ents);
    % end

    % function nEnts = getNumberOfEntities(target,mesh,tags,source)
    %   if nargin < 3
    %     ents = getEntities(obj,mesh);
    %   else
    %     ents = getEntities(obj,mesh,tags);
    %   end
    %   nEnts = numel(ents);
    % end

    % function ents = getEntityFromElement(outEntity,sourceEntity,mesh,el,nc)
    % 
    %   if sourceEntity == entityField.cell
    % 
    %     switch outEntity
    %       case entityField.node
    %         ents = mesh.cells(el,1:mesh.cellNumVerts(el));
    %       case entityField.cell
    %         ents = el;
    %       otherwise
    %         error("Method not yet supported for entity of type %s");
    %     end
    % 
    %   elseif sourceEntity == entityField.surface
    % 
    %     switch outEntity
    %       case entityField.node
    %         ents = mesh.surfaces(el,1:mesh.surfaceNumVerts(el));
    %       case entityField.surface
    %         ents = el;
    %       otherwise
    %         error("Method not yet supported for entity of type %s",outEntity);
    %     end
    % 
    %   else
    % 
    %     error("Method not yet supported for source entity %s",sourceEntity)
    % 
    %   end
    % 
    %   % if nargin > 3
    %   %   ents = DoFManager.dofExpand(ents,nc);
    %   % end
    % 
    % end


    % function map = getEntityMapping(targetEntity,sourceEntity,grid,varargin)
    %   % return a sparse matrix that perform geometric interpolation from a
    %   % sourceEntity to a targetEntitiy.
    %   % mesh (grid object handling all geometric information on the
    %   % entities)
    %   % since this is is an interpolation operator, rows sum to one
    %   % optional input: a list for a subset of source entities
    % 
    %   nS = getNumberOfEntities(sourceEntity);
    % 
    %   if targetEntity == sourceEntity
    %     map = speye(nS);
    %     return
    %   end
    % 
    % 
    %   if nargin > 3
    %     list = varargin{1};
    %   else
    %     list = reshape(1:nS,[],1);
    %   end
    % 
    %   if targetEntity ~= entityField.node
    %     error(['Mapping between different entities is currently ' ...
    %       'available only for nodal target'])
    %   end
    % 
    %   switch sourceEntity
    %     case entityField.surface
    %       getInf = @(el) 
    %     case entityField.cell
    %   end
    %   end






      % use getEntityFromElement to map the entities and preallocate the
      % matrix

      % then use findNodeVolume and findNodeAea


    %end


  %end

end


