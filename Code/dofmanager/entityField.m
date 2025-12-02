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

    
    function ents = getEntities(obj, mesh, tags)
      % GETENTITIES returns the entities of this type from the mesh
      % obj      : the enum instance (entityField.node, etc.)
      % mesh     : mesh object containing cells, faces, etc.
      % tags     : optional tags to filter entities

      if nargin < 3
        tags = [];
      end

      if mesh.nCells == 0
        assert(obj == entityField.node || obj == entityField.surface,...
          "Cannot retrieve entity %s from a 2D mesh object",obj);
        % 2D mesh as input
        source = mesh.surfaces;
        sourceTag = mesh.surfaceTag;
      else
        source = mesh.cells;
        sourceTag = mesh.cellTag;
      end

      switch obj
        case entityField.node
          if isempty(tags)
            ents = unique(source(:));
          else
            cellsID = ismember(sourceTag, tags);
            ents = unique(source(cellsID,:));
          end

        case entityField.face
          mesh = mesh.topology;
          faces = mesh.faces;
          fN = faces.faceNeighbors;
          if isempty(tags)
            ents = 1:size(fN,1);  % all faces
          else
            cellsID = ismember(sourceTag, tags);
            ents = find(any(ismember(fN, find(cellsID)),2));
          end

        case entityField.cell
          if isempty(tags)
            ents = 1:size(source,1);
          else
            ents = find(ismember(sourceTag, tags));
          end

        case entityField.surface

          if isempty(tags)
            ents = 1:size(source,1);
          else
            ents = find(ismember(sourceTag, tags));
          end
        otherwise
          error('Unknown entityField type.');
      end
    end

    function nEnts = getNumberOfEntities(obj,mesh,tags)
      if nargin < 3
        ents = getEntities(obj,mesh);
      else
        ents = getEntities(obj,mesh,tags);
      end
      nEnts = numel(ents);
    end

    function ents = getEntityFromElement(outEntity,sourceEntity,mesh,el,nc)

      if sourceEntity == entityField.cell

        switch outEntity
          case entityField.node
            ents = mesh.cells(el,1:mesh.cellNumVerts(el));
          case entityField.cell
            ents = el;
          otherwise
            error("Method not yet supported for entity of type %s");
        end

      elseif sourceEntity == entityField.surface

        switch outEntity
          case entityField.node
            ents = mesh.surfaces(el,1:mesh.surfaceNumVerts(el));
          case entityField.surface
            ents = el;
          otherwise
            error("Method not yet supported for entity of type %s",outEntity);
        end

      else 
        
        error("Method not yet supported for source entity %s",sourceEntity)

      end

      if nargin > 3
        ents = DoFManager.dofExpand(ents,nc);
      end

    end

  end


end


