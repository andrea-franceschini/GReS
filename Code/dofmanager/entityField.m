classdef entityField
  enumeration
    node
    face
    cell
    surface
  end


  methods
    function ents = getEntities(obj, grid, tags)
      % GETENTITIES returns the entities of this type from the mesh
      % obj      : the enum instance (entityField.node, etc.)
      % mesh     : mesh object containing cells, faces, etc.
      % tags     : optional tags to filter entities

      if nargin < 3
        tags = [];
      end

      switch obj
        case entityField.node
          if isempty(tags)
            ents = unique(grid.topology.cells(:));
          else
            cellsID = ismember(grid.topology.cellTag, tags);
            ents = unique(grid.cells(cellsID,:));
          end

        case entityField.face
          fN = grid.faces.faceNeighbors;
          if isempty(tags)
            ents = 1:size(fN,1);  % all faces
          else
            cellsID = ismember(grid.topology.cellTag, tags);
            ents = find(any(ismember(fN, find(cellsID)),2));
          end

        case entityField.cell
          if isempty(tags)
            ents = 1:size(grid.cells,1);
          else
            ents = find(ismember(grid.topology.cellTag, tags));
          end

        case entityField.surface
          % Replace with your interface extraction logic
          if isempty(tags)
            ents = 1:size(grid.surfaces,1);
          else
            ents = find(ismember(grid.topology.surfaceTag, tags));
          end
        otherwise
          error('Unknown entityField type.');
      end
    end

    function nEnts = getNumberOfEntities(obj,mesh,tags)
      ents = getEntities(obj,mesh,tags);
      nEnts = numel(ents);
    end


  end


end


