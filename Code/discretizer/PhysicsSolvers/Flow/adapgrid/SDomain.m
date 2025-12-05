classdef SDomain < handle
   % STRUCTUREDGRID Class to control all information related to the
   % structured grid

   properties (Access = public)
      grid         % ijk grid representing each cell
      boundaries   % the boundaries for this sub-domain.
   end

   properties (Access = private)
      numCells    % Number of exist cells
      numNodes    % Number of exist nodes
      numFaces    % Number of exist faces
   end

   methods(Access = public)
      function obj = SDomain(data)
         %ADAPGRID Construct an instance of this class
         obj.constructIJK(data.faces,1,[0 0 0]);
         obj.numCells = data.mesh.nCells;
         obj.numNodes = data.mesh.nNodes;
         obj.numFaces = data.faces.nFaces;

         % Saving boundary condition information.
         bkeys = data.bcs.db.keys;
         obj.boundaries = containers.Map('KeyType','char','ValueType','any');
         for i=1:length(bkeys)
           it = data.bcs.db(bkeys{i});
           physics = data.model.findPhysicsFromID(data.model.findIDPhysics(it.physics));
           if strcmp(physics,'Flow')
             obj.boundaries(bkeys{i})=GrowingBoundary(it.data,'Constant');
           end
         end
      end

      function [n_cells, n_ijk] = findNeighborhod(obj,cell_ijk)
         %FINDNEIGHBORHOD Find the neighboring information of a reference cell
         n_ijk = [-1  0  0;
                   1  0  0;
                   0 -1  0;
                   0  1  0;
                   0  0 -1;
                   0  0  1];
         n_ijk=repelem(cell_ijk,6,1)+n_ijk;
         n_cells = zeros(6,1);
         for i=1:6
            pos = find(ismember(obj.grid, n_ijk(i,:), 'rows'));
            if ~isempty(pos)
               n_cells(i) = pos;
            end
         end
      end

      function [n_cells, n_ijk] = findNeighborhodOfCell(obj,cellID)
         %FINDNEIGHBORHOD Find the neighboring information of a reference cell
         cell_ijk = obj.grid(cellID,:);
         [n_cells, n_ijk] = findNeighborhod(obj,cell_ijk);
      end

      function updateBorder(obj,data,cell,ijk_ref,ijk_new)
        keys = obj.boundaries.keys;
        dataRef = obj.findFacesCells(data.faces,ijk_ref);
        dataNew = obj.findFacesCells(data.faces,ijk_new);
        for i=1:length(keys)
          ref = obj.boundaries(keys{i});
          ref.update(dataNew,dataRef,cell.faceGrow);
        end
      end

      function out = findFacesCells(obj,data,ijk)
        [cells, ~] = obj.findNeighborhod(ijk);
        out.cell=cells';
        cellRef = find(ismember(obj.grid,ijk,'rows'));
        map = repelem(1:length(data.mapF2E)-1,diff(data.mapF2E))';
        map = (map==cellRef);
        out.faces = data.faces2Elements(map,1)';
        out.faces(data.faces2Elements(map,2)) = out.faces;
      end

      % Functions to update the structure.
      function flag = updateIJK(obj,ijk,ncells,nnodes,nfaces)
         % Update this class, increasing one cell.
         obj.grid(end+1,:)=ijk;
         obj.numCells=ncells;
         obj.numNodes=nnodes;
         obj.numFaces=nfaces;
         flag = true;
      end


   end

   methods (Access = private)
      % Function related to construct this class.
      function constructIJK(obj,face,cell,number)
         %CONSTRUCTIJK Construct the ijk instance for each cell in the grid
         ncells = length(face.mapF2E)-1;
         facesBelement = diff(face.mapF2E);
         face2Element = repelem(1:ncells,facesBelement)';

          % Test to check if all elements is hexahedra
         assert(isempty(find(facesBelement~=6, 1))==1);

         obj.grid = zeros(ncells,3);
         checkList = cell;         % FIFO list with cell to check  
         checkedList = [];         % FIFO list with the cell checked
         obj.grid(cell,:) = number; % for a generic numbering.
         count = 1;
         while (~isempty(checkList)) && (count<ncells)
            ref_cell = checkList(1);
            faces = face.faces2Elements(find(face2Element==ref_cell),:);
            ref_neigh = sum(face.faceNeighbors(faces(:,1),:),2)-ref_cell;
            for i=1:6
               if (~ismember(ref_neigh(i), checkedList)) && (~ismember(ref_neigh(i), checkList)) && (ref_neigh(i) ~= 0)
                  checkList = [checkList, ref_neigh(i)];
                  obj.grid(ref_neigh(i),:) = obj.grid(ref_cell,:);
                  switch faces(i,2)
                     case 3
                        obj.grid(ref_neigh(i),2) = obj.grid(ref_neigh(i),2) - 1;
                     case 4
                        obj.grid(ref_neigh(i),2) = obj.grid(ref_neigh(i),2) + 1;
                     case 1
                        obj.grid(ref_neigh(i),1) = obj.grid(ref_neigh(i),1) - 1;
                     case 2
                        obj.grid(ref_neigh(i),1) = obj.grid(ref_neigh(i),1) + 1;
                     case 5
                        obj.grid(ref_neigh(i),3) = obj.grid(ref_neigh(i),3) - 1;
                     case 6
                        obj.grid(ref_neigh(i),3) = obj.grid(ref_neigh(i),3) + 1;
                  end
               end
            end

            % Update the 2 list.
            count=count+1;
            checkedList = [checkedList, ref_cell];
            if length(checkList)>1
               checkList = checkList(2:end);
            else
               checkList =[];
            end
         end
      end
      
   end

   methods (Static)

   end

end