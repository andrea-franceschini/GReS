classdef adapgrid < handle
   % ADAPGRID Class to control the mesh grow

   properties (Access = public)
      ijk         % ijk information for each cell
      % fFacesCells % number of free faces for each cell
   end

   properties (Access = private)
      nodes       % local nodes for the new cell
      faces       % local faces for the new cell
      coordinates % local coordinates for the new cell
      numNewNodes % Number of nodes to be created
      numNewFaces % Number of faces to be created
      numCells    % Number of exist cells
      numNodes    % Number of exist nodes
      numFaces    % Number of exist faces
   end

   methods(Access = public)
      function obj = adapgrid(grid)
         %ADAPGRID Construct an instance of this class
         obj.constructIJK(grid.faces,1,[0 0 0]);
         obj.numCells = grid.topology.nCells;
         obj.numNodes = grid.topology.nNodes;
         obj.numFaces = grid.faces.nFaces;

         % obj.constructFreeFaces(grid.faces);
      end

      function addCell(obj,grid,cell,dir)
         %ADDCELL Add a cell in the grid and update the mesh
         
         % Find the position to where to grow the mesh.
         [n_cells,n_ijk] = obj.findNeighborhod(obj.ijk(cell,:));
         if obj.canGrow(dir,n_cells)
            % Find the information about the cell to be created.
            cell_ijk_grow=n_ijk(dir,:);
            obj.findExistData(grid,cell_ijk_grow);
            obj.fillExistData(dir,0.5);

            % Update the mesh
            obj.updateFaces(grid.faces,cell_ijk_grow);
            obj.updateElements(grid.cells,cell_ijk_grow);
         end         
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

         obj.ijk = zeros(ncells,3);
         checkList = cell;         % FIFO list with cell to check  
         checkedList = [];         % FIFO list with the cell checked
         obj.ijk(cell,:) = number; % for a generic numbering.
         count = 1;
         while (~isempty(checkList)) && (count<ncells)
            ref_cell = checkList(1);
            faces = face.faces2Elements(find(face2Element==ref_cell),:);
            ref_neigh = sum(face.faceNeighbors(faces(:,1),:),2)-ref_cell;
            for i=1:6
               if (~ismember(ref_neigh(i), checkedList)) && (~ismember(ref_neigh(i), checkList)) && (ref_neigh(i) ~= 0)
                  checkList = [checkList, ref_neigh(i)];
                  obj.ijk(ref_neigh(i),:) = obj.ijk(ref_cell,:);
                  switch faces(i,2)
                     case 1
                        obj.ijk(ref_neigh(i),1) = obj.ijk(ref_neigh(i),1) - 1;
                     case 2
                        obj.ijk(ref_neigh(i),1) = obj.ijk(ref_neigh(i),1) + 1;
                     case 3
                        obj.ijk(ref_neigh(i),2) = obj.ijk(ref_neigh(i),2) - 1;
                     case 4
                        obj.ijk(ref_neigh(i),2) = obj.ijk(ref_neigh(i),2) + 1;
                     case 5
                        obj.ijk(ref_neigh(i),3) = obj.ijk(ref_neigh(i),3) - 1;
                     case 6
                        obj.ijk(ref_neigh(i),3) = obj.ijk(ref_neigh(i),3) + 1;
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

      function constructFreeFaces(obj,face)
         %CONSTRUCTFREEFACES Construct the number of free faces by cell.
         ncells = length(face.mapF2E)-1;

         faceA = face.faceNeighbors(:,1)==0;
         faceB = face.faceNeighbors(:,2)==0;
         cells = [face.faceNeighbors(faceA,2);face.faceNeighbors(faceB,1)];
         obj.fFacesCells = accumarray(cells,ones(length(cells),1),[ncells 1]);
      end



      % Functions to identify information in the ijk cell.
      function [n_cells, n_ijk] = findNeighborhod(obj,cell_ijk)
         %FINDNEIGHBORHOD Find the neighboring information of a reference cell
         % cell_ijk = obj.ijk(cellID,:);
         n_ijk = [ -1  0  0;
                    1  0  0;
                    0 -1  0;
                    0  1  0;
                    0  0 -1;
                    0  0  1];
         n_ijk=repelem(cell_ijk,6,1)+n_ijk;
         n_cells = zeros(6,1);
         for i=1:6
            pos = find(ismember(obj.ijk, n_ijk(i,:), 'rows'));
            if ~isempty(pos)
               n_cells(i) = pos;
            end
         end
      end

      function findExistData(obj,grid,ijk_ref)
         %FINDEXISTDATA Find the existent data for the new cell.

         % Some important variables.
         local_nodes(1:8) = 0;
         mapEFaces = repelem(1:obj.numCells,diff(grid.faces.mapF2E))';
         mapNFaces = repelem(1:obj.numFaces,diff(grid.faces.mapN2F))';
         ref_pos = [6 5 1 2       % sth (1)
                    3 4 8 7       % nth (2)
                    5 8 4 1       % wst (3)
                    2 3 7 6       % est (4)
                    1 4 3 2       % bot (5)
                    6 7 8 5       % top (6)
                    ];

         % Find the cells around.
         [n_cells,~] = obj.findNeighborhod(ijk_ref);

         % Build the information already provided.
         local_faces(1:6)=0;
         local_coord(1:8,1:3)=0;
         for i=1:6
            if(n_cells(i)>0)
               % Finding the face.
               dir = i+2*mod(i,2)-1;
               faceNum = grid.faces.faces2Elements(mapEFaces==n_cells(i),:);
               faceNum = faceNum(faceNum(:,2)==dir,1);

               % Copy information.
               local_faces(i)=faceNum;
               local_nodes(ref_pos(i,:))=grid.faces.nodes2Faces(mapNFaces==faceNum);
               local_coord(ref_pos(i,:),:) = grid.topology.coordinates(local_nodes(local_nodes~=0),:);
            end
         end

         % store the information in this object.
         obj.faces = local_faces;
         obj.nodes = local_nodes;
         obj.coordinates = local_coord;
         obj.numNewNodes = sum(local_nodes==0);
         obj.numNewFaces = sum(local_faces==0);
      end

      function fillExistData(obj,dir,len)
         %FILLEXISTDATA Complement the information necessary to create the
         % new cell.

         % Mask to some Variables.
         nodes2Add = obj.numNewNodes;
         faces2Add = obj.numNewFaces;
         lastNode = obj.numNodes;
         lastFace = obj.numFaces;

         % Direction of the Cell Grow.         
         ijk_loc = mod(ceil(dir/2)-1, 3) + 1;       

         % Creating Nodes and Coordenates.
         if nodes2Add~=0
            % Nodes that where copy from the others in the direction.
            ref = [2 1 4 3 6 5 8 7;
                4 3 2 1 8 7 6 5;
                5 6 7 8 1 2 3 4];
            ref = ref(ijk_loc,:);

            % Indication of the new nodes
            ind = obj.nodes==0;

            % Create the new nodes numbers
            obj.nodes(ind)=lastNode+1:lastNode+nodes2Add;

            % Create the nodes positions
            obj.coordinates(ind,:) = obj.coordinates(ref(ind),:);
            obj.coordinates(ind,ijk_loc) = obj.coordinates(ind,ijk_loc) + (1-2*mod(dir,2))*len;
         end

         % Creating Faces.
         if faces2Add~=0
            % Indication of the new faces
            ind = obj.faces==0;

            % Create the new faces numbers
            obj.faces(ind)=lastFace+1:lastFace+faces2Add;
         end

         % 
         % 
         % % Find the cells around.
         % infFaces = obj.faces;
         % facesNOTLook = find(obj.faces~=0);
         % facesLook = facesNOTLook+2*mod(facesNOTLook,2)-1;
         % infFaces(facesLook)=1;
         % facesLook = [facesLook find(infFaces==0)];
         % 
         % % Build the information for the other cells
         % lastFace = grid.faces.nFaces;
         % lastCell = grid.topology.nCells;
         % lastNode = grid.topology.nNodes;
         % 
         % 
         % 
         % for i=1:length(facesLook)
         %    obj.nodes==0
         %    % Finding the face.
         %    % dir = i+2*mod(i,2)-1;
         %    faceNum = grid.faces.faces2Elements(mapEFaces==n_cells(i),:);
         %    faceNum = faceNum(faceNum(:,2)==dir,1);
         % 
         %    % Copy information.
         %    faces(i)=faceNum;
         %    nodes(ref_pos(i,:))=grid.faces.nodes2Faces(mapNFaces==faceNum);
         % end
      end



      % Functions to check the mesh
      function flag = canGrow(obj,dir,n_cells)
         %CANGROW Function to check if a cell can be add in the mesh, for
         %the direction specified.
         flag = n_cells(dir)==0;
      end



      % Functions to update the mesh
      function updateFaces(obj,faces,refCell)
         %UPDATEFACES update the face information to add the new cell

         % Mask to some Variables.
         faces2Add = obj.numNewFaces;
         lastFace = obj.numFaces;
         % nodes2Add = obj.numNewNodes;
         % lastNode = obj.numNodes;

         % Find the cells around.
         [neighcells,~] = obj.findNeighborhod(refCell);

         % Some Important Variables.
         % mapEFaces = repelem(1:obj.numCells,diff(faces.mapF2E))';
         % mapNFaces = repelem(1:obj.numFaces,diff(faces.mapN2F))';
         ref_pos = [6 5 1 2       % sth (1)
                    3 4 8 7       % nth (2)
                    5 8 4 1       % wst (3)
                    2 3 7 6       % est (4)
                    1 4 3 2       % bot (5)
                    6 7 8 5       % top (6)
                    ];

         faceCentroid(1:faces2Add,1:3)=0;
         faceNormal(1:faces2Add,1:3)=0;

         % Variables to add with simple construction.
         mapN2F(1:faces2Add)=max(faces.mapN2F)+4*(1:faces2Add);
         faces2Elements=[obj.faces' (1:6)'];
         mapF2E=max(faces.mapF2E)+6;

         % creating information for each face
         nodes2Faces(1:4*faces2Add)=0;
         faceNeighbors(1:faces2Add,1:2)=obj.numCells+1;
         for i=1:faces2Add
            faceNum = obj.numFaces+i;
            faceNodes = ref_pos(obj.faces == faceNum,:);
            faceNodeGL = obj.nodes(faceNodes);
            faceCentroid(i,:)=sum(obj.coordinates(faceNodes,:))/4;
            direction = faces2Elements(faces2Elements(:,1)==faceNum,2);

            [~,nodeC]=min(faceNodeGL);
            nodeL = mod(nodeC - 2, 4) + 1;
            nodeR = mod(nodeC, 4) + 1;
            nodeA = mod(nodeC + 1, 4) + 1;
            if (faceNodeGL(nodeL)<faceNodeGL(nodeR))
               nodes2Faces(4*(i-1)+1:4*i)=[
                  faceNodeGL(nodeC) faceNodeGL(nodeL) ...
                  faceNodeGL(nodeA) faceNodeGL(nodeR)];
               faceNeighbors(i,1)=neighcells(direction);   % <---(i,1) or (i,2) to check
            else
               nodes2Faces(4*(i-1)+1:4*i)=[
                  faceNodeGL(nodeC) faceNodeGL(nodeR) ...
                  faceNodeGL(nodeA) faceNodeGL(nodeL)];
               faceNeighbors(i,2)=neighcells(direction);   % <---(i,1) or (i,2) to check
            end
            v1 = obj.coordinates(faceNodes(nodeL),:)-obj.coordinates(faceNodes(nodeC),:);
            v2 = obj.coordinates(faceNodes(nodeA),:)-obj.coordinates(faceNodes(nodeC),:);
            % v3 = obj.coordinates(faceNodes(nodeR),:)-obj.coordinates(faceNodes(nodeC),:);
            % area = 0.5*norm(cross(v1,v2)) + 0.5*norm(cross(v2,v3));
            faceNormal(i,:)=cross(v1,v2);
         end

         % UPDATE THE FACE STRUCT.
         faces.nodes2Faces = [faces.nodes2Faces; nodes2Faces'];
         faces.mapN2F = [faces.mapN2F; mapN2F'];
         faces.faces2Elements = [faces.faces2Elements; faces2Elements];
         faces.mapF2E = [faces.mapF2E; mapF2E];
         faces.faceNeighbors = [faces.faceNeighbors; faceNeighbors];
         faces.faceCentroid = [faces.faceCentroid; faceCentroid];
         faces.faceNormal = [faces.faceNormal; faceNormal];
         faces.nFaces = faces.nFaces+faces2Add;

         facesNot2Add = 6 - faces2Add;
         existFaces = faces2Elements(faces2Elements(:,1)<=lastFace,1);
         for i=1:facesNot2Add
            loc=faces.faceNeighbors(existFaces(i),:)==0;
            faces.faceNeighbors(existFaces(i),loc)=obj.numCells+1;
         end
      end

      function updateElements(obj,grid,refCell)
         %UPDATEFACES update the face information to add the new cell
         center = mean(obj.coordinates);
         V = obj.coordinates;

         tets = [ 1 2 4 5;
            2 3 4 7;
            2 5 6 7;
            4 5 7 8;
            2 4 5 7 ];

         vol = 0;
         for i = 1:size(tets,1)
            A = V(tets(i,1),:);
            B = V(tets(i,2),:);
            C = V(tets(i,3),:);
            D = V(tets(i,4),:);

            % Compute volume of tetrahedron
            tet_vol = abs(dot(B - A, cross(C - A, D - A))) / 6;
            vol = vol + tet_vol;
         end



         g=1;
      end

      function updateMesh(obj)
      end


   end

   methods (Static)
      % Some moviment in the grid
      function ijk_loc = cell_IJK_from_Ref(ijk_neig,dir)
         %IJK_LOC Function to find cell ijk to where to grow
         switch dir
            case 1
               ijk_loc = ijk_neig(dir,:)+[1 0 0];
            case 2
               ijk_loc = ijk_neig(dir,:)-[1 0 0];
            case 3
               ijk_loc = ijk_neig(dir,:)+[0 1 0];
            case 4
               ijk_loc = ijk_neig(dir,:)-[0 1 0];
            case 5
               ijk_loc = ijk_neig(dir,:)+[0 0 1];
            case 6
               ijk_loc = ijk_neig(dir,:)-[0 0 1];
         end
      end

      function coord = growDirection(pos,dir,len)
         coord = pos;
         switch dir
            case 1
               coord(1)=coord(1)-len;
            case 2
               coord(1)=coord(1)+len;
            case 3
               coord(2)=coord(2)-len;
            case 4
               coord(2)=coord(2)+len;
            case 5
               coord(3)=coord(3)-len;
            case 6
               coord(3)=coord(3)+len;
         end
      end



   end



end