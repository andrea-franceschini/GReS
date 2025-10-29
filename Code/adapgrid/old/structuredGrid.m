classdef structuredGrid < handle
   % STRUCTUREDGRID Class to control all information related to the
   % structured grid

   properties (Access = public)
      ijk         % ijk information for each cell
      % fFacesCells % number of free faces for each cell
   end

   properties (Access = private)
      numCells    % Number of exist cells
      numNodes    % Number of exist nodes
      numFaces    % Number of exist faces
   end

   methods(Access = public)
      function obj = structuredGrid(mesh,faces)
         %ADAPGRID Construct an instance of this class
         obj.constructIJK(faces,1,[0 0 0]);
         obj.numCells = mesh.nCells;
         obj.numNodes = mesh.nNodes;
         obj.numFaces = faces.nFaces;
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

      function [newCell] = newCellData(obj,mesh,cells,ijk_ref,dir,len,refCell)
         %NEWCELLDATA Create the faces information for the new cell.

         % Some important variables.
         mapEFaces = repelem(1:obj.numCells,diff(cells.mapF2E))';
         mapNFaces = repelem(1:obj.numFaces,diff(cells.mapN2F))';
         ref_nodes = [ 1 2 6 5       % sth (1)
                       4 3 7 8       % nth (2)
                       1 5 8 4       % wst (3)
                       2 6 7 3       % est (4)
                       1 2 3 4       % bot (5)
                       5 6 7 8       % top (6)
                    ];

         % Find the cells around.
         [n_cells,~] = obj.findNeighborhod(ijk_ref);

         % Build the information already provided.
         newCell.face(1:6) = 0;
         newCell.node(1:8) = 0;
         newCell.coord(1:8,1:3) = 0;
         for face=1:6
            if(n_cells(face)>0)
               % Finding the face.
               faceMirror = face+2*mod(face,2)-1;
               faceNode = mesh.cells(n_cells(face),ref_nodes(faceMirror,:));
               newCell.node(ref_nodes(face,:)) = faceNode;

               faceNum = cells.faces2Elements(mapEFaces==n_cells(face),:);
               faceNum = faceNum(faceNum(:,2)==faceMirror,1);
               newCell.face(face)=faceNum;
            end
         end
         pos = find((newCell.node~=0)==1);
         newCell.coord(pos,:)=mesh.coordinates(newCell.node(pos),:);
         newCell.nNodes = sum(newCell.node==0);
         newCell.nFaces = sum(newCell.face==0);

         % Direction of the Cell Grow.         
         ijk_loc = mod(ceil(dir/2)-1, 3) + 1;       

         % Creating Nodes and Coordenates.
         if newCell.nNodes~=0
            % Nodes that where copy from the others in the direction.
            % ref = [2 1 4 3 6 5 8 7;
            %     4 3 2 1 8 7 6 5;
            %     5 6 7 8 1 4 3 2];
            % ref = [ 1 2 5 6 4 8 7 3;
            %         2 6 7 3 1 5 8 4;
            %         5 6 7 8 1 4 3 2 ];

            ref = [ 4 3 7 8 1 2 6 5;
                    1 5 8 4 2 6 7 3;
                    5 6 7 8 1 2 3 4 ];

            ref = ref(ijk_loc,:);

            % Indication of the new nodes
            ind = newCell.node==0;

            % Create the new nodes numbers
            newCell.node(ind)=mesh.nNodes+1:mesh.nNodes+newCell.nNodes;

            % Create the nodes positions
            newCell.coord(ind,:) = newCell.coord(ref(ind),:);
            newCell.coord(ind,ijk_loc) = newCell.coord(ind,ijk_loc) + (1-2*mod(dir,2))*len;
         end

         % Creating Faces.
         if newCell.nFaces~=0
            % Indication of the new faces
            ind = newCell.face==0;

            % Create the new faces numbers
            newCell.face(ind)=cells.nFaces+1:cells.nFaces+newCell.nFaces;
         end
         newCell.gtCell = mesh.nCells;
         newCell.gtNode = mesh.nNodes;
         newCell.gtFace = cells.nFaces;

         % Some information from reference cell
         newCell.refBound(1:6)=0;
         refHexa = mesh.cells(refCell,:);
         nodefaces = [6 5 1 2       % sth (1)
                    3 4 8 7       % nth (2)
                    5 8 4 1       % wst (3)
                    2 3 7 6       % est (4)
                    1 4 3 2       % bot (5)
                    6 7 8 5       % top (6)
                    ];
         for face=1:6
            nodes = refHexa(nodefaces(face,:));
            locs=sum(ismember(mesh.surfaces,nodes),2);
            loc=find(locs==4);
            if ~isempty(loc)
               newCell.refBound(face)=loc;
            end
         end
         newCell.ijk = ijk_ref;
         newCell.faceGrow = face+2*mod(face,2)-1;
      end

      % Functions to update the structure.
      function flag = updateIJK(obj,ijk,ncells,nnodes,nfaces)
         % Update this class, incresing one cell.
         obj.ijk(end+1,:)=ijk;
         obj.numCells=ncells;
         obj.numNodes=nnodes;
         obj.numFaces=nfaces;
         flag = true;
      end





      function [newCell] = findExistData(obj,mesh,cells,ijk_ref)
         %FINDEXISTDATA Find the existent data for the new cell.

         % Some important variables.
         local_nodes(1:8) = 0;
         mapEFaces = repelem(1:obj.numCells,diff(cells.mapF2E))';
         mapNFaces = repelem(1:obj.numFaces,diff(cells.mapN2F))';
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
         % struct('faces', 1, 'nodes', 1)
         local_faces(1:6)=0;
         local_coord(1:8,1:3)=0;
         for i=1:6
            if(n_cells(i)>0)
               % Finding the face.
               dir = i+2*mod(i,2)-1;
               faceNum = cells.faces2Elements(mapEFaces==n_cells(i),:);
               faceNum = faceNum(faceNum(:,2)==dir,1);

               % Copy information.
               local_faces(i)=faceNum;
               local_nodes(ref_pos(i,:))=cells.nodes2Faces(mapNFaces==faceNum);
               local_coord(ref_pos(i,:),:) = mesh.coordinates(local_nodes(local_nodes~=0),:);
            end
         end

         % store the information in this object.
         obj.locFaces = local_faces;
         obj.locNodes = local_nodes;
         obj.coordinates = local_coord;
         obj.numNewNodes = sum(local_nodes==0);
         obj.numNewFaces = sum(local_faces==0);

         newCell.face = local_faces;
         newCell.node = local_nodes;
         newCell.coord = local_coord;
         newCell.nNodes = obj.numNewNodes;
         newCell.nFaces = obj.numNewFaces;
      end


      function [newCell] = fillExistData(obj,newCell,dir,len)
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
            ind = obj.locNodes==0;

            % Create the new nodes numbers
            obj.locNodes(ind)=lastNode+1:lastNode+nodes2Add;

            % Create the nodes positions
            obj.coordinates(ind,:) = obj.coordinates(ref(ind),:);
            obj.coordinates(ind,ijk_loc) = obj.coordinates(ind,ijk_loc) + (1-2*mod(dir,2))*len;
         end

         % Creating Faces.
         if faces2Add~=0
            % Indication of the new faces
            ind = obj.locFaces==0;

            % Create the new faces numbers
            obj.locFaces(ind)=lastFace+1:lastFace+faces2Add;
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
                        obj.ijk(ref_neigh(i),2) = obj.ijk(ref_neigh(i),2) - 1;
                     case 2
                        obj.ijk(ref_neigh(i),2) = obj.ijk(ref_neigh(i),2) + 1;
                     case 3
                        obj.ijk(ref_neigh(i),1) = obj.ijk(ref_neigh(i),1) - 1;
                     case 4
                        obj.ijk(ref_neigh(i),1) = obj.ijk(ref_neigh(i),1) + 1;
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

      % Unnecessary function.
      % % % % function constructFreeFaces(obj,face)
      % % % %    %CONSTRUCTFREEFACES Construct the number of free faces by cell.
      % % % %    ncells = length(face.mapF2E)-1;
      % % % % 
      % % % %    faceA = face.faceNeighbors(:,1)==0;
      % % % %    faceB = face.faceNeighbors(:,2)==0;
      % % % %    cells = [face.faceNeighbors(faceA,2);face.faceNeighbors(faceB,1)];
      % % % %    obj.fFacesCells = accumarray(cells,ones(length(cells),1),[ncells 1]);
      % % % % end




      % Functions to update the mesh - OLD WRONG
      % % % % function updateFaces(obj,faces,refCell)
      % % % %    %UPDATEFACES update the face information to add the new cell
      % % % % 
      % % % %    % Mask to some Variables.
      % % % %    faces2Add = obj.numNewFaces;
      % % % %    lastFace = obj.numFaces;
      % % % %    % nodes2Add = obj.numNewNodes;
      % % % %    % lastNode = obj.numNodes;
      % % % % 
      % % % %    % Find the cells around.
      % % % %    [neighcells,~] = obj.findNeighborhod(refCell);
      % % % % 
      % % % %    % Some Important Variables.
      % % % %    % mapEFaces = repelem(1:obj.numCells,diff(faces.mapF2E))';
      % % % %    % mapNFaces = repelem(1:obj.numFaces,diff(faces.mapN2F))';
      % % % %    ref_pos = [6 5 1 2       % sth (1)
      % % % %               3 4 8 7       % nth (2)
      % % % %               5 8 4 1       % wst (3)
      % % % %               2 3 7 6       % est (4)
      % % % %               1 4 3 2       % bot (5)
      % % % %               6 7 8 5       % top (6)
      % % % %               ];
      % % % % 
      % % % %    faceCentroid(1:faces2Add,1:3)=0;
      % % % %    faceNormal(1:faces2Add,1:3)=0;
      % % % % 
      % % % %    % Variables to add with simple construction.
      % % % %    mapN2F(1:faces2Add)=max(faces.mapN2F)+4*(1:faces2Add);
      % % % %    faces2Elements=[obj.faces' (1:6)'];
      % % % %    mapF2E=max(faces.mapF2E)+6;
      % % % % 
      % % % %    % creating information for each face
      % % % %    nodes2Faces(1:4*faces2Add)=0;
      % % % %    faceNeighbors(1:faces2Add,1:2)=obj.numCells+1;
      % % % %    for i=1:faces2Add
      % % % %       faceNum = obj.numFaces+i;
      % % % %       faceNodes = ref_pos(obj.faces == faceNum,:);
      % % % %       faceNodeGL = obj.nodes(faceNodes);
      % % % %       faceCentroid(i,:)=sum(obj.coordinates(faceNodes,:))/4;
      % % % %       direction = faces2Elements(faces2Elements(:,1)==faceNum,2);
      % % % % 
      % % % %       [~,nodeC]=min(faceNodeGL);
      % % % %       nodeL = mod(nodeC - 2, 4) + 1;
      % % % %       nodeR = mod(nodeC, 4) + 1;
      % % % %       nodeA = mod(nodeC + 1, 4) + 1;
      % % % %       if (faceNodeGL(nodeL)<faceNodeGL(nodeR))
      % % % %          nodes2Faces(4*(i-1)+1:4*i)=[
      % % % %             faceNodeGL(nodeC) faceNodeGL(nodeL) ...
      % % % %             faceNodeGL(nodeA) faceNodeGL(nodeR)];
      % % % %          faceNeighbors(i,1)=neighcells(direction);   % <---(i,1) or (i,2) to check
      % % % %       else
      % % % %          nodes2Faces(4*(i-1)+1:4*i)=[
      % % % %             faceNodeGL(nodeC) faceNodeGL(nodeR) ...
      % % % %             faceNodeGL(nodeA) faceNodeGL(nodeL)];
      % % % %          faceNeighbors(i,2)=neighcells(direction);   % <---(i,1) or (i,2) to check
      % % % %       end
      % % % %       v1 = obj.coordinates(faceNodes(nodeL),:)-obj.coordinates(faceNodes(nodeC),:);
      % % % %       v2 = obj.coordinates(faceNodes(nodeA),:)-obj.coordinates(faceNodes(nodeC),:);
      % % % %       % v3 = obj.coordinates(faceNodes(nodeR),:)-obj.coordinates(faceNodes(nodeC),:);
      % % % %       % area = 0.5*norm(cross(v1,v2)) + 0.5*norm(cross(v2,v3));
      % % % %       faceNormal(i,:)=cross(v1,v2);
      % % % %    end
      % % % % 
      % % % %    % UPDATE THE FACE STRUCT.
      % % % %    faces.nodes2Faces = [faces.nodes2Faces; nodes2Faces'];
      % % % %    faces.mapN2F = [faces.mapN2F; mapN2F'];
      % % % %    faces.faces2Elements = [faces.faces2Elements; faces2Elements];
      % % % %    faces.mapF2E = [faces.mapF2E; mapF2E];
      % % % %    faces.faceNeighbors = [faces.faceNeighbors; faceNeighbors];
      % % % %    faces.faceCentroid = [faces.faceCentroid; faceCentroid];
      % % % %    faces.faceNormal = [faces.faceNormal; faceNormal];
      % % % %    faces.nFaces = faces.nFaces+faces2Add;
      % % % % 
      % % % %    facesNot2Add = 6 - faces2Add;
      % % % %    existFaces = faces2Elements(faces2Elements(:,1)<=lastFace,1);
      % % % %    for i=1:facesNot2Add
      % % % %       loc=faces.faceNeighbors(existFaces(i),:)==0;
      % % % %       faces.faceNeighbors(existFaces(i),loc)=obj.numCells+1;
      % % % %    end
      % % % % end
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