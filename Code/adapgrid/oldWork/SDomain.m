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


      function [cell] = newCellData(obj,mesh,cells,ijk_ref,dir,len,refCell)
        %NEWCELLDATA Create the faces information for the new cell.
        %
        % output:
        %    cell ->  segregate the information 
        %    cell.new -> the information necessary to 
        %  cell.point -> relate to 
        %   cell.surf ->
        %   cell.volu ->
        %       newCell.face(6) -> face number
        %       newCell.node(8) -> node number
        %    newCell.coord(8,3) -> nodes coordenates
        %     newCell.nNodes(1) -> number of new nodes
        %     newCell.nFaces(1) -> number of new faces
        %     newCell.gtCell(1) -> number of cell in the domain
        %     newCell.gtNode(1) -> number of nodes in the domain
        %     newCell.gtFace(1) -> number of faces in the domain
        %   newCell.refBound(6) ->
        %        newCell.ijk(3) -> ijk of the new cell
        %   newCell.faceGrow(1) -> direction of grow in reference of the face
        %  newCell.bordRef(2,6) -> cell and direction to be used as
        %  reference. The first line refering to the cell and the second
        %  to the face axis.

        % NODE ORDERING ASSUMPTION:
        % Hexahedron (8 nodes):
        %
        %        Z
        % 5----------8
        % |\     ^   |\
        % | \    |   | \
        % |  \   |   |  \
        % |   6------+---7
        % |   |  +-- |-- | -> Y
        % 1---+---\--4   |
        %  \  |    \  \  |
        %   \ |     \  \ |
        %    \|      X  \|
        %     2----------3
        % map face
        % 1(sth): 6-5-1-2
        % 2(nth): 3-4-8-7
        % 3(est): 5-8-4-1
        % 4(wst): 2-3-7-6
        % 5(bot): 1-4-3-2
        % 6(top): 5-6-7-8


        % cell.new.ijk(3) -> ijk of the new cell
        % cell.new.nNodes(1)-> number of new nodes
        % cell.new.nFaces(1)-> number of new faces
        % cell.new.faceGrow(1) -> direction of grow in reference of the face

        % cell.point.node(8)    -> node number
        % cell.point.coord(8,3) -> nodes coordenates
        % cell.point.gtNode(1)  -> number of nodes in the domain

        % cell.surf.face(6) -> face number
        % cell.surf.gtFace(1) -> number of faces in the domain
        % cell.surf.bordRef(2,6) -> cell and direction to be used as
        %  reference. The first line refering to the cell and the second
        %  to the face axis.

        % cell.volu.gtCell(1) -> number of cell in the domain
        % 
        




        cell = struct('new',[],'point',[],'surf',[],'volu',[]);
        cell.new.ijk= ijk_ref;
        cell.new.nNodes = [];
        cell.new.nFaces = [];
        cell.new.faceGrow = [];
        cell.new.cellTag = 1;

        cell.point.node(1:8) = 0;
        cell.point.coord(1:8,1:3) = 0.;
        cell.point.gtNode = mesh.nNodes;

        cell.surf.face(1:6) = 0;
        cell.surf.gtFace = cells.nFaces;
        cell.surf.bordRef = zeros(2,6);
        cell.surf.area(1:6) = 0;
        cell.surf.normal(1:6,1:3) = 0;
        cell.surf.centroid(1:6,1:3) = 0;

        cell.volu.gtCell = mesh.nCells;
        cell.volu.volume = 0.;
        cell.volu.centroid(1:3) = 0.;

        % Some important variables.
        mapEFaces = repelem(1:obj.numCells,diff(cells.mapF2E))';

        % bot = [1 2 3 4];
        % top = [5 6 7 8];
        % est = [2 6 7 3];
        % wst = [1 5 8 4];
        % sth = [1 2 6 5];
        % nth = [4 3 7 8];
        % ref_nodes = [wst; est; sth; nth; bot; top];

        % Original
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
        for face=1:6
          if(n_cells(face)>0)
            % Finding the face.
            faceMirror = face+2*mod(face,2)-1;
            faceNode = mesh.cells(n_cells(face),ref_nodes(faceMirror,:));
            cell.point.node(ref_nodes(face,:)) = faceNode;

            faceNum = cells.faces2Elements(mapEFaces==n_cells(face),:);
            faceNum = faceNum(faceNum(:,2)==faceMirror,1);
            cell.surf.face(face)=faceNum;
          end
        end
        pos = find((cell.point.node~=0)==1);
        cell.point.coord(pos,:)=mesh.coordinates(cell.point.node(pos),:);
        cell.new.nNodes = sum(cell.point.node==0);
        cell.new.nFaces = sum(cell.surf.face==0);

        % Direction of the Cell Grow.
        ijk_loc = mod(ceil(dir/2)-1, 3) + 1;

        % Creating Nodes and Coordenates.
        if cell.new.nNodes~=0
          % Nodes that where copy from the others in the direction.
          % ref = [
          %   1 5 8 4 2 6 7 3;
          %   4 3 7 8 1 2 6 5;
          %   5 6 7 8 1 2 3 4 ];

          % Original
          ref = [ 4 3 7 8 1 2 6 5;
            1 5 8 4 2 6 7 3;
            5 6 7 8 1 2 3 4 ];
          ref = ref(ijk_loc,:);

          % Indication of the new nodes
          ind = cell.point.node==0;

          % Create the new nodes numbers
          cell.point.node(ind)=mesh.nNodes+1:mesh.nNodes+cell.new.nNodes;

          % Create the nodes positions
          cell.point.coord(ind,:) = cell.point.coord(ref(ind),:);
          cell.point.coord(ind,ijk_loc) = cell.point.coord(ind,ijk_loc) + (1-2*mod(dir,2))*len;
        end
        cell.volu.centroid = mean(cell.point.coord);

        % Creating Faces.
        if cell.new.nFaces~=0
          % Indication of the new faces
          ind = cell.surf.face==0;

          % Create the new faces numbers
          cell.surf.face(ind)=cells.nFaces+1:cells.nFaces+cell.new.nFaces;
        end

        % -----------------------------------------------------------------
        % Some information from reference cell
        refHexa = mesh.cells(refCell,:);

        % top = [5 6 7 8];
        % est = [2 6 7 3];
        % wst = [1 5 8 4];
        % sth = [1 2 6 5];
        % nth = [4 3 7 8];
        % bot = [1 2 3 4];

        % bot = [1 4 3 2];
        % top = [6 7 8 5];
        % est = [2 3 7 6];
        % wst = [5 8 4 1];
        % sth = [6 5 1 2];
        % nth = [3 4 8 7];

        % nodefaces = [wst; est; sth; nth; bot; top];
        
        % Original
        bot = [1 4 3 2];
        top = [6 7 8 5];
        est = [2 3 7 6];
        wst = [5 8 4 1];
        sth = [6 5 1 2];
        nth = [3 4 8 7];
        % nodefaces = [6 5 1 2       % sth (1)
        %   3 4 8 7       % nth (2)
        %   5 8 4 1       % wst (3)
        %   2 3 7 6       % est (4)
        %   1 4 3 2       % bot (5)
        %   6 7 8 5       % top (6)
        %   ];

        nodefaces = [sth; nth; wst; est; bot; top];

        cell.new.faceGrow = face+2*mod(face,2)-1;

        % -----------------------------------------------------------------
        % Compute faces areas, centroid and normal
        for face=1:6
          nodes = nodefaces(face,:);
          coords = cell.point.coord(nodes,:);
          centerPt = sum(coords)/4;
          centerTri(1:4,1:3)=0;
          surfTri(1:4)=0;
          for j=1:4
            centerTri(j,:)=sum([ centerPt; coords(j,:); coords(mod(j,4)+1,:) ])/3;
            surfTri(j)=norm(0.5*cross( (coords(j,:)-centerPt),...
              (coords(mod(j,4)+1,:)-centerPt)));
          end
          cell.surf.area(face)=sum(surfTri);
          cell.surf.centroid(face,:)=sum(surfTri'.*centerTri)/cell.surf.area(face);

          v1 = coords(2,:)-coords(1,:);
          v2 = coords(3,:)-coords(1,:);
          v2xv1 = cross(v1,v2);
          cell.surf.normal(face,:) = v2xv1/norm(v2xv1);
        end

        % -----------------------------------------------------------------
        % Compute the Volume
        tets = [ 1 2 4 5;
          2 3 4 7;
          2 5 6 7;
          4 5 7 8;
          2 4 5 7 ];
        for i = 1:size(tets,1)
          coords = cell.point.coord(tets(i,:),:);
          A = coords(1,:);
          B = coords(2,:);
          C = coords(3,:);
          D = coords(4,:);

          % Compute volume of tetrahedron
          tet_vol = abs(dot(B - A, cross(C - A, D - A))) / 6;
          cell.volu.volume = cell.volu.volume + tet_vol;
        end

        % -----------------------------------------------------------------

        cell.surf.bordRef(1,:)=refCell;
        cell.surf.bordRef(2,:)=dir;
        [n_cells,~] = obj.findNeighborhod(ijk_ref);
        for i=1:6
          if n_cells(i) ~= 0
            cell.surf.bordRef(1,i)=n_cells(i);
            cell.surf.bordRef(2,i)=i+2*mod(i,2)-1;
          end
        end
        [n_cells,~] = obj.findNeighborhodOfCell(refCell);
        for i=1:6
          if n_cells(i) ~= 0
            cell.surf.bordRef(1,i)=n_cells(i);
            cell.surf.bordRef(2,i)=i;
          end
        end
        cell.surf.bordRef(:,dir+2*mod(dir,2)-1)=0;


        [n_cells,~] = obj.findNeighborhod(ijk_ref);
        tl = HexaCell(mesh,cells,n_cells,dir,len);
      end


      function updateBorder(obj,data,cell,cellRef)
        keys = obj.boundaries.keys;
        dataRef = obj.findFacesCells(data.faces,obj.grid(cellRef,:));
        dataNew = obj.findFacesCells(data.faces,cell.new.ijk);
        for i=1:length(keys)
          ref = obj.boundaries(keys{i});
          ref.update(dataNew,dataRef,cell);
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
         % Update this class, incresing one cell.
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
                     case 1
                        obj.grid(ref_neigh(i),2) = obj.grid(ref_neigh(i),2) - 1;
                     case 2
                        obj.grid(ref_neigh(i),2) = obj.grid(ref_neigh(i),2) + 1;
                     case 3
                        obj.grid(ref_neigh(i),1) = obj.grid(ref_neigh(i),1) - 1;
                     case 4
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
      % Some moviment in the grid
      % function ijk_loc = cell_IJK_from_Ref(ijk_neig,dir)
      %    %IJK_LOC Function to find cell ijk to where to grow
      %    switch dir
      %       case 1
      %          ijk_loc = ijk_neig(dir,:)+[1 0 0];
      %       case 2
      %          ijk_loc = ijk_neig(dir,:)-[1 0 0];
      %       case 3
      %          ijk_loc = ijk_neig(dir,:)+[0 1 0];
      %       case 4
      %          ijk_loc = ijk_neig(dir,:)-[0 1 0];
      %       case 5
      %          ijk_loc = ijk_neig(dir,:)+[0 0 1];
      %       case 6
      %          ijk_loc = ijk_neig(dir,:)-[0 0 1];
      %    end
      % end

      % function coord = growDirection(pos,dir,len)
      %    coord = pos;
      %    switch dir
      %       case 1
      %          coord(1)=coord(1)-len;
      %       case 2
      %          coord(1)=coord(1)+len;
      %       case 3
      %          coord(2)=coord(2)-len;
      %       case 4
      %          coord(2)=coord(2)+len;
      %       case 5
      %          coord(3)=coord(3)-len;
      %       case 6
      %          coord(3)=coord(3)+len;
      %    end
      % end

   end



end