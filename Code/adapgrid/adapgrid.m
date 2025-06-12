classdef adapgrid < handle
   % ADAPGRID Class to control the mesh grow

   properties
      ijk         % ijk information for each cell
      fFacesCells % number of free faces for each cell
   end

   methods(Access = public)
      function obj = adapgrid(grid)
         %ADAPGRID Construct an instance of this class
         obj.constructIJK(grid.faces,1,[0 0 0]);
         obj.constructFreeFaces(grid.faces);
      end

   end

   methods (Access = private)
      % 
      % function outputArg = method1(obj,inputArg)
      %    %METHOD1 Summary of this method goes here
      %    %   Detailed explanation goes here
      %    outputArg = obj.Property1 + inputArg;
      % end

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


   end
end