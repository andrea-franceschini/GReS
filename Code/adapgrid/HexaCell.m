classdef HexaCell < handle
  % HEXACELL cell class to construct the adaptative mesh.
  %
  % NODE ORDERING ASSUMPTION (same as Gmsh output):
  % Hexahedron (8 nodes):
  %        v
  % 4----------3
  % |\     ^   |\
  % | \    |   | \
  % |  \   |   |  \
  % |   8------+---7
  % |   |  +-- |-- | -> u
  % 1---+---\--2   |
  %  \  |    \  \  |
  %   \ |     \  \ |
  %    \|      w  \|
  %     5----------6

  % map face
  % 1: 1-2-3-4
  % 2: 1-4-5-8
  % 3: 1-2-5-6
  % 4: 2-3-6-7
  % 5: 3-4-7-8
  % 6: 5-6-7-8

  properties(Access = public)
    coord (8,3) = 0
    nodes (1,8) = 0
    faces (1,6) = 0

    centroid (1,3) = 0.
    faceArea (1,6) = 0.
    faceNormal (6,3) = 0.
    faceCentroid (6,3) = 0.
    volume (1,1) = 0.

    lastNode = 0
    lastFace = 0
    lastCell = 0

    newNodes = 0
    newFaces = 0

    cellNeigh (6,1) = zeros(6,1)
    facesNew (6,1) = true(6,1)

    tag
    faceGrow
  end

  properties (Constant)
    % FaceOrder = ["wst"; "est"; "sth"; "nth"; "bot"; "top"]
    FacesToNodes = [
      1 5 8 4
      2 3 7 6
      1 2 6 5
      4 8 7 3
      1 4 3 2
      5 6 7 8 ]
    swapNodesBMirrorFaces = [1 4 3 2]
  end


  methods (Access = public)
    function obj = HexaCell(mesh,faces,neigh,direction,len,tag)
      %HEXACELL Construct an instance of this class
      obj.lastCell = mesh.nCells;
      obj.lastFace = faces.nFaces;
      obj.lastNode = mesh.nNodes;
      obj.faceGrow = direction+2*mod(direction,2)-1;

      % Obtains new cell information from neighboring cells.
      obj.copyExistData(mesh,faces,neigh); % --- Parece OK.

      % Create the coordinates for the missing nodes.
      obj.createCoords(direction,len);      

      % Creating Faces.
      if obj.newFaces~=0
        % Indication of the new faces
        ind = obj.faces==0;

        % Create the new faces numbers
        obj.faces(ind)=obj.lastFace+1:obj.lastFace+obj.newFaces;
      end

      % Compute faces areas, centroid and normal
      obj.createFaceData()

      % Compute the volume of the cell.
      obj.computeVolume()

      obj.cellNeigh = neigh;
      obj.facesNew = obj.faces > obj.lastFace;
      obj.tag = tag;
    end

    function out = faces2Elements(obj)
      out = [obj.faces' (1:6)'];
    end

    function out = newNodes2Faces(obj)
      out=zeros(4*obj.newFaces,1);
      pos=1;
      for faceId=1:6
        if obj.facesNew(faceId)
          out(pos:pos+3) = obj.nodes(obj.FacesToNodes(faceId,:));
          pos=pos+4;
        end
      end
    end

    function out = mapN2F(obj)
      out = 4*(1:obj.newFaces);
    end

    function out = mapF2E(~)
      out = 6;
    end

    function out = faceNeighbors(obj,lnkFace)
      tmpNeigh(1:6,1:2)=obj.lastCell+1;
      tmpNeigh(obj.facesNew,2)=0;
      tmpNeigh=tmpNeigh(obj.facesNew,:);      
      out=[lnkFace.faceNeighbors; tmpNeigh];
      for i=1:6
        if ~obj.facesNew(i)
          loc=out(obj.faces(i),:)==0;
          out(obj.faces(i),loc)=obj.lastCell+1;
        end
      end
    end

    function out = coordinates(obj)
      newPts = obj.nodes>obj.lastNode;
      out = obj.coord(newPts,:);
    end

    function out = cellTag(obj)
      out = obj.tag;
    end 

    function out = cellNumVerts(~)
      out = 8;
    end

    function out = cellVTKType(~)
      out = 12;
    end

    function out = computeTrans(obj,KMat)
      % Compute the half transmissibility of a hexacell
      r = [1, 1, 1, 2, 2, 2, 3, 3, 3];
      c = [1, 2, 3, 1, 2, 3, 1, 2, 3];
      L = obj.faceCentroid-obj.centroid;
      N = obj.faceArea'.*obj.faceNormal;      
      hT = zeros(6,1);
      for k=1:9
        hT = hT + L(:,r(k)) .* KMat(k) .* N(:,c(k));
      end
      out = hT./sum(L.*L,2);
    end

    function out = newNodesMirrorNodes(obj)
      ijk_dir = mod(ceil(obj.faceGrow/2)-1, 3) + 1;

      out = [];
      % Creating Nodes and Coordinates.
      if obj.newNodes~=0
        % New nodes
        ind = obj.nodes>obj.lastNode;
        % Nodes that where copy from the others in the direction.
        nodesMirrorByAxis = [
          2 1 4 3 6 5 8 7;
          4 3 2 1 8 7 6 5;
          5 6 7 8 1 2 3 4];
        ref = nodesMirrorByAxis(ijk_dir,:);

        % Copy the nodes
        out = obj.nodes(ref(ind))';
      end
    end

  end

  methods (Access = private)

    function copyExistData(obj,mesh,faces,neigh)
      for faceId=1:6
        if (neigh(faceId)>0)
          loopCellNeigh = neigh(faceId);
          loopFaceNeigh = faceId+2*mod(faceId,2)-1;
          loopNeighFacesDir = faces.faces2Elements(faces.mapF2E(loopCellNeigh):faces.mapF2E(loopCellNeigh+1)-1,:);
          loopCommonFace = loopNeighFacesDir(loopNeighFacesDir(:,2)==loopFaceNeigh,1);
          loopNodePos = obj.FacesToNodes(faceId,:);
          loopNodePos = loopNodePos(obj.swapNodesBMirrorFaces);
          loopNeighNodePos = obj.FacesToNodes(loopFaceNeigh,:);
          locNode = mesh.cells(loopCellNeigh,loopNeighNodePos);
          locCoord = mesh.coordinates(locNode,:);

          % copy the information
          obj.nodes(loopNodePos)=locNode;
          obj.faces(faceId)=loopCommonFace;
          obj.coord(loopNodePos,:)=locCoord;
        end
      end
      obj.newNodes = sum(obj.nodes==0);
      obj.newFaces = sum(obj.faces==0);
    end

    function createCoords(obj,direction,len)
      % ijk_dir = mod(ceil(direction/2)-1, 3) + 1;
      ijk_dir = mod(ceil(obj.faceGrow/2)-1, 3) + 1;

      % Creating Nodes and Coordinates.
      if obj.newNodes~=0
        % Create the new nodes
        ind = obj.nodes==0;
        obj.nodes(ind)=obj.lastNode+1:obj.lastNode+obj.newNodes;

        % Nodes that where copy from the others in the direction.
        nodesMirrorByAxis = [
          2 1 4 3 6 5 8 7;
          4 3 2 1 8 7 6 5;
          5 6 7 8 1 2 3 4];
        ref = nodesMirrorByAxis(ijk_dir,:);

        % Create the nodes positions
        obj.coord(ind,:) = obj.coord(ref(ind),:);
        % obj.coord(ind,ijk_dir) = obj.coord(ind,ijk_dir) + (1-2*mod(direction,2))*len;
        obj.coord(ind,ijk_dir) = obj.coord(ind,ijk_dir) + (2*mod(obj.faceGrow,2)-1)*len;
      end
      obj.centroid = mean(obj.coord);
    end

    function createFaceData(obj)
      for faceId=1:6
        locNode = obj.FacesToNodes(faceId,:);
        locCoord = obj.coord(locNode,:);
        locCenter = sum(locCoord)/4;
        locFaceCenter(1:4,1:3)=0;
        locFaceArea(1:4)=0;
        for j=1:4
          locFaceCenter(j,:)=sum([ locCenter; locCoord(j,:); locCoord(mod(j,4)+1,:) ])/3;
          locFaceArea(j)=norm(0.5*cross( (locCoord(j,:)-locCenter),...
            (locCoord(mod(j,4)+1,:)-locCenter)));
        end
        obj.faceArea(faceId)=sum(locFaceArea);
        obj.faceCentroid(faceId,:)=sum(locFaceArea'.*locFaceCenter)/obj.faceArea(faceId);

        v1 = locCoord(2,:)-locCoord(1,:);
        v2 = locCoord(3,:)-locCoord(1,:);
        v2xv1 = cross(v1,v2);
        obj.faceNormal(faceId,:) = v2xv1/norm(v2xv1);
      end
    end

    function computeVolume(obj)
      tets = [ 1 2 4 5;
        2 3 4 7;
        2 5 6 7;
        4 5 7 8;
        2 4 5 7 ];
      for i=1:5
        locCoord = obj.coord(tets(i,:),:);
        A = locCoord(1,:);
        B = locCoord(2,:);
        C = locCoord(3,:);
        D = locCoord(4,:);

        % Compute volume of tetrahedron
        loopVol = abs(dot(B - A, cross(C - A, D - A))) / 6;
        obj.volume = obj.volume + loopVol;
      end
    end



  end

  methods (Static)
  end


end

