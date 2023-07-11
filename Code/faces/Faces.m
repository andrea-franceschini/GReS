classdef Faces < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    nodes2Faces      % List of nodes of each face. The data is arranged vectorwise
    mapN2F           % Indirection map for nodes2Faces    
    faces2Elements   % List of faces of each element. The data is arranged vectorwise
    mapF2E           % Indirection map for faces2Elements
    faceNeighbors    % IDs of the cells to the left and right of each face (nf x 2 matrix)
    faceCentroid     % Centroid coordinates of each face
    faceNormal       % Normal to each face (the magnitude of the vector is equal to the area)
    nFaces           % # of faces
  end
  
  properties (Access = private)
    mesh
    model
  end
  
  methods (Access = public)
    function obj = Faces(simmod,msh)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      obj.setFaces(simmod,msh);
    end
    %     Triangle:
% 
%     v
%     ^
%     |
%     2                 
%     |`\               
%     |  `\             
%     |    `\           
%     |      `\         
%     |        `\       
%     0----------1--> u

%       Quadrangle:
% 
%             v
%             ^
%             |
%       3-----------2       
%       |     |     |       
%       |     |     |       
%       |     +---- | --> u 
%       |           |       
%       |           |      
%       0-----------1     
    
    function areaTri = computeAreaTri(obj,idTri)
      triNod = obj.mesh.surfaces(idTri,:);
%       ptr1 = triNod(2,:);
%       ptr2 = triNod(1,:);
%       v1 = obj.mesh.coordinates(ptr1,:) - obj.mesh.coordinates(ptr2,:);
      v1 = obj.mesh.coordinates(triNod(:,2),:) - obj.mesh.coordinates(triNod(:,1),:);
      v2 = obj.mesh.coordinates(triNod(:,3),:) - obj.mesh.coordinates(triNod(:,1),:);
      areaTri = 0.5*vecnorm(cross(v1,v2,2),2,2);
    end
    
    function areaQuad = computeAreaQuad(obj,idQuad)
      nQuad = length(idQuad);
      quadNod = obj.mesh.surfaces(idQuad,1:4);
      quadNod = reshape(quadNod',[],1);
      surfCentroid = sparse(repelem(1:nQuad,4)', ...
        quadNod,0.25*ones(length(quadNod),1), ...
        nQuad,obj.mesh.nNodes)*obj.mesh.coordinates;
      d = obj.mesh.coordinates(quadNod,:) - repelem(surfCentroid,4,1);
      % Works only for 4-point faces
      ptr = repmat(mod(1:4,4)+1,1,nQuad) + repelem(cumsum([0 4*ones(1,nQuad-1)]),4);
      normalTri = cross(d,d(ptr,:),2);
      areaTri = vecnorm(normalTri,2,2)/2;
%       faceArea = sum(reshape(areaTri,4,[]));
      areaQuad = accumarray(repelem(1:nQuad,4)',areaTri);
    end
  end
    
  methods (Access = private)
    function setFaces(obj,simmod,msh)
      obj.mesh = msh;
      obj.model = simmod;
      %
% % %       if obj.mesh.nSurfaces > 0
% % %         obj.nFaces = obj.mesh.nSurfaces;
% % %         areaSurf = zeros(obj.mesh.nSurfaces,1);
% % %         idTri = (obj.mesh.surfaceVTKType == 5);
% % %         idQuad = (obj.mesh.surfaceVTKType == 9);
% % %         if any(idTri)
% % %           areaSurf(idTri) = obj.computeAreaTri(idTri);
% % %         end
% % %         if any(idQuad)
% % %           areaSurf(idQuad) = obj.computeAreaQuad(idQuad);
% % %         end
% % %         tmpMat = obj.mesh.surfaces';
% % %         topolVec = tmpMat(tmpMat ~= 0);
% % %         [obj.loadedNodesSurf,~,rowID] = unique(topolVec);
% % %         clear matTmp
% % %         colID = repelem(1:obj.mesh.nSurfaces,obj.mesh.surfaceNumVerts);
% % %         areaSurf = areaSurf./obj.mesh.surfaceNumVerts;
% % %         aNod = repelem(areaSurf,obj.mesh.surfaceNumVerts); %./obj.mesh.surfaceNumVerts;
% % %         obj.nodeArea = sparse(rowID,colID,aNod,length(obj.loadedNodesSurf),obj.mesh.nSurfaces);
% % %       end
      %
      if obj.model.isFVTPFABased('Flow')
        obj.setupFaceTopology;
        obj.computeFaceProperties;
      end
    end
        
%         id = 1:obj.mesh.nSurfaces;
%         vTmp = double(obj.mesh.surfaceVTKType == 5) + 4*double(obj.mesh.surfaceVTKType == 9);
%         obj.areaNod = zeros(sum(vTmp),1);
%         obj.mapAreaNod = cumsum([1 vTmp']);
%         idTri = id(vTmp == 1);
%         idQuad = id(vTmp == 4);
%         if ~isempty(idTri)
%           obj.areaNod(obj.mapAreaNod(idTri)) = obj.computeAreaTri(idTri);
%           obj.areaNod(obj.mapAreaNod(idTri)) = (1/3)*obj.areaNod(obj.mapAreaNod(idTri));
%         end
%         if ~isempty(idQuad)
%           vTmp1 = obj.computeAreaQuad(idQuad);
%           obj.areaNod(find(repelem(vTmp,vTmp)==4))=0.25*repelem(vTmp1,4);
%         end
%         obj.nFaces = obj.mesh.nSurfaces;
%       end 
%       %
%       if obj.model.isFVTPFABased('Flow')
%         obj.setupFaceTopology;
%         obj.computeFaceProperties;
%       end
%     end
    
    function computeFaceProperties(obj)
%       % Compute faces center point by averaging
%       nNF = diff(obj.mapN2F);
%       obj.faceCentroid = sparse(repelem(1:obj.nFaces,nNF), ...
%          obj.nodes2Faces,(repelem(nNF,nNF)).^(-1),...
%          obj.nFaces,obj.mesh.nNodes)*obj.mesh.coordinates;
%       %
%       d = obj.mesh.coordinates(obj.nodes2Faces,:) - repelem(obj.faceCentroid,nNF,1);
%       % Works only for 4-point faces
%       ptr = repmat(mod(1:4,4)+1,1,obj.nFaces) + repelem(cumsum([0 nNF(1:end-1)']),4);
%       normalTri = cross(d,d(ptr,:),2);
%       areaTri = vecnorm(normalTri,2,2)/2;
%       faceArea = sum(reshape(areaTri,4,[]));
%       obj.faceNormal = reshape(sum(reshape(normalTri',3,4,[]),2),3,[])';
%       centroidTri = (obj.mesh.coordinates(obj.nodes2Faces,:) + ...
%         obj.mesh.coordinates(obj.nodes2Faces(ptr),:) + repelem(obj.faceCentroid,nNF,1))/3;
%       obj.faceCentroid = reshape(sum(reshape((areaTri.*centroidTri)',3,4,[]),2),3,[])';
%       obj.faceCentroid = obj.faceCentroid.*(1./faceArea)';
%       %
%       obj.faceNormal = obj.faceNormal./(vecnorm(obj.faceNormal,2,2)).*(faceArea');
      % Compute faces center point by averaging
      nNF = diff(obj.mapN2F);
      obj.faceCentroid = sparse(repelem(1:obj.nFaces,nNF), ...
         obj.nodes2Faces,(repelem(nNF,nNF)).^(-1),...
         obj.nFaces,obj.mesh.nNodes)*obj.mesh.coordinates;
      %
      d = obj.mesh.coordinates(obj.nodes2Faces,:) - repelem(obj.faceCentroid,nNF,1);
      % Works only for 4-point faces
      ptr = repmat(mod(1:4,4)+1,1,obj.nFaces) + repelem(cumsum([0 nNF(1:end-1)']),4);
      normalTri = cross(d,d(ptr,:),2);
      areaTri = vecnorm(normalTri,2,2)/2;
%       faceArea = sum(reshape(areaTri,4,[]));
      faceArea = accumarray(repelem(1:obj.nFaces,nNF)',areaTri);
      obj.faceNormal = reshape(sum(reshape(normalTri',3,4,[]),2),3,[])';
      centroidTri = (obj.mesh.coordinates(obj.nodes2Faces,:) + ...
        obj.mesh.coordinates(obj.nodes2Faces(ptr),:) + repelem(obj.faceCentroid,nNF,1))/3;
      obj.faceCentroid = reshape(sum(reshape((areaTri.*centroidTri)',3,4,[]),2),3,[])';
      obj.faceCentroid = obj.faceCentroid.*(1./faceArea);
      %
      obj.faceNormal = obj.faceNormal./(vecnorm(obj.faceNormal,2,2)).*faceArea;

%       nNF = diff(obj.mapN2F);
%       obj.faceCentroid = sparse(repelem(1:obj.nFaces,nNF), ...
%          obj.nodes2Faces,(repelem(nNF,nNF)).^(-1),...
%          obj.nFaces,obj.mesh.nNodes)*obj.mesh.coordinates;
%       %
%       d = obj.mesh.coordinates(obj.nodes2Faces,:) - repelem(obj.faceCentroid,nNF,1);
%       ptr = repmat(mod(1:4,4)+1,1,obj.nFaces) + repelem(cumsum([0 nNF(1:end-1)']),4);
%       normalTri = cross(d,d(ptr,:),2);
%       areaTri = vecnorm(normalTri,2,2)/2;
%       faceArea = sum(reshape(areaTri,4,[]));
%       tmpMat = reshape(normalTri',3,4,[]);
%       tmpMat2 = sum(tmpMat,2);
%       obj.faceNormal = reshape(tmpMat2,3,[])';
%       centroidTri = (obj.mesh.coordinates(obj.nodes2Faces,:) + obj.mesh.coordinates(obj.nodes2Faces(ptr),:) + repelem(obj.faceCentroid,nNF,1))/3;
%       centroidTri = areaTri.*centroidTri;
%       tmpMat = reshape(centroidTri',3,4,[]);
%       tmpMat2 = sum(tmpMat,2);
%       obj.faceCentroid = reshape(tmpMat2,3,[])';
%       obj.faceCentroid = (1./faceArea)'.*obj.faceCentroid;
         
    end
    
    function setupFaceTopology(obj)
      % DISCLAIMER: Thus far hexas are the only allowed elements
      % Implemented following the piece of code in MRST
      % Definition of faces on hex
      bot = obj.mesh.cells(:, [1, 4, 3, 2]);
      top = obj.mesh.cells(:, [6, 7, 8, 5]);
      est = obj.mesh.cells(:, [2, 3, 7, 6]);
      wst = obj.mesh.cells(:, [5, 8, 4, 1]);
      sth = obj.mesh.cells(:, [6, 5, 1, 2]);
      nth = obj.mesh.cells(:, [3, 4, 8, 7]);
      %
      % Rearrange such that for each row, the smalles node index is in the
      % first column
      hfaces   = [wst; est; sth; nth; bot; top];
      [m, i]   = min(hfaces, [], 2); %#ok
      I        = bsxfun(@mod, bsxfun(@plus, (1:4)-1, i-1), 4)+1;
      k        = sub2ind(size(hfaces), repmat((1:obj.mesh.nCells*6)', [1, 4]), I);
      hfaces   = hfaces(k);
      %
      % Rearrange such that for each row, the second smallest node is in the
      % second column.  We use 'i' as the 'sign' of each half-face.
      [m, i] = min(hfaces(:, [2,end]), [], 2); %#ok
      hfaces(i==2, :) = hfaces(i==2, [1,4,3,2]);
      %
      % Better safe than sorry
      [m, k] = min(hfaces(:,[2,end]), [], 2); %#ok
      assert(all(k==1));
      clear k
      %
      % id is [cellnumber, half-face tag]
      id       = [(1:obj.mesh.nCells)', repmat(1, [obj.mesh.nCells, 1]);...
                  (1:obj.mesh.nCells)', repmat(2, [obj.mesh.nCells, 1]);...
                  (1:obj.mesh.nCells)', repmat(3, [obj.mesh.nCells, 1]);...
                  (1:obj.mesh.nCells)', repmat(4, [obj.mesh.nCells, 1]);...
                  (1:obj.mesh.nCells)', repmat(5, [obj.mesh.nCells, 1]);...
                  (1:obj.mesh.nCells)', repmat(6, [obj.mesh.nCells, 1]) ];

      % Sort rows to find pairs of cells sharing a face
      [hfaces, j] = sortrows(hfaces);
      
      % Run length compression to obtain unique face numbering
      [fnodes,n]  = Faces.rlencode(hfaces);
      obj.nodes2Faces = reshape(fnodes', [], 1);
      obj.nFaces          = numel(n);

      % Expand face numbers to each half-face
      N           = Faces.rldecode(1:obj.nFaces,n,2)';

      % Accumulate results in face-2-cell map
      obj.faceNeighbors  = accumarray([N, i(j)], id(j,1), [obj.nFaces, 2]);

      % ... and a map of local directions
      hftag  = accumarray([N, i(j)], id(j,2), [obj.nFaces, 2]);

      % ... and explicitly give number of nodes for each face.  So far, this
      % code is oblivious to pinch which may render faces triangular or
      % collapsed.
      nnodes = repmat(4, [obj.nFaces, 1]);
      obj.mapN2F = [1; cumsum(nnodes)+1];
      
      addFaces(obj, nnodes, hftag);
    end
    
    function addFaces (obj, nnodes, tags)
      %Add faces F from grid structure
      %
      % SYNOPSIS:
      %   H             = addFaces(G, facenodes, numnodes, neighbors, tags)
      %   [H, newfaces] = addFaces(...)
      %
      % PARAMETERS:
      %   G       - Grid structure as described by grid_structure.
      %
      %   fnodes  - face nodes
      %
      %   nnodes  - number of nodes for each face
      %   neigh   - cell neighbors of each face
      %   tags    - cellFaces tags for each half-face associated with face.
      %
      % RETURNS:
      %   G       - Grid structure where the following fields have been modified:
      %
      %                  G.faces.nodes
      %
      %                  G.cells.faces
      %                  G.cells.facePos
      %                  G.cells.numFaces
      %
      %                  G.faces.neighbors
      %                  G.faces.numNodes
      %                  G.faces.nodePos
      %                  G.faces.tag
      %                  G.faces.num
      %
      %  newfaces - Numbering of new faces.
      %
      % COMMENTS:
      %   Cells may not be closed polygons/polyhedra after this call.

      % BUG
      %
      % o Does not preserve orientation in 2D

      %{
      Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

      This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

      MRST is free software: you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation, either version 3 of the License, or
      (at your option) any later version.

      MRST is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      GNU General Public License for more details.

      You should have received a copy of the GNU General Public License
      along with MRST.  If not, see <http://www.gnu.org/licenses/>.
      %}

         assert(numel(nnodes) == size(obj.faceNeighbors, 1));
         assert(sum(nnodes) == numel(obj.nodes2Faces));

%          G.faces.neighbors = [G.faces.neighbors; neighbors];
         %G.faces.tag       = [G.faces.tag; tags];
%          G.faces.nodes     = [G.faces.nodes; facenodes];

%          numf        = numel(nnodes);
%          newfaces    = G.faces.num + (1 : obj.nFaces)';
%          G.faces.num = double(G.faces.num + obj.nFaces);
%          clear numf

%          if isfield(G.faces, 'numNodes')
%             G.faces.numNodes = [G.faces.numNodes; nnodes];
%             G.faces.nodePos = cumsum([1; double(G.faces.numNodes)]);
%          else
%             G.faces.nodePos = cumsum([1; double(diff(G.faces.nodePos)); ...
%                                       double(nnodes)]);
%          end

         faces = repmat((1 : obj.nFaces)', [2,1]);
         k     = obj.faceNeighbors(:) ~= 0;

         if nargin == 3
            tags = tags(k);
         end

         faces     = faces(k);
         neighbors = obj.faceNeighbors;
         neighbors= neighbors(k);
%          neigh = obj.faceNeighbors;

         obj.mapF2E = zeros(obj.mesh.nCells+1,1);
         obj.faces2Elements = zeros(0,2);
         if nargin == 2
           insertInPackedData(obj,neighbors, faces(:));
         else
           insertInPackedData(obj,neighbors, [faces(:), tags(:)]);
         end

         % keep for the time being
%          if isfield(G.faces, 'tag')
%             G.faces.tag = [G.faces.tag; zeros(numel(nnodes), size(G.faces.tag, 2))];
%          end
% 
%          H = G;
    end
    
    function insertInPackedData(obj,r,c)
    %Insert c into row r of packed array (data, pos)
    %
    % SYNOPSIS:
    %   [data, pos] = insertInPackedData(pos, data, r, c)
    %
    % PARAMETERS:
    %   pos, data - Packed data array.
    %   r         - Row indices of new data entries.
    %   c         - New data entries.
    %
    % RETURNS:
    %   data, pos - Updated packed data array.  Note reverse order compared to
    %               input.
    %
    % EXAMPLE:
    %   [G.cells.faces, G.cells.facePos] = ...
    %       insertInPackedData(G.cells.facePos, G.cells.faces, ...
    %                          cells, [faceIDs, hfTags])
    %
    %   [G.faces.nodes, G.faces.nodePos] = ...
    %       insertInPackedData(G.faces.nodePos, G.faces.nodes, faces, newnodes)

    %{
    Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

    This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

    MRST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRST.  If not, see <http://www.gnu.org/licenses/>.
    %}

       n   = numel(obj.mapF2E)-1;
       tmp = sortrows([r, c]);
       new = tmp(:,2:end);
       t   = accumarray(r, 1, [n, 1]);
       num = diff(obj.mapF2E);
       obj.mapF2E = cumsum([1;double(num)+t]);

       r   = unique(r);
       ix  = Faces.mcolon(obj.mapF2E(r)+num(r), obj.mapF2E(r+1)-1);
       newData = zeros(size(obj.faces2Elements, 1)+size(new, 1), size(obj.faces2Elements, 2));
       i    = (1:size(newData))';
       i(ix)=[];
       newData(i, :) = obj.faces2Elements;
       newData(ix,1:size(new,2)) = new;
       obj.faces2Elements = newData;
    end
  end
  
  methods (Static = true)
    function [A, n] = rlencode(A, dim)
    %Compute run length encoding of array A along dimension dim.
    %
    % SYNOPSIS:
    %   [A, n] = rlencode(A)
    %   [A, n] = rlencode(A, dim)
    %
    % PARAMETERS:
    %   A         - Array.  Numeric or cell array of string ("cellstring").
    %   dim       - dimension of `A` where run length encoding is done.
    %               `dim > 0`, `dim <= ndims(A)`.
    %               OPTIONAL.  Default value: `dim=1`.
    %
    % RETURNS:
    %   A         - Compressed `A` where repeated layers are removed.
    %   n         - repetition count of repeated layers in original `A`.
    %
    % EXAMPLE:
    %   % 1) Regular numeric matrix
    %   A      = [ 1, 2, 3, 4 ; ...
    %              1, 2, 3, 4 ; ...
    %              3, 4, 5, 6 ; ...
    %              3, 3, 3, 3 ; ...
    %              3, 3, 4, 5 ; ...
    %              3, 3, 4, 5 ];
    %
    %   [A, n] = rlencode(A, 1);
    %
    %   assert (isequal(A, [ 1, 2, 3, 4 ; ...
    %                        3, 4, 5, 6 ; ...
    %                        3, 3, 3, 3 ; ...
    %                        3, 3, 4, 5 ]))
    %
    %   assert (isequal(n, [2, 1, 1, 2] .'))
    %
    %   % 2) Cell array of string
    %   S = { 'aaaa' ; 'aaaa' ; 'bb' ; 'cccc' ; 'cccc' ; 'd'; 'd'; 'd' ; 'd' };
    %   [s, n] = rlencode(S);
    %
    %   assert (isequal(s, { 'aaaa' ; 'bb'; 'cccc'; 'd' }))
    %   assert (isequal(n, [ 2; 1; 2; 4 ]))
    %
    % SEE ALSO:
    %   `rldecode`.

    %{
    Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

    This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

    MRST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRST.  If not, see <http://www.gnu.org/licenses/>.
    %}

     if isempty(A)
        n = 0;
        return
     end

     if nargin < 2
        dim = 1;
     end

     % Pick out (1:end - 1, :, ...) and (2:end, :, ...) in a multidimensional
     % way.  In particular, swap dimensions 1 and 'dim' so that we only have
     % to consider compressing along the first dimension.
     d = 1 : ndims(A) ;   d([1, dim]) = [dim, 1];
     B = permute(A, d);

     % Find positions where layers differ.
     i = [find(any(Faces.different_layers(B), 2)) ; size(B, 1)];

     % Compare differences in position to find run length.
     n = diff([0; i]);

     % Re-swap dimensions 1 and 'dim' to form return value.
     sz = size(B); sz(1) = numel(i);
     A  = permute(reshape(B(i,:), sz), d);
    end
    %----------------------------------------------------%
    function d = different_layers(B)
     top = B(1 : (end - 1), :);
     bot = B(2 : (end - 0), :);

     if isnumeric(B)
        d = xor(top ~= bot, isnan(top) & isnan(bot));
     else
        % Handles case of cellstring and, in MATLAB >= R2016b, STRING too.
        d = ~strcmp(top, bot);
     end
    end
    
    function A = rldecode(A, n, dim)
    %Decompress run length encoding of array `A` along dimension `dim`.
    %
    % SYNOPSIS:
    %   B = rldecode(A, n, dim)
    %   B = rldecode(A, n) % dim assumed to be 1
    %
    % PARAMETERS:
    %   A         - encoded array
    %   n         - repetition of each layer along dimension `dim`. `n` can be
    %               either a scalar or one repetition number for each layer.
    %   dim       - dimension of `A` where run length encoding is done.
    %               dim > 0.
    %
    % RETURNS:
    %   B         - uncompressed `A`
    %
    % EXAMPLE:
    %   % 1) Numerical example:
    %   A = [1,2,3,4;1,2,3,4;3,4,5,6;3,3,3,3;3,3,4,5;3,3,4,5]
    %   [B,n] = rlencode(A,1)
    %   C = rldecode(B,n,1)
    %
    %   % 2) Retrive 'first' column of G.cells.faces (see grid_structure):
    %   G = cartGrid([10, 10, 2]);
    %   cellNo  = rldecode(1:G.cells.num, double(G.cells.numFaces), 2) .';
    %   disp(['CellFace nr. 10 belongs to cell nr: ', num2str(cellNo(10))]);
    %
    % SEE ALSO:
    %   `rlencode`

    %{
    Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

    This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

    MRST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRST.  If not, see <http://www.gnu.org/licenses/>.
    %}


    if nargin < 3
      dim = 1;
    end

    assert(dim > 0, 'Third argument DIM must be positive');

    if numel(n) == 1
       n = repmat(n, [size(A, dim), 1]);
    end

    assert (all( n(:)>=0 ), 'All repeat counts should be nonnegative.');
    if nargin < 3
       assert (numel(n) == size(A, dim), ...
       sprintf(['There should be a repeat count for each value along dimension dim.\n',...
        'The default value of dim is 1. Did you forget to specify dim?']));
    else
       assert (numel(n) == size(A, dim), ...
       'There should be a repeat count for each value along dimension dim.');
    end

    % Take dimension we compress along to be first dimension,
    % i.e., swap dimensions 1 and dim.
    d      = 1:max(dim, ndims(A));
    d([1, dim])   = [dim, 1];
    B      = permute(A,d);

    r      = n(:)~=0;
    B      = reshape(B(r, :), sum(r), []);

    % Insert repeated layers and permute back
    i      = cumsum([1; double(reshape(n(r), [], 1))]);
    j      = zeros(i(end)-1,1);
    j(i(1:end-1)) = 1;

    szA    = [size(A), ones(1, dim-ndims(A))];
    A      = permute(reshape(B(cumsum(j),:), [sum(n(:)), szA(d(2:end))]), d);
    end
    
    function x = mcolon(lo, hi, s)
%Compute vector of consecutive indices from separate lower/upper bounds.
%
% SYNOPSIS:
%   ind = mcolon(lo, hi)
%   ind = mcolon(lo, hi, stride)
%
% PARAMETERS:
%   lo  - Vector of start values (lower bounds).
%   hi  - Vector of end values (upper bounds).
%   s   - Vector of strides.
%         Optional.  Default value: s = ones([numel(lo), 1]) (unit stride).
%
% RETURNS:
%   ind - `[lo(1):hi(1)     , lo(2):hi(2)     ,..., lo(end):hi(end)]`
%   ind - `[lo(1):s(1):hi(1), lo(2):s(2):hi(2),...,lo(end):s(end):hi(end)]`
%
%
% EXAMPLE:
%   lo  = [1 1 1 1]; hi = [2 3 4 5];
%   ind = mcolon(lo, hi)
%
% NOTES:
%   Note that `ind` has type `double` irrespective of the type of its input
%   parameters.
%
%   `mcolon` may be implemented in terms of `arrayfun` and `horzcat`, e.g., ::
%
%      ind = arrayfun(@colon, lo, hi, 'UniformOutput', false);
%      ind = [ind{:}];
%   or ::
%
%      ind = arrayfun(@colon, lo, s, hi, 'UniformOutput', false);
%      ind = [ind{:}];
%   but the current implementation is faster.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


% Support general input class.
lo = double(lo);
hi = double(hi);

if numel(lo) == 1
   lo = repmat(lo, size(hi));
end
if numel(hi)==1
   hi = repmat(hi, size(lo));
end

if numel(lo) ~= numel(hi)
  error('In mcolon: lo and hi must have same number of elements!');
elseif numel(lo) == 0
  % Empty index boundaries -> empty result.
  x = [];
elseif nargin < 3
  % x = mcolon(lo, hi), stride = 1.

  % Remove lo-hi pairs where numel(lo:hi)==0
  i    = hi >= lo;

  if ~ any(i), x = []; return, end

  hi   = hi(i);
  lo   = lo(i);
  m    = numel(lo);           % Number of bins
  d    = double(hi - lo + 1); % Number of elements per bin (hi/lo included)
  n    = sum(d);              % Total number of (expanded) elements

  % Preallocate result.  Recall: CUMSUM(ONES([1,n])) == 1:n
  x    = ones([1, n]);
  x(1) = lo(1); % Starting index of first bin in first position.

  % Prepare for running sum accumulation.  The following properties are
  % used
  %
  %   1) 1 + cumsum(d(1:end-1)) is starting positions of bins 2:m
  %   2) lo(2:m) - hi(1:m-1) is offset from hi(1:m-1) to reach lo(2:m)
  %
  x(1 + cumsum(d(1:end-1))) = lo(2:m) - hi(1:m-1);

  % Alternative setups to previous two statements.
  %   x(cumsum([1 , d(1:end-1)])) = [ lo(1) , lo(2:m) - hi(1:m-1) ]
  %   x(cumsum([1 , d(1:end-1)])) = lo - [ 0 , hi(1:m-1) ]

  % Compute result by running sum of index accumulators (mostly ONES).
  %
  % Note: Running sum is adjusted (i.e., reset to lo(i)) at each bin
  % boundary 'i' by 2) above.
  x = cumsum(x);
else
  % x = mcolon(lo,hi,s), stride = s.

  s = double(s);
  if numel(s) == 1
     s = repmat(s, size(lo));
  end

  % Remove lo-hi-s triplets where numel(lo:s:hi)==0
  i    = ((hi >= lo) & (s > 0)) | ((hi <= lo) & (s < 0));
  hi   = hi(i);
  lo   = lo(i);
  s    = s(i);

  if sum(i) == 0, x=[];return; end

  % Compute lo + (0:(hi-lo)/stride)*stride
  % Fix or hack: avoid roundoff error in floor((hi-lo)./s) when hi-lo = N*s
  % for some natural number N.
  e    =  (1- 2*(hi<lo)) * eps;
  hi   = fix((e+hi-lo)./s);


  m    = numel(lo);
  d    = double(hi + 1);
  n    = sum(d);

  assert (all(d > 0), ...
         ['Internal error in ''%s'': Bins with non-positive ', ...
          'number of elements detected'], mfilename);

  ind = 1+cumsum(d(1:end-1));

  % Expand lo  to [lo(1) lo(1) ... lo(2) lo(2) ... lo(end)]
  LO   = zeros(1,n);
  LO(1)=lo(1);
  LO(ind) = lo(2:m) - lo(1:m-1);
  LO   = cumsum(LO);

  % Expand stride
  S    = zeros(1,n);
  S(1) = s(1);
  S(ind) = s(2:m) - s(1:m-1);
  S    = cumsum(S);

  x    = ones(1,n);
  x(1) = 0;
  x(ind) = -hi(1:m-1);

  x    = cumsum(x).*S + LO;
end
end
  end
end

