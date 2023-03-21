classdef Faces < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    mapN2F
    nodes2Faces
    faceNeighbors
    mapF2E
    faces2Elements
    faceCentroid
    faceNormal
    nFaces
  end
  
  properties (Access = private)
    mesh
  end
  
  methods (Access = public)
    function obj = Faces(msh)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      obj.setFaces(msh);
      obj.setupFaceTopology;
    end
  end
    
  methods (Access = private)
    function setFaces(obj,msh)
      obj.mesh = msh;
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

         faces = repmat(1 : obj.nFaces, [2,1]);
         k     = obj.faceNeighbors(:) ~= 0;

         if nargin == 3
            tags = tags(k);
         end

         faces     = faces(k);
         neigh = obj.faceNeighbors(k);

         obj.mapF2E = zeros(obj.mesh.nCells+1,1);
         obj.faces2Elements = zeros(0,2);
         if nargin == 2
               insertInPackedData(obj, neigh(:), faces(:));
         else
               insertInPackedData(obj, neigh(:), [faces(:), tags(:)]);
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

%        n   = numel(obj.mapF2E)-1;
       tmp = sortrows([r, c]);
       new = tmp(:,2:end);
%        t   = accumarray(r, 1, [n, 1]);
%        num = diff(obj.mapF2E);
%        pos = cumsum([1;double(num)+t]);

%        r   = unique(r);
%        ix  = mcolon(pos(r)+num(r), pos(r+1)-1);
%        newData = zeros(size(obj.faces2Elem, 1)+size(new, 1), size(obj.faces2Elem, 2));
%        i    = (1:size(newData))';
%        i(ix)=[];
%        newData(i, :) = data;
       obj.faces2Elem = new;
%        data = newData;
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

    if numel(n) == 1,
       n = repmat(n, [size(A, dim), 1]);
    end

    assert (all( n(:)>=0 ), 'All repeat counts should be nonnegative.');
    if nargin < 3,
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
  end
end

