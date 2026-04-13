classdef ArrayOfArrays < handle
  % class implementing management for Array of arrays of different sizes
  % used for connectivities

  % under development

  properties (Access = private)
    data      % (1 x N)
    ptr       % (1,nr+1)
    nrows = 0;
  end

  % Construction
  %
  %   A = ArrayOfArrays(cellArray)
  %       Build from a cell array of row vectors, e.g.
  %           ArrayOfArrays({[1 2 3], [4 5], [6 7 8 9]})
  %
  %   A = ArrayOfArrays(matrix)
  %       Build from a dense numeric matrix where each row is zero-padded on
  %       the right to the width of the longest row.  Trailing zeros are
  %       stripped; leading/interior zeros are kept as real values.
  %           conn = [3 7 9 0;
  %                   1 4 6 8;
  %                   2 5 0 0];
  %           A = ArrayOfArrays(conn)
  %
  %
  %   A = ArrayOfArrays(data, sizes)
  %       Build from a pre-flattened data vector and a sizes vector
  %       (length of each row). Sizes and ptr are distinguished by length:
  %           numel(sizes) == nrows   → treated as sizes
  %           numel(ptr)   == nrows+1 → treated as ptr
  %
  % Indexing (operator overloading)
  %   A(r, :)    – returns the full r-th row as a 1-D array
  %   A(r, j)    – returns the j-th element of row r
  %   A(r, j:k)  – returns a slice of row r
  %
  % Other methods
  %   len  = A.arraySize(r)   – length of row r
  %   lens = A.arraySize()    – lengths of all rows in a flat array


  methods

    function obj = ArrayOfArrays(varargin)

      if nargin == 0
        return
      end

      if nargin == 1

        if strcmp(class(varargin{1}),"ArrayOfArrays")
          obj = varargin{1};
          return
        end

        mat = varargin{1};

        assert(isnumeric(mat) && ismatrix(mat),"Input must be a matrix")
        % remove continguous zeros from the right
        id = mat > 0;
        rl = sum(id,2);
        mat = mat';
        id = id';
        flat = mat(id(:));
      else
        flat = varargin{1};
        rl = varargin{2};
      end

      if length(flat)~=sum(rl)
        error("Lenght of flat data and row length is not matching")
      end

      obj.nrows = numel(rl);
      rl = reshape(rl,1,[]);
      obj.ptr   = [1, cumsum(rl) + 1];
      obj.data  = flat;

    end


    function [data,ptr] = getData(obj)

      % return flat data and indirection map

      data = obj.data(:);

      if nargout > 1
        ptr = obj.ptr(:);
      end


    end



    function out = get(obj, r, c)


      if ~isscalar(r)
        m = obj.toMatrix(r);
        out = m(:,c);
        return
      end

      if r < 1 || r > obj.nrows
        error('ArrayOfArrays: row index must be a scalar in [1, %d].', obj.nrows);
      end

      % get row slice
      row_data = obj.getArray(r);

      if ischar(c) && c == ':'
        out = row_data;
      else
        out = row_data(c);
      end
    end


    function mat = toMatrix(obj, rows)
      % TOMATRIX  Reconstruct a zero-padded dense matrix.
      %
      %   M = obj.toMatrix()       – all rows         (fast path)
      %   M = obj.toMatrix(r1:r2)  – contiguous chunk (fast path)
      %
      %
      %   colId(k) = global_position(k) - ptr(row_of_k) + 1
      %
      % Full matrix:  chunk = obj.data  (no copy, used as-is)
      %               offsets from obj.ptr directly
      %
      % Chunk:        chunk = obj.data(a:b)  (one contiguous slice)
      %               offsets from sptr = ptr(rows) shifted to 1-base

      if nargin == 1
        % entire matrix
        lRows = diff(obj.ptr);
        nRows = obj.nrows;
        nCols = max(lRows);
        mat   = zeros(nRows, nCols);

        rowId = repelem(1:nRows, lRows);
        %  ptr(rowId) is the start of each element's row in data
        %  subtracting it gives the 1-based column within that row
        colId = (1:numel(obj.data)) - obj.ptr(rowId) + 1;

        mat(sub2ind([nRows, nCols], rowId, colId)) = obj.data;

      else
        % contiguous chunk
        rows  = rows(:)';
        nRows = numel(rows);

        if rows(end) - rows(1) ~= nRows - 1
          error('ArrayOfArrays: toMatrix only supports contiguous row ranges.');
        end
        if rows(1) < 1 || rows(end) > obj.nrows
          error('ArrayOfArrays: row range [%d,%d] out of bounds [1,%d].', ...
            rows(1), rows(end), obj.nrows);
        end

        lRows = diff(obj.ptr(rows(1) : rows(end)+1));
        nCols = max(lRows);
        mat   = zeros(nRows, nCols);

        % one contiguous slice — no gather needed
        chunk = obj.data(obj.ptr(rows(1)) : obj.ptr(rows(end)+1)-1);

        % shift ptr so the chunk starts at 1
        sptr  = obj.ptr(rows) - obj.ptr(rows(1)) + 1;

        rowId = repelem(1:nRows, lRows);
        colId = (1:numel(chunk)) - sptr(rowId) + 1;

        mat(sub2ind([nRows, nCols], rowId, colId)) = chunk;
      end
    end



    function sz = size(obj)

      sz = obj.nrows;

    end

    function out = getRows(obj, r)
      %SUBSETROWS Extract arbitrary rows as a new ArrayOfArrays.

      r = r(:);

      if isempty(r)
        out = ArrayOfArrays([], []);
        return
      end

      starts = obj.ptr(r);
      ends   = obj.ptr(r+1) - 1;
      rl     = ends - starts + 1;   % row lengths

      % Fast path: contiguous rows
      if all(diff(r) == 1)
        flat = obj.data(starts(1):ends(end));
      else
        % General case: gather
        nTot = sum(rl);
        flat = zeros(nTot, 1, class(obj.data));

        p = 1;
        for k = 1:numel(r)
          nk = rl(k);
          flat(p:p+nk-1) = obj.data(starts(k):ends(k));
          p = p + nk;
        end
      end
      
      out = ArrayOfArrays(flat, rl);
    end


    function mat = getRowsMatrix(obj, r)
      %GETROWSMATRIX Return rows as a dense matrix.
      % Requires all selected rows to have the same length L.
      % Output size: [numel(r) x L]

      r = r(:)';

      starts = obj.ptr(r);
      ends   = obj.ptr(r+1) - 1;
      lens   = ends - starts + 1;

      if isscalar(r)
        mat = obj.data(starts:ends);
        return
      end

      if any(lens ~= lens(1))
        error('ArrayOfArrays:getRowsMatrix', ...
          'Selected rows must all have the same length.');
      end

      L = lens(1);

      idx = starts + (0:L-1)';   % [L x nRows]
      mat = obj.data(idx)';      % [nRows x L]
    end

    function setArray(obj, r, mat)
      % SETARRAY  Set data of multiple rows with the same length
      %
      %   obj.setArray(r, mat)
      %
      % r   : vector of row indices (not necessarily contiguous)
      % mat : matrix [numel(r) x L] where L matches the row length

      r = r(:)';  % row vector

      if isempty(r)
        return
      end

      % Start and end indices of the rows
      starts = obj.ptr(r);
      ends   = obj.ptr(r+1) - 1;

      % Length of each row
      lens = ends - starts + 1;

      if any(lens ~= lens(1))
        error('All selected rows must have the same length.');
      end

      L = lens(1);

      % Check input matrix dimensions
      if size(mat,1) ~= numel(r)
        error('Number of rows in input matrix must match number of rows selected.');
      end
      if size(mat,2) ~= L
        error('Number of columns in input matrix must match row length.');
      end

      % Compute the linear indices to overwrite
      % Each column of idx is a row's position in obj.data
      idx = starts + (0:L-1)';  % L x numel(r)

      % Write row by row
      obj.data(idx) = mat';
    end



    function out = vertcat(varargin)
      % VERTCAT  Vertical concatenation for ArrayOfArrays
      %
      %   C = [A; B; D; ...]
      %
      % Concatenates rows of multiple ArrayOfArrays objects by simply
      % appending their flat data and their row lengths.

      if nargin == 0
        out = ArrayOfArrays([], []);
        return
      end

      % Validate inputs
      for k = 1:nargin
        if ~isa(varargin{k}, 'ArrayOfArrays')
          error('ArrayOfArrays:vertcat', ...
            'All inputs must be ArrayOfArrays objects.');
        end
      end

      % Total allocation sizes
      totalElems = 0;
      totalRows  = 0;
      for k = 1:nargin
        objk = varargin{k};
        totalElems = totalElems + numel(objk.data);
        totalRows  = totalRows  + objk.nrows;
      end

      % Preallocate
      flat = zeros(1, totalElems, class(varargin{1}.data));
      rl   = zeros(1, totalRows);

      % Fill
      pData = 1;
      pRow  = 1;
      for k = 1:nargin
        objk = varargin{k};

        nk = numel(objk.data);
        nr = objk.nrows;

        flat(pData:pData+nk-1) = objk.data;
        rl(pRow:pRow+nr-1)     = diff(objk.ptr);

        pData = pData + nk;
        pRow  = pRow  + nr;
      end

      out = ArrayOfArrays(flat, rl);
    end

  end


  methods

    function n = arraySize(obj, r)
      % Length of row r.
      if nargin == 1
        n = diff(obj.ptr);
      else
        n = obj.ptr(r+1) - obj.ptr(r);
      end
    end

    function n = totalSize(obj)
      % Length of row r.
      n = numel(obj.data);
    end

    function disp(obj)

      w   = numel(num2str(max(obj.data)));
      fmt = ['%' num2str(w) 'g  '];

      fprintf('Array of Arrays  --  %d rows,  %d total elements\n', ...
        obj.nrows, numel(obj.data));

      for r = 1:obj.nrows
        row = obj.data(obj.ptr(r) : obj.ptr(r+1)-1);
        fprintf('[%s]\n',sprintf(fmt, row));
      end

    end


  end







end

