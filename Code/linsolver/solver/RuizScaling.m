classdef RuizScaling < handle
% RUIZSCALING Manages iterative Ruiz equilibration (symmetric diagonal scaling)
% for block-structured matrices, sparse matrices, and vectors.
%
% Transformations:
%   Matrix:  A_scaled = D * A * D'  (or D_i * A_{ij} * D_j')
%   RHS:     b_scaled = D * b
%   Sol:     x_orig   = D * x_scaled
%   Init:    x0_scale = D^{-1} * x0
%   TestSp:  TV_scale = D^{-1} * TV

   properties (SetAccess = private, GetAccess = public)
      % Diagonal scaling factor vectors for each block (cell array of column vectors)
      D = []

      % Concatenated column vector of all diagonal scaling factors
      fullD = []

      % Number of diagonal blocks
      nBlocks = 0

      % Size of each diagonal block
      blockSizes = []

      % Verbosity flag
      verb = false

      % Maximum iterations for equilibration
      maxit = 10

      % Convergence tolerance
      tol = 1e-3
   end

   properties (Access = public)
      % Flag indicating whether Ruiz scaling is active
      scalingFlag = false
   end

   methods
      % Constructor
      function obj = RuizScaling(debugFlag, maxitin, tolin)
         if nargin >= 1 && ~isempty(debugFlag)
            obj.verb = logical(debugFlag);
         end
         if nargin >= 2 && ~isempty(maxitin)
            obj.maxit = maxitin;
         end
         if nargin >= 3 && ~isempty(tolin)
            obj.tol = tolin;
         end
      end

      % Compute the Ruiz-scaled matrix and store diagonal scaling factors.
      % Supports optional key-value pair 'freezeBlocks', vec.
      function [Abal] = Compute(obj, A, varargin)
         if ~obj.scalingFlag
            Abal = A;
            return;
         end

         if ~iscell(A)
            warning('RuizScaling:Compute:NonCell', ...
               'RuizScaling Compute requires a cell array of block matrices. Scaling disabled.');
            obj.scalingFlag = false;
            Abal = A;
            return;
         end

         % Parse optional key-value pairs
         freezeBlocks = [];
         for k = 1:2:numel(varargin)
            if k+1 <= numel(varargin)
               if strcmpi(varargin{k}, 'freezeBlocks')
                  freezeBlocks = varargin{k+1};
               end
            end
         end

         if isempty(freezeBlocks) || isempty(obj.D)
            % Compute full block equilibration
            [Abal, obj.D] = ruiz_block_symmetric(A, obj.maxit, obj.tol, obj.verb);
            obj.nBlocks = numel(obj.D);
            obj.fullD = vertcat(obj.D{:});
            obj.blockSizes = cellfun(@numel, obj.D);
         else
            % Compute equilibration with specified blocks frozen
            Abal = obj.computeFreezeBlocks(A, freezeBlocks);
         end
      end

      % Internal method to compute Ruiz scaling while keeping specified blocks frozen
      function [Abal] = computeFreezeBlocks(obj, A, frozenBlocks)
         nb = size(A, 1);
         blksize = zeros(nb, 1);

         % Infer block sizes from first non-empty block in each row
         for i = 1:nb
            for j = 1:nb
               if ~isempty(A{i,j})
                  blksize(i) = size(A{i,j}, 1);
                  break;
               end
            end
            if blksize(i) == 0
               error('RuizScaling:computeFreezeBlocks:EmptyRow', ...
                  'Cannot determine size of block-row %d (entire row is empty).', i);
            end
         end

         D_old = obj.D;
         D_init = cell(nb, 1);
         for i = 1:nb
            if ismember(i, frozenBlocks) && ~isempty(D_old) && numel(D_old) >= i && ~isempty(D_old{i})
               D_init{i} = D_old{i};
            else
               D_init{i} = ones(blksize(i), 1);
            end
         end

         % Initialize Abal and scale by initial diagonal factors
         Abal = A;
         for i = 1:nb
            di = D_init{i};
            for j = 1:nb
               if ~isempty(Abal{i,j})
                  Abal{i,j} = di .* Abal{i,j} .* (D_init{j}.');
               end
            end
         end

         % Iteratively equilibrate non-frozen blocks
         activeBlocks = setdiff(1:nb, frozenBlocks);
         D_act = cell(nb, 1);
         for i = 1:nb
            D_act{i} = ones(blksize(i), 1);
         end

         for k = 1:obj.maxit
            s = cell(nb, 1);
            conv = 0;
            for i = activeBlocks
               ri = zeros(blksize(i), 1);
               for j = 1:nb
                  if ~isempty(Abal{i,j})
                     ri = ri + sum(abs(Abal{i,j}).^2, 2);
                  end
               end
               ri = max(sqrt(ri), 1e-15);
               s{i} = 1 ./ sqrt(ri);
               conv = max(conv, max(abs(s{i} - 1)));
            end

            if conv < obj.tol, break; end

            for i = activeBlocks
               D_act{i} = D_act{i} .* s{i};
            end

            for j = activeBlocks
               sj = s{j}.';
               for i = 1:nb
                  if ~isempty(Abal{i,j})
                     Abal{i,j} = Abal{i,j} .* sj;
                  end
               end
            end
            for i = activeBlocks
               si = s{i};
               for j = 1:nb
                  if ~isempty(Abal{i,j})
                     Abal{i,j} = si .* Abal{i,j};
                  end
               end
            end
         end

         % Assemble final scaling vectors
         for i = 1:nb
            obj.D{i} = D_init{i} .* D_act{i};
         end
         obj.nBlocks = numel(obj.D);
         obj.fullD = vertcat(obj.D{:});
         obj.blockSizes = cellfun(@numel, obj.D);
      end

      % Scale a 2D matrix or block-cell matrix.
      %
      % Syntax:
      %   Abal = scaleMat(obj, A)                    -> full cell matrix or full double matrix
      %   Abal = scaleMat(obj, A, iblock)            -> diagonal block i: D_i * A * D_i'
      %   Abal = scaleMat(obj, A, iblock, jblock)    -> off-diagonal block: D_i * A * D_j'
      %   Abal = scaleMat(obj, A, startBlk, endBlk)  -> subrange of cell matrix
      function [Abal] = scaleMat(obj, A, varargin)
         if ~obj.scalingFlag || isempty(obj.D) || isempty(obj.fullD)
            Abal = A;
            return;
         end

         % Guard against vectors to prevent accidental outer-product dense matrix generation
         if ~iscell(A) && isvector(A) && ~isscalar(A)
            error('RuizScaling:scaleMat:VectorInput', ...
               'scaleMat expects a 2D matrix, not a vector. Use applyD for vectors or scaleDiag for diagonal entries.');
         end

         if iscell(A)
            % Determine block range
            if numel(varargin) == 2
               startBlock = varargin{1};
               endBlock = varargin{2};
               range = startBlock:endBlock;
               % If A is an extracted submatrix whose dimensions match numel(range)
               if size(A, 1) == numel(range) && size(A, 2) == numel(range)
                  Abal = A;
                  for i = 1:numel(range)
                     di = obj.D{range(i)};
                     for j = 1:numel(range)
                        if ~isempty(A{i,j})
                           Abal{i,j} = di .* A{i,j} .* (obj.D{range(j)}.');
                        end
                     end
                  end
                  return;
               end
            else
               startBlock = 1;
               endBlock = min(size(A, 1), obj.nBlocks);
            end

            Abal = A;
            for i = startBlock:endBlock
               di = obj.D{i};
               for j = startBlock:endBlock
                  if ~isempty(A{i,j})
                     Abal{i,j} = di .* A{i,j} .* (obj.D{j}.');
                  end
               end
            end

         else
            % Numeric matrix case (sparse or full double)
            if isempty(varargin)
               % Full unpartitioned matrix
               Abal = obj.fullD .* A .* (obj.fullD.');
            elseif isscalar(varargin)
               % Diagonal block i: D_i * A * D_i'
               iblock = varargin{1};
               Abal = obj.D{iblock} .* A .* (obj.D{iblock}.');
            elseif numel(varargin) == 2
               % Off-diagonal block (i, j): D_i * A * D_j'
               iblock = varargin{1};
               jblock = varargin{2};
               Abal = obj.D{iblock} .* A .* (obj.D{jblock}.');
            else
               error('RuizScaling:scaleMat:TooManyArgs', ...
                  'Too many optional arguments passed to scaleMat.');
            end
         end
      end

      % Scale diagonal values from a vector and optionally construct a sparse diagonal matrix.
      %
      % Syntax:
      %   D_mat = scaleDiag(obj, d, iblock)
      %   D_mat = scaleDiag(obj, d, iblock, jblock)
      %   D_mat = scaleDiag(obj, d, iblock, jblock, 'asMatrix', false) -> returns scaled vector
      %
      % Mathematical transformation:
      %   scaled_d = D_i .* d .* D_j
      %   If jblock is omitted, jblock = iblock (D_i^2 .* d).
      function [dOut] = scaleDiag(obj, d, iblock, varargin)
         % Ensure d is a column vector
         dCol = d(:);
         n = numel(dCol);

         % Parse jblock and options
         jblock = iblock;
         asMatrix = true;

         k = 1;
         while k <= numel(varargin)
            arg = varargin{k};
            if isnumeric(arg) && isscalar(arg)
               jblock = arg;
               k = k + 1;
            elseif (ischar(arg) || isstring(arg)) && strcmpi(arg, 'asMatrix')
               if k+1 <= numel(varargin)
                  asMatrix = logical(varargin{k+1});
                  k = k + 2;
               else
                  k = k + 1;
               end
            else
               k = k + 1;
            end
         end

         % Apply scaling if active
         rows = n;
         cols = n;
         if obj.scalingFlag && ~isempty(obj.D)
            di = obj.D{iblock};
            dj = obj.D{jblock};
            rows = numel(di);
            cols = numel(dj);
            minDim = min(rows, cols);
            assert(n == minDim, ...
               'RuizScaling:scaleDiag:DimensionMismatch', ...
               'Diagonal vector length (%d) does not match block dimensions (%d, %d).', ...
               n, rows, cols);
            dCol = di(1:minDim) .* dCol .* dj(1:minDim);
         elseif ~isempty(obj.blockSizes) && numel(obj.blockSizes) >= max(iblock, jblock)
            rows = obj.blockSizes(iblock);
            cols = obj.blockSizes(jblock);
         end

         % Format output
         if asMatrix
            dOut = sparse(1:n, 1:n, dCol, rows, cols);
         else
            if isrow(d)
               dOut = dCol.';
            else
               dOut = dCol;
            end
         end
      end

      % Scale a vector: vec_scaled = D * vec.
      % Preserves row/column orientation and handles cell arrays of vectors.
      %
      % Syntax:
      %   vecBal = applyD(obj, vec)            -> full vector or cell array of vectors
      %   vecBal = applyD(obj, vec, iblock)    -> single block vector
      function [vecBal] = applyD(obj, vec, varargin)
         if ~obj.scalingFlag || isempty(obj.D) || isempty(obj.fullD)
            vecBal = vec;
            return;
         end

         if iscell(vec)
            if numel(varargin) == 2
               startBlock = varargin{1};
               endBlock = varargin{2};
               range = startBlock:endBlock;
               if numel(vec) == numel(range)
                  vecBal = vec;
                  for k = 1:numel(range)
                     vecBal{k} = obj.applyDVector(vec{k}, obj.D{range(k)});
                  end
                  return;
               end
            else
               startBlock = 1;
               endBlock = min(numel(vec), obj.nBlocks);
            end

            vecBal = vec;
            for b = startBlock:endBlock
               vecBal{b} = obj.applyDVector(vec{b}, obj.D{b});
            end

         else
            % Single vector
            if isempty(varargin)
               vecBal = obj.applyDVector(vec, obj.fullD);
            elseif isscalar(varargin)
               iblock = varargin{1};
               vecBal = obj.applyDVector(vec, obj.D{iblock});
            else
               error('RuizScaling:applyD:InvalidArgs', 'Invalid arguments to applyD.');
            end
         end
      end

      % Inverse scale a vector: vec = D^{-1} * vec_scaled.
      % Preserves row/column orientation and handles cell arrays of vectors.
      %
      % Syntax:
      %   vec = applyDinv(obj, vecBal)         -> full vector or cell array of vectors
      %   vec = applyDinv(obj, vecBal, iblock) -> single block vector
      function [vec] = applyDinv(obj, vecBal, varargin)
         if ~obj.scalingFlag || isempty(obj.D) || isempty(obj.fullD)
            vec = vecBal;
            return;
         end

         if iscell(vecBal)
            if numel(varargin) == 2
               startBlock = varargin{1};
               endBlock = varargin{2};
               range = startBlock:endBlock;
               if numel(vecBal) == numel(range)
                  vec = vecBal;
                  for k = 1:numel(range)
                     vec{k} = obj.applyDinvVector(vecBal{k}, obj.D{range(k)});
                  end
                  return;
               end
            else
               startBlock = 1;
               endBlock = min(numel(vecBal), obj.nBlocks);
            end

            vec = vecBal;
            for b = startBlock:endBlock
               vec{b} = obj.applyDinvVector(vecBal{b}, obj.D{b});
            end

         else
            % Single vector
            if isempty(varargin)
               vec = obj.applyDinvVector(vecBal, obj.fullD);
            elseif isscalar(varargin)
               iblock = varargin{1};
               vec = obj.applyDinvVector(vecBal, obj.D{iblock});
            else
               error('RuizScaling:applyDinv:InvalidArgs', 'Invalid arguments to applyDinv.');
            end
         end
      end
   end

   methods (Access = private, Static)
      % Safe vector scaling preserving row/column orientation
      function vOut = applyDVector(v, d)
         if isrow(v)
            vOut = (d.' .* v);
         else
            vOut = (d .* v);
         end
      end

      % Safe vector inverse scaling preserving row/column orientation
      function vOut = applyDinvVector(v, d)
         if isrow(v)
            vOut = (v ./ (d.'));
         else
            vOut = (v ./ d);
         end
      end
   end
end
