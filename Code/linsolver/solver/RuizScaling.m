classdef RuizScaling < handle
   
   properties (SetAccess = private, GetAccess=public)
      % Diagonal scaling matrix
      D = []

      % Sparse double matrix
      fullD = []

      % Number of blocks in D
      nBlocks

      % Verbosity flag
      verb = false
      
      % Max number of iter to compute D
      maxit = 10

      % Tolerance to compute D to
      tol = 1e-3

   end

   properties (Access = public)
      % Flag if to use this scaling
      scalingFlag = false;
   end

   methods
      % Initialize the class
      function obj = RuizScaling(debugFlag, maxitin, tolin)
         if nargin > 1
            obj.verb = debugFlag;
         end
         if nargin > 2
            obj.maxit = maxitin;
         end
         if nargin > 3
            obj.tol = tolin;
         end
      end

      % Compute the Ruiz-scaled matrix and store the diagonal scaling
      % matrix. Supports optional key-value pair 'freezeBlocks', vec
      function [Abal] = Compute(obj, A, varargin) 
         if ~obj.scalingFlag
            Abal = A;
            return;
         end

         if ~iscell(A)
            warning('Non cell A is not supported');
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
            % Compute full Ruiz equilibration
            [Abal, obj.D] = ruiz_block_symmetric(A, obj.maxit, obj.tol, obj.verb);
            obj.nBlocks = numel(obj.D);
            obj.fullD = vertcat(obj.D{:});
         else
            % Compute Ruiz equilibration with specified blocks frozen
            Abal = obj.computeFreezeBlocks(A, freezeBlocks);
         end
      end

      % Internal method to compute Ruiz scaling while keeping specified blocks frozen
      function [Abal] = computeFreezeBlocks(obj, A, frozenBlocks)
         nb = size(A, 1);
         blksize = zeros(nb, 1);

         % Get the size of the matrices
         for i = 1:nb
            blksize(i) = size(A{i,1}, 1);
         end

         D_old = obj.D;
         D_init = cell(nb, 1);
         % Check if frozen values
         for i = 1:nb
            if ismember(i, frozenBlocks)
               D_init{i} = D_old{i};
            else
               D_init{i} = ones(blksize(i), 1);
            end
         end

         % Initialize Abal
         Abal = A;

         % Compute the frozen balanced matrix
         for i = 1:nb
            di = D_init{i};
            for j = 1:nb
               Abal{i,j} = di .* Abal{i,j} .* (D_init{j}.');
            end
         end

         % Get the non frozen diagonal entries
         activeBlocks = setdiff(1:nb, frozenBlocks);
         D_act = cell(nb, 1);
         for i = 1:nb
            D_act{i} = ones(blksize(i), 1);
         end

         % Compute the new diagonal entries
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

         % Finish setting the diagonal scaling with both frozen and new
         % entries
         for i = 1:nb
            if ismember(i, frozenBlocks)
               obj.D{i} = D_old{i};
            else
               obj.D{i} = D_act{i};
            end
         end
         obj.nBlocks = numel(obj.D);
         obj.fullD = vertcat(obj.D{:});
      end

      % Scale a block matrix
      function [Abal] = scaleMat(obj,A,varargin)
         if obj.scalingFlag
            if isempty(obj.D) || isempty(obj.fullD)
               Abal = A;
               return;
            end
            if iscell(A)
               % Validate optional block range
               if ~isempty(varargin)
                  if numel(varargin) == 2
                     startBlock = varargin{1};
                     endBlock = varargin{2};
                  else
                     warning(['scaleMat with cell matrix requires two optional' ...
                        'inputs if any are passed. Falling back to using all blocks']);
                     startBlock = 1;
                     endBlock = obj.nBlocks;
                  end
               else
                  startBlock = 1;
                  endBlock = obj.nBlocks;
               end
      
               % Allocate
               Abal = A;
      
               % Actually scale the matrix
               for i = startBlock:endBlock
                  di = obj.D{i};
                  for j = startBlock:endBlock
                     if ~isempty(A{i,j})
                        Abal{i,j} = di .* A{i,j} .* (obj.D{j}.');
                     end
                  end
               end
            else
               % Single blocks case
               if isequal(obj.nBlocks, 1) || isempty(varargin)
                  Abal = obj.fullD .* A .* (obj.fullD.');
               else
                  if length(varargin) ~= 1
                     warning(['scaleMat with sparse double matrix requires one optional' ...
                        'input if any are passed. Falling back to using the first block']);
                     % Actually scale the matrix
                     Abal = obj.D{1} .* A .* (obj.D{1}.');
                  else
                     iblock = varargin{1};
                     % Actually scale the matrix
                     Abal = obj.D{iblock} .* A .* (obj.D{iblock}.');
                  end
               end
            end
         else
            Abal = A;
         end
      end

      % Scale a vector 
      function [vecBal] = applyD(obj,vec,varargin)
         if obj.scalingFlag
            if isempty(obj.D) || isempty(obj.fullD)
               vecBal = vec;
               return;
            end
            if iscell(vec)
               % Select which blocks are to be scaled or used for the scaling
               if ~isempty(varargin) 
                  if length(varargin) == 2
                     startBlock = varargin{1};
                     endBlock = varargin{2};
                  else
                     warning(['applyD requires two optional inputs if any ' ...
                        ' are passed. Falling back to using all blocks']);
                     startBlock = 1;
                     endBlock = obj.nBlocks;
                  end
               else
                  startBlock = 1;
                  endBlock = obj.nBlocks;
               end
      
               % Actually scale the vector
               range = startBlock:endBlock;
               vecBal = vec;
   
               for k = 1:numel(range)
                  b = range(k);
                  vecBal{b} = obj.D{b} .* vec{b};
               end
            else
               % Single blocks case
               if isequal(obj.nBlocks, 1) || isempty(varargin)
                  vecBal = obj.fullD .* vec;
               else
                  if length(varargin) ~= 1
                     warning(['applyD with sparse double matrix requires one optional' ...
                        ' input if any are passed. Falling back to using the first block']);
                     % Actually scale the vector
                     vecBal = obj.D{1} .* vec;
                  else
                     iblock = varargin{1};
                     % Actually scale the vector
                     vecBal = obj.D{iblock} .* vec;
                  end
               end
            end
         else
            vecBal = vec;
         end
      end

      % Inverse scale a vector 
      function [vec] = applyDinv(obj,vecBal,varargin)
         if obj.scalingFlag
            if isempty(obj.D) || isempty(obj.fullD)
               vec = vecBal;
               return;
            end
            if iscell(vecBal)
               % Select which blocks are to be scaled or used for the scaling
               if ~isempty(varargin) 
                  if length(varargin) == 2
                     startBlock = varargin{1};
                     endBlock = varargin{2};
                  else
                     warning(['applyDinv requires two optional inputs if any ' ...
                        'are passed. Falling back to using all blocks']);
                  end
               else
                  startBlock = 1;
                  endBlock = obj.nBlocks;
               end
      
               % Actually scale the vector
               range = startBlock:endBlock;
               vec = vecBal;
   
               for k = 1:numel(range)
                  b = range(k);
                  vec{b} = obj.D{b} .\ vecBal{b};
               end
            else
               % Single blocks case
               if isequal(obj.nBlocks, 1) || isempty(varargin)
                  vec = obj.fullD .\ vecBal;
               else
                  if length(varargin) ~= 1
                     warning(['applyDinv with sparse double matrix requires one optional' ...
                        'input if any are passed. Falling back to using the first block']);
                     % Actually scale the vector
                     vec = obj.D{1} .\ vecBal;
                  else
                     iblock = varargin{1};
                     % Actually scale the vector
                     vec = obj.D{iblock} .\ vecBal;
                  end
               end
            end
         else
            vec = vecBal;
         end
      end
   end
end
