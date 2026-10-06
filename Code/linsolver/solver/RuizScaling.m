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
      maxit

      % Tolerance to compute D to
      tol

   end

   properties (Access = public)
      % Flag if to use this scaling
      scalingFlag = false;
   end

   methods
      % Initialize the class
      function obj = RuizScaling(debugFlag, maxitin, tolin)
         obj.maxit = maxitin;
         obj.tol = tolin;
         obj.verb = debugFlag;
      end

      % Compute the Ruiz-scaled matrix and store the diagonal scaling
      % matrix
      function [Abal] = Compute(obj, A) 
         if obj.scalingFlag
            if ~iscell(A)
               warning('Non cell A is not supported');
               obj.scalingFlag = false;
               return;
            end
   
            [Abal, obj.D] = ruiz_block_symmetric(A, obj.maxit, obj.tol, obj.verb);
            
            % Get the number of blocks
            obj.nBlocks = numel(obj.D);
  
            % Get the full sparse matrix
            obj.fullD = vertcat(obj.D{:});
         else
            Abal = A;
         end
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
