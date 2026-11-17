classdef (Abstract) preconditioner < handle
   
%   PRECONDITIONER Abstract base class for matrix preconditioners.
%
%   This abstract handle class defines the interface and common properties
%   for preconditioner objects used to accelerate iterative linear solvers.
%   Subclasses must implement the Compute(obj,A,sym,varargin) method to
%   assemble or update the preconditioner for a given matrix A. The class
%   provides read only function handles Apply_L and Apply_R for applying the
%   preconditioner on the left and right, respectively. A DEBUGflag property
%   is available for internal debugging control.
%
%   Typical usage:
%      % Create subclass instance (e.g., myPrecond) and compute
%      p = myPrecond();
%      p.Compute(A,true);
%      % Apply preconditioner via p.Apply_L and/or p.Apply_R


   properties (Access = protected)
      % Flag for debug
      DEBUGflag = false
   end

   properties (GetAccess = public, SetAccess = protected)
      % Preconditioner application
      Apply_L = []
      Apply_R = []

      % Define the positive definiteness of the problem
      posDef = true;

      % Global diagonal block indices owned by this preconditioner (1..nBlocks)
      precBlocks = []

      % Diagonal block metadata slice matching precBlocks (from describeSystem)
      blockInfo = struct([])

      % Cached unscaled test space struct array matching blockInfo
      testSpace = struct([])
   end

   properties (Access = public)
      % Ruiz scaling object
      Ruiz = []
   end

   methods (Abstract)
      % Function to compute the preconditioner
      A = Compute(obj,A,sym,varargin)

      % Getter for the function handle to apply the left preconditioner
      x = ApplyLeft(obj,b,varargin)

      % Getter for the function handle to apply the right preconditioner
      x = ApplyRight(obj,b,varargin)
   end
   
   methods
      % Constructor Function
      function obj = preconditioner()
         % Instantiate default Ruiz scaling object (scalingFlag = false by default)
         obj.Ruiz = RuizScaling();
      end

      % Growth checking function default (do nothing), overridden in
      % growing preconditioner
      function x0 = checkGrowth(obj, linsolver, b)

         % Set default output  
         if ~isempty(linsolver.x0)
            x0 = linsolver.x0;
         else
            x0 = zeros(size(b));
         end
      end

      % Update the growing preconditioner default (do nothing), overridden in
      % growing preconditioner
      function updateGrowingPrec(obj,Amat)

      end

      % Update the state of the blocks default (do nothing), overridden in
      % specific preconditioners
      function updateStateBlocks(obj, A)

      end

      % Condenses the domains and interfaces in a 2x2 matrix, default does
      % nothing. Overridden in specific preconditioners
      function A = condenseDomains(obj,A)

      end

      % Scale the system matrix and RHS before linear solve (default: do nothing)
      function [A, b] = scaleSystem(obj, A, b)

      end

      % Set diagonal block scope from systemInfo
      function setBlockScope(obj, sysInfo, blockIds)
         obj.precBlocks = reshape(blockIds, 1, []);
         obj.blockInfo = sysInfo.blocks(obj.precBlocks);
      end

      % Get local indices matching kind ("variable" or "interface")
      function idx = getKindBlocks(obj, kind)
         idx = find([obj.blockInfo.kind] == kind);
      end

      % Get local indices matching variable names
      function idx = getVarBlocks(obj, varNames)
         idx = find([obj.blockInfo.kind] == "variable" & ismember([obj.blockInfo.varName], varNames));
      end

      % Get local indices matching physics/variable names across variables and interfaces
      function idx = getPhysicsBlocks(obj, varNames)
         idx = find(ismember([obj.blockInfo.varName], varNames));
      end

      % Set sub-preconditioners (no-op default for leaf preconditioners)
      function setSubPrecs(obj, varargin)

      end

      % Determine adaptive block ordering for activeBlockIds
      function blockOrder = getBlockOrdering(obj, blocks, activeBlockIds)
         % Default for leaf preconditioners (aAMG, aFSAI, direct):
         % Retains activeBlockIds in their natural domain sequence.
         blockOrder = activeBlockIds;
      end

      % Retrieve or compute cached unscaled test space struct array
      function TV0 = getBlockTestSpace(obj, generalsolver)
         if isempty(obj.testSpace)
            obj.testSpace = obj.buildTestSpace(generalsolver);
         end
         TV0 = obj.testSpace;
      end

      % Build test space struct array for specified local diagonal blocks
      function TV0 = buildTestSpace(obj, generalsolver, localIdx)
         if nargin < 3 || isempty(localIdx)
            localIdx = 1:numel(obj.blockInfo);
         end
         TV0 = struct('TV', cell(1, numel(localIdx)), ...
                      'varName', cell(1, numel(localIdx)), ...
                      'domainId', cell(1, numel(localIdx)));
         for k = 1:numel(localIdx)
            b = localIdx(k);
            blk = obj.blockInfo(b);
            TV0(k).varName = blk.varName;
            TV0(k).domainId = blk.domainId;
            if blk.kind == "interface"
               TV0(k).TV = [];
            elseif blk.varName == "displacements"
               TV0(k).TV = mk_rbm_3d(generalsolver.domains(blk.domainId(1)).grid.coordinates);
            else
               TV0(k).TV = ones(blk.nDof, 1);
            end
         end
      end
   end

   methods (Static)
      % Extract assembled numeric test space matrix from TV0 struct or matrix
      function matTV = extractTestSpaceMatrix(TV0)
         if isempty(TV0)
            matTV = [];
         elseif isstruct(TV0)
            validTV = {TV0(~cellfun(@isempty, {TV0.TV})).TV};
            if isempty(validTV)
               matTV = [];
            else
               matTV = vertcat(validTV{:});
            end
         else
            matTV = TV0;
         end
      end
   end
end


