classdef (Abstract) preconditioner < handle
   
%   PRECONDITIONER Abstract base class for matrix preconditioners.
%
%   This abstract handle class defines the interface and common properties
%   for preconditioner objects used to accelerate iterative linear solvers.
%   Subclasses must implement the Compute(obj,A,sym,varargin) method to
%   assemble or update the preconditioner for a given matrix A. The class
%   provides readonly function handles Apply_L and Apply_R for applying the
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
   end
   
   methods (Abstract)
      % Function to compute the preconditioner
      Compute(obj,A,sym,varargin)

      % Getter for the function handle to apply the left preconditioner
      x = ApplyLeft(obj,b,varargin)

      % Getter for the function handle to apply the right preconditioner
      x = ApplyRight(obj,b,varargin)
   end
   
   methods
      % Constructor Function
      function obj = preconditioner()

      end

      % Growth checking function default (do nothing), overridden in
      % growing preconditioner
      function x0 = checkGrowth(obj, linsolver, b)

         % Set default output  
         x0 = linsolver.x0;
      end

      % Update the growing preconditioner default (do nothing), overridden in
      % growing preconditioner
      function updateGrowingPrec(obj,Amat)

      end
   end
end


