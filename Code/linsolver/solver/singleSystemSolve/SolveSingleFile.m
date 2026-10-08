function [x, flag, relres, iter, resvec] = SolveSingleFile(matFile, ruizFlag, preconditioner, rhsFlag, tol, verb, isCoords, physics, solverType, maxit)
% SolveSingleFile - Direct file-based linear solver for dumped .mat files.
%
% This function loads a linear system dumped from SolveLin (.mat file) and solves it.
%
% Syntax:
%   [x, flag, relres, iter, resvec] = SolveSingleFile(matFile, ...
%       ruizFlag, preconditioner, rhsFlag, tol, verb, isCoords, physics, solverType, maxit)
%
% Inputs:
%   matFile        - String or char path to .mat file containing 'A', 'b', and 'coordinates'.
%   ruizFlag       - (Optional) Overrides ruizFlag from file (if empty, uses file value).
%   preconditioner - (Optional) Overrides preconditioner from file (if empty, uses file value).
%   rhsFlag        - (Optional) 0: input b [default], 1: ones, 2: rand, 3: A*ones.
%   tol            - (Optional) Relative tolerance [default: 1e-6].
%   verb           - (Optional) Verbosity flag [default: true].
%   isCoords       - (Optional) True if coordinates in file [default: true].
%   physics        - (Optional) Physics string [default: "displacements"].
%   solverType     - (Optional) 'auto' [default], 'sqmr', 'gmres', or 'direct'.
%   maxit          - (Optional) Maximum iterations [default: 1000].
%
% Outputs:
%   x      - Solution vector.
%   flag   - Convergence flag (0 = converged, 1 = no convergence).
%   relres - Final relative residual.
%   iter   - Total iteration count.
%   resvec - Residual history vector.
%
% See also SolveSingle.

   if nargin < 1 || isempty(matFile)
      error('SolveSingleFile: MAT file path must be provided.');
   end

   % Delegate directly to SolveSingle's file-based overload
   [x, flag, relres, iter, resvec] = SolveSingle(matFile, ...
      varargin_or_default(2, nargin, ruizFlag, []), ...
      varargin_or_default(3, nargin, preconditioner, []), ...
      varargin_or_default(4, nargin, rhsFlag, 0), ...
      varargin_or_default(5, nargin, tol, 1e-6), ...
      varargin_or_default(6, nargin, verb, true), ...
      varargin_or_default(7, nargin, isCoords, []), ...
      varargin_or_default(8, nargin, physics, "displacements"), ...
      varargin_or_default(9, nargin, solverType, 'auto'), ...
      varargin_or_default(10, nargin, maxit, 1000));

end

function val = varargin_or_default(pos, totalNargin, argVal, defaultVal)
   if pos <= totalNargin && ~isempty(argVal)
      val = argVal;
   else
      val = defaultVal;
   end
end
