function prec = createInnerPrec(type,debugflag,generalsolver,physname)
% CREATEINNERPREC Factory function to instantiate leaf preconditioners ('amg' or 'fsai').
%
% Usage:
%   prec = createInnerPrec(type, debugflag, generalsolver, physname)
%
% Inputs:
%   type          - string or char: 'amg' (default) or 'fsai'
%   debugflag     - boolean: debug output flag
%   generalsolver - general solver object or struct
%   physname      - string with physics name (e.g. 'pressure', 'displacements')
%
% Outputs:
%   prec          - instance of aAMG or aFSAI preconditioner

   if nargin < 1 || isempty(type)
      type = 'amg';
   end
   if nargin < 2 || isempty(debugflag)
      debugflag = false;
   end
   if nargin < 3
      generalsolver = [];
   end
   if nargin < 4
      physname = [];
   end

   switch lower(string(type))
      case "amg"
         prec = aAMG(debugflag, generalsolver, physname);
      case "fsai"
         prec = aFSAI(debugflag, generalsolver, physname);
      otherwise
         error('createInnerPrec:unknownType', ...
               'Unknown preconditioner type: %s. Supported types are "amg" and "fsai".', type);
   end
end
