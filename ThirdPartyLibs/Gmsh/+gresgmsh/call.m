function varargout = call(command, varargin)
%CALL Execute a supported Gmsh operation in the current MATLAB process.
%   Start with call("initialize") and close with call("finalize").
%   Use onCleanup to finalize the session if an operation fails.
%   Supported commands include model.add, option.setNumber,
%   model.occ.{addBox,addPoint,addLine,addCurveLoop,addPlaneSurface,
%   fragment,synchronize}, and model.mesh.{generate,getNodes,
%   getElementType,getElementsByType}. See meshFracturedBox for usage.
assert(ischar(command) && isrow(command) || isstring(command) && isscalar(command), ...
  'GReS:GmshCommand', 'Specify a Gmsh command as a character vector or string scalar.');
assert(gresgmsh.isAvailable(), 'GReS:GmshUnavailable', ...
  'Run setupGmsh and compileAll(''gmsh'') to enable Gmsh.');

[varargout{1:nargout}] = gresGmshGateway(char(command), varargin{:});
end
