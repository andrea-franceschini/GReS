function initializeActiveSet(contactSolver,N,varargin)
% initialize an Active Set of size N and write it on input obj
% read xml field ActiveSet

assert(isprop(contactSolver,"activeSet"),"A property named 'activeSet' " + ...
  "is required to initialize the activeSet");

assert(nargin < 4,"Too many input for initializeActiveSet() function")

contactSolver.activeSet.curr = repmat(ContactMode.stick,N,1);
contactSolver.activeSet.prev = contactSolver.activeSet.curr;
% count how many times an element has changed state during iteration
contactSolver.activeSet.stateChange = zeros(N,1);


% extract xml field
if isempty(varargin)  || ismissing(varargin{1})
  activeSetParams = [];
else
  activeSetParams = varargin{1};
end


params = readInput(...
  struct('resetActiveSet',1,...
  'forceStickBoundary',missing,'Tolerances',missing),activeSetParams);

if ismissing(params.forceStickBoundary)
  params.forceStickBoundary = false;
end

contactSolver.activeSet = params;

contactSolver.activeSet.curr = repmat(ContactMode.stick,N,1);
contactSolver.activeSet.prev = contactSolver.activeSet.curr;
% count how many times an element has changed state during iteration
contactSolver.activeSet.stateChange = zeros(N,1);

if ismissing(params.Tolerances)
  tol = [];
else
  tol = params.Tolerances;
end

default = struct('sliding',1e-4,...
  'normalGap',1e-6,...
  'normalTraction',1e-3,...
  'tangentialViolation',1e-4,...
  'areaChange',1e-2,...
  'minLimitTraction',1e-4,...
  'maxStateChange',100);

contactSolver.activeSet.tol = readInput(default,tol);


end

