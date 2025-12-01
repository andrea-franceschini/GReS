function varargout = reshapeBasisFunctions(nc,varargin)
assert(numel(varargin)==nargout);
varargout = cell(1,nargout);
for i = 1:numel(varargin)
  varargout{i} = reshapeBasisF(varargin{i},nc);
end
end

function Nout = reshapeBasisF(basis,nc)

% expand and reshape basis function matrix to match number of components
% input: nEnts x nGP
% output: nc x nEnts x nG

[ng,nn,nt] = size(basis);
Nout = zeros(nc,nc*nn,ng,nt);

for i = 1:nt
  % index pattern for matrix reshaping
  ind = repmat(linspace(1,nc^2,nc),1,nn)+(nc^2)*repelem(0:nn-1,1,nc);
  N = zeros(nc*nn,ng); % initialize single page
  N(ind(1:nc*nn),:) = repelem(basis(:,:,i)',nc,1);
  Nout(:,:,:,i) = reshape(N,[nc,nn*nc,ng]); % reshaped 3D matrix
end
end