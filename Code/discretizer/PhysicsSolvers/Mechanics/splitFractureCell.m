function cellModel = splitFractureCell(mesh, boundaryShape, C)
% Assemble a P1 elastic cut cell and condense its interior displacement traces.
% mesh.points: XY or XYZ coordinates; mesh.simplices: triangles or tetrahedra.
% mesh.facets: physical fracture edges or triangles (one row per simplex face).
% mesh.pressure: scalar or one value per physical facet; positive opens.
% Unmarked faces are intact, including any artificial continuation of a tip.
% boundaryShape evaluates the coarse displacement lifting at exterior vertices.
% C uses GReS engineering Voigt order (xx, yy, zz, yz, xz, xy).

points = mesh.points;
simplices = mesh.simplices;
facets = mesh.facets;
dim = size(points,2);
nVertices = dim + 1;
nSimplex = size(simplices,1);
assert(ismember(dim,[2,3]) && size(simplices,2) == nVertices && nSimplex > 0, ...
  'Expected a nonempty triangle or tetrahedron subdivision.');
assert(all(isfinite(points(:))) && all(simplices(:) >= 1) && ...
  all(simplices(:) <= size(points,1)) && all(simplices(:) == fix(simplices(:))), ...
  'Invalid cut-cell vertices or connectivity.');
assert(size(unique(points,'rows'),1) == size(points,1), ...
  'Use shared geometric vertex IDs. The solver duplicates displacement traces.');
if isempty(facets)
  facets = zeros(0,dim);
end
assert(size(facets,2) == dim && size(unique(sort(facets,2),'rows'),1) == size(facets,1), ...
  'List each physical facet once using its geometric vertex IDs.');

% Collect the owner(s) of each geometric face. Each simplex starts with
% independent copies of its vertices. Intact faces identify matching copies.
localFaces = nchoosek(1:nVertices,dim);
faceVertices = zeros(nSimplex*nVertices,dim);
for i = 1:nVertices
  rows = (i-1)*nSimplex + (1:nSimplex);
  faceVertices(rows,:) = sort(simplices(:,localFaces(i,:)),2);
end
[uniqueFaces,~,faceId] = unique(faceVertices,'rows');
owners = repmat((1:nSimplex)',nVertices,1);
ownerCount = accumarray(faceId,1);
assert(all(ownerCount <= 2),'A geometric face has more than two owners.');
firstOwner = accumarray(faceId,owners,[],@min);
lastOwner = accumarray(faceId,owners,[],@max);
[found, fractureId] = ismember(sort(facets,2),uniqueFaces,'rows');
assert(all(found) && all(ownerCount(fractureId) == 2), ...
  'A fracture facet must be shared by two simplices in the host cell.');

copies = reshape(1:nSimplex*nVertices,nVertices,nSimplex).';
intact = ownerCount == 2;
intact(fractureId) = false;
links = zeros(nnz(intact)*dim,2);
row = 0;
for face = find(intact).'
  [~,i] = ismember(uniqueFaces(face,:),simplices(firstOwner(face),:));
  [~,j] = ismember(uniqueFaces(face,:),simplices(lastOwner(face),:));
  links(row+(1:dim),:) = [copies(firstOwner(face),i).',copies(lastOwner(face),j).'];
  row = row + dim;
end
nCopies = numel(copies);
adjacency = sparse(links(:,1),links(:,2),ones(size(links,1),1),nCopies,nCopies);
[order,~,starts] = dmperm(speye(nCopies) + adjacency + adjacency.');
traceId = zeros(nCopies,1);
for i = 1:numel(starts)-1
  traceId(order(starts(i):starts(i+1)-1)) = i;
end
traces = reshape(traceId(copies),size(copies));
nTrace = numel(starts)-1;
traceVertex = accumarray(traces(:),simplices(:),[nTrace,1],@min);

% Anchor only traces incident on an exterior facet. Interior traces that
% happen to share a geometric vertex with an exterior trace remain independent.
fixed = false(nTrace,1);
for face = find(ownerCount == 1).'
  owner = firstOwner(face);
  [~,local] = ismember(uniqueFaces(face,:),simplices(owner,:));
  fixed(traces(owner,local)) = true;
end
boundary = traceDofs(find(fixed),dim);
interior = traceDofs(find(~fixed),dim);
P = kron(boundaryShape(points(traceVertex(fixed),:)),eye(dim));

nDof = dim*nTrace;
nElementDof = dim*nVertices;
rowIds = zeros(nElementDof^2,nSimplex);
colIds = rowIds;
values = rowIds;
for i = 1:nSimplex
  x = points(simplices(i,:),:);
  jac = x(2:end,:) - x(1,:);
  assert(rcond(jac) > 100*eps,'Degenerate cut-cell simplex.');
  measure = abs(det(jac))/factorial(dim);
  grad = jac \ [-ones(dim,1),eye(dim)];
  B = zeros(6,nElementDof);
  B(1,1:dim:end) = grad(1,:);
  B(2,2:dim:end) = grad(2,:);
  B(6,1:dim:end) = grad(2,:);
  B(6,2:dim:end) = grad(1,:);
  if dim == 3
    B(3,3:dim:end) = grad(3,:);
    B(4,2:dim:end) = grad(3,:);
    B(4,3:dim:end) = grad(2,:);
    B(5,1:dim:end) = grad(3,:);
    B(5,3:dim:end) = grad(1,:);
  end
  dofs = traceDofs(traces(i,:),dim);
  [rows,cols] = ndgrid(dofs);
  rowIds(:,i) = rows(:);
  colIds(:,i) = cols(:);
  stiffness = measure*(B.'*C*B);
  values(:,i) = stiffness(:);
end
K = sparse(rowIds(:),colIds(:),values(:),nDof,nDof);
K = (K+K.')/2;

% The two trace IDs at a facet vertex can coincide, for example at a crack tip.
pressure = mesh.pressure(:);
if isscalar(pressure)
  pressure = repmat(pressure,numel(fractureId),1);
end
assert(numel(pressure) == numel(fractureId) && all(isfinite(pressure)), ...
  'Provide a finite pressure for each physical facet.');
plus = zeros(size(facets));
minus = plus;
normal = zeros(size(facets,1),dim);
force = zeros(nDof,1);
for i = 1:numel(fractureId)
  face = fractureId(i);
  ownerMinus = firstOwner(face);
  ownerPlus = lastOwner(face);
  [~,im] = ismember(facets(i,:),simplices(ownerMinus,:));
  [~,ip] = ismember(facets(i,:),simplices(ownerPlus,:));
  minus(i,:) = traces(ownerMinus,im);
  plus(i,:) = traces(ownerPlus,ip);
  x = points(facets(i,:),:);
  edges = x(2:end,:) - x(1,:);
  if dim == 2
    n = [-edges(2),edges(1)];
  else
    n = cross(edges(1,:),edges(2,:))/2;
  end
  area = norm(n);
  direction = mean(points(simplices(ownerPlus,:),:),1) - ...
    mean(points(simplices(ownerMinus,:),:),1);
  normal(i,:) = sign(dot(n,direction))*n/area;
  nodalLoad = repmat(pressure(i)*area/dim*normal(i,:).',dim,1);
  pDofs = traceDofs(plus(i,:),dim);
  mDofs = traceDofs(minus(i,:),dim);
  force(pDofs) = force(pDofs) + nodalLoad;
  force(mDofs) = force(mDofs) - nodalLoad;
end

Kbar = P.'*K(boundary,boundary)*P;
fbar = P.'*force(boundary);
recovery = zeros(0,size(P,2));
particular = zeros(0,1);
if ~isempty(interior)
  [R,flag] = chol(K(interior,interior));
  pivots = abs(diag(R));
  assert(flag == 0 && min(pivots) > sqrt(eps)*max(pivots), ...
    'SplitFractureMechanics:UnanchoredTrace', ...
    'The local stiffness is singular. Check the anchoring of every rock piece.');
  Kib = K(interior,boundary)*P;
  recovery = -(R \ (R.' \ Kib));
  particular = R \ (R.' \ force(interior));
  Kbar = Kbar + Kib.'*recovery;
  fbar = fbar - Kib.'*particular;
end
cellModel = struct('stiffness',(Kbar+Kbar.')/2,'force',fbar, ...
  'K',K,'load',force,'boundary',boundary,'interior',interior,'lifting',P, ...
  'recovery',recovery,'particular',particular,'traceVertex',traceVertex, ...
  'simplexTraces',traces,'plus',plus,'minus',minus,'normal',normal);
end

function dofs = traceDofs(vertices,dim)
dofs = reshape((1:dim).' + dim*(vertices(:).'-1),[],1);
end
