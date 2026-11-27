function tests = testMortarPoisson
% Manufactured Poisson solution on two nonmatching Cartesian subdomains.
  tests = functiontests(localfunctions);
end


function testLinearSegment(~)
  checkRates(1,"SegmentBasedQuadrature");
end


function testLinearElement(~)
  checkRates(1,"ElementBasedQuadrature");
end


function testQuadraticSegment(~)
  checkRates(2,"SegmentBasedQuadrature");
end


function testQuadraticElement(~)
  checkRates(2,"ElementBasedQuadrature");
end


function checkRates(degree,quadrature)
  errors = zeros(2,2);
  for level = 1:2
    errors(level,:) = solveMortarPoisson(degree,quadrature,2^(level-1));
  end
  rates = log2(errors(1,:)./errors(2,:));
  fprintf('Q%d %s: L2 rate %.3f, H1 rate %.3f\n',degree,quadrature,rates);
  expected = [degree+1,degree];
  assert(all(rates > expected-0.2 & rates < expected+0.5), ...
    'Unexpected Q%d convergence for %s: L2 = %.3f, H1 = %.3f.', ...
    degree,quadrature,rates);
end


function errors = solveMortarPoisson(degree,quadrature,refinement)
  u = @(x,y,z) cos(pi*y).*cos(pi*z).*(2*x-x.^2+sin(pi*x));
  ux = @(x,y,z) cos(pi*y).*cos(pi*z).*(2-2*x+pi*cos(pi*x));
  uy = @(x,y,z) -pi*sin(pi*y).*cos(pi*z).*(2*x-x.^2+sin(pi*x));
  uz = @(x,y,z) -pi*cos(pi*y).*sin(pi*z).*(2*x-x.^2+sin(pi*x));
  f = @(x,y,z) cos(pi*y).*cos(pi*z).*(-2-4*pi^2*x+2*pi^2*x.^2-3*pi^2*sin(pi*x));

  left = makeGrid([3,4,4]*refinement,[0,1],degree);
  right = makeGrid([4,6,6]*refinement,[1,2],degree);
  slave = makeDomain(left,0,u,f,ux,uy,uz);
  master = makeDomain(right,2,u,f,ux,uy,uz);
  interface = Interface([slave;master],'slaveSurface',6,'masterSurface',5);
  order = 6;
  if quadrature == "ElementBasedQuadrature"
    order = 10;
  end
  interface.addConstraint("MortarContinuity",'variable',"u", ...
    'multiplierType',"standard",'quadrature',struct('type',quadrature,'gaussOrder',order));

  params = SimulationParameters(struct('Start',0,'End',1,'DtInit',1,'DtMin',1,'DtMax',1));
  model = Model([slave;master],interface);
  solver = NonLinearImplicit('simulationparameters',params,'model',model);
  solver.simulationLoop();

  [l2s,h1s] = slave.getPhysicsSolver("Poisson").computeError();
  [l2m,h1m] = master.getPhysicsSolver("Poisson").computeError();
  errors = [hypot(l2s,l2m),hypot(h1s,h1m)];
end


function domain = makeDomain(grid,xBoundary,u,f,ux,uy,uz)
  % The manufactured solution has homogeneous flux on the y and z walls.
  xyz = grid.coordinates;
  nodes = find(xyz(:,1) == xBoundary);
  bc = Boundaries(grid);
  bc.addBC('name',"exact",'type',"dirichlet",'field',"node", ...
    'entityListType',"bcList",'entityList',nodes,'variable',"u");
  bc.addBCEvent("exact",'time',0,'value',u(xyz(nodes,1),xyz(nodes,2),xyz(nodes,3)));
  domain = Domain('grid',grid,'boundaries',bc);
  domain.addPhysicsSolver("Poisson",'gaussOrder',4);
  domain.getPhysicsSolver("Poisson").setAnalSolution(u,f,ux,uy,uz);
end


function grid = makeGrid(n,xRange,degree)
  grid = structuredMesh(n(1),n(2),n(3),xRange,[0,1],[0,1]);
  if degree == 2
    % Use the twice-refined lattice for shared vertices, mid-edge nodes,
    % face centers, and cell centers in VTK connectivity order.
    fine = structuredMesh(2*n(1),2*n(2),2*n(3),xRange,[0,1],[0,1]);
    lattice = round(2*n.*(grid.coordinates - [xRange(1),0,0]));
    corners = grid.cells.connectivity;
    offsets = HexahedronQuadratic.coordLoc + 1;
    grid.cells.connectivity = nodeIds(lattice(corners(:,1),:),offsets,n);
    grid.cells.numVerts(:) = 27;
    grid.cells.VTKType(:) = double(VTKType.Hexa27);

    corners = grid.surfaces.connectivity;
    faceOrigin = lattice(corners(:,1),:);
    first = (lattice(corners(:,2),:) - faceOrigin)/2;
    second = (lattice(corners(:,4),:) - faceOrigin)/2;
    offsets = QuadrilateralQuadratic.coordLoc + 1;
    shape = 2*n + 1;
    for k = 1:size(offsets,1)
      xyz = faceOrigin + offsets(k,1)*first + offsets(k,2)*second + 1;
      grid.surfaces.connectivity(:,k) = sub2ind(shape,xyz(:,1),xyz(:,2),xyz(:,3));
    end
    grid.surfaces.numVerts(:) = 9;
    grid.surfaces.VTKType(:) = double(VTKType.Quad9);
    grid.coordinates = fine.coordinates;
  end
  grid.processGeometry();
end


function ids = nodeIds(origin,offsets,n)
  shape = 2*n + 1;
  ids = sub2ind(shape,origin(:,1) + offsets(:,1)' + 1, ...
    origin(:,2) + offsets(:,2)' + 1,origin(:,3) + offsets(:,3)' + 1);
end
