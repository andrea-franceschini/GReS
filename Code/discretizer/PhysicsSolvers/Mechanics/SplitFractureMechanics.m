classdef SplitFractureMechanics < PhysicsSolver
  % Linear elastic proof of concept with element-local split traces.
  % Register with 'cutCells', a cell array indexed by host cell ID. An entry
  % contains points, simplices, facets and pressure (see splitFractureCell).
  % Supply a conforming subdivision that tiles the host cell. Fracture
  % pressure is prescribed per facet. Intact faces tie their traces exactly.
  % Exterior traces follow a P1 interpolation of the coarse lifting on
  % exterior subfacets. Use matching face subdivisions at shared 3D faces.
  % Exit jumps are zero. This prototype models prescribed open faces.

  properties (SetAccess = private)
    cellModels
  end

  properties (Access = private)
    fieldId
    stiffness
    force
  end

  methods
    function obj = SplitFractureMechanics(domain)
      obj@PhysicsSolver(domain);
    end

    function registerSolver(obj,varargin)
      defaults = struct('targetRegions',1:obj.grid.cells.nTag, ...
        'cutCells',{cell(obj.grid.cells.num,1)});
      params = readInput(defaults,varargin{:});
      dim = obj.grid.nDim;
      assert(ismember(dim,[2,3]),'Expected a 2D or 3D solid grid.');
      assert(iscell(params.cutCells) && numel(params.cutCells) == obj.grid.cells.num, ...
        'cutCells must contain one entry per host cell.');

      cells = find(ismember(obj.grid.cells.tag,params.targetRegions));
      assert(~isempty(cells) && all(ismember(find(~cellfun(@isempty,params.cutCells)),cells)), ...
        'Supply cut-cell subdivisions only for active mechanical cells.');
      obj.cellModels = cell(obj.grid.cells.num,1);
      for id = cells(:).'
        nodes = obj.grid.getCellNodes(id);
        coords = obj.grid.coordinates(nodes,1:dim);
        [shape,subdivision] = hostInterpolation(coords,obj.grid.cells.VTKType(id));
        mesh = params.cutCells{id};
        if isempty(mesh)
          mesh = struct('points',coords,'simplices',subdivision, ...
            'facets',zeros(0,dim),'pressure',0);
        end
        assert(size(mesh.points,2) == dim,'Subdivision and host dimensions differ.');
        law = obj.domain.materials.getConstitutiveLaw(obj.grid.cells.tag(id));
        assert(isa(law,'Elastic'),'SplitFractureMechanics requires Elastic materials.');
        obj.cellModels{id} = splitFractureCell(mesh,shape,law.getElasticTensor());
      end

      dofm = obj.domain.dofm;
      field = obj.getField();
      dofm.registerVariable(field,entityField.node,dim,params.targetRegions);
      obj.fieldId = dofm.getVariableId(field);
      nDof = dofm.getNumbDoF(obj.fieldId);
      entries = sum((dim*obj.grid.cells.numVerts(cells)).^2);
      assembly = assembler(entries,nDof,nDof);
      obj.force = zeros(nDof,1);
      for id = cells(:).'
        nodes = obj.grid.getCellNodes(id);
        dofs = dofm.getLocalDoF(obj.fieldId,nodes);
        model = obj.cellModels{id};
        assembly.localAssembly(dofs,dofs,model.stiffness);
        obj.force(dofs) = obj.force(dofs) + model.force;
      end
      obj.stiffness = assembly.sparseAssembly();
      state = obj.getState();
      state.(field) = zeros(dim*obj.grid.nNodes,1);
      obj.setState(state);
    end

    function assembleSystem(obj,~)
      field = obj.getField();
      active = obj.domain.dofm.getActiveEntities(obj.fieldId,true);
      u = obj.getState(field);
      model = obj.domain.model;
      model.setJacobian(obj.domain,obj.stiffness,obj.fieldId);
      model.setResidual(obj.domain,obj.stiffness*u(active) - obj.force,obj.fieldId);
    end

    function updateState(obj,solution)
      field = obj.getField();
      dofm = obj.domain.dofm;
      active = dofm.getActiveEntities(obj.fieldId,true);
      u = obj.getState(field);
      u(active) = u(active) + solution(dofm.getDoF(obj.fieldId));
      obj.setState(u,field);
    end

    function [jump,trace] = recoverCell(obj,id,u)
      % jump(vertex, component, facet), oriented from minus to plus.
      model = obj.cellModels{id};
      assert(~isempty(model),'Cell %d has no split-fracture mechanics model.',id);
      dim = obj.grid.nDim;
      nodes = obj.grid.getCellNodes(id);
      if nargin < 3
        u = obj.getState(obj.getField());
      end
      coarse = u(DoFManager.dofExpand(nodes,dim));
      trace = zeros(size(model.K,1),1);
      trace(model.boundary) = model.lifting*coarse;
      trace(model.interior) = model.particular + model.recovery*coarse;
      trace = reshape(trace,dim,[]).';
      jump = zeros(dim,dim,size(model.plus,1));
      for face = 1:size(model.plus,1)
        jump(:,:,face) = trace(model.plus(face,:),:) - trace(model.minus(face,:),:);
      end
    end

    function [cellData,pointData] = writeVTK(obj,fac,time)
      dim = obj.grid.nDim;
      u = obj.domain.state.interpolate(fac,obj.getField());
      xyz = zeros(obj.grid.nNodes,3);
      xyz(:,1:dim) = reshape(u,dim,[]).';
      isFractured = cellfun(@(m) ~isempty(m) && ~isempty(m.facets),obj.cellModels);
      cellData = struct('name','isCellFractured','data',double(isFractured));
      pointData = struct('name','displacements','data',xyz);

      prefix = [obj.domain.getOutName(),'_'];
      out = obj.domain.outstate;
      block = obj.domain.vtmBlock;
      for id = find(isFractured).'
        model = obj.cellModels{id};
        [jump,trace] = obj.recoverCell(id,u);
        [bulk,bulkData,fractures,fractureData] = splitFractureVTKData(model,jump,trace);
        out.writeVTKfile(block,sprintf('%sSplitRock_%d',prefix,id), ...
          makeOutputGrid(bulk),time,bulkData,[],[],[]);
        out.writeVTKfile(block,sprintf('%sFractures_%d',prefix,id), ...
          makeOutputGrid(fractures),time,fractureData,[],[],[]);
      end
    end

    function writeSolution(obj,fac,tID)
      field = obj.getField();
      obj.domain.outstate.results(tID).(field) = obj.domain.state.interpolate(fac,field);
    end

    function out = isSymmetric(~)
      out = true;
    end

    function out = isLinear(~)
      out = true;
    end
  end

  methods (Static)
    function out = getField()
      out = "displacements";
    end
  end
end

function [shape,subdivision] = hostInterpolation(coords,vtkType)
dim = size(coords,2);
if (dim == 2 && vtkType == 5) || (dim == 3 && vtkType == 10)
  origin = coords(1,:);
  affine = [ones(dim+1,1),coords-origin];
  shape = @(x) [ones(size(x,1),1),x-origin]/affine;
  subdivision = 1:dim+1;
  return
end

assert((dim == 2 && vtkType == 9) || (dim == 3 && vtkType == 12), ...
  'Supported hosts are linear simplices and Cartesian quads/hexahedra.');
low = min(coords,[],1);
width = max(coords,[],1)-low;
corners = (coords-low)./width;
assert(all(width > 0) && all(corners(:) == 0 | corners(:) == 1) && ...
  size(unique(corners,'rows'),1) == 2^dim,'Tensor hosts must be Cartesian.');
shape = @(x) tensorShape((x-low)./width,corners);
subdivision = [1,2,3;1,3,4];
if dim == 3
  subdivision = [1,2,3,7;1,3,4,7;1,4,8,7;1,8,5,7;1,5,6,7;1,6,2,7];
end
end

function N = tensorShape(x,corners)
N = ones(size(x,1),size(corners,1));
for i = 1:size(x,2)
  N = N.*(x(:,i)*corners(:,i).' + (1-x(:,i))*(1-corners(:,i)).');
end
end

function grid = makeOutputGrid(mesh)
% Minimal Grid for OutState; keep distinct point IDs at coincident traces.
grid = Grid();
grid.coordinates = mesh.points;
grid.nNodes = size(mesh.points,1);
nCells = size(mesh.cells,1);
grid.cells.connectivity = mesh.cells;
grid.cells.num = nCells;
grid.cells.VTKType = repmat(mesh.vtkType,nCells,1);
grid.cells.numVerts = repmat(size(mesh.cells,2),nCells,1);
end
