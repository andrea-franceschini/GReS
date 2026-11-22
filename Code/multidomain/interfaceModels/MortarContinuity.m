classdef MortarContinuity < Constraint

  % Standard mesh tying between non conforming interfaces
  
  properties
    Dmat                      % scalar coupling to slave trace nodes
    Mmat                      % scalar coupling to master trace nodes
    stabilizationMat
    stabilizationScale = 1.0

    % the constraint can be enforced to more than one fields at the same
    % time
    coupledFields

  end

  properties (Access = private)
    multiplierType
  end


  methods

    function obj = MortarContinuity(interface)
      
      obj@Constraint(interface);


    end

    function registerConstraint(obj,varargin)

      params = readInput(struct('stabilizationScale',1.0,...
                                'variable',missing,...
                                'multiplierType',"P0",...
                                'quadrature',missing),varargin{:});

      obj.multiplierType = validatestring(params.multiplierType,{'P0','standard','dual'});
      obj.interface.setMortarQuadrature(obj.multiplierType,params.quadrature);

      multLocation = entityField.node;
      if strcmp(obj.multiplierType,'P0')
        multLocation = entityField.cell;
      end


      obj.coupledFields = getCoupledField(obj,params.variable);
      fld = getField(obj);

      % scaling for jump stabilization of P0 multipliers
      obj.stabilizationScale = params.stabilizationScale;

      dofmD = obj.interface.domains(MortarSide.slave).dofm;
      dofmM = obj.interface.domains(MortarSide.master).dofm;
      dofmI = obj.interface.dofm;

      state = getState(obj);

      for i = 1:numel(obj.coupledFields)
        field = obj.coupledFields(i);
        nComp = dofmD.getNumberOfComponents(field);
        assert(nComp == dofmM.getNumberOfComponents(field),...
          'Mortar field %s must have the same component count on both domains',field);
        dofmI.registerVariable(fld(i),multLocation,nComp,unique(obj.grid.cells.tag));
        state.(fld(i)) = zeros(getNumbDoF(dofmI,fld(i)),1);
      end
      obj.setState(state);

    end


    function updateState(obj,du)

      dofm = obj.interface.dofm;
      for field = reshape(obj.getField(),1,[])
        dofs = dofm.getDoF(dofm.getVariableId(field));
        obj.setState(obj.getState(field) + du(dofs),field);
      end

    end


    function assembleConstraint(obj,varargin)

      interf = obj.interface;
      model = interf.model;
      master = interf.domains(MortarSide.master);
      slave = interf.domains(MortarSide.slave);
      fields = obj.getCoupledField();
      multFields = obj.getField();

      for i = 1:numel(fields)
        field = fields(i);
        multField = multFields(i);
        [D,M] = obj.computeConstraintMatrices(field);

        model.addJacobian(master,M',field,multField,interf);
        model.addJacobian(interf,M,multField,field,master);
        model.addJacobian(slave,D',field,multField,interf);
        model.addJacobian(interf,D,multField,field,slave);

        increment = obj.getState(multField) - obj.getStateInit(multField);
        resMaster = M' * increment;
        resSlave = D' * increment;
        masterState = obj.getVariableField(MortarSide.master,field);
        slaveState = obj.getVariableField(MortarSide.slave,field);
        resMult = M * masterState(:) + D * slaveState(:);

        if strcmp(obj.multiplierType,'P0')
          obj.computeStabilizationMatrix(field,multField);
          model.setJacobian(interf,-obj.stabilizationMat,multField);
          resMult = resMult - obj.stabilizationMat*increment;
        else
          obj.stabilizationMat = [];
        end

        model.addResidualAtDoFs(model.getGlobalDoFs(master,field),resMaster);
        model.addResidualAtDoFs(model.getGlobalDoFs(slave,field),resSlave);
        model.setResidual(interf,resMult,multField);
      end

    end


    function [D,M] = computeConstraintMatrices(obj,field)

      if isempty(obj.Dmat)
        obj.computeScalarMatrices();
      end
      [D,M] = obj.getFieldMatrices(field,"multiplier_" + field);

    end


    function [surfaceStr,pointStr] = writeVTK(obj,fac,varargin)

      surfaceStr = struct('name',{},'data',{});
      pointStr = struct('name',{},'data',{});
      dofm = obj.interface.dofm;
      for field = reshape(obj.getField(),1,[])
        mult = obj.interface.state.interpolate(fac,field);
        nComp = dofm.getNumberOfComponents(field);
        data = struct('name',char(field),'data',reshape(mult,nComp,[])');
        if dofm.getFieldLocation(field) == entityField.cell
          surfaceStr(end+1) = data;
        else
          pointStr(end+1) = data;
        end
      end

    end


    function writeSolution(obj,fac,tID)

      interf = obj.interface;
      for field = reshape(obj.getField(),1,[])
        interf.outstate.results(tID).(field) = interf.state.interpolate(fac,field);
      end

    end


    function [out] = isLinear(obj)
       out = true;
    end


    function fld = getCoupledField(obj,fld)

      if nargin == 1
        if isempty(obj.coupledFields)
          error('Coupled field not yet specified')
        end
        fld = obj.coupledFields;
        return
      end

      dom = obj.interface.domains;

      var1 = dom(1).dofm.getVariableNames;
      var2 = dom(2).dofm.getVariableNames;
      sharedVars = intersect(var1,var2);

      if isscalar(fld) && ismissing(fld)
        fld = sharedVars;
      else
        fld = string(fld);
        if ~all(ismember(fld,sharedVars))
          error('Specified variable is not present in connected domains')
        end
      end

    end




    function fld = getField(obj)

      if isempty(obj.coupledFields)
        error("Coupled fields not yet specified")
      end

      fld = "multiplier_" + obj.coupledFields;

    end


  end


  methods (Access = private)

    function computeScalarMatrices(obj)

      quad = obj.interface.quadrature;
      mesh = quad.grid;
      slaveTrace = mesh.traceGrids(MortarSide.slave);
      masterTrace = mesh.traceGrids(MortarSide.master);
      interfaceElements = obj.createTraceElements(mesh,quad.gaussOrder);
      slaveElements = obj.createTraceElements(slaveTrace,quad.gaussOrder);
      masterElements = obj.createTraceElements(masterTrace,quad.gaussOrder);

      neighbors = mesh.neighbors;
      nodal = ~strcmp(obj.multiplierType,'P0');
      nRows = ones(size(neighbors,1),1);
      nMult = mesh.cells.num;
      if nodal
        nRows = mesh.cells.numVerts(neighbors(:,1));
        nMult = mesh.nNodes;
      end

      nSlave = sum(nRows .* slaveTrace.cells.numVerts(neighbors(:,2)));
      nMaster = sum(nRows .* masterTrace.cells.numVerts(neighbors(:,3)));
      asbD = assembler(nSlave,nMult,slaveTrace.nNodes);
      asbM = assembler(nMaster,nMult,masterTrace.nNodes);

      for triplet = reshape(quad.getActiveTriplets(),1,[])
        ids = neighbors(triplet,:);
        elemI = interfaceElements{mesh.cells.VTKType(ids(1))};
        elemS = slaveElements{slaveTrace.cells.VTKType(ids(2))};
        elemM = masterElements{masterTrace.cells.VTKType(ids(3))};
        [Ns,Nm,Nmult] = quad.getMortarBasisFunctions(...
          triplet,elemI,elemM,elemS,obj.multiplierType);

        rows = ids(1);
        if nodal
          rows = mesh.getCellNodes(ids(1));
          rows = reshape(rows(1:size(Nmult,2)),[],1);
        end
        nodesS = slaveTrace.getCellNodes(ids(2));
        nodesM = masterTrace.getCellNodes(ids(3));
        weights = quad.getIntegrationWeights(triplet);
        weighted = Nmult' .* weights(:)';
        asbD.localAssembly(rows,reshape(nodesS(1:size(Ns,2)),[],1),-weighted*Ns);
        asbM.localAssembly(rows,reshape(nodesM(1:size(Nm,2)),[],1),weighted*Nm);
      end

      obj.Dmat = asbD.sparseAssembly();
      obj.Mmat = asbM.sparseAssembly();

    end


    function [D,M] = getFieldMatrices(obj,field,multiplier)

      interf = obj.interface;
      mesh = interf.grid;
      slave = interf.domains(MortarSide.slave);
      master = interf.domains(MortarSide.master);
      dofm = interf.dofm;
      multId = dofm.getVariableId(multiplier);
      nComp = dofm.getNumberOfComponents(multId);

      D = obj.projectTraceField(slave,mesh.traceGrids(MortarSide.slave),field,obj.Dmat);
      M = obj.projectTraceField(master,mesh.traceGrids(MortarSide.master),field,obj.Mmat);

      % Map rows to the multiplier manager's scalar entity ordering.
      entities = dofm.getActiveEntities(multId);
      dofs = dofm.getLocalDoF(multId,entities);
      [~,order] = sort(dofs(1:nComp:end));
      D = D(entities(order),:);
      M = M(entities(order),:);

      if nComp > 1
        D = kron(D,speye(nComp));
        M = kron(M,speye(nComp));
      end
      if nComp == 3
        R = obj.getMultiplierRotation(mesh,dofm,multId);
        D = R' * D;
        M = R' * M;
      end

    end


    function matrix = projectTraceField(~,domain,trace,field,matrix)

      dofm = domain.dofm;
      id = dofm.getVariableId(field);
      entities = trace.cells.loc2glob;
      nComp = dofm.getNumberOfComponents(id);
      dofs = dofm.getLocalDoF(id,entities);
      firstDofs = dofs(1:nComp:end);
      active = find(firstDofs > 0);
      cols = (firstDofs(active)-1)/nComp + 1;
      transfer = sparse(active,cols,ones(numel(active),1),...
        numel(entities),dofm.getNumbDoF(id)/nComp);
      matrix = matrix * transfer;

    end


    function elements = createTraceElements(~,mesh,gaussOrder)

      % Collect element objects for the input mesh.
      types = double(mesh.cells.vtkTypes);
      elements = cell(1,max(types));
      for vtk = reshape(types,1,[])
        elements{vtk} = FiniteElementType.create(vtk,mesh,'gaussOrder',gaussOrder);
      end

    end


    function R = getMultiplierRotation(~,mesh,dofm,field)

      entities = dofm.getActiveEntities(field);
      normals = mesh.cells.normal;
      if dofm.getFieldLocation(field) == entityField.node
        nodalNormals = zeros(mesh.nNodes,3);
        for cellId = 1:mesh.cells.num
          nodes = mesh.getCellNodes(cellId);
          nodes = nodes(1:mesh.cells.numVerts(cellId));
          nodalNormals(nodes,:) = nodalNormals(nodes,:) + ...
            normals(cellId,:)*mesh.cells.area(cellId);
        end
        normals = nodalNormals;
      end

      nDofs = dofm.getNumbDoF(field);
      blocks = assembler(9*numel(entities),nDofs,nDofs);
      for entity = reshape(entities,1,[])
        dofs = dofm.getLocalDoF(field,entity);
        blocks.localAssembly(dofs,dofs,mxComputeRotationMat(normals(entity,:)));
      end
      R = blocks.sparseAssembly();

    end


    function computeStabilizationMatrix(obj,field,multiplier)

      interf = obj.interface;
      mesh = interf.grid;
      masterTrace = mesh.traceGrids(MortarSide.master);
      master = interf.domains(MortarSide.master);
      slave = interf.domains(MortarSide.slave);
      dofm = interf.dofm;
      multId = dofm.getVariableId(multiplier);
      nComp = dofm.getNumberOfComponents(multId);
      nMult = dofm.getNumbDoF(multId);
      edges = mesh.faces.neighbors(~mesh.faces.isBoundary,:);
      masterEdges = masterTrace.faces.neighbors(~masterTrace.faces.isBoundary,:);
      triplets = mesh.neighbors;

      model = interf.model;
      slaveDiag = diag(cell2mat(model.getJacobian(slave,field)));
      slaveCoupling = cell2mat(model.getJacobian(interf,multiplier,field,slave));
      sameDomain = master == slave;
      if ~sameDomain
        masterDiag = diag(cell2mat(model.getJacobian(master,field)));
        masterCoupling = cell2mat(model.getJacobian(interf,multiplier,field,master));
      end
      stabilization = assembler(size(edges,1)*(2*nComp)^2,nMult,nMult);

      % Master-edge patches are mapped to interface cells through the triplets.
      for edge = 1:size(masterEdges,1)
        inPatch = ismember(triplets(:,3),masterEdges(edge,:));
        cells = unique(triplets(inPatch,1));
        patchEdges = edges(all(ismember(edges,cells),2),:);
        if isempty(patchEdges)
          continue
        end

        patchDofs = dofm.getLocalDoF(multId,cells);
        S = obj.computeSchurPatch(slaveCoupling(patchDofs,:),slaveDiag);
        S = S + obj.computeSchurPatch(masterCoupling(patchDofs,:),masterDiag);

        % D and M already use the local multiplier frames.
        [~,positions] = ismember(patchEdges,cells);
        for i = 1:size(patchEdges,1)
          first = DoFManager.dofExpand(positions(i,1),nComp);
          second = DoFManager.dofExpand(positions(i,2),nComp);
          weight = 0.5*obj.stabilizationScale*(S(first,first) + S(second,second));
          rows = dofm.getLocalDoF(multId,patchEdges(i,:));
          stabilization.localAssembly(rows,rows,[weight,-weight;-weight,weight]);
        end
      end

      obj.stabilizationMat = stabilization.sparseAssembly();

    end


    function S = computeSchurPatch(~,coupling,diagonal)

      cols = find(any(coupling,1));
      coupling = coupling(:,cols);
      S = coupling * (spdiags(1./diagonal(cols),0,numel(cols),numel(cols))*coupling');

    end
  end
end
