classdef MortarContinuity < Constraint

  % Standard mesh tying between non conforming interfaces
  
  properties
    Dmat                      % scalar coupling to slave trace nodes
    Mmat                      % scalar coupling to master trace nodes
    stabilizationMat
    stabilizationScale = 1.0

  end

  properties (Access = protected)
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

      % preprocess data for cross grid quadrature
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

      % Build frames after quadrature has pruned and renumbered the grids.
      mesh = obj.grid;
      obj.precomputeMultiplierRotations(mesh,multLocation);
      for side = [MortarSide.slave,MortarSide.master]
        obj.precomputeMultiplierRotations(mesh.traceGrids(side),multLocation);
      end

    end


    function R = getMultiplierRotation(obj,varargin)
      % (mesh,dofm,field) assembles the global multiplier rotation.
      % (cellId) uses an interface cell, and (cellId,side) uses a trace cell.
      if nargin == 4
        [mesh,dofm,field] = varargin{:};
        entities = dofm.getActiveEntities(field);
        location = dofm.getFieldLocation(field);
        nDofs = dofm.getNumbDoF(field);
        dofs = reshape(dofm.getLocalDoF(field,entities),3,[]);
      else
        cellId = varargin{1};
        mesh = obj.grid;
        if nargin == 3
          mesh = mesh.traceGrids(varargin{2});
        end
        location = obj.multiplierLocation;
        if location == entityField.cell
          entities = cellId;
        else
          entities = mesh.getCellNodes(cellId);
          entities = entities(1:mesh.cells.numVerts(cellId));
        end
        nDofs = 3*numel(entities);
        dofs = reshape(1:nDofs,3,[]);
      end

      if location == entityField.cell
        rotations = mesh.cells.rotationMatrices;
      else
        rotations = mesh.cells.nodalRotationMatrices;
      end
      if nargin < 4 && location == entityField.cell
        R = reshape(rotations(entities,:),3,3);
        return
      end

      blocks = assembler(9*numel(entities),nDofs,nDofs);
      for i = 1:numel(entities)
        rows = dofs(:,i);
        blocks.localAssembly(rows,rows,reshape(rotations(entities(i),:),3,3));
      end
      R = blocks.sparseAssembly();
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
        % assemble mortar matrices for each input field

        field = fields(i);
        multField = multFields(i);

        % get mortar matrices (with right size and local
        % orientation)
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

        % compute and add stabilization matrix
        if obj.multiplierLocation == entityField.cell
          if isempty(obj.stabilizationMat)
            obj.computeStabilizationMatrix(field,multField);
          end
          model.setJacobian(interf,-obj.stabilizationMat,multField);
          resMult = resMult - obj.stabilizationMat*increment;
        end

        model.addResidual(master,resMaster,field);
        model.addResidual(slave,resSlave,field);
        model.setResidual(interf,resMult,multField);
      end

    end


    function varargout = computeConstraintMatrices(obj,field)

      if isempty(obj.Dmat)
        obj.computeScalarMatrices();
      end

      fld = getField(obj);

      if numel(fld) > 1
        fld = fld(contains(fld,field));
      end

      if nargout == 2
       [varargout{1},varargout{2}] = obj.getMortarMatrices(field,fld);
      end

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

    function initialize(obj)


      removeSlaveBCents(obj)


    end

    function removeSlaveBCents(obj)

      nodSlave = obj.grid.traceGrids(MortarSide.slave).cells.loc2glob;

      bc = obj.interface.domains(MortarSide.slave).bcs;
      bcList = bc.getBCList();


      for bcId = bcList
        if getType(bc,bcId) == BCtype.dirichlet
          bcNodes = intersect(bc.getTargetEntities(bcId),nodSlave);

          if obj.multiplierLocation == entityField.node
            bc.removeTargetEntities(bcId,bcNodes);
          end
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


  methods (Access = protected)

    function loc = multiplierLocation(obj)

      switch obj.multiplierType
        case "P0"
          loc = entityField.cell;
        case {"standard","dual"}
          loc = entityField.node;
      end
    end

    function elements = createTraceElements(~,mesh,gaussOrder)

      % Collect element objects for the input mesh.
      types = double(mesh.cells.vtkTypes);
      elements = cell(1,max(types));
      for vtk = reshape(types,1,[])
        elements{vtk} = FiniteElementType.create(vtk,mesh,'gaussOrder',gaussOrder);
      end

    end


    function computeStabilizationMatrix(obj,field,multiplier)

      % get ingredients
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

      % diagonal approximation of leading jacobian blocks
      model = interf.model;
      Ks = diag(cell2mat(model.getJacobian(slave,field)));
      Km = diag(cell2mat(model.getJacobian(master,field)));
      [D,M] = getMortarMatrices(obj,field,multiplier);

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
        % compute the local Schur complement approximation for the patch
        S = obj.computeSchurPatch(D(patchDofs,:),Ks);
        S = S + obj.computeSchurPatch(M(patchDofs,:),Km);

        % Remember: D and M already use the local multiplier frames.
        [~,positions] = ismember(patchEdges,cells);
        for i = 1:size(patchEdges,1)
          % assemble the jump-jump matrix
          f1 = DoFManager.dofExpand(positions(i,1),nComp);
          f2 = DoFManager.dofExpand(positions(i,2),nComp);
          Sloc = 0.5*obj.stabilizationScale*(S(f1,f1) + S(f2,f2));
          rows = dofm.getLocalDoF(multId,patchEdges(i,:));
          stabilization.localAssembly(rows,rows,[Sloc,-Sloc;-Sloc,Sloc]);
        end
      end

      obj.stabilizationMat = stabilization.sparseAssembly();

    end


  end


  methods (Access = private)

    function computeScalarMatrices(obj)

      quad = obj.interface.quadrature;
      mesh = quad.grid;
      slaveTrace = mesh.traceGrids(MortarSide.slave);
      masterTrace = mesh.traceGrids(MortarSide.master);

      % get handle to all possible finite element involved in quadrature
      % this avoids nested loops over aailable vtk indices
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

      % prepare mortar matrices assemblers
      nSlave = sum(nRows .* slaveTrace.cells.numVerts(neighbors(:,2)));
      nMaster = sum(nRows .* masterTrace.cells.numVerts(neighbors(:,3)));
      asbD = assembler(nSlave,nMult,slaveTrace.nNodes);
      asbM = assembler(nMaster,nMult,masterTrace.nNodes);

      % loop over triplets (I,S,M)
      for triplet = reshape(quad.getActiveTriplets(),1,[])
        ids = neighbors(triplet,:);
        % interface element (I)
        elemI = interfaceElements{mesh.cells.VTKType(ids(1))};
        % slave element (S)
        elemS = slaveElements{slaveTrace.cells.VTKType(ids(2))};
        % master element (M)
        elemM = masterElements{masterTrace.cells.VTKType(ids(3))};

        % use the quadrature object to get the basis functions of (I,S,M)
        [Ns,Nm,Nmult] = quad.getMortarBasisFunctions(...
          triplet,elemI,elemM,elemS,obj.multiplierType);

        % select multiplier indices
        rows = ids(1);
        if nodal
          rows = mesh.getCellNodes(ids(1));
          rows = reshape(rows(1:size(Nmult,2)),[],1);
        end

        % select trace nodes indices
        nodesS = slaveTrace.getCellNodes(ids(2));
        nodesM = masterTrace.getCellNodes(ids(3));

        % unified mortar integration
        weights = quad.getIntegrationWeights(triplet);
        weighted = Nmult' .* weights(:)';
        asbD.localAssembly(rows,reshape(nodesS(1:size(Ns,2)),[],1),-weighted*Ns);
        asbM.localAssembly(rows,reshape(nodesM(1:size(Nm,2)),[],1),weighted*Nm);
      end

      obj.Dmat = asbD.sparseAssembly();
      obj.Mmat = asbM.sparseAssembly();

    end


    function [D,M] = getMortarMatrices(obj,field,multiplier)

      if nargin == 1
        field = obj.coupledFields;
        multiplier = obj.getField;
      end

      interf = obj.interface;
      mesh = interf.grid;
      slave = interf.domains(MortarSide.slave);
      master = interf.domains(MortarSide.master);
      dofm = interf.dofm;
      multId = dofm.getVariableId(multiplier);
      nComp = dofm.getNumberOfComponents(multId);

      % Map columns from local trace grid node id to global bulk grid node
      % id
      D = obj.projectTraceField(slave,mesh.traceGrids(MortarSide.slave),field,obj.Dmat);
      M = obj.projectTraceField(master,mesh.traceGrids(MortarSide.master),field,obj.Mmat);

      % Map rows to the multiplier manager's scalar entity ordering.
      entities = dofm.getActiveEntities(multId);
      dofs = dofm.getLocalDoF(multId,entities);
      [~,order] = sort(dofs(1:nComp:end));
      D = D(entities(order),:);
      M = M(entities(order),:);

      % expand the scalar mortar matrices for multi-valued fields
      if nComp > 1
        D = kron(D,speye(nComp));
        M = kron(M,speye(nComp));
      end

      % apply rotation matrix for multipliers computed in the local frame
      if nComp == 3
        R = obj.getMultiplierRotation(mesh,dofm,multId);
        D = R' * D;
        M = R' * M;
      end

    end


    function matrix = projectTraceField(~,domain,trace,field,matrix)

      % map the column of interface matrices from local trace grid
      % numbering to global domain dof numbering

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


    function precomputeMultiplierRotations(~,mesh,location)
      normals = mesh.cells.normal;
      if location == entityField.node
        field = 'nodalRotationMatrices';
        nodalNormals = zeros(mesh.nNodes,3);
        for cellId = 1:mesh.cells.num
          nodes = mesh.getCellNodes(cellId);
          nodes = nodes(1:mesh.cells.numVerts(cellId));
          nodalNormals(nodes,:) = nodalNormals(nodes,:) + ...
            normals(cellId,:)*mesh.cells.area(cellId);
        end
        normals = nodalNormals;
      else
        field = 'rotationMatrices';
      end

      rotations = zeros(size(normals,1),9);
      for i = 1:size(normals,1)
        R = mxComputeRotationMat(normals(i,:));
        rotations(i,:) = R(:).';
      end
      mesh.cells.(field) = rotations;
    end


    function S = computeSchurPatch(~,B,K)

      cols = find(any(B,1));  % id of coupled dofs 

      invK = 1./K(cols);
      B = B(:,cols);
      S = B * (invK.*B');

    end
  end
end
