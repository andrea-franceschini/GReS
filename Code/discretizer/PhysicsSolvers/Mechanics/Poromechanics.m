classdef Poromechanics < PhysicsSolver

  properties
    K               % the stiffness matrix free of boundary conditions
    fInt            % internal forces
    flOut = true

    % stress and strain tensor use engineering voigt notation
    % s_xx,s_yy,s_zz,tau_yz,tau_xz,tau_xy
  end

  properties (Access = protected)
    gaussOrder      % (0 means the minimum required by the fem type)
  end

  properties (Access = private)
    fieldId
  end

  methods (Access = public)

    function obj = Poromechanics(domain)

      % call physicsSolver constructor
      obj@PhysicsSolver(domain);

    end

    function registerSolver(obj,varargin)

      cells = obj.grid.cells;

      default = struct('targetRegions',1:cells.nTag,...
        'gaussOrder',0);


      params = readInput(default,varargin{:});

      dofm = obj.domain.dofm;

      % register nodal displacements on target regions
      dofm.registerVariable(obj.getField(),entityField.node,3,params.targetRegions);

      % gauss integration order for fem
      obj.gaussOrder = params.gaussOrder;

      % store the id of the field in the degree of freedom manager
      obj.fieldId = dofm.getVariableId(obj.getField());

      % initialize the state object
      initState(obj);

    end

    function assembleSystem(obj,dt)


      obj.fInt = zeros(obj.domain.dofm.getNumbDoF(obj.fieldId),1);
      % declare the assembler function
      localAssembler = @(cellList,elem) ...
        assembleMechanics(obj,cellList,elem,dt);

      % utility to assemble the finite element matrix

      obj.K = FEMassembly(obj,localAssembler,'chunkSize',1e4);

      obj.domain.J{obj.fieldId,obj.fieldId} = obj.K;
      obj.domain.rhs{obj.fieldId} = obj.fInt;

    end



    % function assembleSystem(obj,dt)
    %   % compute the displacements matrices and rhs in the domain
    %   assembleSystemMechanics(obj,dt);
    %   obj.domain.J{obj.fieldId,obj.fieldId} = obj.K;
    %   obj.domain.rhs{obj.fieldId} = obj.fInt;
    % end




    function assembleSystemMechanics(obj,dt)
      % general sparse assembly loop over elements for Poromechanics

      % define local assembler
      %assembleKloc = @(elemId,counter) computeLocalStiff(obj,elemId,dt,counter);

      % shortcuts
      dofm = obj.domain.dofm;
      coordinates = obj.grid.coordinates;
      cells = obj.grid.cells;
      s = getState(obj);
      sOld = getStateOld(obj);
      iniStress = getStateInit(obj,'stress');

      t = s.time;

      % allocate
      subCells = dofm.getFieldCells(obj.fieldId);
      n = sum((obj.grid.nDim^2)*(obj.grid.cells.numVerts(subCells)).^2);
      Ndof = dofm.getNumbDoF(obj.fieldId);
      obj.fInt = zeros(Ndof,1);
      assembleK = assembler(n,Ndof,Ndof);

      % get state variables
      du = s.displacements - sOld.displacements;

      gpMap = obj.domain.gpMap;


      if isempty(obj.K) || ~isLinear(obj)
        computeK = true;
      else
        computeK = false;
      end


      for cTag = 1:cells.nTag

        % extract the constitutive law
        constLaw = obj.domain.materials.getConstitutiveLaw(cTag);

        % extract cells belonging to subregion
        subRegionCells = subCells(cells.tag(subCells) == cTag);


        for vtkId = cells.vtkTypes

          % extract cells of the subregion of homogeneous vtk type
          cellId = cells.VTKType(subRegionCells) == vtkId;
          subCellsLoc = subRegionCells(cellId);

          cellList = find(cellId);

          elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);

          % get node topology for given vtk type
          topol = obj.grid.getCellNodes(subCellsLoc);

          nG = elem.getNumbGaussPts;

          for i = 1:numel(subCellsLoc)

            % assembly loop for homogeneous element type

            el = subCellsLoc(i);

            l = gpMap(el,1);

            nodes = topol(i,:);
            dof = dofm.getLocalDoF(obj.fieldId,nodes);
            coords = coordinates(nodes,:);

            % compute strain
            [gradN,dJWeighed] = getDerBasisFAndDet(elem,coords);
            B = elem.getStrainMatrix(gradN);
            s.strain(l:l+nG-1,:) = reshape(pagemtimes(B,du(dof)),6,nG)';


            [sigma,D] = constLaw.constitutiveUpdate(cellList(i),...
              sOld.stress(l:l+nG-1,:),...
              s.strain(l:l+nG-1,:),...
              dt,...
              t);


            % update stress map and gp counter
            s.stress(l:(l+nG-1),:) = sigma;

            % internal forces (initial stresses do not contribute!)
            sz = sigma - iniStress(l:l+nG-1,:);
            sz = reshape(sz',6,1,nG);
            fTmp = pagemtimes(B,'ctranspose',sz,'none');
            fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
            fLoc = sum(fTmp,3);

            % assemble internal forces
            obj.fInt(dof) = obj.fInt(dof)+fLoc;

            if computeK
              % local stiffnes Kloc = B^T * D * B
              KLoc = obj.computeKloc(B,D,B,dJWeighed);
              assembleK.localAssembly(dof,dof,KLoc);
            end


          end % end sub cells loop

        end % end vtk loop

      end % end cell tag region loop

      % assemble stiffness matrix
      if computeK
        obj.K = assembleK.sparseAssembly();
      end


      % update modified state object with new stress and strain
      setState(obj,s);

    end



    function K = assembleMechanics(obj,cellList,elem,dt)

      % assemble local stiffness matrix for a given element type and region
      % tag

      dofm = obj.domain.dofm;
      s = getState(obj);
      sOld = getStateOld(obj);
      iniStress = getStateInit(obj,'stress');
      du = s.displacements - sOld.displacements;
      topol = obj.grid.getCellNodes(cellList);
      nG = elem.getNumbGaussPts;


      cTag = obj.grid.cells.tag(cellList(1));
      constLaw = obj.domain.materials.getConstitutiveLaw(cTag);

      % define assembler for matrix
      nDoF = dofm.getNumbDoF(obj.fieldId);
      nEntries = (obj.grid.nDim^2)*numel(cellList)*elem.nNode.^2;
      assembleK = assembler(nEntries,nDoF,nDoF);

      % vectorized cell processing

      nodes = topol';
      dof = dofm.getLocalDoF(obj.fieldId,nodes(:));
      dof = reshape(dof,3*elem.nNode,numel(cellList));

      coords = obj.grid.coordinates(nodes(:),:);
      coords = reshape(coords,elem.nNode,numel(cellList),3);
      coords = permute(coords,[1,3,2]);

      % batched computation of basis function gradients
      [gradN,dJWeighed] = getDerBasisFAndDet(elem,coords);
      B = elem.getStrainMatrix(gradN);

      gpId = getGaussPointIds(obj.domain.gpMap,cellList);

      % B*du for every Gauss point and cell
      duLoc = reshape(du(dof),size(dof,1),1,1,numel(cellList));
      strainIncrement = pagemtimes(B,duLoc);

      % Preserve state storage:
      % one Gauss point per row, with Gauss points flattened cell by cell.
      s.strain(gpId,:) = reshape(strainIncrement,6,[])';

      % batched constitutive update
      locCells = obj.domain.materials.getMaterialCells(cellList);
      [sigma,D] = constLaw.constitutiveUpdate(locCells,...
        sOld.stress(gpId,:),s.strain(gpId,:),dt,s.time);

      % sigma remains stored as (nG*numel(cellList)) x 6
      s.stress(gpId,:) = sigma;

      stressIncrement = sigma - iniStress(gpId,:);
      stressIncrement = reshape(stressIncrement', ...
        6,1,nG,numel(cellList));

      fLoc = pagemtimes(B,'ctranspose',stressIncrement,'none');
      fLoc = sum(fLoc.*reshape(dJWeighed, ...
        1,1,nG,numel(cellList)),3);
      fLoc = reshape(fLoc,size(dof,1),numel(cellList));

      % Accumulate repeated global DOFs shared by different cells.
      obj.fInt = obj.fInt + accumarray( ...
        dof(:),fLoc(:),[nDoF,1],@sum,0);

      Kloc = obj.computeKloc(B,D,B,dJWeighed);

      % dof is nDofPerCell x nCells and Kloc is
      % nDofPerCell x nDofPerCell x nCells.

      assembleK.localAssembly(dof,dof,Kloc);

      K = assembleK.sparseAssembly;
      
      setState(obj,s);

    end


    function initialize(obj)

      % initial stress - assumed balanced with external forces
      computeAvgStressAndStrain(obj);
      setStateInit(obj,getState(obj));
      setStateOld(obj,getState(obj));

      mat = obj.domain.materials;
      matNames = mat.getMaterialNames;

      for matID = matNames

        tags = mat.getMaterialTags(matID);

        getConstitutiveLaw(obj.domain.materials,matID).initialize(tags,obj.domain);

      end

    end



    function initState(obj)
      % add poromechanics fields to state structure
      nGP = FiniteElementType.getTotGPinGrid(obj.grid,obj.gaussOrder);
      cells = obj.grid.cells;
      state = getState(obj);
      state.stress = zeros(nGP,6);

      state.status = zeros(nGP,6);
      state.strain = zeros(nGP,6);
      state.displacements = zeros(obj.grid.nDim*obj.grid.nNodes,1);
      state.avgStress = zeros(cells.num,6);
      state.avgStrain = zeros(cells.num,6);
      % load current state
      setState(obj,state);

    end


    function goBackState(obj)

      setState(obj,getStateOld(obj));

      for tag = 1:obj.grid.cells.nTag
        getConstitutiveLaw(obj.domain.materials,tag).goBackStatus();
      end

    end

    function advanceState(obj)
      % Set converged state to current state after newton convergence
      stateCurr = getState(obj);
      setStateOld(obj,stateCurr);
      % set strain to zero
      stateCurr.strain(:) = 0.0;

      for tag = 1:obj.grid.cells.nTag
        getConstitutiveLaw(obj.domain.materials,tag).advanceStatus();
      end



      obj.flOut = true;
    end

    function updateState(obj,solution)
      % Update state structure with last solution increment
      dofm = obj.domain.dofm;
      ents = dofm.getActiveEntities(obj.fieldId,1);

      uCurr = getState(obj,obj.getField);
      %stateOld = obj.getStateOld();

      if nargin > 1
        % apply newton update to current displacements
        uCurr(ents) = uCurr(ents) + solution(getDoF(dofm,obj.fieldId));
      end

      setState(obj,uCurr,obj.getField);

    end

    function computeAvgStressAndStrain(obj)

      kernel = @(cellList,elem) computeAvgStressAndStrainKernel(obj,cellList,elem);
      FEMassembly(obj,kernel,'chunkSize',1e5);

    end

 

    function applyBC(obj,bcId,t)

      bcType = obj.domain.bcs.getType(bcId);

      switch bcType
        case 'volumeforce' % custom bc type
          applyVolumeForceBC(obj,bcId,t);
        otherwise
          applyBC@PhysicsSolver(obj,bcId,t);
      end

    end


    function applyVolumeForceBC(obj,bcId,t)

      % handle volume force boundary condition.
      % poromechanical coupling with a scalar pressure field

      bc = obj.domain.bcs;

      srcVal = bc.getSourceVals(bcId,t);

      assert(getField(bc,bcId)==entityField.cell,"field of 'volumeforce' BC" + ...
        " must be 'cell'")

      cellPress = bc.getSourceEntities(bcId);
      cells = obj.grid.cells;
      coordinates = obj.grid.coordinates;

      dofm = obj.domain.dofm;
      nodeId = bc.getTargetEntities(bcId);

      % preallocate vector for later assembly
      vals = zeros(3*sum(cells.numVerts(cellPress)),1);
      dofs = zeros(3*sum(cells.numVerts(cellPress)),1);
      k = 0;

      for vtkId = cells.vtkTypes

        subCellsLoc = obj.grid.getCellsByVTKId(vtkId,cellPress);
        elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);

        % get node topology for given vtk type
        topol = obj.grid.getCellNodes(subCellsLoc);

        nG = elem.getNumbGaussPts;

        bcvalLoc = srcVal(cells.VTKType(cellPress)==vtkId);

        for i = 1:numel(subCellsLoc)

          % assembly loop for homogeneous element type

          %el = subCellsLoc(i);

          nodes = topol(i,:);
          dof = dofm.getLocalDoF(obj.fieldId,nodes);

          coords = coordinates(nodes,:);

          % compute strain
          [gradN,dJWeighed] = getDerBasisFAndDet(elem,coords);
          B = elem.getStrainMatrix(gradN);
          kron = [1;1;1;0;0;0];
          iN = repmat(kron,1,1,nG);
          Qs = pagemtimes(B,'ctranspose',iN,'none');
          Qs = Qs.*reshape(dJWeighed,1,1,[]);
          Qloc = sum(Qs,3);

          % accumulate bc values
          n = size(Qloc,1);
          vals(k+1:k+n) = Qloc*srcVal(i);
          dofs(k+1:k+n) = dof;
          k = k+n;

        end

      end

      % accumulate results
      bcVals = accumarray(dofs,vals,[dofm.getNumbDoF(obj.fieldId) 1]);
      bcDofs = dofm.getLocalDoF(obj.fieldId);

      applyNeuBC(obj,bcId,bcDofs,bcVals);

    end


    % function applyDirVal(obj,bcId,t)
    %
    %   bcVar = obj.domain.bcs.getVariable(bcId);
    %
    %   if ~strcmp(bcVar,obj.getField())
    %     return
    %   end
    %
    %   bcEnts = getBCents(obj,bcId);
    %   bcVals = getBCVals(obj,bcId,t);
    %
    %   obj.getState().data.displacements(bcEnts) = bcVals;
    % end

    % function rhs = computeRhs(obj,varargin)
    %
    %   dofm = obj.domain.dofm;
    %
    %   cells = obj.grid.cells;
    %
    %   rhs = obj.fInt;
    %
    % if isLinear(obj) % linear case
    %
    %   J = getJacobian(obj);
    %
    %   s = obj.getState();
    %   sOld = obj.getStateOld();
    %
    %   % update elastic stress
    %   subCells = dofm.getFieldCells(obj.fieldId);
    %
    %   for el=subCells'
    %
    %     D = getElasticTensor(obj.domain.materials.getMaterial(cells.tag(el)).ConstLaw);
    %
    %     l = obj.cell2stress(el,1);
    %     nG = obj.cell2stress(el,2);
    %
    %     % Get the right material stiffness for each element
    %
    %     s.data.stress((l):(l+nG-1),:) = ...
    %       sOld.data.stress((l):(l+nG-1),:)+...
    %       s.data.strain((l):(l+nG-1),:)*D;
    %
    %   end
    %
    %   ents = dofm.getActiveEntities(obj.fieldId,1);
    %
    %   u = s.data.(obj.getField());
    %   uOld = sOld.data.(obj.getField());
    %
    %   if obj.domain.simparams.isTimeDependent
    %     theta = obj.domain.simparams.theta;
    %     rhs = J*u(ents) + (1/theta-1)*J*uOld(ents);
    %   else
    %     rhs = J*u(ents);
    %   end
    % else % non linear case: rhs computed with internal forces (B^T*sigma)
    %   rhs = obj.fInt; % provisional assuming theta = 1;
    % end

    %end


    function order = getGaussOrder(obj)
      order = obj.gaussOrder;
    end


    function iniStress = getInitialStress(obj)

      iniStress = obj.iniStress;

    end




    function out = isLinear(obj)

      % check if model is pure linear elasticity

      out = false;

      % check if there is not embedded fractures
      if any(contains(obj.domain.solverNames,"EmbeddedFractureMechanics"))
        return
      end

      for i = 1:obj.grid.cells.nTag
        out = obj.domain.materials.getConstitutiveLaw(i).isLinear;
        if ~out
          return;
        end
      end
    end


    function out = isSymmetric(obj)

      % if the problem is linear, then Poromechanics is also symmetric
      out = true;

    end


    function writeSolution(obj,fac,tID)

      uInterp = obj.domain.state.interpolate(fac,obj.getField());

      obj.domain.outstate.results(tID).(obj.getField()) = uInterp;

    end

    function [cellData,pointData] = writeVTK(obj,fac,varargin)

      if obj.flOut
        % finalize output only once
        computeAvgStressAndStrain(obj);
        obj.flOut = false;
      end

      stateInterp = obj.domain.state.interpolate(fac);

      [cellData,pointData] = Poromechanics.buildPrintStruct(stateInterp);

    end

  end

  methods (Access=private)

    % function computeStrain(obj)
    %
    %   stateCurr = obj.getState();
    %   stateOld = obj.getStateOld();
    %   coordinates = obj.grid.coordinates;
    %
    %   % displacement increment at current iteration
    %   du = stateCurr.displacements - stateOld.displacements;
    %
    %   cells = obj.grid.cells;
    %
    %   for vtkId = cells.vtkTypes
    %
    %     cellList = obj.grid.getCellsByVTKId(vtkId);
    %     elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);
    %     nG = elem.getNumbGaussPts;
    %
    %     % get node topology for given vtk type
    %     topol = obj.grid.getCellNodes(cellList);
    %
    %
    %     for i = 1:numel(cellList)
    %
    %       el = cellList(i);
    %
    %       l = obj.cell2stress(el);
    %
    %       nodes = topol(i,:);
    %       coord = coordinates(nodes,:);
    %       dof = dofId(nodes,3);
    %
    %       gradN = getDerBasisFAndDet(elem,coord);
    %       B = elem.getStrainMatrix(gradN);
    %
    %       stateCurr.strain(l:l+nG-1,:) = reshape(pagemtimes(B,du(dof)),6,nG)';
    %
    %     end
    %
    %   end
    %
    %   setState(obj,stateCurr.strain,'strain')
    %
    % end

    function computeAvgStressAndStrainKernel(obj,cellList,elem)
      % compute cell average values of stress and strain for print purpose
      % note: output strain is the total strain = B*u

      stress = getState(obj,'stress');
      strain = getState(obj,'strain');
      avStress = getState(obj,'avgStress');
      avStrain = getState(obj,'avgStrain');

      topol = obj.grid.getCellNodes(cellList);

      nodes = topol';
      coords = obj.grid.coordinates(nodes(:),:);
      coords = reshape(coords,elem.nNode,numel(cellList),3);
      coords = permute(coords,[1,3,2]);

      gpId = getGaussPointIds(obj.domain.gpMap,cellList);

      stresses = stress(gpId,:);
      strains = strain(gpId,:);

      [~,dJWeighed] = getDerBasisFAndDet(elem,coords);

      nG = elem.getNumbGaussPts();
      nC = numel(cellList);

      dJW = reshape(dJWeighed,nG,nC,1);
      vol = sum(dJW,1);

      stresses = reshape(stresses,nG,nC,6);
      stresses = sum(stresses .* dJW,1);
      stresses = reshape(stresses,nC,6);
      avStress(cellList,:) = stresses ./ vol';

      strains = reshape(strains,nG,nC,6);
      strains = sum(strains .* dJW,1);
      strains = reshape(strains,nC,6);
      avStrain(cellList,:) = strains ./ vol';

      setState(obj,avStress,'avgStress');
      setState(obj,avStrain,'avgStrain');

    end


    function map = getCell2StressMap(obj)
      cells = obj.grid.cells;
      vtks = cells.vtkTypes;

      ngCells = zeros(cells.num,1);

      for i = length(vtks)
        el = FiniteElementType.create(vtks(i));
        cId = getCellsByVTKId(obj.grid,vtks(i));
        ngCells(cId) = el.getNumbGaussPts;
      end

      map = [1;1+cumsum(ngCells)];
      map = map(1:end-1);

    end



  end

  methods (Static)

    function [cellStr,pointStr] = buildPrintStruct(state)

      disp = state.displacements;
      stress = state.avgStress;
      strain = state.avgStrain;

      nCellData = 2;
      nPointData = 1;
      pointStr = repmat(struct('name', 1, 'data', 1), nPointData, 1);
      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      % Displacements
      pointStr(1).name = 'displacements';
      pointStr(1).data = [disp(1:3:end) disp(2:3:end) disp(3:3:end)];

      % Permutation needed to be consistent with paraview output
      % Stress
      cellStr(1).name = 'stress';
      cellStr(1).data = stress(:,[1 2 3 6 4 5]);
      cellStr(2).name = 'strain';
      cellStr(2).data = strain(:,[1 2 3 6 4 5]);

    end

    function indB = setStrainMatIndex(np)
      % Preapare indices of strain matrix for direct assignment of shape
      % function derivatives
      % we use standard voigt notation, hence the order is
      % ex,ey,ez,gyz,gxz,gxy

      % np: number of nodes in the element x number of GP
      % indB(:,1) -> index in matrix of basis function derivatives of size
      % 3xnNxnGP

      % indB(:,2) -> index for B of size 6x(3*nN)xnGP for each gauss point
      % where Ni is the basis function of node i

      indB = zeros(9*np,2);
      indB(:,1) = repmat([1, 3, 2, 2, 3, 1, 3, 2, 1],[1,np]);
      indB(:,2) = repmat([1, 5, 6, 8,10,12,15,16,17],[1,np]);
      indB(:,1) = indB(:,1) + repelem(3*(0:(np-1))',9,1);
      indB(:,2) = indB(:,2) + repelem(18*(0:(np-1))',9,1);


    end


    function Kloc = computeKloc(a,b,c,dJW)

      nC = size(b,4);
      ng = size(b,3);

      Ks = pagemtimes(pagemtimes(a,'ctranspose',b,'none'),c);
      Ks = Ks.*reshape(dJW,1,1,ng,nC);
      Kloc = sum(Ks,3);
      Kloc = squeeze(Kloc);
    end

    function out = getField()
      out = "displacements";
    end

  end

end
