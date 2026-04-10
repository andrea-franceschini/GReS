classdef Poromechanics < PhysicsSolver

  properties
    K               % the stiffness matrix free of boundary conditions
    fInt            % internal forces
    cell2stress     % map cell ID to position in stress/strain matrix
    avStress
    avStrain
    avStressOld
    avStrainOld
    flOut = true

    % stress and strain tensor use engineering voigt notation
    % s_xx,s_yy,s_zz,tau_yz,tau_xz,tau_xy
  end

  properties (Access = protected)
    iniStress
  end

  properties (Access = private)
    fieldId
    gaussOrder      % (0 means the minimum required by the fem type)
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

      obj.cell2stress = zeros(cells.num,2);

    end

    function assembleSystem(obj,dt)
      % compute the displacements matrices and rhs in the domain
      obj.domain.J{obj.fieldId,obj.fieldId} = computeJacobian(obj,dt);
      obj.domain.rhs{obj.fieldId} = computeRhs(obj);
    end


    function Jmat = computeJacobian(obj,dt)
      if ~isLinear(obj) || isempty(getJacobian(obj))
        % recompute matrix if the model is non linear
        % define size of output matrix
        if isempty(getJacobian(obj))
          % compute strain due to initial boundary displacements
          computeStrain(obj);
        end
        computeStiffMat(obj,dt);
      end

      if obj.domain.simparams.isTimeDependent
        Jmat = obj.domain.simparams.theta*obj.K;
      else
        Jmat = obj.K;
      end
    end


    function computeStiffMat(obj,dt)
      % general sparse assembly loop over elements for Poromechanics

      % define local assembler
      %assembleKloc = @(elemId,counter) computeLocalStiff(obj,elemId,dt,counter);

      % shortcuts
      dofm = obj.domain.dofm;
      coordinates = obj.grid.coordinates;
      cells = obj.grid.cells;
      s = getState(obj);
      sOld = getStateOld(obj);

      % allocate
      subCells = dofm.getFieldCells(obj.fieldId);
      n = sum((obj.grid.nDim^2)*(obj.grid.cells.numVerts(subCells)).^2);
      Ndof = dofm.getNumbDoF(obj.fieldId);
      obj.fInt = zeros(Ndof,1);
      assembleK = assembler(n,Ndof,Ndof);


      for vtkId = cells.vtkTypes

        tmp = obj.grid.getCellsByVTKId(vtkId);
        subCellsLoc = reshape(intersect(subCells,tmp,'sorted'),1,[]);
        elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);

        % get node topology for given vtk type
        topol = obj.grid.getCellNodes(subCellsLoc);

        for i = 1:numel(subCellsLoc)

          % assembly loop for homogeneous element type

          el = subCellsLoc(i);

          % position in stress/strain matrix
          l = obj.cell2stress(el,1);
          nG = obj.cell2stress(el,2);

          nodes = topol(i,:);
          dof = dofm.getLocalDoF(obj.fieldId,nodes);
          coords = coordinates(nodes,:);

          % get strain matrix B
          [gradN,dJWeighed] = getDerBasisFAndDet(elem,coords);
          B = elem.getStrainMatrix(gradN);

          % constitutive update (will be updated with proper
          % constitutiveLaw class)
          [D, sigma, status] = obj.domain.materials.updateMaterial( ...
            cells.tag(el), ...
            sOld.data.stress(l:l+nG-1,:), ...
            s.data.strain(l:l+nG-1,:), ...
            dt, sOld.data.status(l:l+nG-1,:), el, s.t);

          % local stiffnes Kloc = B^T * D * B
          KLoc = obj.computeKloc(B,D,B,dJWeighed);

          % internal forces (initial stresses do not contribute!)
          sz = sigma - obj.iniStress(l:l+nG-1,:);
          sz = reshape(sz',6,1,nG);
          fTmp = pagemtimes(B,'ctranspose',sz,'none');
          fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
          fLoc = sum(fTmp,3);

          % assemble internal forces
          obj.fInt(dof) = obj.fInt(dof)+fLoc;

          % assemble local stiffness
          assembleK.localAssembly(dof,dof,KLoc);

          % update stress map and gp counter
          stateCurr.data.status(l:l+nG-1,:) = status;
          stateCurr.data.stress((l):(l+nG-1),:) = sigma;

        end % end sub cells loop

      end % end vtk loop

      % assemble stiffness matrix
      obj.K = assembleK.sparseAssembly();

    end

   

    % BUBBLES ARE NOT SUPPORTED ANYMORE
    % function [dofr,dofc,Kub,Kbb,varargout] = computeLocalStiffBubble(obj,el,dt,varargin)
    %   % compute local stiffness matrix contribution due to bubble basis
    %   % functions
    %   % the method return Kub and Kbb for later use
    %   % faceId: local index of face holding bubble dof
    %   % assumption: the grid consist only of tetra or hexa, not mixed.
    %   % the unstabilized stiffness has been already assembled
    % 
    %   s = getState(obj);
    %   sOld = getStateOld(obj);
    % 
    %   vtkId = obj.grid.cellVTKType(el);
    %   elem = getElement(obj.elements,vtkId);
    %   nG = elem.GaussPts.nNode;
    %   l = obj.cell2stress(el);      % get index to access stress and strain matrix
    %   [Nu,dJWeighed] = getDerBasisFAndDet(elem,el,1);
    %   Nb = getDerBubbleBasisFAndDet(elem,el,2);
    %   % get strain matrices
    %   Bu = zeros(6,elem.nNode*obj.grid.nDim,nG);
    %   Bu(elem.indB(:,2)) = Nu(elem.indB(:,1));
    %   Bb = zeros(6,elem.nFace*obj.grid.nDim,nG);
    %   Bb(elem.indBbubble(:,2)) = Nb(elem.indBbubble(:,1));
    %   [D, ~, ~] = obj.domain.materials.updateMaterial( ...
    %     obj.grid.cellTag(el), ...
    %     sOld.data.stress(l:l+nG-1,:), ...
    %     s.data.strain(l:l+nG-1,:), ...
    %     dt,sOld.data.status(l:l+nG-1,:), el, s.t);
    % 
    %   Kub = Poromechanics.computeKloc(Bu,D,Bb,dJWeighed);
    %   Kbb = Poromechanics.computeKloc(Bb,D,Bb,dJWeighed);
    % 
    %   % important: right hand side in the unstabilized block already
    %   % considers the bubble contribution due to enhanced strain
    % 
    %   % get global DoF
    %   nodes = obj.grid.cells(el,1:obj.grid.cellNumVerts(el));
    %   dof = obj.domain.dofm.getLocalDoF(obj.fieldId,nodes);
    %   dofr = dof; dofc = dof;
    % 
    %   % get variable output from matrix
    %   if ~isempty(varargin)
    %     varargout = cell(numel(varargin),1);
    %     for i = 1:numel(varargin)
    %       switch varargin{i}
    %         case 'Bb'
    %           varargout{i} = Bb;
    %         case 'Bu'
    %           varargout{i} = Bu;
    %         case 'D'
    %           varargout{i} = D;
    %       end
    %     end
    %   end
    % end

    function initialize(obj)

      % initial stress - assumed balanced with external forces
      state = getState(obj);
      obj.iniStress = state.data.stress;
      
    end



    function initState(obj)
      % add poromechanics fields to state structure
      nGP = FiniteElementType.getTotGPinGrid(obj.grid,obj.gaussOrder);
      state = getState(obj);
      state.data.stress = zeros(nGP,6);
      % initial stress ( assumed balanced with external forces)
      obj.iniStress = zeros(nGP,6);
      state.data.status = zeros(nGP,6);
      state.data.strain = zeros(nGP,6);
      state.data.(obj.getField()) = zeros(obj.grid.nDim*obj.grid.nNodes,1);
    end

    function advanceState(obj)
      % Set converged state to current state after newton convergence
      stateOld = getStateOld(obj);
      state = getState(obj);
      stateOld.data.displacements = getState(obj,"displacements");
      state.data.strain = 0.0*state.data.strain;
      stateOld.data.stress = getState(obj,"stress");
      stateOld.data.status = getState(obj,"status");
      obj.flOut = true;
    end

    function updateState(obj,solution)
      % Update state structure with last solution increment
      dofm = obj.domain.dofm;
      ents = dofm.getActiveEntities(obj.fieldId,1);

      stateCurr = obj.getState();
      %stateOld = obj.getStateOld();

      if nargin > 1
        % apply newton update to current displacements
        stateCurr.data.displacements(ents) = stateCurr.data.displacements(ents) + ...
          solution(getDoF(dofm,obj.fieldId));
      end

      computeStrain(obj);

    end

    function [avStress,avStrain] = finalizeState(obj,stateIn)
      % compute cell average values of stress and strain for print purpose
      % note: output strain is the total one strain = B*u

      cells = obj.grid.cells;

      avStress = zeros(cells.num,6);
      avStrain = zeros(cells.num,6);

      coordinates = obj.grid.coordinates;

      dofm = obj.domain.dofm;

      for vtkId = cells.vtkTypes

        cellList = obj.grid.getCellsByVTKId(vtkId);
        elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);
        nG = elem.getNumbGaussPts;

        % get node topology for given vtk type
        topol = obj.grid.getCellNodes(cellList);

        for i = 1:numel(cellList)

          el = cellList(i);
          l = obj.cell2stress(el,1);
          if l == 0
            % cell is not active
            continue
          end

          stress = stateIn.data.stress(l:l+nG-1,:);

          nodes = topol(i,:);
          coord = coordinates(nodes,:);
          dofs = dofm.getLocalDoF(obj.fieldId,nodes);
          u = stateIn.data.displacements(dofs);

          [gradN,dJWeighed] = getDerBasisFAndDet(elem,coord);

          vol = sum(dJWeighed);

          B = elem.getStrainMatrix(gradN);

          % compute average stress and strain
          avStress(el,:) = sum(diag(dJWeighed)*stress)/vol;
          dStrain = pagemtimes(B,u);
          dStrain = dStrain.*reshape(dJWeighed,1,1,[]);
          avStrain(el,:) = sum(dStrain,3)/vol;

        end
      end

    end

    function applyBC(obj,bcId,t)

      bcType = obj.domain.bcs.getType(bcId);

      switch bcType
        case 'volumeforce' % custom bc type
          applyVolumeForceBC(obj,bcId,bcDofs);
        otherwise
          applyBC@PhysicsSolver(obj,bcId,t);
      end
      
    end


    function applyVolumeForceBC(obj,bcId,t)

      % handle volume force boundary condition.
      % poromechanical coupling with a scalar pressure field

      bc = obj.domain.bcs;

      srcVal = bc.getSourceVals(id,t);

      assert(getField(bc,bcId)==entityField.cell,"field of 'volumeforce' BC" + ...
        " must be 'cell'")

      cells = bc.getSourceEntities(bcId);

      nodeId = bc.getTargetEntities(bcId);

      % preallocate vector for later assembly
      valsC = zeros(3*sum(obj.grid.cellNumVerts(cells)),1);
      dofs = zeros(3*sum(obj.grid.cellNumVerts(cells)),1);
      k = 0;

      for i = 1:numel(cells)

        % local coupling to map cell pressure to nodal force (as in Biot)
        % assumes unit biot coefficient
        el = cells(i);

        % compute local biot coupling matrix
        elem = getElement(obj.elements,obj.grid.cellVTKType(el));
        nG = elem.GaussPts.nNode;
        n = 3*elem.nNode;
        [N,dJWeighed] = getDerBasisFAndDet(elem,el,1);
        B = zeros(6,elem.nNode*obj.grid.nDim,nG);
        B(elem.indB(:,2)) = N(elem.indB(:,1));
        kron = [1;1;1;0;0;0];
        iN = repmat(kron,1,1,nG);
        Qs = pagemtimes(B,'ctranspose',iN,'none');
        Qs = Qs.*reshape(dJWeighed,1,1,[]);
        Qloc = sum(Qs,3);

        % accumulate bc values
        valsC(k+1:k+n) = Qloc*srcVal(i);
        dofState(k+1:k+n) = dofId(obj.grid.cells(el,:),3);
        k = k+n;

      end

      % accumulate results
      valsC = accumarray(dofs,valsC,[3*obj.grid.nNodes 1]);

      dofState = dofId(nodeId,3);

      % retrieve dof numbering of constrained dofs
      dofs = dofId(obj.dofm.getLocalEnts(obj.fieldId,nodeId),3);
      vals = valsC(dofState);

      applyNeuBC(obj,bcId,dofs,vals);

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

    function rhs = computeRhs(obj,varargin)

      dofm = obj.domain.dofm;

      cells = obj.grid.cells;

      if isLinear(obj) % linear case

        J = getJacobian(obj);

        s = obj.getState();
        sOld = obj.getStateOld();

        % update elastic stress
        subCells = dofm.getFieldCells(obj.fieldId);

        for el=subCells'

          D = getElasticTensor(obj.domain.materials.getMaterial(cells.tag(el)).ConstLaw);

          l = obj.cell2stress(el,1);
          nG = obj.cell2stress(el,2);

          % Get the right material stiffness for each element

          s.data.stress((l):(l+nG-1),:) = ...
            sOld.data.stress((l):(l+nG-1),:)+...
            s.data.strain((l):(l+nG-1),:)*D;

        end

        ents = dofm.getActiveEntities(obj.fieldId,1);

        u = s.data.(obj.getField());
        uOld = sOld.data.(obj.getField());

        if obj.domain.simparams.isTimeDependent
          theta = obj.domain.simparams.theta;
          rhs = J*u(ents) + (1/theta-1)*J*uOld(ents);
        else
          rhs = J*u(ents);
        end
      else % non linear case: rhs computed with internal forces (B^T*sigma)
        rhs = obj.fInt; % provisional assuming theta = 1;
      end
      
    end


    function order = getGaussOrder(obj)
      order = obj.gaussOrder;
    end




    function out = isLinear(obj)
      out = false;

      % check if there is not embedded fractures
      if any(contains(obj.domain.solverNames,"EmbeddedFractureMechanics"))
        return
      end


      for i = 1:obj.grid.cells.nTag
        out = isa(obj.domain.materials.getMaterial(i).ConstLaw,"Elastic");
        if ~out
          return;
        end
      end
    end

    function out = isSymmetric(obj)

      % if the problem is linear, then Poromechanics is symmetric
      out = isLinear(obj);

    end

    function writeSolution(obj,fac,tID)

      uOld = getStateOld(obj,obj.getField());
      uCurr = getState(obj,obj.getField());

      obj.domain.outstate.results(tID).(obj.getField()) = uCurr*fac+uOld*(1-fac);

      % TO DO: optional print of stresses as an input for the solver

    end

    function [cellData,pointData] = writeVTK(obj,fac,varargin)

      % append state variable to output structure
      stateCurr = getState(obj);
      stateOld = getStateOld(obj);

      displ = getState(obj,obj.getField());
      dispOld = getStateOld(obj,obj.getField());

      if obj.flOut
        % perform this computation once per time step
        if isempty(obj.avStress)
          [obj.avStress,obj.avStrain] = finalizeState(obj,stateOld);
        end
        obj.avStressOld = obj.avStress;
        obj.avStrainOld = obj.avStrain;
        [obj.avStress,obj.avStrain] = finalizeState(obj,stateCurr);
        obj.flOut = false;
      end

      avgStress = obj.avStress*fac+obj.avStressOld*(1-fac);
      avgStrain = obj.avStrain*fac+obj.avStrainOld*(1-fac);
      displ = displ*fac+dispOld*(1-fac);

      [cellData,pointData] = Poromechanics.buildPrintStruct(displ,avgStress,avgStrain);

    end

  end

  methods (Access=private)

    function computeStrain(obj)

      stateCurr = obj.getState();
      stateOld = obj.getStateOld();
      coordinates = obj.grid.coordinates;

      % displacement increment at current iteration
      du = stateCurr.data.displacements - stateOld.data.displacements;

      % Update strain
      l = 1;


      cells = obj.grid.cells;

      for vtkId = cells.vtkTypes

        cellList = obj.grid.getCellsByVTKId(vtkId);
        elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);
        nG = elem.getNumbGaussPts;

        % get node topology for given vtk type
        topol = obj.grid.getCellNodes(cellList);


        for i = 1:numel(cellList)

          el = cellList(i);

          nodes = topol(i,:);
          coord = coordinates(nodes,:);
          dof = dofId(nodes,3);

          gradN = getDerBasisFAndDet(elem,coord);
          B = elem.getStrainMatrix(gradN);

          stateCurr.data.strain(l:l+nG-1,:) = reshape(pagemtimes(B,du(dof)),6,nG)';

          obj.cell2stress(el,:) = [l nG];
          l = l + nG;

        end

      end

    end



    % function dof = getBCdofs(obj,bcId)
    % 
    %   bc = obj.domain.bcs;
    % 
    %   % get BC entity
    %   ents = bc.getBCentities(bcId);
    % 
    %   % map to local dof numbering
    %   ents = obj.domain.dofm.getLocalEnts(obj.fieldId,ents);
    % 
    %   dof = bc.getCompEntities(bcId,ents);
    % 
    %   if strcmp(bc.getCond(bcId),'VolumeForce')
    %     dof = dofId(dof,3);
    %   end
    % 
    % end
    % 
    % function ents = getBCents(obj,bcId)
    % 
    %   bc = obj.domain.bcs;
    %   % get BC entity
    %   ents = bc.getBCentities(bcId);
    % 
    %   % get component dof for multi-component bcs
    %   ents = bc.getCompEntities(bcId,ents);
    % 
    %   if strcmp(bc.getCond(bcId),'VolumeForce')
    %     ents = dofId(ents,3);
    %   end
    % 
    % 
    % end


    % function vals = getBCVals(obj,id,t)
    % 
    %   bc = obj.domain.bcs;
    % 
    %   vals = bc.getVals(id,t);
    % 
    %   if strcmp(bc.getCond(id),'surface')
    % 
    %     % nodeArea*bcValue
    %     entInfl = bc.getEntitiesInfluence(id);
    %     vals = entInfl*vals;
    % 
    %   elseif strcmp(bc.getCond(id),'volumeforce')
    % 
    %     % imposed volume pressure
    %     valsCell = vals;
    %     cells = bc.getEntities(id);
    % 
    %     % preallocate vector for later assembly
    %     vals = zeros(3*sum(obj.grid.cellNumVerts),1);
    %     dofs = zeros(3*sum(obj.grid.cellNumVerts),1);
    %     k = 0;
    % 
    %     for i = 1:numel(cells)
    % 
    %       % local coupling to map cell pressure to nodal force (as in Biot)
    %       % assumes unit biot coefficient
    %       el = cells(i);
    %       elem = getElement(obj.elements,obj.grid.cellVTKType(el));
    %       nG = elem.GaussPts.nNode;
    %       n = 3*elem.nNode;
    %       [N,dJWeighed] = getDerBasisFAndDet(elem,el,1);
    %       B = zeros(6,elem.nNode*obj.grid.nDim,nG);
    %       B(elem.indB(:,2)) = N(elem.indB(:,1));
    %       kron = [1;1;1;0;0;0];
    %       iN = repmat(kron,1,1,nG);
    %       Qs = pagemtimes(B,'ctranspose',iN,'none'); 
    %       Qs = Qs.*reshape(dJWeighed,1,1,[]);
    %       Qloc = sum(Qs,3);
    %       vals(k+1:k+n) = Qloc*valsCell(i);
    %       dofs(k+1:k+n) = dofId(obj.grid.cells(el,:),3);
    %       k = k+n;
    % 
    %     end
    % 
    %     % accumulate results
    %     vals = accumarray(dofs,vals,[3*obj.grid.nNodes 1]);
    %     dof = obj.getBCdofs(id);
    %     vals = vals(dof);
    %   end
    % 
    % end

  end

  methods (Static)
    
    function [cellStr,pointStr] = buildPrintStruct(disp,stress,strain)
      
      nCellData = 2;
      nPointData = 1;
      pointStr = repmat(struct('name', 1, 'data', 1), nPointData, 1);
      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      % Displacements
      pointStr(1).name = 'displacements';
      pointStr(1).data = [disp(1:3:end) disp(2:3:end) disp(3:3:end)];
      % pointStr(1).name = 'ux';
      % pointStr(1).data = disp(1:3:end);
      % pointStr(2).name = 'uy';
      % pointStr(2).data = disp(2:3:end);
      % pointStr(3).name = 'uz';
      % pointStr(3).data = disp(3:3:end);
      %

      % Permutation needed to be consistent with paraview output

      % Stress
      cellStr(1).name = 'stress';
      cellStr(1).data = stress(:,[1 2 3 6 4 5]);
      cellStr(2).name = 'strain';
      cellStr(2).data = strain(:,[1 2 3 6 4 5]);
      % cellStr(3).name = 'sz';
      % cellStr(3).data = stress(:,3);
      % cellStr(4).name = 'tyz';
      % cellStr(4).data = stress(:,4);
      % cellStr(5).name = 'txz';
      % cellStr(5).data = stress(:,5);
      % cellStr(6).name = 'txy';
      % cellStr(6).data = stress(:,6);
      % %
      % % Strain
      % cellStr(7).name = 'ex';
      % cellStr(7).data = strain(:,1);
      % cellStr(8).name = 'ey';
      % cellStr(8).data = strain(:,2);
      % cellStr(9).name = 'ez';
      % cellStr(9).data = strain(:,3);
      % cellStr(10).name = 'gyz';
      % cellStr(10).data = strain(:,4);
      % cellStr(11).name = 'gxz';
      % cellStr(11).data = strain(:,5);
      % cellStr(12).name = 'gxy';
      % cellStr(12).data = strain(:,6);
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
      indB(:,1) = indB(:,1) + repelem(3*(0:(np-1))',9);
      indB(:,2) = indB(:,2) + repelem(18*(0:(np-1))',9);


    end


    function Kloc = computeKloc(a,b,c,dJW)
      Ks = pagemtimes(pagemtimes(a,'ctranspose',b,'none'),c);
      Ks = Ks.*reshape(dJW,1,1,[]);
      Kloc = sum(Ks,3);
    end

    function out = getField()
      out = "displacements";
    end

  end

end

