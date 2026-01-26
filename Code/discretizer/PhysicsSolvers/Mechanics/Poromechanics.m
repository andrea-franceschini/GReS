classdef Poromechanics < PhysicsSolver

  properties
    K               % the stiffness matrix free of boundary conditions
    fInt            % internal forces
    cell2stress     % map cell ID to position in stress/strain matrix
  end

  properties (Access = private)
    fieldId
  end

  methods (Access = public)

    function obj = Poromechanics(domain)

      % call physicsSolver constructor
      obj@PhysicsSolver(domain);

    end

    function registerSolver(obj,solverInput)

      nTags = obj.mesh.nCellTag;

      if ~isempty(solverInput)
        targetRegions = getXMLData(solverInput,1:nTags,"targetRegions");
      else
        targetRegions = 1:nTags;
      end

      % register nodal displacements on target regions
      obj.dofm.registerVariable(obj.getField(),entityField.node,3,targetRegions);

      % store the id of the field in the degree of freedom manager
      obj.fieldId = obj.dofm.getVariableId(obj.getField());

      % initialize the state object
      initState(obj);

      obj.cell2stress = zeros(obj.mesh.nCells,1);

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
      assembleKloc = @(elemId,counter) computeLocalStiff(obj,elemId,dt,counter);

      subCells = obj.dofm.getFieldCells(obj.fieldId);
      n = sum((obj.mesh.nDim^2)*(obj.mesh.cellNumVerts(subCells)).^2);
      l = 0;
      Ndof = obj.dofm.getNumbDoF(obj.fieldId);
      obj.fInt = zeros(Ndof,1);
      assembleK = assembler(n,Ndof,Ndof,assembleKloc);

      stateCurr = getState(obj);

      % loop over active mechanics cells
      for el = subCells'
        [sigma,status] = assembleK.localAssembly(el,l);
        ng = size(sigma,1);
        stateCurr.data.status(l+1:l+ng,:) = status;
        stateCurr.data.stress((l+1):(l+ng),:) = sigma;
        obj.cell2stress(el) = l;
        l = l + ng;
      end
      % populate stiffness matrix
      obj.K = assembleK.sparseAssembly();
    end

    function [dofr,dofc,KLoc,sigma,status] = computeLocalStiff(obj,elID,dt,l)

      % get the state object
      s = getState(obj);
      sOld = getStateOld(obj);

      vtkId = obj.mesh.cellVTKType(elID);
      elem = getElement(obj.elements,vtkId);
      nG = elem.GaussPts.nNode;
      [N,dJWeighed] = getDerBasisFAndDet(elem,elID,1);
      B = zeros(6,elem.nNode*obj.mesh.nDim,nG);
      B(elem.indB(:,2)) = N(elem.indB(:,1));

      % constitutive update
      [D, sigma, status] = obj.materials.updateMaterial( ...
        obj.mesh.cellTag(elID), ...
        sOld.data.stress(l+1:l+nG,:), ...
        s.data.strain(l+1:l+nG,:), ...
        dt, sOld.data.status(l+1:l+nG,:), elID, s.t);

      % compute local stiffness and internal forces
      KLoc = obj.computeKloc(B,D,B,dJWeighed);
      sz = sigma - s.data.iniStress(l+1:l+nG,:);
      sz = reshape(sz',6,1,nG);
      fTmp = pagemtimes(B,'ctranspose',sz,'none');
      fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
      fLoc = sum(fTmp,3);

      % get global DoF
      nodes = obj.mesh.cells(elID,1:obj.mesh.cellNumVerts(elID));
      dof = obj.dofm.getLocalDoF(obj.fieldId,nodes);
      dofr = dof; dofc = dof;

      % assemble internal forces
      obj.fInt(dof) = obj.fInt(dof)+fLoc;
    end


    function [dofr,dofc,Kub,Kbb,varargout] = computeLocalStiffBubble(obj,el,dt,varargin)
      % compute local stiffness matrix contribution due to bubble basis
      % functions
      % the method return Kub and Kbb for later use
      % faceId: local index of face holding bubble dof
      % assumption: the grid consist only of tetra or hexa, not mixed.
      % the unstabilized stiffness has been already assembled

      s = getState(obj);
      sOld = getStateOld(obj);

      vtkId = obj.mesh.cellVTKType(el);
      elem = getElement(obj.elements,vtkId);
      nG = elem.GaussPts.nNode;
      l = obj.cell2stress(el);      % get index to access stress and strain matrix
      [Nu,dJWeighed] = getDerBasisFAndDet(elem,el,1);
      Nb = getDerBubbleBasisFAndDet(elem,el,2);
      % get strain matrices
      Bu = zeros(6,elem.nNode*obj.mesh.nDim,nG);
      Bu(elem.indB(:,2)) = Nu(elem.indB(:,1));
      Bb = zeros(6,elem.nFace*obj.mesh.nDim,nG);
      Bb(elem.indBbubble(:,2)) = Nb(elem.indBbubble(:,1));
      [D, ~, ~] = obj.materials.updateMaterial( ...
        obj.mesh.cellTag(el), ...
        sOld.data.stress(l+1:l+nG,:), ...
        s.data.strain(l+1:l+nG,:), ...
        dt,sOld.data.status(l+1:l+nG,:), el, s.t);

      Kub = Poromechanics.computeKloc(Bu,D,Bb,dJWeighed);
      Kbb = Poromechanics.computeKloc(Bb,D,Bb,dJWeighed);

      % important: right hand side in the unstabilized block already
      % considers the bubble contribution due to enhanced strain

      % get global DoF
      nodes = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
      dof = obj.dofm.getLocalDoF(obj.fieldId,nodes);
      dofr = dof; dofc = dof;

      % get variable output from matrix
      if ~isempty(varargin)
        varargout = cell(numel(varargin),1);
        for i = 1:numel(varargin)
          switch varargin{i}
            case 'Bb'
              varargout{i} = Bb;
            case 'Bu'
              varargout{i} = Bu;
            case 'D'
              varargout{i} = D;
          end
        end
      end
    end

    function initState(obj)
      % add poromechanics fields to state structure
      Ndata = getNumbCellData(obj.elements);
      state = getState(obj);
      state.data.stress = zeros(Ndata,6);
      state.data.iniStress = zeros(Ndata,6);
      state.data.status = zeros(Ndata,6);
      state.data.strain = zeros(Ndata,6);
      state.data.(obj.getField()) = zeros(obj.mesh.nDim*obj.mesh.nNodes,1);
    end

    function advanceState(obj)
      % Set converged state to current state after newton convergence
      stateOld = getStateOld(obj);
      state = getState(obj);
      stateOld.data.displacements = getState(obj,"displacements");
      state.data.strain = 0.0*state.data.strain;
      stateOld.data.stress = getState(obj,"stress");
      stateOld.data.status = getState(obj,"status");
    end

    function updateState(obj,solution)
      % Update state structure with last solution increment
      ents = obj.dofm.getActiveEntities(obj.fieldId,1);

      stateCurr = obj.getState();
      %stateOld = obj.getStateOld();

      if nargin > 1
        % apply newton update to current displacements
        stateCurr.data.displacements(ents) = stateCurr.data.displacements(ents) + ...
          solution(getDoF(obj.dofm,obj.fieldId));
      end

      computeStrain(obj);
    end

    function [avStress,avStrain] = finalizeState(obj,stateIn)
      % compute cell average values of stress and strain
      avStress = zeros(obj.mesh.nCells,6);
      avStrain = zeros(obj.mesh.nCells,6);

      l = 0;
      for el = 1:obj.mesh.nCells
        dof = getDoFID(obj.mesh,el);
        vtkId = obj.mesh.cellVTKType(el);
        elem = getElement(obj.elements,vtkId);
        nG = elem.GaussPts.nNode;
        vol = obj.mesh.cellVolume(el);
        [N,dJWeighed] = getDerBasisFAndDet(elem,el,1);
        B = zeros(6,elem.nNode*obj.mesh.nDim,nG);
        B(elem.indB(:,2)) = N(elem.indB(:,1));
        avStress(el,:) = sum(diag(dJWeighed)* ...
          stateIn.data.stress((l+1):(l+nG),:))/vol;
        dStrain = pagemtimes(B,stateIn.data.displacements(dof));
        dStrain = dStrain.*reshape(dJWeighed,1,1,[]);
        avStrain(el,:) = sum(dStrain,3)/vol;
        l = l + nG;
      end
    end

    function applyBC(obj,bcId,t)

      if ~BCapplies(obj,bcId)
        return
      end

      % get bcDofs and bcVals
      [bcDofs,bcVals] = getBC(obj,bcId,t);

      bcType = obj.bcs.getType(bcId);

      switch bcType
        case 'Dirichlet'
          applyDirBC(obj,bcId,bcDofs);
        case {'Neumann','VolumeForce'}
          applyNeuBC(obj,bcId,bcDofs,bcVals);
        otherwise
          error("Error in %s: Boundary condition type '%s' is not " + ...
            "available in %s",class(obj),bcType);
      end
      
    end


    function [dof,vals] = getBC(obj,bcId,t)
      %
      dof = obj.getBCdofs(bcId);
      vals = obj.getBCVals(bcId,t);
    end


    function applyDirVal(obj,bcId,t)

      bcVar = obj.bcs.getVariable(bcId);

      if ~strcmp(bcVar,obj.getField())
        return
      end

      [bcDofs,bcVals] = getBC(obj,bcId,t);

      obj.getState().data.displacements(bcDofs) = bcVals;
    end

    function rhs = computeRhs(obj,varargin)

      if isLinear(obj) % linear case

        J = getJacobian(obj);

        s = obj.getState();
        sOld = obj.getStateOld();

        % update elastic stress
        subCells = obj.dofm.getFieldCells(obj.fieldId);
        l1 = 0;

        for el=subCells'

          D = getElasticTensor(obj.materials.getMaterial(obj.mesh.cellTag(el)).ConstLaw);

          % Get the right material stiffness for each element
          vtk = obj.mesh.cellVTKType(el);
          elem = obj.elements.getElement(vtk);
          nG = elem.GaussPts.nNode;
          s.data.stress((l1+1):(l1+nG),:) = ...
            s.data.stress((l1+1):(l1+nG),:)+...
            s.data.strain((l1+1):(l1+nG),:)*D;
          l1 = l1+nG;
        end

        ents = obj.dofm.getActiveEntities(obj.fieldId,1);

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


    function out = isLinear(obj)
      out = false;

      % check if there is not embedded fractures
      if any(contains(obj.domain.solverNames,"EmbeddedFractureMechanics"))
        return
      end


      for i = 1:obj.mesh.nCellTag
        out = isa(obj.materials.getMaterial(i).ConstLaw,"Elastic");
        if ~out
          return;
        end
      end
    end

    function out = isSymmetric(obj)

      % if the problem is linear, then Poromechanics is symmetric
      out = isLinear(obj);

    end

    function writeMatFile(obj,fac,tID)

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
      [avStress,avStrain] = finalizeState(obj,stateCurr);

      dispOld = getStateOld(obj,obj.getField());
      [avStressOld,avStrainOld] = finalizeState(obj,stateOld);

      avStress = avStress*fac+avStressOld*(1-fac);
      avStrain = avStrain*fac+avStrainOld*(1-fac);
      displ = displ*fac+dispOld*(1-fac);

      [cellData,pointData] = Poromechanics.buildPrintStruct(displ,avStress,avStrain);

    end

  end

  methods (Access=private)

    function computeStrain(obj)

      stateCurr = obj.getState();
      stateOld = obj.getStateOld();

      % displacement increment at current iteration
      du = stateCurr.data.displacements - stateOld.data.displacements;

      % Update strain
      l = 0;

      for el=1:obj.mesh.nCells
        dof = getEntityFromElement(entityField.node,entityField.cell,obj.mesh,el,3);
        vtkId = obj.mesh.cellVTKType(el);
        elem = getElement(obj.elements,vtkId);
        nG = elem.GaussPts.nNode;
        N = getDerBasisFAndDet(elem,el,2);
        B = zeros(6,elem.nNode*obj.mesh.nDim,nG);
        B(elem.indB(:,2)) = N(elem.indB(:,1));
        stateCurr.data.strain(l+1:l+nG,:) = reshape(pagemtimes(B,du(dof)),6,nG)';
        l = l + nG;
      end

    end

    function dof = getBCdofs(obj,bcId)

      % get BC entity
      ents = obj.bcs.getBCentities(bcId);

      % get local entity numbering
      ents = obj.dofm.getLocalEnts(obj.fieldId,ents);

      % get component dof for multi-component bcs
      dof = obj.bcs.getCompEntities(bcId,ents);

      % Volume forces are but isotropically act on all directions
      if strcmp(obj.bcs.getCond(bcId),'VolumeForce')
        dof = dofId(dof,3);
      end

    end

    function vals = getBCVals(obj,id,t)

      vals = obj.bcs.getVals(id,t);

      if strcmp(obj.bcs.getCond(id),'SurfBC')

        % nodeArea*bcValue
        entInfl = obj.bcs.getEntitiesInfluence(id);
        vals = entInfl*vals;

      elseif strcmp(obj.bcs.getCond(id),'VolumeForce')

        % imposed volume pressure
        valsCell = vals;
        cells = obj.bcs.getEntities(id);

        % preallocate vector for later assembly
        vals = zeros(3*sum(obj.mesh.cellNumVerts),1);
        dofs = zeros(3*sum(obj.mesh.cellNumVerts),1);
        k = 0;

        for i = 1:numel(cells)

          % local coupling to map cell pressure to nodal force (as in Biot)
          % assumes unit biot coefficient
          el = cells(i);
          elem = getElement(obj.elements,obj.mesh.cellVTKType(el));
          nG = elem.GaussPts.nNode;
          n = 3*elem.nNode;
          [N,dJWeighed] = getDerBasisFAndDet(elem,el,1);
          B = zeros(6,elem.nNode*obj.mesh.nDim,nG);
          B(elem.indB(:,2)) = N(elem.indB(:,1));
          kron = [1;1;1;0;0;0];
          iN = repmat(kron,1,1,nG);
          Qs = pagemtimes(B,'ctranspose',iN,'none'); 
          Qs = Qs.*reshape(dJWeighed,1,1,[]);
          Qloc = sum(Qs,3);
          vals(k+1:k+n) = Qloc*valsCell(i);
          dofs(k+1:k+n) = dofId(obj.mesh.cells(el,:),3);
          k = k+n;

        end

        % accumulate results
        vals = accumarray(dofs,vals,[3*obj.mesh.nNodes 1]);
        dof = obj.getBCdofs(id);
        vals = vals(dof);
      end

    end

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

      % Stress
      cellStr(1).name = 'stress';
      cellStr(1).data = stress;
      cellStr(2).name = 'strain';
      cellStr(2).data = strain;
      % cellStr(3).name = 'sz';
      % cellStr(3).data = stress(:,3);
      % cellStr(4).name = 'txy';
      % cellStr(4).data = stress(:,4);
      % cellStr(5).name = 'tyz';
      % cellStr(5).data = stress(:,5);
      % cellStr(6).name = 'txz';
      % cellStr(6).data = stress(:,6);
      % %
      % % Strain
      % cellStr(7).name = 'ex';
      % cellStr(7).data = strain(:,1);
      % cellStr(8).name = 'ey';
      % cellStr(8).data = strain(:,2);
      % cellStr(9).name = 'ez';
      % cellStr(9).data = strain(:,3);
      % cellStr(10).name = 'gxy';
      % cellStr(10).data = strain(:,4);
      % cellStr(11).name = 'gyz';
      % cellStr(11).data = strain(:,5);
      % cellStr(12).name = 'gxz';
      % cellStr(12).data = strain(:,6);
    end

    function indB = setStrainMatIndex(N)
      % Preapare indices of strain matrix for direct assignment of shape
      % function derivatives
      indB = zeros(9*N,2);
      indB(:,1) = repmat([1, 2, 3, 2, 1, 3, 3, 2, 1],[1,N]);
      indB(:,2) = repmat([1, 4, 6, 8,10,11,15,17,18],[1,N]);
      indB(:,1) = indB(:,1) + repelem(3*(0:(N-1))',9);
      indB(:,2) = indB(:,2) + repelem(18*(0:(N-1))',9);
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

