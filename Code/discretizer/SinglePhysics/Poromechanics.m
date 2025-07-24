classdef Poromechanics < SinglePhysics

  properties        
    fInt            % internal forces
    cell2stress     % map cell ID to position in stress/strain matrix
  end

  properties (Constant)
    field = 'Poromechanics'
  end

  methods (Access = public)
    function obj = Poromechanics(symmod,params,dofManager,grid,mat,state)
      obj@SinglePhysics(symmod,params,dofManager,grid,mat,state);
      setPoromechanics(obj)
    end

%     function state = computeMatOld(obj,state,~,dt)
%       if ~isLinear(obj) || isempty(obj.J)
%         % recompute matrix if the model is non linear
%         
%         obj.computeStiffMat(dt);
%       end
%       if obj.simParams.isTimeDependent
%         obj.J = obj.simParams.theta*obj.J;
%       end
%     end

    function computeMat(obj,~,dt)
      if ~isLinear(obj) || isempty(obj.J)
        % recompute matrix if the model is non linear
        assembler = @(elemId,counter) computeLocalStiff(obj,elemId,dt,counter);
        % define size of output matrix
        computeStiffMat(obj,assembler);
      end
      if obj.simParams.isTimeDependent
        obj.J = obj.simParams.theta*obj.J;
      end
    end


    function computeStiffMat(obj,assembleKloc)
      % general sparse assembly loop over elements for Poromechanics
      subCells = obj.dofm.getFieldCells(obj.field);
      n = sum((obj.mesh.nDim^2)*(obj.mesh.cellNumVerts(subCells)).^2);
      l = 0;
      Ndof = obj.dofm.getNumDoF(obj.field);
      obj.fInt = zeros(Ndof,1);
      assembleK = assembler(n,assembleKloc,Ndof,Ndof);
      % loop over cells
      for el = subCells'
        % get dof id and local matrix
        [dsigma,status] = assembleK.localAssembly(el,l);
        s = size(dsigma,1);
        obj.state.data.curr.status(l+1:l+s,:) = status;
        obj.state.data.curr.stress((l+1):(l+s),:) = dsigma;
        obj.cell2stress(el) = l;
        l = l + s;
      end
      % populate stiffness matrix
      obj.J = assembleK.sparseAssembly();
    end

    function [dofr,dofc,KLoc,sigma,status] = computeLocalStiff(obj,elID,dt,l)
      vtkId = obj.mesh.cellVTKType(elID);
      elem = getElement(obj.elements,vtkId);
      nG = elem.GaussPts.nNode;
      [N,dJWeighed] = getDerBasisFAndDet(elem,elID,1);
      B = zeros(6,elem.nNode*obj.mesh.nDim,nG);
      B(elem.indB(:,2)) = N(elem.indB(:,1));
      [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(elID), ...
        obj.state.data.conv.stress(l+1:l+nG,:), ...
        obj.state.data.curr.strain(l+1:l+nG,:), ...
        dt,obj.state.data.conv.status(l+1:l+nG,:), elID, obj.state.t);
      %           obj.state.data.curr.status(l+1:l+obj.GaussPts.nNode,:) = status;
      %           obj.state.data.curr.stress(l+1:l+obj.GaussPts.nNode,:) = ones(8,6);
      KLoc = obj.computeKloc(B,D,B,dJWeighed);
      sz = sigma - obj.state.data.iniStress(l+1:l+nG,:);
      sz = reshape(sz',6,1,nG);
      fTmp = pagemtimes(B,'ctranspose',sz,'none');
      fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
      fLoc = sum(fTmp,3);
      % get global DoF
      nodes = obj.mesh.cells(elID,1:obj.mesh.cellNumVerts(elID));
      dof = obj.dofm.getLocalDoF(nodes,obj.fldId);
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
      [D, ~, ~] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
        obj.state.data.conv.stress(l+1:l+nG,:), ...
        obj.state.data.curr.strain(l+1:l+nG,:), ...
        dt,obj.state.data.conv.status(l+1:l+nG,:), el, obj.state.t);
      Kub = Poromechanics.computeKloc(Bu,D,Bb,dJWeighed);
      Kbb = Poromechanics.computeKloc(Bb,D,Bb,dJWeighed);

      % important: right hand side in the unstabilized block already
      % considers the bubble contribution due to enhanced strain
      % get global DoF
      nodes = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
      dof = obj.dofm.getLocalDoF(nodes,obj.fldId);
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
      obj.state.data.conv = struct('strain', [], 'stress', [], 'status', []);
      obj.state.data.curr = struct('strain', [], 'stress', [], 'status', []);
      obj.state.data.curr.stress = zeros(Ndata,6);
      obj.state.data.curr.strain = zeros(Ndata,6);
      obj.state.data.curr.status = zeros(Ndata,2);
      obj.state.data.conv = obj.state.data.curr;
      obj.state.data.iniStress = zeros(Ndata,6);
      obj.state.data.dispConv = zeros(obj.mesh.nDim*obj.mesh.nNodes,1);
      obj.state.data.dispCurr = zeros(obj.mesh.nDim*obj.mesh.nNodes,1);
    end

    function advanceState(obj)
      % Set converged state to current state after newton convergence
      obj.state.data.dispConv = obj.state.data.dispCurr;
      % we store the incremental strain during the time step
      obj.state.data.curr.strain = 0.0*obj.state.data.curr.strain;
      obj.state.data.conv.stress = obj.state.data.curr.stress;
      obj.state.data.conv.status = obj.state.data.curr.status;
    end

    function updateState(obj,dSol)
      % Update state structure with last solution increment
      ents = obj.dofm.getActiveEnts(obj.field);
      if nargin > 1
        % apply newton update to current displacements
        obj.state.data.dispCurr(ents) = obj.state.data.dispCurr(ents) + dSol(getDoF(obj.dofm,obj.field));
      end
      du = obj.state.data.dispCurr - obj.state.data.dispConv;
      % Update strain
      l = 0;
      for el=1:obj.mesh.nCells
        dof = getDoFID(obj.mesh,el);
        vtkId = obj.mesh.cellVTKType(el);
        elem = getElement(obj.elements,vtkId);
        nG = elem.GaussPts.nNode;
        N = getDerBasisFAndDet(elem,el,2);
        B = zeros(6,elem.nNode*obj.mesh.nDim,nG);
        B(elem.indB(:,2)) = N(elem.indB(:,1));
        obj.state.data.curr.strain((l+1):(l+nG),:) = ...
          reshape(pagemtimes(B,du(dof)),6,nG)';
        l = l + nG;
      end
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
          stateIn.data.curr.stress((l+1):(l+nG),:))/vol;
        dStrain = pagemtimes(B,stateIn.data.dispCurr(dof));
        dStrain = dStrain.*reshape(dJWeighed,1,1,[]);
        avStrain(el,:) = sum(dStrain,3)/vol;
        l = l + nG;
      end
    end

    function var = getState(obj,varargin)
      % input: state structure
      % output: current primary variable
      if isempty(varargin)
        var = obj.state.data.dispCurr;
      else
        stateIn = varargin{1};
        var = stateIn.data.dispCurr;
      end
    end

    function setState(obj,id,vals)
      % set values of the primary variable  
      obj.state.data.dispCurr(id) = vals; 
    end

    function [dof,vals] = getBC(obj,bc,id,t,~)
      dof = obj.getBCdofs(bc,id);
      vals = obj.getBCVals(bc,id,t);
    end

    function applyDirVal(obj,dof,vals)
      obj.state.data.dispConv(dof) = vals;
      obj.state.data.dispCurr(dof) = vals;
    end

    function computeRhs(obj,varargin)
      %computing rhs for poromechanics field
      if isLinear(obj) % linear case
        % update elastic stress
        subCells = obj.dofm.getFieldCells(obj.field);
        l1 = 0;
        for el=subCells'
          D = getElasticTensor(obj.material.getMaterial(obj.mesh.cellTag(el)).ConstLaw);
          % Get the right material stiffness for each element
          vtk = obj.mesh.cellVTKType(el);
          elem = obj.elements.getElement(vtk);
          nG = elem.GaussPts.nNode;
          obj.state.data.curr.stress((l1+1):(l1+nG),:) = ...
            obj.state.data.curr.stress((l1+1):(l1+nG),:)+...
            obj.state.data.curr.strain((l1+1):(l1+nG),:)*D;
          l1 = l1+nG;
        end
        ents = obj.dofm.getActiveEnts(obj.getField());
        if obj.simParams.isTimeDependent
          theta = obj.simParams.theta;
          obj.rhs = obj.J*obj.state.data.dispCurr(ents) + ...
            (1/theta-1)*obj.J*obj.state.data.dispConv(ents);
        else
          obj.rhs = obj.J*obj.state.data.dispCurr(ents);
        end
      else % non linear case: rhs computed with internal forces (B^T*sigma)
        obj.rhs = obj.fInt; % provisional assuming theta = 1;
      end
    end
% 
%     function blk = blockJacobian(obj,varargin)
%       fRow = varargin{1};
%       fCol = varargin{2};
%       locRow = obj.dofm.field2block(fRow);
%       locCol = obj.dofm.field2block(fCol);
%       blk = obj.simParams.theta*obj.K(locRow,locCol);
%     end
% 
%     function blk = blockRhs(obj, fld)
%       if ~strcmp(obj.dofm.subPhysics(fld), 'Poro')
%         % no contribution to non poro fields
%         blk = 0;
%       else
%         dofs = obj.dofm.field2block(fld);
%         blk = obj.rhs(dofs);
%       end
%     end

    function out = isLinear(obj)
      out = false;
      for i = 1:obj.mesh.nCellTag
        out = isa(obj.material.getMaterial(i).ConstLaw,"Elastic");
        if ~out
          return;
        end
      end
    end

    function [cellData,pointData] = printState(obj,sOld,sNew,t)
      % append state variable to output structure
      switch nargin
        case 2
          [stress,strain] = finalizeState(obj,sOld);
          displ = sOld.data.dispConv;
        case 4
          % linearly interpolate state variables containing print time
          fac = (t - sOld.t)/(sNew.t - sOld.t);
          [avStressOld,avStrainOld] = finalizeState(obj,sOld);
          [avStressNew,avStrainNew] = finalizeState(obj,sNew);
          stress = avStressNew*fac+avStressOld*(1-fac);
          strain = avStrainNew*fac+avStrainOld*(1-fac);
          displ = sNew.data.dispConv*fac+sOld.data.dispConv*(1-fac);
        otherwise
          error('Wrong number of input arguments');
      end
      [cellData,pointData] = Poromechanics.buildPrintStruct(displ,stress,strain);
    end
  end

  methods (Access=private)
    function dof = getBCdofs(obj,bc,id)
      switch bc.getCond(id)
        case 'NodeBC'
          ents = bc.getEntities(id);
        case 'SurfBC'
          ents = bc.getLoadedEntities(id);
          % node id contained by constrained surface
        otherwise
          error('BC type %s is not available for %s field',cond,obj.field);
      end
      % map entities dof to local dof numbering
      dof = obj.dofm.getLocalEnts(ents,obj.fldId);
      dof = bc.getCompEntities(id,dof);
    end

    function vals = getBCVals(obj,bc,id,t)
      vals = bc.getVals(id,t);
      if strcmp(bc.getCond(id),'SurfBC')
        entInfl = bc.getEntitiesInfluence(id);
        vals = entInfl*vals;
      end
    end

    function setPoromechanics(obj)
      obj.cell2stress = zeros(obj.mesh.nCells,1);
    end
  end

  methods (Static)
    function [cellStr,pointStr] = buildPrintStruct(disp,stress,strain)
      nCellData = 12;
      nPointData = 3;
      pointStr = repmat(struct('name', 1, 'data', 1), nPointData, 1);
      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      % Displacement
      pointStr(1).name = 'ux';
      pointStr(1).data = disp(1:3:end);
      pointStr(2).name = 'uy';
      pointStr(2).data = disp(2:3:end);
      pointStr(3).name = 'uz';
      pointStr(3).data = disp(3:3:end);
      %
      % Stress
      cellStr(1).name = 'sx';
      cellStr(1).data = stress(:,1);
      cellStr(2).name = 'sy';
      cellStr(2).data = stress(:,2);
      cellStr(3).name = 'sz';
      cellStr(3).data = stress(:,3);
      cellStr(4).name = 'txy';
      cellStr(4).data = stress(:,4);
      cellStr(5).name = 'tyz';
      cellStr(5).data = stress(:,5);
      cellStr(6).name = 'txz';
      cellStr(6).data = stress(:,6);
      %
      % Strain
      cellStr(7).name = 'ex';
      cellStr(7).data = strain(:,1);
      cellStr(8).name = 'ey';
      cellStr(8).data = strain(:,2);
      cellStr(9).name = 'ez';
      cellStr(9).data = strain(:,3);
      cellStr(10).name = 'gxy';
      cellStr(10).data = strain(:,4);
      cellStr(11).name = 'gyz';
      cellStr(11).data = strain(:,5);
      cellStr(12).name = 'gxz';
      cellStr(12).data = strain(:,6);
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
      out = Poromechanics.field;
    end
  end
end

