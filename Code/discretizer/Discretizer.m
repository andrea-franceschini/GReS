classdef Discretizer < handle           
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    H
    P
    J
    rhs
    blockJ
    blockRhs
    KPoro
    Q
  end
  
  properties (Access = public)
    model
    simParams
    dofm
    mesh
    elements
    faces
    material
%     bound
%     BCName
%     state
    GaussPts
%     probType  %either linear (lin) or nonlinear (nonlin)
%     flCompRHS = false
%     nE   % nE = [#tetra, #hexa, #wed, #pyr]
%     nEntryKLoc
%     fConst
    trans
    RHSGravTerm
    upElem
    isIntFaces
%     Sw
%     dSw
%     lw
%     dlw
    nEntryKLoc
  end
  
  methods (Access = public)
    function obj = Discretizer(symmod,simParams,dofManager,grid,mat,varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      obj.setDiscretizer(symmod,simParams,dofManager,grid,mat,varargin);
    end
    
    function trans = getFaceTransmissibilities(obj,faceID)
      trans = obj.trans(faceID);
    end
    
    
    
    
    
    function computeSPFMatrices(obj)
      if obj.model.isFEMBased('Flow')
        computeFlowStiffAndCapMatFEM_Test(obj);
      elseif obj.model.isFVTPFABased('Flow')
        mu = obj.material.getFluid().getDynViscosity();
        computeFlowStiffMatFV_Test(obj,1/mu);
        computeFlowCapMatFV_Test(obj);
      end
    end
    
    function computeVSFMatricesAndRHS(obj,statek,stateTmp,dt)
      pkpt = obj.simParams.theta*stateTmp.pressure + ...
        (1 - obj.simParams.theta)*statek.pressure;
      [Swkpt,dSwkpt,lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt);
      computeFlowStiffMatFV(obj,lwkpt);
      computeFlowCapMatFV(obj,Swkpt,dSwkpt);
      computeFlowJacobian(obj,dt,statek,stateTmp,pkpt,dSwkpt,dlwkpt);
      computeFlowRHS(obj,statek,stateTmp,dt,lwkpt);
    end
    
    function computeFlowStiffAndCapMatFEM(obj)
      % Compute the stiffness (H) and mass (P) matrices for the flow problem by FEM
      iiVec = zeros((obj.elements.nNodesElem.^2)*obj.elements.nCellsByType,1);
      jjVec = zeros((obj.elements.nNodesElem.^2)*obj.elements.nCellsByType,1);
      HVec = zeros((obj.elements.nNodesElem.^2)*obj.elements.nCellsByType,1);
      PVec = zeros((obj.elements.nNodesElem.^2)*obj.elements.nCellsByType,1);
      % Get the fluid compressibility
      beta = obj.material.getFluid().getFluidCompressibility();
      if obj.elements.nCellsByType(2) > 0
        N1 = obj.elements.hexa.getBasisFinGPoints();
      end
      % Get the fluid dynamic viscosity
      mu = obj.material.getFluid().getDynViscosity();
      %
      l1 = 0;
      for el=1:obj.mesh.nCells
        % Get the rock permeability, porosity and compressibility
        permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
        poro = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPorosity();
        if isPoromechanics(obj.model)
            alpha = 0; %this term is not needed in coupled formulation
        else
            alpha = obj.material.getMaterial(obj.mesh.cellTag(el)).ConstLaw.getRockCompressibility();
            %solid skeleton contribution to storage term as oedometric compressibility . 
        end
        % Compute the element matrices based on the element type
        % (tetrahedra vs. hexahedra)
        switch obj.mesh.cellVTKType(el)
          case 10 % Tetrahedra
            % Computing the H matrix contribution
            N = obj.elements.tetra.getDerBasisF(el);
%               vol = getVolume(obj.elements,el);
            HLoc = N'*permMat*N*obj.elements.vol(el)/mu;
            s1 = obj.elements.nNodesElem(1)^2;
            % Computing the P matrix contribution
            PLoc = ((alpha + poro*beta)*obj.elements.vol(el)/20)*(ones(obj.elements.nNodesElem(1))...
                    + eye(obj.elements.nNodesElem(1)));
          case 12 % Hexa
            [N,dJWeighed] = obj.elements.hexa.getDerBasisFAndDet(el,1);
            permMat = permMat/mu;
            Hs = pagemtimes(pagemtimes(N,'ctranspose',permMat,'none'),N);
            Hs = Hs.*reshape(dJWeighed,1,1,[]);
            HLoc = sum(Hs,3);
            clear Hs;
            s1 = obj.elements.nNodesElem(2)^2;
            % Computing the P matrix contribution
            PLoc = (alpha+poro*beta)*(N1'*diag(dJWeighed)*N1);
        end
        %
        dof = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
        [jjLoc,iiLoc] = meshgrid(dof,dof);
        iiVec(l1+1:l1+s1) = iiLoc(:);
        jjVec(l1+1:l1+s1) = jjLoc(:);
        HVec(l1+1:l1+s1) = HLoc(:);
        PVec(l1+1:l1+s1) = PLoc(:);
        l1 = l1 + s1;
      end
      % Assemble H and P matrices
      obj.H = sparse(iiVec,jjVec,HVec,obj.mesh.nNodes,obj.mesh.nNodes);
      obj.P = sparse(iiVec,jjVec,PVec,obj.mesh.nNodes,obj.mesh.nNodes);
    end

    function computeFlowStiffAndCapMatFEM_Test(obj) %provisional method exploiting dof manager workflow
      nSubPh = length(obj.dofm.subList);
      nSub = length(obj.dofm.subDomains);
      %nBlock = size(obj.blockJ,1);
      for i = 1:nSub
         if any(strcmp(obj.dofm.subDomains(i).physics,'Flow'))  
              subReg = obj.dofm.subDomains(i).regions; %region Tag for block's subdomain
              subCells = find(ismember(obj.mesh.cellTag,subReg)); %id of cells in subdomain
              %nSubCells = length(subCells); %number of cells in subdomain
              nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]); 
              % Compute the stiffness (H) and mass (P) matrices for the flow problem by FEM
              iiVec = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
              jjVec = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
              phVeci = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
              phVecj = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
              HVec = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
              PVec = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
              % Get the fluid compressibility
              beta = obj.material.getFluid().getFluidCompressibility();
              if nSubCellsByType(2) > 0
                N1 = obj.elements.hexa.getBasisFinGPoints();
              end
              % Get the fluid dynamic viscosity
              mu = obj.material.getFluid().getDynViscosity();
              %
              l1 = 0;
              for el = subCells'
                % Get the rock permeability, porosity and compressibility
                permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
                poro = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPorosity();
                if  any(strcmp(obj.dofm.subDomains(i).physics,'Poro'))
                    alpha = 0; %this term is not needed in coupled formulation
                else
                    alpha = obj.material.getMaterial(obj.mesh.cellTag(el)).ConstLaw.getRockCompressibility();
                    %solid skeleton contribution to storage term as oedometric compressibility . 
                end
                % Compute the element matrices based on the element type
                % (tetrahedra vs. hexahedra)
                switch obj.mesh.cellVTKType(el)
                  case 10 % Tetrahedra
                    % Computing the H matrix contribution
                    N = obj.elements.tetra.getDerBasisF(el);
        %               vol = getVolume(obj.elements,el);
                    HLoc = N'*permMat*N*obj.elements.vol(el)/mu;
                    s1 = obj.elements.nNodesElem(1)^2;
                    % Computing the P matrix contribution
                    PLoc = ((alpha + poro*beta)*obj.elements.vol(el)/20)*(ones(obj.elements.nNodesElem(1))...
                            + eye(obj.elements.nNodesElem(1)));
                  case 12 % Hexa
                    [N,dJWeighed] = obj.elements.hexa.getDerBasisFAndDet(el,1);
                    permMat = permMat/mu;
                    Hs = pagemtimes(pagemtimes(N,'ctranspose',permMat,'none'),N);
                    Hs = Hs.*reshape(dJWeighed,1,1,[]);
                    HLoc = sum(Hs,3);
                    clear Hs;
                    s1 = obj.elements.nNodesElem(2)^2;
                    % Computing the P matrix contribution
                    PLoc = (alpha+poro*beta)*(N1'*diag(dJWeighed)*N1);
                end
                %Getting dof associated to Flow subphysic
                dof = obj.dofm.getLocDoF('Flow',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                ph = obj.dofm.getSubPhysic('Flow',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                %dof = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
                [jjLoc,iiLoc] = meshgrid(dof,dof);
                [phLocj,phLoci] = meshgrid(ph,ph);
                iiVec(l1+1:l1+s1) = iiLoc(:);
                jjVec(l1+1:l1+s1) = jjLoc(:);
                phVeci(l1+1:l1+s1) = phLoci(:);
                phVecj(l1+1:l1+s1) = phLocj(:);
                HVec(l1+1:l1+s1) = HLoc(:);
                PVec(l1+1:l1+s1) = PLoc(:);
                l1 = l1 + s1;
              end
              % Assemble H and P matrices defined as new fields of
              % block Jacobian matrix
              k=1;
              % Populate the Jacobian
              for ii = 1:nSubPh
                  for jj = k:nSubPh
                      if all(strcmp(obj.blockJ(ii,jj).physics,'Flow'))
                      %find pairs of row,col indices belonging to ii,jj
                      %block of Jacobian
                      tmp = intersect(find(phVeci == ii),find(phVecj == jj));
                      obj.blockJ(ii,jj).H = obj.blockJ(ii,jj).H + sparse(iiVec(tmp),jjVec(tmp),HVec(tmp),obj.dofm.numDof(ii),obj.dofm.numDof(jj)); 
                      obj.blockJ(ii,jj).P = obj.blockJ(ii,jj).P + sparse(iiVec(tmp),jjVec(tmp),PVec(tmp),obj.dofm.numDof(ii),obj.dofm.numDof(jj)); 
                      %consider a transpose option for the following lines
                      if ii ~= jj
                        obj.blockJ(jj,ii).H = (obj.blockJ(ii,jj).H)';
                        obj.blockJ(jj,ii).P = (obj.blockJ(ii,jj).P)';
                      end
                      end
                  end
                  k = k+1;  %loop only trough upper triangular
              end
          end
      end
    end



    function computeFlowStiffMatFV(obj,lw)
      % Inspired by MRST
      % Pairs of cells sharing internal faces
      neigh1 = obj.faces.faceNeighbors(obj.isIntFaces,1);
      neigh2 = obj.faces.faceNeighbors(obj.isIntFaces,2);
      % Transmissibility of internal faces
      tmpVec = lw.*obj.trans(obj.isIntFaces);
%       tmpVec = lw.*tmpVec;
      sumDiagTrans = accumarray([neigh1; neigh2], ...
        repmat(tmpVec,[2,1]),[obj.mesh.nCells,1]);
      obj.H = sparse([neigh1; neigh2; (1:obj.mesh.nCells)'], ...
                     [neigh2; neigh1; (1:obj.mesh.nCells)'], ...
                     [-tmpVec; -tmpVec; ...
                      sumDiagTrans],obj.mesh.nCells,obj.mesh.nCells);
    end
    
      


    function computeFlowStiffMatFV_Test(obj,lw)
      % Inspired by MRST
      nSubPh = length(obj.dofm.subList);
      for i = 1:length(obj.dofm.subDomains)
          if any(strcmp(obj.dofm.subDomains(i).physics,'Flow'))
          subCells = find(obj.dofm.subCells(:,i));
          nSubCells = length(subCells); 
          %get pairs of faces that contribute to the subdomain
          neigh1 = obj.faces.faceNeighbors(obj.isIntFaces(:,i),1);
          neigh2 = obj.faces.faceNeighbors(obj.isIntFaces(:,i),2);
          neigh1dof = obj.dofm.getLocDoF('Flow',neigh1);
          neigh2dof = obj.dofm.getLocDoF('Flow',neigh2);
          neigh1ph =  obj.dofm.getSubPhysic('Flow',neigh1);
          neigh2ph =  obj.dofm.getSubPhysic('Flow',neigh2);
          % Transmissibility of internal faces
          tmpVec = lw.*obj.trans(obj.isIntFaces(:,i));
          % tmpVec = lw.*tmpVec;
          sumDiagTrans = accumarray([neigh1dof; neigh2dof], ...
            repmat(tmpVec,[2,1]),[nSubCells,1]);
          % Assemble H matrices defined as new fields of block Jacobian matrix
          k=1;
          % Populate the Jacobian
          for ii = 1:nSubPh
              for jj = k:nSubPh
                  if all(strcmp(obj.blockJ(ii,jj).physics,'Flow'))
                  %find pairs of row,col indices belonging to ii,jj
                  %block of Jacobian
                  tmp = intersect(find(neigh1ph == ii),find(neigh2ph == jj));
                  if ii == jj
                      %diagonal block involving diagonal terms
                      obj.blockJ(ii,jj).H  = obj.blockJ(ii,jj).H + sparse([neigh1dof(tmp); neigh2dof(tmp); (1:nSubCells)'], ...
                         [neigh2dof(tmp); neigh1dof(tmp); (1:nSubCells)'], ...
                         [-tmpVec(tmp); -tmpVec(tmp); ...
                          sumDiagTrans],obj.dofm.numDof(ii),obj.dofm.numDof(jj));
                  else
                      %extra-diagonal block
                      obj.blockJ(ii,jj).H  = obj.blockJ(ii,jj).H + sparse([neigh1dof(tmp); neigh2dof(tmp)], ...
                         [neigh2dof(tmp); neigh1dof(tmp)], ...
                         [-tmpVec(tmp); -tmpVec(tmp)],obj.dofm.numDof(ii),obj.dofm.numDof(jj));     
                      %exploiting symmetry of flow stiffness matrix
                      obj.blockJ(jj,ii).H = (obj.blockJ(ii,jj).H)';
                  end
                  end
              end
          end

          end
      end
    end
    
    function computeFlowCapMatFV(obj,varargin)
      % if isVariabSatFlow
      % varargin{1} -> Sw
      % varargin{2} -> dSw
      poroMat = zeros(obj.mesh.nCellTag,1);
      alphaMat = zeros(obj.mesh.nCellTag,1);
      beta = obj.material.getFluid().getFluidCompressibility();
      for m = 1:obj.mesh.nCellTag
        poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
        alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
      end
      % (alpha+poro*beta)
      PVal = alphaMat(obj.mesh.cellTag) + beta*poroMat(obj.mesh.cellTag);
      if obj.model.isVariabSatFlow()
        PVal = PVal.*varargin{1} + poroMat(obj.mesh.cellTag).*varargin{2};
      end
      PVal = PVal.*obj.elements.vol;
      obj.P = sparse(1:obj.mesh.nCells,1:obj.mesh.nCells,PVal,obj.mesh.nCells,obj.mesh.nCells);
    end

    function computeFlowCapMatFV_Test(obj,varargin)
      % if isVariabSatFlow
      % varargin{1} -> Sw
      % varargin{2} -> dSw
      nSubPh = length(obj.dofm.subList);
      for i = 1:length(obj.dofm.subDomains)
          if any(strcmp(obj.dofm.subDomains(i).physics,'Flow'))
              subCells = find(obj.dofm.subCells(:,i));
              nSubCells = length(subCells); 
              poroMat = zeros(nSubCells,1);
              alphaMat = zeros(nSubCells,1);
              beta = obj.material.getFluid().getFluidCompressibility();
              for m = 1:obj.mesh.nCellTag
                if ~any(strcmp(obj.dofm.subDomains(i).physics,'Poro'))
                    % compute alpha only if there's no coupling in the
                    % subdomain 
                    alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
                end
                poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
              end
              % (alpha+poro*beta)
              PVal = alphaMat(obj.mesh.cellTag(subCells)) + beta*poroMat(obj.mesh.cellTag(subCells));
              if obj.model.isVariabSatFlow()
                PVal = PVal.*varargin{1} + poroMat(obj.mesh.cellTag(subCells)).*varargin{2};
              end
              PVal = PVal.*obj.elements.vol(subCells);
              for ii = 1:nSubPh
                if all(strcmp(obj.blockJ(ii,ii).physics,'Flow'))
                 obj.blockJ(ii,ii).P = obj.blockJ(ii,ii).P + sparse(1:nSubCells,1:nSubCells,PVal,nSubCells,nSubCells);
                end
              end
          end
      end
    end

    function computeFlowJacobian(obj,dt,varargin)
      % varargin -> statek,stateTmp,pkpt,dSwkpt,dlwkpt
      % IF SINGLE PHASE
      obj.J = obj.simParams.theta*obj.H + obj.P/dt;
      if obj.model.isVariabSatFlow() && obj.simParams.isNewtonNLSolver()
        % Compute the additional terms of the Jacobian
        [JNewt] = computeNewtPartOfJacobian(obj,dt,varargin{1}, ...
          varargin{2},varargin{3},varargin{4},varargin{5});
        obj.J = obj.J + JNewt;
      end
    end

    function computeFlowJacobian_Test(obj,dt,varargin)
      nBlock = size(obj.blockJ,1);
          for i=1:nBlock^2
              if all(strcmp(obj.blockJ(i).physics,'Flow')) 
                  % varargin -> statek,stateTmp,pkpt,dSwkpt,dlwkpt
                  % IF SINGLE PHASE
                  obj.blockJ(i).block = obj.simParams.theta*obj.blockJ(i).H + obj.blockJ(i).P/dt;
                  if obj.model.isVariabSatFlow() && obj.simParams.isNewtonNLSolver()
                    % Compute the additional terms of the Jacobian
                    [JNewt] = computeNewtPartOfJacobian(obj,dt,varargin{1}, ...
                      varargin{2},varargin{3},varargin{4},varargin{5});
                    obj.J = obj.J + JNewt;
                  end
              end
          end
    end

    function computeFlowRHS_test(obj,statek,stateTmp,dt,lw)
      % Compute the residual of the flow problem
%       obj.rhs = (obj.simParams.theta*obj.H + obj.P/dt)*stateTmp.pressure ...
%         - (obj.P/dt - (1 - obj.simParams.theta)*obj.H)*statek.pressure;
     %loop trough rhs blocks to find coupled poromechanics residual 
      nRhs = length(obj.blockRhs);
      theta = obj.simParams.theta;
      for i = 1:nRhs
          if strcmp(obj.blockRhs(i).physic,'Flow')
              for j = 1:nRhs
                  if all(strcmp(obj.blockJ(i,j).physics,'Flow'))
                    dofs = obj.dofm.getLocDoF(j);  
                    obj.blockRhs(i).block = obj.blockRhs(i).block + ...
                        theta*obj.blockJ(i,j).H*stateTmp.pressure(dofs) + (1-theta)*obj.blockJ(i,j).H*statek.pressure(dofs) + ...
                        (obj.blockJ(i,j).P/dt)*(stateTmp.pressure(dofs) - statek.pressure(dofs));
                  end
              end
              gamma = obj.material.getFluid().getFluidSpecWeight();
              %adding gravity rhs contribute
              if gamma > 0
                if isFEMBased(obj.model,'Flow')
                 obj.blockRhs(i).block = obj.blockRhs(i).block + obj.blockRhs(i).rhsGrav;
                elseif isFVTPFABased(obj.model,'Flow')
                  obj.blockRhs(i).block = obj.blockRhs(i).block + finalizeRHSGravTerm(obj,i,lw);
                end
              end
          end
      end
      % obj.rhs = obj.simParams.theta*obj.H*stateTmp.pressure + 1/dt*obj.P*stateTmp.pressure ...
      %   - 1/dt*obj.P*statek.pressure + (1 - obj.simParams.theta)*obj.H*statek.pressure;
    end

    function computeFlowRHS(obj,statek,stateTmp,dt,lw)
      % Compute the residual of the flow problem
%       obj.rhs = (obj.simParams.theta*obj.H + obj.P/dt)*stateTmp.pressure ...
%         - (obj.P/dt - (1 - obj.simParams.theta)*obj.H)*statek.pressure;
      obj.rhs = obj.simParams.theta*obj.H*stateTmp.pressure + 1/dt*obj.P*stateTmp.pressure ...
        - 1/dt*obj.P*statek.pressure + (1 - obj.simParams.theta)*obj.H*statek.pressure;
      gamma = obj.material.getFluid().getFluidSpecWeight();
      if gamma > 0
        if isFEMBased(obj.model,'Flow')
          obj.rhs = obj.rhs + obj.RHSGravTerm.*lw;
        elseif isFVTPFABased(obj.model,'Flow')
          obj.rhs = obj.rhs + finalizeRHSGravTerm(obj,lw);
        end
      end
    end

    function computePoroRhs(obj,stateTmp) 
      %loop trough rhs blocks to accumulate coupled poromechanics residual 
      nRhs = length(obj.blockRhs);
      theta = obj.simParams.theta;
      for i = 1:nRhs
          if strcmp(obj.blockRhs(i).physic,'Poro')
              for j = 1:nRhs
                  if all(strcmp(obj.blockJ(i,j).physics,'Poro'))
                    dofs = obj.dofm.getLocDoF(j);  
                    obj.blockRhs(i).block = obj.blockRhs(i).block + ...
                        obj.blockJ(i,j).block*stateTmp.dispCurr(dofs) + ...
                        ((1-theta)/theta)*obj.blockJ(i,j).block*stateTmp.dispConv(dofs);
                  end
              end
          end
      end
    end


    
    function computePoroSyst(obj,state,dt)
      % Compute the Jacobian and residual of the geomechanical problem
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
%       KLocSize = obj.mesh.cellNumVerts(1)*obj.mesh.nDim;
      obj.rhs = zeros(obj.mesh.nNodes*obj.mesh.nDim,1);
      iiVec = zeros(obj.nEntryKLoc*obj.elements.nCellsByType,1);
      jjVec = zeros(obj.nEntryKLoc*obj.elements.nCellsByType,1);
      JVec = zeros(obj.nEntryKLoc*obj.elements.nCellsByType,1);
      %
      l1 = 0;
      l2 = 0;
      for el=1:obj.mesh.nCells
        % Get the right material stiffness for each element
        switch obj.mesh.cellVTKType(el)
          case 10 % Tetrahedra
            N = getDerBasisF(obj.elements.tetra,el);
            vol = findVolume(obj.elements.tetra,el);
            B = zeros(6,4*obj.mesh.nDim);
            B(obj.elements.indB(1:36,2)) = N(obj.elements.indB(1:36,1));
            [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
                 state.conv.stress(l2+1,:), ...
                 state.curr.strain(l2+1,:), ...
                 dt,state.conv.status(l2+1,:), el, state.t);
             state.curr.status(l2+1,:) = status;
             state.curr.stress(l2+1,:) = sigma;
           % D = obj.preP.getStiffMatrix(el,state.stress(l2+1,3)+state.iniStress(l2+1,3));
            JLoc = B'*D*B*vol;
            s1 = obj.nEntryKLoc(1);
            %
%             if obj.flCompRHS
              sz = sigma - state.iniStress(l2+1,:);
              fLoc = (B')*sz'*vol;
              s2 = 1;
%             end
          case 12 % Hexahedra
              [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
              B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
              B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
              [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
                 state.conv.stress(l2+1:l2+obj.GaussPts.nNode,:), ...
                 state.curr.strain(l2+1:l2+obj.GaussPts.nNode,:), ...
                 dt,state.conv.status(l2+1:l2+obj.GaussPts.nNode,:), el, state.t);
              state.curr.status(l2+1:l2+obj.GaussPts.nNode,:) = status;
              state.curr.stress((l2+1):(l2+obj.GaussPts.nNode),:) = sigma;
              Js = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
              Js = Js.*reshape(dJWeighed,1,1,[]);
              JLoc = sum(Js,3);
              clear Js;
              s1 = obj.nEntryKLoc(2);
              sz = sigma - state.iniStress(l2+1:l2+obj.GaussPts.nNode,:);
              sz = reshape(sz',6,1,obj.GaussPts.nNode);
              fTmp = pagemtimes(B,'ctranspose',sz,'none');
              fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
              fLoc = sum(fTmp,3);
              s2 = obj.GaussPts.nNode;
%               D = obj.preP.getStiffMatrix(el,state.stress(l2+1:l2+obj.GaussPts.nNode,3)+state.iniStress(l2+1:l2+obj.GaussPts.nNode,3));
%               Js = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
%               Js = Js.*reshape(dJWeighed,1,1,[]);
%               JLoc = sum(Js,3);
%               clear Ks;
%               s1 = obj.nEntryKLoc(2);
%               %
% %               if obj.flCompRHS
%                 sz = state.stress(l2+1:l2+obj.GaussPts.nNode,:);
%                 sz = reshape(sz',6,1,obj.GaussPts.nNode);
%                 fTmp = pagemtimes(B,'ctranspose',sz,'none');
%                 fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
%                 fLoc = sum(fTmp,3);
%                 s2 = obj.GaussPts.nNode;
% %               end
        end
        %
        dof = getDoFID(obj.mesh,el);
        [jjLoc,iiLoc] = meshgrid(dof,dof);
        iiVec(l1+1:l1+s1) = iiLoc(:);
        jjVec(l1+1:l1+s1) = jjLoc(:);
        JVec(l1+1:l1+s1) = JLoc(:);
        % Accumulate the residual contributions
          obj.rhs(dof) = obj.rhs(dof) + fLoc;
%         end
        l1 = l1 + s1;
        l2 = l2 + s2;
      end
      % Populate the Jacobian
      obj.J = sparse(iiVec,jjVec,JVec,obj.mesh.nNodes*obj.mesh.nDim, ...
                     obj.mesh.nNodes*obj.mesh.nDim);
    end

    function computePoroSyst_Test(obj,state,dt)
    % Compute the Jacobian and residual of the geomechanical problem
      %nBlock = size(obj.blockJ,1);
      nSubPh = length(obj.dofm.subList);
      % Loop trough subdomains
      for i = 1:length(obj.dofm.subDomains)
          if any(strcmp(obj.dofm.subDomains(i).physics,'Poro')) && ~obj.dofm.subDomains(i).coupling   
              subReg = obj.dofm.subDomains(i).regions; %region for block's subdomain
              subCells = find(ismember(obj.mesh.cellTag,subReg)); %id of cells in subdomain
              nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]); 
              iiVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
              jjVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
              phVeci = zeros(obj.nEntryKLoc*nSubCellsByType,1);
              phVecj = zeros(obj.nEntryKLoc*nSubCellsByType,1);
              KVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
              %
              l1 = 0;
              l2 = 0;
              for el=subCells'
                % Get the right material stiffness for each element
                switch obj.mesh.cellVTKType(el)
                  case 10 % Tetrahedra
                    N = getDerBasisF(obj.elements.tetra,el);
                    vol = findVolume(obj.elements.tetra,el);
                    B = zeros(6,4*obj.mesh.nDim);
                    B(obj.elements.indB(1:36,2)) = N(obj.elements.indB(1:36,1));
                    [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
                         state.conv.stress(l2+1,:), ...
                         state.curr.strain(l2+1,:), ...
                         dt,state.conv.status(l2+1,:), el, state.t);
                     state.curr.status(l2+1,:) = status;
                     state.curr.stress(l2+1,:) = sigma;
                   % D = obj.preP.getStiffMatrix(el,state.stress(l2+1,3)+state.iniStress(l2+1,3));
                    KLoc = B'*D*B*vol;
                    s1 = obj.nEntryKLoc(1);
                    %
        %             if obj.flCompRHS
                      sz = sigma - state.iniStress(l2+1,:);
                      fLoc = (B')*sz'*vol;
                      s2 = 1;
        %             end
                  case 12 % Hexahedra
                      [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                      B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
                      B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                      [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
                         state.conv.stress(l2+1:l2+obj.GaussPts.nNode,:), ...
                         state.curr.strain(l2+1:l2+obj.GaussPts.nNode,:), ...
                         dt,state.conv.status(l2+1:l2+obj.GaussPts.nNode,:), el, state.t);
                      state.curr.status(l2+1:l2+obj.GaussPts.nNode,:) = status;
                      state.curr.stress((l2+1):(l2+obj.GaussPts.nNode),:) = sigma;
                      Ks = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
                      Ks = Ks.*reshape(dJWeighed,1,1,[]);
                      KLoc = sum(Ks,3);
                      clear Ks;
                      s1 = obj.nEntryKLoc(2);
                      sz = sigma - state.iniStress(l2+1:l2+obj.GaussPts.nNode,:);
                      sz = reshape(sz',6,1,obj.GaussPts.nNode);
                      fTmp = pagemtimes(B,'ctranspose',sz,'none');
                      fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
                      fLoc = sum(fTmp,3);
                      s2 = obj.GaussPts.nNode;
                end
                % get local dof (refered to local blocks of Jacobian and
                % rhs)
                dof = obj.dofm.getLocDoF('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                ph = obj.dofm.getSubPhysic('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                [jjLoc,iiLoc] = meshgrid(dof,dof);
                [phLocj,phLoci] = meshgrid(ph,ph);
                iiVec(l1+1:l1+s1) = iiLoc(:);
                jjVec(l1+1:l1+s1) = jjLoc(:);
                phVeci(l1+1:l1+s1) = phLoci(:);
                phVecj(l1+1:l1+s1) = phLocj(:);
                KVec(l1+1:l1+s1) = KLoc(:);
                % Accumulate the residual contributions trough all residual blocks 
                % theta method not required for uncoupled Poromechanics
                for j = 1:nSubPh
                    if strcmp(obj.blockRhs(j).physic,'Poro')
                    %get Rhs ID corresponding to Poro subphysics of subdomain i
                    tmp = ph == j;
                    obj.blockRhs(j).block(dof(tmp)) = obj.blockRhs(j).block(dof(tmp)) + fLoc(tmp);
                    end
                end
        %         end
                l1 = l1 + s1;
                l2 = l2 + s2;
              end
              k=1;
              % Populate the Jacobian
              for ii = 1:nSubPh
                  for jj = k:nSubPh
                      if all(strcmp(obj.blockJ(ii,jj).physics,'Poro'))
                      %find pairs of row,col indices belonging to ii,jj
                      %block of Jacobian
                      tmp = intersect(find(phVeci == ii),find(phVecj == jj));
                      obj.blockJ(ii,jj).block = sparse(iiVec(tmp),jjVec(tmp),KVec(tmp), ...
                          obj.dofm.numDof(ii),obj.dofm.numDof(jj)); 
                      %consider a transpose option for the following line
                      obj.blockJ(jj,ii).block = (obj.blockJ(ii,jj).block)';
                      end
                  end
                  k = k+1;  %loop only trough upper triangular
              end
          end
      end
    end

    function computePoroCoupled_Test(obj,state,dt)
    % Compute the Jacobian and residual of the geomechanical problem
      %nBlock = size(obj.blockJ,1);
      nSubPh = length(obj.dofm.subList);
      nSub = length(obj.dofm.subDomains);
      % Loop trough elements of blockJ matrix to find pairs PORO-PORO
      for i = 1:nSub
          if any(strcmp(obj.dofm.subDomains(i).physics,'Poro'))   
              %subReg = obj.dofm.subDomains(i).regions; %region for block's subdomain
              subCells = find(obj.dofm.subCells(:,i));
              %subCells = find(ismember(obj.mesh.cellTag,subReg)); %id of cells in subdomain
              nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]); 
              iiVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
              jjVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
              phVeci = zeros(obj.nEntryKLoc*nSubCellsByType,1);
              phVecj = zeros(obj.nEntryKLoc*nSubCellsByType,1);
              KVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
              %
              l1 = 0;
              l2 = 0;
              for el=subCells'
                % Get the right material stiffness for each element
                switch obj.mesh.cellVTKType(el)
                  case 10 % Tetrahedra
                    N = getDerBasisF(obj.elements.tetra,el);
                    vol = findVolume(obj.elements.tetra,el);
                    B = zeros(6,4*obj.mesh.nDim);
                    B(obj.elements.indB(1:36,2)) = N(obj.elements.indB(1:36,1));
                    [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
                         state.conv.stress(l2+1,:), ...
                         state.curr.strain(l2+1,:), ...
                         dt,state.conv.status(l2+1,:), el, state.t);
                     state.curr.status(l2+1,:) = status;
                     state.curr.stress(l2+1,:) = sigma;
                   % D = obj.preP.getStiffMatrix(el,state.stress(l2+1,3)+state.iniStress(l2+1,3));
                    KLoc = B'*D*B*vol;
                    s1 = obj.nEntryKLoc(1);
                    %
        %             if obj.flCompRHS
                      sz = sigma - state.iniStress(l2+1,:);
                      fLoc = (B')*sz'*vol;
                      s2 = 1;
        %             end
                  case 12 % Hexahedra
                      [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                      B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
                      B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                      [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
                         state.conv.stress(l2+1:l2+obj.GaussPts.nNode,:), ...
                         state.curr.strain(l2+1:l2+obj.GaussPts.nNode,:), ...
                         dt,state.conv.status(l2+1:l2+obj.GaussPts.nNode,:), el, state.t);
                      state.curr.status(l2+1:l2+obj.GaussPts.nNode,:) = status;
                      state.curr.stress((l2+1):(l2+obj.GaussPts.nNode),:) = sigma;
                      Ks = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
                      Ks = Ks.*reshape(dJWeighed,1,1,[]);
                      KLoc = sum(Ks,3);
                      clear Ks;
                      s1 = obj.nEntryKLoc(2);
                      sz = sigma - state.iniStress(l2+1:l2+obj.GaussPts.nNode,:);
                      sz = reshape(sz',6,1,obj.GaussPts.nNode);
                      % fTmp = pagemtimes(B,'ctranspose',sz,'none');
                      % fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
                      % fLoc = sum(fTmp,3);
                      s2 = obj.GaussPts.nNode;
                end
                % get local dof (refered to local blocks of Jacobian and
                % rhs)
                dof = obj.dofm.getLocDoF('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                ph = obj.dofm.getSubPhysic('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                [jjLoc,iiLoc] = meshgrid(dof,dof);
                [phLocj,phLoci] = meshgrid(ph,ph);
                iiVec(l1+1:l1+s1) = iiLoc(:);
                jjVec(l1+1:l1+s1) = jjLoc(:);
                phVeci(l1+1:l1+s1) = phLoci(:);
                phVecj(l1+1:l1+s1) = phLocj(:);
                KVec(l1+1:l1+s1) = KLoc(:);
                l1 = l1 + s1;
                l2 = l2 + s2;
              end
              k=1;
              % Populate the Jacobian
              for ii = 1:nSubPh
                  for jj = k:nSubPh
                      if all(strcmp(obj.blockJ(ii,jj).physics,'Poro'))
                      %find pairs of row,col indices belonging to ii,jj
                      %block of Jacobian
                      tmp = intersect(find(phVeci == ii),find(phVecj == jj));
                      obj.blockJ(ii,jj).block = obj.blockJ(ii,jj).block + obj.simParams.theta*sparse(iiVec(tmp),jjVec(tmp),KVec(tmp), ...
                          obj.dofm.numDof(ii),obj.dofm.numDof(jj)); 
                      %consider a transpose option for the following line
                      if ii ~= jj
                        obj.blockJ(jj,ii).block = (obj.blockJ(ii,jj).block)';
                      end
                      end
                  end
                  k = k+1;  %loop only trough upper triangular
              end
          end
      end
    end

    % function computePoroCoupled_Test(obj,state,dt)
    % % Compute the Jacobian and residual of the geomechanical problem
    %   nBlock = size(obj.blockJ,1);
    %   % Loop trough elements of blockJ matrix to find pairs PORO-PORO
    %   for i=1:nBlock
    %       for j=1:nBlock
    %           if all(strcmp(obj.blockJ(i,j).physics,'Poro')) && obj.blockJ(i,j).coupling   
    %               subReg = obj.dofm.subDomains(obj.blockJ(i,j).subID).regions; %region for block's subdomain
    %               subCells = find(ismember(obj.mesh.cellTag,subReg)); %id of cells in subdomain
    %               nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]); 
    %               % Compute the stiffness (H) and mass (P) matrices for the flow problem by FEM
    %               iiVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
    %               jjVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
    %               KVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
    %               %
    %               l1 = 0;
    %               l2 = 0;
    %               for el=subCells'
    %                 % Get the right material stiffness for each element
    %                 switch obj.mesh.cellVTKType(el)
    %                   case 10 % Tetrahedra
    %                     N = getDerBasisF(obj.elements.tetra,el);
    %                     vol = findVolume(obj.elements.tetra,el);
    %                     B = zeros(6,4*obj.mesh.nDim);
    %                     B(obj.elements.indB(1:36,2)) = N(obj.elements.indB(1:36,1));
    %                     [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
    %                          state.conv.stress(l2+1,:), ...
    %                          state.curr.strain(l2+1,:), ...
    %                          dt,state.conv.status(l2+1,:), el, state.t);
    %                      state.curr.status(l2+1,:) = status;
    %                      state.curr.stress(l2+1,:) = sigma;
    %                    % D = obj.preP.getStiffMatrix(el,state.stress(l2+1,3)+state.iniStress(l2+1,3));
    %                     KLoc = B'*D*B*vol;
    %                     s1 = obj.nEntryKLoc(1);
    %                     %
    %         %             if obj.flCompRHS
    %                       sz = sigma - state.iniStress(l2+1,:);
    %                       fLoc = (B')*sz'*vol;
    %                       s2 = 1;
    %         %             end
    %                   case 12 % Hexahedra
    %                       [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
    %                       B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
    %                       B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
    %                       [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
    %                          state.conv.stress(l2+1:l2+obj.GaussPts.nNode,:), ...
    %                          state.curr.strain(l2+1:l2+obj.GaussPts.nNode,:), ...
    %                          dt,state.conv.status(l2+1:l2+obj.GaussPts.nNode,:), el, state.t);
    %                       state.curr.status(l2+1:l2+obj.GaussPts.nNode,:) = status;
    %                       state.curr.stress((l2+1):(l2+obj.GaussPts.nNode),:) = sigma;
    %                       Ks = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
    %                       Ks = Ks.*reshape(dJWeighed,1,1,[]);
    %                       KLoc = sum(Ks,3);
    %                       clear Ks;
    %                       s1 = obj.nEntryKLoc(2);
    %                       sz = sigma - state.iniStress(l2+1:l2+obj.GaussPts.nNode,:);
    %                       sz = reshape(sz',6,1,obj.GaussPts.nNode);
    %                       fTmp = pagemtimes(B,'ctranspose',sz,'none');
    %                       fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
    %                       fLoc = sum(fTmp,3);
    %                       s2 = obj.GaussPts.nNode;
    %                 end
    %                 %
    %                 dof = obj.dofm.getLocDoF('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
    %                 [jjLoc,iiLoc] = meshgrid(dof,dof);
    %                 iiVec(l1+1:l1+s1) = iiLoc(:);
    %                 jjVec(l1+1:l1+s1) = jjLoc(:);
    %                 KVec(l1+1:l1+s1) = KLoc(:);
    %                 % Accumulate the residual contributions - 
    %                 % theta method not required for uncoupled Poromechanics
    %                 % obj.blockRhs(i).block(dof) = obj.blockRhs(i).block(dof) + fLoc;
    %         %         end
    %                 l1 = l1 + s1;
    %                 l2 = l2 + s2;
    %               end
    %               % Populate the Jacobian 
    %               % Coupled Poro Block is already multiplied by theta
    %               obj.blockJ(i,j).block = obj.simParams.theta*sparse(iiVec,jjVec,KVec,obj.mesh.nNodes*obj.mesh.nDim, ...
    %                              obj.mesh.nNodes*obj.mesh.nDim);
    %           end
    %       end
    %   end
    % end
    % 
    
     function computePoroCoupled(obj,state,dt)
      % Compute only the Jacobian of geomechanical problem inside Coupled
      % HydroMechanics. At the moment the rhs is calculated separately
      % with the product KPoro*stateTmp.displ. Will be fixed in future
      %obj.rhs = zeros(obj.mesh.nNodes*obj.mesh.nDim,1);
      iiVec = zeros(obj.nEntryKLoc*obj.elements.nCellsByType,1);
      jjVec = zeros(obj.nEntryKLoc*obj.elements.nCellsByType,1);
      KVec = zeros(obj.nEntryKLoc*obj.elements.nCellsByType,1);
      %
      l1 = 0;
      l2 = 0;
      for el=1:obj.mesh.nCells
        % Get the right material stiffness for each element
        switch obj.mesh.cellVTKType(el)
          case 10 % Tetrahedraobj.dof
           N = getDerBasisF(obj.elements.tetra,el);
            vol = findVolume(obj.elements.tetra,el);
            B = zeros(6,4*obj.mesh.nDim);
            B(obj.elements.indB(1:36,2)) = N(obj.elements.indB(1:36,1));
            [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
                 state.conv.stress(l2+1,:), ...
                 state.curr.strain(l2+1,:), ...
                 dt,state.conv.status(l2+1,:), el, state.t);
             state.curr.status(l2+1,:) = status;
             state.curr.stress(l2+1,:) = sigma;
           % D = obj.preP.getStiffMatrix(el,state.stress(l2+1,3)+state.iniStress(l2+1,3));
            KLoc = B'*D*B*vol;
            s1 = obj.nEntryKLoc(1);
            %
%             if obj.flCompRHS
              sz = sigma - state.iniStress(l2+1,:);
              fLoc = (B')*sz'*vol;
              s2 = 1;
%             end
          case 12 % Hexahedra
              [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
              B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
              B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
              [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
                 state.conv.stress(l2+1:l2+obj.GaussPts.nNode,:), ...
                 state.curr.strain(l2+1:l2+obj.GaussPts.nNode,:), ...
                 dt,state.conv.status(l2+1:l2+obj.GaussPts.nNode,:), el, state.t);
              state.curr.status(l2+1:l2+obj.GaussPts.nNode,:) = status;
              state.curr.stress((l2+1):(l2+obj.GaussPts.nNode),:) = sigma;
              Ks = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
              Ks = Ks.*reshape(dJWeighed,1,1,[]);
              KLoc = sum(Ks,3);
              clear Ks;
              s1 = obj.nEntryKLoc(2);
              %sz = sigma - state.iniStress(l2+1:l2+obj.GaussPts.nNode,:);
              %sz = reshape(sz',6,1,obj.GaussPts.nNode);
              %fTmp = pagemtimes(B,'ctranspose',sz,'none');
              %fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
              %fLoc = sum(fTmp,3);
              s2 = obj.GaussPts.nNode;
%               D = obj.preP.getStiffMatrix(el,state.stress(l2+1:l2+obj.GaussPts.nNode,3)+state.iniStress(l2+1:l2+obj.GaussPts.nNode,3));
%               Js = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
%               Js = Js.*reshape(dJWeighed,1,1,[]);
%               JLoc = sum(Js,3);
%               clear Ks;
%               s1 = obj.nEntryKLoc(2);
%               %
% %               if obj.flCompRHS
%                 sz = state.stress(l2+1:l2+obj.GaussPts.nNode,:);
%                 sz = reshape(sz',6,1,obj.GaussPts.nNode);
%                 fTmp = pagemtimes(B,'ctranspose',sz,'none');
%                 fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
%                 fLoc = sum(fTmp,3);
%                 s2 = obj.GaussPts.nNode;
% %               end
        end
        %
        dof = getDoFID(obj.mesh,el);
        [jjLoc,iiLoc] = meshgrid(dof,dof);
        iiVec(l1+1:l1+s1) = iiLoc(:);
        jjVec(l1+1:l1+s1) = jjLoc(:);
        KVec(l1+1:l1+s1) = KLoc(:);
        % Accumulate the residual contributions
          %obj.rhs(dof) = obj.rhs(dof) + fLoc;
%         end
        l1 = l1 + s1;
        l2 = l2 + s2;
      end
      % Populate the Jacobian
      obj.KPoro = sparse(iiVec,jjVec,KVec,obj.mesh.nNodes*obj.mesh.nDim, ...
                     obj.mesh.nNodes*obj.mesh.nDim);
     end

     


    function computeBiotMat(obj,dt)
      %loop search inside Jacobian block matrix
      nBlock = size(obj.blockJ,1);
      k = 1;
      % Loop trough elements of blockJ matrix to find pairs PORO-FLOW
      for i=1:nBlock
          for j=k:nBlock
              if any(strcmp(obj.blockJ(i,j).physics,'Poro')) && any(strcmp(obj.blockJ(i,j).physics,'Flow')) && obj.blockJ(i,j).coupling 
                  subReg = obj.dofm.subDomains(obj.blockJ(i,j).subID).regions; %region for block's subdomain
                  subCells = find(ismember(obj.mesh.cellTag,subReg)); %id of cells in subdomain
                  nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]); 
                  % Compute the stiffness (H) and mass (P) matrices for the flow problem by FEM
                  iivec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
                  jjvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
                  Qvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
                  %
                  if nSubCellsByType(2) > 0
                    N1 = obj.elements.hexa.getBasisFinGPoints();
                  end
                  l1 = 0;
                  for el=subCells'
                    % Get the right material stiffness for each element
                   biot = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getBiotCoefficient();
                   switch obj.mesh.cellVTKType(el)
                       case 10 %Tetrahedrons, direct integration 
                           vol = findVolume(obj.elements.tetra,el);
                           der = getDerBasisF(obj.elements.tetra,el);
                           Qloc = biot*0.25*repelem(der(:),1,4)*vol;
                           s1 = obj.elements.nNodesElem(1).^2*obj.mesh.nDim;
                       case 12 %Hexahedrons, Gauss integration
                           [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                           nG = obj.GaussPts.nNode;
                           iN = zeros(6,8,nG); %matrix product i*N
                           B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
                           B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                           iN(1,:,:) = reshape(N1',1,8,nG);
                           iN(2:3,:,:) = repmat(iN(1,:,:),2,1,1);
                           Qs = biot*pagemtimes(B,'ctranspose',iN,'none');
                           Qs = Qs.*reshape(dJWeighed,1,1,[]);
                           Qloc = sum(Qs,3);
                           clear Qs;
                           s1 = obj.elements.nNodesElem(2)^2*obj.mesh.nDim;   
                   end
                   %
                   %assembly Coupling Matrix
                   if strcmp(obj.blockRhs(i).physic,'Poro')
                       dofrow = obj.dofm.getLocDoF('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                       dofcol = obj.dofm.getLocDoF('Flow',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el))); 
                   elseif strcmp(obj.blockRhs.physic,'Flow')
                       dofrow = obj.dofm.getLocDoF('Flow',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                       dofcol = obj.dofm.getLocDoF('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el))); 
                   end
                   [jjloc,iiloc] = meshgrid(dofcol,dofrow);
                   iivec(l1+1:l1+s1) = iiloc(:);
                   jjvec(l1+1:l1+s1) = jjloc(:);
                   Qvec(l1+1:l1+s1) = Qloc(:); 
                   l1 = l1+s1;
                  end
                  % Populate the Jacobian
                  Qtmp = sparse(iivec,jjvec,Qvec,obj.dofm.numDof(i),obj.dofm.numDof(j));
                  if strcmp(obj.blockJ(i,j).physics(1),'Poro')
                      obj.blockJ(i,j).block = -obj.simParams.theta*Qtmp;
                      obj.blockJ(j,i).block = Qtmp'/dt;
                  else
                      obj.blockJ(i,j).block = Qtmp'/dt;
                      obj.blockJ(j,i).block = -obj.simParams.theta*Qtmp;
                  end
              end
          end
          k = k+1; %look for upper triangular only;
      end   
    end

    function computeBiotMat_FEM_FEM_Test(obj,dt)
      %loop search inside Jacobian block matrix
      nSubPh = size(obj.blockJ,1);
      nSub = length(obj.dofm.subDomains);
      % Loop trough elements of blockJ matrix to find pairs PORO-FLOW
      for i = 1:nSub
          if obj.dofm.subDomains(i).coupling 
              subReg = obj.dofm.subDomains(i).regions; %region for block's subdomain
              subCells = find(ismember(obj.mesh.cellTag,subReg)); %id of cells in subdomain
              nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]); 
              % Compute the stiffness (H) and mass (P) matrices for the flow problem by FEM
              iivec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
              jjvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
              phVeci = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
              phVecj = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
              Qvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
              %
              l1 = 0;
              if nSubCellsByType(2) > 0
                N1 = obj.elements.hexa.getBasisFinGPoints();
              end              
              for el=subCells'
                % Get the right material stiffness for each element
               biot = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getBiotCoefficient();
               switch obj.mesh.cellVTKType(el)
                   case 10 %Tetrahedrons, direct integration 
                       vol = findVolume(obj.elements.tetra,el);
                       der = getDerBasisF(obj.elements.tetra,el);
                       Qloc = biot*0.25*repelem(der(:),1,4)*vol;
                       s1 = obj.elements.nNodesElem(1).^2*obj.mesh.nDim;
                   case 12 %Hexahedrons, Gauss integration
                       [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                       nG = obj.GaussPts.nNode;
                       iN = zeros(6,8,nG); %matrix product i*N
                       B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
                       B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                       iN(1,:,:) = reshape(N1',1,8,nG);
                       iN(2:3,:,:) = repmat(iN(1,:,:),2,1,1);
                       Qs = biot*pagemtimes(B,'ctranspose',iN,'none');
                       Qs = Qs.*reshape(dJWeighed,1,1,[]);
                       Qloc = sum(Qs,3);
                       clear Qs;
                       s1 = obj.elements.nNodesElem(2)^2*obj.mesh.nDim;   
               end
               %
               %assembly Coupling Matrix
               dofrow = obj.dofm.getLocDoF('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
               phrow = obj.dofm.getSubPhysic('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
               dofcol = obj.dofm.getLocDoF('Flow',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
               phcol = obj.dofm.getSubPhysic('Flow',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));

               [jjloc,iiloc] = meshgrid(dofcol,dofrow);
               [phLocj,phLoci] = meshgrid(phcol,phrow);
               iivec(l1+1:l1+s1) = iiloc(:);
               jjvec(l1+1:l1+s1) = jjloc(:);
               phVeci(l1+1:l1+s1) = phLoci(:);
               phVecj(l1+1:l1+s1) = phLocj(:);
               Qvec(l1+1:l1+s1) = Qloc(:); 
               l1 = l1+s1;
              end
              k = 1;
              % Populate the Jacobian
              for ii = 1:nSubPh
                  for jj = 1:nSubPh
                      if strcmp(obj.dofm.subPhysics(ii),'Poro') && any(strcmp(obj.blockJ(ii,jj).physics,'Poro')) && any(strcmp(obj.blockJ(ii,jj).physics,'Flow')) 
                         %find pairs of row,col indices belonging to ii,jj
                        tmp = intersect(find(phVeci == ii),find(phVecj == jj));
                        Qtmp = sparse(iivec(tmp),jjvec(tmp),Qvec(tmp),obj.dofm.numDof(ii),obj.dofm.numDof(jj));
                        obj.blockJ(ii,jj).block = obj.blockJ(ii,jj).block - obj.simParams.theta*Qtmp;
                        obj.blockJ(jj,ii).block = obj.blockJ(jj,ii).block + Qtmp'/dt;
                      end
                  end
              end
          end
      end
    end

    function computeBiotMat_FEM_FV_Test(obj,dt)
      %loop search inside Jacobian block matrix
      nSubPh = size(obj.blockJ,1);
      nSub = length(obj.dofm.subDomains);
      % Loop trough elements of blockJ matrix to find pairs PORO-FLOW
      for i = 1:nSub
          if obj.dofm.subDomains(i).coupling 
              subReg = obj.dofm.subDomains(i).regions; %region for block's subdomain
              subCells = find(ismember(obj.mesh.cellTag,subReg)); %id of cells in subdomain
              nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]); 
              iivec = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
              jjvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
              phVeci = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
              phVecj = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
              Qvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
              %
              if nSubCellsByType(2) > 0
                  N1 = obj.elements.hexa.getBasisFinGPoints();
              end
              l1 = 0;
              for el=subCells'
                % Get the right material stiffness for each element
               biot = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getBiotCoefficient();
               switch obj.mesh.cellVTKType(el)
                   case 10 %Tetrahedrons, direct integration 
                       error('TPFA-FV not implemented for Tetrahedrons')
                       vol = findVolume(obj.elements.tetra,el);
                       der = getDerBasisF(obj.elements.tetra,el);
                       Qloc = biot*0.25*repelem(der(:),1,4)*vol;
                       s1 = obj.elements.nNodesElem(1).^2*obj.mesh.nDim;
                   case 12 %Hexahedrons, Gauss integration
                       [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                       nG = obj.GaussPts.nNode;
                       %iN = zeros(6,8,nG); %matrix product i*N
                       B = zeros(6,8*obj.mesh.nDim,nG);
                       B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                       tmp = [1;1;1;0;0;0];
                       kron = repmat(tmp,1,1,8);
                       % iN(1,:,:) = reshape(N1',1,8,nG);
                       % iN(2:3,:,:) = repmat(iN(1,:,:),2,1,1);
                       Qs = biot*pagemtimes(B,'ctranspose',kron,'none');
                       Qs = Qs.*reshape(dJWeighed,1,1,[]);
                       Qloc = sum(Qs,3);
                       clear Qs;
                       s1 = obj.elements.nNodesElem(2)*obj.mesh.nDim;   
               end
               %
               %assembly Coupling Matrix
               dofrow = obj.dofm.getLocDoF('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
               phrow = obj.dofm.getSubPhysic('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
               dofcol = obj.dofm.getLocDoF('Flow',el);
               phcol = obj.dofm.getSubPhysic('Flow',el);

               [jjloc,iiloc] = meshgrid(dofcol,dofrow);
               [phLocj,phLoci] = meshgrid(phcol,phrow);
               iivec(l1+1:l1+s1) = iiloc(:);
               jjvec(l1+1:l1+s1) = jjloc(:);
               phVeci(l1+1:l1+s1) = phLoci(:);
               phVecj(l1+1:l1+s1) = phLocj(:);
               Qvec(l1+1:l1+s1) = Qloc(:); 
               l1 = l1+s1;
              end
              k = 1;
              % Populate the Jacobian
              for ii = 1:nSubPh
                  for jj = 1:nSubPh
                      if strcmp(obj.dofm.subPhysics(ii),'Poro') && any(strcmp(obj.blockJ(ii,jj).physics,'Poro')) && any(strcmp(obj.blockJ(ii,jj).physics,'Flow')) 
                         %find pairs of row,col indices belonging to ii,jj
                        tmp = intersect(find(phVeci == ii),find(phVecj == jj));
                        Qtmp = sparse(iivec(tmp),jjvec(tmp),Qvec(tmp),obj.dofm.numDof(ii),obj.dofm.numDof(jj));
                        obj.blockJ(ii,jj).block = obj.blockJ(ii,jj).block - obj.simParams.theta*Qtmp;
                        obj.blockJ(jj,ii).block = obj.blockJ(jj,ii).block + Qtmp'/dt;
                      end
                  end
              end
          end
      end
    end

    function computeBiotRhs(obj,stateTmp,statek)
      %loop trough rhs blocks to find Biot coupling residual contribution 
      nRhs = length(obj.blockRhs);
      theta = obj.simParams.theta;
      for i = 1:nRhs
          if strcmp(obj.blockRhs(i).physic,'Poro')
              for j = 1:nRhs
                  if any(strcmp(obj.blockJ(i,j).physics,'Poro')) && any(strcmp(obj.blockJ(i,j).physics,'Flow'))
                    dofs = obj.dofm.getLocDoF(j);  
                    obj.blockRhs(i).block = obj.blockRhs(i).block + ...
                        obj.blockJ(i,j).block*stateTmp.pressure(dofs) + ...
                        (1/theta-1)*(obj.blockJ(i,j).block*statek.pressure(dofs));
                  end
              end
          elseif strcmp(obj.blockRhs(i).physic,'Flow') 
             for j = 1:nRhs
              if any(strcmp(obj.blockJ(i,j).physics,'Poro')) && any(strcmp(obj.blockJ(i,j).physics,'Flow'))
                dofs = obj.dofm.getLocDoF(j);  
                obj.blockRhs(i).block = obj.blockRhs(i).block+...
                obj.blockJ(i,j).block*(stateTmp.dispCurr(dofs)-stateTmp.dispConv(dofs));
              end
             end

          end
      end
    end
    
    
    function computeCoupleMat(obj)
       %initializing index vectors for assembly 
       iivec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*obj.elements.nCellsByType,1);
       jjvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*obj.elements.nCellsByType,1);
       Qvec =  zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*obj.elements.nCellsByType,1);
       
       l1=0;
       
             
       if obj.elements.nCellsByType(2) > 0 %at least one Hexahedron
           N1 = getBasisFinGPoints(obj.elements.hexa); %compute Basis Function matrix for further calculations
       end
       
       for el=1:obj.mesh.nCells
           biot = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getBiotCoefficient();
           switch obj.mesh.cellVTKType(el)
               case 10 %Tetrahedrons, direct integration 
                   vol = findVolume(obj.elements.tetra,el);
                   der = getDerBasisF(obj.elements.tetra,el);
                   Qloc = biot*0.25*repelem(der(:),1,4)*vol;
                   s1 = obj.elements.nNodesElem(1).^2*obj.mesh.nDim;
               case 12 %Hexahedrons, Gauss integration
                   [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                   nG = obj.GaussPts.nNode;
                   iN = zeros(6,8,nG); %matrix product i*N
                   B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
                   B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                   iN(1,:,:) = reshape(N1',1,8,nG);
                   iN(2:3,:,:) = repmat(iN(1,:,:),2,1,1);
                   Qs = biot*pagemtimes(B,'ctranspose',iN,'none');
                   Qs = Qs.*reshape(dJWeighed,1,1,[]);
                   Qloc = sum(Qs,3);
                   clear Qs;
                   s1 = obj.elements.nNodesElem(2)^2*obj.mesh.nDim;   
           end
       
       %assembly Coupling Matrix
       dofporo = getDoFID(obj.mesh,el);
       dofflow = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
       [jjloc,iiloc] = meshgrid(dofflow,dofporo);
       iivec(l1+1:l1+s1) = iiloc(:);
       jjvec(l1+1:l1+s1) = jjloc(:);
       Qvec(l1+1:l1+s1) = Qloc(:);
       
       l1 = l1+s1;
       end  
       obj.Q = sparse(iivec,jjvec,Qvec,obj.mesh.nNodes*obj.mesh.nDim, ...
                     obj.mesh.nNodes);
    end

     function computeBiotSyst(obj,dt,statek,stateTmp)
        %This methods performs matrix assembly and residual computation
        %of all terms dealing with coupled poromechanics 
        %Flow matrices are computed before time loop
        %compute coupled Poromechanics blocks
        computePoroCoupled_Test(obj,stateTmp,dt)
        if isFEMBased(obj.model,'Flow')
            computeBiotMat_FEM_FEM_Test(obj,dt);
        elseif isFVTPFABased(obj.model,'Flow')
            computeBiotMat_FEM_FV_Test(obj,dt);
        end
        computePoroRhs(obj,stateTmp);
        computeBiotRhs(obj,stateTmp,statek);
     end

    function buildGlobalJacobianAndRhs(obj)
        %get Jacobian and rhs vector from block struct
        nBlock = length(obj.blockRhs); 
        cellJ = struct2cell(obj.blockJ);
        cellJ = reshape(cellJ(1,:,:),nBlock,nBlock);
        obj.J = cell2mat(cellJ);
        cellRhs = struct2cell(obj.blockRhs);
        cellRhs = reshape(cellRhs(1,:),[],1);
        obj.rhs = cell2mat(cellRhs);
    end

    function resetJacobianAndRhs(obj)
        %reset Jacobian and Rhs to zero
        dim = length(obj.dofm.subList);
        for i = 1:dim
            obj.blockRhs(i).block = zeros(obj.dofm.numDof(i),1);
            for j=1:dim
                obj.blockJ(i,j).block = sparse(obj.dofm.numDof(i),obj.dofm.numDof(j));
                % flow matrices H and P stay the same
            end
        end
    end  


     
    function computeCoupleSyst(obj,theta,dt,statek,stateTmp)
        %Not necessary since they are already computed before the NR loop
        %         obj.computeFlowMat(); 
        %         obj.computeFlowRHSGravContribute();

        %initializing rhs vector
        nDofPoro = obj.mesh.nDim*obj.mesh.nNodes;
        obj.rhs = zeros(nDofPoro + obj.mesh.nNodes,1);
        

        obj.computeCoupleMat(); %computing coupling matrix
        %obj.rhsFlow = obj.RHSGravTerm; %initializing flow RHS, not needed
        %here, directly added after
        obj.computePoroCoupled(stateTmp, dt); %computing Poromechanics stiffness and initializing Poromechanics rhs with internal forces
        %obj.computeCoupleFlowSystMat(theta,dt);
        %computing rhs contributes (neumann and volume forces still
        %missing! they are added with applyBCandForces function)
        
        %temporarily poromechanics internal forces are computed with
        %displacement. Update with stress computation!
        obj.rhs(1:nDofPoro) = theta*obj.KPoro*stateTmp.dispCurr - theta*obj.Q*stateTmp.pressure + ...
            (1-theta)*(obj.KPoro*stateTmp.dispConv-obj.Q*statek.pressure);
        
        obj.rhs(nDofPoro+1:end) = ((obj.Q)'/dt)*stateTmp.dispCurr+(theta*obj.H+obj.P/dt)*stateTmp.pressure +...
            (1-theta)*obj.H*statek.pressure -(1/dt)*((obj.Q)'*stateTmp.dispConv+obj.P*statek.pressure);
        
       %construct Jacobian and rhs for Coupled System
        obj.J = [theta*obj.KPoro, -theta*obj.Q; (obj.Q)'/dt, theta*obj.H+obj.P/dt];
    end

    
    
    
    
    
    
    
    
%     function applyBC(obj)
%       l = length(obj.BCName);
%       if l == 0
%         error('Warning: No boundary conditions will be applied.');
%       end
%       for i=1:l
%         cond = getBC(obj.bound,obj.BCName(i));
%         if isequal(cond.boundType,'neu')  % Apply Neumann conditions,if any
%           obj.rhs(cond.boundDof) = obj.rhs(cond.boundDof) + cond.boundVal;
%         end
%         %
%         if isequal(cond.boundType,'dir')  % Apply Dirichlet conditions
%           maxVal = max(abs(obj.K), [], 'all');
%           obj.rhs(cond.boundDof) = cond.boundVal*(maxVal*10^10);
%           obj.K(obj.mesh.nNodes*obj.mesh.nDim*(cond.boundDof-1)+cond.boundDof) = maxVal*10^10;
%         end
%       end
%     end
  end
  
  methods(Access = private)
    function setDiscretizer(obj,symmod,params,dofManager,grid,mat,data)
      obj.model = symmod;
      obj.simParams = params;
      obj.dofm = dofManager;
      obj.mesh = grid.topology;
      obj.elements = grid.cells;
      obj.faces = grid.faces;
      obj.material = mat;
%       obj.probType = pType;
%       obj.bound = bc;
%       obj.BCName = BCName;
%       obj.state = stat;
      if ~isempty(data)
        obj.GaussPts = data{1};
      end
      %
      obj.nEntryKLoc = (obj.mesh.nDim^2)*(obj.elements.nNodesElem).^2;
      %
      initializeBlockJacobianAndRhs(obj)
      %%% richards model -  will be checked in the future %%%%%%%%%%%%%%%%
      if obj.model.isFVTPFABased('Flow')
        obj.computeTrans;
        %get cells with active flow model
        
        % internal faces inside each subdomain
        nSub = length(obj.dofm.subDomains);
        IntFaces = zeros(obj.faces.nFaces,nSub);
        flowCells = [];
        for i = 1:nSub
            if any(strcmp(obj.dofm.subDomains(i).physics,"Flow"))
            flowCells = [flowCells; find(obj.dofm.subCells(:,i))];
            end
        end

        for i = 1:nSub
            intcells = find(obj.dofm.subCells(:,i));
            tmp1 = any(ismember(obj.faces.faceNeighbors, intcells),2);
            tmp2 = all(ismember(obj.faces.faceNeighbors, flowCells),2);
            IntFaces(:,i) = all([tmp1 tmp2],2);
        end
        obj.isIntFaces = logical(IntFaces);
        if obj.model.isVariabSatFlow()
          obj.upElem = zeros(nnz(obj.isIntFaces),1);
        end
      end
      %
      if obj.model.isFlow()
        computeFlowRHSGravTerm_Test(obj);
      end
%       if isSinglePhaseFlow(obj.model)
%         obj.fOld = zeros(obj.mesh.nNodes,1);
%         obj.fNew = zeros(obj.mesh.nNodes,1);
%       end
%       if strcmpi(obj.probType,'lin')
%         obj.flCompRHS = false;
%       elseif strcmpi(obj.probType,'nonlin')
%         obj.flCompRHS = true;
%       end
    end
    

    function [JNewt] = computeNewtPartOfJacobian(obj,dt,statek,stateTmp,pkpt,dSwkpt,dlwkpt)
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      zVec = obj.elements.cellCentroid(:,3);
      zNeigh = zVec(neigh);
      gamma = obj.material.getFluid().getFluidSpecWeight();
      tmpVec1 = (dlwkpt.*obj.trans(obj.isIntFaces)).*(pkpt(neigh(:,1)) - pkpt(neigh(:,2)) + gamma*(zNeigh(:,1) - zNeigh(:,2)));
      %
      poroMat = zeros(obj.mesh.nCellTag,1);
      alphaMat = zeros(obj.mesh.nCellTag,1);
      beta = obj.material.getFluid().getFluidCompressibility();
      for m = 1:obj.mesh.nCellTag
        poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
        alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
      end
      % (alpha+poro*beta)
      tmpVec2 = alphaMat(obj.mesh.cellTag) + beta*poroMat(obj.mesh.cellTag);
      tmpVec2 = 1/dt*((tmpVec2.*dSwkpt).*(stateTmp.pressure - statek.pressure)).*obj.elements.vol;
%       tmpVec = lw.*tmpVec;
%       sumDiagTrans = accumarray([neigh1; neigh2], ...
%         repmat(tmpVec,[2,1]),[obj.mesh.nCells,1]);
      JNewt = sparse([neigh(:,1); neigh(:,2); (1:obj.mesh.nCells)'], ...
                     [repmat(obj.upElem,[2,1]); (1:obj.mesh.nCells)'], ...
                     [tmpVec1; -tmpVec1; tmpVec2],obj.mesh.nCells,obj.mesh.nCells);
      JNewt = obj.simParams.theta*JNewt;
    end
    
    function computeTrans(obj)   % Inspired by MRST
      % Compute first the vector connecting each cell centroid to the
      % half-face
      r = [1, 1, 1, 2, 2, 2, 3, 3, 3];
      c = [1, 2, 3, 1, 2, 3, 1, 2, 3];
      hf2Cell = repelem((1:obj.mesh.nCells)',diff(obj.faces.mapF2E));
      L = obj.faces.faceCentroid(obj.faces.faces2Elements(:,1),:) - obj.elements.cellCentroid(hf2Cell,:);
      sgn = 2*(hf2Cell == obj.faces.faceNeighbors(obj.faces.faces2Elements(:,1))) - 1;
      N = bsxfun(@times,sgn,obj.faces.faceNormal(obj.faces.faces2Elements(:,1),:));
      KMat = zeros(obj.mesh.nCellTag,9);
      for i=1:obj.mesh.nCellTag
        KMat(i,:) = obj.material.getMaterial(i).PorousRock.getPermVector();
      end
      hT = zeros(length(hf2Cell),1);
      for k=1:length(r)
        hT = hT + L(:,r(k)) .* KMat(obj.mesh.cellTag(hf2Cell),k) .* N(:,c(k));
      end
      hT = hT./sum(L.*L,2);
%       mu = obj.material.getMaterial(obj.mesh.nCellTag+1).getDynViscosity();
%       hT = hT/mu;
      %
      obj.trans = 1 ./ accumarray(obj.faces.faces2Elements(:,1),1 ./ hT,[obj.faces.nFaces,1]);
    end
    
    function computeFlowRHSGravTerm(obj)
      % Compute the gravity contribution
      % Get the fluid specific weight and viscosity
      gamma = obj.material.getFluid().getFluidSpecWeight();
      if gamma > 0
%         mu = obj.material.getMaterial(obj.mesh.nCellTag+1).getDynViscosity();
        if isFEMBased(obj.model,'Flow')
          nr = obj.mesh.nNodes;
          obj.RHSGravTerm = zeros(nr,1);
          for el=1:obj.mesh.nCells
            % Get the material permeability
            permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
%             permMat = permMat/mu;
            switch obj.mesh.cellVTKType(el)
              case 10 % Tetrahedra
                N = obj.elements.tetra.getDerBasisF(el);
  %               volSign = getVolumeSign(obj.elements,el);
%                 vol = getVolume(obj.elements,el);
  %               fLoc = (N'*permMat(:,3))*volSign*vol*gamma;
                fLoc = (N'*permMat(:,3))*obj.elements.vol(el)*gamma;
              case 12 % Hexa
                [N,dJWeighed] = obj.elements.hexa.getDerBasisFAndDet(el,1);
                fs = pagemtimes(N,'ctranspose',permMat(:,3),'none');
                fs = fs.*reshape(dJWeighed,1,1,[]);
                fLoc = sum(fs,3)*gamma;
            end
            %
            dof = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
            obj.RHSGravTerm(dof) = obj.RHSGravTerm(dof) + fLoc;
          end
        elseif isFVTPFABased(obj.model,'Flow')
          neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
%           nr = nnz(obj.isIntFaces);
          zVec = obj.elements.cellCentroid(:,3);
          zNeigh = zVec(neigh);
          obj.RHSGravTerm = gamma*obj.trans(obj.isIntFaces).*(zNeigh(:,1) - zNeigh(:,2));
        end
      end
      % Initializing fOld
      % SOURCE/SINK CONTRIBUTION IS MISSING AT THE MOMENT
%       obj.fOld = obj.fConst;
    end


    function computeFlowRHSGravTerm_Test(obj)
      % Compute the gravity contribution
      % Get the fluid specific weight and viscosity
      gamma = obj.material.getFluid().getFluidSpecWeight();
      if gamma > 0
          nBlock = length(obj.blockRhs);
%         mu = obj.material.getMaterial(obj.mesh.nCellTag+1).getDynViscosity();
                for i = 1:nBlock
                    if strcmp(obj.blockRhs(i).physic,'Flow')
                      obj.blockRhs(i).rhsGrav = obj.blockRhs(i).block;
                      subReg = obj.dofm.subDomains(obj.dofm.subList(i)).regions; %region for block's subdomain
                      subCells = find(ismember(obj.mesh.cellTag,subReg)); %id of cells in subdomain
                      %nSubCells = length(subCells); %number of cells in
                      %subdomain
                      if isFEMBased(obj.model,'Flow')
                      l1 = 0;
                      for el = subCells'
                        % Get the material permeability
                        permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
            %             permMat = permMat/mu;
                        switch obj.mesh.cellVTKType(el)
                          case 10 % Tetrahedra
                            N = obj.elements.tetra.getDerBasisF(el);
              %               volSign = getVolumeSign(obj.elements,el);
            %                 vol = getVolume(obj.elements,el);
              %               fLoc = (N'*permMat(:,3))*volSign*vol*gamma;
                            fLoc = (N'*permMat(:,3))*obj.elements.vol(el)*gamma;
                            s1 = obj.elements.nNodesElem(1);
                          case 12 % Hexa
                            [N,dJWeighed] = obj.elements.hexa.getDerBasisFAndDet(el,1);
                            fs = pagemtimes(N,'ctranspose',permMat(:,3),'none');
                            fs = fs.*reshape(dJWeighed,1,1,[]);
                            fLoc = sum(fs,3)*gamma;
                            s1 = obj.elements.nNodesElem(2);
                        end
                        %
                        dof = obj.dofm.getLocDoF('Flow',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                        ph = obj.dofm.getSubPhysic('Flow',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                        dofVec(l1+1:l1+s1) = dof;
                        phVec(l1+1:l1+s1) = ph;
                        fVec(l1+1:l1+s1) = fLoc;
                        l1 = l1 + s1;
                      end
                      for ii = 1:nBlock
                          tmp = phVec == ii;
                          obj.blockRhs(ii).rhsGrav(dofVec(tmp))  = obj.blockRhs(ii).rhsGrav(dofVec(tmp)) + fVec(tmp);
                      end
                      elseif isFVTPFABased(obj.model,'Flow')
                          subID = obj.dofm.subList(i);
                          neigh = obj.faces.faceNeighbors(obj.isIntFaces(:,subID),:);
                %           nr = nnz(obj.isIntFaces);
                          zVec = obj.elements.cellCentroid(:,3);
                          zNeigh = zVec(neigh);
                          obj.blockRhs(i).rhsGrav = gamma*obj.trans(obj.isIntFaces(:,subID)).*(zNeigh(:,1) - zNeigh(:,2));
                      end
                    end
                end
      end
      % Initializing fOld
      % SOURCE/SINK CONTRIBUTION IS MISSING AT THE MOMENT
%       obj.fOld = obj.fConst;
    end
    
    function gTerm = finalizeRHSGravTerm(obj,idblock,lw)
      % varargin{1} -> lw
      subID = obj.dofm.subList(idblock);
      neigh = obj.faces.faceNeighbors(obj.isIntFaces(:,subID),:);
      dof = obj.dofm.getLocDoF('Flow',neigh(:));
      if obj.model.isSinglePhaseFlow()
        gTerm = accumarray(dof,[obj.blockRhs(idblock).rhsGrav; -obj.blockRhs(idblock).rhsGrav],[obj.dofm.numDof(idblock),1]);
      elseif obj.model.isVariabSatFlow()
        %tmpVec = lw.*obj.RHSGravTerm;
      gTerm = accumarray(dof,[lw.*obj.blockRhs.rhsGrav; ...
        -lw.*obj.blockRhs.rhsGrav],[obj.dofm.numDof(idblock),1]);
%       end
      end
    end
    
    function initializeBlockJacobianAndRhs(obj)
        %block dimension of the Jacobian matrix
        dim = length(obj.dofm.subList);
        obj.blockJ = repmat(struct('block',[],'subID',[],'physics',[],'coupling',[]),dim);
        obj.blockRhs = repmat(struct('block',[],'subID',[],'physic',[]),dim,1);
        for i = 1:dim
            subID = obj.dofm.subList(i);
            obj.blockRhs(i).subID = subID;
            obj.blockRhs(i).physic = obj.dofm.subPhysics(i);
            obj.blockRhs(i).block = zeros(obj.dofm.numDof(i),1);
            %obj.blockRhs(i).coupling = obj.dofm.subDomains(subID).coupling;
            for j=1:dim
                subID = obj.dofm.subList(i);
                obj.blockJ(i,j).subID = subID;
                obj.blockJ(i,j).physics = obj.dofm.subPhysics([i j]);
                obj.blockJ(i,j).block = sparse(obj.dofm.numDof(i),obj.dofm.numDof(j));
                if all(strcmp(obj.blockJ(i,j).physics,'Flow'))
                    %Initialize flow elementary matrix
                    obj.blockJ(i,j).H = sparse(obj.dofm.numDof(i),obj.dofm.numDof(j));
                    obj.blockJ(i,j).P = sparse(obj.dofm.numDof(i),obj.dofm.numDof(j));
                end
            end
        end  
    end

              
    function [Swkpt,dSwkpt,lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt)
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      gamma = obj.material.getFluid().getFluidSpecWeight();
      if gamma > 0
        zVec = obj.elements.cellCentroid(:,3);
        zNeigh = zVec(neigh);
        lElemIsUp = pkpt(neigh(:,1)) - pkpt(neigh(:,2)) + gamma*(zNeigh(:,1) - zNeigh(:,2)) >= 0;
      else
        lElemIsUp = pkpt(neigh(:,1)) >= pkpt(neigh(:,2));
      end
      obj.upElem(lElemIsUp) = neigh(lElemIsUp,1);
      obj.upElem(~lElemIsUp) = neigh(~lElemIsUp,2);
      [Swkpt,dSwkpt] = obj.material.computeSwAnddSw(obj.mesh,pkpt);
%       dSwkpt = zeros(length(dSwkpt),1);
      dSwkpt = - dSwkpt;
      [lwkpt,dlwkpt] = obj.material.computeLwAnddLw(obj.mesh,obj.upElem,pkpt);
      dlwkpt = - dlwkpt;
    end
  end
end