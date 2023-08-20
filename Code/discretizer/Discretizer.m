classdef Discretizer < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    H
    P
    J
    rhs
    KPoro
    Q
  end
  
  properties (Access = public)
    model
    simParams
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
    function obj = Discretizer(symmod,simParams,grid,mat,varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      obj.setDiscretizer(symmod,simParams,grid,mat,varargin);
    end
    
    function trans = getFaceTransmissibilities(obj,faceID)
      trans = obj.trans(faceID);
    end
    
    function computeSPFMatrices(obj)
      if obj.model.isFEMBased('Flow')
        computeFlowStiffAndCapMatFEM(obj);
      elseif obj.model.isFVTPFABased('Flow')
        mu = obj.material.getMaterial(obj.mesh.nCellTag+1).getDynViscosity();
        computeFlowStiffMatFV(obj,1/mu);
        computeFlowCapMatFV(obj);
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
      beta = obj.material.getMaterial(obj.mesh.nCellTag+1).getFluidCompressibility();
      if obj.elements.nCellsByType(2) > 0
        N1 = obj.elements.hexa.getBasisFinGPoints();
      end
      % Get the fluid dynamic viscosity
      mu = obj.material.getMaterial(obj.mesh.nCellTag+1).getDynViscosity();
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
    
    function computeFlowStiffMatFV(obj,lw)
      % Inspired by MRST
      neigh1 = obj.faces.faceNeighbors(obj.isIntFaces,1);
      neigh2 = obj.faces.faceNeighbors(obj.isIntFaces,2);
      tmpVec = lw.*obj.trans(obj.isIntFaces);
%       tmpVec = lw.*tmpVec;
      sumDiagTrans = accumarray([neigh1; neigh2], ...
        repmat(tmpVec,[2,1]),[obj.mesh.nCells,1]);
      obj.H = sparse([neigh1; neigh2; (1:obj.mesh.nCells)'], ...
                     [neigh2; neigh1; (1:obj.mesh.nCells)'], ...
                     [-tmpVec; -tmpVec; ...
                      sumDiagTrans],obj.mesh.nCells,obj.mesh.nCells);
    end
    
    function computeFlowCapMatFV(obj,varargin)
      % if isVariabSatFlow
      % varargin{1} -> Sw
      % varargin{2} -> dSw
      poroMat = zeros(obj.mesh.nCellTag,1);
      alphaMat = zeros(obj.mesh.nCellTag,1);
      beta = obj.material.getMaterial(obj.mesh.nCellTag+1).getFluidCompressibility();
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

    function computeFlowRHS(obj,statek,stateTmp,dt,lw)
      % Compute the residual of the flow problem
%       obj.rhs = (obj.simParams.theta*obj.H + obj.P/dt)*stateTmp.pressure ...
%         - (obj.P/dt - (1 - obj.simParams.theta)*obj.H)*statek.pressure;
      obj.rhs = obj.simParams.theta*obj.H*stateTmp.pressure + 1/dt*obj.P*stateTmp.pressure ...
        - 1/dt*obj.P*statek.pressure + (1 - obj.simParams.theta)*obj.H*statek.pressure;
      gamma = obj.material.getMaterial(obj.mesh.nCellTag+1).getFluidSpecWeight();
      if gamma > 0
        if isFEMBased(obj.model,'Flow')
          obj.rhs = obj.rhs + obj.RHSGravTerm.*lw;
        elseif isFVTPFABased(obj.model,'Flow')
          obj.rhs = obj.rhs + finalizeRHSGravTerm(obj,lw);
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
                 dt,state.conv.status(l2+1,:));
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
                 dt,state.conv.status(l2+1:l2+obj.GaussPts.nNode,:));
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
          case 10 % Tetrahedra
           N = getDerBasisF(obj.elements.tetra,el);
            vol = findVolume(obj.elements.tetra,el);
            B = zeros(6,4*obj.mesh.nDim);
            B(obj.elements.indB(1:36,2)) = N(obj.elements.indB(1:36,1));
            [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
                 state.conv.stress(l2+1,:), ...
                 state.curr.strain(l2+1,:), ...
                 dt,state.conv.status(l2+1,:));
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
                 dt,state.conv.status(l2+1:l2+obj.GaussPts.nNode,:));
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
    
    
    function computeCoupleMat(obj)
       %initializing index vectors for assembly 
       iivec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*obj.elements.nCellsByType,1);
       jjvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*obj.elements.nCellsByType,1);
       Qvec =  zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*obj.elements.nCellsByType,1);
       
       l1=0;
       biot = 1; %TO DO: call to Biot Coefficient in PorousRock class
       
       if obj.elements.nCellsByType(2) > 0 %at least one Hexahedron
           N1 = getBasisFinGPoints(obj.elements.hexa); %compute Basis Function matrix for further calculations
       end
       
       for el=1:obj.mesh.nCells
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
    function setDiscretizer(obj,symmod,params,grid,mat,data)
      obj.model = symmod;
      obj.simParams = params;
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
      if obj.model.isFVTPFABased('Flow')
        obj.computeTrans;
        obj.isIntFaces = all(obj.faces.faceNeighbors ~= 0,2);
        if obj.model.isVariabSatFlow()
          obj.upElem = zeros(nnz(obj.isIntFaces),1);
        end
      end
      %
      if obj.model.isFlow()
        computeFlowRHSGravTerm(obj);
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
      gamma = obj.material.getMaterial(obj.mesh.nCellTag+1).getFluidSpecWeight();
      tmpVec1 = (dlwkpt.*obj.trans(obj.isIntFaces)).*(pkpt(neigh(:,1)) - pkpt(neigh(:,2)) + gamma*(zNeigh(:,1) - zNeigh(:,2)));
      %
      poroMat = zeros(obj.mesh.nCellTag,1);
      alphaMat = zeros(obj.mesh.nCellTag,1);
      beta = obj.material.getMaterial(obj.mesh.nCellTag+1).getFluidCompressibility();
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
      gamma = obj.material.getMaterial(obj.mesh.nCellTag+1).getFluidSpecWeight();
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
    
    function gTerm = finalizeRHSGravTerm(obj,lw)
      % varargin{1} -> lw
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
%       if obj.model.isSinglePhaseFlow()
%         gTerm = accumarray(neigh(:),[obj.RHSGravTerm; -obj.RHSGravTerm],[obj.mesh.nCells,1]);
%       elseif obj.model.isVariabSatFlow()
%         tmpVec = lw.*obj.RHSGravTerm;
      gTerm = accumarray(neigh(:),[lw.*obj.RHSGravTerm; ...
        -lw.*obj.RHSGravTerm],[obj.mesh.nCells,1]);
%       end
    end
    
    function [Swkpt,dSwkpt,lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt)
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      gamma = obj.material.getMaterial(obj.mesh.nCellTag+1).getFluidSpecWeight();
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