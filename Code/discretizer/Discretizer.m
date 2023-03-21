classdef Discretizer < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    H
    P
    K
    K2
    rhs
  end
  
  properties (Access = private)
    model
    mesh
    elements
    material
%     bound
%     BCName
    preP
%     state
    GaussPts
%     probType  %either linear (lin) or nonlinear (nonlin)
%     flCompRHS = false
%     nE   % nE = [#tetra, #hexa, #wed, #pyr]
%     nEntryKLoc
    fConst
  end
  
  methods (Access = public)
    function obj = Discretizer(symmod,msh,elem,mat,pre,varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      nIn = nargin;
      data = varargin;
      obj.setDiscretizer(nIn,symmod,msh,elem,mat,pre,data);
    end
    
    function computeFlowMat(obj)
      % Compute the stiffness (H) and mass (P) matrices for the flow problem by FEM
      iiVec = zeros((obj.preP.nNodesElem.^2)*obj.preP.nE,1);
      jjVec = zeros((obj.preP.nNodesElem.^2)*obj.preP.nE,1);
      HVec = zeros((obj.preP.nNodesElem.^2)*obj.preP.nE,1);
      PVec = zeros((obj.preP.nNodesElem.^2)*obj.preP.nE,1);
      % Get the fluid compressibility
      beta = obj.material.getMaterial(2*obj.preP.nMat+1).getFluidCompressibility();
      if obj.preP.nE(2) > 0
        N1 = getBasisFinGPoints(obj.elements);
      end
      % Get the fluid dynamic viscosity
      mu = obj.material.getMaterial(2*obj.preP.nMat+1).getDynViscosity();
      %
      l1 = 0;
      for el=1:obj.mesh.nCells
        % Get the rock permeability, porosity and compressibility
        permMat = obj.material.getMaterial(obj.preP.nMat+obj.mesh.cellTag(el)).getPermMatrix();
        poro = obj.material.getMaterial(obj.preP.nMat+obj.mesh.cellTag(el)).getPorosity();
        alpha = obj.material.getMaterial(obj.mesh.cellTag(el)).getRockCompressibility();
        % Compute the element matrices based on the element type
        % (tetrahedra vs. hexahedra)
        switch obj.mesh.cellVTKType(el)
          case 10 % Tetrahedra
            % Computing the H matrix contribution
            N = getDerBasisF(obj.elements,el);
            vol = getVolume(obj.elements,el);
            HLoc = N'*permMat*N*vol/mu;
            s1 = obj.preP.nNodesElem(1)^2;
            % Computing the P matrix contribution
            PLoc = ((alpha + poro*beta)*vol/20)*(ones(obj.preP.nNodesElem(1)) + eye(obj.preP.nNodesElem(1)));
          case 12 % Hexa
            [N,dJWeighed] = getDerBasisFAndDet(obj.elements,el);
            permMat = permMat/mu;
            Hs = pagemtimes(pagemtimes(N,'ctranspose',permMat,'none'),N);
            Hs = Hs.*reshape(dJWeighed,1,1,[]);
            HLoc = sum(Hs,3);
            clear Hs;
            s1 = obj.preP.nNodesElem(2)^2;
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
    
    function computeFlowSystMat(obj,theta,dt)
      % Compute matrices K and K2
      obj.K = theta*obj.H + obj.P/dt;
      obj.K2 = obj.P/dt - (1-theta)*obj.H;
    end
    
    function computeFlowRHSGravContribute(obj)
      % Compute the gravity contribution
      obj.fConst = zeros(obj.mesh.nNodes,1);
      % Get the fluid specific weight and viscosity
      gamma = obj.material.getMaterial(2*obj.preP.nMat+1).getFluidSpecWeight();
      mu = obj.material.getMaterial(2*obj.preP.nMat+1).getDynViscosity();
      if gamma > 0
        for el=1:obj.mesh.nCells
          % Get the material permeability
          permMat = obj.material.getMaterial(obj.preP.nMat+obj.mesh.cellTag(el)).getPermMatrix();
          permMat = permMat/mu;
          switch obj.mesh.cellVTKType(el)
            case 10 % Tetrahedra
              N = getDerBasisF(obj.elements,el);
%               volSign = getVolumeSign(obj.elements,el);
              vol = getVolume(obj.elements,el);
%               fLoc = (N'*permMat(:,3))*volSign*vol*gamma;
              fLoc = (N'*permMat(:,3))*vol*gamma;
            case 12 % Hexa
              [N,dJWeighed] = getDerBasisFAndDet(obj.elements,el);
              fs = pagemtimes(N,'ctranspose',permMat(:,3),'none');
              fs = fs.*reshape(dJWeighed,1,1,[]);
              fLoc = sum(fs,3)*gamma;
          end
          %
          dof = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
          obj.fConst(dof) = obj.fConst(dof) + fLoc;
        end
      end
      % Initializing fOld
      % SOURCE/SINK CONTRIBUTION IS MISSING AT THE MOMENT
%       obj.fOld = obj.fConst;
    end
    
    function computeFlowRHS(obj,statek,stateTmp)
      % Compute the residual of the flow problem
      obj.rhs = obj.fConst + obj.K*stateTmp.pressure - obj.K2*statek.pressure;
    end
    
    function computePoroSyst(obj,state)
      % Compute the Jacobian and residual of the geomechanical problem
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
%       KLocSize = obj.mesh.cellNumVerts(1)*obj.mesh.nDim;
      obj.rhs = zeros(obj.mesh.nNodes*obj.mesh.nDim,1);
      iiVec = zeros(obj.preP.nEntryKLoc*obj.preP.nE,1);
      jjVec = zeros(obj.preP.nEntryKLoc*obj.preP.nE,1);
      KVec = zeros(obj.preP.nEntryKLoc*obj.preP.nE,1);
      %
      l1 = 0;
      l2 = 0;
      for el=1:obj.mesh.nCells
        % Get the right material stiffness for each element
        switch obj.mesh.cellVTKType(el)
          case 10 % Tetrahedra
            N = getDerBasisF(obj.elements,el);
            vol = getVolume(obj.elements,el);
            B = zeros(6,4*obj.mesh.nDim);
            B(obj.preP.indB(1:36,2)) = N(obj.preP.indB(1:36,1));
            D = obj.preP.getStiffMatrix(el,state.stress(l2+1,3)+state.iniStress(l2+1,3));
            KLoc = B'*D*B*vol;
            s1 = obj.preP.nEntryKLoc(1);
            %
%             if obj.flCompRHS
              fLoc = (B')*(state.stress(l2+1,:))'*vol;
              s2 = 1;
%             end
          case 12 % Hexahedra
%             [N,dJWeighed] = getDerBasisFAndDet(obj.elements,el);
%             sh1 = 0;
%             KLoc = zeros(obj.mesh.cellNumVerts(1)*obj.mesh.nDim);
%             for i=1:obj.GaussPts.nNode
%               B = zeros(6,obj.mesh.cellNumVerts(1)*obj.mesh.nDim);
%               B(obj.i2) = N(obj.i1+sh1);
%               KLoc = KLoc + B'*D*B*dJWeighed(i);
%               sh1 = sh1 + obj.mesh.nDim*obj.mesh.cellNumVerts(1);
%             end
              [N,dJWeighed] = getDerBasisFAndDet(obj.elements,el);
              B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
              B(obj.preP.indB(:,2)) = N(obj.preP.indB(:,1));
%               KTmp = pagemtimes(B,'ctranspose',D,'none');
%       end
              D = obj.preP.getStiffMatrix(el,state.stress(l2+1:l2+obj.GaussPts.nNode,3)+state.iniStress(l2+1:l2+obj.GaussPts.nNode,3));
              Ks = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
              Ks = Ks.*reshape(dJWeighed,1,1,[]);
              KLoc = sum(Ks,3);
              clear Ks;
              s1 = obj.preP.nEntryKLoc(2);
              %
%               if obj.flCompRHS
                sz = state.stress(l2+1:l2+obj.GaussPts.nNode,:);
                sz = reshape(sz',6,1,obj.GaussPts.nNode);
                fTmp = pagemtimes(B,'ctranspose',sz,'none');
                fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
                fLoc = sum(fTmp,3);
                s2 = obj.GaussPts.nNode;
%               end
        end
        %
        dof = obj.preP.getDoFID(el);
        [jjLoc,iiLoc] = meshgrid(dof,dof);
        iiVec(l1+1:l1+s1) = iiLoc(:);
        jjVec(l1+1:l1+s1) = jjLoc(:);
        KVec(l1+1:l1+s1) = KLoc(:);
        % Accumulate the residual contributions
          obj.rhs(dof) = obj.rhs(dof) + fLoc;
%         end
        l1 = l1 + s1;
        l2 = l2 + s2;
      end
      % Populate the Jacobian
      obj.K = sparse(iiVec,jjVec,KVec,obj.mesh.nNodes*obj.mesh.nDim, ...
                     obj.mesh.nNodes*obj.mesh.nDim);
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
    function setDiscretizer(obj,nIn,symmod,msh,elem,mat,pre,data)
      obj.model = symmod;
      obj.mesh = msh;
      obj.elements = elem;
      obj.material = mat;
%       obj.probType = pType;
%       obj.bound = bc;
%       obj.BCName = BCName;
      obj.preP = pre;
%       obj.state = stat;
      if nIn > 5
        obj.GaussPts = data{1};
      end
      %
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
  end
end