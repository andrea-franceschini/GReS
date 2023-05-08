classdef Discretizer < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    H
    P
    K
    K2
    Q
    rhs
    rhsPoro
    rhsFlow
    KPoro
  end
  
  properties (Access = private)
    model
    mesh
    elements
    faces
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
    trans
  end
  
  methods (Access = public)
    function obj = Discretizer(symmod,grid,mat,pre,varargin)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      nIn = nargin;
      data = varargin;
      obj.setDiscretizer(nIn,symmod,grid,mat,pre,data);
    end
    
    function computeFlowMat(obj)
      if obj.model.isFEMBased('Flow')
        % Compute the stiffness (H) and mass (P) matrices for the flow problem by FEM
        iiVec = zeros((obj.preP.nNodesElem.^2)*obj.preP.nE,1);
        jjVec = zeros((obj.preP.nNodesElem.^2)*obj.preP.nE,1);
        HVec = zeros((obj.preP.nNodesElem.^2)*obj.preP.nE,1);
        PVec = zeros((obj.preP.nNodesElem.^2)*obj.preP.nE,1);
        % Get the fluid compressibility
        beta = obj.material.getMaterial(2*obj.preP.nMat+1).getFluidCompressibility();
        if obj.preP.nE(2) > 0
          N1 = obj.elements.hexa.getBasisFinGPoints();
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
              N = obj.elements.tetra.getDerBasisF(el);
%               vol = getVolume(obj.elements,el);
              HLoc = N'*permMat*N*obj.elements.vol(el)/mu;
              s1 = obj.preP.nNodesElem(1)^2;
              % Computing the P matrix contribution
              PLoc = ((alpha + poro*beta)*obj.elements.vol(el)/20)*(ones(obj.preP.nNodesElem(1))...
                      + eye(obj.preP.nNodesElem(1)));
            case 12 % Hexa
              [N,dJWeighed] = obj.elements.hexa.getDerBasisFAndDet(el,1);
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
      elseif obj.model.isFVTPFABased('Flow')  % Inspired by MRST
        idIntFaces = all(obj.faces.faceNeighbors ~= 0,2);
        neigh1 = obj.faces.faceNeighbors(idIntFaces,1);
        neigh2 = obj.faces.faceNeighbors(idIntFaces,2);
        sumDiagTrans = accumarray([neigh1; neigh2], ...
          repmat(obj.trans(idIntFaces),[2,1]),[obj.mesh.nCells,1]);
        obj.H = sparse([neigh1; neigh2; (1:obj.mesh.nCells)'], ...
                       [neigh2; neigh1; (1:obj.mesh.nCells)'], ...
                       [-obj.trans(idIntFaces); -obj.trans(idIntFaces); ...
                        sumDiagTrans],obj.mesh.nCells,obj.mesh.nCells);
        poroMat = zeros(obj.preP.nMat,1);
        alphaMat = zeros(obj.preP.nMat,1);
        beta = obj.material.getMaterial(2*obj.preP.nMat+1).getFluidCompressibility();
        for m = 1:obj.preP.nMat
          poroMat(m) = obj.material.getMaterial(obj.preP.nMat+m).getPorosity();
          alphaMat(m) = obj.material.getMaterial(m).getRockCompressibility();
        end
        % (alpha+poro*beta)
        PVal = (alphaMat(obj.mesh.cellTag) + beta*poroMat(obj.mesh.cellTag)).*obj.elements.vol;
        obj.P = sparse(1:obj.mesh.nCells,1:obj.mesh.nCells,PVal,obj.mesh.nCells,obj.mesh.nCells);
      end
    end
    
    function computeFlowSystMat(obj,theta,dt)
      % Compute matrices K and K2
      obj.K = theta*obj.H + obj.P/dt;
      obj.K2 = obj.P/dt - (1-theta)*obj.H;
    end
    
    function computeFlowRHSGravContribute(obj)
      % Compute the gravity contribution
      if isFEMBased(obj.model,'Flow')
        nr = obj.mesh.nNodes;
      elseif isFVTPFABased(obj.model,'Flow')
        nr = obj.mesh.nCells;
      end
      obj.fConst = zeros(nr,1); clear nr;
      % Get the fluid specific weight and viscosity
      gamma = obj.material.getMaterial(2*obj.preP.nMat+1).getFluidSpecWeight();
      if gamma > 0
        if isFEMBased(obj.model,'Flow')
          mu = obj.material.getMaterial(2*obj.preP.nMat+1).getDynViscosity();
          for el=1:obj.mesh.nCells
            % Get the material permeability
            permMat = obj.material.getMaterial(obj.preP.nMat+obj.mesh.cellTag(el)).getPermMatrix();
            permMat = permMat/mu;
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
            obj.fConst(dof) = obj.fConst(dof) + fLoc;
          end
        elseif isFVTPFABased(obj.model,'Flow')
          idIntFaces = all(obj.faces.faceNeighbors ~= 0,2);
          neigh = obj.faces.faceNeighbors(idIntFaces,:);
%           zNeigh = obj.elements.cellCentroid(neigh,3);
          zVec = obj.elements.cellCentroid(:,3);
          zNeigh = zVec(neigh);
          dz = zNeigh(:,1) - zNeigh(:,2);
          dz = gamma*obj.trans(idIntFaces).*dz;
          obj.fConst = accumarray(neigh(:),[dz; -dz],[obj.mesh.nCells,1]);
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
            N = getDerBasisF(obj.elements.tetra,el);
            vol = findVolume(obj.elements.tetra,el);
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
              [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
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
    function  rhs = computePoroCouple(obj)
       % calculate poromechanics stiffness matrix in coupled model
       % output: right-hand side (internal forces) waiting for external
       % forces and Neumann conditions
       
       %same instructions for computePoroSyst, different Output variables
      rhs = zeros(obj.mesh.nNodes*obj.mesh.nDim,1);
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
            N = getDerBasisF(obj.elements.tetra,el);
            vol = findVolume(obj.elements.tetra,el);
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
              [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
              B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
              B(obj.preP.indB(:,2)) = N(obj.preP.indB(:,1));
              D = obj.preP.getStiffMatrix(el,state.stress(l2+1:l2+obj.GaussPts.nNode,3)+state.iniStress(l2+1:l2+obj.GaussPts.nNode,3));
              Ks = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
              Ks = Ks.*reshape(dJWeighed,1,1,[]);
              KLoc = sum(Ks,3);
              clear Ks;
              s1 = obj.preP.nEntryKLoc(2);
              sz = state.stress(l2+1:l2+obj.GaussPts.nNode,:);
              sz = reshape(sz',6,1,obj.GaussPts.nNode);
              fTmp = pagemtimes(B,'ctranspose',sz,'none'); %internal forces local calculation
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
        rhs(dof) = rhs(dof) + fLoc;
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
       iivec = zeros((obj.preP.nNodesElem.^2*obj.mesh.nDim)*obj.preP.ne,1);
       jjvec = zeros((obj.preP.nNodesElem.^2*obj.mesh.nDim)*obj.preP.ne,1);
       Qvec =  zeros((obj.preP.nNodesElem.^2*obj.mesh.nDim)*obj.preP.ne,1);
       
       l1=0;
       biot = 1; %TO DO: call to Biot Coefficient in PorousRock class
       
       if obj.preP.ne(2) > 0 %at least one Hexahedron
           N1 = getBasisFinGPoints(obj.elements); %compute Basis Function matrix for further calculations
       end
       
       for el=1:obj.mesh.nCells
           switch obj.mesh.cellVTKType(el)
               case 10 %Tetrahedrons, direct integration 
                   vol = findVolume(obj.elements.tetra,el);
                   der = getDerBasisF(obj.elements.tetra,el);
                   der = der';
                   Qloc = biot*0.25*repelem(der(:),1,4)*vol;
                   s1 = obj.preP.nNodesElem(1).^2*obj.mesh.ndim;
               case 12 %Hexahedrons, Gauss integration
                   [N,dJWeighed] = getDerBasisFandDet(obj.elements.hexa,el,1);
                   nG = obj.GaussPts.nNode;
                   iN = zeros(6,8,nG); %matrix product i*N
                   B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
                   B(obj.preP.indB(:,2)) = N(obj.preP.indB(:,1));
                   iN(1,:,:) = reshape(N1',1,8,nG);
                   iN(2:3,:,:) = repmat(iN(1,:,:),2,1,1);
                   Qs = biot*pagemtimes(B,'ctranspose',iN,'none');
                   Qs = Qs.*reshape(dJWeighed,1,1,[]);
                   Qloc = sum(Qs,3);
                   clear Qs;
                   s1 = obj.preP.nNodesElem(2)^2*obj.mesh.nDim;   
           end
       
       %assembly Coupling Matrix
       dofporo = obj.preP.getDofID(el);
       dofflow = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
       [jjloc,iiloc] = meshgrid(dofflow,dofporo);
       iivec(l1+1:l1+s1) = iiloc(:);
       jjvec(l1+1:l1+s1) = jjloc(:);
       Qvec(l1+1:l1+s1) = Qloc(:);
       
       l1 = l1+s1;
       end  
       obj.Q = sparse(iiVec,jjVec,QVec,obj.mesh.nNodes*obj.mesh.nDim, ...
                     obj.mesh.nNodes);
    end

    function computeCoupleSyst(obj,theta,dt,statek,stateTmp)
%Not necessary since they are already computed before the NR loop
%         obj.computeFlowMat(); 
%         obj.computeFlowRHSGravContribute();
        obj.computeCoupleMat;
        obj.rhsFlow = obj.fConst; %initializing flow RHS
        obj.rhsPoro = obj.computePoroCouple(); 

        %computing rhs contributes (neumann and volume forces still
        %missing!)
        
        obj.rhsPoro = obj.rhsPoro + (-obj.Q*stateTmp.pressure) - (((1-theta)/theta)*...
            (obj.Q*statek.pressure-obj.KPoro*statek.displ));
        
        obj.rhsFlow = obj.rhsFlow + (-(obj.Q)'*stateTmp.displ- ...
            (theta*dt*obj.H+obj.P)*stateTmp.pressure) - ...
            ((-obj.Q)'*statek.displ+((1-theta)*dt*obj.H-obj.P)*statek.pressure);
        
    end
        

  end
  
  
  
  
  
  methods(Access = private)
    function setDiscretizer(obj,nIn,symmod,grid,mat,pre,data)
      obj.model = symmod;
      obj.mesh = grid.topology;
      obj.elements = grid.cells;
      obj.faces = grid.faces;
      obj.material = mat;
%       obj.probType = pType;
%       obj.bound = bc;
%       obj.BCName = BCName;
      obj.preP = pre;
%       obj.state = stat;
      if nIn > 4
        obj.GaussPts = data{1};
      end
      %
      if obj.model.isFVTPFABased('Flow')
        obj.computeTrans;
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
    
%     function computeTrans(obj)
%       
%     end
    
    function computeTrans(obj)   % Inspired by MRST
      % Compute first the vector connecting each cell centroid to the
      % half-face
      r = [1, 1, 1, 2, 2, 2, 3, 3, 3];
      c = [1, 2, 3, 1, 2, 3, 1, 2, 3];
      hf2Cell = repelem((1:obj.mesh.nCells)',diff(obj.faces.mapF2E));
      L = obj.faces.faceCentroid(obj.faces.faces2Elements(:,1),:) - obj.elements.cellCentroid(hf2Cell,:);
      sgn = 2*(hf2Cell == obj.faces.faceNeighbors(obj.faces.faces2Elements(:,1))) - 1;
      N = bsxfun(@times,sgn,obj.faces.faceNormal(obj.faces.faces2Elements(:,1),:));
      KMat = zeros(obj.preP.nMat,9);
      for i=1:obj.preP.nMat
        KMat(i,:) = obj.material.getMaterial(obj.preP.nMat+i).getPermVector();
      end
      hT = zeros(length(hf2Cell),1);
      for k=1:length(r)
        hT = hT + L(:,r(k)) .* KMat(obj.mesh.cellTag(hf2Cell),k) .* N(:,c(k));
      end
      hT = hT./sum(L.*L,2);
      mu = obj.material.getMaterial(2*obj.preP.nMat+1).getDynViscosity();
      hT = hT/mu;
      %
      obj.trans = 1 ./ accumarray(obj.faces.faces2Elements(:,1),1 ./ hT,[obj.faces.nFaces,1]);
    end
  end
end