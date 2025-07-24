classdef SinglePhaseFlow < SinglePhysics
   %POROMECHANICS
   % Subclass of Discretizer
   % Implement Poromechanics methods to assemble the stiffness matrix and
   % the residual contribution

   properties
      trans
      isIntFaces
      rhsGrav         % gravity contribution to rhs
      H
      P
   end

   properties (Constant)
     field = 'SinglePhaseFlow'
   end

   methods (Access = public)
      function obj = SinglePhaseFlow(symmod,params,dofManager,grid,mat,state)
         obj@SinglePhysics(symmod,params,dofManager,grid,mat,state);
         if obj.model.isFVTPFABased('Flow')
            obj.computeTrans;
            %get cells with active flow model
            flowCells = obj.dofm.getFieldCells(obj.field);
            % Find internal faces (i.e. shared by two cells with active
            % flow)
            obj.isIntFaces = all(ismember(obj.faces.faceNeighbors, flowCells), 2);
         end
         computeRHSGravTerm(obj);
      end

      function initState(obj)
         if isFEMBased(obj.model,'Flow')
            n = obj.mesh.nNodes;
         elseif isFVTPFABased(obj.model,'Flow')
            n = obj.mesh.nCells;
         end
         obj.state.data.pressure = zeros(n,1);
      end

      function updateState(obj,dSol)
         ents = obj.dofm.getActiveEnts(obj.field);
         obj.state.data.pressure(ents) = obj.state.data.pressure(ents) + dSol(obj.dofm.getDoF(obj.field));
      end

      function fluidPot = finalizeState(obj,stateIn)
         fluidPot = stateIn.data.pressure;
         gamma = obj.material.getFluid().getFluidSpecWeight();
         if gamma > 0
            if isFEMBased(obj.model,'Flow')
               fluidPot = fluidPot + gamma*obj.mesh.coordinates(:,3);
            elseif isFVTPFABased(obj.model,'Flow')
               fluidPot = fluidPot + gamma*obj.mesh.cellCentroid(:,3);
            end
         end
      end

      function var = getState(obj,varargin)
        % input: state structure
        % output: current primary variable
        if isempty(varargin)
          var = obj.state.data.pressure;
        else
          stateIn = varargin{1};
          var = stateIn.data.pressure;
        end
      end

      function [cellData,pointData] = printState(obj,sOld,sNew,t)
         % append state variable to output structure
         switch nargin
            case 2
               fluidPot = finalizeState(obj,sOld);
               pressure = sOld.data.pressure;
            case 4
               % linearly interpolate state variables containing print time
               fac = (t - sOld.t)/(sNew.t - sOld.t);
               fluidPotOld = finalizeState(obj,sOld);
               fluidPotNew = finalizeState(obj,sNew);
               fluidPot = fluidPotNew*fac+fluidPotOld*(1-fac);
               pressure = sNew.data.pressure*fac+sOld.data.pressure*(1-fac);
            otherwise
               error('Wrong number of input arguments');
         end
         [cellData,pointData] = SinglePhaseFlow.buildPrintStruct(obj.model,pressure,fluidPot);
      end

      function computeMat(obj,~,dt)
         % recompute elementary matrices only if the model is linear
         if ~isLinear(obj) || isempty(obj.J)
            if obj.model.isFEMBased('Flow')
               computeMatFEM(obj);
            elseif obj.model.isFVTPFABased('Flow')
               mu = obj.material.getFluid().getDynViscosity();
               computeStiffMatFV(obj,1/mu);
               computeCapMatFV(obj);
            end
         end
         if obj.simParams.isTimeDependent
            obj.J = obj.simParams.theta*obj.H + obj.P/dt;
         else
            obj.J = obj.H;
         end
      end

      function computeMatFEM(obj)

        subCells = obj.dofm.getFieldCells(obj.field);
        nEntries = sum(obj.mesh.cellNumVerts(subCells).^2);

        [iiVec,jjVec,HVec,PVec] = deal(zeros(nEntries,1));

        % Get the fluid compressibility
        beta = obj.material.getFluid().getFluidCompressibility();

        % Get the fluid dynamic viscosity
        mu = obj.material.getFluid().getDynViscosity();
        %
        l1 = 0;
        for el = subCells'
          permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
          poro = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPorosity();
          alpha = getRockCompressibility(obj,el);
          % Compute the element matrices based on the element type
          % (tetrahedra vs. hexahedra)
          elem = getElement(obj.elements,obj.mesh.cellVTKType(el));
          [gradN,dJWeighed] = getDerBasisFAndDet(elem,el,1);
          N = getBasisFinGPoints(elem);
          permMat = permMat/mu;
          Hs = pagemtimes(pagemtimes(gradN,'ctranspose',permMat,'none'),gradN);
          Hs = Hs.*reshape(dJWeighed,1,1,[]);
          HLoc = sum(Hs,3);
          clear Hs;
          s1 = numel(HLoc);
          % Computing the P matrix contribution
          PLoc = (alpha+poro*beta)*(N'*diag(dJWeighed)*N);
          %Getting dof associated to Flow subphysic
          nodes = (obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
          dof = obj.dofm.getLocalDoF(nodes,obj.fldId);
          [jjLoc,iiLoc] = meshgrid(dof,dof);
          iiVec(l1+1:l1+s1) = iiLoc(:);
          jjVec(l1+1:l1+s1) = jjLoc(:);
          HVec(l1+1:l1+s1) = HLoc(:);
          PVec(l1+1:l1+s1) = PLoc(:);
          l1 = l1 + s1;
        end
        % renumber indices according to active nodes
        nDoF = obj.dofm.getNumDoF(obj.field);
        % Assemble H and P matrices defined as new fields of
        obj.H = sparse(iiVec, jjVec, HVec, nDoF, nDoF);
        obj.P = sparse(iiVec, jjVec, PVec, nDoF, nDoF);
      end


      function computeStiffMatFV(obj,lw)
         % Inspired by MRST
         subCells = obj.dofm.getFieldCells(obj.field);
         nSubCells = length(subCells);
         %get pairs of faces that contribute to the subdomain
         neigh1 = obj.faces.faceNeighbors(obj.isIntFaces,1);
         neigh2 = obj.faces.faceNeighbors(obj.isIntFaces,2);
         % Transmissibility of internal faces
         tmpVec = lw.*obj.trans(obj.isIntFaces);
         % tmpVec = lw.*tmpVec;
         [~,~,neigh1] = unique(neigh1);
         [~,~,neigh2] = unique(neigh2);
         sumDiagTrans = accumarray([neigh1; neigh2], ...
            repmat(tmpVec,[2,1]),[nSubCells,1]);
         % Assemble H matrix
         nDoF = obj.dofm.getNumDoF(obj.field);
         obj.H = sparse([neigh1; neigh2; (1:nSubCells)'], ...
            [neigh2; neigh1; (1:nSubCells)'], ...
            [-tmpVec; -tmpVec; ...
            sumDiagTrans],nDoF,nDoF);
      end



      function computeCapMatFV(obj,varargin)
         subCells = obj.dofm.getFieldCells(obj.field);
         nSubCells = length(subCells);
         poroMat = zeros(nSubCells,1);
         alphaMat = zeros(nSubCells,1);
         beta = obj.material.getFluid().getFluidCompressibility();
         for m = 1:obj.mesh.nCellTag
            if ~ismember(m,obj.dofm.getFieldCellTags({obj.field,'Poromechanics'}))
               % compute alpha only if there's no coupling in the
               % subdomain
               alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
            end
            poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
         end
         % (alpha+poro*beta)
         PVal = alphaMat(obj.mesh.cellTag(subCells)) + beta*poroMat(obj.mesh.cellTag(subCells));
         if ~isempty(varargin)
             % variably saturated flow model
             PVal = PVal.*varargin{1} + poroMat(obj.mesh.cellTag(subCells)).*varargin{2};
         end
         PVal = PVal.*obj.mesh.cellVolume(subCells);
         nDoF = obj.dofm.getNumDoF(obj.field);
         [~,~,dof] = unique(subCells);
         obj.P = sparse(dof,dof,PVal,nDoF,nDoF);
      end


      function computeRhs(obj,stateOld,dt)
         % Compute the residual of the flow problem
         lw = obj.material.getFluid().getDynViscosity();
         ents = obj.dofm.getActiveEnts(obj.field);
         if ~obj.simParams.isTimeDependent
            obj.rhs = obj.H*obj.state.data.pressure(ents);
         else
            theta = obj.simParams.theta;
            rhsStiff = theta*obj.H*obj.state.data.pressure(ents) + (1-theta)*obj.H*stateOld.data.pressure(ents);
            rhsCap = (obj.P/dt)*(obj.state.data.pressure(ents) - stateOld.data.pressure(ents));
            obj.rhs = rhsStiff + rhsCap;
         end
         gamma = obj.material.getFluid().getFluidSpecWeight();
         %adding gravity rhs contribute
         if gamma > 0
            if isFEMBased(obj.model,'Flow')
               obj.rhs = obj.rhs + obj.rhsGrav;
            elseif isFVTPFABased(obj.model,'Flow')
               obj.rhs = obj.rhs + finalizeRHSGravTerm(obj,lw);
            end
         end
      end

      function computeRHSGravTerm(obj)
         % Compute the gravity contribution
         % Get the fluid specific weight and viscosity'
         rhsTmp = zeros(obj.dofm.getNumDoF(obj.field),1);
         gamma = obj.material.getFluid().getFluidSpecWeight();
         if gamma > 0
            subCells = obj.dofm.getFieldCells(obj.field);
            if isFEMBased(obj.model,'Flow')
               for el = subCells'
                  % Get the material permeability
                  permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
                  %             permMat = permMat/mu;
                  switch obj.mesh.cellVTKType(el)
                     case 10 % Tetrahedra
                        N = obj.elements.tetra.getDerBasisF(el);
                        rhsLoc = (N'*permMat(:,3))*obj.mesh.cellVolume(el)*gamma;
                     case 12 % Hexa
                        [N,dJWeighed] = obj.elements.hexa.getDerBasisFAndDet(el,1);
                        fs = pagemtimes(N,'ctranspose',permMat(:,3),'none');
                        fs = fs.*reshape(dJWeighed,1,1,[]);
                        rhsLoc = sum(fs,3)*gamma;
                  end
                  %
                  entsId = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
                  rhsTmp = rhsTmp(entsId) + rhsLoc;
               end
               obj.rhsGrav = rhsTmp(obj.dofm.getActiveEnts(obj.field));
            elseif isFVTPFABased(obj.model,'Flow')
               neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
               zVec = obj.mesh.cellCentroid(:,3);
               zNeigh = zVec(neigh);
               obj.rhsGrav = gamma*obj.trans(obj.isIntFaces).*(zNeigh(:,1) - zNeigh(:,2));
            end
         end
         % remove inactive components of rhs vector
      end

      function gTerm = finalizeRHSGravTerm(obj,lw)
         nCells = obj.dofm.getNumDoF(obj.field);
         neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
         gTerm = accumarray(neigh(:),[lw.*obj.rhsGrav; ...
            -lw.*obj.rhsGrav],[nCells,1]);
         gTerm = gTerm(obj.dofm.getActiveEnts(obj.field));
      end

      function [dof,vals] = getBC(obj,bc,id,t)
         switch bc.getCond(id)
            case {'NodeBC','ElementBC'}
               ents = bc.getEntities(id);
               vals = bc.getVals(id,t);
            case 'SurfBC'
               v = bc.getVals(id,t);
               if isFVTPFABased(obj.model,'Flow')
                  faceID = bc.getEntities(id);
                  ents = sum(obj.faces.faceNeighbors(faceID,:),2);
                  [ents,~,ind] = unique(ents);
                  switch bc.getType(id)
                     case 'Neu'
                        area = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;
                        vals = accumarray(ind, area);
                     case 'Dir'
                        gamma = obj.material.getFluid().getFluidSpecWeight();
                        mu = obj.material.getFluid().getDynViscosity();
                        tr = obj.getFaceTransmissibilities(faceID);
                        q = 1/mu*tr.*((obj.state.data.pressure(ents) - v)...
                           + gamma*(obj.mesh.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID,3)));
                        vals = [1/mu*tr,accumarray(ind,q)]; % {JacobianVal,rhsVal]
                  end
               elseif isFEMBased(obj.model,'Flow')
                  ents = bc.getLoadedEntities(id);
                  entitiesInfl = bc.getEntitiesInfluence(id);
                  vals = entitiesInfl*v;
               end
           case 'VolumeForce'
             v = bc.getVals(id,t);
             if isFVTPFABased(obj.model,'Flow')
               ents = bc.getEntities(id);
               vals = v.*obj.mesh.cellVolume(ents);
             elseif isFEMBased(obj.model,'Flow')
               ents = bc.getLoadedEntities(id);
               entitiesInfl = bc.getEntitiesInfluence(id);
               vals = entitiesInfl*v;
             end
         end
         % get local dof numbering
         dof = obj.dofm.getLocalDoF(ents,obj.fldId);
      end

      function applyDirVal(obj,dof,vals)
         if isFVTPFABased(obj.model,'Flow')
            % Dirichlet BCs cannot be directly applied to the solution
            % vector
            return
         end
         obj.state.data.pressure(dof) = vals;
      end

      function applyDirBC(obj,~,ents,vals)
         % apply Dirichlet BCs
         % ents: id of constrained faces without any dof mapping applied
         % vals(:,1): Jacobian BC contrib vals(:,2): rhs BC contrib
         if isFVTPFABased(obj.model,'Flow')
            % BCs imposition for finite volumes
            assert(size(vals,2)==2,'Invalid matrix size for BC values');
            nDoF = obj.dofm.getNumDoF(obj.field);
            obj.J(nDoF*(ents-1) + ents) = obj.J(nDoF*(ents-1) + ents) + vals(:,1);
            obj.rhs(ents) = obj.rhs(ents) + vals(:,2);
         else
            % strong nodal BCs imposition
            applyDirBC@SinglePhysics(obj,obj.field,ents)
         end
      end

      function out = isLinear(obj)
         out = true;
      end


      function trans = getFaceTransmissibilities(obj,faceID)
         trans = obj.trans(faceID);
      end

      function computeTrans(obj)   % Inspired by MRST
         % Compute first the vector connecting each cell centroid to the
         % half-face
         r = [1, 1, 1, 2, 2, 2, 3, 3, 3];
         c = [1, 2, 3, 1, 2, 3, 1, 2, 3];
         hf2Cell = repelem((1:obj.mesh.nCells)',diff(obj.faces.mapF2E));
         L = obj.faces.faceCentroid(obj.faces.faces2Elements(:,1),:) - obj.mesh.cellCentroid(hf2Cell,:);
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

      function alpha = getRockCompressibility(obj,el)
        if ismember(obj.mesh.cellTag(el),getFieldCellTags(obj.dofm,{obj.field,'Poromechanics'}))
          alpha = 0; %this term is not needed in coupled formulation
        else
          if isfield(obj.material.getMaterial(obj.mesh.cellTag(el)),"ConstLaw")
            %solid skeleton contribution to storage term as oedometric compressibility .
            alpha = obj.material.getMaterial(obj.mesh.cellTag(el)).ConstLaw.getRockCompressibility();
          else
            alpha = 0;
          end
        end
      end

   end


   methods (Static)
      function [cellStr,pointStr] = buildPrintStruct(mod,press,pot)
         if isFEMBased(mod,'Flow')
           cellStr = [];
           pointStr = repmat(struct('name', 1, 'data', 1), 2, 1);
           pointStr(1).name = 'pressure';
           pointStr(1).data = press;
           pointStr(2).name = 'potential';
           pointStr(2).data = pot;
         elseif isFVTPFABased(mod,'Flow')
           pointStr = [];
           cellStr = repmat(struct('name', 1, 'data', 1), 2, 1);
           cellStr(1).name = 'pressure';
           cellStr(1).data = press;
           cellStr(2).name = 'potential';
           cellStr(2).data = pot;
         end
      end

      function out = getField()
        out = SinglePhaseFlow.field;
      end
   end
end

