classdef SPFlow < SinglePhysics
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

   methods (Access = public)
      function obj = SPFlow(symmod,params,dofManager,grid,mat,data,varargin)
         if isempty(varargin)
            fieldName = 'SPFlow';
         else
            fieldName = varargin{1};
         end
         obj@SinglePhysics(fieldName,symmod,params,dofManager,grid,mat,data);
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

      function state = setState(obj,state)
         if isFEMBased(obj.model,'Flow')
            n = obj.mesh.nNodes;
         elseif isFVTPFABased(obj.model,'Flow')
            n = obj.mesh.nCells;
            state.flux = zeros(n,3);
         end
         state.pressure = zeros(n,1);
      end

      function state = updateState(obj,state,dSol)
         ents = obj.dofm.getActiveEnts(obj.field);
         state.pressure(ents) = state.pressure(ents) + dSol(obj.dofm.getDoF(obj.field));
      end

      function state = finalizeState(obj,bound,state)
         % Compute the posprocessing variables for the module.
         state.potential = computePotential(obj,state.pressure);

         mu = (1/obj.material.getFluid().getDynViscosity())*ones(obj.faces.nFaces,1);
         mu = mu(obj.isIntFaces);
         state.perm = printPermeab(obj);
         state.flux = computeFlux(obj,mu,bound,state.potential,state.t);
      end

      function [cellData,pointData] = printState(obj,bound,sOld,sNew,t)
         % append state variable to output structure
         state = SPFlow.buildState([]);
         switch nargin
            case 3
               state.pressure = sOld.pressure;
               state.t = sOld.t;
            case 5
               % linearly interpolate the pressure containing print time
               fac = (t - sOld.t)/(sNew.t - sOld.t);
               state.pressure = sNew.pressure*fac+sOld.pressure*(1-fac);
               state.t = sNew.t*fac+sOld.t*(1-fac);
            otherwise
               error('Wrong number of input arguments');
         end
         % posprocessing the structure of VSFlow.
         state = finalizeState(obj,bound,state);
         [cellData,pointData] = SPFlow.buildPrintStruct(obj.model,state);
      end

      function state = computeMat(obj,state,~,dt)
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
         % dealing with input params
         subCells = obj.dofm.getFieldCells(obj.field);
         nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]);
         % Compute the stiffness (H) and mass (P) matrices for the flow problem by FEM
         iiVec = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
         jjVec = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
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
            permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
            poro = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPorosity();
            if ismember(obj.mesh.cellTag(el),getFieldCellTags(obj.dofm,{obj.field,'Poromechanics'}))
               alpha = 0; %this term is not needed in coupled formulation
            else
               alpha = obj.material.getMaterial(obj.mesh.cellTag(el)).ConstLaw.getRockCompressibility();
               %solid skeleton contribution to storage term as oedometric compressibility .
            end             % Compute the element matrices based on the element type
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
            dof = dofId(obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)),1);
            [jjLoc,iiLoc] = meshgrid(dof,dof);
            iiVec(l1+1:l1+s1) = iiLoc(:);
            jjVec(l1+1:l1+s1) = jjLoc(:);
            HVec(l1+1:l1+s1) = HLoc(:);
            PVec(l1+1:l1+s1) = PLoc(:);
            l1 = l1 + s1;
         end
         % renumber indices according to active nodes
         [~,~,iiVec] = unique(iiVec);
         [~,~,jjVec] = unique(jjVec);
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
         neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
         % Transmissibility of internal faces
         tmpVec = lw.*obj.trans(obj.isIntFaces);
         nneigh = length(tmpVec);
         % [~,~,reorder] = unique([neigh(:,1); neigh(:,2); subCells]);
         [~,~,reorder] = unique([neigh(:,1); neigh(:,2)]);
         neigh1 = reorder(1:nneigh);
         neigh2 = reorder(nneigh+1:2*nneigh);
         sumDiagTrans = accumarray( [neigh1;neigh2], repmat(tmpVec,[2,1]), ...
             [nSubCells,1]);
         % Assemble H matrix
         nDoF = obj.dofm.getNumDoF(obj.field);
         obj.H = sparse([neigh1; neigh2; (1:nSubCells)'],...
             [neigh2; neigh1; (1:nSubCells)'],...
             [-tmpVec; -tmpVec; sumDiagTrans], nDoF, nDoF);
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
         PVal = PVal.*obj.elements.vol(subCells);
         nDoF = obj.dofm.getNumDoF(obj.field);
         [~,~,dof] = unique(subCells);
         obj.P = sparse(dof,dof,PVal,nDoF,nDoF);
      end

      function stateTmp = computeRhs(obj,stateTmp,statek,dt)
         % Compute the residual of the flow problem
         lw = obj.material.getFluid().getDynViscosity();
         ents = obj.dofm.getActiveEnts(obj.field);
         if ~obj.simParams.isTimeDependent
            obj.rhs = obj.H*stateTmp.pressure(ents);
         else
            theta = obj.simParams.theta;
            rhsStiff = theta*obj.H*stateTmp.pressure(ents) + (1-theta)*obj.H*statek.pressure(ents);
            rhsCap = (obj.P/dt)*(stateTmp.pressure(ents) - statek.pressure(ents));
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
                        rhsLoc = (N'*permMat(:,3))*obj.elements.vol(el)*gamma;
                     case 12 % Hexa
                        [N,dJWeighed] = obj.elements.hexa.getDerBasisFAndDet(el,1);
                        fs = pagemtimes(N,'ctranspose',permMat(:,3),'none');
                        fs = fs.*reshape(dJWeighed,1,1,[]);
                        rhsLoc = sum(fs,3)*gamma;
                  end
                  %
                  entsId = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
                  rhsTmp(entsId) = rhsTmp(entsId) + rhsLoc;
               end
               obj.rhsGrav = rhsTmp(obj.dofm.getActiveEnts(obj.field));
            elseif isFVTPFABased(obj.model,'Flow')
               neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
               zVec = obj.elements.cellCentroid(:,3);
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

      function [dof,vals] = getBC(obj,bc,id,t,state)
         switch bc.getCond(id)
            case {'NodeBC','ElementBC'}
               ents = bc.getEntities(id);
               vals = bound.getVals(id,t);
            case 'SurfBC'
               v = bc.getVals(id,t);
               if isFVTPFABased(obj.model,'Flow')
                  faceID = bc.getEntities(id);
                  ents = sum(obj.faces.faceNeighbors(faceID,:),2);
                  %[ents,~,ind] = unique(ents1);
                  switch bc.getType(id)
                     case 'Neu'
                        vals = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;
                     case 'Dir'
                       gamma = obj.material.getFluid().getFluidSpecWeight();
                       mu = obj.material.getFluid().getDynViscosity();
                       tr = obj.getFaceTransmissibilities(faceID);
                       % q(ind) = 1/mu*tr.*((state.pressure(ents) - v)...
                       % + gamma*(obj.elements.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID(ind),3)));
                       q = 1/mu*tr.*((state.pressure(ents) - v)...
                         + gamma*(obj.elements.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID,3)));
                       %vals = [1/mu*tr,accumarray(ind,q)]; % {JacobianVal,rhsVal]
                       vals = [1/mu*tr,q];
                  end
               elseif isFEMBased(obj.model,'Flow')
                  ents = bc.getLoadedEntities(id);
                  entitiesInfl = bc.getEntitiesInfluence(id);
                  vals = entitiesInfl*v;
               end
            case 'VolumeForce'
               v = bc.getVals(id,t);
               ents = bc.getEntities(id);
               if isFVTPFABased(obj.model,'Flow')
                  vals = v.*obj.elements.vol(ents);
               elseif isFEMBased(obj.model,'Flow')
                  entitiesInfl = bc.getEntitiesInfluence(id);
                  vals = entitiesInfl*v;
               end
         end
         % get local dof numbering
         dof = obj.dofm.getLocalDoF(ents,obj.field);
      end

      function state = applyDirVal(obj,bc,id,t,state)
         if isFVTPFABased(obj.model,'Flow')
            % Dirichlet BCs cannot be directly applied to the solution
            % vector
            return
         end
         switch bc.getCond(id)
            case 'NodeBC'
               ents = bc.getEntities(id);
               vals = bc.getVals(id,t);
            case 'SurfBC'
               ents = bc.getLoadedEntities(id);
               s2n = bc.getEntitiesInfluence(id);
               vals = s2n*bc.getVals(id,t);
               % node id contained by constrained surface
            otherwise
               error('BC type %s is not available for %s field',cond,obj.field);
         end
         dof = bc.getCompEntities(id,ents);
         state.pressure(dof) = vals;
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

      function potential = computePotential(obj,pressure)
          %COMPUTEFLUX - compute the potential for the cell or element.
          potential = pressure;
          gamma = obj.material.getFluid().getFluidSpecWeight();
          if gamma > 0
              if isFEMBased(obj.model,'Flow')
                  potential = potential + gamma*obj.mesh.coordinates(:,3);
              elseif isFVTPFABased(obj.model,'Flow')
                  potential = potential + gamma*obj.elements.cellCentroid(:,3);
              end
          end
      end

      function flux = computeFlux(obj,mob,bound,pot,t)
          %COMPUTEFLUX - compute the flux at the faces, than accumulate
          %the value at the nodes.
          flux = zeros(obj.mesh.nNodes,3);
          if isFEMBased(obj.model,'Flow')
              % flux = fluidPot + gamma*obj.mesh.coordinates(:,3);
          elseif isFVTPFABased(obj.model,'Flow')
              nnodesBfaces = diff(obj.faces.mapN2F);
              neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
              Node2Face = repelem((1:obj.faces.nFaces)',nnodesBfaces);
              sgn = 2*(obj.faces.faceNeighbors(:,1)==0) - 1;

              fluxFaces = zeros(obj.faces.nFaces,1);              
              fluxFaces(obj.isIntFaces) = pot(neigh(:,1))-pot(neigh(:,2));
              fluxFaces(obj.isIntFaces) = mob.*obj.trans(obj.isIntFaces).*fluxFaces(obj.isIntFaces);

              faceUnit = obj.faces.faceNormal./vecnorm(obj.faces.faceNormal,2,2);
              % faceUnit = obj.faces.neighNormal;
              areaSq = vecnorm(obj.faces.faceNormal,2,2);
              areaSq = areaSq.*nnodesBfaces;
              fluxFaces = fluxFaces./areaSq;
              fluxFaces = fluxFaces.*faceUnit;  % flux at the faces

              % Contribution
              fluxFaces = repelem(fluxFaces,nnodesBfaces,1);
              axis = ones(length(obj.faces.nodes2Faces),1);
              flux = accumarray([[obj.faces.nodes2Faces axis]; ...
                  [obj.faces.nodes2Faces 2*axis]; ...
                  [obj.faces.nodes2Faces 3*axis]], fluxFaces(:));

              % add boundary condition              
              bcList = bound.db.keys;
              for bc = string(bcList)
                  field = translatePhysic(bound.getPhysics(bc),obj.model);
                  for f = field
                      if (field == "SPFlow") || (field == "VSFlow")
                          v = bound.getVals(bc,t);
                          switch bound.getCond(bc)
                              case {'NodeBC','ElementBC'}
                              case 'SurfBC'
                                  faceID = bound.getEntities(bc);
                                  switch bound.getType(bc)
                                      case 'Neu'
                                          vals = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;
                                      case 'Dir'
                                          ents = sum(obj.faces.faceNeighbors(faceID,:),2);
                                          gamma = obj.material.getFluid().getFluidSpecWeight();
                                          potBd = pot(ents)-(v+gamma*obj.faces.faceCentroid(faceID,3));
                                          vals = mob(ents).*obj.trans(faceID).*potBd;
                                  end
                                  vals = sgn(faceID).*vals./areaSq(faceID).*faceUnit(faceID,:);
                                  % sml = nnodesBfaces(faceID);
                                  vals = -repelem(vals,nnodesBfaces(faceID),1);
                                  nodes = obj.faces.nodes2Faces(ismember(Node2Face,faceID));
                                  [loc,~,pos] = unique(nodes);
                                  axis = ones(length(nodes),1);
                                  fluxB = accumarray([[pos axis]; [pos 2*axis]; ...
                                      [pos 3*axis]], vals(:));
                                  flux(loc,:)=fluxB;
                              case 'VolumeForce'
                                  % Find the cell to apply the boundary condition
                                  cellID = bound.getEntities(bc);
                                  vals = v.*obj.elements.vol(cellID);
                                  facesBcell = diff(obj.faces.mapF2E);

                                  % Find the faces to distribute the contribution.
                                  vals = vals./facesBcell(cellID);
                                  vals = repelem(vals,facesBcell(cellID),1);

                                  hf2Cell = repelem((1:obj.mesh.nCells)',facesBcell);                                  
                                  faceID = obj.faces.faces2Elements(hf2Cell == cellID,1);

                                  vals = sgn(faceID).*vals./areaSq(faceID).*faceUnit(faceID,:);
                                  vals = -repelem(vals,nnodesBfaces(faceID),1);
                                  nodes = obj.faces.nodes2Faces(ismember(Node2Face,faceID));
                                  [loc,~,pos] = unique(nodes);
                                  axis = ones(length(nodes),1);
                                  fluxB = accumarray([[pos axis]; [pos 2*axis]; ...
                                      [pos 3*axis]], vals(:));
                                  flux(loc,:)=flux(loc,:)+fluxB;
                          end
                      end
                  end
              end
          end
      end

      function perm = printPermeab(obj)
          %printPropState - print the potential for the cell or element.
          perm = zeros(obj.mesh.nCells,6);
          for el=1:obj.mesh.nCells
              ktmp = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
              perm(el,1)=ktmp(1,1);
              perm(el,2)=ktmp(2,2);
              perm(el,3)=ktmp(3,3);
              perm(el,4)=ktmp(1,2);
              perm(el,5)=ktmp(2,3);
              perm(el,6)=ktmp(1,3);
          end
      end

   end

   methods (Access=private)
      % function dof = getBCdofs(obj,bc,id)
      %    switch bc.getCond(id)
      %       case 'NodeBC'
      %          ents = bc.getEntities(id,obj.field);
      %       case 'SurfBC'
      %          ents = bc.getLoadedEntities(id);
      %          % node id contained by constrained surface
      %       otherwise
      %          error('BC type %s is not available for %s field',cond,obj.field);
      %    end
      %    % map entities dof to local dof numbering
      %    dof = obj.dofm.getLocalDoF(ents,obj.field);
      %    switch bc.getType(id)
      %       case 'Dir'
      %          % component multiplication of BC dofs
      %          dof = bc.getCompEntities(id,dof);
      %       case 'Neu'
      %          dir = obj.getDirection(identifier);
      %          c = find(strcmp(['x','y','z'],dir));
      %          dof = c*dof;
      %    end
      % end
      %
      % function vals = getBCVals(obj,bc,id,t)
      %    if strcmp(bc.getType(id),'Dir')
      %       vals = [];
      %       return
      %    end
      %    vals = bc.getVals(id,t);
      %    if strcmp(bc.getCond(id),'SurfBC')
      %       entInfl = bc.getEntitiesInfluence(id);
      %       vals = entInfl*vals;
      %    end
      % end
   end

   methods (Static)
      function state = buildState(state)
         % Constructor for the all variables avaible in this module
         state.t = [];
         state.pressure = [];
         state.potential = [];
         state.flux = [];
         state.perm = [];
      end

      function [cellStr,pointStr] = buildPrintStruct(mod,state)
         if isFEMBased(mod,'Flow')
            cellStr = repmat(struct('name', 1, 'data', 1), 1, 1);
            cellStr(1).name = 'permeability';
            cellStr(1).data = state.perm;

            pointStr = repmat(struct('name', 1, 'data', 1), 2, 1);
            pointStr(1).name = 'pressure';
            pointStr(1).data = state.pressure;
            pointStr(2).name = 'potential';
            pointStr(2).data = state.potential;
         elseif isFVTPFABased(mod,'Flow')
            pointStr = repmat(struct('name', 1, 'data', 1), 1, 1);
            pointStr(1).name = 'flux';
            pointStr(1).data = state.flux;
            % pointStr = [];

            cellStr = repmat(struct('name', 1, 'data', 1), 3, 1);
            cellStr(1).name = 'pressure';
            cellStr(1).data = state.pressure;
            cellStr(2).name = 'potential';
            cellStr(2).data = state.potential;
            cellStr(3).name = 'permeability';
            cellStr(3).data = state.perm;
            % cellStr(4).name = 'flux';
            % cellStr(4).data = state.flux;
         end
      end
   end
end

