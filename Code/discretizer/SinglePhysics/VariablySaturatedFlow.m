classdef VariablySaturatedFlow < SinglePhaseFlow
    % Variably Saturated flow
    % Subclass of SPFlow
    % Implements Richards equations for unsaturated flow in vadose region
    % It is subclass of SPFlow since most of the methods are the same

    properties
        lwkpt           % mobility
        JNewt = []      % newton jacobian contribution
        upElem          % upstream elements array for each face
    end

    methods (Access = public)
        function obj = VariablySaturatedFlow(symmod,params,dofManager,grid,mat,bc,state)
            % initialize as SPFlow class
            obj@SinglePhaseFlow(symmod,params,dofManager,grid,mat,bc,state);
        end

        function computeMat(obj,stateOld,dt)
          theta = obj.simParams.theta;
          pkpt = theta*obj.state.data.pressure + ...
            (1 - theta)*stateOld.data.pressure;
          % pkpt = obj.simParams.theta*obj.state.data.pressure + ...
          %   (1 - obj.simParams.theta)*stateOld.data.pressure;
          [Swkpt,dSwkpt,d2Swkpt, obj.lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt);
          computeStiffMatFV(obj,obj.lwkpt);
          computeCapMatFV(obj,Swkpt,dSwkpt);
          obj.J = theta*obj.H + obj.P/dt;
          if isNewtonNLSolver(obj.simParams)
            % computeNewtPartOfJacobian(obj,dt,stateOld,pkpt,dSwkpt,d2Swkpt,dlwkpt)
            Jh = computeJacobianPartJhNewton(obj,pkpt,dlwkpt);
            Jp = computeJacobianPartJpNewton(obj,dt,obj.state.data.pressure,stateOld.data.pressure,Swkpt,dSwkpt,d2Swkpt);
            obj.JNewt = Jp+Jh;
            obj.J = obj.J + obj.JNewt;
          end
        end

        function computeRhs(obj,stateOld,dt)
            % Compute the residual of the flow problem
            theta = obj.simParams.theta;
            ents = obj.dofm.getActiveEnts(obj.getField());

            % rhsStiff = theta*obj.H*obj.state.data.pressure(ents) + (1-theta)*obj.H*stateOld.data.pressure(ents);
            % rhsCap = (obj.P/dt)*(obj.state.data.pressure(ents) - stateOld.data.pressure(ents));
            % obj.rhs = rhsStiff + rhsCap;
            pkpt = theta*obj.state.data.pressure(ents) + (1-theta)*stateOld.data.pressure(ents);
            pkdiff = obj.state.data.pressure(ents) - stateOld.data.pressure(ents);
            obj.rhs = obj.H*pkpt + (obj.P/dt)*pkdiff;

            gamma = obj.material.getFluid().getFluidSpecWeight();
            % adding gravity rhs contribute
            if gamma > 0
               obj.rhs = obj.rhs + finalizeRHSGravTerm(obj,obj.lwkpt);
               % mu = obj.material.getFluid().getDynViscosity();
               % lwkpt = mu*ones(length(obj.lwkpt),1);
               % lwkpt = 1./obj.lwkpt;
               % lwkpt = obj.lwkpt;
               % obj.rhs = obj.rhs + finalizeRHSGravTerm(obj,lwkpt);
            end
        end

        function initState(obj)
            n = obj.mesh.nCells;
            obj.state.data.pressure = zeros(n,1);
            obj.state.data.saturation = zeros(n,1);
        end

        function updateState(obj,dSol)
          if nargin > 1
            ents = obj.dofm.getActiveEnts(obj.getField());
            obj.state.data.pressure(ents) = obj.state.data.pressure(ents) + dSol(obj.dofm.getDoF(obj.getField()));
            obj.state.data.saturation = computeSaturation(obj,obj.state.data.pressure);
          end
        end

        function states = finalizeState(obj,states)
            % Compute the posprocessing variables for the module.
            pressure = states.pressure;
            states.potential = computePotential(obj,pressure);
            states.saturation = computeSaturation(obj,pressure);
            [mob ,~] = computeMobility(obj,pressure);
            % states.flux = computeFlux(obj,mob,pressure);
            % flux = computeFluxBound(obj,flux,bound,potential,pressure,t);
            % mass = checkMassCons(obj,mob,potential);
            states.perm = printPermeab(obj);
        end

        function [cellData,pointData] = printState(obj,sOld,sNew,t)
            % append state variable to output structure
            outPrint = [];
            switch nargin
                case 2
                    outPrint.pressure = sOld.data.pressure;
                case 4
                    % linearly interpolate state variables containing print time
                    fac = (t - sOld.t)/(sNew.t - sOld.t);
                    outPrint.pressure = sNew.data.pressure*fac+sOld.data.pressure*(1-fac);
                otherwise
                    error('Wrong number of input arguments');
            end
            % posprocessing the structure of VSFlow.
            outPrint = finalizeState(obj,outPrint);
            [cellData,pointData] = VariablySaturatedFlow.buildPrintStruct(outPrint);
        end

        function out = isLinear(obj)
            out = false;
        end

        function [dof,vals] = getBC(obj,id,t)
           switch obj.bcs.getCond(id)
              case {'NodeBC','ElementBC'}
                 ents = obj.bcs.getEntities(id);
                 vals = obj.bcs.getVals(id,t);
              case 'SurfBC'
                 % faceID = obj.bcs.getEntities(id);
                 % ents = sum(obj.faces.faceNeighbors(faceID,:),2);
                 % v = obj.bcs.getVals(id,t);
                 % [ents,~,ind] = unique(ents);

                 [faceID, faceOrder] = sort(obj.bcs.getEntities(id));
                 ents = sum(obj.faces.faceNeighbors(faceID,:),2);
                 v(faceOrder,1) = obj.bcs.getVals(id,t);
                 switch obj.bcs.getType(id)
                    case 'Neu'
                       vals = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;
                       % area = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;
                       % vals = accumarray(ind, area);
                    case 'Dir'
                       % theta = obj.simParams.theta;
                       gamma = obj.material.getFluid().getFluidSpecWeight();
                       [mob, dmob] = obj.computeMobilityBoundary(obj.state.data.pressure(ents),v,faceID);
                       tr = obj.getFaceTransmissibilities(faceID);
                       press = obj.state.data.pressure(ents) - v;
                       gravT =  gamma*(obj.mesh.cellCentroid(ents,3) ...
                          - obj.faces.faceCentroid(faceID,3));
                       dirJ = mob.*tr;
                       q = dirJ.*(press+gravT);
                       if isNewtonNLSolver(obj.simParams)
                          % Contibution of Jh part in the boundary
                          dirJh = dmob.*tr;
                          dirJ = dirJ + dirJh.*(press+gravT);
                       end
                       vals = [dirJ,q]; % {JacobianVal,rhsVal]
                       % vals = [dirJ,accumarray(ind,q)]; % {JacobianVal,rhsVal]
                    case 'Spg'
                       gamma = obj.material.getFluid().getFluidSpecWeight();
                       assert(gamma>0.,'To impose Seepage boundary condition is necessary the fluid specify weight be bigger than zero!');

                       % theta = obj.simParams.theta;
                       zbc = obj.faces.faceCentroid(faceID,3);
                       href = v;                       
                       v = gamma*(href(1)-zbc);

                       v(v<=0)=0.;
                       [mob, dmob] = obj.computeMobilityBoundary(obj.state.data.pressure(ents),v,faceID);
                       tr = obj.getFaceTransmissibilities(faceID);
                       press = obj.state.data.pressure(ents) - v;
                       gravT = gamma*(obj.mesh.cellCentroid(ents,3) ...
                          - obj.faces.faceCentroid(faceID,3));
                       dirJ = mob.*tr;
                       q = dirJ.*(press+gravT);
                       if isNewtonNLSolver(obj.simParams)
                          % Contibution of Jh part in the boundary
                          dirJh = dmob.*tr;
                          dirJ = dirJ + dirJh.*(press+gravT);
                       end
                       vals = [dirJ,q]; % {JacobianVal,rhsVal]

                       % % pos=v>=0;
                       % % % v(v<=0)=0;  % Atmosferic pressure
                       % % % resize the number of boundary condition.
                       % % [ents,~,ind] = unique(ents(pos));
                       % % v=v(pos); faceID = faceID(pos);
                       % % [mob, dmob] = obj.computeMobilityBoundary(state.pressure(ents),v,faceID);
                       % % tr = obj.getFaceTransmissibilities(faceID);
                       % % press = state.pressure(ents) - v;
                       % % gravT = gamma*(obj.elements.cellCentroid(ents,3) ...
                       % %    - obj.faces.faceCentroid(faceID,3));
                       % % dirJ = mob.*tr;
                       % % q = dirJ.*(press+gravT);
                       % % if isNewtonNLSolver(obj.simParams)
                       % %    % Contibution of Jh part in the boundary                          
                       % %    dirJh = dmob.*tr.*press;
                       % %    dirJ = dirJ + dirJh;
                       % % end
                       % % vals = [dirJ,accumarray(ind,q)]; % {JacobianVal,rhsVal]
                       % % % vals = [dirJ,q]; % {JacobianVal,rhsVal]
                 end

              case 'VolumeForce'
                v = obj.bcs.getVals(id,t);
                ents = obj.bcs.getEntities(id);
                vals = v.*obj.mesh.cellVolume(ents);
           end
           % get local dof numbering
           dof = obj.dofm.getLocalDoF(ents,obj.fldId);
        end

        function flux = computeFluxBound(obj,flux,bound,pot,pres,t)
            %COMPUTEFLUX - compute the flux at the faces, than accumulate
            %the value at the nodes.
            % flux = zeros(obj.mesh.nNodes,3);
            if isFEMBased(obj.model,'Flow')
                % flux = fluidPot + gamma*obj.mesh.coordinates(:,3);
            elseif isFVTPFABased(obj.model,'Flow')
                nnodesBfaces = diff(obj.faces.mapN2F);
                neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
                Node2Face = repelem((1:obj.faces.nFaces)',nnodesBfaces);
                sgn = 2*(obj.faces.faceNeighbors(:,1)==0) - 1;

                areaSq = vecnorm(obj.faces.faceNormal,2,2);
                faceUnit = obj.faces.faceNormal./areaSq;
                areaSq = areaSq.*nnodesBfaces;

                % add boundary condition
                bcList = bound.db.keys;
                for bc = string(bcList)
                    field = translatePhysic(bound.getPhysics(bc),obj.model);
                    for f = field
                        if field == "VSFlow"
                            v = bound.getVals(bc,t);
                            switch bound.getCond(bc)
                                case {'NodeBC','ElementBC'}
                                case 'SurfBC'
                                    faceID = sort(bound.getEntities(bc));
                                    gamma = obj.material.getFluid().getFluidSpecWeight();
                                    switch bound.getType(bc)
                                        case 'Neu'
                                            vals = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;
                                        case 'Dir'
                                            ents = sum(obj.faces.faceNeighbors(faceID,:),2);
                                            potBd = pot(ents)-(v+gamma*obj.faces.faceCentroid(faceID,3));
                                            mob = computeMobilityBoundary(obj,pres(ents),v,faceID);
                                            vals = -mob.*obj.trans(faceID).*potBd;
                                        case 'Spg'
                                            ents = sum(obj.faces.faceNeighbors(faceID,:),2);
                                            % zbc = obj.faces.faceCentroid(faceID,3);
                                            % v= gamma*(v(1)-zbc);

                                            Datum = max(obj.mesh.coordinates);
                                            zbc = Datum(3)-obj.faces.faceCentroid(faceID,3);
                                            href = Datum(3)-v(1);
                                            v = gamma*(zbc-href(1));

                                            v(v<=0)=0.;                                            
                                            potBd = pot(ents)-(v+gamma*obj.faces.faceCentroid(faceID,3));
                                            mob = computeMobilityBoundary(obj,pres(ents),v,faceID);
                                            vals = mob.*obj.trans(faceID).*potBd;
                                            vals(:) = 0.; % after find the error, delete this line.
                                    end
                                    dir = sgn(faceID).*faceUnit(faceID,:);
                                    vals = vals./areaSq(faceID).*dir;
                                    vals = repelem(vals,nnodesBfaces(faceID),1);
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
    
    end

    methods (Access = private)
        function [Swkpt,dSwkpt,d2Swkpt] = computeSaturation(obj,pkpt)
            % COMPUTESATURATION compute the saturation and it's derivatives
            Swkpt = zeros(obj.mesh.nCells,1);
            dSwkpt = zeros(obj.mesh.nCells,1);
            d2Swkpt = zeros(obj.mesh.nCells,1);
            for m = 1:obj.mesh.nCellTag
                isElMat = obj.mesh.cellTag == m;
                p = pkpt(isElMat);
                Sws = obj.material.getMaterial(m).PorousRock.getMaxSaturation();
                Swr = obj.material.getMaterial(m).PorousRock.getResidualSaturation();
                [Swkpt(isElMat), dSwkpt(isElMat), d2Swkpt(isElMat)] = obj.material.getMaterial(m).Curves.computeSwAnddSw(p);
                Swkpt(isElMat) = Swr + (Sws-Swr)*Swkpt(isElMat);
                dSwkpt(isElMat) = (Sws-Swr)*dSwkpt(isElMat);
                d2Swkpt(isElMat) = (Sws-Swr)*d2Swkpt(isElMat);
            end
        end

        function [lwkpt,dlwkpt] = computeMobility(obj,pkpt)
            % COMPUTEMOBILITY compute the mobility and it's derivatives
            % for the upstream elements for each face
            nIntFaces = length(obj.upElem);
            lwkpt = zeros(nIntFaces,1);
            dlwkpt = zeros(nIntFaces,1);
            matUpElem = obj.mesh.cellTag(obj.upElem);
            for m = 1:obj.mesh.nCellTag
                isElMat = matUpElem == m;
                p = pkpt(obj.upElem(isElMat));
                [lwkpt(isElMat), dlwkpt(isElMat)] = obj.material.getMaterial(m).Curves.computeRelativePermeability(p);
                % [lwkpt(isElMat), dlwkpt(isElMat)] = obj.material.getMaterial(m).RelativePermCurve.interpTable(p);
            end
            mu = obj.material.getFluid().getDynViscosity();
            lwkpt = lwkpt/mu;
            dlwkpt = dlwkpt/mu;
        end

        function [lpt, dlpt] = computeMobilityBoundary(obj,pcells,pface,faceID)
            % COMPUTEMOBILITYBOUNDARY compute the mobility for the
            % upstream elements in the boundary
            elms = obj.faces.faceNeighbors(faceID,:);
            elms = elms(elms~=0);
            materials = obj.mesh.cellTag(elms);

            % Find the direction of the flux;
            gamma = obj.material.getFluid().getFluidSpecWeight();
            if gamma > 0
                zfaces = obj.faces.faceCentroid(faceID,3);
                cellz = obj.mesh.cellCentroid(elms,3);
                lElemIsUp = (pcells - pface) + gamma*(cellz- zfaces) >= 0;
            else
                lElemIsUp = pcells >= pface;
            end

            % Find the upstream pressure.
            pres = zeros(length(lElemIsUp),1);
            pres(lElemIsUp) = pcells(lElemIsUp);
            pres(~lElemIsUp) = pface(~lElemIsUp);

            mu = obj.material.getFluid().getDynViscosity();
            lpt = zeros(length(pcells),1);
            dlpt = zeros(length(pcells),1);
            for i=1:length(materials)
                % [krel,dkrel] = obj.material.computeRelativePermeability(pres(i),materials(i));
                [krel,dkrel] = obj.material.getMaterial(materials(i)).Curves.computeRelativePermeability(pres(i));
                lpt(i) = krel/mu;
                dlpt(i) = -dkrel/mu;
            end
            dlpt = dlpt.*lElemIsUp;
        end

        function [Swkpt,dSwkpt,d2Swkpt,lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt)
            % compute upstream elements for each face
            % interpolate effective saturation and relative permeability
            % and first derivatives
            % compute also second derivative for saturation
            neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
            gamma = obj.material.getFluid().getFluidSpecWeight();
            if gamma > 0
                zVec = obj.elements.mesh.cellCentroid(:,3);
                zNeigh = zVec(neigh);
                lElemIsUp = (pkpt(neigh(:,1)) - pkpt(neigh(:,2))) + gamma*(zNeigh(:,1) - zNeigh(:,2)) >= 0;
            else
                lElemIsUp = pkpt(neigh(:,1)) >= pkpt(neigh(:,2));
            end
            obj.upElem(lElemIsUp) = neigh(lElemIsUp,1);
            obj.upElem(~lElemIsUp) = neigh(~lElemIsUp,2);
            [Swkpt,dSwkpt,d2Swkpt] = computeSaturation(obj,pkpt);
            dSwkpt = - dSwkpt;
            [lwkpt,dlwkpt] = computeMobility(obj,pkpt);
            dlwkpt = - dlwkpt;
        end

        function Jh = computeJacobianPartJhNewton(obj,pTau,dMob)
          % COMPUTEJACOBIANFORNEWTONJFPART Method to compute the Jh part
          % of the jacobian used in the Newton-Raphson.
          %% Equation for the part.
          % $$ J_{h} = \theta p_{\tau}^{m}
          % \frac{\partial\lambda_{u}^{t+dt,m}}{\partial p_{t+dt}} \mathbf{T}$$
          % $$ J_{f} = \theta \gamma
          % \frac{\partial\lambda_{u}^{t+dt,m}}{\partial p_{t+dt}} \mathbf{T}$$
          % $$ J_{hf} = J_{h} + J_{f}$$

          % Define some values.
          % theta = obj.simParams.theta;
          gamma = obj.material.getFluid.getFluidSpecWeight();
          % subCells = obj.dofm.getFieldCells(obj.getField());

          % Get pairs of faces that contribute to the subdomain
          neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
          zVec = obj.mesh.cellCentroid(:,3);
          zNeigh = zVec(neigh);
          pTau = pTau(neigh);
          nneigh = length(dMob);
          % [~,~,reorder] = unique([neigh(:,1); neigh(:,2); obj.upElem'; subCells]);
          [~,~,reorder] = unique([neigh(:,1); neigh(:,2); obj.upElem']);
          neigh1 = reorder(1:nneigh);
          neigh2 = reorder(nneigh+1:2*nneigh);
          neighU = reorder(2*nneigh+1:3*nneigh);

          % Transmissibility of internal faces
          % Jh = theta*dMob.*obj.trans(obj.isIntFaces);
          Jh = dMob.*obj.trans(obj.isIntFaces);
          Jh = Jh.*(pTau(:,1) - pTau(:,2) + gamma*(zNeigh(:,1) - zNeigh(:,2)));

          % Computing the Jacobian Part.
          nDoF = obj.dofm.getNumDoF(obj.getField());
          Jh = sparse( [neigh1; neigh2], [neighU; neighU], ...
            [Jh; -Jh], nDoF, nDoF);
        end

        function Jp = computeJacobianPartJpNewton(obj,dt,pTmp,pOld,Stau,dStau,ddStau)
          % COMPUTEJACOBIANFORNEWTONPPART Method to compute the Jp part
          % of the jacobian used in the Newton-Raphson.
          %% Equation for the part.
          % $$ J_{p} = \left(\theta/\delta t\right) \left[
          % \alpha^{e}\beta S_{\tau}^{e,m}
          % + \left(2\alpha^{e}+\beta\phi_{\tau}^{e}\right)\frac{\partial S_{\tau}^{e,m}}{\partial p}
          % + \phi_{\tau}^{e} \frac{\partial^{2} S_{t+dt}^{m,e}}{\partial p_{t+dt}^{2}} \right]
          % \left(\frac{p_{t+dt}^{m}-p_{t}}{dt}\right)
          % \left|\Omega\right|^{e}$$

          % Define some constant.
          % theta = obj.simParams.theta;
          beta = obj.material.getFluid().getFluidCompressibility();
          subCells = obj.dofm.getFieldCells(obj.getField());
          nSubCells = length(subCells);
          poroMat = zeros(nSubCells,1);
          alphaMat = zeros(nSubCells,1);
          for m = 1:obj.mesh.nCellTag
            if ~ismember(m,obj.dofm.getFieldCellTags({obj.getField(),'Poromechanics'}))
              % compute alpha only if there's no coupling in the subdomain
              alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
            end
            poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
          end

          alphaMat=alphaMat(obj.mesh.cellTag(subCells));
          poroMat=poroMat(obj.mesh.cellTag(subCells));

          % Define some parameters.
          pTmp = pTmp(subCells);
          pOld = pOld(subCells);
          pdiff = pTmp-pOld;

          % Computing the Jacobian Part.
          % Jp = poroMat.*ddStau;
          % Jp = (theta/dt)*obj.elements.vol(subCells).*Jp.*pdiff;
          Jp = beta*alphaMat.*Stau + (2*alphaMat+beta*poroMat).*dStau + poroMat.*ddStau;
          Jp = (1./dt)*obj.elements.mesh.cellVolume(subCells).*Jp.*pdiff;

          nDoF = obj.dofm.getNumDoF(obj.getField());
          [~,~,dof] = unique(subCells);
          Jp = sparse(dof,dof,Jp,nDoF,nDoF);
        end
    end

    methods (Static)
        function [cellStr,pointStr] = buildPrintStruct(state)
            % pointStr(1).name = 'flux';
            % pointStr(1).data = state.flux;
            pointStr = [];
            cellStr = repmat(struct('name', 1, 'data', 1), 3, 1);
            cellStr(1).name = 'pressure';
            cellStr(1).data = state.pressure;
            cellStr(2).name = 'potential';
            cellStr(2).data = state.potential;
            cellStr(3).name = 'saturation';
            cellStr(3).data = state.saturation;
            cellStr(4).name = 'permeability';
            cellStr(4).data = state.perm;
            % cellStr(5).name = 'mass_cons';
            % cellStr(5).data = state.mass;
        end

        function out = getField()
          out = 'VariablySaturatedFlow';
        end
    end
end

