classdef VSFlow < SPFlow
    % Variably Saturated flow
    % Subclass of SPFlow
    % Implements Richards equations for unsaturated flow in vadose region
    % It is subclass of SPFlow since most of the method are the same

    properties
        lwkpt           % mobility
        JNewt = []      % newton jacobian contribution
        upElem          % upstream elements array for each face
    end

    methods (Access = public)
        function obj = VSFlow(symmod,params,dofManager,grid,mat,data)
            % initialize as SPFlow class
            obj@SPFlow(symmod,params,dofManager,grid,mat,data,'VSFlow');
        end

        function stateTmp = computeMat(obj,stateTmp,statek,dt)
            theta = obj.simParams.theta;
            pkpt = theta*stateTmp.pressure + (1 - theta)*statek.pressure;
            % pkpt = obj.simParams.theta*stateTmp.pressure + ...
            %     (1 - obj.simParams.theta)*statek.pressure;
            [Swkpt,dSwkpt,d2Swkpt, obj.lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt);
            computeStiffMatFV(obj,obj.lwkpt);
            computeCapMatFV(obj,Swkpt,dSwkpt);
            obj.J = theta*obj.H + obj.P/dt;
            if isNewtonNLSolver(obj.simParams)
                Jh = computeJacobianPartJhNewton(obj,pkpt,dlwkpt);
                Jp = computeJacobianPartJpNewton(obj,dt,stateTmp.pressure,statek.pressure,Swkpt,dSwkpt,d2Swkpt);
                obj.JNewt = Jp+Jh;
                obj.J = obj.J + obj.JNewt;
                % JNewt2 = computeNewtPartOfJacobian(obj,dt,statek,stateTmp,pkpt,dSwkpt,d2Swkpt,dlwkpt);
                % [Jh,Jp] = computeNewtPartOfJacobian2(obj,dt,stateTmp.pressure,statek.pressure,dSwkpt,d2Swkpt,dlwkpt);
                % l1 = normest(Jh);
                % l2 = normest(Jp);
            end
        end

        function stateTmp = computeRhs(obj,stateTmp,statek,dt)
            % Compute the residual of the flow problem
            theta = obj.simParams.theta;
            ents = obj.dofm.getActiveEnts(obj.field);
            rhsStiff = theta*obj.H*stateTmp.pressure(ents) + (1-theta)*obj.H*statek.pressure(ents);
            rhsCap = (obj.P/dt)*(stateTmp.pressure(ents) - statek.pressure(ents));
            obj.rhs = rhsStiff + rhsCap;
            gamma = obj.material.getFluid().getFluidSpecWeight();
            %adding gravity rhs contribute
            if gamma > 0
                if isFEMBased(obj.model,'Flow')
                    obj.rhs = obj.rhs + obj.rhsGrav;
                elseif isFVTPFABased(obj.model,'Flow')
                    obj.rhs = obj.rhs + finalizeRHSGravTerm(obj,obj.lwkpt);
                end
            end
        end

        function state = setState(obj,state)
            n = obj.mesh.nCells;
            state.pressure = zeros(n,1);
            state.saturation = zeros(n,1);
        end

        function state = updateState(obj,state,dSol)
            ents = obj.dofm.getActiveEnts(obj.field);
            state.pressure(ents) = state.pressure(ents) + dSol(obj.dofm.getDoF(obj.field));
            state.saturation = obj.material.computeSwAnddSw(obj.mesh,state.pressure);
        end

        function [cellData,pointData] = printState(obj,sOld,sNew,t)
            % append state variable to output structure
            switch nargin
                case 2
                    fluidPot = finalizeState(obj,sOld);
                    pressure = sOld.pressure;
                case 4
                    % linearly interpolate state variables containing print time
                    fac = (t - sOld.t)/(sNew.t - sOld.t);
                    fluidPotOld = finalizeState(obj,sOld);
                    fluidPotNew = finalizeState(obj,sNew);
                    fluidPot = fluidPotNew*fac+fluidPotOld*(1-fac);
                    pressure = sNew.pressure*fac+sOld.pressure*(1-fac);
                otherwise
                    error('Wrong number of input arguments');
            end
            saturation = obj.material.computeSwAnddSw(obj.mesh,pressure);
            [cellData,pointData] = VSFlow.buildPrintStruct(pressure,fluidPot,saturation);
        end

        function out = isLinear(obj)
            out = false;
        end

        function [dof,vals] = getBC(obj,bc,id,t,state)
            switch bc.getCond(id)
                case {'NodeBC','ElementBC'}
                    ents = bc.getEntities(id);
                    vals = bound.getVals(id,t);
                case 'SurfBC'
                    if isFVTPFABased(obj.model,'Flow')
                        faceID = bc.getEntities(id);
                        ents = sum(obj.faces.faceNeighbors(faceID,:),2);
                        [ents,~,ind] = unique(ents);
                        switch bc.getType(id)
                            case 'Neu'
                                v = bc.getVals(id,t);
                                area = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;
                                vals = accumarray(ind, area);
                            case 'Dir'
                                v = bc.getVals(id,t);
                                gamma = obj.material.getFluid().getFluidSpecWeight();
                                [mob, dmob] = obj.computeMobilityBoundary(state.pressure(ents),v,faceID);
                                tr = obj.getFaceTransmissibilities(faceID);
                                potential = state.pressure(ents) - v ...
                                    + gamma*(obj.elements.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID,3));
                                dirJ = mob.*tr;
                                q = dirJ.*potential;
                                if isNewtonNLSolver(obj.simParams)
                                    % Contibution of Jh part in the boundary
                                    theta = obj.simParams.theta;
                                    dirJh = theta*dmob.*tr.*potential;
                                    dirJ = dirJ + dirJh;
                                end
                                vals = [dirJ,accumarray(ind,q)]; % {JacobianVal,rhsVal]
                            case 'Spg'
                                zbc = obj.elements.cellCentroid(ents,:);
                                href = bc.getVals(id,t);
                                v=zeros(length(zbc));
                                for jj=1:length(zbc)
                                    if zbc(:,3)<href(1,3)
                                        v(jj)=0.;
                                    else
                                        v(jj)=zbc(jj,3)-href(jj,3);
                                    end
                                end
                                gamma = obj.material.getFluid().getFluidSpecWeight();
                                mu = obj.material.getFluid().getDynViscosity();
                                tr = obj.getFaceTransmissibilities(faceID);
                                q = 1/mu*tr.*((state.pressure(ents) - v)...
                                    + gamma*(obj.elements.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID,3)));
                                vals = [1/mu*tr,accumarray(ind,q)]; % {JacobianVal,rhsVal]
                        end
                    elseif isFEMBased(obj.model,'Flow')
                        v = bc.getVals(id,t);
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
    end

    methods (Access = private)
        function [Swkpt,dSwkpt,d2Swkpt,lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt)
            % compute upstream elements for each face
            % interpolate effective saturation and relative permeability
            % and first derivatives
            % compute also second derivative for saturation
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
            [Swkpt,dSwkpt,d2Swkpt] = obj.material.computeSwAnddSw(obj.mesh,pkpt);
            dSwkpt = - dSwkpt;
            [lwkpt,dlwkpt] = obj.material.computeLwAnddLw(obj.mesh,obj.upElem,pkpt);
            dlwkpt = - dlwkpt;
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

                cellz = obj.elements.cellCentroid(elms,3);
                lElemIsUp = pcells - pface + gamma*(cellz- zfaces) >= 0;
            else
                lElemIsUp = pcells >= pface;
            end

            % Find the upstream pressure.
            pres = zeros(length(lElemIsUp),1);
            pres(lElemIsUp) = pcells(lElemIsUp);
            pres(~lElemIsUp) = pface(~lElemIsUp);
            % pres = pcells;
            mu = obj.material.getFluid().getDynViscosity();
            lpt = zeros(length(pcells),1);
            dlpt = zeros(length(pcells),1);
            for i=1:length(materials)
                [krel,dkrel] = obj.material.computeRelativePermeability(pres(i),materials(i));
                lpt(i) = krel/mu;
                dlpt(i) = -dkrel/mu;
            end
            dlpt = dlpt.*lElemIsUp;
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
            theta = obj.simParams.theta;
            gamma = obj.material.getFluid.getFluidSpecWeight();

            % subCells = obj.dofm.getFieldCells(obj.field);

            % Get pairs of faces that contribute to the subdomain
            neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
            zVec = obj.elements.cellCentroid(:,3);
            zNeigh = zVec(neigh);
            pTau = pTau(neigh);
            nneigh = length(dMob);
            % [~,~,reorder] = unique([neigh(:,1); neigh(:,2); obj.upElem'; subCells]);
            [~,~,reorder] = unique([neigh(:,1); neigh(:,2); obj.upElem']);
            neigh1 = reorder(1:nneigh);
            neigh2 = reorder(nneigh+1:2*nneigh);
            neighU = reorder(2*nneigh+1:3*nneigh);

            % Transmissibility of internal faces            
            Jh = theta*dMob.*obj.trans(obj.isIntFaces);
            Jh = Jh.*(pTau(:,1) - pTau(:,2) + gamma*(zNeigh(:,1) - zNeigh(:,2)));
            
            % Computing the Jacobian Part.
            nDoF = obj.dofm.getNumDoF(obj.field);
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
            theta = obj.simParams.theta;
            beta = obj.material.getFluid().getFluidCompressibility();
            subCells = obj.dofm.getFieldCells(obj.field);
            nSubCells = length(subCells);
            poroMat = zeros(nSubCells,1);
            alphaMat = zeros(nSubCells,1);            
            for m = 1:obj.mesh.nCellTag
                if ~ismember(m,obj.dofm.getFieldCellTags({obj.field,'Poromechanics'}))
                    % compute alpha only if there's no coupling in the
                    % subdomain
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
            Jp = beta*alphaMat.*Stau + (2*alphaMat+beta*poroMat).*dStau + poroMat.*ddStau;
            Jp = (theta/dt)*obj.elements.vol(subCells).*Jp.*pdiff;
            nDoF = obj.dofm.getNumDoF(obj.field);
            [~,~,dof] = unique(subCells);
            Jp = sparse(dof,dof,Jp,nDoF,nDoF);
        end
    end

    methods (Static)
        function [cellStr,pointStr] = buildPrintStruct(press,pot,sat)
            pointStr = [];
            cellStr = repmat(struct('name', 1, 'data', 1), 3, 1);
            cellStr(1).name = 'pressure';
            cellStr(1).data = press;
            cellStr(2).name = 'potential';
            cellStr(2).data = pot;
            cellStr(3).name = 'saturation';
            cellStr(3).data = sat;
        end
    end

    methods (Access = private)
        % function [Swkpt,dSwkpt,d2Swkpt] = computeSaturation(obj,pkpt)
        %     % COMPUTESATURATION compute the effective saturation and it's
        %     % derivatives for the upstream elements for each face
        %     neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
        %     gamma = obj.material.getFluid().getFluidSpecWeight();
        %     if gamma > 0
        %         zVec = obj.elements.cellCentroid(:,3);
        %         zNeigh = zVec(neigh);
        %         lElemIsUp = pkpt(neigh(:,1)) - pkpt(neigh(:,2)) + gamma*(zNeigh(:,1) - zNeigh(:,2)) >= 0;
        %     else
        %         lElemIsUp = pkpt(neigh(:,1)) >= pkpt(neigh(:,2));
        %     end
        %     obj.upElem(lElemIsUp) = neigh(lElemIsUp,1);
        %     obj.upElem(~lElemIsUp) = neigh(~lElemIsUp,2);
        %     [Swkpt,dSwkpt,d2Swkpt] = obj.material.computeSwAnddSw(obj.mesh,pkpt);
        %     dSwkpt = - dSwkpt;
        % end

        % function [lpt,dlpt] = computeMobility(obj,pkpt)
        %     % COMPUTEMOBILITY compute the mobility and it's derivatives
        %     % for the upstream elements for each face
        %     neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
        %     gamma = obj.material.getFluid().getFluidSpecWeight();
        %     if gamma > 0
        %         zVec = obj.elements.cellCentroid(:,3);
        %         zNeigh = zVec(neigh);
        %         lElemIsUp = pkpt(neigh(:,1)) - pkpt(neigh(:,2)) + gamma*(zNeigh(:,1) - zNeigh(:,2)) >= 0;
        %     else
        %         lElemIsUp = pkpt(neigh(:,1)) >= pkpt(neigh(:,2));
        %     end
        %     obj.upElem(lElemIsUp) = neigh(lElemIsUp,1);
        %     obj.upElem(~lElemIsUp) = neigh(~lElemIsUp,2);
        %     [lpt,dlpt] = obj.material.computeLwAnddLw(obj.mesh,obj.upElem,pkpt);
        %     dlpt = - dlpt;
        % end

        % function JNewt = computeNewtPartOfJacobian(obj,dt,statek,stateTmp,pkpt,dSwkpt,d2Swkpt,dlwkpt)
        %     subCells = obj.dofm.getFieldCells(obj.field);
        %     nSubCells = length(subCells);
        %     % compute matrices J1 and J2 (gathering non linear terms)
        %     neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
        %     zVec = obj.elements.cellCentroid(:,3);
        %     zNeigh = zVec(neigh);
        %     gamma = obj.material.getFluid().getFluidSpecWeight();
        %     tmpVec1 = (dlwkpt.*obj.trans(obj.isIntFaces)).*(pkpt(neigh(:,1)) - pkpt(neigh(:,2)) + gamma*(zNeigh(:,1) - zNeigh(:,2)));
        %     %
        %     poroMat = zeros(obj.mesh.nCellTag,1);
        %     alphaMat = zeros(obj.mesh.nCellTag,1);
        %     beta = obj.material.getFluid().getFluidCompressibility();
        %     for m = 1:obj.mesh.nCellTag
        %         poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
        %         alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
        %     end
        %     tmpVec2 = alphaMat(obj.mesh.cellTag) + beta*poroMat(obj.mesh.cellTag);
        %     tmpVec2 = 1/dt*((tmpVec2(subCells).*dSwkpt(subCells) + poroMat(obj.mesh.cellTag).*d2Swkpt(subCells)).*(stateTmp.pressure(subCells) - statek.pressure(subCells))).*obj.elements.vol(subCells);
        %     [~,~,reordercells] = unique([neigh(:,1);neigh(:,2);obj.upElem';subCells]);
        % 
        %     nneigh = length(tmpVec1);
        %     neigh1 = reordercells(1:nneigh);
        %     neigh2 = reordercells(nneigh+1:2*nneigh);
        %     neighU = reordercells(2*nneigh+1:3*nneigh);
        %     JNewt = sparse([neigh1; neigh2; SubCells], ...
        %         [neighU; neighU; SubCells], ...
        %         [tmpVec1; -tmpVec1; tmpVec2],nSubCells,nSubCells);
        % 
        %     % [~,~,neigh1] = unique(neigh(:,1));
        %     % [~,~,neigh2] = unique(neigh(:,2));
        %     % [~,~,upElemdof] = unique(obj.upElem);
        %     % JNewt = sparse([neigh1; neigh2; (1:nSubCells)'], ...
        %     %     [repmat(upElemdof,[2,1]);  (1:nSubCells)'], ...
        %     %     [tmpVec1; -tmpVec1; tmpVec2],nSubCells,nSubCells);
        % 
        %     JNewt = obj.simParams.theta*JNewt;
        % end

        % function [Jh, Jp] = computeNewtPartOfJacobian2(obj,dt,pTmp,pOld,dSwkpt,d2Swkpt,dlwkpt)
        %     % Initial considerations.
        %     subCells = obj.dofm.getFieldCells(obj.field);
        %     nSubCells = length(subCells);
        %     cellTag = obj.mesh.cellTag(subCells);
        %     nCellTag = obj.mesh.nCellTag;
        %     neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
        % 
        %     theta = obj.simParams.theta;
        %     gamma = obj.material.getFluid().getFluidSpecWeight();
        %     beta = obj.material.getFluid().getFluidCompressibility();
        % 
        %     pTau = theta*pTmp+(1-theta)*pOld;
        %     pTmp = pTmp(subCells);
        %     pOld = pOld(subCells);
        % 
        %     % compute matrices J1 and J2 (gathering non linear terms)
        %     % Computing Jh term
        %     zVec = obj.elements.cellCentroid(:,3);
        %     zNeigh = zVec(neigh);
        %     pTau = pTau(neigh);
        %     Jh = theta*dlwkpt.*obj.trans(obj.isIntFaces);
        %     Jh = Jh.*(pTau(:,1) - pTau(:,2) + gamma*(zNeigh(:,1) - zNeigh(:,2)));
        % 
        %     % Computing Jp term
        %     poroMat = zeros(nCellTag,1);
        %     alphaMat = zeros(nCellTag,1);            
        %     for m = 1:nCellTag
        %         poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
        %         alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
        %     end
        %     Jp = alphaMat(cellTag) + beta*poroMat(cellTag);
        %     Jp = 1/dt*((Jp(subCells).*dSwkpt(subCells) + poroMat(cellTag).*d2Swkpt(subCells)).*(pTmp - pOld)).*obj.elements.vol(subCells);
        % 
        %     % Making the matrix's.
        %     [~,~,neigh] = unique([neigh(:,1);neigh(:,2)]);
        %     sumDiagTrans = accumarray(neigh,repmat(Jh,[2,1]),[nSubCells,1]);
        %     Jh = sparse([neigh; subCells], ...
        %         [neigh; subCells], ...
        %         [-Jh;-Jh;sumDiagTrans], nSubCells, nSubCells);
        % 
        %     Jp = sparse(subCells,subCells,Jp,nSubCells,nSubCells);
        % end
    end
end

