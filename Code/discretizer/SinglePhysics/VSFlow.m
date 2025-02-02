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
            pkpt = obj.simParams.theta*stateTmp.pressure + ...
                (1 - obj.simParams.theta)*statek.pressure;
            [Swkpt,dSwkpt,d2Swkpt, obj.lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt);
            computeStiffMatFV(obj,obj.lwkpt);
            computeCapMatFV(obj,Swkpt,dSwkpt);
            obj.J = obj.simParams.theta*obj.H + obj.P/dt;
            if isNewtonNLSolver(obj.simParams)
                computeNewtPartOfJacobian(obj,dt,statek,stateTmp,pkpt,dSwkpt,d2Swkpt,dlwkpt)
                obj.J = obj.J + obj.JNewt;
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

        function computeNewtPartOfJacobian(obj,dt,statek,stateTmp,pkpt,dSwkpt,d2Swkpt,dlwkpt)
            subCells = obj.dofm.getFieldCells(obj.field);
            nSubCells = length(subCells);
            % compute matrices J1 and J2 (gathering non linear terms)
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
            tmpVec2 = alphaMat(obj.mesh.cellTag) + beta*poroMat(obj.mesh.cellTag);
            tmpVec2 = 1/dt*((tmpVec2(subCells).*dSwkpt(subCells) + poroMat(obj.mesh.cellTag).*d2Swkpt(subCells)).*(stateTmp.pressure(subCells) - statek.pressure(subCells))).*obj.elements.vol(subCells);
            [~,~,neigh1] = unique(neigh(:,1));
            [~,~,neigh2] = unique(neigh(:,2));
            [~,~,upElemdof] = unique(obj.upElem);
            obj.JNewt = sparse([neigh1; neigh2; (1:nSubCells)'], ...
                [repmat(upElemdof,[2,1]);  (1:nSubCells)'], ...
                [tmpVec1; -tmpVec1; tmpVec2],nSubCells,nSubCells);
            obj.JNewt = obj.simParams.theta*obj.JNewt;
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
end

