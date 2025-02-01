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

        function computeMat(obj,stateTmp,statek,dt)
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


        function computeRhs(obj,stateTmp,statek,dt)
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

        function gTerm = finalizeRHSGravTerm(obj,lw)            
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'VSFlow'));
            nSubCells = nnz(obj.dofm.subCells(:,subInd));
            intFaces = any(obj.isIntFaces(:,subInd),2);
            neigh = obj.faces.faceNeighbors(intFaces,:);
            neigh = obj.dofm.ent2field('VSFlow',neigh(:));
            gTerm = accumarray(neigh,[lw.*obj.rhsGrav; ...
                 -lw.*obj.rhsGrav],[nSubCells,1]);
        end

        function blk = blockJacobian(obj,varargin)
            fRow = varargin{1}; 
            fCol = varargin{2};
            dt = varargin{3};
            locRow = obj.dofm.field2block(fRow);
            locCol = obj.dofm.field2block(fCol);
            if isNewtonNLSolver(obj.simParams)
                blk = obj.simParams.theta*obj.H(locRow,locCol) + obj.P(locRow,locCol)/dt + ...
                    obj.JNewt(locRow,locCol);
            else
                blk = obj.simParams.theta*obj.H(locRow,locCol) + obj.P(locRow,locCol)/dt;

            end
        end

        function blk = blockRhs(obj, fld)
            if ~strcmp(obj.dofm.subPhysics(fld), 'VSFlow')
                % no contribution to non poro fields
                blk = 0;
            else
                dofs = obj.dofm.field2block(fld);
                blk = obj.rhs(dofs);
            end
        end

        function out = isLinear(obj)
            out = false;
        end


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
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'VSFlow'));
            [subCells, ~] = find(obj.dofm.subCells(:,subInd));
            nSubCells = length(subCells);
            intFaces = any(obj.isIntFaces(:,subInd),2);
            % compute matrices J1 and J2 (gathering non linear terms)
            neigh = obj.faces.faceNeighbors(intFaces,:);
            zVec = obj.elements.cellCentroid(:,3);
            zNeigh = zVec(neigh);
            gamma = obj.material.getFluid().getFluidSpecWeight();
            tmpVec1 = (dlwkpt.*obj.trans(intFaces)).*(pkpt(neigh(:,1)) - pkpt(neigh(:,2)) + gamma*(zNeigh(:,1) - zNeigh(:,2)));
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
            neigh1dof = obj.dofm.ent2field('VSFlow',neigh(:,1));
            neigh2dof = obj.dofm.ent2field('VSFlow',neigh(:,2));
            upElemdof = obj.dofm.ent2field('VSFlow',obj.upElem);
            obj.JNewt = sparse([neigh1dof; neigh2dof; (1:nSubCells)'], ...
                [repmat(upElemdof,[2,1]);  (1:nSubCells)'], ...
                [tmpVec1; -tmpVec1; tmpVec2],nSubCells,nSubCells);
            obj.JNewt = obj.simParams.theta*obj.JNewt;
        end
    end
end

