function [rhsP,rhsU] = solvePoro(model,interfF2S,interfS2F,K,poro,alpha)
% Sequential algorithm to solve the Pore scale problem
% 1) solve flow
% LOOP UNTIL REL.PRESSURE CHANGE < TOL
% 2) use mortar to transfer pressures to mechanics nodes
% 3) solve mechanics
% 4) use mortar to compute displacement in the flow grid to update channel
% aperture
% 4) update K
% 5) solve pressure

tol = 1e-12;
itMax = 6;

% this function works with 2 domains
% domain 1: flow matrix
% domain 2: grains
t = 0;
tStep = 0;
dt = model(1).SimParams.dtIni;
% create statek object for flow domain
statFlow = copy(model(1).State);
rhsP = zeros(itMax,1);
rhsU = zeros(itMax,1);
rhsP2 = zeros(itMax,1);
rhsU2 = zeros(itMax,1);
while t < model(1).SimParams.tMax
    t = t + dt;
    tStep = tStep + 1;
    % compute flow matrix with optional parameters
    model(1).Discretizer.getField('Flow').computeMat(K,poro,alpha);
    % apply dirichlet values to incremental system
    applyDirVal(model(1).ModelType,model(1).BoundaryConditions,0,...
        model(1).State);
    % compute Rhs
    model(1).Discretizer.getField('Flow').computeRhs(model(1).State,statFlow,dt);
    % Retrieve block jacobian and rhs
    model(1).Discretizer.computeBlockJacobianAndRhs(dt);
    % apply BCs
    s1 = struct2cell(model(1));
    applyBCandForces(s1{[3,5,8,6]},0,s1{[11,9]})
    rhsNormP = norm(model(1).Discretizer.rhs{1});
    % solve system
    dp = model(1).Discretizer.solve();
    model(1).Discretizer.resetJacobianAndRhs();
    % Update tmpState
    model(1).State.updateState(dp,model(1).DoFManager);
    pnew = model(1).State.pressure;
    iter = 0;
    %rhsNormU = 2*tol;
    rhsP(iter+1) = rhsNormP;
    resU = 2*tol;
    resP = 2*tol;
    while (iter < itMax) && (resP > tol && resU > tol)
        iter = iter+1;
        % Mechanical domain
        %statPoro = copy(model(2).State);
        % compute mechanics matrix
        if isempty(model(2).Discretizer.getField('Poro').K)
        model(2).Discretizer.getField('Poro').computeMat(model(2).State,dt);
        end
        % apply dirichlet values to incremental system
        applyDirVal(model(2).ModelType,model(2).BoundaryConditions,0,...
            model(2).State);
        % compute Rhs
        model(2).Discretizer.getField('Poro').computeRhs(model(2).State);
        % Retrieve block jacobian and rhs
        model(2).Discretizer.computeBlockJacobianAndRhs(dt);

        % employ interpolation operators to transfer pressure to the grains
        % compute forcing term due to fluid pressure
        for i = 1:numel(interfF2S)
            Ep = interfF2S(i).InterpOperator;
            %Eu = interfaces(i).InterpS2M;
            p = model(1).State.pressure(interfF2S(i).masterSet);
            % interpolation
            pInt = Ep*p;
            % compute nodal forces
            F = (interfF2S(i).nodeNormal.*pInt)';
            dofs = model(2).DoFManager.getDoF('Poro',interfF2S(i).slaveSet);
            dofs = glob2blockDoF(model(2).DoFManager,dofs);
            % update rhs with forcing terms
            model(2).Discretizer.rhs{1}(dofs) = model(2).Discretizer.rhs{1}(dofs) + F(:);
        end

        % apply BCs
        s2 = struct2cell(model(2));
        applyBCandForces(s2{[3,5,8,6]},0,s2{[11,9]})
        
        rhsNormU = norm(model(2).Discretizer.rhs{1});
        rhsU(iter) = rhsNormU;

        du = model(2).Discretizer.solve();
        model(2).Discretizer.resetJacobianAndRhs();
        % Update tmpState
        model(2).State.updateState(du,model(2).DoFManager);
        % Interpolate displacements on the interface 
        % and update coordinates of flow mesh to update channel size
        %du_mat = reshape(du,3,[])';
        for i = 1:numel(interfS2F)
            Eu = interfS2F(i).InterpOperator;
            u = du(model(2).DoFManager.getDoF('Poro',interfS2F(i).masterSet));
            for j = 1:3
            % interpolate each coordinate separately
            uInt = Eu*u(j:3:end);
            interfF2S(i).mortar.intMaster.coordinates(:,j) = ...
                interfF2S(i).mortar.intMaster.coordinates(:,j) + uInt;
            end
        end
        Kold = K;
        d = zeros(model(1).Grid.topology.nCells,1);
        K = zeros(model(1).Grid.topology.nCells,1);
        for i = 1:model(1).Grid.topology.nCells
            % compute channel size
            centroid = model(1).Grid.cells.cellCentroid(i,:);
            d(i) = computeChannelSize(centroid,interfF2S(1).mortar.intMaster,interfF2S(2).mortar.intMaster);
            K(i) = ((d(i))^3/12);
        end
        normK = (K-Kold)./Kold;
        % Solve pressure equation
        model(1).Discretizer.getField('Flow').computeMat(K,poro,alpha);
        % compute Rhs
        model(1).Discretizer.getField('Flow').computeRhs(model(1).State,statFlow,dt);
        % Retrieve block jacobian and rhs
        model(1).Discretizer.computeBlockJacobianAndRhs(dt);
        % apply BCs
        s1 = struct2cell(model(1));
        applyBCandForces(s1{[3,5,8,6]},0,s1{[11,9]})
        rhsNormP = norm(model(1).Discretizer.rhs{1});
        rhsP(iter+1) = rhsNormP;
        % solve system
        dp = model(1).Discretizer.solve();
        model(1).Discretizer.resetJacobianAndRhs();
        % Update tmpState
        model(1).State.updateState(dp,model(1).DoFManager);
        pold = pnew;
        pnew = model(1).State.pressure;

        % control convergence of pressure solution
        resP = norm(pnew-pold)/norm(pold);
        if iter > 1
            uold = unew;
            unew = model(2).State.dispCurr;
            resU = norm(unew-uold)/norm(uold);
        else
            unew = model(2).State.dispCurr;
        end
        rhsP2(iter+1) = resP;
        rhsU2(iter+1) = resU;
    end
    model(2).State.dispConv = model(2).State.dispCurr;
    % print results
    printState(model(1).OutState,model(1).State);
    printState(model(2).OutState,model(2).State);
end
end

