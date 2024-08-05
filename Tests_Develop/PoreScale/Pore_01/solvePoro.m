function solvePoro(model,interfF2S,interfS2F,K,poro,alpha)
% Sequential algorithm to solve the Pore scale problem
% 1) solve flow
% LOOP UNTIL REL.PRESSURE CHANGE < TOL
% 2) use mortar to transfer surface pressures
% 3) solve mechanics
% 4) update K
% 5) solve pressure

tol = 1e-6;
itMax = 5;

% this function works with 2 domains
% domain 1: flow matrix
% domain 2: grains
t = 0;
tStep = 0;
dt = model(1).SimParams.dtIni;
% create statek object for flow domain
statFlow = copy(model(1).State);
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
    resP = norm(model(1).Discretizer.rhs{1});
    % solve system
    dp = model(1).Discretizer.solve();
    model(1).Discretizer.resetJacobianAndRhs();
    % Update tmpState
    model(1).State.updateState(dp,model(1).DoFManager);
    iter = 0;
    while resP > tol && iter < itMax
        iter = iter+1;
        % Mechanical domain
        %statPoro = copy(model(2).State);
        % compute mechanics matrix
        model(2).Discretizer.getField('Poro').computeMat(model(2).State,dt);
        % apply dirichlet values to incremental system
        applyDirVal(model(2).ModelType,model(2).BoundaryConditions,0,...
            model(2).State);
        % compute Rhs
        model(2).Discretizer.getField('Poro').computeRhs(model(2).State);
        % Retrieve block jacobian and rhs
        model(2).Discretizer.computeBlockJacobianAndRhs(dt);

        % apply BCs
        s2 = struct2cell(model(2));
        applyBCandForces(s2{[3,5,8,6]},0,s2{[11,9]})

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
        du = model(2).Discretizer.solve();
        model(2).Discretizer.resetJacobianAndRhs();
        % Update tmpState
        model(2).State.updateState(du,model(2).DoFManager);
        % Update coordinates of grains
        du_mat = reshape(du,3,[])';
        for i = 1:numel(interfF2S)
            interfF2S(i).mortar.intSlave.coordinates = ...
                interfF2S(i).mortar.intSlave.coordinates + du_mat(interfF2S(i).slaveSet,:);
        end
        Kold = K;
        K = zeros(model(1).Grid.topology.nCells,1);
        for i = 1:model(1).Grid.topology.nCells
            % compute channel size
            centroid = model(1).Grid.cells.cellCentroid(i,:);
            d(i) = computeChannelSize(centroid,interfF2S(1).mortar.intSlave,interfF2S(2).mortar.intSlave);
            K(i) = ((d(i))^3/12);
        end
        normK = norm((K-Kold)./Kold);
        % Solve pressure equation
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
        % solve system
        dp = model(1).Discretizer.solve();
        model(1).Discretizer.resetJacobianAndRhs();
        % Update tmpState
        model(1).State.updateState(dp,model(1).DoFManager);
    end
    % print results
    printState(model(1).OutState,model(1).State);
    printState(model(2).OutState,model(2).State);
end
end

