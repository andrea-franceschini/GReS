function [J,rhs] = applyBCAndForces_MD(domID, mG, t, state, J, rhs)
% Apply Boundary conditions to multidomain global incremental
% linear system
% NOTE: BCs at the interface are applied only to master dofs
fields = find(~strcmp(extractfield(mG.MD_struct,'type'),'slave'));

% Get list of BCs for input domain
bc = mG.model(domID).BoundaryConditions;
keys = bc.db.keys;
model = mG.model(domID).ModelType;
% Recall: Dirichlet BCs add null rhs value to incremental system

for i = 1 : length(keys)
    % Get BC dofs and values
    dirVal = []; % if stays empty Penalty method is used
    rhsVal = []; % if stays empty Penalty method is used
    type = 'tmp';
    cond = bc.getCond(keys{i});
    ph = bc.getPhysics(keys{i});
    ph_mod = translatePhysic(ph,model);
    if ~strcmp(cond,'VolumeForce')
        type = bc.getType(keys{i});
    end
    switch cond
        case 'NodeBC'
            switch bc.getType(keys{i})
                case 'Neu'
                    rhsVal = - bc.getVals(keys{i}, t);
            end
        case 'SurfBC'
            if isFEMBased(model,ph)
                switch type
                    case 'Neu'
                        entitiesInfl = bc.getEntitiesInfluence(keys{i});
                        q = bc.getVals(keys{i}, t);
                        rhsVal = - entitiesInfl*q;
                end
            elseif isFVTPFABased(model,ph)
                faceID = bc.getEntities(keys{i},1);
                neigh = sum(grid.faces.faceNeighbors(faceID,:),2); % possibly repeated cells
                doftmp = bc.dof.getDoF(ph_mod,neigh);
                [bcDofs,~,ind] = unique(doftmp); % cell dof with no repetitions (corner cells have more than one face on the bcary!)
                switch bc.getType(keys{i})
                    case 'Neu'
                        area = -vecnorm(grid.faces.faceNormal(faceID,:),2,2).*bc.getVals(keys{i}, t);
                        rhsVal = accumarray(ind, area);
                    case 'Dir'
                        gamma = material.getFluid().getFluidSpecWeight();
                        mu = material.getFluid().getDynViscosity();
                        trans = getFaceTransmissibilities(syst.getField(ph),faceID);
                        q = 1/mu*trans.*((state.pressure(neigh) - bc.getVals(keys{i}, t))...
                            + gamma*(grid.cells.cellCentroid(neigh,3) - grid.faces.faceCentroid(faceID,3)));
                        rhsVal = accumarray(ind,q);
                        dirVal = 1/mu*trans;
                end
            end
        case 'VolumeForce'
            if isFEMBased(model, ph)
                entitiesInfl = bc.getEntitiesInfluence(keys{i});
                q = bc.getVals(keys{i}, t);
                rhsVal = - entitiesInfl*q;
            elseif isFVTPFABased(model, ph)
                rhsVal = - bc.getVals(keys{i}, t).*grid.cells.vol(bc.getEntities(keys{i}));
            end
    end

    bcDofs = [];
    for j = fields
        % get constrained DoF in global multidomain system
        if mG.MD_struct(j).dom == domID && strcmp(ph_mod,mG.MD_struct(j).physic)
        bcDofs = bc.getDofs_MD(keys{i},mG,j);
        bcDofs = bcDofs + mG.countDoF(j);
        end

        % ----------------------------- APPLY BC ----------------------------------
        if ~isempty(bcDofs) 
            switch type
                case 'Dir' % Dirichlet BC
                    if isempty(rhsVal)
                        vals = zeros(numel(bcDofs),1);
                        [J,rhs] = applyDir(bcDofs,vals,J,rhs);
                    end
                case 'Neu'
                    rhs(bcDofs) = rhs(bcDofs) + rhsVal;
                otherwise
                    rhs(bcDofs) = rhs(bcDofs) + rhsVal;
            end
        end
        bcDofs = [];
    end
end
end