function old_ApplyBCAndForces(model, grid, bound, material, t, syst, state)
% old version of the function for the application of boundary conditions
% this applied bcs to the glboal system, instead of doing it block wise
keys = bound.db.keys;
maxVal = max(syst.J, [], "all");
for i = 1:length(keys)
    % --------------------------NODE BC ----------------------------------- 
    if strcmp(bound.getCond(keys{i}),'NodeBC') 
      if strcmp(bound.getType(keys{i}), 'Neu')  % Apply Neumann conditions,if any
        syst.rhs(bound.getDofs(keys{i})) = syst.rhs(bound.getDofs(keys{i})) - bound.getVals(keys{i}, t);
      elseif strcmp(bound.getType(keys{i}), 'Dir')  % Apply Dirichlet conditions
        nrows = size(syst.J,1);
        syst.J(nrows*(bound.getDofs(keys{i})-1) + bound.getDofs(keys{i})) = maxVal*1.e10;
        syst.rhs(bound.getDofs(keys{i})) = 0;
      end
    % --------------------------ELEMENT BC ---------------------------------   
    elseif strcmp(bound.getCond(keys{i}),'ElementBC') %apply dirichlet condition assigned on cells (FVTPFA scheme)
      nrows = size(syst.J,1);
      syst.J(nrows*(bound.getDofs(keys{i})-1) + bound.getDofs(keys{i})) = maxVal*1.e10;
      syst.rhs(bound.getDofs(keys{i})) = 0;
    % --------------------------SURFACE BC --------------------------------- 
    elseif strcmp(bound.getCond(keys{i}),'SurfBC')
      %FEM surface BC  
      if isFEMBased(model,bound.getPhysics(keys{i}))
        if strcmp(bound.getType(keys{i}), 'Neu')
          entitiesInfl = bound.getEntitiesInfluence(keys{i});
          q = bound.getVals(keys{i}, t);
          syst.rhs(bound.getDofs(keys{i})) = syst.rhs(bound.getDofs(keys{i})) - entitiesInfl*q;
        elseif strcmp(bound.getType(keys{i}), 'Dir')
            nrows = size(syst.J,1);
            syst.J(nrows*(bound.getDofs(keys{i})-1) + bound.getDofs(keys{i})) = maxVal*1.e10;
            syst.rhs(bound.getDofs(keys{i})) = 0;
        end
      elseif isFVTPFABased(model,bound.getPhysics(keys{i}))
        if strcmp(bound.getType(keys{i}), 'Neu')
          faceID = bound.getEntities(keys{i});
%           neigh = grid.faces.faceNeighbors(faceID,:);
%           assert(all(sum(neigh~=0,2) == 1),'Corrupted face numbering in %s',bound.getName(keys{i}));
          neigh = sum(grid.faces.faceNeighbors(faceID,:),2);
          syst.rhs(neigh) = syst.rhs(neigh) - vecnorm(grid.faces.faceNormal(faceID,:),2,2).*bound.getVals(keys{i}, t);
        elseif strcmp(bound.getType(keys{i}), 'Dir')
%           error('Dirichlet cond. for FVTPFA on surf has not been implemented yet (ref. %s)', ...
%           bound.getName(keys{i}));
          nrows = size(syst.J,1);
          faceID = bound.getEntities(keys{i});
          neighEl = sum(grid.faces.faceNeighbors(faceID,:),2);
          neighElLoc = bound.dof.getLocDoF('Flow',neighEl);
          neighEl = bound.dof.getDoF('Flow',neighEl);
          gamma = material.getFluid().getFluidSpecWeight();
          mu = material.getFluid().getDynViscosity();
          trans = getFaceTransmissibilities(syst,faceID);
          q = 1/mu*trans.*((state.pressure(neighElLoc) - bound.getVals(keys{i}, t))...
            + gamma*(grid.cells.cellCentroid(neighElLoc,3) - grid.faces.faceCentroid(faceID,3)));
          syst.rhs(neighEl) = syst.rhs(neighEl) + q;
          syst.J(nrows*(neighEl-1) + neighEl) = syst.J(nrows*(neighEl-1) + neighEl) + 1/mu*trans;
        end
      end
    % --------------------------VOLUME FORCE ---------------------------------   
    elseif strcmp(bound.getCond(keys{i}), 'VolumeForce')
      if isFEMBased(model,bound.getPhysics(keys{i}))
        entitiesInfl = bound.getEntitiesInfluence(keys{i});
        q = bound.getVals(keys{i}, t);
        entitiesForce = entitiesInfl*q;
        syst.rhs(bound.getDofs(keys{i})) = syst.rhs(bound.getDofs(keys{i})) - entitiesForce;
      elseif isFVTPFABased(model,bound.getPhysics(keys{i}))
        syst.rhs(bound.getDofs(keys{i})) = syst.rhs(bound.getDofs(keys{i})) - bound.getVals(keys{i}, t).*grid.cells.vol(bound.getEntities(keys{i}));
      end
    end
end
end

