function applyBCandForces_test(model, grid, bound, material, t, syst, state)
  % Apply Boundary condition to the block partitioned system. The local 
  % dof indexing is used. BC assignment is done by ApplyBC function.
  % Impose BC to the linearized system (Jacobian matrix + RHS)
  % The Penalty method is used for the Dirichlet conditions
  keys = bound.db.keys;
  for i = 1 : length(keys)
      dirVal = []; % if stays empty Penalty method is used
      rhsVal = []; % if stays empty Penalty method is used
      type = 'tmp';
      cond = bound.getCond(keys{i});
      if ~strcmp(cond,"VolumeForce")
        type = bound.getType(keys{i});
      end
      ph = bound.getPhysics(keys{i});
      switch cond
          case 'NodeBC'
              bcDofs = bound.getDofs(keys{i});
              switch bound.getType(keys{i})
                  case 'Neu'
                      rhsVal = - bound.getVals(keys{i}, t);
              end
          case 'ElementBC'
              bcDofs = bound.getDofs(keys{i});
          case 'SurfBC'
              if isFEMBased(model,bound.getPhysics(keys{i}))
                  bcDofs = bound.getDofs(keys{i});
                  switch type
                      case 'Neu'
                          entitiesInfl = bound.getEntitiesInfluence(keys{i});
                          q = bound.getVals(keys{i}, t);
                          rhsVal = - entitiesInfl*q;
                  end    
              elseif isFVTPFABased(model,bound.getPhysics(keys{i}))
                  faceID = bound.getEntities(keys{i});
                  neigh = sum(grid.faces.faceNeighbors(faceID,:),2); % possibly repeated cells
                  doftmp = bound.dof.getDoF('Flow',neigh); 
                  [bcDofs,~,ind] = unique(doftmp); % cell dof with no repetitions (corner cells have more than one face on the boundary!)
                  switch bound.getType(keys{i})                     
                      case 'Neu'
                          area = -vecnorm(grid.faces.faceNormal(faceID,:),2,2).*bound.getVals(keys{i}, t);
                          rhsVal = accumarray(ind, area);
                      case 'Dir'
                          gamma = material.getFluid().getFluidSpecWeight();
                          mu = material.getFluid().getDynViscosity();
                          trans = getFaceTransmissibilities(syst,faceID);
                          q = 1/mu*trans.*((state.pressure(neigh) - bound.getVals(keys{i}, t))...
                            + gamma*(grid.cells.cellCentroid(neigh,3) - grid.faces.faceCentroid(faceID,3)));
                          rhsVal = accumrray(ind,q);
                          dirVal = 1/mu*trans;
                  end    
              end
          case 'VolumeForce'
              bcDofs = bound.getDofs(keys{i});
              if isFEMBased(model,bound.getPhysics(keys{i}))
                entitiesInfl = bound.getEntitiesInfluence(keys{i});
                q = bound.getVals(keys{i}, t);
                rhsVal = - entitiesInfl*q;
              elseif isFVTPFABased(model,bound.getPhysics(keys{i}))
                rhsVal = - bound.getVals(keys{i}, t).*grid.cells.vol(bound.getEntities(keys{i}));
              end
      end

% ----------------------------- APPLY BC ----------------------------------
    % get entity type associated to BC
    if isFEMBased(model, ph)
        ent = 'node';
    elseif isFVTPFABased(model, ph)
        ent = 'elem';
    end

    
    locDofs = bound.dof.glob2loc(bcDofs,ent); % get local dofs 
    phDofs = bound.dof.glob2sub(bcDofs, ent);  % get block ID associated to given bc global dofs
    switch type
        case 'Dir' % Dirichlet BC
            for j = unique(phDofs)'
                dof = locDofs(phDofs == j);
                nrows = size(syst.blockJ(j,j).block,1);
                if isempty(rhsVal) && isempty(dirVal) % penalty method
                    maxVal = max(syst.blockJ(j,j).block, [], "all");
                    syst.blockJ(j,j).block(nrows*(dof-1) + dof) = maxVal*1.e10;
                    syst.blockRhs(j).block(dof) = 0;
                else
                    syst.blockJ(j,j).block(nrows*(dof-1) + dof) = syst.blockJ(j,j).block(nrows*(dof-1) + dof) + dirVal;
                    syst.blockRhs(j).block(dof) = syst.blockRhs(j).block(dof) + rhsVal;
                end
            end
        case 'Neu'
            for j = unique(phDofs)'
                dof = locDofs(phDofs == j);
                syst.blockRhs(j).block(dof) = syst.blockRhs(j).block(dof) + rhsVal;
            end
        otherwise
            for j = unique(phDofs)'
                dof = locDofs(phDofs == j);
                syst.blockRhs(j).block(dof) = syst.blockRhs(j).block(dof) + rhsVal;
            end

    end
  end

end