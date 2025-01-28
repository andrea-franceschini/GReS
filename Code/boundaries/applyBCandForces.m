function applyBCandForces(model, grid, bound, material, t, syst, state)
% Apply Boundary condition to linear system
% BCs is imposed to each physics solver separately
keys = bound.db.keys;
for i = 1 : length(keys)
   dirVal = [];
   rhsVal = [];
   type = 'tmp';
   cond = bound.getCond(keys{i});
   ph = bound.getPhysics(keys{i});
   ph_mod = translatePhysic(ph,model);
   if ~strcmp(cond,'VolumeForce')
      type = bound.getType(keys{i});
   end
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
         if isFEMBased(model,ph)
            bcDofs = bound.getDofs(keys{i});
            switch type
               case 'Neu'
                  entitiesInfl = bound.getEntitiesInfluence(keys{i});
                  q = bound.getVals(keys{i}, t);
                  rhsVal = - entitiesInfl*q;
            end
         elseif isFVTPFABased(model,ph)
            faceID = bound.getEntities(keys{i});
            neigh = sum(grid.faces.faceNeighbors(faceID,:),2); % possibly repeated cells
            doftmp = bound.dof.getDoF(ph_mod,neigh);
            [bcDofs,~,ind] = unique(doftmp); % cell dof with no repetitions (corner cells have more than one face on the boundary!)
            switch bound.getType(keys{i})
               case 'Neu'
                  area = -vecnorm(grid.faces.faceNormal(faceID,:),2,2).*bound.getVals(keys{i}, t);
                  rhsVal = accumarray(ind, area);
               case 'Dir'
                  gamma = material.getFluid().getFluidSpecWeight();
                  mu = material.getFluid().getDynViscosity();
                  trans = getFaceTransmissibilities(syst.getField(ph),faceID);
                  q = 1/mu*trans.*((state.pressure(neigh) - bound.getVals(keys{i}, t))...
                     + gamma*(grid.cells.cellCentroid(neigh,3) - grid.faces.faceCentroid(faceID,3)));
                  rhsVal = accumarray(ind,q);
                  dirVal = 1/mu*trans;
            end
         end
      case 'VolumeForce'
         bcDofs = bound.getDofs(keys{i});
         if isFEMBased(model, ph)
            entitiesInfl = bound.getEntitiesInfluence(keys{i});
            q = bound.getVals(keys{i}, t);
            rhsVal = - entitiesInfl*q;
         elseif isFVTPFABased(model, ph)
            rhsVal = - bound.getVals(keys{i}, t).*grid.cells.vol(bound.getEntities(keys{i}));
         end
   end

   % ----------------------------- APPLY BC ----------------------------------

   locDofs = bound.dof.glob2blockDoF(bcDofs); % get local block dofs
   phDofs = bound.dof.glob2block(bcDofs);  % get block ID associated to given bc global dofs
   switch type
      case 'Dir' % Dirichlet BC
         for j = unique(phDofs)'
            dof = locDofs(phDofs == j);
            nrows = size(syst.J{j,j},1);
            if isempty(rhsVal) && isempty(dirVal)
               %                     maxVal = max(syst.J{j,j}, [], "all");
               %                     syst.J{j,j}(nrows*(dof-1) + dof) = maxVal*1.e10;
               %                     syst.rhs{j}(dof) = 0;
               vals = zeros(numel(dof),1);
               [syst.J,syst.rhs] = applyDirBlock(dof,vals,syst.J,syst.rhs,j);
            else
               syst.J{j,j}(nrows*(dof-1) + dof) = syst.J{j,j}(nrows*(dof-1) + dof) + dirVal(ind);
               syst.rhs{j}(dof) = syst.rhs{j}(dof) + rhsVal(ind);
            end
         end
      case 'Neu'
         for j = unique(phDofs)'
            ind = phDofs == j;
            dof = locDofs(ind);
            syst.rhs{j}(dof) = syst.rhs{j}(dof) + rhsVal(ind);
         end
      otherwise
         for j = unique(phDofs)'
            ind = phDofs == j;
            dof = locDofs(ind);
            syst.rhs{j}(dof) = syst.rhs{j}(dof) + rhsVal(ind);
         end
   end
end

end