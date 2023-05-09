 function applyBCandForces(model, grid, bound, t, syst,theta,dt)
  % Impose BC to the linearized system (Jacobian matrix + RHS)
  % The Penalty method is used for the Dirichlet conditions
  %
  maxVal = max(abs(syst.K), [], 'all');
  %
  keys = bound.db.keys;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%IF COUPLED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %BCs added according to formulation n.1 in IJNAMG2009
  if isCoupled(model)
    for i = 1 : length(keys)
    if strcmp(bound.getCond(keys{i}),'NodeBC')
      if strcmp(bound.getPhysics(keys{i}),'Poro')
          if strcmp(bound.getType(keys{i}), 'Neu')  % Apply Neumann conditions,if any
            syst.rhsPoro(bound.getDofs(keys{i})) = syst.rhsPoro(bound.getDofs(keys{i})) - (1/theta)*bound.getVals(keys{i}, t);
          elseif strcmp(bound.getType(keys{i}), 'Dir')  % Apply Dirichlet conditions
            nrows = size(syst.KPoro,1);
            syst.KPoro(nrows*(bound.getDofs(keys{i})-1) + bound.getDofs(keys{i})) = maxVal*1.e10;
            syst.rhsPoro(bound.getDofs(keys{i})) = 0;
          end
      elseif strcmp(bound.getPhysics(keys{i}),'Flow')
          if strcmp(bound.getType(keys{i}), 'Neu')  % Apply Neumann conditions,if any
            syst.rhsFlow(bound.getDofs(keys{i})) = syst.rhsFlow(bound.getDofs(keys{i})) + dt*bound.getVals(keys{i}, t);
          elseif strcmp(bound.getType(keys{i}), 'Dir')  % Apply Dirichlet conditions
            nrows = size(syst.KFlow,1);
            syst.KFlow(nrows*(bound.getDofs(keys{i})-1) + bound.getDofs(keys{i})) = maxVal*1.e10;
            syst.rhsFlow(bound.getDofs(keys{i})) = 0;
          end
      end
    elseif strcmp(bound.getCond(keys{i}),'ElementBC')
      nrows = size(syst.K,1);
      syst.K(nrows*(bound.getDofs(keys{i})-1) + bound.getDofs(keys{i})) = maxVal*1.e10;
      syst.rhs(bound.getDofs(keys{i})) = 0;
    elseif strcmp(bound.getCond(keys{i}),'SurfBC')
      if isFEMBased(model,bound.getPhysics(keys{i}))
        if strcmp(bound.getType(keys{i}), 'Neu')
%           nod2FaceCond = obj.grid.topology.surfaces(bound.getDofs(keys{i}),:);
%           nod2FaceCond = nod2FaceCond';
%           nod2faceVec = nod2FaceCond(nod2FaceCond ~= 0);
%           fval = repelem(bound.getVals(keys{i}, t),obj.grid.topology.surfaceNumVerts(bound.getDofs(keys{i})));
%           s = obj.grid.faces.areaNod
%           syst.rhs(nod2faceVec) = syst.rhs(nod2faceVec) +
          %
          % TO DO: Vectorize the following lines
          vals = bound.getVals(keys{i}, t);
          id = 0;
          for s = bound.getDofs(keys{i})'
            id = id + 1;
            nodes = grid.topology.surfaces(s,1:grid.topology.surfaceNumVerts(s));
            fval = vals(id)*ones(length(nodes),1);
            aNod = grid.faces.areaNod(grid.faces.mapAreaNod(s):grid.faces.mapAreaNod(s+1)-1);
            if length(aNod) == 1
              aNod = aNod*ones(3,1);
            end
            syst.rhsFlow(nodes) = syst.rhsFlow(nodes) + dt*fval.*aNod;
          end
        elseif strcmp(bound.getType(keys{i}), 'Dir')
          error('Dirichlet cond. for FEM on surf is not available (ref. %s)', ...
            bound.getToken(keys{i}));
        end
      elseif isFVTPFABased(model,bound.getPhysics(keys{i}))
        if strcmp(bound.getType(keys{i}), 'Neu')
          faceID = bound.getDofs(keys{i});
          neigh = grid.faces.faceNeighbors(faceID,:);
          assert(all(sum(neigh~=0,2) == 1),'Corrupted face numbering in %s',bound.getName(keys{i}));
          neigh = sum(neigh,2);
          syst.rhs(neigh) = syst.rhs(neigh) - bound.getVals(keys{i}, t);
        elseif strcmp(bound.getType(keys{i}), 'Dir')
          error('Dirichlet cond. for FVTPFA on surf has not been implemented yet (ref. %s)', ...
          bound.getName(keys{i}));
        end
      end
    elseif strcmp(bound.getCond(keys{i}), 'VolumeForce')
      if isFEMBased(model,bound.getPhysics(keys{i}))
        val = bound.getVals(keys{i}, t);
        ptr = 0;
        if any(grid.topology.cellVTKType(bound.getDofs(keys{i}))== 12)
          N1 = getBasisFinGPoints(grid.cells.hexa);
        end
        for el=(bound.getDofs(keys{i}))'
          ptr = ptr + 1;
          nodes = grid.topology.cells(el,1:grid.topology.cellNumVerts(el));
          switch grid.topology.cellVTKType(el)
            case 10 % Tetrahedron
              addRhs = (grid.cells.vol(el)/grid.topology.cellNumVerts(el))*repmat(val(ptr),[grid.topology.cellNumVerts(el),1]);
            case 12 % Hexahedron
              dJWeighed = getDerBasisFAndDet(grid.cells.hexa,el,3);
              addRhs = val(ptr)*N1'*dJWeighed';
          end
          syst.rhs(nodes) = syst.rhs(nodes) - addRhs;
        end
      elseif isFVTPFABased(model,bound.getPhysics(keys{i}))
        syst.rhs(bound.getDofs(keys{i})) = syst.rhs(bound.getDofs(keys{i})) - bound.getVals(keys{i}, t).*grid.cells.vol(bound.getDofs(keys{i}));
      end
    end
    end 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IF NOT COUPLED%%%%%%%%%%%%%%%%%%%%%%
  else
  for i = 1 : length(keys)
    if strcmp(bound.getCond(keys{i}),'NodeBC')
      if strcmp(bound.getType(keys{i}), 'Neu')  % Apply Neumann conditions,if any
        syst.rhs(bound.getDofs(keys{i})) = syst.rhs(bound.getDofs(keys{i})) - bound.getVals(keys{i}, t);
      elseif strcmp(bound.getType(keys{i}), 'Dir')  % Apply Dirichlet conditions
        nrows = size(syst.K,1);
        syst.K(nrows*(bound.getDofs(keys{i})-1) + bound.getDofs(keys{i})) = maxVal*1.e10;
        syst.rhs(bound.getDofs(keys{i})) = 0;
      end
    elseif strcmp(bound.getCond(keys{i}),'ElementBC')
      nrows = size(syst.K,1);
      syst.K(nrows*(bound.getDofs(keys{i})-1) + bound.getDofs(keys{i})) = maxVal*1.e10;
      syst.rhs(bound.getDofs(keys{i})) = 0;
    elseif strcmp(bound.getCond(keys{i}),'SurfBC')
      if isFEMBased(model,bound.getPhysics(keys{i}))
        if strcmp(bound.getType(keys{i}), 'Neu')
%           nod2FaceCond = obj.grid.topology.surfaces(bound.getDofs(keys{i}),:);
%           nod2FaceCond = nod2FaceCond';
%           nod2faceVec = nod2FaceCond(nod2FaceCond ~= 0);
%           fval = repelem(bound.getVals(keys{i}, t),obj.grid.topology.surfaceNumVerts(bound.getDofs(keys{i})));
%           s = obj.grid.faces.areaNod
%           syst.rhs(nod2faceVec) = syst.rhs(nod2faceVec) +
          %
          % TO DO: Vectorize the following lines
          vals = bound.getVals(keys{i}, t);
          id = 0;
          for s = bound.getDofs(keys{i})'
            id = id + 1;
            nodes = grid.topology.surfaces(s,1:grid.topology.surfaceNumVerts(s));
            fval = vals(id)*ones(length(nodes),1);
            aNod = grid.faces.areaNod(grid.faces.mapAreaNod(s):grid.faces.mapAreaNod(s+1)-1);
            if length(aNod) == 1
              aNod = aNod*ones(3,1);
            end
            syst.rhs(nodes) = syst.rhs(nodes) - fval.*aNod;
          end
        elseif strcmp(bound.getType(keys{i}), 'Dir')
          error('Dirichlet cond. for FEM on surf is not available (ref. %s)', ...
            bound.getToken(keys{i}));
        end
      elseif isFVTPFABased(model,bound.getPhysics(keys{i}))
        if strcmp(bound.getType(keys{i}), 'Neu')
          faceID = bound.getDofs(keys{i});
          neigh = grid.faces.faceNeighbors(faceID,:);
          assert(all(sum(neigh~=0,2) == 1),'Corrupted face numbering in %s',bound.getName(keys{i}));
          neigh = sum(neigh,2);
          syst.rhs(neigh) = syst.rhs(neigh) - bound.getVals(keys{i}, t);
        elseif strcmp(bound.getType(keys{i}), 'Dir')
          error('Dirichlet cond. for FVTPFA on surf has not been implemented yet (ref. %s)', ...
          bound.getName(keys{i}));
        end
      end
    elseif strcmp(bound.getCond(keys{i}), 'VolumeForce')
      if isFEMBased(model,bound.getPhysics(keys{i}))
        val = bound.getVals(keys{i}, t);
        ptr = 0;
        if any(grid.topology.cellVTKType(bound.getDofs(keys{i}))== 12)
          N1 = getBasisFinGPoints(grid.cells.hexa);
        end
        for el=(bound.getDofs(keys{i}))'
          ptr = ptr + 1;
          nodes = grid.topology.cells(el,1:grid.topology.cellNumVerts(el));
          switch grid.topology.cellVTKType(el)
            case 10 % Tetrahedron
              addRhs = (grid.cells.vol(el)/grid.topology.cellNumVerts(el))*repmat(val(ptr),[grid.topology.cellNumVerts(el),1]);
            case 12 % Hexahedron
              dJWeighed = getDerBasisFAndDet(grid.cells.hexa,el,3);
              addRhs = val(ptr)*N1'*dJWeighed';
          end
          syst.rhs(nodes) = syst.rhs(nodes) - addRhs;
        end
      elseif isFVTPFABased(model,bound.getPhysics(keys{i}))
        syst.rhs(bound.getDofs(keys{i})) = syst.rhs(bound.getDofs(keys{i})) - bound.getVals(keys{i}, t).*grid.cells.vol(bound.getDofs(keys{i}));
      end
    end
  end
  end
end

%   if isFEMBased(obj.model)
%             for el=cond.boundDof'
%               nodes = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
%               switch obj.mesh.cellVTKType(el)
%                 case 10 % Tetrahedron
%                   addRhs = (obj.elements.vol(el)/obj.mesh.cellNumVerts(el))*repmat(tmpVec(i),[obj.mesh.cellNumVerts(el),1]);
%                 case 12 % Hexahedron
%                   N1 = getBasisFinGPoints(obj.elements); % For better performance, N1 can be sought for only once since it is constant
%                   dJWeighed = getDerBasisFAndDet(obj.elements,el);
%                   addRhs = tmpVec(i)*N1'*dJWeighed';
%               end
%               syst.rhs(nodes) = syst.rhs(nodes) - addRhs;
%             end
%           elseif isFVTPFABased(obj.model)
%             syst.rhs(cond.boundDof) = syst.rhs(cond.boundDof) - tmpVec.*obj.elements.vol(cond.boundDof);
%           end
