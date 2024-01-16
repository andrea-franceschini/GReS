classdef Boundaries < handle
  % BOUNDARY CONDITIONS - General boundary conditions class

  properties (Access = public)
    % Creation of a Map object for the boundary conditions
    db
    dof
  end
  
    properties (Access = private)
    model
    grid
  end

  methods (Access = public)
    % Class constructor method
    function obj = Boundaries(fileNames,model,grid,dof) %,model,grid
      % MATLAB evaluates the assignment expression for each instance, which
      % ensures that each instance has a unique value
      obj.db = containers.Map('KeyType','char','ValueType','any');
      obj.dof = dof;
      obj.model = model;
      obj.grid = grid;
      % Calling the function to read input data from file
      obj.readInputFiles(fileNames);
      obj.computeBoundaryProperties(model,grid);
      linkBoundSurf2TPFAFace(model,obj,grid);
    end
    
    function delete(obj)
      remove(obj.db,keys(obj.db));
      obj.db = [];
    end

    % Check if the identifier defined by the user is a key of the Map object
    function bc = getData(obj,identifier)
      if (obj.db.isKey(identifier))
        bc = obj.db(identifier);
      else
        % Displaying error message if the identifier does not refer
        % to any existing class
        error('Boundary condition %s not present', identifier);
      end
    end
    
    function vals = getVals(obj, identifier, t)
      vals = obj.getData(identifier).data.getValues(t);
    end

    function list = getDofs(obj, identifier) 
      %return loaded DOFS for the specified BC in global indexing
      %list = obj.getData(identifier).data.entities; OLD VERSION
      %%%%update to getDofs method
      col = obj.dof.getColTable(obj.getPhysics(identifier));
      if strcmp(obj.getCond(identifier),'NodeBC') | strcmp(obj.getCond(identifier),'ElementBC') 
            nEnts = obj.getData(identifier).data.nEntities;
            entities = obj.getData(identifier).data.entities;
            i1 = 1;
            for i = 1:length(nEnts)
              i2 = i1 + nEnts(i);
              if strcmp(obj.getCond(identifier),'NodeBC') % Node BC ---> Node Dof
                list(i1:i2-1) = obj.dof.nodeDofTable(entities(i1:i2-1),col(i));
              elseif strcmp(obj.getCond(identifier),'ElementBC') % Element BC ---> Element Dof
                list(i1:i2-1) = obj.dof.elemDofTable(entities(i1:i2-1),col(i));  
              end
              i1 = i2;
            end
      elseif strcmp(obj.getCond(identifier),'SurfBC')
          % SurfBC ---> Node Dof
          if isFEMBased(obj.model,obj.getPhysics(identifier))
              if strcmp(obj.getType(identifier),'Neu')
                  loadedEnts = obj.getLoadedEntities(identifier);
                  if strcmp(obj.getPhysics(identifier),'Poro')
                      direction = obj.getDirection(identifier);
                      switch direction
                          case 'x'
                            list = obj.dof.nodeDofTable(loadedEnts,col(1));
                          case 'y'
                            list = obj.dof.nodeDofTable(loadedEnts,col(2));
                          case 'z'
                            list = obj.dof.nodeDofTable(loadedEnts,col(3));
                      end
                  elseif strcmp(obj.getPhysics(identifier),'Flow')
                      list = obj.dof.nodeDofTable(loadedEnts,col);
                  end
              elseif strcmp(obj.getType(identifier),'Dir')           
                nEnts = obj.getNumbLoadedEntities(identifier);
                ents = obj.getLoadedEntities(identifier);
                i1 = 1;
                for i = 1:length(nEnts)
                  i2 = i1 + nEnts(i);
                  list(i1:i2-1) = obj.dof.nodeDofTable(ents(i1:i2-1),col(i));
                  i1 = i2;
                end
              end
            elseif isFVTPFA(obj.model,bound.getPhysics(identifier))
          % SurfBC ---> Element Dof
                if  strcmp(obj.getPhysics(identifier),'Poro')
                    direction = obj.getDirection(identifier);
                    switch direction
                        case 'x'
                            list = obj.dof.elemDofTable(loadedEnts,col(1));
                        case 'y'
                            list = obj.dof.elemDofTable(loadedEnts,col(2));
                        case 'z'
                            list = obj.dof.elemDofTable(loadedEnts,col(3));
                    end
                elseif strcmp(obj.getPhysics(identifier),'Flow')
                    list = obj.dof.elemDofTable(loadedEnts,col);
                end
          end
      elseif strcmp(obj.getCond(identifier),'VolumeForce')
          % Volume Force ---> Node Dof
          %VolumeForce are available only for flow model
          if isFEMBased(obj.model,obj.getPhysics(identifier))
              loadedEnts = obj.getLoadedEntities(identifier);
              list = obj.dof.nodeDofTable(loadedEnts,col);
          elseif isFVTPFABased(obj.model,obj.getPhysics(identifier))
          % Volume Force ---> Element Dof
              ents = obj.getData(identifier).data.entities;
              list = obj.dof.elemDofTable(ents,col);
          end
      end
      if any(list == 0)
          error('Boundary conditions %s not supported from subdomain',identifier)
      end
               
    end

    function list = getLocDofs(obj, identifier)
      %return loaded DOFS for the specified BC in global indexing
      %list = obj.getData(identifier).data.entities; OLD VERSION
      %%%%update to getDofs method
      col = obj.dof.getColTable(obj.getPhysics(identifier));
      if strcmp(obj.getCond(identifier),'NodeBC') | strcmp(obj.getCond(identifier),'ElementBC') 
            nEnts = obj.getData(identifier).data.nEntities;
            entities = obj.getData(identifier).data.entities;
            i1 = 1;
            for i = 1:length(nEnts)
              i2 = i1 + nEnts(i);
              if strcmp(obj.getCond(identifier),'NodeBC') % Node BC ---> Node Dof
                list(i1:i2-1) = obj.dof.nodeDofTable(entities(i1:i2-1),col(i),2);
              elseif strcmp(obj.getCond(identifier),'ElementBC') % Element BC ---> Element Dof
                list(i1:i2-1) = obj.dof.elemDofTable(entities(i1:i2-1),col(i));  
              end
              i1 = i2;
            end
      elseif strcmp(obj.getCond(identifier),'SurfBC')
          % SurfBC ---> Node Dof
          if isFEMBased(obj.model,obj.getPhysics(identifier))
              if strcmp(obj.getType(identifier),'Neu')
                  loadedEnts = obj.getLoadedEntities(identifier);
                  if strcmp(obj.getPhysics(identifier),'Poro')
                      direction = obj.getDirection(identifier);
                      switch direction
                          case 'x'
                            list = obj.dof.nodeDofTable(loadedEnts,col(1),2);
                          case 'y'
                            list = obj.dof.nodeDofTable(loadedEnts,col(2),2);
                          case 'z'
                            list = obj.dof.nodeDofTable(loadedEnts,col(3),2);
                      end
                  elseif strcmp(obj.getPhysics(identifier),'Flow')
                      list = obj.dof.nodeDofTable(loadedEnts,col);
                  end
              elseif strcmp(obj.getType(identifier),'Dir')           
                nEnts = obj.getNumbLoadedEntities(identifier);
                ents = obj.getLoadedEntities(identifier);
                i1 = 1;
                for i = 1:length(nEnts)
                  i2 = i1 + nEnts(i);
                  list(i1:i2-1) = obj.dof.nodeDofTable(ents(i1:i2-1),col(i),2);
                  i1 = i2;
                end
              end
            elseif isFVTPFA(obj.model,bound.getPhysics(identifier))
          % SurfBC ---> Element Dof
                if  strcmp(obj.getPhysics(identifier),'Poro')
                    direction = obj.getDirection(identifier);
                    switch direction
                        case 'x'
                            list = obj.dof.elemDofTable(loadedEnts,col(1),2);
                        case 'y'
                            list = obj.dof.elemDofTable(loadedEnts,col(2),2);
                        case 'z'
                            list = obj.dof.elemDofTable(loadedEnts,col(3),2);
                    end
                elseif strcmp(obj.getPhysics(identifier),'Flow')
                    list = obj.dof.elemDofTable(loadedEnts,col,2);
                end
          end
      elseif strcmp(obj.getCond(identifier),'VolumeForce')
          % Volume Force ---> Node Dof
          %VolumeForce are available only for flow model
          if isFEMBased(obj.model,obj.getPhysics(identifier))
              loadedEnts = obj.getLoadedEntities(identifier);
              list = obj.dof.nodeDofTable(loadedEnts,col,2);
          elseif isFVTPFABased(obj.model,obj.getPhysics(identifier))
          % Volume Force ---> Element Dof
              ents = obj.getData(identifier).data.entities;
              list = obj.dof.elemDofTable(ents,col,2);
          end
      end
      if any(list == 0)
          error('Boundary conditions %s not supported from subdomain',identifier)
      end
               
    end
    
    function cond = getCond(obj, identifier)
      cond = obj.getData(identifier).cond;
    end
    
    function name = getName(obj, identifier)
      name = obj.getData(identifier).data.name;
    end

    function type = getType(obj, identifier)
      type = obj.getData(identifier).type;
    end

    function physics = getPhysics(obj, identifier)
      physics = obj.getData(identifier).physics;
    end
    
    function ents = getEntities(obj, identifier, varargin)
      % notempty varargin ---> return entities with component multiplication
      %this method return the index of constrained entities inside solution vectors
      %needed for applying Dirichlet BCs
      %(pressure, displacement)...  
      ents = obj.getData(identifier).data.entities; 
      % consider solution components
      if ~isempty(varargin)
          nEnts = obj.getData(identifier).data.nEntities;
          comps = length(nEnts);
          %return consistent entitity with discretization DoF
          %Example FEM - SurfBC ---> NodeBC 
          %handle multicomponents dof
          i1 = 1;
          for i = 1:comps
              i2 = i1 + nEnts(i);
              ents(i1:i2-1) = comps*(ents(i1:i2-1)-1)+i;
              i1 = i2;
          end
      end
     end
    
     function ents = getLoadedEntities(obj, identifier,varargin)
      % return loaded entities for Volume or Surface BCs
      % if varargin not empty, return loaded ents with component
      % multiplication
      ents = obj.getData(identifier).loadedEnts;
      if ~isempty(varargin)
          nEnts = obj.getNumbLoadedEntities(identifier);
          comps = length(nEnts);
          i1 = 1;
          for i = 1:comps
              i2 = i1 + nEnts(i);
              ents(i1:i2-1) = comps*(ents(i1:i2-1)-1)+i;
              i1 = i2;
          end
      end
          

    end

    function ent = getNumbLoadedEntities(obj, identifier)
      ent = obj.getData(identifier).nloadedEnts;
    end
    
    function infl = getEntitiesInfluence(obj, identifier)
      infl = obj.getData(identifier).entitiesInfl;
    end
    
    function direction = getDirection(obj, identifier)
      direction = obj.getData(identifier).direction;
    end
    
    function setDofs(obj, identifier, list)
      obj.getData(identifier).data.entities = list;
    end
    
    function computeBoundaryProperties(obj,model,grid)
      keys = obj.db.keys;
      for i = 1 : length(keys)
        if strcmp(obj.getCond(keys{i}), 'VolumeForce') && ...
            isFEMBased(model,obj.getPhysics(keys{i})) 
          dofs = obj.getEntities(keys{i});
          tmpMat = grid.topology.cells(dofs,:)';
          [loadedEnts] = unique(tmpMat(tmpMat ~= 0));
          ptrHexa = grid.topology.cellVTKType(dofs) == 12;
          ptrTetra = grid.topology.cellVTKType(dofs) == 10;
          nHexa = nnz(ptrHexa);
          nTetra = nnz(ptrTetra);
          rowID = zeros(8*nHexa+4*nTetra,1);
          colID = zeros(8*nHexa+4*nTetra,1);
          nodeVol = zeros(8*nHexa+4*nTetra,1);
          ptr = 0;
          if any(ptrHexa)
            nodeVol(ptr+1:ptr+8*nHexa) = grid.cells.hexa.findNodeVolume(dofs(ptrHexa)');
            topolHexa = tmpMat(:,ptrHexa);
            [~,rowID(ptr+1:ptr+8*nHexa)] = ismember(topolHexa(topolHexa ~= 0),loadedEnts);
            colID(ptr+1:ptr+8*nHexa) = repelem(find(ptrHexa),8);
            ptr = ptr + 8*nHexa;
            clear topolHexa
          end
          if any(ptrTetra)
            nodeVol(ptr+1:ptr+4*nTetra) = 0.25*repelem(grid.cells.vol(dofs(ptrTetra)),4);
            topolTetra = tmpMat(:,ptrTetra);
            [~,rowID(ptr+1:ptr+4*nTetra)] = ismember(topolTetra(topolTetra ~= 0),loadedEnts);
            colID(ptr+1:ptr+4*nTetra) = repelem(find(ptrTetra),4);
            clear topolTetra
          end
          nodeVolume = sparse(rowID,colID,nodeVol,length(loadedEnts),length(dofs));
          tmpDbEntry = obj.getData(keys{i});
          tmpDbEntry.entitiesInfl = nodeVolume;
          tmpDbEntry.loadedEnts = loadedEnts;
          obj.db(keys{i}) = tmpDbEntry;
        elseif strcmp(obj.getCond(keys{i}), 'SurfBC') && ...
          isFEMBased(model,obj.getPhysics(keys{i})) && ...
          strcmp(obj.getType(keys{i}),'Neu')
          dofs = obj.getEntities(keys{i});
          tmpMat = grid.topology.surfaces(dofs,:)';
          [loadedEnts] = unique(tmpMat(tmpMat ~= 0));
          ptrTri = grid.topology.surfaceVTKType(dofs) == 5;
          ptrQuad = grid.topology.surfaceVTKType(dofs) == 9;
          nQuad = nnz(ptrQuad);
          nTri = nnz(ptrTri);
          rowID = zeros(4*nQuad+3*nTri,1);
          colID = zeros(4*nQuad+3*nTri,1);
          areaSurf = zeros(4*nQuad+3*nTri,1);
          ptr = 0;
          if any(ptrQuad)
            areaSurf(ptr+1:ptr+4*nQuad) = 0.25*repelem(grid.faces.computeAreaQuad(dofs(ptrQuad)),4);
            topolQuad = tmpMat(:,ptrQuad);
            [~,rowID(ptr+1:ptr+4*nQuad)] = ismember(topolQuad(topolQuad ~= 0),loadedEnts);
            colID(ptr+1:ptr+4*nQuad) = repelem(find(ptrQuad),4);
            ptr = ptr + 4*nQuad;
            clear topolQuad
          end
          if any(ptrTri)
            areaSurf(ptr+1:ptr+3*nTri) = 1/3*repelem(grid.faces.computeAreaTri(dofs(ptrTri)),3);
            topolTri = tmpMat(:,ptrTri);
            [~,rowID(ptr+1:ptr+3*nTri)] = ismember(topolTri(topolTri ~= 0),loadedEnts);
            colID(ptr+1:ptr+3*nTri) = repelem(find(ptrTri),3);
            clear topolTri
          end
          nodeArea = sparse(rowID,colID,areaSurf,length(loadedEnts),length(dofs));
          %
          %this instructions will be later moved in the general
          %getDofs method
%           if strcmp(obj.getPhysics(keys{i}),'Poro')
%             direction = obj.getDirection(keys{i});
%             switch direction
%               case 'x'
%                 loadedEnts = 3*loadedEnts - 2;
%               case 'y'
%                 loadedEnts = 3*loadedEnts - 1;
%               case 'z'
%                 loadedEnts = 3*loadedEnts;
%             end
%           end
          tmpDbEntry = obj.getData(keys{i});
          tmpDbEntry.entitiesInfl = nodeArea;
          tmpDbEntry.loadedEnts = loadedEnts;
          obj.db(keys{i}) = tmpDbEntry;
        elseif strcmp(obj.getCond(keys{i}), 'SurfBC') && ...
          isFEMBased(model,obj.getPhysics(keys{i})) && ...
          strcmp(obj.getType(keys{i}),'Dir')
          dofs = obj.getEntities(keys{i});
          nEnts = obj.getData(keys{i}).data.nEntities;
          comps = length(nEnts);
          nLoadEnts = zeros(comps,1);
          i1 = 1; j1 = 1;
          nodes = []; sid = [];
          for ic = 1: comps
              i2 = i1+nEnts(ic);
              surfs = (grid.topology.surfaces(dofs(i1:i2-1),:))';
              [nod, ind, ~] = unique(surfs,'last');
              nLoadEnts(ic) = length(nod);
              j2 = j1 + nLoadEnts(ic);
              [~, s] = ind2sub(size(surfs),ind);
              sid(j1:j2-1,1) = s + i1-1;
              nodes(j1:j2-1,1) = nod;
              i1 = i2;
              j1 = j2;
              %               i2 = i1 + nEnts(i);
              % ents(i1:i2-1) = comps*(ents(i1:i2-1)-1)+i;
              % i1 = i2;
          end
          nnod = length(nodes);
          mapNodSurf = sparse(1:nnod,sid,ones(nnod,1));
          tmpDbEntry = obj.getData(keys{i});
          tmpDbEntry.entitiesInfl = mapNodSurf;
          tmpDbEntry.loadedEnts = nodes;
          tmpDbEntry.nloadedEnts = nLoadEnts;
          obj.db(keys{i}) = tmpDbEntry;
        end
      end
    end
    
%     function iniBC(obj,BCList,state)
%       % Initialize BC
%       l = length(BCList);
%       assert(l>0,'Warning: No boundary conditions will be applied.');
%       obj.sizeM = max(length(state.displ),length(state.pressure));
%       obj.BCName = BCList;
%       obj.resState = state;
%     end
    
%     function du = applyDirVal(obj)
%       du = zeros(obj.sizeM,1);
%       for i=1:length(obj.BCName)
%         cond = getBC(obj,obj.BCName(i));
%         %
%         if strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
%           du(cond.boundDof) = cond.boundVal;
%         end
%       end
%     end
    %
%     function applyDirVal(obj,t)
%       % Apply Dirichlet conditions to the solution vector
%       for i=1:length(obj.BCName)
%         cond = getBC(obj,obj.BCName(i));
%         %
%         % Check if the BC need to be updated (we include this check here to
%         % limit the number of times the db is accessed)
%         %
%         Boundaries.checkAndUpdateCond(cond,t);
%         if isa(cond,'NodeBC')
%           if strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
%             % Interpolate in the step interval
%             % fac = (t-t_old)/(t_new-t_old)
%             fac = (t - cond.timeInt(1,cond.timeInt(2,1)))/ ...
%               (cond.timeInt(1,cond.timeInt(2,2)) - ...
%                cond.timeInt(1,cond.timeInt(2,1)));
%             tmpVec = (1-fac)*cond.boundVal(:,cond.timeInt(2,1)) + ...
%               fac*cond.boundVal(:,cond.timeInt(2,2));
%             if strcmp(cond.boundPhysics,'flow')
%               obj.resState.pressure(cond.boundDof) = tmpVec;
%             elseif strcmp(cond.boundPhysics,'poro')
%               obj.resState.displ(cond.boundDof) = tmpVec;
%             end
%           end
%         end
%       end
%     end
    
%     function applyBCandForces(obj,syst,t)
%       % Impose BC to the linearized system (Jacobian matrix + RHS)
%       % The Penalty method is used for the Dirichlet conditions
%       %
% %       l = length(obj.BCName);
% %       assert(l>0,'Warning: No boundary conditions will be applied.');
% %       if strcmp(probType,'lin'); fac = 1; else; fac = -1; end
%       maxVal = max(abs(syst.K), [], 'all');
%       %
%       for i=1:length(obj.BCName)
%         cond = getBC(obj,obj.BCName(i));
%         % fac = (t-t_old)/(t_new-t_old)
%         fac = (t - cond.timeInt(1,cond.timeInt(2,1)))/ ...
%           (cond.timeInt(1,cond.timeInt(2,2)) - ...
%            cond.timeInt(1,cond.timeInt(2,1)));
%         tmpVec = (1-fac)*cond.boundVal(:,cond.timeInt(2,1)) + ...
%           fac*cond.boundVal(:,cond.timeInt(2,2));
%         if isa(cond,'NodeBC')
%           if strcmp(cond.boundType,'neu')  % Apply Neumann conditions,if any
%             syst.rhs(cond.boundDof) = syst.rhs(cond.boundDof) - tmpVec;
%           elseif strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
%             syst.K(obj.sizeM*(cond.boundDof-1)+cond.boundDof) = maxVal*10^10;
%             syst.rhs(cond.boundDof) = 0;
%   %           if strcmp(probType,'lin')
%   %             syst.rhs(cond.boundDof) = cond.boundVal*(maxVal*10^10);
%   %           else
%   %             varargin{1}.displ(cond.boundDof) = cond.boundVal;
%   %             syst.rhs(cond.boundDof) = 0;
%   %           end
%           end
%         elseif isa(cond,'VolumeForce')
%           if isFEMBased(obj.model)
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
%         end
%       end
%     end
    
%     function applyBCNeu(obj,syst)
% %       l = length(obj.BCName);
% %       assert(l>0,'Warning: No boundary conditions will be applied.');
% %       if strcmp(probType,'lin'); fac = 1; else; fac = -1; end
% %       maxVal = max(abs(syst.K), [], 'all');
%       %
%       for i=1:length(obj.BCName)
%         cond = getBC(obj,obj.BCName(i));
%         if strcmp(cond.boundType,'neu')  % Apply Neumann conditions,if any
%           syst.rhs(cond.boundDof) = syst.rhs(cond.boundDof) - cond.boundVal;
%         end
%         %
% %         if strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
% %           syst.K(obj.sizeM*(cond.boundDof-1)+cond.boundDof) = maxVal*10^10;
% %           syst.rhs(cond.boundDof) = 0;
% %           if strcmp(probType,'lin')
% %             syst.rhs(cond.boundDof) = cond.boundVal*(maxVal*10^10);
% %           else
% %             varargin{1}.displ(cond.boundDof) = cond.boundVal;
% %             syst.rhs(cond.boundDof) = 0;
% %           end
% %         end
%       end
%     end
%     
%     function applyBCDir(obj,syst)
% %       l = length(obj.BCName);
% %       assert(l>0,'Warning: No boundary conditions will be applied.');
% %       if strcmp(probType,'lin'); fac = 1; else; fac = -1; end
%       maxVal = max(abs(syst.K), [], 'all');
%       %
%       for i=1:length(obj.BCName)
%         cond = getBC(obj,obj.BCName(i));
% %         if strcmp(cond.boundType,'neu')  % Apply Neumann conditions,if any
% %           syst.rhs(cond.boundDof) = syst.rhs(cond.boundDof) - cond.boundVal;
% %         end
%         %
%         if strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
%           syst.K(obj.sizeM*(cond.boundDof-1)+cond.boundDof) = maxVal*10^10;
%           syst.rhs(cond.boundDof) = 0;
% %           if strcmp(probType,'lin')
% %             syst.rhs(cond.boundDof) = cond.boundVal*(maxVal*10^10);
% %           else
% %             varargin{1}.displ(cond.boundDof) = cond.boundVal;
% %             syst.rhs(cond.boundDof) = 0;
% %           end
%         end
%       end
%     end
%     function checkBCStep(obj,t)
%       % NOTE: In case of backstep with simulation time falling in the
%       % previous BC step we need to go back in the file and read the
%       % old condition again. Use the fseek command (TO BE IMPLEMENTED)
%       for i=1:length(obj.BCName)
%         cond = getBC(obj,obj.BCName(i));
%         while max(cond.timeInt(1,:) < t)
%           cond.updateBC();
%         end
%       end
%     end
  end

  methods (Access = private)
    % Reading boundary input file
    function readInputFile(obj,fileName)
      fid = fopen(fileName, 'r');
      if (fid == -1)
        error('File %s not opened correctly',fileName);
      end
      token = Boundaries.readToken(fid);
      if (~ismember(convertCharsToStrings(token), ["NodeBC", "SurfBC", "VolumeForce","ElementBC"]))
        error(['%s condition is unknown\n', ...
          'Accepted types are: NodeBC   -> Boundary cond. on nodes\n',...
          '                    SurfBC   -> Boundary cond. on surfaces\n',...
          '                    ElementBC   -> Boundary cond. on elements\n',...                      
          '                    VolumeForce -> Volume force on elements'], token);
      end
      if ismember(convertCharsToStrings(token), ["NodeBC", "SurfBC", "ElementBC"])
        type = Boundaries.readToken(fid);
        if (~ismember(type, ['Dir', 'Neu']))
          error(['%s boundary condition is not admitted\n', ...
            'Accepted types are: Dir -> Dirichlet, Neu -> Neumann'], type);
        end
      end
      physics = Boundaries.readToken(fid);
      
      if strcmp(physics,'Poro') && strcmp(token,'SurfBC') && strcmp(type,'Neu') 
        direction = Boundaries.readToken(fid);
        if ~ismember(direction,['x','y','z'])
          error(['%s is an invalid direction of the distributed load\n', ...
            'Accepted directions are: x, y, and z'],direction);
        end
      end
      name = Boundaries.readToken(fid);
      setFile = Boundaries.readToken(fid);
      [times, dataFiles] = Boundaries.readDataFiles(fid);
      fclose(fid);
      if obj.db.isKey(name)
        error('%s boundary condition name already defined', name);
      end
      switch token
        case 'NodeBC'
          obj.db(name) = struct('data', BoundaryEntities(name, setFile, times, dataFiles), ...
            'cond',token,'type', type, 'physics', physics);
        case 'SurfBC'
          switch physics
            case 'Flow'
              obj.db(name) = struct('data', BoundaryEntities(name, setFile, times, dataFiles), ...
                'cond',token,'type', type, 'physics', physics);
            case 'Poro'
                switch type
                    case 'Neu'
                      obj.db(name) = struct('data', BoundaryEntities(name, setFile, times, dataFiles), ...
                        'cond', token,'direction', direction, 'type', type, 'physics', physics);
                    case 'Dir'
                      obj.db(name) = struct('data', BoundaryEntities(name, setFile, times, dataFiles), ...
                        'cond', token,'type', type, 'physics', physics);                        
                end
          end
        case {'VolumeForce', 'ElementBC'}
          obj.db(name) = struct('data', BoundaryEntities(name, setFile, times, dataFiles), ...
            'cond',token, 'physics', physics);
      end
      %
    end
    
    % Reading boundary input file
    function readInputFiles(obj,fileNames)
      n = length(fileNames);
      assert(n > 0,'No boundary conditions are to be imposed');
      for i = 1 : n
        readInputFile(obj,fileNames(i));
      end
    end
    
%     % Set Boundaries class
%     function setBoundaries(obj,symmod,grid)
%       obj.model = symmod;
%       obj.mesh = grid.topology;
%       obj.elements = grid.cells;
%     end
    % Reading boundary input file
%     function readInputFile(obj,fileName)
%       nBC = length(fileName);
%       assert(nBC > 0,'No boundary conditions are to be imposed');
% %       BCType = lower(BCType);
% %       if ~ismember(BCType,["dir", "neu"])
% %         error('%s boundary condition is not admitted\n Accepted types are: dir -> Dirichlet, neu -> Neumann',BCType);
% %       end
%       %
%       for i=1:nBC
%         fid = fopen(fileName(i), 'r');
%         if fid == -1
%           error('File %s not opened correctly',fileName(i));
%         end
%         block = '';
%         [flEof,line] = Boundaries.readLine(fid);
%         % Reading the BC
%         while (~isempty(line))
%           if flEof == 1
%             error('End of file while reading BC in file %s',fileName(i));
%           end
%           line = strtrim(strtok(line,'%'));
%           block = [block, line, ' '];
%           [flEof,line] = Boundaries.readLine(fid);
%         end
%         blockSplt = strsplit(string(deblank(block)));
%         % Calling the specific BC class based on the BC name
%         %
%         % Structure of blockSplt (see also setBC):
%         %       blockSplt(1) -> label defining the discrete item on which 
%         %                       the BC is applied
%         %       blockSplt(2) -> type of BC (Dir = Dirichlet, Neu = Neumann)
%         %       blockSplt(3) -> physical process for which the BC is
%         %                       specified (Poro = poromechanics, 1PFlow = 
%         %                       single-phase flow)
%         %       blockSplt(4) -> BC name
%         if ~isnan(str2double(blockSplt(1)))
%           error('The first entry in %s must be a string',fileName(i));
%         end
%         switch lower(blockSplt(1))
%           case 'nodebc'
%             pos = 1:4;
%           case 'volumeforce'
%             pos = 1:3;
%           otherwise
%             error('Condition %s not available', blockSplt(1));
%         end
%         %
%         if ~all(isnan(str2double(blockSplt(pos(2:end)))))
%           error('The first %d entries in %s must be strings',pos(end),fileName(i));
%         end
%         blockSplt(pos) = lower(blockSplt(pos));
%         if obj.db.isKey(blockSplt(pos(end)))
%           error('%s condition already defined',blockSplt(pos(end)));
%         end
%         %
%         switch blockSplt(1)
%           case 'nodebc'
%             obj.db(blockSplt(pos(end))) = NodeBC(fid,blockSplt);
%           case 'volumeforce'
%             obj.db(blockSplt(pos(end))) = VolumeForce(fid,blockSplt);
%         end
%         cond = getBC(obj,blockSplt(pos(end)));
%         % Read the BC values at t=0 and at the end of the event
%         Boundaries.readStep(cond,1);
%         Boundaries.readStep(cond,2);
% %         fclose(fid);
%       end
%     end
  end
  
  methods(Static = true)
    % Read the next token and check for eof
    function [token] = readToken(fid)
      flEof = feof(fid);   % end-of-file flag
      if flEof == 1
        error('No token available in boundary condition file.');
      else
        token = sscanf(fgetl(fid), '%s', 1);
      end
    end
    
    function [times, data] = readDataFiles(fid)
      nDataMax = 100;
      data = repmat(struct('time', 0, 'fileName', []), nDataMax, 1);
      times = zeros(nDataMax,1);
      id = 0;
      while (~feof(fid))
        line = fgetl(fid);
        if (strcmp(line, 'End'))
          break;
        end
        word = sscanf(line, '%s', 1);
        if (~strcmp(word(1), '%'))
          [time, ~, ~, pos] = sscanf(line, '%e', 1);
          id = id + 1;
          if (id > nDataMax)
            nDataMax = 2*nDataMax;
            data(nDataMax) = data(1);
            times(nDataMax) = 0.0;
          end
          times(id) = time;
          data(id).time = time;
          data(id).fileName = strtrim(line(pos:end));
        end
      end
      data = data(1:id);
      times = times(1:id);
    end
    
%     function checkAndUpdateCond(cond,t)
%       % Check if we move to the next event
%       while max(cond.timeInt(1,:)) < t
%         if cond.timeInt(1,1) < cond.timeInt(1,2)
% %         pos = 1;
%           cond.timeInt(2,1) = 2;
%           cond.timeInt(2,2) = 1;
%         else
% %         pos = 2;
%           cond.timeInt(2,1) = 1;
%           cond.timeInt(2,2) = 2;
%         end
%         Boundaries.readStep(cond,cond.timeInt(2,2));
%       end
%     end
 
%     function readStep(cond,pos)
%       % Read the BC values in the next event
%       try
%       cond.timeInt(1,pos) = fscanf(cond.fID,['TIME' '%e'],[1 1]);
%       catch
%         error('Wrong format of the TIME header in %s %s BC at time %f', ...
%           cond.boundPhysics,cond.boundType,cond.timeInt(1,pos))
%       end
%       try
%         cond.boundVal(:,pos) = fscanf(cond.fID,'%e\n',[1 length(cond.boundDof)]);
%       catch
%         error('Wrong number of values in %s %s BC at time %f', ...
%         cond.boundPhysics,cond.boundType,cond.timeInt(1,pos))
%       end
% %       if ~feof(obj.fID); tmpLine = fgetl(obj.fID); end
%     end
    
%     % Read the next line and check for eof
%     function [flEof,line] = readLine(fid)
%       flEof = feof(fid);   % end-of-file flag
%       if flEof == 1
%         line = '';
%       else
%         line = deblank(fgetl(fid));
% %         if isempty(line)
% %           error('No blank lines are admitted in BC file')
% %         end
%       end  
%     end
  end
end