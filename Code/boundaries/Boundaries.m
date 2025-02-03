classdef Boundaries < handle
  % BOUNDARY CONDITIONS - General boundary conditions class

  properties (Access = public)
    % Creation of a Map object for the boundary conditions
    db
  end
  
    properties (Access = private)
    model
    grid
  end

  methods (Access = public)
    % Class constructor method
    function obj = Boundaries(fileNames,model,grid) %,model,grid
      % MATLAB evaluates the assignment expression for each instance, which
      % ensures that each instance has a unique value
      obj.db = containers.Map('KeyType','char','ValueType','any');
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
 
    function cond = getCond(obj, identifier)
      cond = obj.getData(identifier).cond;
    end
    
    function name = getName(obj, identifier)
      name = obj.getData(identifier).data.name;
    end

    function type = getType(obj, identifier)
       if ~strcmp(obj.getCond(identifier),'VolumeForce')
          type = obj.getData(identifier).type;
       else
          type = 'VolumeForce';
       end
    end

    function physics = getPhysics(obj, identifier)
      physics = obj.getData(identifier).physics;
    end
   
    function dofs = getCompEntities(obj,identifier,ents)
       % get component dof of Dirichlet BC loaded entities
       nEnts = getNumbLoadedEntities(obj,identifier);
       % component multiplication of BC entities
       dim = length(nEnts);
       i1 = 1;
       dofs = zeros(numel(ents),1);
       for i = 1 : dim
          i2 = i1 + nEnts(i);
          dofs(i1:i2-1) = dim*(ents(i1:i2-1)-1) + i;
          i1 = i2;
       end
    end

    function ents = getEntities(obj,identifier)
       ents = obj.getData(identifier).data.entities;
    end
    
     function ents = getLoadedEntities(obj, identifier)
        % return loaded entities for Volume or Surface BCs
        ents = obj.getData(identifier).loadedEnts;
     end

     function nEnts = getNumbEntities(obj,identifier)
        nEnts = obj.getData(identifier).data.nEntities;
     end

     function ents = getNumbLoadedEntities(obj, identifier)
        bc = obj.getData(identifier);
        if isfield(bc,'nloadedEnts')
           % Surface BC
           ents = bc.nloadedEnts;
        else
           % Node BC
           ents = bc.data.nEntities;
        end
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
          case {'NodeBC', 'ElementBC'}
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
        case 'VolumeForce'
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
        if (strcmpi(line, 'End'))
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
  end
end