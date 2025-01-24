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

    function list = getDofs(obj, identifier)
        %return loaded DOFS for the specified BC in global indexing
        %list = obj.getData(identifier).data.entities; OLD VERSION
        %%%%update to getDofs method
        ph = obj.getPhysics(identifier);
        col = obj.dof.getColTable(translatePhysic(ph, obj.model));
        if strcmp(obj.getCond(identifier),'NodeBC') || strcmp(obj.getCond(identifier),'ElementBC')
            nEnts = obj.getData(identifier).data.nEntities;
            entities = obj.getData(identifier).data.entities;
            i1 = 1;
            for i = 1:length(nEnts)
                i2 = i1 + nEnts(i);
                if strcmp(obj.getCond(identifier),'NodeBC') % Node BC ---> Node Dof
                    list(i1:i2-1) = obj.dof.nod2dof(entities(i1:i2-1),col(i));
                elseif strcmp(obj.getCond(identifier),'ElementBC') % Element BC ---> Element Dof
                    list(i1:i2-1) = obj.dof.elem2dof(entities(i1:i2-1),col(i));
                end
                i1 = i2;
            end
        elseif strcmp(obj.getCond(identifier),'SurfBC')
            % SurfBC ---> Node Dof
            if isFEMBased(obj.model,ph)
                if strcmp(obj.getType(identifier),'Neu')
                    loadedEnts = obj.getLoadedEntities(identifier);
                    if strcmp(ph,'Poro')
                        direction = obj.getDirection(identifier);
                        switch direction
                            case 'x'
                                list = obj.dof.nod2dof(loadedEnts,col(1));
                            case 'y'
                                list = obj.dof.nod2dof(loadedEnts,col(2));
                            case 'z'
                                list = obj.dof.nod2dof(loadedEnts,col(3));
                        end
                    elseif strcmp(ph,'Flow')
                        list = obj.dof.nod2dof(loadedEnts,col);
                    end
                elseif strcmp(obj.getType(identifier),'Dir')
                    nEnts = obj.getNumbLoadedEntities(identifier);
                    ents = obj.getLoadedEntities(identifier);
                    i1 = 1;
                    for i = 1:length(nEnts)
                        i2 = i1 + nEnts(i);
                        list(i1:i2-1) = obj.dof.nod2dof(ents(i1:i2-1),col(i));
                        i1 = i2;
                    end
                end
            elseif isFVTPFA(obj.model, ph)
                % SurfBC ---> Element Dof
                loadedEnts = obj.getLoadedEntities(identifier);
                if  strcmp(ph, 'Poro')
                    direction = obj.getDirection(identifier);
                    switch direction
                        case 'x'
                            list = obj.dof.elem2dof(loadedEnts,col(1));
                        case 'y'
                            list = obj.dof.elem2dof(loadedEnts,col(2));
                        case 'z'
                            list = obj.dof.elem2dof(loadedEnts,col(3));
                    end
                elseif strcmp(ph,'Flow')
                    list = obj.dof.elem2dof(loadedEnts,col);
                end
            end
        elseif strcmp(obj.getCond(identifier),'VolumeForce')
            % Volume Force ---> Node Dof
            %VolumeForce are available only for flow model
            if isFEMBased(obj.model,ph)
                loadedEnts = obj.getLoadedEntities(identifier);
                list = obj.dof.nod2dof(loadedEnts,col);
            elseif isFVTPFABased(obj.model,ph)
                % Volume Force ---> Element Dof
                ents = obj.getData(identifier).data.entities;
                list = obj.dof.elem2dof(ents,col);
            end
        end
        if any(list == 0)
            error('Boundary conditions %s not supported from subdomain',identifier)
        end
    end


    function list = getDofs_MD(obj,identifier,mG,d)
        % get BC Dofs in multidomain numbering
        ph = obj.getPhysics(identifier);
        list = [];
        nEnts = 1;
        % retrieve entities
        switch obj.getCond(identifier)
            case {'NodeBC','ElementBC'}
                nEnts = obj.getData(identifier).data.nEntities;
                ents = obj.getData(identifier).data.entities;
            case 'SurfBC'
                ents = obj.getLoadedEntities(identifier);
                if strcmp(obj.getType(identifier),'Dir')
                    nEnts = obj.getNumbLoadedEntities(identifier);
                end
            case 'VolumeForce'
                if isFEMBased(obj.model,ph)
                    ents = obj.getLoadedEntities(identifier);
                elseif isFVTPFABased(obj.model,ph)
                    % Volume Force ---> Element Dof
                    ents = obj.getData(identifier).data.entities;
                end
        end

        % check if any entity belong to requested field
        %if ~any(ismembc(ents,))

        
        % get component-wise numbering
        if nEnts > 1
            comps = length(nEnts);
            %handle multicomponents dof
            i1 = 1;
            for ii = 1:comps
                i2 = i1 + nEnts(ii);
                add_ents = find(ismember(mG.MD_struct(d).entities,ents(i1:i2-1)));
                add_ents = comps*(add_ents-1)+ii;
                list = [list;add_ents];
                i1 = i2;
            end
        else
            if strcmp(ph,'Poro')
                ents = find(ismembc(mG.MD_struct(d).entities,ents));
                c = find(strcmp(["x","y","z"],obj.getDirection(identifier)));
                list = 3*(ents-1)+c;
            else
               list = ents;
            end
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

    function ents = getCompEntities(obj,identifier,ents)
       % map entities number to component dof
       % consider solution components
       nEnts = obj.getNumbLoadedEntities(identifier);
       comps = length(nEnts);
       i1 = 1;
       for i = 1:comps
          i2 = i1 + nEnts(i);
          ents(i1:i2-1) = comps*(ents(i1:i2-1)-1)+i;
          i1 = i2;
       end
    end

    function ents = getEntities(obj, identifier)
       % notempty varargin ---> return entities with component multiplication
       %this method return the index of constrained entities inside solution vectors
       %needed for applying Dirichlet BCs
       %(pressure, displacement)...
       ents = obj.getData(identifier).data.entities;
    end
    
     function ents = getLoadedEntities(obj, identifier)
        % return loaded entities for Volume or Surface BCs
        % if varargin not empty, return loaded ents with component
        % multiplication (according to specified direction)
        ents = obj.getData(identifier).loadedEnts;
     end

     function ents = getCompLoadedEntities(obj, identifier)
        % return loaded entities for Volume or Surface BCs
        % if varargin not empty, return loaded ents with component
        % multiplication (according to specified direction)
        ents = obj.getData(identifier).loadedEnts;
        if strcmp(obj.getType(identifier),'Dir')
           ents = obj.getCompEntities(identifier,ents);
        else
           c = find(strcmp(["x","y","z"],obj.getDirection(identifier)));
           ents = 3*(ents-1)+c;
        end
     end

     % function getEntitiesDof
     % end

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