classdef Boundaries < handle
  % BOUNDARY CONDITIONS - General boundary conditions class

  properties (Access = public)
    % Creation of a Map object for the boundary conditions
    db
  end

  methods (Access = public)
    % Class constructor method
    function obj = Boundaries(fileName)
      % MATLAB evaluates the assignment expression for each instance, which
      % ensures that each instance has a unique value
      obj.db = containers.Map('KeyType','char','ValueType','any');
      % Calling the function to read input data from file
      obj.readInputFiles(fileName);
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
      list = obj.getData(identifier).data.entities;
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
    
    function setDofs(obj, identifier, list)
      obj.getData(identifier).data.entities = list;
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
      if (~ismember(token, ['NodeBC', 'SurfBC', 'VolForce']))
        error(['%s condition is not admitted\n', ...
          'Accepted types are: NodeBC   -> Boundary cond. on nodes\n',...
          '                    SurfBC   -> Boundary cond. on surfaces\n',...
          '                    VolForce -> Volume force on elements'], token);
      end
      if ismember(token, ['NodeBC', 'SurfBC'])
        type = Boundaries.readToken(fid);
        if (~ismember(type, ['Dir', 'Neu']))
          error(['%s boundary condition is not admitted\n', ...
            'Accepted types are: Dir -> Dirichlet, Neu -> Neumann'], type);
        end
      end
      physics = Boundaries.readToken(fid);
      name = Boundaries.readToken(fid);
      setFile = Boundaries.readToken(fid);
      [times, dataFiles] = Boundaries.readDataFiles(fid);
      fclose(fid);
      if obj.db.isKey(name)
        error('%s boundary condition name already defined', name);
      end
      switch (token)
        case {'NodeBC', 'SurfBC'}
          obj.db(name) = struct('data', BoundaryEntities(name, setFile, times, dataFiles), ...
            'cond',token,'type', type, 'physics', physics);
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