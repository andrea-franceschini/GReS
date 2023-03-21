classdef Boundaries < handle
  % BOUNDARY CONDITIONS - General boundary conditions class

  properties (Access = public)
    % Creation of a Map object for the boundary conditions
    db = containers.Map('KeyType','char','ValueType','any')
    %
    sizeM
    BCName
    resState
  end

  methods (Access = public)
    % Class constructor method
    function obj = Boundaries(fileName)
      % Calling the function to read input data from file
      obj.readInputFile(fileName);
    end
    
    function delete(obj)
      % Close all BC files
      for i=1:length(obj.BCName)
        cond = getBC(obj,obj.BCName(i));
        fclose(cond.fID);
      end
      remove(obj.db,keys(obj.db));
      obj.db = [];
      obj.sizeM = [];
      obj.BCName = [];
    end

    % Check if the BCIdentifier defined by the user is a key of the Map object
    function bcType = getBC(obj,BCIdentifier)
      if (obj.db.isKey(lower(BCIdentifier)))
        bcType = obj.db(lower(BCIdentifier));
      else
        % Displaying error message if the identifier does not refer
        % to any existing class
        error('Boundary condition %s not present', BCIdentifier);
      end
    end
    
    function iniBC(obj,BCList,state)
      % Initialize BC
      l = length(BCList);
      assert(l>0,'Warning: No boundary conditions will be applied.');
      obj.sizeM = max(length(state.displ),length(state.pressure));
      obj.BCName = BCList;
      obj.resState = state;
    end
    
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
    function applyDirVal(obj,t)
      % Apply Dirichlet conditions to the solution vector
      for i=1:length(obj.BCName)
        cond = getBC(obj,obj.BCName(i));
        %
        % Check if the BC need to be updated (we include this check here to
        % limit the number of times the db is accessed)
        %
        Boundaries.checkAndUpdateCond(cond,t);
        if strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
          % Interpolate in the step interval
          % fac = (t-t_old)/(t_new-t_old)
          fac = (t - cond.timeInt(1,cond.timeInt(2,1)))/ ...
            (cond.timeInt(1,cond.timeInt(2,2)) - ...
             cond.timeInt(1,cond.timeInt(2,1)));
          tmpVec = (1-fac)*cond.boundVal(:,cond.timeInt(2,1)) + ...
            fac*cond.boundVal(:,cond.timeInt(2,2));
          if strcmp(cond.boundPhysics,'flow')
            obj.resState.pressure(cond.boundDof) = tmpVec;
          elseif strcmp(cond.boundPhysics,'poro')
            obj.resState.displ(cond.boundDof) = tmpVec;
          end
        end
      end
    end
    
    function applyBC(obj,syst,t)
      % Impose BC to the linearized system (Jacobian matrix + RHS)
      % The Penalty method is used for the Dirichlet conditions
      %
%       l = length(obj.BCName);
%       assert(l>0,'Warning: No boundary conditions will be applied.');
%       if strcmp(probType,'lin'); fac = 1; else; fac = -1; end
      maxVal = max(abs(syst.K), [], 'all');
      %
      for i=1:length(obj.BCName)
        cond = getBC(obj,obj.BCName(i));
        % fac = (t-t_old)/(t_new-t_old)
        fac = (t - cond.timeInt(1,cond.timeInt(2,1)))/ ...
          (cond.timeInt(1,cond.timeInt(2,2)) - ...
           cond.timeInt(1,cond.timeInt(2,1)));
        tmpVec = (1-fac)*cond.boundVal(:,cond.timeInt(2,1)) + ...
          fac*cond.boundVal(:,cond.timeInt(2,2));
        if strcmp(cond.boundType,'neu')  % Apply Neumann conditions,if any
          syst.rhs(cond.boundDof) = syst.rhs(cond.boundDof) - tmpVec;
        elseif strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
          syst.K(obj.sizeM*(cond.boundDof-1)+cond.boundDof) = maxVal*10^10;
          syst.rhs(cond.boundDof) = 0;
%           if strcmp(probType,'lin')
%             syst.rhs(cond.boundDof) = cond.boundVal*(maxVal*10^10);
%           else
%             varargin{1}.displ(cond.boundDof) = cond.boundVal;
%             syst.rhs(cond.boundDof) = 0;
%           end
        end
      end
    end
    
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
      nBC = length(fileName);
      assert(nBC > 0,'No boundary conditions are to be imposed');
%       BCType = lower(BCType);
%       if ~ismember(BCType,["dir", "neu"])
%         error('%s boundary condition is not admitted\n Accepted types are: dir -> Dirichlet, neu -> Neumann',BCType);
%       end
      %
      for i=1:nBC
        fid = fopen(fileName(i), 'r');
        if fid == -1
          error('File %s not opened correctly',fileName(i));
        end
        block = '';
        [flEof,line] = Boundaries.readLine(fid);
        % Reading the BC
        while (~isempty(line))
          if flEof == 1
            error('End of file while reading BC in file %s',fileName);
          end
          line = strtrim(strtok(line,'%'));
          block = [block, line, ' '];
          [flEof,line] = Boundaries.readLine(fid);
        end
        blockSplt = strsplit(string(deblank(block)));
        % Calling the specific BC class based on the BC name
        %
        % Structure of blockSplt (see also setBC):
        %       blockSplt(1) -> label defining the discrete item on which 
        %                       the BC is applied
        %       blockSplt(2) -> type of BC (Dir = Dirichlet, Neu = Neumann)
        %       blockSplt(3) -> physical process for which the BC is
        %                       specified (Poro = poromechanics, 1PFlow = 
        %                       single-phase flow)
        %       blockSplt(4) -> BC name
        if ~all(isnan(str2double(blockSplt(1:4))))
          error('The first four entries in %s must be strings',fileName);
        end
        blockSplt(1:4) = lower(blockSplt(1:4));
        switch blockSplt(1)
          case 'nodebc'
            if obj.db.isKey(blockSplt(4))
              error('%s BC already defined',blockSplt(4));
            else
              obj.db(blockSplt(4)) = NodeBC(fid,blockSplt);
            end
          otherwise
            error('Boundary condition %s not available', blockSplt(1));
        end
        cond = getBC(obj,blockSplt(4));
        % Read the BC values at t=0 and at the end of the event
        Boundaries.readStep(cond,1);
        Boundaries.readStep(cond,2);
%         fclose(fid);
      end
    end
  end
  
  methods(Static = true)
    function checkAndUpdateCond(cond,t)
      % Check if we move to the next event
      while max(cond.timeInt(1,:)) < t
        if cond.timeInt(1,1) < cond.timeInt(1,2)
%         pos = 1;
          cond.timeInt(2,1) = 2;
          cond.timeInt(2,2) = 1;
        else
%         pos = 2;
          cond.timeInt(2,1) = 1;
          cond.timeInt(2,2) = 2;
        end
        Boundaries.readStep(cond,cond.timeInt(2,2));
      end
    end
 
    function readStep(cond,pos)
      % Read the BC values in the next event
      try
      cond.timeInt(1,pos) = fscanf(cond.fID,['TIME' '%e'],[1 1]);
      catch
        error('Wrong format of the TIME header in %s %s BC at time %f', ...
          cond.boundPhysics,cond.boundType,cond.timeInt(1,pos))
      end
      try
        cond.boundVal(:,pos) = fscanf(cond.fID,'%e\n',[1 length(cond.boundDof)]);
      catch
        error('Wrong number of values in %s %s BC at time %f', ...
        cond.boundPhysics,cond.boundType,cond.timeInt(1,pos))
      end
%       if ~feof(obj.fID); tmpLine = fgetl(obj.fID); end
    end
    
    % Read the next line and check for eof
    function [flEof,line] = readLine(fid)
      flEof = feof(fid);   % end-of-file flag
      if flEof == 1
        line = '';
      else
        line = deblank(fgetl(fid));
%         if isempty(line)
%           error('No blank lines are admitted in BC file')
%         end
      end  
    end
  end
end