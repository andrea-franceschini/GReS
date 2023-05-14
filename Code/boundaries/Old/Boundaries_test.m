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
%       disp('I am deleting')
      remove(obj.db,keys(obj.db));
      obj.db = [];
      obj.sizeM = [];
      obj.BCName = [];
      obj.resState = [];
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
    function applyDirVal(obj)
      for i=1:length(obj.BCName)
        cond = getBC(obj,obj.BCName(i));
        %
        if strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
          if strcmp(cond.boundPhysics,'flow')
            obj.resState.pressure(cond.boundDof) = cond.boundVal;
          elseif strcmp(cond.boundPhysics,'poro')
            obj.resState.displ(cond.boundDof) = cond.boundVal;
          end
        end
      end
    end
    
    function applyBC(obj,syst)
%       l = length(obj.BCName);
%       assert(l>0,'Warning: No boundary conditions will be applied.');
%       if strcmp(probType,'lin'); fac = 1; else; fac = -1; end
      maxVal = max(abs(syst.K), [], 'all');
      %
      for i=1:length(obj.BCName)
        cond = getBC(obj,obj.BCName(i));
        if strcmp(cond.boundType,'neu')  % Apply Neumann conditions,if any
          syst.rhs(cond.boundDof) = syst.rhs(cond.boundDof) - cond.boundVal;
        elseif strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
          syst.K(obj.sizeM*(cond.boundDof-1)+cond.boundDof) = maxVal*10^10;
          syst.rhs(cond.boundDof) = maxVal*10^10*cond.boundVal;
%           if strcmp(probType,'lin')
%             syst.rhs(cond.boundDof) = cond.boundVal*(maxVal*10^10);
%           else
%             varargin{1}.displ(cond.boundDof) = cond.boundVal;
%             syst.rhs(cond.boundDof) = 0;
%           end
        end
      end
    end
    
    function applyBCNeu(obj,syst)
%       l = length(obj.BCName);
%       assert(l>0,'Warning: No boundary conditions will be applied.');
%       if strcmp(probType,'lin'); fac = 1; else; fac = -1; end
%       maxVal = max(abs(syst.K), [], 'all');
      %
      for i=1:length(obj.BCName)
        cond = getBC(obj,obj.BCName(i));
        if strcmp(cond.boundType,'neu')  % Apply Neumann conditions,if any
          syst.rhs(cond.boundDof) = syst.rhs(cond.boundDof) - cond.boundVal;
        end
        %
%         if strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
%           syst.K(obj.sizeM*(cond.boundDof-1)+cond.boundDof) = maxVal*10^10;
%           syst.rhs(cond.boundDof) = 0;
%           if strcmp(probType,'lin')
%             syst.rhs(cond.boundDof) = cond.boundVal*(maxVal*10^10);
%           else
%             varargin{1}.displ(cond.boundDof) = cond.boundVal;
%             syst.rhs(cond.boundDof) = 0;
%           end
%         end
      end
    end
    
    function applyBCDir(obj,syst)
%       l = length(obj.BCName);
%       assert(l>0,'Warning: No boundary conditions will be applied.');
%       if strcmp(probType,'lin'); fac = 1; else; fac = -1; end
      maxVal = max(abs(syst.K), [], 'all');
      %
      for i=1:length(obj.BCName)
        cond = getBC(obj,obj.BCName(i));
%         if strcmp(cond.boundType,'neu')  % Apply Neumann conditions,if any
%           syst.rhs(cond.boundDof) = syst.rhs(cond.boundDof) - cond.boundVal;
%         end
        %
        if strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
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
    
  end

  methods (Access = private)
    % Reading boundary input file
    function readInputFile(obj,fileName)
      nBC = length(fileName);
      assert(nBC > 0,'No boundary conditions are imposed');
%       BCType = lower(BCType);
%       if ~ismember(BCType,["dir", "neu"])
%         error('%s boundary condition is not admitted\n Accepted types are: dir -> Dirichlet, neu -> Neumann',BCType);
%       end
      %
      for i=1:nBC
        fid = fopen(fileName(i), 'r');
        block = '';
        [flEof,line] = Boundaries.readLine(fid);
        % Reading the BC
        while (~strcmpi(line,'End'))
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
              obj.db(blockSplt(4)) = NodeBC(blockSplt);
            end
          otherwise
            error('Boundary condition %s not available', blockSplt(1));
        end
        fclose(fid);
      end
    end
  end

  methods(Static = true)
    % Read the next line and check for eof
    function [flEof,line] = readLine(fid)
      flEof = feof(fid);   % end-of-file flag
      if flEof == 1
        line = '';
      else
        line = deblank(fgetl(fid));
        if isempty(line)
          error('No blank lines are admitted in BC file')
        end
      end  
    end
  end
end