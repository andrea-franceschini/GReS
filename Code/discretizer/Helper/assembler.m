classdef assembler < handle
  % General class for assembling sparse matrices in FEM

  properties
    computeLocal          % function handle of local matrix computation
    iVec                  % array of row indices for sparse assembly
    jVec                  % array of column indices for sparse assembly
    valsVec               % array of values for sparse assembly
    count = 0             % counter for sparse assembly
    Nrows                 % number of rows 
    Ncols                 % number of columns 
  end

  methods
    function obj = assembler(n,nrows,ncols,varargin)
      [obj.iVec,obj.jVec,obj.valsVec] = deal(zeros(n,1));
      
      if ~isempty(nrows) && ~isempty(ncols)
        obj.Nrows = nrows;
        obj.Ncols = ncols;
      end

      if ~isempty(varargin)
        obj.computeLocal = varargin{1};
      else
        obj.computeLocal = @(dofR,dofC,mat)...
        assembler.computeLocalBase(dofR,dofC,mat);
      end
    end
    
    function varargout = localAssembly(obj,varargin)
      % call back to local matrix computation
      nOut = nargout;

      % Call local assembler and manage output
      allOut = cell(3 + nOut, 1);
      [allOut{:}] = obj.computeLocal(varargin{:});

      assert(numel(allOut) > 2,['At least 3 output are required: \n' ...
        'out(1): list of local row dofs \n' ...
        'out(2): list of local column dofs \n' ...
        'out(3): list of local values.']);

      dofr      = allOut{1};
      dofc      = allOut{2};
      localMat  = allOut{3};
      varargout = allOut(4:end);

      n = numel(localMat);
      [jLoc,iLoc] = meshgrid(dofc,dofr);
      obj.iVec(obj.count+1:obj.count+n) = iLoc(:);
      obj.jVec(obj.count+1:obj.count+n) = jLoc(:);
      obj.valsVec(obj.count+1:obj.count+n) = localMat(:);
      obj.count = obj.count + n;
    end

    function mat = sparseAssembly(obj)
      % trim row col indices if needed
      if length(obj.iVec) > obj.count
        obj.iVec = obj.iVec(1:obj.count);
        obj.jVec = obj.jVec(1:obj.count);
        obj.valsVec = obj.valsVec(1:obj.count);
      end
      %
      if ~isempty(obj.Nrows)
        mat = sparse(obj.iVec,obj.jVec,obj.valsVec,obj.Nrows,obj.Ncols);
      else
        mat = sparse(obj.iVec,obj.jVec,obj.valsVec);
      end

    end
  end

  methods (Static)
    function [dofr,dofc,mat] = computeLocalBase(dofr,dofc,mat)
      % placeholder for internal assembly kernel when local matrix 
      % and degree of freedom are already available
    end

  end
end

