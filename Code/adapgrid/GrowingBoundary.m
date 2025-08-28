classdef GrowingBoundary < handle
  % STRUCTUREDGRID Class to control all information related to the
  % structured grid

  properties (Access = public)
    cond        % Boundary condition applied 

    % totEnts     % Total number of constrained entities    
    nEntities   % Number of constrained entities for each degree of freedom
    entities    % Indices of constrained degrees of freedom
    % nTimes      % Number of input times
    % times       % Set of input times
    % availVals   % Values of currently-stored boundary conditions
    % availSteps  % Time id of currently-stored boundary conditions

    npoints     % number of points used to interpolate the new boundaries
    refindex    % Index of the boundary  
    weight      % weight 
  end

  properties (Access = private)
  end

  methods(Access = public)
    function obj = GrowingBoundary(cond,data,ty)
      %ADAPGRID Construct an instance of this class
      % obj.totEnts = data.totEnts;
      obj.nEntities = data.nEntities;
      obj.entities = data.entities;
      % obj.nTimes = data.nTimes;
      % obj.times = data.times;
      % obj.availVals = data.availVals;
      % obj.availSteps = data.availSteps;

      % How the boundary condition grow.
      switch ty
        case 'Constant'
          obj.cond = cond;
          obj.npoints = 1;
          obj.refindex = zeros(obj.nEntities,1);
          obj.weight = ones(obj.nEntities,1);
          obj.refindex = 1:obj.nEntities;          
        % TODO Create other ways to define the values for the boundary
        % conditions in newly created cells.
      end
    end

  end

  methods (Access = private)
    % Function related to construct this class.
    function constructIJK(obj,face,cell,number)
      %CONSTRUCTIJK Construct the ijk instance for each cell in the grid
    end

  end

  methods (Static)
    % Some moviment in the grid
    function ijk_loc = cell_IJK_from_Ref(ijk_neig,dir)
      %IJK_LOC Function to find cell ijk to where to grow
    end

  end



end