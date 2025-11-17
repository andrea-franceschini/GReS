classdef (Abstract) SolutionScheme < handle
  % General interface for creating a custom solution strategy

  properties
    solverStat
    %linearSolver -> for the future
    simparams
    domains
    interfaces
  end


  methods (Access = public)
    function obj = SolutionScheme(model,simparams,varargin)

      % input: additional input for the solution strategy
      obj.setSolutionScheme(model,simparams,varargin{:});
    end

    function setSolutionScheme(obj,model,simparams,varargin)
      obj.simparams = simparams;
      obj.domains = model.domains;
      obj.interfaces  = model.interfaces;
    end

  end


  methods (Abstract)
    flConv = solveTimeStep(obj,t,dt);
  end


end

