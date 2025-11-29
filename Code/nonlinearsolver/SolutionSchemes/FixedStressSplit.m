classdef  FullyCoupledGeneral < SolutionScheme
  % Fully coupled fully implicit newton solver

  % This is a general method working with any number of domains and
  % interfaces

  methods (Access = public)
    function obj = FullyCoupledGeneral(simparams)
      obj.setSolutionScheme(simparams);
    end

  end


  methods 
    function flConv = solveTimeStep(obj,t,dt)
      % solve linear system with newton raphson strategy

    end
  end


end

