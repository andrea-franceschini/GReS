classdef SinglePhysics < PhysicsSolver
  % Biot model subclass
  % Coupled Poromechanics with SinglePhaseFlow

  properties

  end

  properties (Access = private)
    Solver
  end

  methods (Access = public)
    function obj = SinglePhysics(domain)
      obj@PhysicsSolver(domain);
    end

    function registerSolver(obj,input)
      % setup the solver with custom input
      if ~isfield(input,"physics")
        error("physics field not declared, it is necessary to " + ...
          " declare one of the physics models for the simulation");
      end
      obj.Solver=feval(input.physics,obj.domain);
      registerSolver(obj.Solver,[]);
    end

    function assembleSystem(obj,dt)
      obj.Solver.assembleSystem(dt);
    end

    function applyBC(obj,bcId,t)
      obj.Solver.applyBC(bcId,t);
    end

    function applyDirVal(obj,bcId,t)
      obj.Solver.applyDirVal(bcId,t);
    end

    function updateState(obj,solution)
      obj.Solver.updateState(solution);
    end

    function [cellData,pointData] = writeVTK(obj,fac,t)
      [cellData,pointData] = obj.Solver.writeVTK(fac,t);
    end

    function writeMatFile(obj,fac,tID)
      obj.Solver.writeMatFile(fac,tID);
    end

    function out = isLinear(obj)
      out = obj.Solver.isLinear();
    end
  end

  methods (Static)


  end

end


