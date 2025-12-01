classdef Sedimentation < PhysicsSolver
  % Sedimentation model subclass
  % Coupled Poromechanics with SinglePhaseFlow

  properties
     structuredGrid
  end

  properties (Access = private)
    Solver
  end

  methods (Access = public)
    function obj = Sedimentation(domain)
      obj@PhysicsSolver(domain);
    end

    function registerSolver(obj,input)
      % setup the solver with custom input
      % str2num(input.domain.division)
      switch input.physics
        case "VariablySaturatedFlow"
          obj.Solver = VariablySaturatedFlow(obj.domain);
        otherwise
          obj.Solver = SinglePhaseFlowFVTPFA(obj.domain);
      end
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

    function out = getField()
      out = SinglePhaseFlow.getField();
    end

    function out = isSymmetric()
      out = obj.Solver.isSymmetric();
    end

  end

end


