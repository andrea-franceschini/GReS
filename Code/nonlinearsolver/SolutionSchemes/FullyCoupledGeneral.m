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

      applyDirVal(obj);
      computeMatricesAndRhs(obj);
      applyBC(obj);
      rhs = assembleRhs(obj);
      rhsNorm = norm(cell2mat(rhs),2);

      tolWeigh = obj.simParameters.relTol*rhsNorm;
      obj.iter = 0;
      %
      
      if obj.simParameters.verbosity > 1
        fprintf('0     %e\n',rhsNorm);
      end

      while ((rhsNorm > tolWeigh) && (obj.iter < obj.simParameters.itMaxNR) ...
          && (rhsNorm > absTol)) || obj.iter == 0

        obj.iter = obj.iter + 1;

        J = assembleJacobian(obj);

        du = solve(obj,J,rhs);

        % update primary variables and multipliers
        updateState(obj,du);

        computeMatricesAndRhs(obj);
        applyBC(obj);
        rhs = assembleRhs(obj);
        rhsNorm = norm(cell2mat(rhs),2);


        if obj.simParameters.verbosity > 1
          fprintf('%d     %e\n',obj.iter,rhsNorm);
        end
      end
      %
      % Check for convergence
      flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);

      if flConv % Convergence
        % Advance state of non linear models
        for i = 1:obj.nDom
          obj.state(i).curr.t = obj.t;
          if isPoromechanics(obj.domains(i).model)
            obj.domains(i).getSolver('Poromechanics').advanceState();
          end
        end

        printState(obj);
      end

    end
  end


end

