classdef NonLinearImplicit < SolutionScheme
  % Class for solving non linear implicit problems with changing
  % configuration

  % The time loop is implemented in the base class SolutionScheme

  properties (SetAccess = protected,GetAccess=public)
    iterNL = 0          % nonlinear iteration number
    iterConfig = 0      % configuration iteration number
    targetVariables     % variables currently solved for
    totIter = 0
  end



  methods (Access = public)
    function obj = NonLinearImplicit(varargin)

      obj@SolutionScheme(varargin{:});

    end

    %

    function converged = solveStep(obj,varargin)

      % set target variable for each domain and interface
      obj.targetVariables = [varargin{:}];
      if ~isempty(varargin)
        error("solveStep with target variables is not yet implemented")
      end

      hasConfigurationChanged = true;
      newtonConv = false;
      absTol = obj.simparams.absTol;
      obj.iterConfig = 0;
      obj.iterNL = 0;

      while (hasConfigurationChanged) && (obj.iterConfig < obj.simparams.itMaxConfig)

        gresLog().log(0,'\nConfiguration iteration n. %i \n', obj.iterConfig);
        obj.iterConfig = obj.iterConfig + 1;

        obj.applyDirVal();

        gresLog().log(1,'Iter     ||res||     ||res||/||res_0||\n');

        obj.model.assembleSystem(obj.dt);
        obj.applyBC();
        rhsNorm = norm(cell2mat(obj.model.res),2);
        obj.rhsNormIt0 = rhsNorm;
        rhsNormIt0 = rhsNorm;

        tolWeigh = obj.simparams.relTol*rhsNorm;

        gresLog().log(1,'0     %e     %e\n',rhsNorm,rhsNorm/rhsNormIt0);

        newtonConv = false;

        % reset non linear iteration counter
        newtonIter = 0;

        %%% NEWTON LOOP %%%
        while (~newtonConv) && (newtonIter < obj.simparams.itMaxNR)

          newtonIter = newtonIter + 1;
          obj.iterNL = obj.iterNL + 1;
          obj.totIter = obj.totIter + 1;

          du = solve(obj,obj.model.J,obj.model.res,newtonIter);
          obj.model.updateState(du);

          obj.model.assembleSystem(obj.dt);
          obj.applyBC();
          rhsNorm = norm(cell2mat(obj.model.res),2);

          gresLog().log(1,'%d     %e     %e\n',newtonIter,rhsNorm,rhsNorm/rhsNormIt0);

          % Check for convergence
          newtonConv = (rhsNorm < tolWeigh || rhsNorm < absTol);

        end % end newton loop

        if newtonConv % Newton Convergence

          hasConfigurationChanged = obj.model.updateConfiguration();

        else

          break

        end

      end

      converged = newtonConv && ~hasConfigurationChanged;

    end
  end




  methods (Access = protected)


    function setLinearSolver(obj,varargin)

      % The Chronos preconditioners still require scheme-owned domains.
      options = struct('useMatlab',true);
      context = struct('simparams',struct('linSolverParams',options));
      obj.linsolver = linearSolver(context,[]);

    end


    function out = computeRhsNorm(obj)

      out = norm(cell2mat(obj.model.res),obj.simparams.pNorm);

    end
  end
end
