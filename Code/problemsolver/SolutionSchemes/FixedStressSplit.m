classdef FixedStressSplit < SolutionScheme
  
% Fixed-stress split algorithm
% This scheme is valid only when used together with the
% BiotFixedStressSplit physics solver.

  properties
    % instance of linearSolver for singlePhysics
    solverMech
    solverFlow
    maxIterFSS = 10;
    tolFSS = 1e-5;
    iterFSS
    iters
  end


  % Chronos still reads these from the solution scheme.
  properties (Dependent, SetAccess = private)
    domains
    nDom
    nInterf
  end

  properties (SetAccess = private)
    iterConfig = 1
  end


  methods
    function obj = FixedStressSplit(varargin)
      obj@SolutionScheme(varargin{:});
    end

    function dom = get.domains(obj)
      dom = obj.model.domains;
    end

    function n = get.nDom(obj)
      n = numel(obj.model.domains);
    end

    function n = get.nInterf(obj)
      n = numel(obj.model.interfaces);
    end


    function fSSplitConv = solveStep(obj,varargin)

      % set target variable for each domain and interface
      if ~isempty(varargin)
        error("solveStep with target variables is not yet implemented")
      end

      dom = obj.model.domains;
      physSolv = dom.getPhysicsSolver("BiotFixedStressSplit");

      % apply dirichlet value at the beginning of the time step
      dom.applyDirVal(obj.t);

      obj.iterFSS = 0;

      fSSplitConv = false;

      % critical: initialize the converged state to match the current state
      physSolv.conv.pressure = getState(dom,"pressure");
      physSolv.conv.displacements = getState(dom,"displacements");


      while (~fSSplitConv) && (obj.iterFSS < obj.maxIterFSS)

        gresLog().log(0,'\nFixed Stress Split iteration n. %i \n', obj.iterFSS);
        obj.iterFSS = obj.iterFSS + 1;

        % solve pressure
        gresLog().log(1,"\nSolving pressure\n")
        newtConvFlow = obj.nonLinearSolve("pressure");

        if ~newtConvFlow, break, end

        % check convergence of fixed stress split
        currPress = getState(physSolv,"pressure");
        prevPress = physSolv.conv.pressure;
        relPressChange = norm(currPress-prevPress);
        if relPressChange > 0
          relPressChange = relPressChange/norm(prevPress);
        end

        physSolv.advanceState("pressure");

        % solve mechanics
        gresLog().log(1,"\nSolving mechanics\n")
        newtConvMech = obj.nonLinearSolve("displacements");

        if ~newtConvMech, break, end

        physSolv.advanceState("displacements");

        % check convergence of the scheme
        gresLog().log(1,"\nRelative pressure change = %1.4e\n",relPressChange)
        if relPressChange < obj.tolFSS
          fSSplitConv = true;
          obj.iters(obj.tStep) = obj.iterFSS;
        end
      end
    end


    function iters = getFixedStressIters(obj)
      iters = obj.iters(1:obj.tStep);
    end

  end




  methods (Access = protected)

    function setSolutionScheme(obj,varargin)

      setSolutionScheme@SolutionScheme(obj,varargin{:})

      % barrier for this solver
      assert(obj.nDom == 1,"Only one domain is admitted when using fixedStressSplit");
      assert(isequal(obj.domains.solverNames,"BiotFixedStressSplit"),...
        "FixedStressSplit algorithm only requires 'BiotFixedStressSplit' physics solver");
      assert(obj.nInterf == 0,"FixedStressSplit with internal interfaces is not yet implemented");

      default = struct('maxiterations',10,'reltolerance',1e-5);
      parm = readInput(default,varargin{:});
      obj.maxIterFSS = parm.maxiterations;
      obj.tolFSS = parm.reltolerance;

    end


    function newtonConv = nonLinearSolve(obj,varName)

      % nonlinear loop for single physics model
      % consider replacing this with a call to
      % NonLinearImplicit.solveStep(varName)

      setLinearSolver(obj,varName);
      
      model = obj.model;
      dom = model.domains(1);
      physSolv = dom.getPhysicsSolver("BiotFixedStressSplit");

      physSolv.assembleSystem(obj.dt,varName);

      dom.applyBC(obj.t,varName);

      gresLog().log(1,'Iter     ||rhs||     ||rhs||/||rhs_0||\n');

      rhs = model.getResidual(dom,varName);
      rhsNorm = norm(rhs,2);
      rhsNormIt0 = rhsNorm;

      tolWeigh = obj.simparams.relTol*rhsNorm;

      gresLog().log(1,'0     %e     %e\n',rhsNorm,rhsNorm/rhsNormIt0);

      newtonConv = (rhsNorm < tolWeigh || rhsNorm < obj.simparams.absTol);

      % reset non linear iteration counter
      iter = 0;

      %%% NEWTON LOOP FOR FLOW %%%
      while (~newtonConv) && (iter < obj.simparams.itMaxNR)

        iter= iter + 1;

        % Extract only the current subproblem; its update is field-local.
        J = model.getJacobian(dom,varName);
        du = solve(obj,J,{rhs},iter);

        % update state variable calling directly the BiotFixedStressSplit
        % solver
        physSolv.updateState(du,varName);

        % reassemble system
        physSolv.assembleSystem(obj.dt,varName);
        dom.applyBC(obj.t,varName);

        rhs = model.getResidual(dom,varName);
        rhsNorm = norm(rhs,2);

        gresLog().log(1,'%d     %e     %e\n',iter,rhsNorm,rhsNorm/rhsNormIt0);

        % Check for convergence
        newtonConv = (rhsNorm < tolWeigh || rhsNorm < obj.simparams.absTol);

      end % end newton loop
    end


    

    function setLinearSolver(obj,varName)

      if nargin <= 1
        obj.solverFlow = linearSolver(obj,"pressure");
        obj.solverMech = linearSolver(obj,"displacements");
      else
        if strcmp(varName,"pressure")
          obj.linsolver = obj.solverFlow;
        elseif strcmp(varName,"displacements")
          obj.linsolver = obj.solverMech;
        else
          error("Invalid variable name for linear solver in " + ...
            "FixedStressSplit algorithm")
        end

      end

    end
  end
end
