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


  methods (Access = public)

    function fSSplitConv = solveStep(obj,varargin)

      % set target variable for each domain and interface
      if ~isempty(varargin)
        error("solveStep with target variables is not yet implemented")
      end

      dom = obj.domains;
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
        currPress = dom.state.data.pressure;
        prevPress = physSolv.conv.pressure;
        relPressChange = norm(currPress-prevPress)/norm(prevPress);

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

      % Check that we have an even number of inputs
      if mod(length(varargin), 2) ~= 0
        error('Arguments must come in key-value pairs.');
      end

      % Loop through the key-value pairs
      for k = 1:2:length(varargin)
        key = varargin{k};
        value = varargin{k+1};

        if isempty(value)
          continue
        end

        if ~ischar(key) && ~isstring(key)
          error('Keys must be strings');
        end

        switch lower(key)
          % case 'simulationparameters'
          %   assert(isa(value, 'SimulationParameters')|| isempty(value),msg)
          %   obj.simparams = value;
          case 'simulationparameters'
            obj.simparams = value;
          case 'output'
            obj.output = value;
          case {'domain','domains'}
            obj.domains = value;
          case 'maxiterations'
            obj.maxIterFSS = value;
          case 'reltolerance'
            obj.tolFSS = value;
          otherwise
            error('Unknown input key %s for SolutionScheme \n', key);
        end
      end

      obj.nDom = numel(obj.domains);
      obj.nInterf = numel(obj.interfaces);

      assert(~isempty(obj.simparams),"Input 'simulationParameters'" + ...
        " is required for SolutionScheme")
      assert(obj.nDom == 1,"Only one domain is admitted when using fixedStressSplit");
      assert(obj.nInterf == 0,"FixedStressSplit with internal interfaces is not yet implmented");

      obj.nVars = 0;

      for i = 1:obj.nDom
        obj.domains(i).domainId = i;
        obj.domains(i).simparams = obj.simparams;
        obj.domains(i).outstate = obj.output;
        obj.domains(i).stateOld = copy(obj.domains(i).getState());
        obj.nVars = obj.nVars + obj.domains(i).dofm.getNumberOfVariables();
      end

      assert(obj.domains.solverNames == "BiotFixedStressSplit",...
        "FixedStressSplit algorithm only requires 'BiotFixedStressSplit' physics solver")

    end


    function newtonConv = nonLinearSolve(obj,varName)

      % nonlinear loop for single physics model
      % consider replacing this with a call to
      % NonLinearImplicit.solveStep(varName)

      setLinearSolver(obj,[],varName);

      
      varId = obj.domains.dofm.getVariableId(varName);
      dom = obj.domains(1);
      physSolv = dom.getPhysicsSolver("BiotFixedStressSplit");

      physSolv.assembleSystem(obj.dt,varName);

      dom.applyBC(obj.t,varName);

      gresLog().log(1,'Iter     ||rhs||     ||rhs||/||rhs_0||\n');

      rhs = dom.rhs{varId};
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

        J = dom.J(varId,varId);
        rhs = dom.rhs(varId);

        % solve linear system
        du = solve(obj,J,rhs);

        % update state variable calling directly the BiotFixedStressSplit
        % solver
        physSolv.updateState(du,varName);

        % reassemble system
        physSolv.assembleSystem(obj.dt,varName);
        obj.domains.applyBC(obj.t,varName);

        rhs = dom.rhs{varId};
        rhsNorm = norm(rhs,2);

        gresLog().log(1,'%d     %e     %e\n',iter,rhsNorm,rhsNorm/rhsNormIt0);

        % Check for convergence
        newtonConv = (rhsNorm < tolWeigh || rhsNorm < obj.simparams.absTol);

      end % end newton loop
    end



    function setLinearSolver(obj,xmlInput,varName)

      if nargin == 1
        xmlInput = [];
      end
      if nargin <= 2
        obj.solverFlow = linearSolver(obj,xmlInput,"pressure");
        obj.solverMech = linearSolver(obj,xmlInput,"displacements");
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

