classdef FCSolver < handle
  % Built in fully coupled solver
  % All equations are solved once using Newton Raphson

  properties (Access = private)
    %
    simparams
    domain
    %
    t = 0
    tStep = 0
    iter
    dt
  end

  properties (Access = public)
     solStatistics
     % to solve
     linsolver
  end

  methods (Access = public)
    function obj = FCSolver(simparams,domain,varargin)
      obj.setNonLinearSolver(simparams,domain);
      saveStasticts = false(3,1);
      for k = 1:(nargin-1)/2
        pos = 2*(k-1)+1;
        switch varargin{pos}
          case 'SaveRelError'
            saveStasticts(1) = varargin{pos+1};
          case 'SaveAbsError'
            saveStasticts(2) = varargin{pos+1};
          case 'SaveBStepInf'
            saveStasticts(3) = varargin{pos+1};
          case 'SaveConverg'
            saveStasticts(:) = varargin{pos+1};
          otherwise
        end
      end

      obj.solStatistics = SolverStatistics(obj.simparams.itMaxNR,...
                                           obj.simparams.relTol,...
                                           obj.simparams.absTol,...
                                           saveStasticts);

      % Check if there is manual input from the user, if not use defaults
      start_dir = pwd;
      chronos_xml = fullfile(start_dir,'linsolver.xml');
      if(isfile(chronos_xml))
         obj.linsolver = linearSolver(obj.domain,[],chronos_xml);
      else
         if gresLog().getVerbosity > 2
            fprintf('Using default values for linsolver\n');
         end
         obj.linsolver = linearSolver(obj.domain,[]);
      end
    end

    function [simStat] = NonLinearLoop(obj)
      simStat = 1;
      % Initialize the time step increment
      obj.dt = obj.simparams.dtIni;

      % initialize the state object
      applyDirVal(obj.domain,obj.t);

      %state = getState(obj.domain);


      % Loop over time
      while obj.t < obj.simparams.tMax

        absTol = obj.simparams.absTol;
        residual = zeros(obj.simparams.itMaxNR+1,2);

        % Update the simulation time and time step ID
        obj.tStep = obj.tStep + 1;
        obj.t = obj.t + obj.dt;

        % Apply the Dirichlet condition value to the solution vector
        applyDirVal(obj.domain,obj.t);
        %
        gresLog().log(0,'\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
        gresLog().log(0,'-----------------------------------------------------------\n');
        gresLog().log(1,'Iter     ||rhs||     ||rhs||/||rhs_0||\n');

        % Compute Rhs and matrices of NonLinear models
        % Loop over available linear models and compute the jacobian
        assembleSystem(obj.domain,obj.dt);

        % Apply BCs to the blocks of the linear system
        applyBC(obj.domain,obj.t);

        % compute Rhs norm
        rhs = getRhs(obj.domain);
        rhsNorm = norm(cell2matrix(rhs),2);
        rhsNormIt0 = rhsNorm;
        residual(1,1) = rhsNormIt0;
        residual(1,2) = 1.;

        % consider output of local field rhs contribution
        tolWeigh = obj.simparams.relTol*rhsNorm;
        obj.iter = 0;
        %
        gresLog().log(1,'0     %e     %e\n',rhsNorm,rhsNorm/rhsNormIt0);

        while ((rhsNorm > tolWeigh) && (obj.iter < obj.simparams.itMaxNR) ...
            && (rhsNorm > absTol)) || obj.iter == 0
          obj.iter = obj.iter + 1;

          % Get system Jacobian
          J = getJacobian(obj.domain);

          % Solve linear system
          du = FCSolver.solve(obj,J,rhs);

          % Update current model state
          updateState(obj.domain,du);

          % Compute Rhs and Matrices of NonLinear models
          assembleSystem(obj.domain,obj.dt);

          % Apply BCs to the blocks of the linear system
          applyBC(obj.domain,obj.t);

          rhs = getRhs(obj.domain);
          % compute Rhs norm
          rhsNorm = norm(cell2matrix(rhs),2);
          residual(obj.iter+1,1) = rhsNorm;
          residual(obj.iter+1,2) = rhsNorm/rhsNormIt0;

          gresLog().log(1,'%d     %e     %e\n',obj.iter,residual(obj.iter+1,1),residual(obj.iter+1,2));
        end
        %

        % Check for convergence
        flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);

        % save solution statistics
        if flConv
          obj.solStatistics.saveIt(obj.t,residual(1:obj.iter+1,1),residual(1:obj.iter+1,2));
        else
          % save backstep info
          obj.solStatistics.saveBackIt();
        end

        % Manage next time step
        manageNextTimeStep(obj,flConv);
      end

    end
  end

  methods (Access = private)
    function setNonLinearSolver(obj,simparams,domain)

      obj.simparams = simparams;
      assert(isscalar(domain),"FCSolver works only with one single domain")

      obj.domain = domain;
      obj.domain.stateOld = copy(domain.getState());
      obj.domain.simparams = simparams;

      assert(~isempty(obj.domain.solverNames),"Missing PhysicsSolvers" + ...
        " in Discretizer. Use the method addPhysicsSolver()")

    end


    % function [t, dt] = updateTime(obj,conv,dt)
    % 
    %   if obj.domain.outstate.modTime
    %     tmp = find(obj.t<obj.domain.outstate.timeList(),1,'first');
    %     if ~conv
    %       t = min([obj.t + obj.dt, obj.t + dt, obj.domain.outstate.timeList(tmp)]);
    %     else
    %       t = min([obj.t + obj.dt, obj.domain.outstate.timeList(tmp)]);
    %     end
    %   else
    %     t = obj.t + obj.dt;
    %   end
    %   dt = t - obj.t;
    % end

    function manageNextTimeStep(obj,flConv)

      if flConv % Convergence

        obj.domain.state.t = obj.t;

        % print results at current time
        printState(obj.domain);

        % advance the domain state after convergence
        % here copy current state to old state
        advanceState(obj.domain);

        % go to next time step
        tmpVec = obj.simparams.multFac;
        obj.dt = min([obj.dt * min(tmpVec), obj.simparams.dtMax]);
        obj.dt = max([obj.dt obj.simparams.dtMin]);


        % limit time step to end of simulation time
        if ((obj.t + obj.dt) > obj.simparams.tMax)
          obj.dt = obj.simparams.tMax - obj.t;
        end

      else

        % backstep
        goBackState(obj.domain);

        obj.t = obj.t - obj.dt;
        obj.tStep = obj.tStep - 1;

        obj.dt = obj.dt/obj.simparams.divFac;

        if obj.dt < obj.simparams.dtMin
          error('Minimum time step reached')
        else 
          gresLog().log(0,'\n %s \n','BACKSTEP')
        end
      end

    end

  end

  methods (Static)
    function sol = solve(obj,J,rhs)
      rhs = cell2matrix(rhs);

      % Actual solution of the system
      [sol,~] = obj.linsolver.Solve(J,-rhs,obj.t);

    end
  end
end
