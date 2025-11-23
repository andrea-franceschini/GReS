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
    % toGrow
  end

  properties (Access = public)
    solStatistics
  end

  methods (Access = public)
    function obj = FCSolver(simparams,linSyst,varargin)
      obj.setNonLinearSolver(simparams,linSyst);
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
    end

    function [simStat] = NonLinearLoop(obj)
      simStat = 1;
      % Initialize the time step increment
      obj.dt = obj.simparams.dtIni;
      delta_t = obj.dt; % dynamic time step

      % initialize the state object
      applyDirVal(obj.domain,obj.t);

      state = getState(obj.domain);
      stateOld = getStateOld(obj.domain);

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
        if obj.simparams.verbosity > 0
          fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,delta_t);
          fprintf('-----------------------------------------------------------\n');
        end
        if obj.simparams.verbosity > 1
          % fprintf('Iter     ||rhs||\n');
          fprintf('Iter     ||rhs||     ||rhs||/||rhs_0||\n');
        end

        % Compute Rhs and matrices of NonLinear models
        % Loop over available linear models and compute the jacobian
        assembleSystem(obj.domain,obj.dt);

        % Apply BCs to the blocks of the linear system
        applyBC(obj.domain,obj.t);

        % compute Rhs norm
        rhs = getRhs(obj.domain);
        rhsNorm = norm(rhs,2);
        rhsNormIt0 = rhsNorm;
        residual(1,1) = rhsNormIt0;
        residual(1,2) = 1.;

        % consider output of local field rhs contribution
        tolWeigh = obj.simparams.relTol*rhsNorm;
        obj.iter = 0;
        %
        if obj.simparams.verbosity > 1
          % fprintf('0     %e\n',rhsNorm);
          fprintf('0     %e     %e\n',rhsNorm,rhsNorm/rhsNormIt0);
        end

        while ((rhsNorm > tolWeigh) && (obj.iter < obj.simparams.itMaxNR) ...
            && (rhsNorm > absTol)) || obj.iter == 0
          obj.iter = obj.iter + 1;

          % Get system Jacobian
          J = getJacobian(obj.domain);

          % Solve linear system
          du = FCSolver.solve(J,rhs);

          % Update current model state
          updateState(obj.domain,du);

          % Compute Rhs and Matrices of NonLinear models
          assembleSystem(obj.domain,obj.dt);

          % Apply BCs to the blocks of the linear system
          applyBC(obj.domain,obj.t);

          rhs = getRhs(obj.domain);
          % compute Rhs norm
          rhsNorm = norm(rhs,2);
          residual(obj.iter+1,1) = rhsNorm;
          residual(obj.iter+1,2) = rhsNorm/rhsNormIt0;

          if obj.simparams.verbosity > 1
            fprintf('%d     %e     %e\n',obj.iter,residual(obj.iter+1,1),residual(obj.iter+1,2));
          end
        end
        %

        % Check for convergence
        flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);
        % flConv = (rhsNorm < tolWeigh);

        if flConv % Convergence

          state.t = obj.t;

          % print results at current time
          printState(obj.domain);

          % advance the state after convergence
          advanceState(obj.domain);

          % save solution statistics
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

      if ~flConv

        % backstep
        obj.domain.state = copy(obj.domain.stateOld);
        obj.t = obj.t - obj.dt;
        obj.tStep = obj.tStep - 1;

        obj.dt = obj.dt/obj.simparams.divFac; 

        if obj.dt < obj.simparams.dtMin
            error('Minimum time step reached')
        elseif obj.simparams.verbosity > 0
          fprintf('\n %s \n','BACKSTEP');
        end

      else

        % go to next time step
        tmpVec = obj.simparams.multFac;
        obj.dt = min([obj.dt * min(tmpVec), obj.simparams.dtMax]);
        obj.dt = max([obj.dt obj.simparams.dtMin]);
        obj.domain.stateOld = copy(obj.domain.state);
        
        % limit time step to end of simulation time
        if ((obj.t + obj.dt) > obj.simparams.tMax)
          obj.dt = obj.simparams.tMax - obj.t;
        end
      end
    end

  end

  methods (Static)
    function sol = solve(J,rhs)
      sol = J\(-rhs);
    end

  end
end