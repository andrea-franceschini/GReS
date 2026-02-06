classdef SolutionScheme < handle
  % General solution scheme class

  properties (Access = protected)
    %
    nDom                % number of domains in the model
    nInterf             % number of interfaces in the model
    %
    tOld                % tOld: previous converged time instant
    t = 0               % simulation time
    tStep = 0           % simulation time step
    dt                  % current time step size
    nVars               % total number of inner variable fields in the model
    attemptedReset      % flag for attempting a configuration reset
  end


  properties (Access = public)
    linsolver             % instance of linear solver object
    output                % object handling the output of the simulation
    simparams             % parameters of the simulations (shared)
    domains               % array of Discretizer objects
    interfaces            % cell array of interfaces objects
  end


  methods (Access = public)
    function obj = SolutionScheme(varargin)

      assert(nargin > 1 && nargin < 4,"Wrong number of input arguments " + ...
        "for general solver")

      obj.setSolutionScheme(varargin{:});

    end

    function simulationLoop(obj)

      % Initialize the time step increment
      obj.dt = obj.simparams.dtIni;

      obj.linsolver = setLinearSolver(obj);

      while obj.t < obj.simparams.tMax

        initializeTimeStep(obj)

        gresLog().log(-1,'\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
        gresLog().log(-1,'-----------------------------------------------------------\n');

        % solve current time step
        % flConv: flag for convergence
        % dtOut: requested time step from the physics solver
        [flConv,dtOut] = solveStep(obj);

        % move to the next time step
        manageNextTimeStep(obj,flConv,dtOut)

      end

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
          case {'interfaces','interface'}
            obj.interfaces = value;
          otherwise
            error('Unknown input key %s for SolutionScheme \n', key);
        end
      end

      obj.nDom = numel(obj.domains);
      obj.nInterf = numel(obj.interfaces);

      assert(~isempty(obj.simparams),"Input 'simulationParameters'" + ...
        " is required for SolutionScheme")
      assert(obj.nDom > 0,"Input 'domains'" + ...
        " is required for SolutionScheme")

      obj.nVars = 0;

      for i = 1:obj.nDom
        obj.domains(i).domainId = i;
        obj.domains(i).simparams = obj.simparams;
        obj.domains(i).stateOld = copy(obj.domains(iD).getState());
        obj.nVars = obj.nVars + obj.domains(iD).dofm.getNumberOfVariables();
      end

      obj.attemptedReset = ~obj.simparams.attemptSimplestConfiguration || obj.nInterf == 0;

    end
    

    function manageNextTimeStep(obj,flConv,dtPhysics)

      if ~flConv && ~obj.attemptedReset

        % allow a configuration reset to attempt saving the simulation

        for i = 1:obj.nDom
          resetConfiguration(obj.domains(i));
        end

        for i = 1:obj.nInterf
          resetConfiguration(obj.interfaces{i});
        end

        obj.tStep = obj.tStep - 1;
        obj.t = obj.t - obj.dt;

        obj.attemptedReset = true;

        gresLog().log(1,"Reset to simplest configuration \n")

        return

      end


      if ~flConv
        % BACKSTEP
        % newton did not converge or configuration changed too many times

        obj.t = obj.tOld;
        obj.tStep = obj.tStep - 1;
        obj.dt = obj.dt/obj.simparams.divFac;  % Time increment chop

        goBackState(obj)

        if obj.dt < obj.simparams.dtMin
          error('Minimum time step reached')
        else
          gresLog().log(0,'\n %s \n','BACKSTEP')
        end

        return

      else

        % TIME STEP CONVERGED - advance to the next time step
        printState(obj);
        advanceState(obj);

        % go to next time step
        tmpVec = obj.simparams.multFac;
        obj.dt = min([obj.dt * min(tmpVec), dtPhysics, obj.simparams.dtMax]);
        obj.dt = max([obj.dt obj.simparams.dtMin]);

        % limit time step to end of simulation time
        if ((obj.t + obj.dt) > obj.simparams.tMax)
          obj.dt = obj.simparams.tMax - obj.t;
        end

        % allow new survival attempts on new time steps
        if obj.simparams.attemptSimplestConfiguration
          obj.attemptedReset = false;
        end

      end

    end

  
    function setLinearSolver(obj)
      % Check if there is manual input from the user, if not use defaults
      start_dir = pwd;
      chronos_xml = fullfile(start_dir,'linsolver.xml');
      if(isfile(chronos_xml))
        obj.linsolver = linearSolver(obj.domains,obj.interfaces,chronos_xml);
      else
        if gresLog().getVerbosity > 2
          fprintf('Using default values for linsolver\n');
        end
        obj.linsolver = linearSolver(obj.domains,obj.interfaces);
      end
    end


    function sol = solve(obj,J,rhs)

      rhs = cell2matrix(rhs);

      % Actual solution of the system
      [sol,~] = obj.linsolver.Solve(J,-rhs,obj.t);
    end



    function applyDirVal(obj)
      for i = 1:obj.nDom
        discretizer = obj.domains(i);

        % Check if boundary conditions are defined for the i-th domain
        if ~isempty(obj.domains(i).bcs)

          % Apply Dirichlet boundary values to i-th domain
          applyDirVal(discretizer,obj.t);
        end
      end
    end


    function applyBC(obj)
      for i = 1:obj.nDom
        discretizer = obj.domains(i);
        % Apply BCs to the blocks of the linear system
        applyBC(discretizer, obj.t);

        % Apply BC to domain coupling matrices
        for j = discretizer.interfaceList
          applyBC(obj.interfaces{j},i,discretizer.bcs,obj.t);
        end
      end
    end


    function updateState(obj,dSol)
      % update domain and interface state using incremental solution
      dSol_fix = dSol;
      for i = 1:obj.nDom
        N = obj.domains(i).dofm.totDoF;
        du = dSol(1:N);
        updateState(obj.domains(i),du);
        dSol = dSol(N+1:end);
      end

      % update interface state
      for j = 1:obj.nInterf
        N = obj.interfaces{j}.totMult;
        if N == 0
          du = dSol_fix;
        else
          du = dSol(1:N);
        end
        obj.interfaces{j}.updateState(du);
        dSol = dSol(N+1:end);
      end
    end

    function initializeTimeStep(obj)

      obj.tStep = obj.tStep + 1;
      obj.tOld = obj.t;
      obj.t = obj.t + obj.dt;

      % set current time into state objects
      for i = 1:obj.nDom
        dom = obj.domains(i);
        dom.state.t = obj.t;
      end

      for i = 1:obj.nInterf
        interf = obj.interfaces{i};
        interf.state.t = obj.t;
      end

    end

    function printState(obj)

      if isempty(obj.output)
        return
      end

      if obj.output.timeID <= length(obj.output.timeList)

        outTime = obj.output.timeList(obj.output.timeID);

        % loop over print times contained in the current time step
        
        while outTime <= obj.t

          assert(outTime >= obj.tOld, 'Print time %f out of range (%f - %f)',...
            outTime, obj.tOld, obj.t);

          assert(obj.t - obj.tOld > eps('double'),...
            'Time step is too small for printing purposes');

          % compute factor to interpolate current and old state variables
          fac = (outTime - obj.t)/(obj.t - obj.tOld);

          if isnan(fac) || isinf(fac)
            fac = 1;
          end

          % print vtk and matlab files
          for i = 1:obj.nDom
            obj.domains(i).writeVTK(obj,fac,outTime);
            obj.domains(i).writeMatFile(fac,obj.output.timeID);
          end

          for i = 1:obj.nInterf
            obj.interfaces{i}.writeVTK(obj,fac,outTime);
            obj.interfaces{i}.writeMatFile(fac,obj.output.timeID);
          end

          obj.output.timeID = obj.output.timeID + 1;

          if obj.output.timeID > length(obj.output.timeList)
            break
          else
            outTime = obj.output.timeList(obj.output.timeID);
          end

        end

      end

    end


    function advanceState(obj)

      for i = 1:obj.nDom
        dom = obj.domains(i);
        advanceState(dom);
      end

      for i = 1:obj.nInterf
        interf = obj.interfaces{i};
        advanceState(interf);
      end

    end

    function goBackState(obj)

      for i = 1:obj.nDom
        dom = obj.domains(i);
        goBackState(dom);
      end

      for i = 1:obj.nInterf
        interf = obj.interfaces{i};
        goBackState(interf);
      end


    end

  end
end

