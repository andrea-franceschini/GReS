classdef ModelManager < handle
  % General model class
  % This class stores all the information related to the problem at hand
  % and store a map to existing domain solvers and



  properties (GetAccess=private, SetAccess=private)
    domains             % db for physics solver in the model
    interfaces          % db for interface solvers in the model
    solutionScheme
    simparams           % simulation parameters
    t = 0               % simulation time
    tStep = 0           % simulation time step
    dt
    nDom
    nInterf
  end


  methods (Access = public)
    function obj = ModelManager(varargin)

      obj.domains = containers.Map('KeyType','double','ValueType','any');
      obj.interfaces = containers.Map('KeyType','double','ValueType','any');

      if ~isempty(varargin)
        obj.simparams = varargin{1};
      end
    end


    function runProblem(obj)

      obj.nDom = numel(keys(obj.domains));
      obj.nInterf = numel(keys(obj.domains));


      % this method implement the simulation time loop
      % solutionScheme is an handle to a function implementing a specific type of solution
      % strategy for each time step

      while obj.t < obj.simparams.tMax

        % Update the simulation time and time step ID
        obj.tStep = obj.tStep + 1;

        obj.t = obj.t + obj.dt;

        % solve time step with chosen solution scheme
        isConverged = obj.solutionScheme.solveTimeStep(obj.t,obj.dt);

        % manage the next time step depending on convergence
        delta_t = manageNextTimeStep(obj,obj.dt,isConverged);

      end
      %
    end

    % function printSolution()
    %   % replaces finalize state
    % end

    function createModel(obj,fileName)
      % utility to create an entire GReS model from a single xml file

      str = readstruct(fileName,AttributeSuffix="");

      obj.simparams = SimulationParameters(fileName);

      for i = 1:numel(str.Domain)
        addDomain(obj,str.Domain(i));
      end

      if isfield(str,"SolutionStrategy")
        strategy = fieldnames(str.SolutionStrategy);
        obj.solutionScheme = feval(strategy{1},obj.simparams,str.SolutionStrategy);
      else
        % use general fully coupled solution strategy as the default one
        obj.solutionScheme = FullyCoupled(obj,obj.simparams);
      end

      if isfield(str,"Interface")
        for iI = 1:numel(str.Interface)
          addInterface(obj,str.Interface(iD));
        end
      else
        if numel(str.Domain) > 1
          warning("Running a problem with multiple domain " + ...
            "without any connecting interfaces");
        end
      end
    end


    function addDomain(obj,varargin)

      % add a domain
      n = numel(keys(obj.domains));
      obj.domains(n+1) = Domain(varargin{:});
    end

    function addInterface(obj,varargin)

      % TO DO
    end

    function out = getDomain(obj,id)
      out = obj.domains(id);
    end

    function out = getInterface(obj,id)
      out = obj.interfaces(id);
    end

  end


  methods (Access = private)
    function dt = manageNextTimeStep(obj, dt, flConv)

      % Time step did not converge - perform a backstep

      if ~flConv
        goBackState(obj);

        % Roll back time
        obj.t     = obj.t     - obj.dt;
        obj.tStep = obj.tStep - 1;

        % Reduce time step
        dt    = dt    / obj.simparams.divFac;
        obj.dt = obj.dt / obj.simparams.divFac;

        % Check for minimum time step
        if min(dt, obj.dt) < obj.simparams.dtMin
          if obj.simparams.goOnBackstep == 1
            flConv = 1;     % force convergence and move on
          elseif obj.simparams.goOnBackstep == 0
            error('Minimum time step reached');
          end

        % Optional verbosity
        elseif obj.simparams.verbosity > 0
          fprintf('\n %s \n', 'BACKSTEP');
        end
      end

      % Time step converged - increase dt and advance

      if flConv
        tmpVec = obj.simparams.multFac;

        % Increase dt but respect bounds
        obj.dt = min([obj.dt * min(tmpVec), obj.simparams.dtMax]);
        obj.dt = max(obj.dt, obj.simparams.dtMin);

        goOnState(obj);

        % Do not exceed final time
        if obj.t + obj.dt > obj.simparams.tMax
          obj.dt = obj.simparams.tMax - obj.t;
        end
      end
    end


    function goOnState(obj)
      % transfer current state into previous state
      for i = 1:obj.nDom
        obj.domains(i).goOnState();
      end

      for i = 1:obj.nInterf
        obj.interfaces(i).goOnState();
      end
    end

    function goBackState(obj)
      % transfer previous state into current state (backstep)
      for i = 1:nueml
        obj.domains(i).goBackState();
      end

      for i = 1:obj.nInterf
        obj.interfaces(i).goBackState();
      end
    end

  end

end
