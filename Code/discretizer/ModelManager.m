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

      % this method implement the simulation time loop
      % solutionScheme is an handle to a function implementing a specific type of solution
      % strategy for each time step

      while obj.t < obj.simparams.tMax

        % Update the simulation time and time step ID
        obj.tStep = obj.tStep + 1;

        [obj.t, dt] = obj.updateTime(flConv, delta_t);

        % solve time step with chosen solution scheme
        isConverged = obj.solutionScheme.solveTimeStep(obj.t,dt);

        % manage the next time step depending on convergence
        delta_t = manageNextTimeStep(obj,dt,isConverged);
        
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

end
