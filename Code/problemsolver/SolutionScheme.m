classdef (Abstract) SolutionScheme < handle
  % General solution scheme class
  % Implement the time loop and basic operations that are common to all
  % solution schemes
  

  properties (SetAccess = protected, GetAccess = public)

    tOld                % tOld: previous converged time instant
    t                   % simulation time
    tStep               % simulation time step
    dt                  % current time step size
    dtSave              % stored time step size 
    attemptedReset      % flag for attempting a configuration reset
    iniState            % initial state of the simulation for solver reset
    isFirstRun = true   % flag if the simulation is first ever or first after a reset   
    targetTimeID        % id for printTimes list for which solution is required
    restartControl      % variable to control the save point for restart
    rhsNormIt0 = 0      % initial residual norm of current step
    
  end


  properties (Access = public)

    model
    simparams
    output
    linsolver             

  end


  methods(Abstract,Access=public)
    % every solution scheme must implement the logic to solve a time step
    % and give feedback on convergence
    conv = solveStep(obj)
  end


  methods (Access = public)
    function obj = SolutionScheme(varargin)
      % assert(nargin > 1 && nargin < 9,"Wrong number of input arguments " + ...
      %   "for general solver")
      obj.setSolutionScheme(varargin{:});
    end

    function simulationLoop(obj,varargin)

      % Check restart option

      if nargin < 2 || ~strcmpi(varargin{1},"restart")

        % The simulation is assumed to start from the beginning
        obj.tStep = 0;
        obj.t = obj.simparams.tIni;
        obj.tOld = obj.t;
        obj.dt = obj.simparams.dtIni;
        obj.dtSave = obj.dt;

        initialize(obj);
        setLinearSolver(obj,varargin{:});
      end

      obj.t = obj.t + obj.dt;

      while obj.t <= obj.simparams.tMax
        initializeTimeStep(obj)

        gresLog().log(-1,'\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
        gresLog().log(-1,'-----------------------------------------------------------\n');

        % solve current time step
        % flConv: flag for convergence
        % dtOut: requested time step from the physics solver

        conv = solveStep(obj);

        % move to the next time step
        manageNextTimeStep(obj,conv)
      end

      obj.finalize;

      gresLog().log(-1,"Simulation completed successfully \n")
    end

    function saveHistory(obj)
      obj.output.saveHistory();
    end

    function saveRestart(restartObj,varargin)

      if restartObj.output.restart.overwrite
        restartFile = restartObj.output.restart.filename;
      else
        % If the file already exists, it will be overwritten
        restartFile = sprintf('%s_%05d',restartObj.output.restart.filename,...
          round(restartObj.t));
      end

      if any(strcmpi(varargin,'fail'))
        restartFile = sprintf('%s_FailAt%05d',restartObj.output.restart.filename,...
          round(restartObj.t));
      end

      save(restartFile,"restartObj","-v7.3");

      gresLog().log(2,'\n Restart file saved: %s\nSimulation time: %g\n',...
        restartFile,restartObj.t);
    end
  end


  methods (Access = protected)

    function setSolutionScheme(obj,varargin)
      default = struct('simulationparameters',SimulationParameters.empty,...
                       'output',missing,...
                       'model',Model.empty);

      params = readInput(default,varargin{:});
      obj.model = params.model;
      obj.simparams = params.simulationparameters;

      if ~ismissing(params.output)
        obj.output = params.output;
      else
        obj.output = OutState();
      end

      assert(~isempty(obj.simparams),"Input 'simulationParameters'" + ...
        " is required for SolutionScheme")
      
    end


    function initialize(obj)

      % restore the solution scheme object at its initial state
      if ~obj.isFirstRun
        obj.reset();
      end

      % store initial state and setup simulation

      for i = 1:numel(obj.model.domains)
        dom = obj.model.domains(i);
        state = dom.getState;
        state.time = obj.simparams.tIni;
        dom.setState(state);
        initialize(dom);
        state = dom.getState();
        obj.iniState.domains{i} = state;
        % set old and initial state 
        dom.setStateInit(state);
        dom.setStateOld(state);
        dom.outstate = obj.output;
        dom.simparams = obj.simparams;
        dom.domainId = i;
      end

      for i = 1:numel(obj.model.interfaces)
        interf = obj.model.interfaces(i);
        state = interf.getState;
        state.time = obj.simparams.tIni;
        interf.setState(state);
        initialize(interf);
        state = interf.getState();
        obj.iniState.interfaces{i} = state;
        % set old and initial state
        interf.setStateInit(state);
        interf.setStateOld(state);
        interf.id = i;
        interf.outstate = obj.output;
      end

      obj.isFirstRun = false;

      obj.attemptedReset = ~obj.simparams.attemptSimplestConfiguration || isempty(obj.model.interfaces);

      % Prepare folder to save the restart files.
      obj.output = obj.output.prepareRestart();

      if ismissing(obj.output.restart.dt)
        obj.restartControl = 1;
      else
        obj.restartControl = 0.;
      end
      
    end

    function reset(obj)
      % reset the simulation  at its initial state

      for i = 1:numel(obj.model.domains)
        setState(obj.model.domains(i),obj.iniState.domains{i});
      end

      for i = 1:numel(obj.model.interfaces)
        obj.model.interfaces(i).setState(obj.iniState.interfaces{i});
      end

      if ~isempty(obj.output)
        obj.output.reset();
      end

      obj.isFirstRun = false;
    end

    function finalize(obj)
      if ~isempty(obj.output)
        obj.output.savePvd();

        obj.output.saveHistory();
      end
    end

    function manageNextTimeStep(obj,flConv)
      if ~flConv && ~obj.attemptedReset
        % allow a configuration reset to attempt saving the simulation

        obj.model.resetConfiguration();

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
        dtNew = obj.dt/obj.simparams.divFac;  % Time increment chop
        goBackState(obj);
        
        if dtNew < obj.simparams.dtMin
          if obj.output.restart.on
            obj.saveRestart('fail');
          end
          error('Minimum time step reached')
        else
          obj.dt = dtNew;
          obj.dtSave = dtNew;
          obj.t = obj.t + dtNew;
          gresLog().log(0,'\n %s \n','BACKSTEP')
        end

        return

      else
        % TIME STEP CONVERGED - advance to the next time step
        printState(obj);
        advanceState(obj);

        if obj.output.restart.on
          saveState(obj);
        end

        % go to next time step
        obj.dt = getNextDt(obj);
        obj.tOld = obj.t;
        obj.t = obj.t + obj.dt;

        % allow new survival attempts on new time steps
        if obj.simparams.attemptSimplestConfiguration
          obj.attemptedReset = false;
        end

      end
    end

    function sol = solve(obj,J,res,varargin)
      
      % Extract the Newton iteration count when provided
      if ~isempty(varargin)
         newtonIter = varargin{1};
      else
         newtonIter = 1;
      end

      % Actual solution of the system
      [sol,~] = obj.linsolver.SolveLin(J,res,obj.t,newtonIter);
    end

    function setLinearSolver(obj,varargin)

      input = readInput(struct('variableName',missing),varargin{:});

      if ismissing(input.variableName)
        input.variableName = [];
      end

      obj.linsolver = linearSolver(obj,input.variableName);

    end

    function applyDirVal(obj)
      for i = 1:numel(obj.model.domains)
        discretizer = obj.model.domains(i);

        % Check if boundary conditions are defined for the i-th domain
        if ~isempty(obj.model.domains(i).bcs)

          % Apply Dirichlet boundary values to i-th domain
          applyDirVal(discretizer,obj.t);
        end
      end
    end

    function applyBC(obj)
      for i = 1:numel(obj.model.domains)
        domain = obj.model.domains(i);
        % Apply BCs to the blocks of the linear system
        applyBC(domain, obj.t);
      end
    end

    function updateState(obj,dSol)

      obj.model.updateState(dSol);

    end


    function initializeTimeStep(obj)

      obj.tStep = obj.tStep + 1;

      % set current time into state objects
      for i = 1:numel(obj.model.domains)
        dom = obj.model.domains(i);
        dom.setState(obj.t,'time');
        timeStepSetup(obj.model.domains(i));
      end

      for i = 1:numel(obj.model.interfaces)
        interf = obj.model.interfaces(i);
        interf.setState(obj.t,'time');
        timeStepSetup(obj.model.interfaces(i));
      end
    end

    function printState(obj)

      if obj.output.timeID <= length(obj.output.timeList)
        outTime = obj.output.timeList(obj.output.timeID);

        % loop over print times contained in the current time step
        while outTime <= obj.t
          assert(outTime >= obj.tOld, 'Print time %f out of range (%f - %f)',...
            outTime, obj.tOld, obj.t);
          assert(obj.t - obj.tOld > eps('double'),...
            'Time step is too small for printing purposes');

          % compute factor to interpolate current and old state variables
          fac = (outTime - obj.tOld)/(obj.t - obj.tOld);

          if isnan(fac) || isinf(fac)
            fac = 1;
          end

          % write results to vtk file
          obj.printVTK(fac,outTime,obj.output.timeID);

          % write results to MAT-file
          obj.printMAT(fac,obj.output.timeID);

          % move to next print time
          obj.output.timeID = obj.output.timeID + 1;

          if obj.output.timeID > length(obj.output.timeList)
            break
          else
            outTime = obj.output.timeList(obj.output.timeID);
          end
        end
      end
    end


    function saveState(obj)

      % save restart point

      if ~ismissing(obj.output.restart.time)

        % save object at specified times
        if obj.restartControl > length(obj.output.restart.time)
          obj.output.restart.on = false;
        else
          if obj.t >= obj.output.restart.time(obj.restartControl)
            obj.restartControl = obj.restartControl+1;
            obj.saveRestart();
          end
        end        

      elseif ~ismissing(obj.output.restart.dt)

        % restart after a specified time interval
        obj.restartControl = obj.restartControl+obj.dt;

        if obj.restartControl >= obj.output.restart.dt
          obj.restartControl = 0;
          obj.saveRestart();
        end

      end


    end

    function printVTK(obj,fac,outTime,tID)
      if obj.output.writeVtk
        % set folders        
        obj.output.prepareOutputFolders(tID);

        obj.output.vtkFile = com.mathworks.xml.XMLUtils.createDocument('VTKFile');
        toc = obj.output.vtkFile.getDocumentElement;
        toc.setAttribute('type', 'vtkMultiBlockDataSet');
        toc.setAttribute('version', '1.0');
        blocks = obj.output.vtkFile.createElement('vtkMultiBlockDataSet');

        % append blocks looping into domains and interfaces
        for i = 1:numel(obj.model.domains)
          vtmBlock = obj.model.domains(i).writeVTK(fac,outTime);
          if ~isempty(vtmBlock)
            blocks.appendChild(vtmBlock);
          end
        end
        %
        for i = 1:numel(obj.model.interfaces)
          vtmBlock = obj.model.interfaces(i).writeVTKfile(fac,outTime);
          if ~isempty(vtmBlock)
            blocks.appendChild(vtmBlock);
          end
        end

        toc.appendChild(blocks);
        obj.output.writeVTMFile();
      end
    end

    function printMAT(obj,fac,timeID)
      % write results into an output structure
      
      if obj.output.writeSolution
        for i = 1:numel(obj.model.domains)
          obj.model.domains(i).writeSolution(fac,timeID);
        end

        for i = 1:numel(obj.model.interfaces)
          obj.model.interfaces(i).writeSolution(fac,timeID);
        end
      end
    end

    function advanceState(obj)

      obj.model.advanceState();

    end


    function goBackState(obj)

      obj.model.goBackState();

    end


    function dt = getNextDt(obj)

      targetTimes = obj.output.getTargetTimes;
      obj.targetTimeID = 1 + sum(targetTimes <= obj.t + obj.simparams.dtMin);
      if obj.targetTimeID > length(targetTimes)
        targetTime = Inf;
      else
        targetTime = targetTimes(obj.targetTimeID);
      end


      [dt,m] = min([obj.dtSave * obj.simparams.multFac, ...
                    obj.simparams.dtMax,...
                    targetTime - obj.t,...
                    obj.simparams.tMax - obj.t,...
                    ]);

      dt = max([dt obj.simparams.dtMin]);

      if m == 3
        % advance target time ID 
        obj.targetTimeID = obj.targetTimeID + 1;
      else
        % save the current step
        obj.dtSave = dt;
      end

      % limit time step to end of simulation time
      if m == 4 && dt == obj.simparams.dtMin
        % force exit from time simulation loop
        obj.t = obj.simparams.tMax;
        % print here the end of simulation if requested
        printState(obj);
      end
    end

    
    

  end

  methods (Static)

    function obj = loadRestart(restartFile)
      
      obj = load(restartFile).restartObj;

    end

    

  end
end

