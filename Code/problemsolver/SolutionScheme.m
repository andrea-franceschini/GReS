classdef (Abstract) SolutionScheme < handle
  % General solution scheme class
  % Implement the time loop and basic operations that are common to all
  % solution schemes
  

  properties (SetAccess = protected, GetAccess = public)
    %
    %
    tOld                % tOld: previous converged time instant
    t                   % simulation time
    tStep               % simulation time step
    dt                  % current time step size
    dtSave              % stored time step size 
    nVars               % total number of inner variable fields in the model
    attemptedReset      % flag for attempting a configuration reset
    iniState            % initial state of the simulation for solver reset
    isFirstRun = true   % flag if the simulation is first ever or first after a reset   
    targetTimeID        % id for printTimes list for which solution is required
    restartControl      % variable to control the save point for restart
  end


  properties (Access = public)
    nDom                % number of domains in the model
    nInterf             % number of interfaces in the model
    linsolver             % instance of linear solver object
    output                % object handling the output of the simulation
    simparams             % parameters of the simulations (shared)
    domains               % array of Discretizer objects
    interfaces            % cell array of interfaces objects
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
      flagRestart = any(strcmpi(varargin,"restart"));

      % Remove restart from varargin
      varargin(strcmpi(varargin,"restart")) = [];

      if ~flagRestart
        % Initialize time
        obj.tStep = 0;
        obj.t = obj.simparams.tIni;
        obj.tOld = obj.t;
        obj.dt = obj.simparams.dtIni;
        obj.dtSave = obj.dt;

        initialize(obj);
        setLinearSolver(obj,varargin{:});
        if ismissing(obj.output.restart.dt)
          obj.restartControl = 1;
        else
          obj.restartControl = 0.;
        end
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
      if obj.output.restart.on
        obj.saveRestart();
      end
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

      gresLog().log(2,'\nRestart file saved: %s\nSimulation time: %g\n',...
        restartFile,restartObj.t);
    end
  end


  methods (Access = protected)

    function setSolutionScheme(obj,varargin)
      default = struct('simulationparameters',SimulationParameters.empty,...
                       'output',missing,...
                       'domains',Discretizer.empty,...
                       'interface',missing);

      params = readInput(default,varargin{:});

      obj.simparams = params.simulationparameters;
      obj.domains = params.domains;

      if ~ismissing(params.interface)
        obj.interfaces = params.interface;
      end

      if ~ismissing(params.output)
        obj.output = params.output;
      else
        obj.output = OutState();
      end

      obj.nDom = numel(obj.domains);
      obj.nInterf = numel(obj.interfaces);

      assert(~isempty(obj.simparams),"Input 'simulationParameters'" + ...
        " is required for SolutionScheme")
      assert(obj.nDom > 0,"Input 'domains'" + ...
        " is required for SolutionScheme")
    end

    function initialize(obj)
      % restore the solution scheme object at its initial state
      if ~obj.isFirstRun
        obj.reset();
      end

      obj.nVars = 0;

      % store initial state and setup simulation

      for i = 1:obj.nDom
        dom = obj.domains(i);
        state = dom.getState;
        state.time = obj.simparams.tIni;
        dom.setState(state);
        obj.iniState.domains{i} = dom.getState();
        % set old and initial state 
        dom.setStateInit(state);
        dom.setStateOld(state);
        dom.outstate = obj.output;
        dom.simparams = obj.simparams;
        dom.domainId = i;
        obj.nVars = obj.nVars + dom.dofm.getNumberOfVariables();
        initialize(dom);
      end

      for i = 1:obj.nInterf
        interf = obj.interfaces{i};
        state = interf.getState;
        state.time = obj.simparams.tIni;
        obj.iniState.interfaces{i} = interf.getState();
        % set old and initial state
        interf.setState(state);
        interf.setStateInit(state)
        interf.setStateOld(state);
        interf.interfId = i;
        interf.outstate = obj.output;
        initialize(interf);
      end

      obj.isFirstRun = false;

      obj.attemptedReset = ~obj.simparams.attemptSimplestConfiguration || obj.nInterf == 0;

      % Prepare folder to save the restart files.
      obj.output = obj.output.prepareRestart();
    end

    function reset(obj)
      % reset the simulation  at its initial state

      for i = 1:obj.nDom
        setState(obj.domains(i),obj.iniState.domains{i});
      end

      for i = 1:obj.nInterf
        setState(get(obj.iniState.interfaces{i}));
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

        % go to next time step
        obj.dt = getNextDt(obj);
        obj.tOld = obj.t;
        obj.t = obj.t + obj.dt;

        % allow new survival attempts on new time steps
        if obj.simparams.attemptSimplestConfiguration
          obj.attemptedReset = false;
        end

        % save restart point
        if obj.output.restart.on
          if ismissing(obj.output.restart.dt)
            if (obj.t-obj.dt)>obj.output.restart.time(obj.restartControl)
              obj.saveRestart();
              obj.restartControl=obj.restartControl+1;
              if obj.restartControl>length(obj.output.restart.time)
                obj.output.restart.on = false;
              end
            end
          else
            obj.restartControl = obj.restartControl+obj.dt;
            if obj.restartControl>=obj.output.restart.dt
              obj.saveRestart();
              obj.restartControl = 0;
            end
          end
        end

      end
    end

    function sol = solve(obj,J,rhs,varargin)
      
      % Extract the Newton iteration count when provided
      if ~isempty(varargin)
         newtonIter = varargin{1};
      else
         newtonIter = 1;
      end

      % Actual solution of the system
      [sol,~] = obj.linsolver.SolveLin(J,rhs,obj.t,newtonIter);
    end

    function setLinearSolver(obj,varargin)
      if isempty(varargin)
         physname = [];
      else
         % check if the user provided the physics
         physname = varargin{1};
      end

      obj.linsolver = linearSolver(obj,physname);
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

      % set current time into state objects
      for i = 1:obj.nDom
        dom = obj.domains(i);
        dom.setState(obj.t,'time');
        timeStepSetup(obj.domains(i));
      end

      for i = 1:obj.nInterf
        interf = obj.interfaces{i};
        interf.setState(obj.t,'time');
        timeStepSetup(obj.interfaces{i});
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
        for i = 1:obj.nDom
          vtmBlock = obj.domains(i).writeVTK(fac,outTime);
          if ~isempty(vtmBlock)
            blocks.appendChild(vtmBlock);
          end
        end
        %
        for i = 1:obj.nInterf
          vtmBlock = obj.interfaces{i}.writeVTKfile(fac,outTime);
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
        for i = 1:obj.nDom
          obj.domains(i).writeSolution(fac,timeID);
        end

        for i = 1:obj.nInterf
          obj.interfaces{i}.writeSolution(fac,timeID);
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
      load(restartFile);
      obj=restartObj;
    end

    

  end
end

