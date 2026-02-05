classdef SolutionScheme < handle
  % General solution scheme class

  properties (Access = protected)
    %
    nDom                % number of domains in the model
    nInterf             % number of interfaces in the model
    %
    t = 0               % simulation time
    tStep = 0           % simulation time step
    dt                  % current time step size
  end


  properties (Access = public)
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



    function NonLinearLoop(obj)

      % Initialize the time step increment
      obj.dt = obj.simparams.dtIni;

      %

      for i = 1:obj.nDom
        obj.domains(i).applyDirVal(obj.t);
      end



      %%% TIME LOOP %%
      while obj.t < obj.simparams.tMax

        % Update the simulation time and time step ID
        absTol = obj.simparams.absTol;

        initializeTimeStep(obj)

        gresLog().log(-1,'\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
        gresLog().log(-1,'-----------------------------------------------------------\n');

        % reset active set iteration counter
        obj.iterConfig = 0;

        hasConfigurationChanged = true;

        %%% CONFIGURATION LOOP %%
        while (hasConfigurationChanged) && (obj.iterConfig < obj.simparams.itMaxConfig)

          gresLog().log(0,'\nConfiguration iteration n. %i \n', obj.iterConfig);
          obj.iterConfig = obj.iterConfig + 1;

          for i = 1:obj.nDom
            obj.domains(i).applyDirVal(obj.t);
          end

          gresLog().log(1,'Iter     ||rhs||     ||rhs||/||rhs_0||\n');

          for i = 1:obj.nDom
            obj.domains(i).assembleSystem(obj.dt);
          end

          for i = 1:obj.nInterf
            obj.interfaces{i}.assembleConstraint();
          end

          for i = 1:obj.nDom
            obj.domains(i).applyBC(obj.t);
          end

          rhs = assembleRhs(obj);
          rhsNorm = norm(cell2mat(rhs),2);
          rhsNormIt0 = rhsNorm;

          tolWeigh = obj.simparams.relTol*rhsNorm;

          gresLog().log(1,'0     %e     %e\n',rhsNorm,rhsNorm/rhsNormIt0);

          flConv = false;

          % reset non linear iteration counter
          obj.iterNL = 0;

          %%% NONLINEAR LOOP %%%
          while (~flConv) && (obj.iterNL < obj.simparams.itMaxNR)

            obj.iterNL = obj.iterNL + 1;

            J = assembleJacobian(obj);

            % solve linear system
            du = solve(obj,J,rhs);

            c = 0;

            % update simulation state with linear system solution
            for i = 1:obj.nDom
              if obj.nDom == 1 && obj.nInterf == 0
                sol = du;
              else
                nDof = obj.domains(i).getNumbDoF();
                sol = du(c+1:c+nDof);
                c = c + nDof;
              end
              obj.domains(i).updateState(sol);
            end

            for i = 1:obj.nInterf
              nDof = obj.interfaces{i}.getNumbDoF();
              sol = du(c+1:c+nDof);
              obj.interfaces{i}.updateState(sol);
              c = c + nDof;
            end

            % reassemble system
            for i = 1:obj.nDom
              obj.domains(i).assembleSystem(obj.dt);
            end

            for i = 1:obj.nInterf
              obj.interfaces{i}.assembleConstraint();
            end

            for i = 1:obj.nDom
              obj.domains(i).applyBC(obj.t);
            end

            rhs = assembleRhs(obj);
            rhsNorm = norm(cell2mat(rhs),2);
            gresLog().log(1,'%d     %e     %e\n',obj.iterNL,rhsNorm,rhsNorm/rhsNormIt0);

            % Check for convergence
            flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);

          end % end newton loop

          if flConv % Newton Convergence

            hasConfigurationChanged = false;

            % update the active set
            for i = 1:obj.nDom
              hasConfigurationChanged = any([hasConfigurationChanged; ...
                obj.domains(i).updateConfiguration()]);
            end

            for i = 1:obj.nInterf
              hasConfigurationChanged = any([hasConfigurationChanged; ...
                obj.interfaces{i}.updateConfiguration()]);
            end

          else

            break

          end

        end


          manageNextTimeStep(obj,flConv,hasConfigurationChanged);

      end % time marching
      %
      finalizeOutput(obj.output);

      gresLog().log(-1,"Simulation completed successfully \n")
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

      for i = 1:obj.nDom
        obj.domains(i).domainId = i;
        obj.domains(i).simparams = obj.simparams;

      end
      
    end

    function setNonLinearSolver(obj,simparams,dom,interf)

      % assumption: same set of simulation parameters for each domain
      obj.simparams = simparams;
      obj.domains = dom;
      obj.nDom = numel(dom);
      obj.interfaces = interf;
      obj.nInterf = numel(interf);
      obj.attemptedReset = ~obj.simparams.attemptSimplestConfiguration;

      obj.nVars = 0;
      for iD = 1:obj.nDom
        obj.domains(iD).stateOld = copy(obj.domains(iD).getState());
        obj.domains(iD).simparams = simparams;
        obj.nVars = obj.nVars + obj.domains(iD).dofm.getNumberOfVariables();
      end

      for iI = 1:obj.nInterf

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



    function out = computeRhsNorm(obj)

      %Return maximum norm of the entire domain
      rhsNorm = zeros(obj.nDom,1);
      for i = 1:obj.nDom
        nRhs = length(obj.domains(i).dofm.subList);
        rhsNorm_loc = zeros(nRhs,1);
        for j = 1:nRhs
          rhsNorm_loc(j) = norm(obj.domains(i).rhs{j}, obj.simparams.pNorm);
        end
        rhsNorm(i) = sqrt(sum(rhsNorm_loc.^2));
      end
      out = norm(rhsNorm);
    end



    function J = assembleJacobian(obj)

      J = cell(obj.nVars + obj.nInterf);

      k = 0;

      for iD = 1:obj.nDom

        dom = obj.domains(iD);
        nV = dom.dofm.getNumberOfVariables;

        % inner domain blocks
        J(k+1:k+nV,k+1:k+nV) = getJacobian(dom);

        for iI = 1:numel(dom.interfaceList)

          q = dom.interfaceList(iI);

          % domain coupling blocks
          [J(k+1:k+nV,obj.nVars+q), J(obj.nVars+q,k+1:k+nV)] = ...
            getInterfaceJacobian(dom,iI);

        end

        k = k + nV;

      end

      for iI = 1:obj.nInterf
        interf = obj.interfaces{iI};
        % constraint blocks
        J{obj.nVars+iI,obj.nVars+iI} = getJacobian(interf);
      end

    end


    function rhs = assembleRhs(obj)
      % assemble blocks of rhs for multidomain system

      % each variable field of each domain represents a cell row
      rhs = cell(obj.nVars + obj.nInterf,1);
      
      k = 0;

      for iD = 1:obj.nDom
        nV = obj.domains(iD).dofm.getNumberOfVariables;
        rhs(k+1:k+nV) = getRhs(obj.domains(iD));
        k = k+nV;
      end

      for iI = 1:obj.nInterf

        rhs{k+1} = getRhs(obj.interfaces{iI});

        % each interface has one only multiplier field!
        k = k+1;
      end
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
      obj.t = obj.t + obj.dt;

      for i = 1:obj.nDom
        dom = obj.domains(i);
        dom.state.t = obj.t;
      end

      for i = 1:obj.nInterf
        interf = obj.interfaces{i};
        interf.state.t = obj.t;
      end

    end

    function setOutput(obj)

      for iD = 1:obj.nDom
        obj.domains(iD).stateOld = copy(obj.domains(iD).getState());
        obj.domains(iD).simparams = simparams;
        obj.nVars = obj.nVars + obj.domains(iD).dofm.getNumberOfVariables();
      end

      for iI = 1:obj.nInterf

      end
    end



    function manageNextTimeStep(obj,NLConv,configurationChanged)

      if ~NLConv && ~obj.attemptedReset

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


      if ~NLConv || configurationChanged

        % BACKSTEP
        % newton did not converge or configuration changed too many times

        obj.t = obj.t - obj.dt;
        obj.tStep = obj.tStep - 1;
        obj.dt = obj.dt/obj.simparams.divFac;  % Time increment chop

        for i = 1:obj.nDom
          goBackState(obj.domains(i));
        end

        for i = 1:obj.nInterf
          goBackState(obj.interfaces{i},obj.dt);
        end

        if obj.dt < obj.simparams.dtMin
          error('Minimum time step reached')
        else
          gresLog().log(0,'\n %s \n','BACKSTEP')
        end

        return

      else 

        % TIME STEP CONVERGED - advance to the next time step

        for i = 1:obj.nDom
          dom = obj.domains(i);
          printState(dom);
          advanceState(dom);
        end

        for i = 1:obj.nInterf
          interf = obj.interfaces{i};
          printState(interf);
          advanceState(interf);
        end

        % go to next time step
        tmpVec = obj.simparams.multFac;
        obj.dt = min([obj.dt * min(tmpVec), obj.simparams.dtMax]);
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

  end
end

