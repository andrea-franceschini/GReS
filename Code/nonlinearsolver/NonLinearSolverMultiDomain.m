classdef NonLinearSolverMultiDomain < handle
  % Multi domain version of standard non linear solver class
  
  
  properties (Access = private)
      %
      nDom
      nInt
      models
      interface
      %
      t = 0
      tStep = 0
      iter
      dt
  end

  properties (Access = public)
      state
  end
  
  methods (Access = public)
      function obj = NonLinearSolverMultiDomain(models,interfaces)
      obj.setNonLinearSolver(models,interfaces);
    end

    function NonLinearLoop(obj)
      % Initialize the time step increment
      obj.dt = obj.simParameters.dtIni;  
      delta_t = obj.dt; % dynamic time step

      % Compute matrices of Linear models (once for the entire simulation)
      computeLinearMatrices(obj.linSyst,obj.stateTmp,obj.statek,obj.dt)
      %
      flConv = true; %convergence flag
      %
      % Loop over time
      while obj.t < obj.simParameters.tMax
        % Update the simulation time and time step ID
        if (obj.t > 4) && (obj.t < 5.5)
            absTol = 1.e-6;
        else 
            absTol = obj.simParameters.absTol;
        end
        obj.tStep = obj.tStep + 1;
        %new time update to fit the outTime list
        [obj.t, delta_t] = obj.updateTime(flConv, delta_t);
        %obj.t = obj.t + obj.dt;

        % Apply the Dirichlet condition value to the solution vector
        applyDirVal(obj.model, obj.bound, obj.t, obj.stateTmp);
        %
        if obj.simParameters.verbosity > 0
        fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,delta_t);
        fprintf('-----------------------------------------------------------\n');
        end
        if obj.simParameters.verbosity > 1
        fprintf('Iter     ||rhs||\n');
        end
        %
        % Compute Rhs and matrices of NonLinear models
        computeNLMatricesAndRhs(obj.linSyst,obj.stateTmp,obj.statek,obj.dt)
    
        % compute block Jacobian and block Rhs
        obj.linSyst.computeBlockJacobianAndRhs(delta_t);

        % Apply BCs to the block-wise system
        applyBCandForces(obj.model, obj.grid, obj.bound, obj.material, ...
          obj.t, obj.linSyst, obj.stateTmp);

        % compute Rhs norm
        [~, rhsNorm] = computeRhsNorm(obj,obj.linSyst);
        % consider output of local field rhs contribution

        tolWeigh = obj.simParameters.relTol*rhsNorm;
        obj.iter = 0;
        %
        if obj.simParameters.verbosity > 1
        fprintf('0     %e\n',rhsNorm);
        end
        while ((rhsNorm > tolWeigh) && (obj.iter < obj.simParameters.itMaxNR) ...
            && (rhsNorm > absTol)) || obj.iter == 0
          obj.iter = obj.iter + 1;
          %
          % Solve system with increment
          du = solve(obj.linSyst);
          % Reset global Jacobian and Rhs
          obj.linSyst.resetJacobianAndRhs(); 
          % Update tmpState
          obj.stateTmp.updateState(du,obj.dofManager);

          % Compute Rhs and Matrices of NonLinear models
          computeNLMatricesAndRhs(obj.linSyst,obj.stateTmp,obj.statek,obj.dt)

          % compute block Jacobian and block Rhs
          obj.linSyst.computeBlockJacobianAndRhs(delta_t);
          %
          applyBCandForces(obj.model, obj.grid, obj.bound, obj.material, ...
            obj.t, obj.linSyst, obj.stateTmp);
          
          % compute residual norm
          [~, rhsNorm] = computeRhsNorm(obj,obj.linSyst);
          if obj.simParameters.verbosity > 1
          fprintf('%d     %e\n',obj.iter,rhsNorm);
          end
        end
        %
        % Check for convergence
        flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);
        if flConv % Convergence
          obj.stateTmp.t = obj.t;
          % Print the solution, if needed
          if isPoromechanics(obj.model)
            obj.stateTmp.advanceState();
          end
          if isVariabSatFlow(obj.model)
            obj.stateTmp.updateSaturation()
          end
          if obj.t > obj.simParameters.tMax   % For Steady State
            printState(obj.printUtil,obj.stateTmp);
          else
            printState(obj.printUtil,obj.statek,obj.stateTmp);
          end
        end
        %
        % Manage next time step
        delta_t = manageNextTimeStep(obj,delta_t,flConv);
      end
        %
    end
  end

  methods (Access = private)
      function setNonLinearSolver(obj,models,interface)
          obj.models = models;
          obj.interface = interface;
          obj.nDom = numel(obj.models);
          obj.nInt = numel(obj.interface);
          obj.state = repmat(struct('prev',{},'curr',{}),obj.nDom,1);
          for i = 1:obj.nDom
              obj.state(i).prev = obj.models(i).State;
              obj.state(i).curr = copy(obj.state(i).prev);
          end
      end

      function [t, dt] = updateTime(obj,conv,dt)
          if obj.printUtil.modTime
              tmp = find(obj.t<obj.printUtil.timeList(),1,'first');
              if ~conv
                  t = min([obj.t + obj.dt, obj.t + dt, obj.printUtil.timeList(tmp)]);
              else
                  t = min([obj.t + obj.dt, obj.printUtil.timeList(tmp)]);
              end
          else
              t = obj.t + obj.dt;
          end
          dt = t - obj.t;
      end

      function [locRhsNorm, globRhsNorm] = computeRhsNorm(obj,syst)
          %Return maximum norm of all Rhs block
          nRhs = length(syst.dofm.subList);
          locRhsNorm = zeros(nRhs,1);
          for i = 1:nRhs
              locRhsNorm(i) = norm(syst.rhs{i}, obj.simParameters.pNorm);
          end
          globRhsNorm = sqrt(sum(locRhsNorm.^2));
      end


      function [dt] = manageNextTimeStep(obj,dt,flConv)
          if ~flConv   % Perform backstep
              transferState(obj.statek,obj.stateTmp);
              obj.t = obj.t - obj.dt;
              obj.tStep = obj.tStep - 1;
              dt = dt/obj.simParameters.divFac;
              obj.dt = obj.dt/obj.simParameters.divFac;  % Time increment chop
              if min(dt,obj.dt) < obj.simParameters.dtMin
                  if obj.simParameters.goOnBackstep == 1
                      flConv = 1;
                  elseif obj.simParameters.goOnBackstep == 0
                      error('Minimum time step reached')
                  end
              elseif obj.simParameters.verbosity > 0
                  fprintf('\n %s \n','BACKSTEP');
              end
          end
          if flConv % Go on if converged
              tmpVec = obj.simParameters.multFac;
              if isFlow(obj.model)
                  dpMax = max(abs(obj.stateTmp.pressure - obj.statek.pressure));
                  tmpVec = [tmpVec, (1+obj.simParameters.relaxFac)* ...
                      obj.simParameters.pTarget/(dpMax + obj.simParameters.relaxFac* ...
                      obj.simParameters.pTarget)];
                  %           if isVariabSatFlow(obj.model)
                  %             dSwMax = max(abs(obj.stateTmp.watSat-obj.statek.watSat));
                  %             tmpVec = [tmpVec, (1+obj.simParameters.relaxFac)* ...
                  %             obj.simParameters.sTarget/(dSwMax + obj.simParameters.relaxFac* ...
                  %             obj.simParameters.sTarget)];
                  %           end
              end
              obj.dt = min([obj.dt * min(tmpVec),obj.simParameters.dtMax]);
              obj.dt = max([obj.dt obj.simParameters.dtMin]);
              %         if obj.dt < obj.simParameters.dtMin
              %           obj.dt = obj.simParameters.dtMin;
              %         end
              %find interval on printUtil's timeList which contains obj.t
              %         tmp = find(obj.t<obj.printUtil.timeList(1),1,'first');
              %         if ((obj.t + obj.dt) > obj.printUtil.timeList(tmp))
              %             obj.dt = obj.printUtil.timeList(tmp) - obj.t;
              %         end
              %

              transferState(obj.stateTmp,obj.statek);
              %
              if ((obj.t + obj.dt) > obj.simParameters.tMax)
                  obj.dt = obj.simParameters.tMax - obj.t;
              end
          end
      end
  end
end