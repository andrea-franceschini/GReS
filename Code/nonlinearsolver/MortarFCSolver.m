classdef MortarFCSolver < handle
   % Multi domain version of standard non linear solver class

   properties (Access = private)
      %
      simParameters
      nDom
      nInt
      %
      t = 0
      tStep = 0
      iter
      dt
   end

   properties (Access = public)
      state
      models
      meshGlue
   end

   methods (Access = public)
      function obj = MortarFCSolver(simParam,models,interfaces)
         obj.setNonLinearSolver(simParam,models,interfaces);
      end

      function NonLinearLoop(obj)
         % Initialize the time step increment
         obj.dt = obj.simParameters.dtIni;
         delta_t = obj.dt; % dynamic time step
         % Compute matrices of Linear models (once for the entire simulation)
         for i = 1:obj.nDom
            computeLinearMatrices(obj.meshGlue.model(i).Discretizer,obj.state(i).curr,obj.state(i).prev,obj.dt)
         end
         %
         flConv = true; %convergence flag
         %
         % Loop over time
         while obj.t < obj.simParameters.tMax
            % Update the simulation time and time step ID
            absTol = obj.simParameters.absTol;
            obj.tStep = obj.tStep + 1;
            %new time update to fit the outTime list

            [obj.t, delta_t] = obj.updateTime(flConv, delta_t);
            for i = 1:obj.nDom
               obj.meshGlue.model(i).Discretizer.resetJacobianAndRhs(); 
               % Apply the Dirichlet condition value to the solution vector
               if ~isempty(obj.meshGlue.model(i).BoundaryConditions)
                  applyDirVal(obj.meshGlue.model(i).ModelType,obj.meshGlue.model(i).BoundaryConditions,...
                     obj.t, obj.state(i).curr);
               end
               %
               % Compute Rhs and matrices of NonLinear models
               computeNLMatricesAndRhs(obj.meshGlue.model(i).Discretizer,...
                  obj.state(i).curr,obj.state(i).prev,obj.dt);
               
               % compute block Jacobian and block Rhs
               obj.meshGlue.model(i).Discretizer.computeBlockJacobianAndRhs(delta_t);
            end
            % Get unique multidomain solution system
            [J,rhs] = obj.meshGlue.getMDlinSyst();
            for i = 1:obj.nDom
               if ~isempty(obj.meshGlue.model(i).BoundaryConditions)
                  [J,rhs] = applyBCAndForces_MD(i, obj.meshGlue, obj.t, obj.state(i).curr, J, rhs);
               end
            end
            % compute Rhs norm
            rhsNorm = norm(rhs,2);

            if obj.simParameters.verbosity > 0
               fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,delta_t);
               fprintf('-----------------------------------------------------------\n');
            end

            if obj.simParameters.verbosity > 1
               fprintf('Iter     ||rhs||\n');
            end

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
               du = J\-rhs;
               % update solution vector for each model
               obj.updateStateMD(du);
               for i = 1:obj.nDom
                  obj.meshGlue.model(i).Discretizer.resetJacobianAndRhs();
                  % Update tmpState
                  computeNLMatricesAndRhs(obj.meshGlue.model(i).Discretizer,...
                     obj.state(i).curr,obj.state(i).prev,obj.dt);
                  % compute block Jacobian and block Rhs
                  obj.meshGlue.model(i).Discretizer.computeBlockJacobianAndRhs(delta_t);
               end


            % Get unique multidomain solution system
            [J,rhs] = obj.meshGlue.getMDlinSyst();

            % Apply BCs to global linear system
            for i = 1:obj.nDom
               if ~isempty(obj.meshGlue.model(i).BoundaryConditions)
                  [J,rhs] = applyBCAndForces_MD(i, obj.meshGlue, obj.t, obj.state(i).curr, J, rhs);
               end
            end

            % compute Rhs norm
            rhsNorm = norm(rhs,2);
               if obj.simParameters.verbosity > 1
                  fprintf('%d     %e\n',obj.iter,rhsNorm);
               end
            end
            %
            % Check for convergence
            flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);
            if flConv % Convergence
               for i = 1:obj.nDom
                  obj.state(i).curr.t = obj.t;
                  % Print the solution, if needed
                  if isPoromechanics(obj.meshGlue.model(i).ModelType)
                     obj.state(i).curr.advanceState();
                  end
                  if isVariabSatFlow(obj.meshGlue.model(i).ModelType)
                     obj.state(i).curr.updateSaturation()
                  end
               end
               if obj.t > obj.simParameters.tMax   % For Steady State
                  for i = 1:obj.nDom
                     printState(obj.meshGlue.model(i).OutState,obj.state(i).curr);
                  end
               else
                  for i = 1:obj.nDom
                     printState(obj.meshGlue.model(i).OutState,obj.state(i).prev,obj.state(i).curr);
                  end
               end
            end

            %
            % Manage next time step
            delta_t = manageNextTimeStep(obj,delta_t,flConv);
         end
         %
      end


      function NonLinearLoop_multiAquifer(obj)
          % Initialize the time step increment
          obj.dt = obj.simParameters.dtIni;
          delta_t = obj.dt; % dynamic time step
          fprintf('Computing Linear matrices \n')
          % Compute matrices of Linear models (once for the entire simulation)
          for i = 1:obj.nDom
              computeLinearMatrices(obj.meshGlue.model(i).Discretizer,obj.state(i).curr,obj.state(i).prev,obj.dt)
          end
          fprintf('Done Linear matrices computation \n')
          %
          flConv = true; %convergence flag
          %
          % Loop over time
          while obj.t < obj.simParameters.tMax
              % Update the simulation time and time step ID
              absTol = obj.simParameters.absTol;
              obj.tStep = obj.tStep + 1;
              %new time update to fit the outTime list
              [obj.t, delta_t] = obj.updateTime(flConv, delta_t);
              fprintf('Applying Dirichlet Boundary Conditions \n')
              for i = 1:obj.nDom
                  % Apply the Dirichlet condition value to the solution vector
                  %obj.meshGlue.model(i).Discretizer.resetJacobianAndRhs();
                  if ~isempty(obj.meshGlue.model(i).BoundaryConditions)
                      applyDirVal(obj.meshGlue.model(i).ModelType,obj.meshGlue.model(i).BoundaryConditions,...
                          obj.t, obj.state(i).curr);
                  end
                  %
                  % Compute Rhs and matrices of NonLinear models
                  computeNLMatricesAndRhs(obj.meshGlue.model(i).Discretizer,...
                      obj.state(i).curr,obj.state(i).prev,obj.dt);
                  %
                  % compute block Jacobian and block Rhs
                  obj.meshGlue.model(i).Discretizer.computeBlockJacobianAndRhs(delta_t);
              end
              fprintf('Done applying Boundary Conditions \n')
              %
              fprintf('Setting up full linear system \n')
              tic
              % Get unique multidomain solution system
              [J,rhs] = obj.meshGlue.getMDlinSyst();
              fprintf('Done setting up full linear system with t = %1.4f \n',toc)
              % Apply BCs to global linear system
              fprintf('Applying BCs to full linear system \n')
              tic
              for i = 1:obj.nDom
                  if ~isempty(obj.meshGlue.model(i).BoundaryConditions)
                      [J,rhs] = applyBCAndForces_MD(i, obj.meshGlue, obj.t, obj.state(i).curr, J, rhs);
                  end
              end
              fprintf('Done applying BCs to full linear system t = %1.4f\n',toc)
              % compute Rhs norm
              rhsNorm = norm(rhs,2);

              if obj.simParameters.verbosity > 0
                  fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,delta_t);
                  fprintf('-----------------------------------------------------------\n');
              end

              if obj.simParameters.verbosity > 1
                  fprintf('Iter     ||rhs||\n');
              end

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
                  fprintf('Solving linear system \n')
                  tic
                  du = J\-rhs;
%                   str = save('solution.mat',"du");
%                   du = str.du;
                  fprintf('Done solving linear system t = %1.4f \n',toc)
                  % update solution vector for each model
                  obj.updateStateMD(du);
                  for i = 1:obj.nDom
                      obj.meshGlue.model(i).Discretizer.resetJacobianAndRhs();
                      % Update tmpState
                      computeNLMatricesAndRhs(obj.meshGlue.model(i).Discretizer,...
                          obj.state(i).curr,obj.state(i).prev,obj.dt);
                      % compute block Jacobian and block Rhs
                      obj.meshGlue.model(i).Discretizer.computeBlockJacobianAndRhs(delta_t);
                  end


                  % Get unique multidomain solution system
                  [J,rhs] = obj.meshGlue.getMDlinSyst();

                  % Apply BCs to global linear system
                  for i = 1:obj.nDom
                      if ~isempty(obj.meshGlue.model(i).BoundaryConditions)
                          [J,rhs] = applyBCAndForces_MD(i, obj.meshGlue, obj.t, obj.state(i).curr, J, rhs);
                      end
                  end

                  % compute Rhs norm
                  rhsNorm = norm(rhs,2);
                  if obj.simParameters.verbosity > 1
                      fprintf('%d     %e\n',obj.iter,rhsNorm);
                  end
              end
              %
              % Check for convergence
              %flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);
              flConv = true;
              if flConv % Convergence
                  for i = 1:obj.nDom
                      obj.state(i).curr.t = obj.t;
                      % Print the solution, if needed
                      if isPoromechanics(obj.meshGlue.model(i).ModelType)
                          obj.state(i).curr.advanceState();
                      end
                      if isVariabSatFlow(obj.meshGlue.model(i).ModelType)
                          obj.state(i).curr.updateSaturation()
                      end
                  end
                  if obj.t > obj.simParameters.tMax || obj.t == obj.simParameters.tMax % For Steady State
                      for i = 1:obj.nDom
                          printState(obj.meshGlue.model(i).OutState,obj.state(i).curr);
                      end
                  else
                      for i = 1:obj.nDom
                          printState(obj.meshGlue.model(i).OutState,obj.state(i).prev,obj.state(i).curr);
                      end
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
       function setNonLinearSolver(obj,simParam,models,mG)
           obj.simParameters = simParam;
           obj.models = models;
           obj.meshGlue = mG;
           obj.nDom = numel(obj.models);
           obj.nInt = numel(obj.meshGlue.interfaces);
           obj.state = repmat(struct('prev',{},'curr',{}),obj.nDom,1);
           for i = 1:obj.nDom
               obj.state(i).prev = obj.meshGlue.model(i).State;
               obj.state(i).curr = copy(obj.state(i).prev);
           end
       end

       function [t, dt] = updateTime(obj,conv,dt)
           t = obj.simParameters.tMax;
           told = t;
           for i = 1:obj.nDom
               if obj.meshGlue.model(i).OutState.modTime
                   tmp = find(obj.t<obj.meshGlue.model(i).outState.timeList(),1,'first');
                   if ~conv
                       t = min([obj.t + obj.dt, obj.t + dt, obj.meshGlue.model(i).OutState.timeList(tmp)]);
                   else
                       t = min([obj.t + obj.dt, obj.meshGlue.model(i).OutState.timeList(tmp)]);
                   end
               else
                   t = obj.t + obj.dt;
               end
               if t > told
                   t = told;
               end
           end
           dt = t - obj.t;
       end

       function out = computeRhsNorm(obj)
           %Return maximum norm of the entire domain
           rhsNorm = zeros(obj.nDom,1);
           for i = 1:obj.nDom
               nRhs = length(obj.meshGlue.model(i).DoFManager.subList);
               rhsNorm_loc = zeros(nRhs,1);
               for j = 1:nRhs
                   rhsNorm_loc(j) = norm(obj.meshGlue.model(i).Discretizer.rhs{j}, obj.simParameters.pNorm);
               end
               rhsNorm(i) = sqrt(sum(rhsNorm_loc.^2));
           end
           out = norm(rhsNorm);
       end

       function updateStateMD(obj,du)
           % update the MD current state with the computed increment
           for i = 1:numel(obj.meshGlue.MD_struct)
               domID = obj.meshGlue.MD_struct(i).dom;
               ph = obj.meshGlue.MD_struct(i).physic;
               ent_dof = obj.models(domID).DoFManager.ent2field(ph, ...
                   obj.meshGlue.MD_struct(i).entities);
               if any(strcmp(obj.meshGlue.MD_struct(i).type,["inner","master"]))
                   switch ph
                       case 'Poro'
                           obj.state(domID).curr.dispCurr(ent_dof) = obj.state(domID).curr.dispCurr(ent_dof) + du(obj.meshGlue.getDofs_MD(i));
                       case 'SPFlow'
                           obj.state(domID).curr.pressure(ent_dof) = obj.state(domID).curr.pressure(ent_dof) + du(obj.meshGlue.getDofs_MD(i));
                   end
               else
                   idM = getMaster(obj.meshGlue,i);
                   switch ph
                       case 'Poro'
                           E = expandMortarOperator(obj.meshGlue,i);
                           obj.state(domID).curr.dispCurr(ent_dof) = obj.state(domID).curr.dispCurr(ent_dof) + E*du(obj.meshGlue.getDofs_MD(idM));
                       case 'SPFlow'
                           obj.state(domID).curr.pressure(ent_dof) = obj.state(domID).curr.pressure(ent_dof) + E*du(obj.meshGlue.getDofs_MD(idM));
                   end
               end
           end
       end


       function [dt] = manageNextTimeStep(obj,dt,flConv)
           if ~flConv   % Perform backstep
               for i = 1:obj.nDom
                   transferState(obj.state(i).curr,obj.state(i).prev);
               end
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
               for i =1:obj.nDom
                   if isFlow(obj.meshGlue.model(i).ModelType)
                       dpMax = max(abs(obj.state(i).curr.pressure - obj.state(i).prev.pressure));
                       tmpVec = [tmpVec, (1+obj.simParameters.relaxFac)* ...
                           obj.simParameters.pTarget/(dpMax + obj.simParameters.relaxFac* ...
                           obj.simParameters.pTarget)];
                   end
               end
               obj.dt = min([obj.dt * min(tmpVec),obj.simParameters.dtMax]);
               obj.dt = max([obj.dt obj.simParameters.dtMin]);
               for i = 1:obj.nDom
                   transferState(obj.state(i).curr,obj.state(i).prev);
               end
               %
               if ((obj.t + obj.dt) > obj.simParameters.tMax)
                   obj.dt = obj.simParameters.tMax - obj.t;
               end
           end
       end
   end
end