classdef ActiveSetContactSolver < MultidomainFCSolver
  % Class for solving non linear contact problem

  properties (Access = private)
    %
    maxActiveSetIters = 10
    contactInterf
  end


  methods (Access = public)
<<<<<<< HEAD
    function obj = ActiveSetContactSolver(domains,interfaces,varargin)
      obj@MultidomainFCSolver(domains,interfaces)
      setContactInterfaces(obj);
      if ~isempty(varargin)
        obj.maxActiveSetIters = varargin{1};
      end
=======
    function obj = ActiveSetContactSolver(simparams,domains,interfaces,varargin)
      
      obj@MultidomainFCSolver(simparams,domains,interfaces)

      % find which of the available interfaces have an active set
      obj.contactInterf = find(cellfun(@(c) isprop(c,'activeSet'), obj.interfaces));
      if isempty(obj.contactInterf)
        error("Call to ActiveSetContactSolver but no activeSet property " + ...
          " is available for interfaces in the model");
      end
      if ~isempty(varargin)
        obj.maxActiveSetIters = varargin{1};
      end

>>>>>>> origin/main
    end



    function NonLinearLoop(obj)

      % Initialize the time step increment
<<<<<<< HEAD
      obj.dt = obj.simParameters.dtIni;
      delta_t = obj.dt; % dynamic time step

      %
      flConv = true; %convergence flag

      % initialize the state object
      applyDirVal(obj);
      for i = 1:obj.nDom
        obj.state(i).curr = obj.domains(i).state;
        obj.state(i).prev =  copy(obj.state(i).curr);
      end


      % Loop over time
      while obj.t < obj.simParameters.tMax


        % Update the simulation time and time step ID
        absTol = obj.simParameters.absTol;
        obj.tStep = obj.tStep + 1;
        %new time update to fit the outTime list

        if obj.simParameters.verbosity > 0
          fprintf('\n-----------------------------------------------------------\n');
          fprintf('TIME STEP %i\n',obj.tStep)
=======
      obj.dt = obj.simparams.dtIni;

      %

      for i = 1:obj.nDom
        obj.domains(i).applyDirVal(obj.t);
      end



      % Loop over time
      while obj.t < obj.simparams.tMax


        % Update the simulation time and time step ID
        absTol = obj.simparams.absTol;
        obj.tStep = obj.tStep + 1;
        obj.t = obj.t + obj.dt;

        gresLog().log(0,'\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
        gresLog().log(0,'-----------------------------------------------------------\n');

        for i = 1:obj.nDom
          obj.domains(i).applyDirVal(obj.t);
>>>>>>> origin/main
        end

        % reset active set iteration counter
        itAS = 0;
<<<<<<< HEAD

        hasActiveSetChanged = true;

        [obj.t, delta_t] = obj.updateTime(flConv, delta_t);


        applyDirVal(obj);

        while hasActiveSetChanged && itAS <= obj.maxActiveSetIters
          % outer active set loop

          if obj.simParameters.verbosity > 0
            fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,delta_t);
          end

          if obj.simParameters.verbosity > 0
            fprintf('Active set iteration n. %i \n', itAS)
          end

          if obj.simParameters.verbosity > 1
            fprintf('Iter     ||rhs||\n');
          end

          computeMatricesAndRhs(obj);
          applyBC(obj);
          rhs = assembleRhs(obj);
          rhsNorm = norm(cell2mat(rhs),2);

          tolWeigh = obj.simParameters.relTol*rhsNorm;
          obj.iter = 0;
          %
          if obj.simParameters.verbosity > 1
            fprintf('0     %e\n',rhsNorm);
          end

          while ((rhsNorm > tolWeigh) && (obj.iter < obj.simParameters.itMaxNR) ...
              && (rhsNorm > absTol)) || obj.iter == 0
=======
        hasActiveSetChanged = true(numel(obj.contactInterf),1);

        while any(hasActiveSetChanged) && itAS <= obj.maxActiveSetIters
          % outer active set loop

          gresLog().log(0,'\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
          gresLog().log(0,'Active set iteration n. %i \n', itAS);
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
          obj.iter = 0;
          %
          gresLog().log(1,'0     %e     %e\n',rhsNorm,rhsNorm/rhsNormIt0);

          flConv = false;

          while (~flConv) && (obj.iter < obj.simparams.itMaxNR)
>>>>>>> origin/main

            obj.iter = obj.iter + 1;

            J = assembleJacobian(obj);

            du = solve(obj,J,rhs);

<<<<<<< HEAD
            % update primary variables and multipliers
            updateState(obj,du);

            computeMatricesAndRhs(obj);
            applyBC(obj);
            rhs = assembleRhs(obj);
            rhsNorm = norm(cell2mat(rhs),2);


            if obj.simParameters.verbosity > 1
              fprintf('%d     %e\n',obj.iter,rhsNorm);
            end
          end
          %
          % Check for convergence
          flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);

          % set active set not changed by default
          hasActiveSetChanged = false;

          if flConv % Newton Convergence
            % Advance state of non linear models

            for i = obj.contactInterf
              hasActiveSetChanged = updateActiveSet(obj.interfaces{i});
=======
            c = 0;

            for i = 1:obj.nDom
              nDof = obj.domains(i).getNumbDoF();
              sol = du(c+1:c+nDof);
              obj.domains(i).updateState(sol);
              c = c + nDof;
            end

            for i = 1:obj.nInterf
              nDof = obj.interfaces{i}.getNumbDoF();
              sol = du(c+1:c+nDof);
              obj.interfaces{i}.updateState(sol);
              c = c + nDof;
            end

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
            gresLog().log(1,'%d     %e     %e\n',obj.iter,rhsNorm,rhsNorm/rhsNormIt0);


            %
            % Check for convergence
            flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);

          end % end newton loop

          if flConv % Newton Convergence
            
            % update the active set

            for i = 1:numel(obj.contactInterf)
              interf = obj.interfaces{obj.contactInterf(i)};
              hasActiveSetChanged(i) = updateActiveSet(interf);
>>>>>>> origin/main
            end

            itAS = itAS + 1;

<<<<<<< HEAD
            if ~hasActiveSetChanged
            for i = 1:obj.nDom
              obj.state(i).curr.t = obj.t;
              if isPoromechanics(obj.domains(i).model)
                obj.domains(i).getSolver('Poromechanics').advanceState();
              end
            end
            end

%             % this happen also if obj.maxActiveSetIters = 0, break loop
%             if itAS >= obj.maxActiveSetIters 
%               obj.state(1).curr.t = obj.t;
%               obj.state(2).curr.t = obj.t;
%               printState(obj);
%               delta_t = manageNextTimeStep(obj,delta_t,flConv,hasActiveSetChanged);
%               break
%             end

            if itAS == obj.maxActiveSetIters
              fprintf('Reached maximum number of active set iterations \n');
              hasActiveSetChanged = false;
=======
            if any(hasActiveSetChanged) && itAS == obj.maxActiveSetIters
              % force backstep
              fprintf('Reached maximum number of active set iterations \n');
              hasActiveSetChanged(:) = false;
>>>>>>> origin/main
              flConv = false;
            end
          end

<<<<<<< HEAD
          if flConv
            printState(obj);
          end

          delta_t = manageNextTimeStep(obj,delta_t,flConv,hasActiveSetChanged);
=======
          manageNextTimeStep(obj,flConv,hasActiveSetChanged);
>>>>>>> origin/main

        end % outer active set loop
      end % time marching
      %
    end 

  end



  methods (Access = protected)

    

    function setContactInterfaces(obj)

      % get contact interfaces
      for i = 1:obj.nInterf
        if isa(obj.interfaces{i},"ContactMortar")
          obj.contactInterf = [obj.contactInterf,i];
        end
      end
    end

%     function J = assembleJacobian(obj)
%       % assemble blocks of jacobian matrix for multidomain system
%       %[N,Nf,Ni] = deal(obj.systSize(1),obj.systSize(2),obj.systSize(3));
%       J = cell(obj.systSize(1));
%       k = 0;
%       % populate jacobian with inner domain blocks
%       for iDom = 1:obj.nDom
%         discr = obj.domains(iDom);
%         J(k+1:k+obj.nfldDom(iDom),k+1:k+obj.nfldDom(iDom)) = ...
%           discr.assembleJacobian();
%         for iFld = 1:obj.nfldDom(iDom)
%           fld = discr.fields(iFld);
%           for iI = discr.interfaceList
%             jj = obj.systSize(2)+iI;
%             [J{iFld+k,jj},J{jj,iFld+k}] = getJacobian(...
%               obj.interfaces{iI},fld,iDom);
%           end
%         end
%         k = k+obj.nfldDom(iDom);
%       end
% 
%       % assembly multiplier blocks and static condensation terms
%       for iI = 1:obj.nInterf
%         %
%         interf = obj.interfaces{iI};
%         if isa(interf,'MeshGlueDual')
%           id = interf.idDomain;
%           if isempty(J{id(1),id(2)})
%             J{id(1),id(2)} = interf.Jcoupling';
%           else
%             J{id(1),id(2)} =  J{id(1),id(2)} + interf.Jcoupling';
%           end
%           if isempty(J{id(2),id(1)})
%             J{id(2),id(1)} = interf.Jcoupling;
%           else
%             J{id(2),id(1)} = J{id(2),id(1)} + interf.Jcoupling;
%           end
%         else
%           jj = obj.systSize(2)+iI;
%           [J{jj,jj}] = getJacobian(...
%             obj.interfaces{iI},fld);
%         end
%       end
% 
%     end

<<<<<<< HEAD
    function [dt] = manageNextTimeStep(obj,dt,newtonConv,activeSetChanged)
      if ~newtonConv    % time step not converged
        dt = dt/obj.simParameters.divFac;
        obj.dt = obj.dt/obj.simParameters.divFac;  % Time increment chop

        if min(dt,obj.dt) < obj.simParameters.dtMin
          if obj.simParameters.goOnBackstep == 1
            newtonConv = true;
          elseif obj.simParameters.goOnBackstep == 0
            error('Minimum time step reached')
          end
        else
          if obj.simParameters.verbosity > 0
            fprintf('\n %s \n','BACKSTEP');
            goBackState(obj);
            obj.t = obj.t - obj.dt*obj.simParameters.divFac;
            obj.tStep = obj.tStep - 1;
          end
        end
      end
      if newtonConv && ~activeSetChanged  % converged time step
        tmpVec = obj.simParameters.multFac;
        for i = 1:obj.nDom
          if isFlow(obj.domains(i).model)
            pnew = obj.state(i).curr.data.pressure;
            pold = obj.state(i).prev.data.pressure;
            dpMax = max(abs(pnew-pold));
            tmpVec = [tmpVec, (1+obj.simParameters.relaxFac)* ...
              obj.simParameters.pTarget/(dpMax + obj.simParameters.relaxFac* ...
              obj.simParameters.pTarget)];
          end
        end
        obj.dt = min([obj.dt * min(tmpVec),obj.simParameters.dtMax]);
        obj.dt = max([obj.dt obj.simParameters.dtMin]);
        goOnState(obj);
        %
        if ((obj.t + obj.dt) > obj.simParameters.tMax)
          obj.dt = obj.simParameters.tMax - obj.t;
        end
      end
=======


    function manageNextTimeStep(obj,newtonConv,activeSetChanged)

      % case 1: newton not converged (perform backstep)
      if ~newtonConv 

        for i = 1:obj.nDom
          goBackState(obj.domains(i));
        end

        for i = 1:obj.nInterf
          goBackState(obj.interfaces{i});
        end

        % time step not converged
        dt = dt/obj.simparams.divFac;
        obj.dt = obj.dt/obj.simparams.divFac;  % Time increment chop

        obj.t = obj.t - obj.dt;
        obj.tStep = obj.tStep - 1;
        obj.dt = obj.dt/obj.simparams.divFac;

        if obj.dt < obj.simparams.dtMin
          error('Minimum time step reached')
        else
          gresLog().log(0,'\n %s \n','BACKSTEP')
        end

        return 
      end

      % case 2: newton convergence and active set not changed -> go to next
      % steo

      if newtonConv && ~any(activeSetChanged)  % converged time step

        for i = 1:obj.nDom
          dom = obj.domains(i);
          dom.state.t = obj.t;
          printState(dom);
          advanceState(dom);
        end

        for i = 1:obj.nInterf
          interf = obj.interfaces{i};
          interf.state.t = obj.t;
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

        return
      end

      % case 3: newton convergence but active set changed
      % do nothing, keep current state() and stateOld()

>>>>>>> origin/main
    end



  end
end
