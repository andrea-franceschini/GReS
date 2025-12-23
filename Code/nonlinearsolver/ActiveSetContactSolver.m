classdef ActiveSetContactSolver < MultidomainFCSolver
  % Class for solving non linear contact problem

  properties (Access = private)
    %
    maxActiveSetIters = 10
    contactInterf
    itAS
    attemptedReset = false        % flag for attempt to save simulation resetting to stick
  end


  methods (Access = public)
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

    end



    function NonLinearLoop(obj)

      % Initialize the time step increment
      obj.dt = obj.simparams.dtIni;

      %

      for i = 1:obj.nDom
        obj.domains(i).applyDirVal(obj.t);
      end



      % Loop over time
      while obj.t < obj.simparams.tMax


        % Update the simulation time and time step ID
        absTol = obj.simparams.absTol;

        % reset active set iteration counter
        obj.itAS = 0;
        hasActiveSetChanged = true(numel(obj.contactInterf),1);

        while any(hasActiveSetChanged) && obj.itAS <= obj.maxActiveSetIters
          % outer active set loop

          obj.tStep = obj.tStep + 1;
          obj.t = obj.t + obj.dt;

          gresLog().log(0,'\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
          gresLog().log(0,'-----------------------------------------------------------\n');

          for i = 1:obj.nDom
            obj.domains(i).applyDirVal(obj.t);
          end


          gresLog().log(0,'\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
          gresLog().log(0,'Active set iteration n. %i \n', obj.itAS);
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

            obj.iter = obj.iter + 1;

            J = assembleJacobian(obj);

            du = solve(obj,J,rhs);

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

            % line search cut to avoid unnecessary iterations
             if obj.iter > 10 && rhsNorm > 1e1*rhsNormIt0
               gresLog().log(1," Line search cut due to non converging rhs \n")
               hasActiveSetChanged(:) = false;
               break
             end

          end % end newton loop

          if flConv % Newton Convergence
            
            % update the active set

            for i = 1:numel(obj.contactInterf)
              interf = obj.interfaces{obj.contactInterf(i)};
              hasActiveSetChanged(i) = updateActiveSet(interf);
            end

            obj.itAS = obj.itAS + 1;

            if any(hasActiveSetChanged) && obj.itAS == obj.maxActiveSetIters
              % force backstep
              fprintf('Reached maximum number of active set iterations \n');
              hasActiveSetChanged(:) = false;
              flConv = false;
            end
          end

          manageNextTimeStep(obj,flConv,hasActiveSetChanged);

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



    function manageNextTimeStep(obj,newtonConv,activeSetChanged)

      if ~newtonConv && ~obj.attemptedReset 
        reset = false(obj.nInterf,1);
        for i = 1:obj.nInterf
          reset(i) = resetConfiguration(obj.interfaces{i});
        end
        newtonConv = any(reset);
        if newtonConv
          obj.attemptedReset = true;
          activeSetChanged(:) = true;
        end
      end

      if ~newtonConv 
      
        obj.itAS = 0;

        % time step not converged
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
      end

      % case 2: newton convergence and active set not changed -> go to next
      % step and reset active set iters

      if newtonConv && ~any(activeSetChanged)  % converged time step

        obj.itAS = 0;

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

        % allow new survival attempts on new time steps
        obj.attemptedReset = false;

        return
      end

      % case 3: newton convergence but active set changed
      % keep current states, just go back in time
      if newtonConv && any(activeSetChanged)
        obj.t = obj.t - obj.dt;
        obj.tStep = obj.tStep - 1;
      end

    end



  end
end
