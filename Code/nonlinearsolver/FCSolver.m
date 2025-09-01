classdef FCSolver < handle
  % Built in fully coupled solver
  % All equations are solved once using Newton Raphson
  
  properties (Access = private)
    %
    statek
    stateTmp
    domain
    %
    t = 0
    tStep = 0
    iter
    dt
    toGrow
  end

  properties (Access = public)
      solStatistics
  end
  
  methods (Access = public)
      function obj = FCSolver(linSyst,varargin)
        obj.setNonLinearSolver(linSyst);
        saveStasticts = false(3,1);
        for k = 1:(nargin-1)/2
            pos = 2*(k-1)+1;
            switch varargin{pos}
                case 'SaveRelError'
                    saveStasticts(1) = varargin{pos+1};
                case 'SaveAbsError'
                    saveStasticts(2) = varargin{pos+1};
                case 'SaveBStepInf'
                    saveStasticts(3) = varargin{pos+1};
                case 'SaveConverg'
                    saveStasticts(:) = varargin{pos+1};
                otherwise
            end
        end
        obj.solStatistics = SolverStatistics(linSyst.simparams.itMaxNR,linSyst.simparams.relTol,linSyst.simparams.absTol,saveStasticts);
        % obj.toGrow = GrowningDomain(obj.linSyst,obj.bound);
    end

    function [simStat] = NonLinearLoop(obj)
      simStat = 1;
      % Initialize the time step increment
      obj.dt = obj.domain.simparams.dtIni;  
      delta_t = obj.dt; % dynamic time step

      flConv = true; % convergence flag

      % [obj.statek,obj.stateTmp]=obj.toGrow.addCell(obj.linSyst,1,2,6,obj.statek,obj.stateTmp);
      % [obj.statek,obj.stateTmp]=obj.toGrow.addCells(obj.linSyst,1,[2 4 6 8],6,obj.statek,obj.stateTmp);

      % Loop over time
      while obj.t < obj.domain.simparams.tMax
         absTol = obj.domain.simparams.absTol;
         residual = zeros(obj.domain.simparams.itMaxNR+1,2);

         % =======
         % % add cells
         % % if obj.t>50 && obj.t<65
         % %    [obj.statek,obj.stateTmp]=obj.toGrow.addCells(obj.linSyst,1,[2 4 6 8],6,obj.statek,obj.stateTmp);
         % % end
         % >>>>>>> 1dfffa00097f21a2e1d34699913ab58ea5431391

         % Update the simulation time and time step ID
         obj.tStep = obj.tStep + 1;
         %new time update to fit the outTime list
         [obj.t, delta_t] = obj.updateTime(flConv, delta_t);

         % Apply the Dirichlet condition value to the solution vector
         applyDirVal(obj.domain,obj.t);
         %
         if obj.domain.simparams.verbosity > 0
            fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,delta_t);
            fprintf('-----------------------------------------------------------\n');
         end
         if obj.domain.simparams.verbosity > 1
            % fprintf('Iter     ||rhs||\n');
            fprintf('Iter     ||rhs||     ||rhs||/||rhs_0||\n');
         end

         % Compute Rhs and matrices of NonLinear models
         % Loop over available linear models and compute the jacobian
         computeMatricesAndRhs(obj.domain,obj.statek,obj.dt);

         % Apply BCs to the blocks of the linear system
         applyBC(obj.domain,obj.t);
         
         rhs = assembleRhs(obj.domain);

         % compute Rhs norm
         rhsNorm = norm(cell2mat(rhs),2);
         rhsNormIt0 = rhsNorm;
         residual(1,1) = rhsNormIt0;
         residual(1,2) = 1.;

         % consider output of local field rhs contribution
         tolWeigh = obj.domain.simparams.relTol*rhsNorm;
         obj.iter = 0;
         %
         if obj.domain.simparams.verbosity > 1
            % fprintf('0     %e\n',rhsNorm);
            fprintf('0     %e     %e\n',rhsNorm,rhsNorm/rhsNormIt0);
         end
         while ((rhsNorm > tolWeigh) && (obj.iter < obj.domain.simparams.itMaxNR) ...
               && (rhsNorm > absTol)) || obj.iter == 0
            obj.iter = obj.iter + 1;
            %
            % Solve system with increment
            J = assembleJacobian(obj.domain);            
            du = FCSolver.solve(J,rhs);

            % =======
            % J = assembleJacobian(obj.linSyst);
            % du = J\-rhs;
            % >>>>>>> 1dfffa00097f21a2e1d34699913ab58ea5431391

            % Update current model state
            updateState(obj.domain,du);

            % Compute Rhs and Matrices of NonLinear models
            computeMatricesAndRhs(obj.domain,obj.statek,obj.dt);

            % Apply BCs to the blocks of the linear system
            applyBC(obj.domain,obj.t);

            rhs = assembleRhs(obj.domain);
            % compute Rhs norm
            rhsNorm = norm(cell2mat(rhs),2);
            residual(obj.iter+1,1)=rhsNorm;
            residual(obj.iter+1,2)=rhsNorm/rhsNormIt0;

            if obj.domain.simparams.verbosity > 1
               fprintf('%d     %e     %e\n',obj.iter,residual(obj.iter+1,1),residual(obj.iter+1,2));
            end
         end
         %
         % Check for convergence
         flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);
         % flConv = (rhsNorm < tolWeigh);
         if flConv % Convergence
            obj.stateTmp.t = obj.t;
            obj.domain.state.t = obj.t;
            % Advance state of non linear models
            if isPoromechanics(obj.domain.model)
               advanceState(getSolver(obj.domain,'Poromechanics'));
            end
            
            if obj.t > obj.domain.simparams.tMax   % For Steady State
               printState(obj.domain);
               % printState(obj.printUtil,obj.bound,obj.stateTmp);
            else
               printState(obj.domain,obj.statek);
               % printState(obj.printUtil,obj.linSyst,obj.bound,obj.statek,obj.stateTmp);
            end
            obj.solStatistics.saveIt(obj.t,residual(1:obj.iter+1,1),residual(1:obj.iter+1,2));
         else
             obj.solStatistics.saveBackIt();
         end
         %
         % Manage next time step
         delta_t = manageNextTimeStep(obj,delta_t,flConv);
      end
      %
    end
  end
  
  methods (Access = private)
     function setNonLinearSolver(obj,linearSyst)
        assert(numel(linearSyst)==1,['Multidomain models is not supported' ...
          'by FCSolver class'])
        obj.stateTmp = linearSyst.state;
        obj.statek = copy(obj.stateTmp);
        obj.domain = linearSyst;
     end

    function [t, dt] = updateTime(obj,conv,dt)
        if obj.domain.outstate.modTime
            tmp = find(obj.t<obj.domain.outstate.timeList(),1,'first');
            if ~conv       
                t = min([obj.t + obj.dt, obj.t + dt, obj.domain.outstate.timeList(tmp)]);
            else 
                t = min([obj.t + obj.dt, obj.domain.outstate.timeList(tmp)]);
            end
        else
        t = obj.t + obj.dt;
        end
        dt = t - obj.t; 
    end

    function [dt] = manageNextTimeStep(obj,dt,flConv)
        if ~flConv   % Perform backstep
            obj.stateTmp = copy(obj.statek);
            obj.t = obj.t - obj.dt;
            obj.tStep = obj.tStep - 1;
            dt = dt/obj.domain.simparams.divFac;
            obj.dt = obj.dt/obj.domain.simparams.divFac;  % Time increment chop
            if min(dt,obj.dt) < obj.domain.simparams.dtMin
                if obj.domain.simparams.goOnBackstep == 1
                    flConv = 1;
                elseif obj.domain.simparams.goOnBackstep == 0
                    error('Minimum time step reached')
                end
            elseif obj.domain.simparams.verbosity > 0
                fprintf('\n %s \n','BACKSTEP');
            end
        end
        if flConv % Go on if converged
            tmpVec = obj.domain.simparams.multFac;
            if isFlow(obj.domain.model)
                dpMax = max(abs(obj.stateTmp.data.pressure - obj.statek.data.pressure));
                tmpVec = [tmpVec, (1+obj.domain.simparams.relaxFac)* ...
                  obj.domain.simparams.pTarget/(dpMax + obj.domain.simparams.relaxFac* ...
                  obj.domain.simparams.pTarget)];
            end
            obj.dt = min([obj.dt * min(tmpVec),obj.domain.simparams.dtMax]);
            obj.dt = max([obj.dt obj.domain.simparams.dtMin]);
            obj.statek = copy(obj.stateTmp);
            %
            if ((obj.t + obj.dt) > obj.domain.simparams.tMax)
              obj.dt = obj.domain.simparams.tMax - obj.t;
            end
        end
    end
  end

  methods (Static)
    function sol = solve(J,rhs)
      J = FCSolver.cell2matJac(J);
      rhs = cell2mat(rhs);
      sol = J\(-rhs);
    end

    function mat = cell2matJac(mat)
      % get size of cell matrix
      n = size(mat,1);
      szr = zeros(n,1);
      szc = zeros(n,1);
      for i = 1:n
        for j = 1:n
          if isempty(mat{i,j})
            continue
          else
            szr(i) = size(mat{i,j},1);
            szc(j) = size(mat{i,j},2);
          end
        end
      end

      % populate empty blocks
      for i = 1:n
        for j = 1:n
          if isempty(mat{i,j})
            mat{i,j} = sparse(szr(i),szc(j));
          end
        end
      end

      mat = cell2mat(mat);
    end
  end
end