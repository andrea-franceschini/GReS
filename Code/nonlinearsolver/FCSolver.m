classdef FCSolver < handle
  % Built in fully coupled solver
  % All equations are solved once using Newton Raphson
  
  properties (Access = private)
    %
    model
    simParameters
    dofManager
    grid
    mesh
    elements
%     faces
    material
    bound
%     BCName
    GaussPts
    printUtil
    statek
    stateTmp
    linSyst
    %
    t = 0
    tStep = 0
    iter
    dt
  end
  
  methods (Access = public)
      function obj = FCSolver(symmod,simParam,dofManager,grid,mat,bc,prtUtil,stateIni,linSyst,varargin)
      obj.setNonLinearSolver(symmod,simParam,dofManager,grid,mat,bc,prtUtil,stateIni,linSyst,varargin);
    end

   function [simStat] = NonLinearLoop(obj)
      simStat = 1;
      % Initialize the time step increment
      obj.dt = obj.simParameters.dtIni;  
      delta_t = obj.dt; % dynamic time step

      flConv = true; %convergence flag

      % Loop over time
      while obj.t < obj.simParameters.tMax
         absTol = obj.simParameters.absTol;
         % Update the simulation time and time step ID
         obj.tStep = obj.tStep + 1;
         %new time update to fit the outTime list
         [obj.t, delta_t] = obj.updateTime(flConv, delta_t);
         %obj.t = obj.t + obj.dt;
         % Apply the Dirichlet condition value to the solution vector
         obj.stateTmp = applyDirVal(obj.linSyst,obj.bound,obj.t,obj.stateTmp);
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
         % Loop over available linear models and compute the jacobian
         obj.stateTmp = computeMatricesAndRhs(obj.linSyst,obj.stateTmp,obj.statek,obj.dt);

         % Apply BCs to the blocks of the linear system
         applyBC(obj.linSyst,obj.bound,obj.t,obj.stateTmp);
         
         rhs = assembleRhs(obj.linSyst);

         % compute Rhs norm
         rhsNorm = norm(rhs,2);
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
            J = assembleJacobian(obj.linSyst);
            
            du = J\-rhs;


            % Update current model state
            obj.stateTmp = updateState(obj.linSyst,obj.stateTmp,du);

            % Compute Rhs and Matrices of NonLinear models
            obj.stateTmp = computeMatricesAndRhs(obj.linSyst,obj.stateTmp,obj.statek,obj.dt);

            % Apply BCs to the blocks of the linear system
            applyBC(obj.linSyst,obj.bound,obj.t,obj.stateTmp);

            rhs = assembleRhs(obj.linSyst);
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
            obj.stateTmp.t = obj.t;
            % Advance state of non linear models
            if isPoromechanics(obj.model)
               obj.stateTmp = Poromechanics.advanceState(obj.stateTmp);
            end
            
            if obj.t > obj.simParameters.tMax   % For Steady State
               printState_new(obj.printUtil,obj.stateTmp);
            else
               printState_new(obj.printUtil,obj.linSyst,obj.statek,obj.stateTmp);
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
     function setNonLinearSolver(obj,symmod,simParam,dofManager,grid,mat,bc,prtUtil,stateIni,linearSyst,data)
        obj.model = symmod;
        obj.simParameters = simParam;
        obj.dofManager = dofManager;
        obj.grid = grid;
        obj.mesh = grid.topology;
        obj.elements = grid.cells;
        obj.material = mat;
        obj.bound = bc;
        obj.printUtil = prtUtil;
        obj.statek = stateIni;
        obj.stateTmp = stateIni;
        obj.linSyst = linearSyst;
        if ~isempty(data)
           obj.GaussPts = data{1};
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

    function [dt] = manageNextTimeStep(obj,dt,flConv)
        if ~flConv   % Perform backstep
            obj.stateTmp = obj.statek;
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
            end
            obj.dt = min([obj.dt * min(tmpVec),obj.simParameters.dtMax]);
            obj.dt = max([obj.dt obj.simParameters.dtMin]);
            obj.statek = obj.stateTmp;
            %
            if ((obj.t + obj.dt) > obj.simParameters.tMax)
                obj.dt = obj.simParameters.tMax - obj.t;
            end
        end
    end
  end
end