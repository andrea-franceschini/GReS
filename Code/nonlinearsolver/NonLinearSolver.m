classdef NonLinearSolver < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
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
    %
    t = 0
    tStep = 0
    iter
    dt
  end
  
  methods (Access = public)
      function obj = NonLinearSolver(symmod,simParam,dofManager,grid,mat,bc,prtUtil,stateIni,varargin)
      obj.setNonLinearSolver(symmod,simParam,dofManager,grid,mat,bc,prtUtil,stateIni,varargin);
    end

    function [simStat] = NonLinearLoop(obj)
      simStat = 1;
      if obj.elements.nCellsByType(2) > 0  % There is at least one hexahedron
        linSyst = Discretizer(obj.model,obj.simParameters,obj.dofManager,obj.grid,obj.material, ...
                  obj.GaussPts);
      else
        linSyst = Discretizer(obj.model,obj.simParameters,obj.dofManager,obj.grid,obj.material);
      end
      % Initialize the time step increment
      obj.dt = obj.simParameters.dtIni;  
      delta_t = obj.dt; % dynamic time step
      % Compute matrices of Linear models (once for the entire simulation)
      for i = 1:length(linSyst.fields)
          fld = linSyst.fields{i};
          if isLinear(linSyst.db(fld))
              linSyst.getField(fld).computeMat(obj.stateTmp,delta_t);
          end
      end
      flConv = true; %convergence flag
      %
      %
      % Loop over time
      while obj.t < obj.simParameters.tMax
        % Update the simulation time and time step ID
        obj.tStep = obj.tStep + 1;
        %new time update to fit the outTime list
        [obj.t, delta_t] = obj.updateTime(flConv, delta_t);
        %obj.t = obj.t + obj.dt;
        % Apply the Dirichlet condition value to the solution vector
        applyDirVal(obj.model, obj.bound, obj.t, obj.stateTmp);
        %
        fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,delta_t);
        fprintf('-----------------------------------------------------------\n');
        fprintf('Iter     ||rhs||\n');
        %
        % Compute Rhs and matrices of NonLinear models
        for i = 1:length(linSyst.fields)
            fld = linSyst.fields{i};
            if ~isLinear(linSyst.db(fld))
                linSyst.getField(fld).computeMat(obj.stateTmp,delta_t);
            end
             linSyst.getField(fld).computeRhs(obj.stateTmp,obj.statek,delta_t);
        end
        % compute block Jacobian and block Rhs
        linSyst.computeBlockJacobianAndRhs(delta_t);

        % Apply BCs to the block-wise system
        applyBCandForces(obj.model, obj.grid, obj.bound, obj.material, ...
          obj.t, linSyst, obj.stateTmp);

        % compute Rhs norm
        [locRhsNorm, rhsNorm] = computeRhsNorm(obj,linSyst);
        % consider output of local field rhs contribution

        tolWeigh = obj.simParameters.relTol*rhsNorm;
        obj.iter = 0;
        %
        fprintf('0     %e\n',rhsNorm);
        while ((rhsNorm > tolWeigh) && (obj.iter < obj.simParameters.itMaxNR) ...
            && (rhsNorm > obj.simParameters.absTol)) || obj.iter == 0
          obj.iter = obj.iter + 1;
          %
          % Solve system with increment
          du = solve(linSyst);
          % Reset global Jacobian and Rhs
          linSyst.resetJacobianAndRhs(); 
          % Update tmpState
          obj.stateTmp.updateState(du,obj.dofManager);

          % Compute Rhs and Matrices of NonLinear models
          % Loop trough problem fields (different compared to number of
          % blocks)
          for i = 1:length(linSyst.fields)
              fld = linSyst.fields{i};
              if ~isLinear(linSyst.db(fld))
                  linSyst.getField(fld).computeMat(obj.stateTmp,delta_t);
              end
              linSyst.getField(fld).computeRhs(obj.stateTmp,obj.statek,delta_t);
          end
          % compute block Jacobian and block Rhs (with theta method)
          linSyst.computeBlockJacobianAndRhs(delta_t);
          %
          applyBCandForces(obj.model, obj.grid, obj.bound, obj.material, ...
            obj.t, linSyst, obj.stateTmp);
          % compute residual norm
          [locRhsNorm, rhsNorm] = computeRhsNorm(obj,linSyst);
          fprintf('%d     %e\n',obj.iter,rhsNorm);
        end
        %
        % Check for convergence
        flConv = (rhsNorm < tolWeigh || rhsNorm < obj.simParameters.absTol);
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
    function setNonLinearSolver(obj,symmod,simParam,dofManager,grid,mat,bc,prtUtil,stateIni,data)
      obj.model = symmod;
      obj.simParameters = simParam;
      obj.dofManager = dofManager;
      obj.grid = grid;
      obj.mesh = grid.topology;
      obj.elements = grid.cells;
      obj.material = mat;
      obj.bound = bc;
%       obj.BCName = BName;
      obj.printUtil = prtUtil;
      obj.statek = stateIni;
      obj.stateTmp = copy(stateIni);
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
        assert(min(dt,obj.dt) >= obj.simParameters.dtMin,'Minimum time step reached');
        fprintf('\n %s \n','BACKSTEP');
      else % Go on if converged
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