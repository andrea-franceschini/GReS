classdef NonLinearSolver < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = private)
%     itMax = 10
%     relTol = 1.e-6
%     absTol = 1.e-10
%     pNorm = 2
%     theta = 1
%     dtIni = 0.5    % 0.5
%     dtMin = 1.e-5
%     dtMax = 1
%     tMax = 0.2
%     multFac = 1.1
%     divFac = 2
%     relaxFac = 0
%     pTarget = 1000000000000
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
      % Compute H and P matrices for flow simulations and the gravity
      % contribution
      if isSinglePhaseFlow(obj.model)
        linSyst.computeSPFMatrices();
%         linSyst.computeFlowRHSGravTerm();
      end
%       t  = 0;
%       tStep = 0;
      % Initialize the time step increment
      obj.dt = obj.simParameters.dtIni;  
      delta_t = obj.dt; % dynamic time step
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
        applyDirVal(obj.bound, obj.t, obj.stateTmp);
        %
        fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,delta_t);
        fprintf('-----------------------------------------------------------\n');
        fprintf('Iter     ||rhs||\n');
        %
        %Call different solvers' methods separately
        % ---> SPFLOW <--- (matrices computed once before entering time stepping scheme)
        if isSinglePhaseFlow(obj.model)
            linSyst.computeFlowJacobian_Test(delta_t);
            mu = obj.material.getFluid().getDynViscosity();
            linSyst.computeFlowRHS_test(obj.statek,obj.stateTmp,delta_t,1/mu);
        end

        % ---> POROMECHANICS <--- (if theta method is not required)
        if  isPoromechanics(obj.model) && ~isSinglePhaseFlow(obj.model)
            linSyst.computePoroSyst_Test(obj.stateTmp,delta_t);
        end

        % ---> BIOT COUPLING <---
        if isPoromechanics(obj.model) && isSinglePhaseFlow(obj.model)
            linSyst.computeBiotSyst(delta_t,obj.statek,obj.stateTmp);
        end
 

        % if isPoromechanics(obj.model) && isSinglePhaseFlow(obj.model)
        %     % compute Jacobian and residual of coupled hydro-mechanics
        %     % problem       
        %   linSyst.computeCoupleSyst(obj.simParameters.theta,delta_t,obj.statek,obj.stateTmp)
        % elseif isPoromechanics(obj.model) 
        %   % Compute Jacobian and residual of the poromechanical problem
        %   linSyst.computePoroSyst(obj.stateTmp,delta_t);
        % elseif isSinglePhaseFlow(obj.model)
        %   % Compute Jacobian and residual of the flow problem 
        %   linSyst.computeFlowJacobian(delta_t);
        %   mu = obj.material.getMaterial(obj.mesh.nCellTag+1).getDynViscosity();
        %   linSyst.computeFlowRHS(obj.statek,obj.stateTmp,delta_t,1/mu);
        % elseif isVariabSatFlow(obj.model)
        %   linSyst.computeVSFMatricesAndRHS(obj.statek,obj.stateTmp,delta_t);
        % end
        %

        % Apply Neu and Dir conditions
        linSyst.buildGlobalJacobianAndRhs();
        applyBCandForces_test(obj.model, obj.grid, obj.bound, obj.material, ...
          obj.t, linSyst, obj.stateTmp);
%         obj.bound.applyBCNeu(linSyst);
        rhsNorm = norm(linSyst.rhs,obj.simParameters.pNorm);
        %rhsNorm = computeRhsNorm(obj,linSyst);
        linSyst.resetJacobianAndRhs();
%         obj.bound.applyBCDir(linSyst);
%         obj.bound.applyBCNeu(linSyst);
        tolWeigh = obj.simParameters.relTol*rhsNorm;
        obj.iter = 0;
        %
        fprintf('0     %e\n',rhsNorm);
        while ((rhsNorm > tolWeigh) && (obj.iter < obj.simParameters.itMaxNR) ...
            && (rhsNorm > obj.simParameters.absTol)) || obj.iter == 0
          obj.iter = obj.iter + 1;
          %
          % Solve system with increment
          du = linSyst.J\(-linSyst.rhs);
          %
          % Update tmpState
          obj.stateTmp.updateState_test(du,obj.dofManager);
          %
          % Compute residual and update Jacobian
          if isSinglePhaseFlow(obj.model)
            linSyst.computeFlowJacobian_Test(delta_t);
            mu = obj.material.getFluid().getDynViscosity();
            linSyst.computeFlowRHS_test(obj.statek,obj.stateTmp,delta_t,1/mu);
          end

        % ---> POROMECHANICS <--- (if theta method is not required)
          if  isPoromechanics(obj.model) && ~isSinglePhaseFlow(obj.model)
            linSyst.computePoroSyst_Test(obj.stateTmp,delta_t);
          end

        % ---> BIOT COUPLING <---
          if isPoromechanics(obj.model) && isSinglePhaseFlow(obj.model)
            linSyst.computeBiotSyst(delta_t,obj.statek,obj.stateTmp);
          end
          %
          linSyst.buildGlobalJacobianAndRhs();
          applyBCandForces_test(obj.model, obj.grid, obj.bound, obj.material, ...
            obj.t, linSyst, obj.stateTmp);
%           obj.bound.applyBCNeu(linSyst);
          %rhsNorm = norm(linSyst.rhs,obj.simParameters.pNorm);
          rhsNorm = computeRhsNorm(obj,linSyst);
          linSyst.resetJacobianAndRhs();
%           obj.bound.applyBCDir(linSyst);
%           obj.bound.applyBCNeu(linSyst);
%         rhsNorm = findNorm(obj,linSyst.rhs);
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
        % Manage next time step
          
          
%           flConv = true;
%           stateTmp.time = time;
%           %
%           % Print the solution if needed
%           if time > obj.tMax
%             printState(obj.printUtil,stateTmp);
%           else
%             printState(obj.printUtil,statek,stateTmp);
%           end
%           transferState(stateTmp,statek);
          %

%       end
%         end
%       end
      %
      %Update dt or perform backstep based on flConv
      % IMPLEMENT!
%       simStat = 0;
%       printState(obj.printUtil,statek);
      
%       tmpState = copy(resState);
%       obj.bound.iniBC(tmpState,obj.BCName);
%       du = obj.bound.applyDirDispl();
%       %
%       % Update tmpState
%       tmpState.updateState(du);
%       clear du
%       fprintf('-----------------------------------------------------------\n');
%       fprintf('Iter     ||rhs||\n');
%       %
%       % Apply Dir Bcond to the initial solution
%       % Assemble Jac and rhs
%       linSyst = Discretizer(obj.mesh,obj.elements,obj.material, ...
%                 obj.preP,obj.GaussPts);
%       if ModelType.isPoromechanics(obj.model) || ModelType.isCoupFlowPoro(obj.model)
%         linSyst.computeMat(tmpState);
%       end
%       if ModelType.isFlow(obj.model) || ModelType.isCoupFlowPoro(obj.model)
%         linSyst.computeFlowMat();
%       end
%       % Apply BC - to the Jacobian and RHS
%       obj.bound.applyBC(linSyst);
%       rhsNorm = norm(linSyst.rhs,obj.pNorm);
% %       rhsNorm = findNorm(obj,linSyst.rhs);
%       tolWeigh = obj.tol*rhsNorm;
%       iter = 0;
%       conv = false;
%       %
%       fprintf('0     %e\n',rhsNorm);
%       while ((rhsNorm > tolWeigh) && (iter < obj.itMax))
%         iter = iter + 1;
%         %
%         % Solve system with increment
%         du = linSyst.K\(-linSyst.rhs);
%         %
%         % Update tmpState
%         tmpState.updateState(du);
%         %
%         % Compute residual and updated Jacobian
%         linSyst.computeMat(tmpState);
%         obj.bound.applyBC(linSyst);
%         rhsNorm = norm(linSyst.rhs,obj.pNorm);
% %         rhsNorm = findNorm(obj,linSyst.rhs);
%         fprintf('%d     %e\n',iter,rhsNorm);
%         % 
%       end
%       % Check for convergence
%       if rhsNorm < tolWeigh
%         finalizeState(tmpState);
%         conv = true;
%         resState.stress = tmpState.stress;
%         resState.avStress = tmpState.avStress;
%         resState.displ = tmpState.displ;
%         resState.avStrain = tmpState.avStrain;
%       end
%     end
%   end
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
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %       flErr = false;
% % %       if isnumeric(obj.pNorm)
% % %         if obj.pNorm ~= 2
% % %           flErr = true;
% % %         end
% % %       elseif ischar(obj.pNorm)
% % %         obj.pNorm = lower(obj.pNorm);
% % %         if ~strcmp(obj.pNorm,'inf')
% % %           flErr = true;
% % %         end
% % %       end
% % %       if flErr
% % %         error('Accepted p-norms are: 2 (numeric variable) and "Inf" (char varibale)');
% % %       end
% % %       obj.pNorm = pNorm;
    end
%     function rhsNorm = findNorm(obj,f)
%       l = length(obj.BCName);
%       dofDir = [];
%       for i=1:l
%         cond = getBC(obj.bound,obj.BCName(i));
%         if strcmp(cond.boundType,'dir')  % Apply Dirichlet conditions
%           dofDir = [dofDir; cond.boundDof];
%         end
%       end
%       rhsNorm = norm(f(setxor(1:length(f),dofDir)),obj.pNorm);
%     end
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

    function rhsNorm = computeRhsNorm(obj,syst)
        %Return maximum norm of all Rhs block
        nRhs = length(syst.dofm.subList);
        rhs = zeros(nRhs,1);
        l1 = 0;
        for i = 1:nRhs
            ncomp = syst.dofm.numDof(i);
            rhs(i) = norm(syst.rhs(l1+1:l1+ncomp), obj.simParameters.pNorm);
            l1 = l1 + ncomp;
        end
        rhsNorm = max(rhs);
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