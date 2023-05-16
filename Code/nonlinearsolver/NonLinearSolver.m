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
    grid
    mesh
    elements
%     faces
    material
    bound
    %BCName
    preP
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
    function obj = NonLinearSolver(symmod,simParam,grid,mat,pre,bc,prtUtil,stateIni,varargin)
      nIn = nargin;
      data = varargin;
      obj.setNonLinearSolver(nIn,symmod,simParam,grid,mat,pre,bc,prtUtil,stateIni,data);
    end

    function [simStat] = NonLinearLoop(obj)
      simStat = 1;
      if obj.preP.nE(2) > 0  % There is at least one hexahedron
        linSyst = Discretizer(obj.model,obj.grid,obj.material, ...
                  obj.preP,obj.GaussPts);
      else
        linSyst = Discretizer(obj.model,obj.grid,obj.material, ...
                obj.preP);
      end
      % Compute H and P matrices for flow simulations and the gravity
      % contribution
      if isSinglePhaseFlow(obj.model)
        linSyst.computeFlowMat();
        linSyst.computeFlowRHSGravContribute();
      end
%       t  = 0;
%       tStep = 0;
      % Initialize the time step increment
      obj.dt = obj.simParameters.dtIni;
      %
      % Setup the boundary conditions
%       obj.bound.iniBC(obj.BCName,obj.stateTmp);
      %
      %
      % Loop over time
      while obj.t < obj.simParameters.tMax
        % Update the simulation time and time step ID
        obj.tStep = obj.tStep + 1;
        obj.t = obj.t + obj.dt;
        % Apply the Dirichlet condition value to the solution vector
        % (changes not needed for Coupled model)
        applyDirVal(obj.bound, obj.t, obj.stateTmp);
        %
        fprintf('\nTSTEP %d   ---  TIME %f\n',obj.tStep,obj.t);
        fprintf('-----------------------------------------------------------\n');
        fprintf('Iter     ||rhs||\n');
        
        
        if isCoupled(obj.model)
%         %compute Jacobian and Residual for coupled system
          linSyst.computeCoupleSyst(obj.simParameters.theta,obj.dt,obj.statek,obj.stateTmp)
        elseif isPoromechanics(obj.model)
          % Compute Jacobian and residual of the poromechanical problem
          linSyst.computePoroSyst(obj.stateTmp);
        elseif isSinglePhaseFlow(obj.model)
          % Compute Jacobian and residual of the flow problem 
          linSyst.computeFlowSystMat(obj.simParameters.theta,obj.dt);
          %
          linSyst.computeFlowRHS(obj.statek,obj.stateTmp);
        end
       
        %
        % Apply Neu and Dir conditions
        applyBCandForces(obj.model, obj.grid, obj.bound, ...
          obj.t, linSyst,obj.simParameters.theta,obj.dt);
%         obj.bound.applyBCNeu(linSyst);
        rhsNorm = norm(linSyst.rhs,obj.simParameters.pNorm);
%         obj.bound.applyBCDir(linSyst);
%         obj.bound.applyBCNeu(linSyst);
        tolWeigh = obj.simParameters.relTol*rhsNorm;
        obj.iter = 0;
        %
        fprintf('0     %e\n',rhsNorm);
        while ((rhsNorm > tolWeigh) && (obj.iter < obj.simParameters.itMaxNR) && (rhsNorm > obj.simParameters.absTol)) || obj.iter == 0
          obj.iter = obj.iter + 1;
          %
          % Solve system with increment
          du = linSyst.K\(-linSyst.rhs);
          %
          % Update tmpState
          obj.stateTmp.updateState(du); 
          %
          % Compute residual and update Jacobian
          if isCoupled(obj.model)
              linSyst.computeCoupleSyst(obj.simParameters.theta,obj.dt,obj.statek,obj.stateTmp);
          elseif isPoromechanics(obj.model)
            linSyst.computePoroSyst(obj.stateTmp);
          %
          elseif isSinglePhaseFlow(obj.model)
            linSyst.computeFlowSystMat(obj.simParameters.theta,obj.dt);   % Not needed
%             since the problem is linear
            linSyst.computeFlowRHS(obj.statek,obj.stateTmp);
          end
          %
          applyBCandForces(obj.model, obj.grid, obj.bound, ...
            obj.t, linSyst,obj.simParameters.theta,obj.dt);
%           obj.bound.applyBCNeu(linSyst);
          rhsNorm = norm(linSyst.rhs,obj.simParameters.pNorm);
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
          if obj.t > obj.simParameters.tMax   % For Steady State
            printState(obj.printUtil,obj.stateTmp);
          else
            printState(obj.printUtil,obj.statek,obj.stateTmp);
          end
          %updating output matrixes
        end
        %
        % Manage next time step
        manageNextTimeStep(obj,flConv);
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
    function setNonLinearSolver(obj,nIn,symmod,simParam,grid,mat,pre,bc,prtUtil,stateIni,data)
      obj.model = symmod;
      obj.simParameters = simParam;
      obj.grid = grid;
      obj.mesh = grid.topology;
      obj.elements = grid.cells;
      obj.material = mat;
      obj.preP = pre;
      obj.bound = bc;
      %obj.BCName = BName;
      obj.printUtil = prtUtil;
      obj.statek = stateIni;
      obj.stateTmp = copy(stateIni);
      if nIn > 9
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
  
    function manageNextTimeStep(obj,flConv)
      if ~flConv   % Perform backstep
        transferState(obj.statek,obj.stateTmp);
        obj.t = obj.t - obj.dt;
        obj.tStep = obj.tStep - 1;
        obj.dt = obj.dt/obj.simParameters.divFac;  % Time increment chop
        assert(obj.dt >= obj.simParameters.dtMin,'Minimum time step reached');
        fprintf('\n %s \n','BACKSTEP');
      else % Go on if converged
        if isSinglePhaseFlow(obj.model)
          dpMax = max(abs(obj.stateTmp.pressure-obj.statek.pressure));
          obj.dt = min([obj.dt * min([obj.simParameters.multFac, ...
            (1+obj.simParameters.relaxFac)*obj.simParameters.pTarget/(dpMax + obj.simParameters.relaxFac*obj.simParameters.pTarget)]), ...
            obj.simParameters.dtMax]);
        elseif isPoromechanics(obj.model)
          obj.dt = min(obj.dt * obj.simParameters.multFac,obj.simParameters.dtMax);
        end
        transferState(obj.stateTmp,obj.statek);
        %
        if ((obj.t + obj.dt) > obj.simParameters.tMax)
          obj.dt = obj.simParameters.tMax - obj.t;
        end
      end
    end
  end
end