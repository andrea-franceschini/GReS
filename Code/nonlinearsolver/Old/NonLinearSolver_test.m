classdef NonLinearSolver < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = private)
    itMax = 10
    tol = 1.e-6
    pNorm = 2
    theta = 0.5
    dtIni = 0.5    % 0.5
    tMax = 100   %100
    %
    model
    mesh
    elements
    material
    bound
    BCName
    preP
    GaussPts
    printUtil
  end
  
  methods (Access = public)
    function obj = NonLinearSolver(symmod,msh,elem,mat,pre,bc,BName,prtUtil,varargin)
      nIn = nargin;
      data = varargin;
      obj.setNonLinearSolver(nIn,symmod,msh,elem,mat,pre,bc,BName,prtUtil,data);
    end

    function [simStat] = NonLinearLoop(obj,statek)
      simStat = 1;
      if obj.preP.nE(2) > 0
        linSyst = Discretizer(obj.model,obj.mesh,obj.elements,obj.material, ...
                  obj.preP,obj.GaussPts);
      else
        linSyst = Discretizer(obj.model,obj.mesh,obj.elements,obj.material, ...
                obj.preP);
      end
      % Compute H and P matrices for flow simulations
      if isSinglePhaseFlow(obj.model)
        linSyst.computeFlowMat();
        linSyst.computeFlowRHSGravContribute();
      end
      time  = 0;
      tStep = 0;
      dt = obj.dtIni;
      %
%       stateTmp = copy(statek);
%       obj.bound.iniBC(obj.BCName,stateTmp);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       rhsReg = zeros(250,1);
%       pNodReg = zeros(250,1);
%       pNodReg(1) = statek.pressure(5714);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       printState(obj.printUtil,statek);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      cond = getBC(obj.bound,obj.BCName(1));
%       statek.pressure(cond.boundDof) = cond.boundVal;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Loop over time
      H = linSyst.H;
      P = linSyst.P;
%       pOld = statek.pressure;
%       pOld = zeros(obj.mesh.nNodes,1);
      max_z = 410;
      gamma = obj.material.getMaterial(2*obj.preP.nMat+1).getFluidSpecWeight();
      pOld = gamma*(max_z-obj.mesh.coordinates(:,3));
      rhsConst = linSyst.rhsConst;
      V = VTKOutput(obj.mesh);
      pointData3D = repmat(struct('name', 1, 'data', 1), 1, 1);
      pointData3D(1).name = 'p';
      pointData3D(1).data = pOld;
%       cellData3D = repmat(struct('name', 1, 'data', 1), 0, 1);
      V.writeVTKFile(-1, pointData3D, [], [], []);
      %%%
      % SS sol
      %%%
      H_SS = H;
%       rhs = zeros(obj.mesh.nNodes,1);
      rhs = -rhsConst;
      H_SS(obj.mesh.nNodes*(cond.boundDof-1)+cond.boundDof) = 1.e12;
      rhs(cond.boundDof) = 1.e12*cond.boundVal;
      pSS = H_SS\rhs;
      T_lin = -1;
      l = 1:5:200;
      pointData3D = repmat(struct('name', 1, 'data', 1), 1, 1);
      pointData3D(1).name = 'p';
      pointData3D(1).data = pSS;
%       cellData3D = repmat(struct('name', 1, 'data', 1), 0, 1);
      V.writeVTKFile(0, pointData3D, [], [], []);
      while time < obj.tMax
        tStep = tStep + 1;
        time = time + dt;
%         obj.bound.applyDirVal();
        %
        % Update tmpState
%         stateTmp.updateState(du);
%         clear du
        fprintf('TSTEP %d   ---  TIME %f\n',tStep,time);
        fprintf('-----------------------------------------------------------\n');
        fprintf('Iter     ||rhs||\n');
        
%         if isSinglePhaseFlow(obj.model)
%           linSyst.computeFlowSystMat(obj.theta,dt);
%           %
%           linSyst.computeFlowRHS(statek,stateTmp);
%         end
        K1 = obj.theta*H + P/dt;
        K2 = P/dt - (1-obj.theta)*H;
%         K1 = H;
%         K2 = P;
%         rhs = K2*pOld;
%         rhs2 = H*pOld;
        rhs = K2*pOld - rhsConst;
        rhs2 = H*pOld - rhsConst;
%         rhs = zeros(obj.mesh.nNodes,1);
%         rhs2=rhs;
        %
        K1(obj.mesh.nNodes*(cond.boundDof-1)+cond.boundDof) = 1.e12;
        if time<T_lin
          rhs(cond.boundDof) = 1.e12*cond.boundVal*time/T_lin;
        else
          rhs(cond.boundDof) = 1.e12*cond.boundVal;
        end
%         K1(obj.mesh.nNodes*(cond.boundDof-1)+cond.boundDof) = 1;
%         rhs(cond.boundDof) = cond.boundVal;
%         for i=1:length(cond.boundDof)
%           [~,c] = find(K1(cond.boundDof(i),:));
%           for j=1:length(c)
%             if c(j) ~= cond.boundDof(i)
%               K1(cond.boundDof(i),c(j)) = 0;
%             end
%           end
%         end
        rhs_cp = rhs;
        rhs2_cp = rhs2;
        rhs_cp(cond.boundDof) = 0;
        rhs2_cp(cond.boundDof) = 0;
        fprintf('0     %e    %e\n',norm(rhs_cp),norm(rhs2_cp));
%         rhs2 = rhs;
%         rhs2(cond.boundDof) = 0;
%         rhsNorm2 = norm(rhs2,2);
        %
        % Apply Neu and Dir conditions
%         obj.bound.applyBC(linSyst);
%         obj.bound.applyBCNeu(linSyst);
%         rhsNorm = norm(linSyst.rhs,obj.pNorm);
%         obj.bound.applyBCDir(linSyst);
%         obj.bound.applyBCNeu(linSyst);
%         tolWeigh = obj.tol*rhsNorm;
%         iter = 0;
%         flConv = false;
        %
%         rhs2 = linSyst.H*statek.pressure;
%         cond = getBC(obj.bound,obj.BCName(1));
        
%         rhsReg(tStep) = rhsNorm;
%         while ((rhsNorm > tolWeigh) && (iter < obj.itMax) && (rhsNorm > 1.e-10))
%           iter = iter + 1;
          %
          % Solve system with increment
          pNew = K1\rhs;
          %%% Norm diff
          normDiff = 0;
          for i=1:obj.mesh.nNodes
            normDiff = normDiff + (pNew(i)-pSS(i))^2*obj.elements.volNod(i);
          end
          normDiff = sqrt(normDiff);
          fprintf('%s  %e\n\n','Norm delta sol to SS',normDiff);
          %
          % Update tmpState
%           stateTmp.updateState(du); 
          %
          % Compute residual and updated Jacobian
%           if isPoromechanics(obj.model)
%             linSyst.computePoroSyst(stateTmp);
%           end
          %
%           if isSinglePhaseFlow(obj.model)
%             linSyst.computeFlowSystMat(obj.theta,dt);
%             linSyst.computeFlowRHS(statek,stateTmp);
%           end
          %
%           obj.bound.applyBC(linSyst);
%           obj.bound.applyBCNeu(linSyst);
%           rhsNorm = norm(linSyst.rhs,obj.pNorm);
%           obj.bound.applyBCDir(linSyst);
%           obj.bound.applyBCNeu(linSyst);
%         rhsNorm = findNorm(obj,linSyst.rhs);
%           fprintf('%d     %e\n',iter,rhsNorm);
%         end
        % Check for convergence
%         if rhsNorm < tolWeigh || rhsNorm < 1.e-10 % Convergence
%           pNodReg(tStep+1) = stateTmp.pressure(5714);
%           fprintf('\n');
%           flConv = true;
%           stateTmp.time = time;
          %
          % Print the solution if needed
%           if time > obj.tMax
%             printState(obj.printUtil,stateTmp);
%           else
%             printState(obj.printUtil,statek,stateTmp);
%           end
%           transferState(stateTmp,statek);
%           if ModelType.isPoromechanics(obj.model)
% %             finalizeState(stateTmp);
%             statek.stress = stateTmp.stress;
% %             statek.avStress = stateTmp.avStress;
%             statek.displ = stateTmp.displ;
% %             statek.avStrain = stateTmp.avStrain;
%             statek.time = stateTmp.time;
%           end
%           %
%           if ModelType.isSinglePhaseFlow(obj.model)
%             statek.pressure = stateTmp.pressure;
%             statek.time = stateTmp.time;
%           end
          %
          pOld = pNew;
          if ismember(tStep,l)
            pointData3D = repmat(struct('name', 1, 'data', 1), 1, 1);
            pointData3D(1).name = 'p';
            pointData3D(1).data = pOld;
      %       cellData3D = repmat(struct('name', 1, 'data', 1), 0, 1);
            V.writeVTKFile(time, pointData3D, [], [], []);
          end
%           pointData3D = repmat(struct('name', 1, 'data', 1), 1, 1);
%           pointData3D(1).name = 'p';
%           pointData3D(1).data = pOld;
%           cellData3D = repmat(struct('name', 1, 'data', 1), 0, 1);
%           V.writeVTKFile(time, pointData3D, [], [], []);
%         pause(2)
      end
      if 5<1
      pointData3D = repmat(struct('name', 1, 'data', 1), 1, 1);
      pointData3D(1).name = 'p';
      pointData3D(1).data = pOld;
%       cellData3D = repmat(struct('name', 1, 'data', 1), 0, 1);
      V.writeVTKFile(time, pointData3D, [], [], []);
      end
%       end
      %
      %Update dt or perform backstep based on flConv
      % IMPLEMENT!
      simStat = 0;
%       printState(obj.printUtil,statek);
    end
      
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
  
  methods (Access = private)
    function setNonLinearSolver(obj,nIn,symmod,msh,elem,mat,pre,bc,BName,prtUtil,data)
      obj.model = symmod;
      obj.mesh = msh;
      obj.elements = elem;
      obj.material = mat;
      obj.preP = pre;
      obj.bound = bc;
      obj.BCName = BName;
      obj.printUtil = prtUtil;
      if nIn > 8
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
  end
end

