classdef MultidomainFCSolver < handle
  % Class for solving non linear problem involving multiple non conforming
  % domains

  properties (Access = protected)
    %
    nDom
    nInterf
    %
    t = 0
    tStep = 0
    iter
    dt
    nVars                 % total number of variable fields in the domains
  end


  properties (Access = public)
    simparams
    domains
    interfaces
  end


  methods (Access = public)
    
    function obj = MultidomainFCSolver(simparams,domains,interfaces)
      obj.setNonLinearSolver(simparams,domains,interfaces);
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

        obj.tStep = obj.tStep + 1;
        obj.t = obj.t + obj.dt;

        gresLog().log(0,'\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,obj.dt);
        gresLog().log(0,'-----------------------------------------------------------\n');
        gresLog().log(1,'Iter     ||rhs||     ||rhs||/||rhs_0||\n');

        for i = 1:obj.nDom
          obj.domains(i).applyDirVal(obj.t);
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
        rhsNorm = norm(cell2matrix(rhs),2);
        rhsNormIt0 = rhsNorm;

        tolWeigh = obj.simparams.relTol*rhsNorm;
        obj.iter = 0;

        flConv = false;
        %

        gresLog().log(1,'0     %e     %e\n',rhsNorm,rhsNorm/rhsNormIt0);

        while ~flConv && (obj.iter < obj.simparams.itMaxNR) || obj.iter == 0

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

        end % end newton loop


        % Check for convergence
        flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);

        manageNextTimeStep(obj,flConv);
      end
      %
    end

    function finalizeOutput(obj)
      % finalize print utils for domains and interfaces
      for i = 1:obj.nDom
        obj.domains(i).outstate.finalize();
      end

      for i = 1:obj.nInterf
        obj.interfaces{i}.outstate.finalize();
      end
    end

  end



  methods (Access = protected)
    function setNonLinearSolver(obj,simparams,dom,interf)

      % assumption: same set of simulation parameters for each domain
      obj.simparams = simparams;
      obj.domains = dom;
      obj.nDom = numel(dom);
      obj.interfaces = interf;
      obj.nInterf = numel(interf);

      obj.nVars = 0;
      for iD = 1:obj.nDom
        obj.domains(iD).stateOld = copy(obj.domains(iD).getState());
        obj.domains(iD).simparams = simparams;
        obj.nVars = obj.nVars + obj.domains(iD).dofm.getNumberOfVariables();
      end
    end



    function sol = solve(obj,J,rhs)

      % TO DO: clean this method
      J = cell2matrix(J);
      rhs = cell2matrix(rhs);
      %tic
      
      if size(J,1)>1e7
        fprintf('Solving linear system...\n')
        if norm(J-J','fro') < 1e-10
          % matrix is practically symmetric
          J = 0.5 * (J + J');
          %           J = J + speye(size(J)) * 1e-10;  % Regularize diagonal
          %           opts.type = 'ict';        % incomplete Cholesky with threshold
          %           opts.droptol = 1e-3;       % drop tolerance
          %           opts.diagcomp = 1e-3;      % diagonal compensation
          L = ichol(J);
          [sol,fl2,rr2,it2,rv2] = pcg(J,-rhs,1e-9,300,L,L');
        else
          p = symamd(J);
          sol = J(p,p) \ -rhs(p);
          sol(p) = sol;
        end
        fprintf('Linear system solved in %.4f s \n',toc)
      else
        %direct solver
        sol = J\(-rhs);
      end
      %       clear Jmat
      %       Jstab =

    end



    function out = computeRhsNorm(obj)
      %Return maximum norm of the entire domain
      rhsNorm = zeros(obj.nDom,1);
      for i = 1:obj.nDom
        nRhs = length(obj.domains(i).dofm.subList);
        rhsNorm_loc = zeros(nRhs,1);
        for j = 1:nRhs
          rhsNorm_loc(j) = norm(obj.domains(i).rhs{j}, obj.simparams.pNorm);
        end
        rhsNorm(i) = sqrt(sum(rhsNorm_loc.^2));
      end
      out = norm(rhsNorm);
    end



    function manageNextTimeStep(obj,flConv)

      if flConv % Convergence

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

      else

        % backstep
        for i = 1:obj.nDom
          goBackState(obj.domains(i));
        end

        for i = 1:obj.nInterf
          goBackState(obj.interfaces{i});
        end

        obj.t = obj.t - obj.dt;
        obj.tStep = obj.tStep - 1;

        obj.dt = obj.dt/obj.simparams.divFac;

        if obj.dt < obj.simparams.dtMin
          error('Minimum time step reached')
        else
          gresLog().log(0,'\n %s \n','BACKSTEP')
        end
      end
    end



    function J = assembleJacobian(obj)

      J = cell(obj.nVars + obj.nInterf);

      k = 0;

      for iD = 1:obj.nDom

        dom = obj.domains(iD);
        nV = dom.dofm.getNumberOfVariables;

        % inner domain blocks
        J(k+1:k+nV,k+1:k+nV) = getJacobian(dom);

        for iI = 1:numel(dom.interfaceList)

          q = dom.interfaceList(iI);

          % domain coupling blocks
          [J(k+1:k+nV,obj.nVars+q), J(obj.nVars+q,k+1:k+nV)] = ...
            getInterfaceJacobian(dom,iI);

          k = k + nV;

        end

      end

      for iI = 1:obj.nInterf
        interf = obj.interfaces{iI};
        % constraint blocks
        J{obj.nVars+iI,obj.nVars+iI} = getJacobian(interf);
      end

    end


    function rhs = assembleRhs(obj)
      % assemble blocks of rhs for multidomain system

      % each variable field of each domain represents a cell row
      rhs = cell(obj.nVars + obj.nInterf,1);
      
      k = 0;

      for iD = 1:obj.nDom
        nV = obj.domains(iD).dofm.getNumberOfVariables;
        rhs(k+1:k+nV) = getRhs(obj.domains(iD));
        k = k+nV;
      end

      for iI = 1:obj.nInterf

        rhs{k+1} = getRhs(obj.interfaces{iI});

        % each interface has one only multiplier field!
        k = k+1;
      end
    end



    function applyDirVal(obj)
      for i = 1:obj.nDom
        discretizer = obj.domains(i);

        % Check if boundary conditions are defined for the i-th domain
        if ~isempty(obj.domains(i).bcs)

          % Apply Dirichlet boundary values to i-th domain
            applyDirVal(discretizer,obj.t);
        end
      end
    end


    function computeMatricesAndRhs(obj)

      % Compute domain matrices
      for i = 1:obj.nDom
        discretizer = obj.domains(i);
        for j = 1:numel(discretizer.solverNames)
          s = discretizer.solverNames{j};
          computeMat(discretizer.solver(s), obj.state(i).prev, obj.dt);
        end
      end

      % Compute domain coupling matrices 
      for j = 1:obj.nInterf
        computeMat(obj.interfaces{j}, obj.dt);
      end

      % Compute domain rhs
      for i = 1:obj.nDom
        discretizer = obj.domains(i);
        for j = 1:numel(discretizer.solverNames)
          s = discretizer.solverNames{j};
          computeRhs(discretizer.solver(s), obj.state(i).prev, obj.dt);
        end
      end

      % Compute domain coupling matrices and rhs
      for j = 1:obj.nInterf
        computeRhs(obj.interfaces{j});
      end
      % Compute matrices and residuals for individual models of
      % the i-th domain
      %         computeMatricesAndRhs(...
      %           discretizer, obj.state(i).prev, obj.dt);
    end



    function applyBC(obj)
      for i = 1:obj.nDom
        discretizer = obj.domains(i);
        % Apply BCs to the blocks of the linear system
        applyBC(discretizer, obj.t);

        % Apply BC to coupling matrices
        for j = discretizer.interfaceList
          applyBC(obj.interfaces{j},i,discretizer.bcs,obj.t);
        end
      end
    end


    function updateState(obj,dSol)
      % update domain and interface state using incremental solution
      dSol_fix = dSol;
      for i = 1:obj.nDom
        N = obj.domains(i).dofm.totDoF;
        du = dSol(1:N);
        updateState(obj.domains(i),du);
        dSol = dSol(N+1:end);
      end

      % update interface state
      for j = 1:obj.nInterf
        N = obj.interfaces{j}.totMult;
        if N == 0
          du = dSol_fix;
        else
          du = dSol(1:N);
        end
        obj.interfaces{j}.updateState(du);
        dSol = dSol(N+1:end);
      end
    end


    function printState(obj)
      if obj.t > obj.simparams.tMax
        for i = 1:obj.nDom
          printState(obj.domains(i));
        end
        for i = 1:obj.nInterf
          id = obj.interfaces{i}.idSlave;
          tOld = obj.state(id).prev.t;
          printState(obj.interfaces{i},tOld)
        end
      else
        for i = 1:obj.nDom
          printState(obj.domains(i),obj.state(i).prev);
        end
        for i = 1:obj.nInterf
          id = obj.interfaces{i}.idDomain(2);
          tOld = obj.state(id).prev.t;
          tNew = obj.state(id).curr.t;
          printState(obj.interfaces{i},tOld,tNew)
        end
      end
    end



    function goOnState(obj)
      % transfer current state into previous state
      for i = 1:obj.nDom
        obj.state(i).prev = copy(obj.state(i).curr);
      end

      for i = 1:obj.nInterf
        obj.interfaces{i}.goOnState();
      end
    end

    function goBackState(obj)
      % transfer previous state into current state (backstep)
      for i = 1:obj.nDom
        obj.state(i).prev = copy(obj.state(i).curr);
      end

      for i = 1:obj.nInterf
        obj.interfaces{i}.goBackState();
      end
    end
  end
end
