classdef MultidomainFCSolver < handle
  % Class for solving non linear problem involving multiple non conforming
  % domains using the mortar method

  properties (Access = private)
    %
    simParameters
    nDom
    nInterf
    %
    t = 0
    tStep = 0
    iter
    dt
    statek
    stateTmp
    systSize % [num_blocks, num_field_domain, num_field_interface]
    nfldInt
    nfldDom
  end


  properties (Access = public)
    state
    domains
    interfaces
    nDof
    results
  end


  methods (Access = public)
    function obj = MultidomainFCSolver(simParam,models,interfaces)
      obj.setNonLinearSolver(simParam,models,interfaces);
    end



    function NonLinearLoop(obj)

      % Initialize the time step increment
      obj.dt = obj.simParameters.dtIni;
      delta_t = obj.dt; % dynamic time step

      %
      flConv = true; %convergence flag


      % Loop over time
      while obj.t < obj.simParameters.tMax
        % Update the simulation time and time step ID
        absTol = obj.simParameters.absTol;
        obj.tStep = obj.tStep + 1;
        %new time update to fit the outTime list

        [obj.t, delta_t] = obj.updateTime(flConv, delta_t);

        if obj.simParameters.verbosity > 0
          fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,delta_t);
          fprintf('-----------------------------------------------------------\n');
        end
        if obj.simParameters.verbosity > 1
          fprintf('Iter     ||rhs||\n');
        end

        applyDirVal(obj);
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

          obj.iter = obj.iter + 1;

          J = assembleJacobian(obj);

          du = solve(obj,J,rhs);

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

        if flConv % Convergence
          % Advance state of non linear models
          for i = 1:obj.nDom
            obj.state(i).curr.t = obj.t;
            if isPoromechanics(obj.domains(i).ModelType)
              obj.domains(i).Discretizer.getSolver('Poromechanics').advanceState();
            end
          end

          printState(obj);
        end
        %
        %updateResults(obj);
        % Manage next time step
        delta_t = manageNextTimeStep(obj,delta_t,flConv);
      end
      %
    end

    function finalizeOutput(obj)
      % finalize print utils for domains and interfaces
      for i =1:obj.nDom
        obj.domains(i).OutState.finalize();
      end

      for i = 1:obj.nInterf
        obj.interfaces{i}.finalizeOutput();
      end
    end
  end



  methods (Access = private)
    function setNonLinearSolver(obj,simParam,dom,interf)
      obj.simParameters = simParam;
      obj.domains = dom;
      obj.nDom = numel(dom);
      obj.interfaces = interf;
      obj.nInterf = numel(interf);
      obj.state = repmat(struct('prev',{},'curr',{}),obj.nDom,1);
      % initialize a state structure for each domains
      for i = 1:obj.nDom
        obj.state(i).curr = obj.domains(i).Discretizer.state;
        obj.state(i).prev =  copy(obj.state(i).curr);
      end
      getSystemSize(obj);
      getNumField(obj);
      setDoFcounter(obj);
    end

    function updateResults(obj)
      id = getSymSurf(obj.interfaces{1});
      mult = obj.interfaces{1}.multipliers(1).curr;
      mult = reshape(mult,3,[]);
      obj.results = [obj.results;...
        struct('sn',mult(1,id)','t_norm',sqrt(mult(2,id).^2+mult(3,id).^2))];
    end

    function [t, dt] = updateTime(obj,conv,dt)
      t = obj.simParameters.tMax;
      told = t;
      for i = 1:obj.nDom
        if obj.domains(i).OutState.modTime
          tmp = find(obj.t<obj.domains(i).outState.timeList(),1,'first');
          if ~conv
            t = min([obj.t + obj.dt, obj.t + dt, obj.domains(i).OutState.timeList(tmp)]);
          else
            t = min([obj.t + obj.dt, obj.domains(i).OutState.timeList(tmp)]);
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



    function sol = solve(obj,J,rhs)
      % solve unstabilized system
      J = FCSolver.cell2matJac(J);
      rhs = cell2mat(rhs);
      tic
      if size(J,1)>1e4
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
        nRhs = length(obj.domains(i).DoFManager.subList);
        rhsNorm_loc = zeros(nRhs,1);
        for j = 1:nRhs
          rhsNorm_loc(j) = norm(obj.domains(i).Discretizer.rhs{j}, obj.simParameters.pNorm);
        end
        rhsNorm(i) = sqrt(sum(rhsNorm_loc.^2));
      end
      out = norm(rhsNorm);
    end



    function [dt] = manageNextTimeStep(obj,dt,flConv)
      if ~flConv   % Perform backstep
        goBackState(obj);
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
        for i = 1:obj.nDom
          if isFlow(obj.domains(i).ModelType)
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
    end



    function getSystemSize(obj)
      % assemble blocks of jacobian matrix for multidomain system
      N = 0;
      Ndof = 0;
      for i = 1:obj.nDom
        nf = numel(obj.domains(i).DoFManager.getFieldList());
        N = N + nf;
        Ndof(1) = Ndof(1) + obj.domains(i).DoFManager.totDoF;
      end
      Nfld = N;
      for i = 1:obj.nInterf
        N = N + obj.interfaces{i}.nFld;
      end
      Ninterf = N - Nfld;
      obj.systSize = [N,Nfld,Ninterf];
      obj.nDof = Ndof; % number of primary variable dofs
    end



    function getNumField(obj)
      for i = 1:obj.nDom
        obj.nfldDom(i) = numel(obj.domains(i).Discretizer.fields);
      end
      for i =1:obj.nInterf
        obj.nfldInt(i) = obj.interfaces{i}.nFld;
      end
    end



    function setDoFcounter(obj)
      N = zeros(obj.nDom,1);
      for i = 1:obj.nDom-1
        N(i) = obj.domains(i).DoFManager.totDoF;
      end
      N = [0; cumsum(N)];
      for i = 1:obj.nInterf
        obj.interfaces{i}.setDoFcount(N);
      end
    end



    function J = assembleJacobian(obj)
      % assemble blocks of jacobian matrix for multidomain system
      %[N,Nf,Ni] = deal(obj.systSize(1),obj.systSize(2),obj.systSize(3));
      J = cell(obj.systSize(1));
      f = 0;
      % populate jacobian with inner domain blocks
      for iD = 1:obj.nDom
        discr = obj.domains(iD).Discretizer;
        J(f+1:f+obj.nfldDom(iD),f+1:f+obj.nfldDom(iD)) = ...
          discr.assembleJacobian();
        for iF = 1:obj.nfldDom(iD)
          for iI = discr.interfaceList
            pos = find(strcmp(obj.interfaces{iI}.physics,discr.fields(iF)));
            jj = obj.systSize(2)+obj.nfldInt(iI)-obj.nfldInt(1)+pos;
            [J{iF+f,jj},J{jj,iF+f}] = getJacobian(...
              obj.interfaces{iF},pos,iD);
            [J{jj,jj}] = getJacobian(...
              obj.interfaces{iF},pos,iD);
          end
        end
        f = f + obj.nfldDom(iD);  % update field counter
      end 

      % provisional assembly of static condensation coupling block
      for iI = 1:obj.nInterf
        %
        interf = obj.interfaces{iI};
        if isa(interf,'MeshGlueDual')
          id = interf.idDomain;
          if isempty(J{id(1),id(2)})
            J{id(1),id(2)} = interf.Jcoupling';
          else
            J{id(1),id(2)} =  J{id(1),id(2)} + interf.Jcoupling';
          end
          if isempty(J{id(2),id(1)})
            J{id(2),id(1)} = interf.Jcoupling;
          else
            J{id(2),id(1)} = J{id(2),id(1)} + interf.Jcoupling;
          end
        end
      end

    end


    function rhs = assembleRhs(obj)
      % assemble blocks of rhs for multidomain system
      rhs = cell(obj.systSize(1),1);
      f = 0;

      for iD = 1:obj.nDom
        discr = obj.domains(iD).Discretizer;
        rhs{f+1:f+obj.nfldDom(iD)} = discr.assembleRhs();
        for iF = 1:obj.nfldDom(iD)
          for iI = discr.interfaceList
            pos = find(strcmp(obj.interfaces{iI}.physics,discr.fields(iF)));
            rhs{f+iF} = rhs{f+iF}{:} + getRhs(...
              obj.interfaces{iF},pos,iD);
            iMult = obj.systSize(2)+sum(obj.nfldInt(2:iI))+pos;
            if isempty(rhs{iMult})
              % dont compute rhsMult twice: 1field -> 1 interface
              rhs{iMult} = getRhs(obj.interfaces{iF},pos);
            end
          end
          f = f + 1;
        end
      end
    end



    function applyDirVal(obj)
      for i = 1:obj.nDom
        discretizer = obj.domains(i).Discretizer;
        bc = obj.domains(i).BoundaryConditions;

        % Check if boundary conditions are defined for the i-th domain
        if ~isempty(obj.domains(i).BoundaryConditions)

          % Apply Dirichlet boundary values to i-th domain
            applyDirVal(discretizer, bc, obj.t,i);
        end
      end
    end


    function computeMatricesAndRhs(obj)

      % Compute domain matrices
      for i = 1:obj.nDom
        discretizer = obj.domains(i).Discretizer;
        for j = 1:discretizer.numSolvers
          computeMat(discretizer.solver(j), obj.state(i).prev, obj.dt);
        end
      end

      % Compute domain coupling matrices 
      for j = 1:obj.nInterf
        computeMat(obj.interfaces{j}, obj.dt);
      end

      % Compute domain rhs
      for i = 1:obj.nDom
        discretizer = obj.domains(i).Discretizer;
        for j = 1:discretizer.numSolvers
          computeRhs(discretizer.solver(j), obj.state(i).prev, obj.dt);
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
        discretizer = obj.domains(i).Discretizer;
        bc = obj.domains(i).BoundaryConditions;
        % Apply BCs to the blocks of the linear system
        applyBC(discretizer, bc, obj.t, i);

        % Apply BC to coupling matrices
        for j = discretizer.interfaceList
          applyBC(obj.interfaces{j},i,bc,obj.t);
        end
      end
    end




    function updateState(obj,dSol)
      % update domain and interface state using incremental solution
      dSol_fix = dSol;
      for i = 1:obj.nDom
        discretizer = obj.domains(i).Discretizer;
        N = obj.domains(i).DoFManager.totDoF;
        du = dSol(1:N);
        updateState(discretizer,du);
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
      if obj.t > obj.simParameters.tMax
        for i = 1:obj.nDom
          printState(obj.domains(i).OutState);
        end
        for i = 1:obj.nInterf
          id = obj.interfaces{i}.idSlave;
          tOld = obj.state(id).prev.t;
          printState(obj.interfaces{i},tOld)
        end
      else
        for i = 1:obj.nDom
          printState(obj.domains(i).OutState,obj.domains(i).Discretizer,obj.state(i).prev);
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
