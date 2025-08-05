classdef MeshGlueDual < MeshGlue
  % subclass of Mortar implementing mesh tying between different domains
  % using dual multipliers

  properties (Access = public)
    E
    Jcoupling
    f2
    Jinterf                     % part of the jacobian of the slave side without bcs applied
  end

  methods

    function obj = MeshGlueDual(id,inputStruct,domains)
      obj@MeshGlue(id,inputStruct,domains);
      assert(isfield(inputStruct,"Physics"), ...
        'Missing Physics field for interface %i',id);
      assert(strcmp(obj.multiplierType,'dual'), 'Multiplier type must be dual');
      %assert(strcmp(obj.physics,'Poisson'),'Static condensation is available only with Poisson problem')
      obj.totMult = 0;        % multipliers are condensed 
    end
    %
  end

  methods (Access=public)
    %
    function varargout = getJacobian(obj,varargin)
      % get jacobian blocks associated to specific field and specific
      % domain
      % if nargout = 2 -> get master/slave pair of jacobian blocks
      % if nargout = 1 -> return multiplier jacobian
      varargout = cell(1,nargout);
      for i = 1:nargout
        varargout{i} = [];
      end
    end

    function d = getDiagSlave(obj)
      d = sum(obj.D,2);
      % d=0 for nodes belonging to slave element that are not really in
      % contact with the master side
      d(d==0) = 1;
    end

    function rhs = getRhs(obj,fldId,varargin)
      % return rhs block associated to master/slave field
      switch nargin
        case 2
          % multiplier rhs block
          rhs = [];
        case 3
          domId = varargin{1};
          isDom = ismember(obj.idDomain,domId);
          if all(isDom)
            rhs = obj.rhsMaster{fldId} + obj.rhsSlave{fldId};
          elseif isDom(1)
            rhs = obj.rhsMaster{fldId};
          elseif isDom(2)
            rhs = obj.rhsSlave{fldId};
          else
            error('Input domain %i is not a valid master/slave',domId)
          end
      end
    end


    function computeMat(obj,~)
      solvSlave = obj.solvers(2).getSolver(obj.physics);
      if ~obj.isMatrixComputed()
        % mesh glue matrices are constant troughout the simulation
        % compute only once!
        computeMortarMatrices(obj);
      end
      obj.Jmaster{1} = obj.M;
      obj.Jslave{1} = obj.D;
      obj.Jinterf = solvSlave.J;
      obj.E = computeMortarOperator(obj);
      getCouplingMat(obj);
    end

    function getCouplingMat(obj)
      % return the coupling matrix for static condensation of dual
      % multipliers
      if ~strcmp(obj.multiplierType,'dual')
        return
      end

%       solvSlave = obj.solvers(2).getSolver(obj.physics);

      % coupling block between one domain and the other
      if ~isempty(obj.Jinterf)
        obj.Jcoupling = obj.Jinterf*obj.E;

        % update master block with condensation term
        solvMaster = obj.solvers(1).getSolver(obj.physics);
        solvMaster.J = solvMaster.J + obj.E'*(obj.Jinterf*obj.E);
      end
    end

    function computeRhs(obj)
      % compute rhs contributions for a specified input field
      for i = 1:obj.nFld
        % reset rhs multiplier
        obj.rhsMult{i} = zeros(getNumbMultipliers(obj),1);
        computeRhsMaster(obj,i);
        computeRhsSlave(obj,i);
      end
    end

    function applyBCmaster(obj,bound,bc,t)
      physic = bound.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(1).getSolver(physic),bound,bc,t);
      i = strcmp(obj.physics,physic);
      obj.rhsMaster{i}(bcEnts) = 0;
      obj.Jcoupling(:,bcEnts) = 0;
    end

    function applyBCslave(obj,bound,bc,t)
      physic = bound.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(2).getSolver(physic),bound,bc,t);
      bcEnts = removeSlaveBCdofs(obj,physic,bcEnts);
      i = strcmp(obj.physics,physic);
      obj.rhsSlave{i}(bcEnts) = 0;
      obj.Jcoupling(bcEnts,:) = 0;
      % remove interface slave dofs from matrix system (force zero) 
      dofSlave = getInterfSlaveDoF(obj);
      obj.Jcoupling(dofSlave,:) = 0;
      solvSlave =  obj.solvers(2).getSolver(obj.physics);
      solvSlave.applyDirBC([],dofSlave);
      % zero out rhs of slave interface dofs
      obj.rhsSlave{i}(dofSlave) = 0;
      solvSlave.rhs(dofSlave) = 0;
    end

    function updateState(obj,du)
      % get interface slave dof using mortar operator
      u_master = getState(obj.solvers(1).getSolver(obj.physics));
      u_slave = obj.E*u_master;
      solvSlave = obj.solvers(2).getSolver(obj.physics);
      dofSlave = getInterfSlaveDoF(obj);
      solvSlave.setState(dofSlave,u_slave(dofSlave));
      % reupdate slave domain with interface slave dofs
      solvSlave.updateState();
      % update multipliers
      var = getState(solvSlave);
      D = getDiagSlave(obj);
      obj.multipliers(1).curr = (1./D).*(obj.f2(dofSlave)-obj.Jinterf(dofSlave,:)*var);
    end
  end

  methods (Access = private)


    function computeRhsMaster(obj,i)
      solvSlave = obj.solvers(2).getSolver(obj.physics);
      dofSlave = getInterfSlaveDoF(obj);
      var = getState(solvSlave);
      ents = obj.dofm(2).getActiveEnts(obj.physics(i));
      var = var(ents);
      var(dofSlave) = 0;       % interface slave dof removal
      obj.rhsMaster{i} = obj.Jcoupling'*var;
    end

    function computeRhsSlave(obj,i)
      var = getState(obj.solvers(1).getSolver(obj.physics(i)));
      ents = obj.dofm(1).getActiveEnts(obj.physics(i));
      obj.rhsSlave{i} = obj.Jcoupling*var(ents);
      % update rhs of master side with condensation contribution
      solvSlave = obj.solvers(2).getSolver(obj.physics);
      entsSlave = obj.dofm(2).getActiveEnts(obj.physics(i));
      varSlave = getState(solvSlave);
      obj.f2 = solvSlave.rhs - solvSlave.J*varSlave(entsSlave); % this is just the forcing term
      %obj.f2 = solvSlave.rhs; % - solvSlave.J*varSlave(entsSlave); % this is just the forcing term
      obj.rhsMaster{i} =  obj.rhsMaster{i} + obj.E'*obj.f2;  % ... + E'*f_Gamma2
      % remove slave rhs contribution of interface slave dofs
      dofInter = getInterfSlaveDoF(obj);
      v = varSlave;
      varSlave(dofInter) = 2*varSlave(dofInter);
      varSlave = varSlave - v;          % keep only interface slave rhs
      % remove rhs contribution of interface dof that will be removed
      solvSlave.rhs = solvSlave.rhs - solvSlave.J*varSlave; 
    end

    function dof = getInterfSlaveDoF(obj)
      nodeInterf = obj.mesh.local2glob{2};
      dof = obj.dofm(2).getLocalDoF(nodeInterf,obj.physics);
    end
  end

  methods (Access = public)

    function out = isMatrixComputed(obj)
       out = all([~isempty(obj.D) ~isempty(obj.M)]);
    end

  end
end

