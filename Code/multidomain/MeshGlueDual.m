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
      assert(strcmp(obj.physic,'Poisson'),'Static condensation is available only with Poisson problem')
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

    function rhs = getRhs(obj,fld,varargin)
      % return rhs block associated to master/slave field

      if ~strcmp(fld,obj.physic)
        rhs = [];
        return
      end

      switch nargin
        case 2
          % multiplier rhs block
          rhs = [];
        case 3
          domId = varargin{1};
          isDom = ismember(obj.idDomain,domId);
          if all(isDom)
            rhs = obj.rhsMaster + obj.rhsSlave;
          elseif isDom(1)
            rhs = obj.rhsMaster;
          elseif isDom(2)
            rhs = obj.rhsSlave;
          else
            error('Input domain %i is not a valid master/slave',domId)
          end
      end
    end


    function computeMat(obj,~)
      solvSlave = obj.solvers(2).getSolver(obj.physic);
      if ~obj.isMatrixComputed()
        % mesh glue matrices are constant troughout the simulation
        % compute only once!
        computeMortarMatrices(obj);
      end
      obj.Jmaster = obj.M;
      obj.Jslave = obj.D;
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

%       solvSlave = obj.solvers(2).getSolver(obj.physic);

      % coupling block between one domain and the other
      obj.Jcoupling = obj.Jinterf*obj.E;

      % update master block with condensation term
      solvMaster = obj.solvers(1).getSolver(obj.physic);
      solvMaster.J = solvMaster.J + obj.E'*(obj.Jinterf*obj.E);
    end

    function computeRhs(obj)
      % compute rhs contributions for a specified input field
      % reset rhs multiplier
      obj.rhsMult = zeros(getNumbMultipliers(obj),1);
      computeRhsMaster(obj);
      computeRhsSlave(obj);
    end

    function applyBCmaster(obj,bc,t)
      physic = obj.solvers(1).bcs.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(1).getSolver(physic),bc,t);
      obj.rhsMaster(bcEnts) = 0;
      obj.Jcoupling(:,bcEnts) = 0;
    end

    function applyBCslave(obj,bc,t)
      physic = obj.solvers(1).bcs.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(2).getSolver(physic),bc,t);
      obj.rhsSlave(bcEnts) = 0;
      obj.Jcoupling(bcEnts,:) = 0;
      % remove interface slave dofs from matrix system (force zero) 
      dofSlave = getInterfSlaveDoF(obj);
      obj.Jcoupling(dofSlave,:) = 0;
      solvSlave =  obj.solvers(2).getSolver(obj.physic);
      solvSlave.applyDirBC([],dofSlave);
      % zero out rhs of slave interface dofs
      obj.rhsSlave(dofSlave) = 0;
      solvSlave.rhs(dofSlave) = 0;
    end

    function d = getDiagSlave(obj)
      d = sum(obj.D,2);
      % d=0 for nodes belonging to slave element that are not really in
      % contact with the master side
      d(d==0) = 1;
    end

    function updateState(obj,du)
      % get interface slave dof using mortar operator
      u_master = getState(obj.solvers(1).getSolver(obj.physic));
      u_slave = obj.E*u_master;
      solvSlave = obj.solvers(2).getSolver(obj.physic);
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


    function computeRhsMaster(obj)
      solvSlave = obj.solvers(2).getSolver(obj.physic);
      dofSlave = getInterfSlaveDoF(obj);
      var = getState(solvSlave);
      ents = obj.dofm(2).getActiveEnts(obj.physic);
      var = var(ents);
      var(dofSlave) = 0;       % interface slave dof removal
      obj.rhsMaster = obj.Jcoupling'*var;
    end

    function computeRhsSlave(obj)
      var = getState(obj.solvers(1).getSolver(obj.physic));
      ents = obj.dofm(1).getActiveEnts(obj.physic);
      obj.rhsSlave = obj.Jcoupling*var(ents);
      % update rhs of master side with condensation contribution
      solvSlave = obj.solvers(2).getSolver(obj.physic);
      entsSlave = obj.dofm(2).getActiveEnts(obj.physic);
      varSlave = getState(solvSlave);
      obj.f2 = solvSlave.rhs - solvSlave.J*varSlave(entsSlave); % this is just the forcing term
      %obj.f2 = solvSlave.rhs; % - solvSlave.J*varSlave(entsSlave); % this is just the forcing term
      obj.rhsMaster =  obj.rhsMaster + obj.E'*obj.f2;  % ... + E'*f_Gamma2
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
      dof = obj.dofm(2).getLocalDoF(nodeInterf,obj.physic);
    end
  end

  methods (Access = public)

    function out = isMatrixComputed(obj)
       out = all([~isempty(obj.D) ~isempty(obj.M)]);
    end

  end
end

