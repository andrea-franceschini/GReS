classdef MeshGlue < Mortar

  % Subclass of Mortar implementing mesh tying between different domains
  % mesh tying enforces continuity between primary variables
  % Current limitation: each interface ties a variable field separately
  % To tie different variable fields, create different interfaces

  properties (Access = public)
    Jmaster               % Jacobian block associated with master
    Jslave                % Jacobian block associated with slave
    Jmult                 % Jacobian block associated with multipliers
    rhsMaster             % master rhs
    rhsSlave              % slave rhs
    rhsMult               % multiplier rhs
    multipliers
    iniMultipliers
  end

  methods

    function obj = MeshGlue(id,inputStruct,domains)
      obj@Mortar(inputStruct,domains);

      obj.id = id;
      assert(isfield(inputStruct,"Physics"), ...
        'Missing Physics field for interface %i',id);
      obj.physic = split(inputStruct.Physics);
      assert(isfield(inputStruct,"Multiplier"),['XML:erorr' ...
        ''])

      obj.multiplierType = inputStruct.Multiplier.typeAttribute;
      isFldMaster = isField(obj.dofm(1),obj.physic);
      isFldSlave = isField(obj.dofm(2),obj.physic);
      assert(isFldMaster,['MeshGlue physic not available for ' ...
        'master domain %i'],obj.idDomain(1));
      assert(isFldSlave,['MeshGlue physic not available for ' ...
        'slave domain %i'],obj.idDomain(2));

      %computing mortar matrices D and M and updating list of slave entitities
      computeMortarMatrices(obj);

      %remove inacive multipliers from D and M
      id = find(~any(obj.D,2));
      obj.D(id,:) = [];
      obj.M(id,:) = [];

      % initialize Jacobian and rhs for the interface (now that active
      % multipliers are known)
      initializeInterface(obj)

      finalizeInterface(obj.mesh,obj.solvers);
      
      % use the updated mesh object for the output of the interface
      if ~isempty(obj.outStruct)
        obj.outStruct.VTK = VTKOutput(obj.mesh.msh(2),obj.outStruct.name);
      end
    end
    %
  end

  methods (Access=public)
    %
    function varargout = getJacobian(obj,field,domId)
      % get jacobian blocks associated to specific field and specific
      % domain
      % if nargout = 2 -> get master/slave pair of jacobian blocks
      % if nargout = 1 -> return multiplier jacobian
      
      if ~strcmp(field,obj.physic)
        % input field is not a mortar field
        varargout = cell(1,nargout);
        return
      end

      switch nargout
        case 1
          varargout{1} = obj.Jmult;
          % 0.5 is needed because the multiplier matrix is retrieved twice
        case 2
          % assign master/slave mortar matrix
          isDom = ismember(obj.idDomain,domId);
          if all(isDom)
            varargout{1} = obj.Jmaster' + obj.Jslave';
            varargout{2} = obj.Jmaster + obj.Jslave;
          elseif isDom(1)
            varargout{1} = obj.Jmaster';
            varargout{2} = obj.Jmaster;
          elseif isDom(2)
            varargout{1} = obj.Jslave';
            varargout{2} = obj.Jslave;
          else
            error('Input domain %i is not a valid master/slave',domId)
          end
      end
    end

    function rhs = getRhs(obj,fld,varargin)
      % return rhs block associated to master/slave/multiplier field
      
      if ~strcmp(fld,obj.physic)
        rhs = [];
        return
      end

      switch nargin
        case 2
          % multiplier rhs block
          rhs = obj.rhsMult;
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
      % return matrices for master and slave side in appropriate field
      if ~obj.isMatrixComputed()
        % mesh glue matrices are constant troughout the simulation
        % compute only once!
         computeMortarMatrices(obj);
      end
      obj.Jmaster = obj.M;
      obj.Jslave = obj.D;
    end

    function computeRhs(obj)
      % compute rhs contributions
      % reset rhs multiplier
      nMult = getNumbMultipliers(obj);
      obj.rhsMult = zeros(nMult,1);
      computeRhsMaster(obj);
      computeRhsSlave(obj);
    end

    function applyBCmaster(obj,bc,t)
      ph = obj.solvers(1).bcs.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(1).getSolver(ph),bc,t);
      obj.rhsMaster(bcEnts) = 0;
      obj.Jmaster(:,bcEnts) = 0;
    end

    function applyBCslave(obj,bc,t)
      ph = obj.solvers(2).bcs.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(2).getSolver(ph),bc,t);
      bcEnts = removeSlaveBCdofs(obj,ph,bcEnts);
      obj.rhsSlave(bcEnts) = 0;
      obj.Jslave(:,bcEnts) = 0;
    end

    function updateState(obj,du)
      actMult = getMultiplierDoF(obj);
      n = numel(actMult);
      obj.multipliers.curr(actMult) = obj.multipliers.curr(actMult) + du(1:n);
      %du = du(n+1:end);
    end
  end

  methods (Access = private)

    function initializeInterface(obj)
      obj.multipliers = struct('prev',[],'curr',[]);
      obj.totMult = 0;

      nDofMaster = getNumDoF(obj.dofm(1),obj.physic);
      nDofSlave = getNumDoF(obj.dofm(2),obj.physic);
      nDofMult = getNumbMultipliers(obj);
      obj.rhsMaster = zeros(nDofMaster,1);
      obj.rhsSlave = zeros(nDofSlave,1);
      obj.rhsMult = zeros(nDofMult,1);
      obj.multipliers.curr = zeros(nDofMult,1);
      obj.multipliers.prev = obj.multipliers.curr;
      obj.iniMultipliers = obj.multipliers.curr;
      obj.totMult = obj.totMult + nDofMult;
    end

    function computeRhsMaster(obj)
      actMult = getMultiplierDoF(obj);
      obj.rhsMaster = ...
        obj.Jmaster'*(obj.multipliers.curr(actMult)-obj.iniMultipliers(actMult));
      var = getState(obj.solvers(1).getSolver(obj.physic));
      ents = obj.dofm(1).getActiveEnts(obj.physic);
      obj.rhsMult = obj.rhsMult + obj.Jmaster*var(ents);
    end

    function computeRhsSlave(obj)
      actMult = getMultiplierDoF(obj);
      multIncrement = obj.multipliers.curr(actMult)-obj.iniMultipliers(actMult);
      obj.rhsSlave = obj.Jslave'*multIncrement;
      var = getState(obj.solvers(2).getSolver(obj.physic));
      ents = obj.dofm(2).getActiveEnts(obj.physic);
      obj.rhsMult = obj.rhsMult + obj.Jslave*var(ents);
      if ~isempty(obj.Jmult)
        obj.rhsMult = obj.rhsMult + ...
          obj.Jmult*multIncrement;
      end
    end
  end

  methods (Access = public)


    function [cellStr,pointStr] = buildPrintStruct(obj,fac)

      nCellData = obj.dofm(1).getDoFperEnt(obj.physic);
      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      if nargin == 1
        outVar = obj.multipliers.prev;
      elseif nargin == 2
        outVar = fac*obj.multipliers.curr + ...
          (1-fac)*obj.multipliers.prev;
      else
        error('Invalid number of input arguments for function buildPrintStruct');
      end

      for i = 1:nCellData
        cellStr(i).name = char(strcat(obj.physic,'_',num2str(i)));
        cellStr(i).data = outVar(i:nCellData:end);
      end

      pointStr = [];        % when using P0 multipliers
    end

    function goOnState(obj)
      % update the value of the multipliers
      obj.multipliers.prev = obj.multipliers.curr;
    end

    function goBackState(obj)
        obj.multipliers.curr = obj.multipliers.prev;
    end

%     function setDoFcount(obj,ndof)
%       obj.dofCount = ndof(obj.idDomain);
%     end


    function out = isMatrixComputed(obj)
      out = all([~isempty(obj.D) ~isempty(obj.M)]);
    end

    function E = computeMortarOperator(obj)
      % provide mortar interpolation operator for nodal multipliers

      assert(strcmp(obj.multiplierType,'dual'),['Mortar operator is avalilable' ...
        'with dual multipliers only']);
      E = sparse(size(obj.Jslave,2),size(obj.Jmaster,2));
      % get global index of slave nodes carrying multipliers
      slaveNodes = obj.mesh.local2glob{2};
      slavedof = obj.dofm(2).getLocalDoF(slaveNodes,obj.physic);
      d = getDiagSlave(obj);
      E(slavedof,:) = -obj.Jmaster.*(1./d);
      % expand E
    end

  end
end

