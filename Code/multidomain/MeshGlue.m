classdef MeshGlue < Mortar
  % subclass of Mortar implementing mesh tying between different domains
  % mesh tying enforces continuity between primary variables
  % different physics can be tyied at once

  properties (Access = private)
    dofCount
  end

  properties (Access = public)
    physics
    multipliers
    iniMultipliers
    totMult
  end

  properties (Access = protected)
    multEnts
  end

  methods

    function obj = MeshGlue(id,inputStruct,domains)
      obj@Mortar(inputStruct,domains);
      assert(isfield(inputStruct,"Physics"), ...
        'Missing Physics field for interface %i',id);
      obj.physics = split(inputStruct.Physics);
      obj.nFld = numel(obj.physics);
            assert(isfield(inputStruct,"Multiplier"),['XML:erorr' ...
              ''])
        obj.multiplierType = inputStruct.Multiplier.typeAttribute;
        actCells = getActiveCells(obj.mesh,2);
        switch obj.multiplierType
          case 'P0'
            obj.multEnts = actCells;
          otherwise
            obj.multEnts = (unique(obj.mesh.msh(2).surfaces(actCells,:)));
        end
        %initializeJacobianAndRhs(obj);
        isFldMaster = isField(obj.dofm(1),obj.physics);
        isFldSlave = isField(obj.dofm(2),obj.physics);
        assert(isFldMaster,['MeshGlue physic not available for ' ...
          'master domain %i'],obj.idDomain(1));
        assert(isFldSlave,['MeshGlue physic not available for ' ...
          'slave domain %i'],obj.idDomain(2));

        % computing mortar matrix and updating list of slave entitities
        computeMortarMatrices(obj);
        id = find(all(obj.D==0,2));
        ncomp = obj.dofm(2).getDoFperEnt(obj.physics);
        id2 = id(ncomp:ncomp:end)/ncomp;
        obj.multEnts(id2) = [];
        obj.D(id,:) = [];
        obj.M(id,:) = [];
        initializeJacobianAndRhs(obj)
    end

    
    %
  end

  methods (Access=public)
    %
    function varargout = getJacobian(obj,fldId,domId)
      % get jacobian blocks associated to specific field and specific
      % domain
      % if nargout = 2 -> get master/slave pair of jacobian blocks
      % if nargout = 1 -> return multiplier jacobian
      switch nargout
        case 1
          varargout{1} = obj.Jmult{fldId};
          % 0.5 is needed because the multiplier matrix is retrieved twice
        case 2
          % assign master/slave mortar matrix
          isDom = ismember(obj.idDomain,domId);
          if all(isDom)
            varargout{1} = (obj.Jmaster{fldId})' + (obj.Jslave{fldId})';
            varargout{2} = obj.Jmaster{fldId} + obj.Jslave{fldId};
          elseif isDom(1)
            varargout{1} = (obj.Jmaster{fldId})';
            varargout{2} = obj.Jmaster{fldId};
          elseif isDom(2)
            varargout{1} = (obj.Jslave{fldId})';
            varargout{2} = obj.Jslave{fldId};
          else
            error('Input domain %i is not a valid master/slave',domId)
          end
      end
    end

    function rhs = getRhs(obj,fldId,varargin)
      % return rhs block associated to master/slave/multiplier field
      switch nargin
        case 2
          % multiplier rhs block
          rhs = obj.rhsMult{fldId};
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
      % return matrices for master and slave side in appropriate field
      if ~obj.isMatrixComputed()
        % mesh glue matrices are constant troughout the simulation
        % compute only once!
         computeMortarMatrices(obj);
         % remove inactive multipliers from D and M
      end

      obj.Jmaster{1} = obj.M;
      obj.Jslave{1} = obj.D;
%      
      %       for i = 1:obj.nFld
      %         % map local mortar matrices to global indices
      %         obj.Jmaster{i} = obj.getMatrix(1,obj.physics(i));
      %         obj.Jslave{i} = obj.getMatrix(2,obj.physics(i));
      %       end
    end

    function computeRhs(obj)
      % compute rhs contributions for a specified input field
      for i = 1:obj.nFld
        % reset rhs multiplier
        [nMult,~] = getNumbMultipliers(obj);
        obj.rhsMult{i} = zeros(nMult,1);
        computeRhsMaster(obj,i);
        computeRhsSlave(obj,i);
      end
    end

    function applyBCmaster(obj,bound,bc,t)
      physic = bound.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(1).getSolver(physic),bound,bc,t);
      i = strcmp(obj.physics,physic);
      obj.rhsMaster{i}(bcEnts) = 0;
      obj.Jmaster{i}(:,bcEnts) = 0;
    end

    function applyBCslave(obj,bound,bc,t)
      physic = bound.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(2).getSolver(physic),bound,bc,t);
      bcEnts = removeSlaveBCdofs(obj,physic,bcEnts);
      i = strcmp(obj.physics,physic);
      obj.rhsSlave{i}(bcEnts) = 0;
      obj.Jslave{i}(:,bcEnts) = 0;
    end

    function updateState(obj,du)
      for i = 1:obj.nFld
        actMult = getMultiplierDoF(obj);
        n = numel(actMult);
        obj.multipliers(i).curr(actMult) = obj.multipliers(i).curr(actMult) + du(1:n);
        du = du(n+1:end);
      end
    end
  end

  methods (Access = private)

    function initializeJacobianAndRhs(obj)
      [obj.Jmaster,obj.Jslave,obj.Jmult] = deal(cell(obj.nFld,1));
      [obj.rhsMaster,obj.rhsSlave, obj.iniMultipliers] = deal(cell(obj.nFld,1));
      obj.multipliers = repmat(struct('prev',[],'curr',[]),obj.nFld,1);
      obj.totMult = 0;

      for i = 1:obj.nFld
        nDofMaster = getNumDoF(obj.dofm(1),obj.physics(i));
        nDofSlave = getNumDoF(obj.dofm(2),obj.physics(i));
        [nDofMultAct,nDofMult] = getNumbMultipliers(obj);
        obj.rhsMaster{i} = zeros(nDofMaster,1);
        obj.rhsSlave{i} = zeros(nDofSlave,1);
        obj.rhsMult{i} = zeros(nDofMultAct,1);
        obj.multipliers(i).curr = zeros(nDofMult,1);
        obj.multipliers(i).prev = obj.multipliers(i).curr;
        obj.iniMultipliers{i} = obj.multipliers(i).curr;
        obj.totMult = obj.totMult + nDofMultAct;
      end
    end

    function computeRhsMaster(obj,i)
      actMult = getMultiplierDoF(obj);
      obj.rhsMaster{i} = ...
        obj.Jmaster{i}'*(obj.multipliers(i).curr(actMult)-obj.iniMultipliers{i}(actMult));
      var = getState(obj.solvers(1).getSolver(obj.physics(i)));
      ents = obj.dofm(1).getActiveEnts(obj.physics(i));
      obj.rhsMult{i} = obj.rhsMult{i} + obj.Jmaster{i}*var(ents);
    end

    function computeRhsSlave(obj,i)
      actMult = getMultiplierDoF(obj);
      multIncrement = obj.multipliers(i).curr(actMult)-obj.iniMultipliers{i}(actMult);
      obj.rhsSlave{i} = ...
        obj.Jslave{i}'*multIncrement;
      var = getState(obj.solvers(2).getSolver(obj.physics(i)));
      ents = obj.dofm(2).getActiveEnts(obj.physics(i));
      obj.rhsMult{i} = obj.rhsMult{i} + obj.Jslave{i}*var(ents);
      if ~isempty(obj.Jmult{i})
        obj.rhsMult{i} = obj.rhsMult{i} + ...
          obj.Jmult{i}*multIncrement;
      end
    end
  end

  methods (Access = public)

    function [dofr,dofc,mat] = computeLocMaster(obj,imult,im,Nmult,Nmaster)
      mat = obj.quadrature.integrate(@(a,b) pagemtimes(a,'ctranspose',b,'none'),...
        Nmult,Nmaster);
      nodeMaster = obj.mesh.local2glob{1}(obj.mesh.msh(1).surfaces(im,:));
      fld = obj.dofm(1).getFieldId(obj.physics);
      dofc = obj.dofm(1).getLocalDoF(nodeMaster,fld);
      dofr = getMultiplierDoF(obj,imult);
    end

    function [dofr,dofc,mat] = computeLocSlave(obj,imult,is,mat)
      if strcmp(obj.multiplierType,'dual')
        % lump local D matrix
        mat = diag(sum(mat,2));
      end
      nodeSlave = obj.mesh.local2glob{2}(obj.mesh.msh(2).surfaces(is,:));
      fld = obj.dofm(2).getFieldId(obj.physics);
      dofc = obj.dofm(2).getLocalDoF(nodeSlave,fld);
      dofr = getMultiplierDoF(obj,imult);
    end


    function [cellStr,pointStr] = buildPrintStruct(obj,fldId,fac)

      fieldName = obj.physics(fldId);
      nCellData = obj.dofm(1).getDoFperEnt(fieldName);
      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      if nargin ==2
        outVar = obj.multipliers(fldId).prev;
      elseif nargin == 3
        outVar = fac*obj.multipliers(fldId).curr + ...
          (1-fac)*obj.multipliers(fldId).prev;
      else
        error('Invalid number of input arguments for function buildPrintStruct');
      end

      for i = 1:nCellData
        cellStr(i).name = char(strcat(fieldName,'_',num2str(i)));
        cellStr(i).data = outVar(i:nCellData:end);
      end

      pointStr = [];        % when using P0 multipliers
    end

    function goOnState(obj)
      % update the value of the multipliers
      for i = 1:obj.nFld
        obj.multipliers(i).prev = obj.multipliers(i).curr;
      end
    end

    function goBackState(obj)
      for i = 1:obj.nFld
        obj.multipliers(i).curr = obj.multipliers(i).curr;
      end
    end

    function setDoFcount(obj,ndof)
      obj.dofCount = ndof(obj.idDomain);
    end


    function out = isMatrixComputed(obj)
      out = all([~isempty(obj.D) ~isempty(obj.M)]);
    end

    function Nmult = computeMultiplierBasisF(obj,el,NslaveIn)
      elem = obj.getElem(2,el);
      switch obj.multiplierType
        case 'P0'
          Nmult = ones(size(NslaveIn,1),1);
        case 'standard'
          Nmult = NslaveIn;
        case 'dual'
          Ns = getBasisFinGPoints(elem);
          gpW = getDerBasisFAndDet(elem,el);
          M = Ns'*(Ns.*gpW');
          D = diag(Ns'*gpW');
          A = M\D;
          Nmult = NslaveIn*A;
      end
    end

    function [nDofAct,nDofTot] = getNumbMultipliers(obj)
      % nDoFAct -> number of currently active multipliers in the interface
      % nDoFTot -> total number of dof in the whole interface

      nc = obj.dofm(2).getDoFperEnt(obj.physics);
      switch obj.multiplierType
        case 'P0' 
          nDofAct = nc*numel(obj.multEnts);
          nDofTot = nc*obj.mesh.msh(2).nSurfaces;
        case {'dual','standard'} % nodal multipliers
          % get number of nodes in active cells
          nDofAct = nc*numel(obj.multEnts);
          nDofTot = nc*obj.mesh.msh(2).nNodes;
      end
    end

    function dofs = getMultiplierDoF(obj,id)
      % return global indices associated to active multipliers to access
      % the multiplier vector and temporarily assemble Jmaster and Jslave
      % id - local index of active cell in slave side
      nc = obj.dofm(2).getDoFperEnt(obj.physics);
      if nargin > 1
        assert(isscalar(id),'Input id must be a scalar integer \n')
        if strcmp(obj.multiplierType,'P0')
          dofs = dofId(id,nc);
        else
          is = obj.mesh.activeCells{2}(id);
          nodes = obj.mesh.msh(2).surfaces(is,:);
          dofs = dofId(nodes,nc);
        end
      else
        if strcmp(obj.multiplierType,'P0')
          dofs = dofId(getActiveCells(obj.mesh,2),nc);
        else
          dofs = dofId(obj.multEnts,nc);
        end
      end
    end

    function E = computeMortarOperator(obj)
      % provide mortar interpolation operator for nodal multipliers

      assert(strcmp(obj.multiplierType,'dual'),['Mortar operator is avalilable' ...
        'with dual multipliers only']);
      E = sparse(size(obj.Jslave{1},2),size(obj.Jmaster{1},2));
      % get global index of slave nodes carrying multipliers
      slaveNodes = obj.mesh.local2glob{2}(obj.multEnts);
      slavedof = obj.dofm(2).getLocalDoF(slaveNodes,obj.physics);
      d = getDiagSlave(obj);
      E(slavedof,:) = -obj.Jmaster{1}.*(1./d);
      % expand E
      

    end

    function [bcDofs,bcVals] = removeSlaveBCdofs(obj,bcPhysics,bcData,domId)
      % this method updates the list of bc dofs and values removing dofs on
      % the slave interface (only for nodal multipliers)
      % this avoid overconstraint and solvability issues

      % bcData: nDofBC x 2 matrix.
      % first column -> dofs
      % second columns -> vals (optional)

      bcDofs = bcData(:,1);
      if size(bcData,2) > 1
        bcVals = bcData(:,2);
      else
        bcVals = [];
      end

      if strcmp(obj.multiplierType,'P0')
        return
      end

      if nargin > 3
        if ~(domId == obj.idDomain(2))
          % not slave side
          return
        end
      end

      % get list of nodal slave dofs in the interface
      nodSlave = obj.mesh.local2glob{2}(obj.multEnts);
      % keep only active multipliers
      fldId = obj.dofm(2).getFieldId(bcPhysics);
      dofSlave = getLocalDoF(obj.dofm(2),nodSlave,fldId);
      isBCdof = ismember(bcData(:,1),dofSlave);
      bcDofs = bcDofs(~isBCdof);
      if ~isempty(bcVals)
        bcVals = bcVals(~isBCdof);
      end
    end
  end
end

