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
      isFldMaster = isField(obj.dofm(1),obj.physic);
      isFldSlave = isField(obj.dofm(2),obj.physic);
      assert(isFldMaster,['MeshGlue physic not available for ' ...
        'master domain %i'],obj.idDomain(1));
      assert(isFldSlave,['MeshGlue physic not available for ' ...
        'slave domain %i'],obj.idDomain(2));

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


%     function computeMortarMatrices(obj)
% 
%       % This method computes the cross grid mortar matrices between
%       % connected interfaces
% 
%       % Also remove inactive slave/master elements after performing first
%       % round of mortar integration
% 
%       % number of components per dof of interpolated physics
%       if ~isempty(obj.dofm)
%         ncomp = obj.dofm(2).getDoFperEnt(obj.physic);
%       else
%         ncomp = 1;
%       end
% 
%       % get number of index entries for sparse matrices
%       % overestimate number of sparse indices assuming all quadrilaterals
%       nNmaster = obj.mesh.msh(1).surfaceNumVerts'*obj.mesh.elemConnectivity;
% 
%       switch obj.multiplierType
%         case 'P0'
%           N1 = sum(nNmaster);
%           N2 = sum(obj.mesh.msh(2).surfaceNumVerts);
%         otherwise
%           N1 = nNmaster*obj.mesh.msh(2).surfaceNumVerts;
%           N2 = sum(obj.mesh.msh(2).surfaceNumVerts.^2);
%       end
% 
%       nm = (ncomp^2)*N1;
%       ns = ncomp^2*N2;
% 
%       if ~isempty(obj.solvers)
%         nDofMaster = obj.dofm(1).getNumDoF(obj.physic);
%         nDofSlave = obj.dofm(2).getNumDoF(obj.physic);
%         nDofMult = getNumbMultipliers(obj);
%       else
%         obj.multiplierType = 'dual';
%         nDofMaster = obj.mesh.msh(1).nNodes;
%         nDofSlave = obj.mesh.msh(2).nNodes;
%         nDofMult = nDofSlave;
%       end
% 
%       % define matrix assemblers
%       locM = @(imult,imaster,Nmult,Nmaster,dJw) ...
%         computeLocMortar(obj,1,imult,imaster,Nmult,Nmaster,dJw);
%       locD = @(imult,islave,Nmult,Nslave,dJw) ...
%         computeLocMortar(obj,2,imult,islave,Nmult,Nslave,dJw);
% 
%       asbM = assembler(nm,locM,nDofMult,nDofMaster);
%       asbD = assembler(ns,locD,nDofMult,nDofSlave);
% 
%       for iPair = 1:obj.quadrature.numbMortarPairs
% 
%         is = obj.quadrature.mortarPairs(iPair,1);
%         im = obj.quadrature.mortarPairs(iPair,2);
% 
%         % retrieve mortar integration data
%         xiMaster = obj.quadrature.getMasterGPCoords(iPair);
%         xiSlave = obj.quadrature.getSlaveGPCoords(iPair);
%         dJw = obj.quadrature.getIntegrationWeights(iPair);
% 
%         [Nslave,Nmaster,Nmult] = ...
%           getMortarBasisFunctions(obj,im,is,xiMaster,xiSlave);
% 
%         [Nslave,Nmaster,Nmult] = ...
%           obj.reshapeBasisFunctions(ncomp,Nslave,Nmaster,Nmult);
% 
%         asbM.localAssembly(is,im,-Nmult,Nmaster,dJw);
%         asbD.localAssembly(is,is,Nmult,Nslave,dJw);
% 
%       end
% 
%       obj.M = asbM.sparseAssembly();
%       obj.D = asbD.sparseAssembly();
% 
%       pu = sum([obj.M obj.D],2);
%       assert(norm(pu)<1e-6,'Partiition of unity violated');
% 
%     end


    function [dofr,dofc,mat] = computeLocMortar(obj,side,imult,iu,Nmult,Nu,dJw)
      mat = pagemtimes(Nmult,'ctranspose',Nu,'none');
      mat = mat.*reshape(dJw,1,1,[]);
      mat = sum(mat,3);
      nodes = obj.mesh.local2glob{side}(obj.mesh.msh(side).surfaces(iu,:));
      fld = obj.dofm(side).getFieldId(obj.physic);
      dofc = obj.dofm(side).getLocalDoF(nodes,fld);
      dofr = getMultiplierDoF(obj,imult);
    end


    % function [rhsMaster,rhsSlave,rhsMult] = computeRhs(obj)
    %   % compute rhs contributions
    %   % reset rhs multiplier
    %   nMult = getNumbMultipliers(obj);
    %   obj.rhsMult = zeros(nMult,1);
    %   computeRhsMaster(obj);
    %   computeRhsSlave(obj);
    % end

    function applyBCmaster(obj,bc,t)
      ph = obj.solvers(1).bcs.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(1).getSolver(ph),bc,t);
      obj.rhsMaster(bcEnts) = 0;
      obj.Jmaster(:,bcEnts) = 0;
    end

    function applyBCslave(obj,bc,t)
      ph = obj.solvers(2).bcs.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(2).getSolver(ph),bc,t);
      obj.rhsSlave(bcEnts) = 0;
      obj.Jslave(:,bcEnts) = 0;
    end

    function updateState(obj,du)
      actMult = getMultiplierDoF(obj);
      n = numel(actMult);
      %du(1:3:end) = du(1:3:end); % fix sign convention
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

