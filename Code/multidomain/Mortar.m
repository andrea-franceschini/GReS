classdef Mortar < handle
  % Base class for non conforming interfaces
  
  properties
    id              % interface identifier
    solvers         % instance of domains connected by the interface
    mesh            % instance of interfaceMesh class
    idDomain        % id of connected domains
    quadrature      % mortar integration object
    dofm            % dof manager instance of connected domains
    elements        % instance of Elements for the 2D interfaces
    multiplierType
    D               % slave mortar matrix
    M               % master mortar matrix
    outStruct       % wrapper for print utilities (mimics OutState)
    physic          % string id of mortar field
    totMult
    dirDofs        % list of nodes with dirichlet bcs enforced
  end

  methods
    function [obj] = Mortar(varargin)
      
      assert(nargin==2,'Wrong number of input arguments for mortar class')
      
      if isstruct(varargin{1})
        % complete mortar definition with input file
        inputStruct = varargin{1};
        domains = varargin{2};
        setMortar(obj,inputStruct,domains);
      else
        % simple mortar definition with 2 surfaces already defined
        % no need to use input file
        assert(isa(varargin{1},'Mesh'))
        mshMaster = varargin{1};
        mshSlave = varargin{2};
        obj.mesh = interfaceMesh(mshMaster,mshSlave);
      end
    end

%     function [r,c,v] = allocateMatrix(obj,sideID)
%       nEntries = nnz(obj.mesh.elemConnectivity)...
%         *obj.mesh.nN(sideID);
%       [r,c,v] = deal(zeros(nEntries,1));
%     end


    function applyBC(obj,idDomain,bound,t)
      bcList = bound.db.keys;

      for bcId = string(bcList)
        field = bound.getPhysics(bcId);
        if ~ismember(field,obj.physic)
          continue
        end

        if ~strcmp(bound.getType(bcId),'Dir')
          continue
          % only dirichlet bc has to be enforced to mortar blocks
        else
          if idDomain == obj.idDomain(1)
            applyBCmaster(obj,bcId,t)
          end
          if idDomain == obj.idDomain(2)
            applyBCslave(obj,bcId,t)
          end
        end
      end
    end


    function [bcDofs,bcVals] = removeSlaveBCdofs(obj,bcPhysics,bcData,domId)
      % this method updates the list of bc dofs and values removing dofs on
      % the slave interface (only for nodal multipliers)
      % this avoid possible overconstraint and solvability issues

      % bcData: nDofBC x 2 matrix.
      % first column -> dofs
      % second columns -> vals (optional)

      bcDofs = bcData(:,1);
      if size(bcData,2) > 1
        bcVals = bcData(:,2);
      else
        bcVals = [];
      end


      if nargin > 3
        if ~(domId == obj.idDomain(2))
          % not slave side
          return
        end
      end

      % get list of nodal slave dofs in the interface
      nodSlave = obj.mesh.local2glob{2};
      % keep only active multipliers
      fldId = obj.dofm(2).getFieldId(bcPhysics);
      dofSlave = getLocalDoF(obj.dofm(2),nodSlave,fldId);
      isBCdof = ismember(bcData(:,1),dofSlave);
      % update list of dirichlet nodes
      obj.dirDofs = unique([obj.dirDofs; bcDofs(isBCdof)]);

      if strcmp(obj.multiplierType,'P0')
        return
      end

      % remove slave dofs only for nodal multipliers on slave side

      bcDofs = bcDofs(~isBCdof);
      if ~isempty(bcVals)
        bcVals = bcVals(~isBCdof);
      end
    end


    function elem = getElem(obj,sideID,id)
      % get instance of element class on one the sides of the interface
      % Assumption: same element type in the entire interface
      % get istance of element class based on cell type
      type = obj.mesh.msh(sideID).surfaceVTKType(id);
      elem = obj.elements(sideID).getElement(type);
    end


    function dofs = getMultiplierDoF(obj,is)
      % return multiplier dofs associated with slave element #is inside
      % obj.mesh.msh(2)

      nc = obj.dofm(2).getDoFperEnt(obj.physic);
      if nargin > 1
        assert(isscalar(is),'Input id must be a scalar integer \n')
        if strcmp(obj.multiplierType,'P0')
          dofs = dofId(is,nc);
        else
          nodes = obj.mesh.msh(2).surfaces(is,:);
          dofs = dofId(nodes,nc);
        end
      else
        if strcmp(obj.multiplierType,'P0')
          dofs = dofId((1:obj.mesh.msh(2).nSurfaces)',nc);
        else
          dofs = dofId((1:obj.mesh.msh(2).nNodes)',nc);
        end
      end
    end


    function nDof = getNumbMultipliers(obj)
      % nDoFAct -> number of currently active multipliers in the interface
      % nDoFTot -> total number of dof in the whole interface

      nc = obj.dofm(2).getDoFperEnt(obj.physic);
      switch obj.multiplierType
        case 'P0'
          nDof = nc*obj.mesh.msh(2).nSurfaces;
        case {'dual','standard'} % nodal multipliers
          % get number of nodes in active cells
          nDof = nc*obj.mesh.msh(2).nNodes;
      end
    end


    function computeMortarInterpolation(obj)

      processMortarPairs(obj.quadrature); 

      inactiveMaster = find(~ismember(1:obj.mesh.msh(1).nSurfaces,...
        obj.quadrature.mortarPairs(:,2)));

      inactiveSlave = find(~ismember(1:obj.mesh.msh(2).nSurfaces,...
        obj.quadrature.mortarPairs(:,1)));

      % remove master elements
      obj.mesh.removeMortarSurface(1,inactiveMaster);

      % remove slave elements
      obj.mesh.removeMortarSurface(2,inactiveSlave);


    end

    %     function mat = getMatrix(obj,sideID,field)
    %       n = obj.dofm(sideID).getDoFperEnt(field);
    %       dofMult = dofId(1:obj.mesh.nEl(2),n);
    %       dof = obj.mesh.local2glob{sideID}(1:size(obj.mortarMatrix{sideID},2));
    %       dof = obj.dofm(sideID).getLocalDoF(dof,field);
    %       [j,i] = meshgrid(dof,dofMult);
    %       nr = n*obj.mesh.nEl(2);
    %       nc = obj.dofm(sideID).getNumDoF(field);
    %       vals = Discretizer.expandMat(obj.mortarMatrix{sideID},n);
    %       mat = sparse(i(:),j(:),vals(:),nr,nc); % minus sign!
    %     end

 
    function computeMortarMatrices(obj,~)

      % This method computes the cross grid mortar matrices between
      % connected interfaces

      % Also remove inactive slave/master elements after performing first
      % round of mortar integration

       % number of components per dof of interpolated physics
      if ~isempty(obj.dofm)
        ncomp = obj.dofm(2).getDoFperEnt(obj.physic);
      else
        ncomp = 1;
      end

      % logical index to track slave and master elements that do not really
      % participate to the mortar surface (fake connectivity)
      isInactiveSlave = false(obj.mesh.msh(2).nSurfaces,1);

      % get number of index entries for sparse matrices
      % overestimate number of sparse indices assuming all quadrilaterals
      nNmaster = obj.mesh.msh(1).surfaceNumVerts'*obj.mesh.elemConnectivity;

      switch obj.multiplierType
        case 'P0'
          N1 = sum(nNmaster);
          N2 = sum(obj.mesh.msh(2).surfaceNumVerts);
        otherwise
          N1 = nNmaster*obj.mesh.msh(2).surfaceNumVerts;
          N2 = sum(obj.mesh.msh(2).surfaceNumVerts.^2);
      end

      nm = (ncomp^2)*N1;
      ns = ncomp^2*N2;

      if ~isempty(obj.solvers)
        nDofMaster = obj.dofm(1).getNumDoF(obj.physic);
        nDofSlave = obj.dofm(2).getNumDoF(obj.physic);
        nDofMult = getNumbMultipliers(obj);
      else
        obj.multiplierType = 'dual'; 
        nDofMaster = obj.mesh.msh(1).nNodes;
        nDofSlave = obj.mesh.msh(2).nNodes;
        nDofMult = nDofSlave;
      end

      % define matrix assemblers
      locM = @(imult,imaster,Nmult,Nmaster) ...
        computeLocMaster(obj,imult,imaster,Nmult,Nmaster);
      locD = @(imult,islave,Dloc) ...
        computeLocSlave(obj,imult,islave,Dloc);

      asbM = assembler(nm,locM,nDofMult,nDofMaster);
      asbD = assembler(ns,locD,nDofMult,nDofSlave);

      for is = 1:obj.mesh.msh(2).nSurfaces
        masterElems = find(obj.mesh.elemConnectivity(:,is));
        if isempty(masterElems)
          continue
        end

        elSlave = getElem(obj,2,is);
        nN = elSlave.nNode;
        switch obj.multiplierType
          case 'P0'
            Dloc = zeros(ncomp,ncomp*nN);
          otherwise
            Dloc = zeros(ncomp*nN,ncomp*nN);
        end

        for im = masterElems'

          [Nslave,Nmaster,Nmult] = ...
            getMortarBasisFunctions(obj.quadrature,is,im);

          if isempty(Nmaster)
            % refine connectivity matrix
            %obj.mesh.elemConnectivity(im,is) = 0;
            continue
          end

          [Nslave,Nmaster,Nmult] = ...
            obj.reshapeBasisFunctions(ncomp,Nslave,Nmaster,Nmult);

          asbM.localAssembly(is,im,-Nmult,Nmaster);

          Dloc = Dloc + ...
            obj.quadrature.integrate(@(a,b) pagemtimes(a,'ctranspose',b,'none'),...
            Nmult,Nslave);
        end

        % if Dloc is empty, the current slave element is inactive remove
        % also corresponding master elements if it is not connected to any
        % other element
        if all(Dloc==0,"all")
          isInactiveSlave(is) = true;
        end

        asbD.localAssembly(is,is,Dloc);
      end

      % remove inactive slave elements
      removeMortarSurface(obj.mesh,2,isInactiveSlave);

      % find unconnected master elements
      isInactiveMaster = ~any(obj.mesh.elemConnectivity, 2);
      % remove unconnected elements
      removeMortarSurface(obj.mesh,1,isInactiveMaster);

      obj.M = asbM.sparseAssembly();
      obj.D = asbD.sparseAssembly();

      % check satisfaction of partition of unity (mortar consistency)
      
      % remove rows of inactive multipliers from Jmaster and Jslave
      if ~isempty(obj.solvers)
        dofMult = getMultiplierDoF(obj);
        obj.M = obj.M(dofMult,:);
        obj.D = obj.D(dofMult,:);
      end

      pu = sum([obj.M obj.D],2);
      assert(norm(pu)<1e-6,'Partiition of unity violated');
%       fprintf('Done computing mortar matrix in %.4f s \n',cputime-tIni)
    end


    function sideStr = getSide(obj,idDomain)
      % get side of the interface 'master' or 'slave' based on the
      % domain input id
      isMaster = obj.idDomain(1) == idDomain;
      isSlave = obj.idDomain(2) == idDomain;
      if isMaster
        sideStr = 'master';
      elseif isSlave
        sideStr = 'slave';
      elseif isMaster && isSlave
        sideStr = 'master'+'slave';
      else
        % consider the case where both sides belong to the same domain,
        % something like 'master_slave'
        error('Input domain not belonging to the interface');
      end
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
          Ml = Ns'*(Ns.*gpW');
          Dl = diag(Ns'*gpW');
          A = Ml\Dl;
          Nmult = NslaveIn*A;
      end
    end

    function finalizeOutput(obj)
      if ~isempty(obj.outStruct)
        obj.outStruct.VTK.finalize();
      end
    end

    function printState(obj,tOld,tNew)
      if isempty(obj.outStruct)
        return
      end
      cellData2D = [];
      pointData2D = [];
      if nargin == 2
        t = tOld;
        [cellData,pointData] = buildPrintStruct(obj);
        cellData2D = OutState.mergeOutFields(cellData2D,cellData);
        pointData2D = OutState.mergeOutFields(pointData2D,pointData);
        obj.VTK.writeVTKFile(t, [], [], pointData2D, cellData2D);
      elseif nargin == 3
        tList = obj.outStruct.tList;
        tID = obj.outStruct.tID;
        if tID <= length(tList)
          while tList(tID) <= tNew
            t = tList(tID);
            %Linear interpolation
            fac = (t - tOld)/(tNew - tOld);
            [cellData,pointData] = buildPrintStruct(obj,fac);
            cellData2D = OutState.mergeOutFields(cellData2D,cellData);
            pointData2D = OutState.mergeOutFields(pointData2D,pointData);
            tID = tID + 1;
            obj.outStruct.VTK.writeVTKFile(t, [], [], pointData2D, cellData2D);
            if tID > length(tList)
              break
            end
          end
          obj.outStruct.tID = tID;
        end
      end
    end
  end


  methods (Access = private)

    function setPrintUtils(obj,str,outState)
      
      if ~isfield(str,'Print')
        return
      else
        out = struct(...
          'name',[],...
          'tID', 1,...
          'tList',[],...
          'VTK',[]);

        out.name = str.Print.nameAttribute;
        out.VTK = VTKOutput(obj.mesh.msh(2),out.name);
        out.tList = outState.timeList;
        obj.outStruct = out;
      end
    end

    function checkInterfaceDisjoint(obj)
      % check that the nodes of mortar and slave side are disjoint
      if obj.idDomain(1)==obj.idDomain(2)
        % interface defined within the same domain 
        out = setdiff(obj.mesh.local2glob{1},obj.mesh.local2glob{2});
        if ~all(ismember(obj.mesh.local2glob{1},out))
          error('Nodes of master and slave side are not disjoint');
        end
      end
    end

    function setMortar(obj,inputStruct,domains)
      
      obj.solvers = domains;
      obj.idDomain = [inputStruct.Master.idAttribute;
        inputStruct.Slave.idAttribute];
      masterSurf = inputStruct.Master.surfaceTagAttribute;
      if isstring(masterSurf)
        masterSurf = str2num(masterSurf);
      end
      slaveSurf = inputStruct.Slave.surfaceTagAttribute;
      if isstring(slaveSurf)
        slaveSurf = str2num(slaveSurf);
      end

      obj.mesh = interfaceMesh(domains(1).grid.topology,domains(2).grid.topology,...
        masterSurf,slaveSurf);

      % check that master and slave node sets are disjoint
      checkInterfaceDisjoint(obj);

      obj.dofm = [domains(1).dofm;
        domains(2).dofm];

      quadType = inputStruct.Quadrature.typeAttribute;
      nG = inputStruct.Quadrature.nGPAttribute;
      if strcmp(quadType,'RBF')
        nInt = inputStruct.Quadrature.nIntAttribute;
      else
        nInt = [];
      end

      obj.setQuadrature(quadType,nG,nInt)

      computeMortarInterpolation(obj)

      setPrintUtils(obj,inputStruct,domains(2).outstate);

    end

    function setQuadrature(obj,quadType,nG,nInt)
      % define quadrature algorithm and element utilities

      switch quadType
        case 'RBF'
          assert(~isempty(nInt),['Missing number of interpolation points for' ...
            'RBF quadrature'])
          obj.quadrature = RBFquadratureNew(obj,nG,nInt);
        case 'SegmentBased'
          obj.quadrature = SegmentBasedQuadratureNew(obj,nG);
        case 'ElementBased'
          obj.quadrature = ElementBasedQuadratureNew(obj,nG);
      end
    end

  end


  methods (Static)
    function [interfaceStruct,domains] = buildInterfaces(fileName,domains)
      if domains(1).simparams.verbosity > 0
        fprintf('Mortar initialization... \n')
      end
      % read interface file and construct array of MeshGlue objects
      interfStr = readstruct(fileName);
      interfStr = interfStr.Interface;
      nInterfaces = numel(interfStr);
      interfaceStruct = cell(nInterfaces,1);
      for i = 1:nInterfaces
        idMaster = interfStr(i).Master.idAttribute;
        idSlave = interfStr(i).Slave.idAttribute;
        type = interfStr(i).Type;
        switch type
          case 'MeshTying'
            if (~isfield(interfStr(i),"Stabilization")) || ismissing(interfStr(i).Stabilization)
              % standard mesh tying with dual multipliers
              interfaceStruct{i} = MeshGlue(i,interfStr(i),domains([idMaster,idSlave]));
            elseif strcmp(interfStr(i).Stabilization.typeAttribute,'Jump')
              interfaceStruct{i} = MeshGlueJumpStabilization(i,interfStr(i), ...
                domains([idMaster,idSlave]));
            elseif strcmp(interfStr(i).Stabilization.typeAttribute,'Bubble')
              interfaceStruct{i} = MeshGlueBubbleStabilization(i,interfStr(i), ...
                domains([idMaster,idSlave]));
            else
              error('Invalid input argument for interface %i',i)
            end
          case 'ContactMortar'
            % not yet implemented!
            interfaceStruct{i} = ContactMortar(i,interfStr(i),domains([idMaster,idSlave]));
          case 'MeshTyingCondensation'
            interfaceStruct{i} = MeshGlueDual(i,interfStr(i),domains([idMaster,idSlave]));
          otherwise
            error(['Invalid interface law type for interface %i in file' ...
              '%s. \nAvailable types are: \nMeshTying \nFault'],i,fileName);
        end
        addInterface(domains(idMaster),i,interfaceStruct{i});
        addInterface(domains(idSlave),i,interfaceStruct{i});
      end
    end

    function varargout = reshapeBasisFunctions(nc,varargin)
      assert(numel(varargin)==nargout);
      varargout = cell(1,nargout);
      for i = 1:numel(varargin)
        varargout{i} = Mortar.reshapeBasisF(varargin{i},nc);
      end
    end


    function Nout = reshapeBasisF(basis,nc)
      % reshape basis functions to match number of components of the
      % physics
      % output: nc x nN x nG
      [ng,nn,nt] = size(basis);
      Nout = zeros(nc,nc*nn,ng,nt);
      for i = 1:nt
        % index pattern for matrix reshaping
        ind = repmat(linspace(1,nc^2,nc),1,nn)+(nc^2)*repelem(0:nn-1,1,nc);
        N = zeros(nc*nn,ng); % initialize single page
        N(ind(1:nc*nn),:) = repelem(basis(:,:,i)',nc,1);
        Nout(:,:,:,i) = reshape(N,[nc,nn*nc,ng]); % reshaped 3D matrix
      end
    end


  end
end

