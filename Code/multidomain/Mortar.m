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
    contactHelper
  end

  methods
    function [obj] = Mortar(varargin)
      
      assert(nargin==2,'Wrong number of input arguments for mortar class')
      
      if isstruct(varargin{1})
        % complete mortar definition with input file
        inputStruct = varargin{1};
        domains = varargin{2};
        if ~isfield(inputStruct,"Multiplier")
          obj.multiplierType = "P0";
        else
          obj.multiplierType = inputStruct.Multiplier.typeAttribute;
        end
        setMortar(obj,inputStruct,domains);
      else
        % simple mortar definition with 2 surfaces already defined
        % no need to use input file
        assert(isa(varargin{1},'Mesh'))
        mshMaster = varargin{1};
        mshSlave = varargin{2};
        obj.mesh = interfaceMesh(mshMaster,mshSlave);
      end
      obj.contactHelper = ContactHelper(obj);
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

        if ~strcmp(bound.getType(bcId),'Dirichlet')
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

    function computeMortarMatrices(obj)

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

      fldM = obj.dofm(1).getFieldId(obj.physic);
      fldS = obj.dofm(2).getFieldId(obj.physic);

      asbM = assembler(nm,nDofMult,nDofMaster);
      asbD = assembler(ns,nDofMult,nDofSlave);

      f = @(a,b) pagemtimes(a,'ctranspose',b,'none');

      for iPair = 1:obj.quadrature.numbMortarPairs

        is = obj.quadrature.mortarPairs(iPair,1);
        im = obj.quadrature.mortarPairs(iPair,2);

        nodeSlave = obj.mesh.local2glob{2}(obj.mesh.msh(2).surfaces(is,:));
        usDof = obj.dofm(2).getLocalDoF(nodeSlave,fldS);
        nodeMaster = obj.mesh.local2glob{1}(obj.mesh.msh(1).surfaces(im,:));
        umDof = obj.dofm(1).getLocalDoF(nodeMaster,fldM);
        tDof = getMultiplierDoF(obj,is);

        % retrieve mortar integration data
        xiMaster = obj.quadrature.getMasterGPCoords(iPair);
        xiSlave = obj.quadrature.getSlaveGPCoords(iPair);
        dJw = obj.quadrature.getIntegrationWeights(iPair);

        [Ns,Nm,Nmult] = ...
          getMortarBasisFunctions(obj,im,is,xiMaster,xiSlave);

        [Ns,Nm,Nmult] = ...
          obj.reshapeBasisFunctions(ncomp,Ns,Nm,Nmult);

        if ncomp > 1
          % vector field, local reference needed
          normIdx = 1:3:length(tDof);
          tangIdx = setdiff(1:length(tDof), 1:3:length(tDof));


          % rotation of multipliers
          R = getRotationMatrix(obj.contactHelper,is);
          NmultR = pagemtimes(Nmult,R);

          % Reduce the dimension of multiplier basis functions exploiting
          % the local definition of degrees of freedom

          % the normal component of the multipliers basis
          Nmult_n = -Nmult(1,normIdx,:);
          % tangential component of multiplier basis
          Nmult_t = NmultR(:,tangIdx,:);

          % get normal at the gauss points (for warped facets...?)
          normalNodes = obj.contactHelper.getNodalNormal(is);
          normal = pagemtimes(Ns,normalNodes);

          % operator selecting only tangential components of the
          % displacements
          T = eye(3) - pagemtimes(normal,'none',normal,'transpose');

          % normal and tangential component of displacement basis functions
          Nm_n = pagemtimes(normal,'transpose',Nm,'none');
          Nm_t = pagemtimes(T,Nm);
          Ns_n = pagemtimes(normal,'transpose',Ns,'none');
          Ns_t = pagemtimes(T,Ns);

          Atm_n =  MortarQuadrature.integrate(f,Nmult_n,Nm_n,dJw);
          Ats_n =  MortarQuadrature.integrate(f,Nmult_n,Ns_n,dJw);

          Atm_t =  MortarQuadrature.integrate(f,Nmult_t,Nm_t,dJw);
          Ats_t =  MortarQuadrature.integrate(f,Nmult_t,Ns_t,dJw);

          asbM.localAssembly(tDof(normIdx),umDof,-Atm_n);
          asbM.localAssembly(tDof(tangIdx),umDof,-Atm_t);
          asbD.localAssembly(tDof(normIdx),usDof,Ats_n);
          asbD.localAssembly(tDof(tangIdx),usDof,Ats_t);

        else
          Mloc = MortarQuadrature.integrate(f,Nmult,Nm,dJw);
          Dloc = MortarQuadrature.integrate(f,Nmult,Ns,dJw);
          asbM.localAssembly(tDof,umDof,-Mloc);
          asbD.localAssembly(tDof,usDof,Dloc);
        end

      end

      obj.M = asbM.sparseAssembly();
      obj.D = asbD.sparseAssembly();

      pu = sum([obj.M obj.D],2);
      assert(norm(pu)<1e-6,'Partiition of unity violated');

    end


    function [dofr,dofc,mat] = computeLocMortar(obj,side,imult,iu,Nmult,Nu,dJw)
      mat = pagemtimes(Nmult,'ctranspose',Nu,'none');
      mat = mat.*reshape(dJw,1,1,[]);
      mat = sum(mat,3);
      nodes = obj.mesh.local2glob{side}(obj.mesh.msh(side).surfaces(iu,:));
      fld = obj.dofm(side).getFieldId(obj.physic);
      dofc = obj.dofm(side).getLocalDoF(nodes,fld);
      dofr = getMultiplierDoF(obj,imult);
    end

    function [Nslave, Nmaster, Nmult, varargout] = getMortarBasisFunctions(obj,im,is,xiMaster,xiSlave)
      elemMaster = obj.getElem(1,im);
      elemSlave = obj.getElem(2,is);
      Nmaster = elemMaster.computeBasisF(xiMaster);
      Nslave = elemSlave.computeBasisF(xiSlave);
      Nmult = obj.computeMultiplierBasisF(is,Nslave);

      if nargout ==4
        % outout the bubble function on the slave side
        varargout{1} = elemSlave.computeBubbleBasisF(xiSlave);
      end
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
          obj.quadrature = RBFquadrature(obj,nG,nInt);
        case 'SegmentBased'
          obj.quadrature = SegmentBasedQuadrature(obj,nG);
        case 'ElementBased'
          obj.quadrature = ElementBasedQuadrature(obj,nG);
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

