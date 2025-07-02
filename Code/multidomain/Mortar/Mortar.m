classdef Mortar < handle
  % This class implement some basic operation for mortar interpolation
  
  properties
    outStruct   
    solvers
    mesh         % instance of interfaceMesh class
    idDomain
    Jmaster
    Jslave
    Jmult
    rhsMaster
    rhsSlave
    rhsMult
    quadrature % integration scheme used for the interface
    quadrature2
    % (RBF or ElementBased)
    dofm
    nFld          
    mortarMatrix
    elements
    multiplierType = 'dual'
  end

  methods
    function [obj] = Mortar(inputStruct,domains)
      obj.solvers = [domains(1).Discretizer,domains(2).Discretizer];
      obj.idDomain = [inputStruct.Master.idAttribute;
                inputStruct.Slave.idAttribute];
      surfId = {inputStruct.Master.surfaceTagAttribute;
        inputStruct.Slave.surfaceTagAttribute};
      obj.mesh = interfaceMesh(domains,surfId);
      checkInterfaceDisjoint(obj);
      obj.dofm = [domains(1).DoFManager;
                        domains(2).DoFManager];
      switch inputStruct.Quadrature.typeAttribute
        case 'RBF'
          nG = inputStruct.Quadrature.nGPAttribute;
          nInt = inputStruct.Quadrature.nIntAttribute;
          obj.elements = [Elements(obj.mesh.msh(1),nG),...
            Elements(obj.mesh.msh(2),nG)];
          obj.quadrature = RBFquadrature(obj,nInt);
        case 'ElementBased'
           nG = inputStruct.Quadrature.nGPAttribute;
          obj.elements = [Elements(obj.mesh.msh(1),nG),...
            Elements(obj.mesh.msh(2),nG)];
          obj.quadrature = ElementBasedQuadrature(obj,nG);
          %obj.quadrature2 = RBFquadrature(obj,6);
        case 'SegmentBased'
          obj.quadrature = SegmentBasedQuadrature(obj,inputStruct.Quadrature.nGPAttribute);
          nG = 3; 
          % nG for dual multipliers inversion
          obj.elements = [Elements(obj.mesh.msh(1),nG),...
            Elements(obj.mesh.msh(2),nG)];
      end
      setPrintUtils(obj,inputStruct,domains(2).OutState);
    end

    function [r,c,v] = allocateMatrix(obj,sideID)
      nEntries = nnz(obj.mesh.elemConnectivity)...
        *obj.mesh.nN(sideID);
      [r,c,v] = deal(zeros(nEntries,1));
    end


    function applyBC(obj,idDomain,bound,t)
      bcList = bound.db.keys;

      for bc = string(bcList)
        field = bound.getPhysics(bc);
        if ~ismember(field,obj.physics)
          continue
        end

        if ~strcmp(bound.getType(bc),'Dir')
          continue
          % only dirichlet bc has to be enforced to mortar blocks
        else
          if idDomain == obj.idDomain(1)
            applyBCmaster(obj,bound,bc,t)
          end
          if idDomain == obj.idDomain(2)
            applyBCslave(obj,bound,bc,t)
          end
        end
      end
    end

    function elem = getElem(obj,sideID,id)
      % get instance of element class on one the sides of the interface
      % Assumption: same element type in the entire interface
      % get istance of element class based on cell type
      type = obj.mesh.msh(sideID).surfaceVTKType(id);
      elem = obj.elements(sideID).getElement(type);
    end

    function mat = getMatrix(obj,sideID,field)
      n = obj.dofm(sideID).getDoFperEnt(field);
      dofMult = dofId(1:obj.mesh.nEl(2),n);
      dof = obj.mesh.local2glob{sideID}(1:size(obj.mortarMatrix{sideID},2));
      fld = obj.dofm(sideID).getFieldId(field);
      dof = obj.dofm(sideID).getLocalDoF(dof,fld);
      [j,i] = meshgrid(dof,dofMult);
      nr = n*obj.mesh.nEl(2);
      nc = obj.dofm(sideID).getNumDoF(field);
      vals = Discretizer.expandMat(obj.mortarMatrix{sideID},n);
      mat = sparse(i(:),j(:),vals(:),nr,nc); % minus sign!
    end

    function computeMortarMatrices(obj,~)
      fprintf('Computing mortar matrices...\n')
      tIni = cputime;
      % assumption: each element has only one bubble face
      % loop over slave faces and:
      % 1) compute Aub, Abu and Abb local matrix from the neighbor cell
      % 2) compute local M, D and Db
      % 3) assemble static condensation blocks and Jmult

      % differently from the base method of the mortar class, here global
      % dof indexing is used for the assembled matrices

      ncomp = obj.dofm(2).getDoFperEnt(obj.physics);

      % get number of index entries for sparse matrices
      % overestimate number of sparse indices assuming all quadrilaterals
      cellsSlave = obj.mesh.getActiveCells(2);
      %cellsMaster = obj.mesh.getActiveCells(1);
      nNmaster = obj.mesh.msh(1).surfaceNumVerts'*obj.mesh.elemConnectivity;
      switch obj.multiplierType
        case 'P0'
          N1 = sum(nNmaster(cellsSlave));
          N2 = sum(obj.mesh(2).surfaceNumVerts(cellsSlave));
        otherwise
          N1 = nNmaster*obj.mesh.msh(2).surfaceNumVerts;
          N2 = sum(obj.mesh.msh(2).surfaceNumVerts(cellsSlave).^2);
      end

      nm = (ncomp^2)*N1;
      ns = ncomp^2*N2;

      % get number of dofs for each block
      nDofMaster = obj.dofm(1).getNumDoF(obj.physics);
      nDofSlave = obj.dofm(2).getNumDoF(obj.physics);
      nDofMult = getNumbMultipliers(obj);

      % define matrix assemblers
      locM = @(imult,imaster,Nmult,Nmaster) ...
        computeLocMaster(obj,imult,imaster,Nmult,Nmaster);
      locD = @(imult,islave,Dloc) ...
        computeLocSlave(obj,imult,islave,Dloc);

      asbM = assembler(nm,locM,nDofMult,nDofMaster);
      asbD = assembler(ns,locD,nDofMult,nDofSlave);

      for i = 1:obj.mesh.nEl(2)
        is = cellsSlave(i);
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

          asbM.localAssembly(i,im,-Nmult,Nmaster);

          Dloc = Dloc + ...
            obj.quadrature.integrate(@(a,b) pagemtimes(a,'ctranspose',b,'none'),...
            Nmult,Nslave);
        end

        asbD.localAssembly(i,is,Dloc);
      end

      obj.Jmaster{1} = asbM.sparseAssembly();
      obj.Jslave{1} = asbD.sparseAssembly();

      % check satisfaction of partition of unity (mortar consistency)
      pu = sum([obj.Jmaster{1} obj.Jslave{1}],2);
      assert(norm(pu)<1e-6,'Partiition of unity violated');
      fprintf('Done computing mortar matrix in %.4f s \n',cputime-tIni)
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

    function finalizeOutput(obj)
      obj.outStruct.VTK.finalize();
    end

    function printState(obj,tOld,tNew)
      if isempty(obj.outStruct)
        return
      end
      cellData2D = [];
      pointData2D = [];
      if nargin == 2
        t = tOld;
        for i = 1:obj.nFld
          [cellData,pointData] = buildPrintStruct(obj,i);
          cellData2D = OutState.mergeOutFields(cellData2D,cellData);
          pointData2D = OutState.mergeOutFields(pointData2D,pointData);
        end
        obj.VTK.writeVTKFile(t, [], [], pointData2D, cellData2D);
      elseif nargin == 3
        tList = obj.outStruct.tList;
        tID = obj.outStruct.tID;
        if tID <= length(tList)
          while tList(tID) <= tNew
            t = tList(tID);
            %Linear interpolation
            fac = (t - tOld)/(tNew - tOld);
            for i = 1:obj.nFld
              [cellData,pointData] = buildPrintStruct(obj,i,fac);
              cellData2D = OutState.mergeOutFields(cellData2D,cellData);
              pointData2D = OutState.mergeOutFields(pointData2D,pointData);
            end
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
      if obj.idDomain(1)==obj.idDomain(2)
        % interface defined within the same domain 
        out = setdiff(obj.mesh.local2glob{1},obj.mesh.local2glob{2});
        if ~all(out==obj.mesh.local2glob{1})
          error('Nodes of master and slave side are not disjoint');
        end
      end
    end
  end


  methods (Static)
    function [interfaceStruct,modelStruct] = buildInterfaceStruct(fileName,modelStruct)
      fprintf('Mortar initialization... \n')
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
            if (~isfield(interfStr,"Stabilization"))
              % standard mesh tying with dual multipliers
              interfaceStruct{i} = MeshGlue(i,interfStr(i),modelStruct([idMaster,idSlave]));
            elseif strcmp(interfStr.Stabilization.typeAttribute,'Jump')
              interfaceStruct{i} = MeshGlueJumpStabilization(i,interfStr(i), ...
                modelStruct([idMaster,idSlave]));
            elseif strcmp(interfStr.Stabilization.typeAttribute,'Bubble')
              interfaceStruct{i} = MeshGlueBubbleStabilization(i,interfStr(i), ...
                modelStruct([idMaster,idSlave]));
            else
              error('Invalid input argument for interface %i',i)
            end
          case 'Fault'
            % not yet implemented!
            interfaceStruct{i} = Fault();
          case 'MeshTyingCondensation'
            interfaceStruct{i} = MeshGlueDual(i,interfStr(i),modelStruct([idMaster,idSlave]));
          otherwise
            error(['Invalid interface law type for interface %i in file' ...
              '%s. \nAvailable types are: \nMeshTying \nFault'],i,fileName);
        end
        addInterface(modelStruct(idMaster).Discretizer,i,interfaceStruct{i});
        addInterface(modelStruct(idSlave).Discretizer,i,interfaceStruct{i});
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
      % reshape basis functions to obtain displacement shape function
      % input: nG x nN matrix
      % output: 3 x nN x nG
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

