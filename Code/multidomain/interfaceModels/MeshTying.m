classdef MeshTying < InterfaceSolver

  % Standard mesh tying between non conforming interfaces
  
  properties

    coupledVariables
    D
    M

  end

  methods

    function obj = MeshTying(id,domains,inputStruct)
      
      obj@InterfaceSolver(id,domains,inputStruct);

    end

    function registerInterface(obj,varargin)

      if ~isscalar(obj.coupledVariables)
        error("A Mesh tying interface can only couple one variable field." + ...
          "If you want to couple more than one variable field, create more than one interface")
      end

      ncomp = obj.domains(2).dofm.getNumberOfComponents(obj.coupledVariables);
      obj.nMult = ncomp * getNumberOfEntities(obj.multiplierLocation,...
        obj.getMesh(MortarSide.slave));

      obj.state.multipliers = zeros(obj.nMult,1);
      obj.state.iniMultipliers = zeros(obj.nMult,1);

      obj.stateOld = obj.state;

    end


    function updateState(obj,du)

      obj.state.multipliers = obj.state.multipliers + du(1:obj.nMult);

    end


    function assembleConstraint(obj)

      computeConstraintMatrices(obj);

      % overwrite current jacobian
      setJum(obj,MortarSide.master,obj.M');
      setJmu(obj,MortarSide.master,obj.M);

      % add slave contribution
      addJum(obj,MortarSide.slave,obj.D');
      addJmu(obj,MortarSide.slave,obj.D);

      [rhsMaster,rhsSlave,rhsMult] = computeRhs(obj);

      addRhs(obj,MortarSide.master,rhsMaster);
      addRhs(obj,MortarSide.slave,rhsSlave);

      obj.rhsConstraint = rhsMult;

    end


    function computeConstraintMatrices(obj)

      % This method computes the cross grid mortar matrices between
      % connected interfaces

      meshMaster = getMesh(obj,MortarSide.master);
      meshSlave = getMesh(obj,MortarSide.slave);

      dofmMaster = getDoFManager(obj,MortarSide.master);
      dofmSlave =  getDoFManager(obj,MortarSide.slave);


      % number of components per dof of interpolated physics
      if ~isempty(obj.domains(2).dofm)
        ncomp = obj.domains(2).dofm.getNumberOfComponents(obj.coupledVariables);
      else
        ncomp = 1;
      end

      % get number of index entries for sparse matrices
      % overestimate number of sparse indices assuming all quadrilaterals
      nNmaster = meshMaster.surfaceNumVerts'*obj.interfMesh.elemConnectivity;

      switch obj.multiplierLocation
        case entityField.cell
          N1 = sum(nNmaster);
          N2 = sum(meshSlave.surfaceNumVerts);
        otherwise
          N1 = nNmaster*meshSlave.surfaceNumVerts;
          N2 = sum(meshSlave.surfaceNumVerts.^2);
      end

      nm = (ncomp^2)*N1;
      ns = ncomp^2*N2;

      nDofMaster = dofmMaster.getNumbDoF(obj.coupledVariables);
      nDofSlave = dofmSlave.getNumbDoF(obj.coupledVariables);
      nDofMult = obj.nMult;


      fldM = dofmMaster.getVariableId(obj.coupledVariables);
      fldS = dofmSlave.getVariableId(obj.coupledVariables);

      asbM = assembler(nm,nDofMult,nDofMaster);
      asbD = assembler(ns,nDofMult,nDofSlave);

      f = @(a,b) pagemtimes(a,'ctranspose',b,'none');

      % loop over pairs of connected master/slave elements
      for iPair = 1:obj.quadrature.numbMortarPairs

        is = obj.quadrature.interfacePairs(iPair,1);
        im = obj.quadrature.interfacePairs(iPair,2);

        % get dofs 
        nodeSlave = obj.interfMesh.local2glob{2}(meshSlave.surfaces(is,:));
        usDof = dofmSlave.getLocalDoF(fldS,nodeSlave);
        nodeMaster = obj.interfMesh.local2glob{1}(meshMaster.surfaces(im,:));
        umDof = dofmMaster.getLocalDoF(fldM,nodeMaster);
        tDof = getMultiplierDoF(obj,is);

        % retrieve mortar integration data
        % gp reference coordinate on master element
        xiMaster = obj.quadrature.getMasterGPCoords(iPair);
        % gp reference coordinate on slave element
        xiSlave = obj.quadrature.getSlaveGPCoords(iPair);
        % determinant of jacobian for each gp
        dJw = obj.quadrature.getIntegrationWeights(iPair);

        % compute slave/master/multipliers basis functions
        [Ns,Nm,Nmult] = ...
          getMortarBasisFunctions(obj.quadrature,im,is,xiMaster,xiSlave);

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
          % move this call to interfaceMesh()
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

          % mortar coupling normal component
          Atm_n =  MortarQuadrature.integrate(f,Nmult_n,Nm_n,dJw);
          Ats_n =  MortarQuadrature.integrate(f,Nmult_n,Ns_n,dJw);

          % mortar coupling tangential component
          Atm_t =  MortarQuadrature.integrate(f,Nmult_t,Nm_t,dJw);
          Ats_t =  MortarQuadrature.integrate(f,Nmult_t,Ns_t,dJw);

          % assemble the local mortar matrix contribution
          asbM.localAssembly(tDof(normIdx),umDof,-Atm_n);
          asbM.localAssembly(tDof(tangIdx),umDof,-Atm_t);
          asbD.localAssembly(tDof(normIdx),usDof,Ats_n);
          asbD.localAssembly(tDof(tangIdx),usDof,Ats_t);

        else

          Mloc = MortarQuadrature.integrate(f,Nmult,Nm,dJw);
          Dloc = MortarQuadrature.integrate(f,Nmult,Ns,dJw);
          asbM.localAssembly(tDof,umDof,-Mloc); % minus sign!
          asbD.localAssembly(tDof,usDof,Dloc);

        end
      end

      obj.M = asbM.sparseAssembly();
      obj.D = asbD.sparseAssembly();

      % verify partition of unity 
      pu = sum([obj.M obj.D],2);
      assert(norm(pu)<1e-7,'Partiition of unity violated');

    end


    function [dofr,dofc,mat] = computeLocMortar(obj,side,imult,iu,Nmult,Nu,dJw)
      mat = pagemtimes(Nmult,'ctranspose',Nu,'none');
      mat = mat.*reshape(dJw,1,1,[]);
      mat = sum(mat,3);
      nodes = obj.interfMesh.local2glob{side}(obj.getMesh(MortarSide.slave).surfaces(iu,:));
      fld = obj.domains(side).dofm.getVariableId(obj.coupledVariables);
      dofc = obj.domains(side).dofm.getLocalDoF(fld,nodes);
      dofr = getMultiplierDoF(obj,imult);
    end


    function [rhsMaster,rhsSlave,rhsMult] = computeRhs(obj)

      % retrieve variable field arrays
      varMaster = getVariableField(obj,MortarSide.master);
      varSlave = getVariableField(obj,MortarSide.slave);

      % retrieve active multipliers
      actMult = getMultiplierDoF(obj);
      iniMult = obj.state.iniMultipliers(actMult);
      mult = obj.state.multipliers(actMult);

      % compute rhs terms
      rhsMaster = obj.M' * (mult - iniMult);
      rhsSlave = obj.D' * (mult - iniMult);

      rhsMult = obj.M * varMaster + obj.D * varSlave;

    end


    function computeMortarInterpolation(obj)

      processMortarPairs(obj.quadrature); 

      inactiveMaster = ~ismember(1:obj.getMesh(MortarSide.master).nSurfaces,...
        obj.quadrature.interfacePairs(:,2));
      [~, ~, obj.quadrature.interfacePairs(:,2)] = ...
        unique(obj.quadrature.interfacePairs(:,2));

      inactiveSlave = ~ismember(1:obj.getMesh(MortarSide.slave).nSurfaces,...
        obj.quadrature.interfacePairs(:,1));
      [~, ~, obj.quadrature.interfacePairs(:,1)] = ...
        unique(obj.quadrature.interfacePairs(:,1));

      % remove master elements
      obj.interfMesh.removeMortarSurface(1,inactiveMaster);

      % remove slave elements
      obj.interfMesh.removeMortarSurface(2,inactiveSlave);


    end


    function sideStr = getSide(obj,idDomain)
      % get side of the interface 'master' or 'slave' based on the
      % domains input id
      isMaster = obj.idDomain(1) == idDomain;
      isSlave = obj.idDomain(2) == idDomain;
      if isMaster
        sideStr = 'master';
      elseif isSlave
        sideStr = 'slave';
      elseif isMaster && isSlave
        sideStr = 'master'+'slave';
      else
        % consider the case where both sides belong to the same domains,
        % something like 'master_slave'
        error('Input domains not belonging to the interface');
      end
    end

    function [surfaceStr,pointStr] = writeVTK(obj,fac,varargin)

      multCurr = obj.state.multipliers;
      multOld = obj.stateOld.multipliers;
      mult = fac*multCurr + (1-fac)*multOld;

      surfaceStr = [];
      pointStr = [];

      % append state variable to output structure
      if obj.multiplierLocation == entityField.surface
        surfaceStr(1).name = 'multiplier';
        surfaceStr(1).data = mult;
      elseif obj.multiplierLocation == entityField.node
        pointStr(1).name = 'multiplier';
        pointStr(1).data = mult;
      end
    end

    function writeMatFile(obj,fac,tID)

      multCurr = obj.state.multipliers;
      multOld = obj.stateOld.multipliers;
      mult = fac*multCurr + (1-fac)*multOld;

      obj.outstate.results(tID).multipliers = mult;

    end

  end


  methods (Static)
 
    function varargout = reshapeBasisFunctions(nc,varargin)
      assert(numel(varargin)==nargout);
      varargout = cell(1,nargout);
      for i = 1:numel(varargin)
        varargout{i} = reshapeBasisF(varargin{i},nc);
      end
    end

  end
end

