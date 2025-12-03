classdef MeshTying < InterfaceSolver

  % Standard mesh tying between non conforming interfaces
  
  properties
    D
    M
    stabilizationMat
    stabilizationScale

  end

  methods

    function obj = MeshTying(id,domains,inputStruct)
      
      obj@InterfaceSolver(id,domains,inputStruct);

    end

    function registerInterface(obj,input)

      if ~isscalar(obj.coupledVariables)
        error("A Mesh tying interface can only couple one variable field." + ...
          "If you want to couple more than one variable field, create more than one interface")
      end

      obj.stabilizationScale = getXMLData(input.Quadrature,1.0,"stabilizationScale");

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

      % reset the jacobian blocks
      obj.setJmu(MortarSide.slave, []);
      obj.setJmu(MortarSide.master, []);
      obj.setJum(MortarSide.slave, []);
      obj.setJum(MortarSide.master, []);

      % overwrite current jacobian
      addJum(obj,MortarSide.master,obj.M');
      addJmu(obj,MortarSide.master,obj.M);

      % add slave contribution
      addJum(obj,MortarSide.slave,obj.D');
      addJmu(obj,MortarSide.slave,obj.D);

      [rhsMaster,rhsSlave,rhsMult] = computeRhs(obj);

      addRhs(obj,MortarSide.master,rhsMaster);
      addRhs(obj,MortarSide.slave,rhsSlave);

      obj.rhsConstraint = rhsMult;

      % apply stabilization for P0 multipliers
      if obj.multiplierLocation == entityField.surface
        
        computeStabilizationMatrix(obj);
        rhsStab = -obj.stabilizationMat*...
          (obj.state.multipliers-obj.state.iniMultipliers);

        obj.Jconstraint = -obj.stabilizationMat;

        rhsMult = rhsMult + rhsStab;
      end

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
      for iPair = 1:obj.quadrature.numbInterfacePairs

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
          reshapeBasisFunctions(ncomp,Ns,Nm,Nmult);

        if ncomp > 1
          % vector field, local reference needed

          normIdx = 1:3:length(tDof);
          tangIdx = setdiff(1:length(tDof), 1:3:length(tDof));

          % rotation of multipliers
          R = getRotationMatrix(obj.interfMesh,is);
          NmultR = pagemtimes(Nmult,R);

          % Reduce the dimension of multiplier basis functions exploiting
          % the local definition of degrees of freedom

          % the normal component of the multipliers basis
          Nmult_n = -Nmult(1,normIdx,:);
          % tangential component of multiplier basis
          Nmult_t = NmultR(:,tangIdx,:);

          % get normal at the gauss points (for warped facets...?)
          % move this call to interfaceMesh()
          normalNodes = obj.interfMesh.getNodalNormal(is);
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
          if strcmp(obj.quadrature.multiplierType,"dual")
            Dloc = diag(sum(Dloc,2));
          end

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


    function [surfaceStr,pointStr] = writeVTK(obj,fac,varargin)

      multCurr = obj.state.multipliers;
      multOld = obj.stateOld.multipliers;
      mult = fac*multCurr + (1-fac)*multOld;

      surfaceStr = [];
      pointStr = [];

      nComp = getDoFManager(obj,MortarSide.slave).getNumberOfComponents(obj.coupledVariables);

      mult = (reshape(mult,nComp,[]))';

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


  methods (Access = protected)


    function computeStabilizationMatrix(obj)

      if ~isempty(obj.stabilizationMat)
        % compute stabilization matrix only once for all edges
        % retrieve row-col needing stabilization at each time step
        return
      end

      nComp = getDoFManager(obj,MortarSide.slave).getNumberOfComponents(obj.coupledVariables);

      % initialize matrix estimating number of entries
      % number of internal slave elements
      nes = sum(all(obj.interfMesh.e2f{2},2));
      nEntries = 2*36*nes; % each cell should contribute at least two times
      nmult = getNumbDoF(obj);
      localKernel = @(S,e1,e2) assembleLocalStabilization(obj,S,e1,e2);
      asbH = assembler(nEntries,nmult,nmult,localKernel);

      % get list of internal master edges
      inEdgeMaster = find(all(obj.interfMesh.e2f{1},2));

      for ieM = inEdgeMaster'
        % loop over internal master edges

        % get 2 master faces sharing internal edge ie
        fM = obj.interfMesh.e2f{1}(ieM,:);
        assert(numel(fM)==2,['Unexpected number of connected faces for' ...
          'master edge %i. Expected 2.'], ieM);

        % get slave faces sharing support with 2 master faces
        fS = unique([find(obj.interfMesh.elemConnectivity(fM(1),:)),...
          find(obj.interfMesh.elemConnectivity(fM(2),:))]);

        if numel(fS) < 2
          continue
        end

        % average master elements area
        Am = mean(obj.interfMesh.msh(1).surfaceArea(fM));

        % get internal edges of slave faces
        eS = unique(obj.interfMesh.f2e{2}(fS,:));
        id = all(ismember(obj.interfMesh.e2f{2}(eS,:),fS),2);
        ieS = eS(id);

        % get master/slave nodes in the local macroelement
        nM = obj.interfMesh.e2n{1}(ieM,:);
        nS = unique(obj.interfMesh.e2n{2}(eS,:));

        % compute local schur complement approximation
        S = computeSchurLocal(obj,nM,nS,fS);

        % assemble stabilization matrix component
        for iesLoc = ieS'
          % loop over internal slave edges in the macroelement

          % get pair of slave faces sharing the edge
          f = obj.interfMesh.e2f{2}(iesLoc,:);

          % mean area of the slave faces
          As = mean(obj.interfMesh.msh(2).surfaceArea(f));
          fLoc1 = dofId(find(fS==f(1)),nComp);
          fLoc2 = dofId(find(fS==f(2)),nComp);

          % local schur complement for macroelement pair of slave faces
          Sloc = 0.5*(Am/As)*(S(fLoc1,fLoc1)+S(fLoc2,fLoc2));
          asbH.localAssembly(Sloc,f(1),f(2));
        end
      end

      obj.stabilizationMat = asbH.sparseAssembly();

      assert(norm(sum(obj.stabilizationMat,2))<1e-8, 'Stabilization matrix is not locally conservative')
    end

    function S = computeSchurLocal(obj,nm,ns,fs)
      % compute approximate schur complement for local nonconforming
      % patch of element (in global reference)
      % input: nm/ns local master/slave node indices
      % fs: local slave faces indices

      % get slave and master dof to access jacobian
      dofMaster = getDoFManager(obj,MortarSide.master);
      dofSlave = getDoFManager(obj,MortarSide.slave);
      fldM = getVariableId(dofMaster,obj.coupledVariables);
      fldS = getVariableId(dofSlave,obj.coupledVariables);
      dofM = dofMaster.getLocalDoF(fldM,obj.interfMesh.local2glob{1}(nm));
      dofS = dofSlave.getLocalDoF(fldS,obj.interfMesh.local2glob{2}(ns));

      % get local mortar matrices
      nc = getDoFManager(obj,MortarSide.slave).getNumberOfComponents(obj.coupledVariables);
      Dloc = obj.D(DoFManager.dofExpand(fs,nc),dofS);
      Mloc = obj.M(DoFManager.dofExpand(fs,nc),dofM);
      V = [Dloc, Mloc];              % minus sign!
      %V = Discretizer.expandMat(V,nc);

      % get local jacobian
      Km = obj.domains(1).J{fldM,fldM}(dofM,dofM);
      Ks = obj.domains(2).J{fldS,fldS}(dofS,dofS);
      Kloc = diag([1./diag(Ks);1./diag(Km)]);

      S = obj.stabilizationScale*V*(Kloc*V');  % compute Schur complement

    end


    function [dofRow,dofCol,mat] = assembleLocalStabilization(obj,S,e1,e2)
      % assemble stabilization matrix S (in global coordinates) for
      % elements e1 and e2.

      nc = getDoFManager(obj,MortarSide.slave).getNumberOfComponents(obj.coupledVariables);
      dof1 = DoFManager.dofExpand(e1,nc);
      dof2 = DoFManager.dofExpand(e2,nc);

      if nc > 1
        % vector field, rotation matrix needed

        % get average rotation matrix
        n1 = getNormal(obj.interfMesh,e1);
        n2 = getNormal(obj.interfMesh,e2);
        if abs(n1'*n2 -1) < 1e4*eps
          avgR = obj.interfMesh.computeRot(n1);
        else
          A1 = obj.interfMesh.msh(2).surfaceArea(e1);
          A2 = obj.interfMesh.msh(2).surfaceArea(e2);
          nAvg = n1*A1 + n2*A2;
          nAvg = nAvg/norm(nAvg);
          avgR = obj.interfMesh.computeRot(nAvg);
        end

        % apply rotation matrix to S
        S = avgR'*S*avgR;

      end

      dofRow = [dof1;dof2];
      dofCol = [dof1;dof2];

      mat = [S,-S;-S,S];

    end


  end

  methods (Static)

    function var = getCoupledVariables()
      var = [];
    end

  end
end



