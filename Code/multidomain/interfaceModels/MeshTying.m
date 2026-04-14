classdef MeshTying < InterfaceSolver

  % Standard mesh tying between non conforming interfaces
  
  properties
    D
    M
    stabilizationMat
    stabilizationScale = 1.0

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

      params = readInput(struct('stabilizationScale',1.0),input);
      obj.stabilizationScale = params.stabilizationScale;

      ncomp = obj.domains(MortarSide.slave).dofm.getNumberOfComponents(obj.coupledVariables);
      obj.nMult = ncomp * getNumberOfEntities(obj.multiplierLocation,...
        obj.grids(MortarSide.slave));

      obj.state.multipliers = zeros(obj.nMult,1);
      obj.state.iniMultipliers = zeros(obj.nMult,1);

      obj.stateOld = obj.state;

    end


    function updateState(obj,du)

      obj.state.multipliers = obj.state.multipliers + du(1:obj.nMult);

    end


    function assembleConstraint(obj,varargin)

      if isempty(obj.D)
        computeConstraintMatrices(obj);
      end

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

      m = MortarSide.master;
      s = MortarSide.slave;

      surfMaster = obj.grids(m).surfaces;
      surfSlave = obj.grids(s).surfaces;

      dofmMaster = getDoFManager(obj,m);
      dofmSlave =  getDoFManager(obj,s);

      elemPairs = obj.quadrature.interfacePairs;

      % number of components per dof of interpolated physics
      if ~isempty(obj.domains(m).dofm)
        ncomp = obj.domains(m).dofm.getNumberOfComponents(obj.coupledVariables);
      else
        ncomp = 1;
      end

      % get number of index entries for sparse matrices
      % overestimate number of sparse indices assuming all quadrilaterals

      % array with number of master nodes per slave
      nv = surfMaster.numVerts(elemPairs(:,m));
      nNMPS = accumarray(elemPairs(:,s),nv,[surfSlave.num,1]); 

      switch obj.multiplierLocation
        case entityField.cell
          N1 = sum(nNMPS);
          N2 = sum(surfSlave.numVerts);
        otherwise
          N1 = nNMPS'*surfSlave.numVerts;
          N2 = sum(surfSlave.numVerts.^2);
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

      topolMaster = getRowsMatrix(surfMaster.connectivity,1:surfMaster.num);
      topolSlave = getRowsMatrix(surfSlave.connectivity,1:surfSlave.num);


      for vtkSlave = surfSlave.vtkTypes

        elSlave = getElement(obj,vtkSlave,s);

        for vtkMaster = surfSlave.vtkTypes

          elMaster = getElement(obj,vtkSlave,m);

          % loop over pairs of connected master/slave elements
          for iPair = 1:obj.quadrature.numbInterfacePairs

            is = elemPairs(iPair,s);
            im = elemPairs(iPair,m);

            if surfSlave.VTKType(is) ~= vtkSlave; continue; end
            if surfMaster.VTKType(im) ~= vtkMaster; continue; end

            % get dofs
            nodeSlave = surfSlave.loc2glob(topolSlave(is,1:elSlave.nNode));
            usDof = dofmSlave.getLocalDoF(fldS,nodeSlave);
            nodeMaster = surfMaster.loc2glob(topolMaster(im,1:elMaster.nNode));
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
              getMortarBasisFunctions(obj.quadrature,im,is,elMaster,elSlave,xiMaster,xiSlave);

            [Ns,Nm,Nmult] = ...
              reshapeBasisFunctions(ncomp,Ns,Nm,Nmult);

            if ncomp > 1
              % vector field, local reference needed

              % rotation of multipliers
              R = getRotationMatrix(obj,s,is);

              % mortar coupling normal component
              Atm =  MortarQuadrature.integrate(f,Nmult,Nm,dJw);
              Ats =  MortarQuadrature.integrate(f,Nmult,Ns,dJw);

              % if strcmp(obj.quadrature.multiplierType,"dual")
              %   Ats = diag(sum(Ats,2));
              % end

              Atm = R'*Atm;
              Ats = R'*Ats;

              % assemble the local mortar matrix contribution
              asbM.localAssembly(tDof,umDof,Atm);
              asbD.localAssembly(tDof,usDof,-Ats);

            else

              Mloc = MortarQuadrature.integrate(f,Nmult,Nm,dJw);
              Dloc = MortarQuadrature.integrate(f,Nmult,Ns,dJw);
              asbM.localAssembly(tDof,umDof,Mloc); % minus sign!
              if strcmp(obj.quadrature.multiplierType,"dual")
                Dloc = diag(sum(Dloc,2));
              end

              asbD.localAssembly(tDof,usDof,-Dloc);

            end
          end
        end
      end

      obj.M = asbM.sparseAssembly();
      obj.D = asbD.sparseAssembly();

      % verify partition of unity
      pu = sum([obj.M obj.D],2);
      assert(norm(pu)<1e-7,'Partiition of unity violated');

    end


    % function [dofr,dofc,mat] = computeLocMortar(obj,side,imult,iu,Nmult,Nu,dJw)
    %   mat = pagemtimes(Nmult,'ctranspose',Nu,'none');
    %   mat = mat.*reshape(dJw,1,1,[]);
    %   mat = sum(mat,3);
    %   nodes = obj.interfMesh.local2glob{side}(obj.getMesh(MortarSide.slave).surfaces(iu,:));
    %   fld = obj.domains(side).dofm.getVariableId(obj.coupledVariables);
    %   dofc = obj.domains(side).dofm.getLocalDoF(fld,nodes);
    %   dofr = getMultiplierDoF(obj,imult);
    % end


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

    function writeSolution(obj,fac,tID)

      multCurr = obj.state.multipliers;
      multOld = obj.stateOld.multipliers;
      mult = fac*multCurr + (1-fac)*multOld;

      obj.outstate.results(tID).multipliers = mult;

    end


  end


  methods (Access = protected)

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
      dofM = dofMaster.getLocalDoF(fldM,nm);
      dofS = dofSlave.getLocalDoF(fldS,ns);

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
    

   function computeStabilizationMatrix(obj)

      m = MortarSide.master;
      s = MortarSide.slave;
      nComp = getDoFManager(obj,MortarSide.slave).getNumberOfComponents(obj.coupledVariables);

      edgeM = obj.grids(m).edges;
      edgeS = obj.grids(s).edges;
      surfM = obj.grids(m).surfaces;
      surfS = obj.grids(s).surfaces;

      % initialize matrix estimating number of entries
      % number of internal slave edges
      nes = sum(~edgeS.isBoundary);
      nEntries = 2*36*nes; % each cell should contribute at least two times
      nmult = getNumbDoF(obj);
      localKernel = @(S,e1,e2) assembleLocalStabilization(obj,S,e1,e2);
      asbH = assembler(nEntries,nmult,nmult,localKernel);

      intEdgeM = (find(~edgeM.isBoundary))';
      surf2edge = getRowsMatrix(surfS.surfaces2edges,1:surfS.num);

      % loop over internal master edges
      for eM = intEdgeM

        % get master node and element in the patch
        nodeM = edgeM.connectivity(eM,:);
        fM = edgeM.neighbors(eM,:);
        
        % get slave faces sharing support with master faces
        ii = ismember(obj.quadrature.interfacePairs(:,m),fM);
        fS = unique(obj.quadrature.interfacePairs(ii,s));

        if numel(fS) < 2
          continue
        end

        % average master elements area
        Am = mean(surfM.area(fM));

        % get edges internal in the macroelement
        eS = surf2edge(fS,:);
        eS = eS(eS > 0);
        [ueS,~,ic] = unique(eS);
        cnt = accumarray(ic,1);
        ieS = ueS(cnt == 2);

        if isempty(ieS)
          continue
        end

        connS = edgeS.connectivity(ieS,:);
        nodeS = unique(connS(:));

        % compute local schur complement approximation
        nodeM = surfM.loc2glob(nodeM);
        nodeS = surfS.loc2glob(nodeS);
        S = computeSchurLocal(obj,nodeM,nodeS,fS);

        % assemble stabilization matrix component
        for iesLoc = ieS'
          % loop over internal slave edges in the macroelement

          % get pair of slave faces sharing the edge
          f = edgeS.neighbors(iesLoc,:);

          % mean area of the slave faces
          As = mean(surfS.area(f));
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

    function [dofRow,dofCol,mat] = assembleLocalStabilization(obj,S,e1,e2)
      % assemble stabilization matrix S (in global coordinates) for
      % elements e1 and e2.

      nc = getDoFManager(obj,MortarSide.slave).getNumberOfComponents(obj.coupledVariables);
      dof1 = DoFManager.dofExpand(e1,nc);
      dof2 = DoFManager.dofExpand(e2,nc);

      % IMPORTANT
      % Rotation is not required since the Schur complement is already
      % rotated in the local frame

      % if nc > 1
      %   % vector field, rotation matrix needed
      % 
      %   % get average rotation matrix
      %   n1 = getNormal(obj.interfMesh,e1);
      %   n2 = getNormal(obj.interfMesh,e2);
      %   if abs(n1'*n2 -1) < 1e4*eps
      %     avgR = obj.interfMesh.computeRot(n1);
      %   else
      %     A1 = obj.interfMesh.msh(2).surfaceArea(e1);
      %     A2 = obj.interfMesh.msh(2).surfaceArea(e2);
      %     nAvg = n1*A1 + n2*A2;
      %     nAvg = nAvg/norm(nAvg);
      %     avgR = obj.interfMesh.computeRot(nAvg);
      %   end

        % apply rotation matrix to S
        %S = avgR'*S*avgR;

      %end

      % prepare matrix for full stick edge

      mat = obj.stabilizationScale*[S,-S;-S,S];
      dofRow = [dof1;dof2];
      dofCol = [dof1;dof2];

    end



  end

  methods (Static)

    function var = getCoupledVariables()
      var = [];
    end

  end
end



