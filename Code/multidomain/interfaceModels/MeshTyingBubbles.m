classdef MeshTyingBubbles < MeshTying

  % Mesh tying between non conforming interfaces using static condensation
  
  properties
    localFaceIndex % slave faces with local index in neighboring cell
    Dbubble

  end

  methods

    function obj = MeshTyingBubbles(id,domains,inputStruct)
      
      obj@MeshTying(id,domains,inputStruct);

    end

    function registerInterface(obj,input)

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

      setBubbleFaces(obj);

    end

    function updateState(obj,dSol)

      solv = obj.domains(2).getPhysicsSolver("Poromechanics");

      % get increment of displacement
      du = solv.domain.state.data.displacements - solv.domain.stateOld.data.displacements;

      % update multipliers
      updateState@MeshTying(obj,dSol);

      dmult = obj.state.multipliers - obj.stateOld.multipliers;

      % Enhance strain with bubble contribution
      for i = 1:getMesh(obj,MortarSide.slave).nSurfaces

        % get index of surface and neighbor cell
        el = obj.interfMesh.f2c{2}(i);

        l = solv.cell2stress(el);

        % get index of cell displacement dofs in du
        dof = getDoFID(solv.mesh,el);

        vtkId = solv.mesh.cellVTKType(el);
        elem = getElement(solv.elements,vtkId);
        nG = elem.GaussPts.nNode;

        % get local Db matrix
        Dbub = diag(repelem(obj.Dbubble(i),3));

        % the last input is the time increment. consider adding it as a
        % property of State.m
        [~,~,Kub,Kbb,Bb] = computeLocalStiffBubble(solv,el,0,'Bb');
        faceId = dofId(obj.localFaceIndex{2}(i),3);
        Kub = Kub(:,faceId);
        Kbb = Kbb(faceId,faceId);
        
        % get bubble dof
        duLoc = du(dof);
        dmultLoc = dmult(dofId(i,3));

        % recover bubble dof value
        duBub = -Kbb\(Kub'*duLoc + Dbub'*dmultLoc);

        % compute and accumulate bubble strain contribution
        strainBubble = pagemtimes(Bb(:,faceId,:),duBub);
        solv.domain.state.data.strain((l+1):(l+nG),:) = ...
          solv.domain.state.data.strain((l+1):(l+nG),:) + ...
          reshape(strainBubble,6,nG)';
      end
    end


    function assembleConstraint(obj,dt)

      % reset the jacobian blocks
      obj.setJmu(MortarSide.slave, []);
      obj.setJmu(MortarSide.master, []);
      obj.setJum(MortarSide.slave, []);
      obj.setJum(MortarSide.master, []);

      %if isempty(obj.D)
      computeConstraintMatrices(obj,dt);
      % end

      % recompute slave rhs with static condensation contribution
      obj.domains(2).rhs{1} = getPhysicsSolver(obj.domains(2),"Poromechanics").computeRhs();

      [rhsMaster,rhsSlave,rhsMult] = computeRhs(obj);
      %
      addRhs(obj,MortarSide.master,rhsMaster);
      addRhs(obj,MortarSide.slave,rhsSlave);
      % 
      obj.rhsConstraint = rhsMult;

    end


    function computeConstraintMatrices(obj,dt)

      % This method computes the cross grid mortar matrices between
      % connected interfaces

      meshMaster = getMesh(obj,MortarSide.master);
      meshSlave = getMesh(obj,MortarSide.slave);

      dofmMaster = getDoFManager(obj,MortarSide.master);
      dofmSlave =  getDoFManager(obj,MortarSide.slave);

      % get number of index entries for sparse matrices
      % overestimate number of sparse indices assuming all quadrilaterals
      nNmaster = meshMaster.surfaceNumVerts'*obj.interfMesh.elemConnectivity;

      ncomp = 3;

      N1 = sum(nNmaster);
      N2 = sum(meshSlave.surfaceNumVerts);

      nm = (ncomp^2)*N1;
      ns = ncomp^2*N2;

      nEl = obj.getMesh(MortarSide.slave).nSurfaces;
      meshSlave3D = obj.domains(2).grid.topology;
      cellsSlave = obj.interfMesh.f2c{2};
      nus = sum(9*(meshSlave3D.cellNumVerts(cellsSlave)).^2);
      nul = 9*sum(meshSlave3D.cellNumVerts(cellsSlave));
      nl = 9*nEl;

      nDofMaster = dofmMaster.getNumbDoF(obj.coupledVariables);
      nDofSlave = dofmSlave.getNumbDoF(obj.coupledVariables);
      nDofMult = obj.nMult;

      fldM = dofmMaster.getVariableId(obj.coupledVariables);
      fldS = dofmSlave.getVariableId(obj.coupledVariables);

      asbM = assembler(nm,nDofMult,nDofMaster);
      asbD = assembler(ns,nDofMult,nDofSlave);
      asbKs = assembler(nus,nDofSlave,nDofSlave);
      asbLag = assembler(nl,nDofMult,nDofMult);
      asbDcond = assembler(nul,nDofMult,nDofSlave);

      % poromechanics object on the slave side
      poroSlave = obj.domains(2).getPhysicsSolver("Poromechanics");

      obj.Dbubble = zeros(nEl,1);

      % keep track if slave element changed
      isOld = 1;

      f = @(a,b) pagemtimes(a,'ctranspose',b,'none');

      Atb = zeros(3,3);

      % loop over pairs of connected master/slave elements
      for iPair = 1:obj.quadrature.numbInterfacePairs

        is = obj.quadrature.interfacePairs(iPair,1);
        im = obj.quadrature.interfacePairs(iPair,2);

        if isOld ~= is

          % assemble bubble related contribution
          cellId = obj.interfMesh.f2c{2}(isOld);
          [dofRow,dofCol,Kub,Kbb] = computeLocalStiffBubble(poroSlave,cellId,dt);
          faceId = dofId(obj.localFaceIndex{2}(is),3);
          Kub = Kub(:,faceId);
          Kbb = Kbb(faceId,faceId);
          invKbb = inv(Kbb);

          Atb = R'*Atb;

          % assemble local stabilization contribution
          asbKs.localAssembly(dofRow,dofCol,-Kub*invKbb*Kub');
          asbLag.localAssembly(dofMult,dofMult,-Atb*invKbb*Atb');
          asbDcond.localAssembly(dofMult,dofRow,-Atb*invKbb*Kub');
          obj.Dbubble(isOld) = Atb(1,1);

          % refresh local Db matrix
          Atb = zeros(3,3);

          isOld = is;

        end

        dofMult = getMultiplierDoF(obj,is);

        % get dofs
        nodeSlave = obj.interfMesh.local2glob{2}(meshSlave.surfaces(is,:));
        usDof = dofmSlave.getLocalDoF(fldS,nodeSlave);
        nodeMaster = obj.interfMesh.local2glob{1}(meshMaster.surfaces(im,:));
        umDof = dofmMaster.getLocalDoF(fldM,nodeMaster);

        % retrieve mortar integration data
        % gp reference coordinate on master element
        xiMaster = obj.quadrature.getMasterGPCoords(iPair);
        % gp reference coordinate on slave element
        xiSlave = obj.quadrature.getSlaveGPCoords(iPair);
        % determinant of jacobian for each gp
        dJw = obj.quadrature.getIntegrationWeights(iPair);

        % compute slave/master/multipliers basis functions
        [Ns,Nm,Nmult,Nbslave] = ...
          getMortarBasisFunctions(obj.quadrature,im,is,xiMaster,xiSlave);

        [Ns,Nm,Nmult,Nbslave] = ...
          reshapeBasisFunctions(ncomp,Ns,Nm,Nmult,Nbslave);

        % rotation of multipliers
        R = getRotationMatrix(obj.interfMesh,is);
        

        % mortar coupling normal component
        Atm =  MortarQuadrature.integrate(f,Nmult,Nm,dJw);
        Ats =  MortarQuadrature.integrate(f,Nmult,Ns,dJw);

        Atm = R'*Atm;
        Ats = R'*Ats;
        Dbtmp = MortarQuadrature.integrate(f,Nmult,Nbslave,dJw);
        Atb = Atb - Dbtmp;
        %Atb = Atb + Dbtmp;

        % assemble the local mortar matrix contribution
        asbM.localAssembly(dofMult,umDof,Atm);
        asbD.localAssembly(dofMult,usDof,-Ats);


      end

      obj.M = asbM.sparseAssembly();
      obj.D = asbD.sparseAssembly() + asbDcond.sparseAssembly();

      obj.addJum(MortarSide.master, obj.M');
      obj.addJum(MortarSide.slave, obj.D');
      obj.addJmu(MortarSide.master, obj.M);
      obj.addJmu(MortarSide.slave, obj.D);

      obj.Jconstraint = asbLag.sparseAssembly();

      poroSlave.domain.J{fldS,fldS} = poroSlave.K + ...
        asbKs.sparseAssembly();

      % verify partition of unity 
      %pu = sum([obj.M obj.D],2);
      %ssert(norm(pu)<1e-7,'Partiition of unity violated');

    end


    function [rhsMaster,rhsSlave,rhsMult] = computeRhs(obj)

      % retrieve variable field arrays
      varMaster = getVariableField(obj,MortarSide.master);
      varSlave = getVariableField(obj,MortarSide.slave);

      Jmu_master = getJmu(obj,MortarSide.master);
      Jmu_slave = getJmu(obj,MortarSide.slave);
      Jum_master = getJum(obj,MortarSide.master);
      Jum_slave = getJum(obj,MortarSide.slave);

      % retrieve active multipliers
      actMult = getMultiplierDoF(obj);
      iniMult = obj.state.iniMultipliers(actMult);
      mult = obj.state.multipliers(actMult);

      % compute rhs terms
      rhsMaster = Jum_master * (mult - iniMult);
      rhsSlave = Jum_slave * (mult - iniMult);

      rhsMult = Jmu_master * varMaster + Jmu_slave * varSlave + ...
        + obj.Jconstraint * (mult - iniMult);

    end


  end


  methods (Access = protected)

      function setBubbleFaces(obj)
      % associate the slave element to the local face index where bubble is
      % located
      
      nElM = obj.domains(1).grid.topology.nSurfaces;
      nElS = obj.domains(2).grid.topology.nSurfaces;
      obj.localFaceIndex = {zeros(nElM,1),...
        zeros(nElS,1)};
      faceNodes = [
        1 2 3 4;
        1 4 5 8;
        1 2 5 6;
        2 3 6 7;
        3 4 7 8;
        5 6 7 8
        ];

      for i = 1:2
        msh = obj.domains(i).grid.topology;
        mapf2c = obj.interfMesh.f2c{i};
        for f = 1:obj.interfMesh.msh(i).nSurfaces
          % get global nodes of the cell
          nodeC = msh.cells(mapf2c(f),:);
          nodeF = obj.interfMesh.msh(i).surfaces(f,:);
          nodeF = obj.interfMesh.local2glob{i}(nodeF);
          cellNodes = sort(find(ismember(nodeC,nodeF)));
          %
          % Face ID lookup using row-wise comparison
          faceID = find(all(bsxfun(@eq, faceNodes, cellNodes), 2), 1);
          obj.localFaceIndex{i}(f) = faceID;
        end
      end
    end


  end

  methods (Static)

    function var = getCoupledVariables()
      var = Poromechanics.getField();
    end

  end
end



