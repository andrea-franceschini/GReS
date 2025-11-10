classdef MeshGlueBubbleStabilization < MeshGlue
  % Mesh glue class implementing bubble stabilization 

  properties
    localFaceIndex % slave faces with local index in neighboring cell
    Dbubble
  end

  methods (Access = public)
    function obj = MeshGlueBubbleStabilization(id,inputStruct,domains)
      obj@MeshGlue(id,inputStruct,domains);
      assert(strcmp(obj.physic,'Poromechanics'),['' ...
        'Bubble stabilization is currently available for Poromechanics only'])
      setFaces(obj);
    end


    function computeMortarMatricesTest(obj,dt)
      % assumption: each element has only one bubble face
      % loop over slave faces and:
      % 1) compute Aub, Abu and Abb local matrix from the neighbor cell
      % 2) compute local M, D and Db
      % 3) assemble static condensation blocks and Jmult

      % differently from the base method of the mortar class, here global
      % dof indexing is used for the assembled matrices

      % get number of index entries for sparse matrix
      nm = nnz(obj.mesh.elemConnectivity)*9*obj.mesh.nN(1);
      nEl = obj.mesh.msh(2).nSurfaces;
      ns = nEl*9*obj.mesh.nN(2);
      mshSlave = obj.solvers(2).grid.topology;
      cellsSlave = obj.mesh.f2c{2};
      nus = sum(9*(mshSlave.cellNumVerts(cellsSlave)).^2);
      nul = 9*sum(mshSlave.cellNumVerts(cellsSlave));
      nl = 9*nEl;

      % get number of dofs for each block
      nDofMaster = obj.dofm(1).getNumDoF(obj.physic);
      nDofSlave = obj.dofm(2).getNumDoF(obj.physic);
      nDofMult = nEl*obj.dofm(2).getDoFperEnt(obj.physic);

      % define matrix assemblers
      locM = @(imult,imaster,Nmult,Nmaster,dJw) ...
        computeLocMortar(obj,1,imult,imaster,Nmult,Nmaster,dJw);
      locD = @(imult,islave,Nmult,Nslave,dJw) ...
        computeLocMortar(obj,2,imult,islave,Nmult,Nslave,dJw);


      asbM = assembler(nm,nDofMult,nDofMaster,locM);
      asbD = assembler(ns,nDofMult,nDofSlave,locD);
      asbKs = assembler(nus,nDofSlave,nDofSlave);
      asbLag = assembler(nl,nDofMult,nDofMult);
      asbDcond = assembler(nul,nDofMult,nDofSlave);

      poroSlave = getSolver(obj.solvers(2),obj.physic);

      obj.Dbubble = zeros(nEl,1);

      is_old = 0;

      for iPair = 1:obj.quadrature.numbMortarPairs

        is = obj.quadrature.mortarPairs(iPair,1);
        im = obj.quadrature.mortarPairs(iPair,2);

        if iPair == 1
          Dbloc = zeros(3,3);
        elseif is_old ~= is
          % assemble bubble related contributions
          % bubble stabilization on the slave side
          cellId = obj.mesh.f2c{2}(is_old);
          [dofRow,dofCol,Kub,Kbb] = computeLocalStiffBubble(poroSlave,cellId,dt);
          faceId = dofId(obj.localFaceIndex{2}(is),3);
          Kub = Kub(:,faceId);
          Kbb = Kbb(faceId,faceId);
          invKbb = inv(Kbb);

          % assemble local stabilization contribution
          asbKs.localAssembly(dofRow,dofCol,-Kub*invKbb*Kub');
          asbLag.localAssembly(dofMult,dofMult,-Dbloc*invKbb*Dbloc');
          asbDcond.localAssembly(dofMult,dofRow,-Dbloc*invKbb*Kub');
          obj.Dbubble(is_old) = Dbloc(1,1);

          % refresh local Db matrix
          Dbloc = zeros(3,3);
        end

        dofMult = dofId(is,3);

        % retrieve mortar integration data
        xiMaster = obj.quadrature.getMasterGPCoords(iPair);
        xiSlave = obj.quadrature.getSlaveGPCoords(iPair);
        dJw = obj.quadrature.getIntegrationWeights(iPair);

        [Nslave,Nmaster,Nmult,Nbslave] = ...
          getMortarBasisFunctions(obj,im,is,xiMaster,xiSlave);

        [Nslave,Nmaster,Nmult,Nbslave] = ...
          obj.reshapeBasisFunctions(3,Nslave,Nmaster,Nmult,Nbslave);

        asbM.localAssembly(is,im,-Nmult,Nmaster,dJw);
        asbD.localAssembly(is,is,Nmult,Nslave,dJw);

        Dbtmp = pagemtimes(Nmult,'ctranspose',Nbslave,'none');
        Dbtmp = Dbtmp.*reshape(dJw,1,1,[]);
        Dbtmp = sum(Dbtmp,3);
        Dbloc = Dbloc + Dbtmp;

        is_old = is;

      end

      obj.Jmaster = asbM.sparseAssembly();
      obj.Jslave = asbD.sparseAssembly() + asbDcond.sparseAssembly();
      obj.Jmult = asbLag.sparseAssembly();
      poroSlave.J = poroSlave.J + asbKs.sparseAssembly();

    end

    %     function computeMortarMatrices(obj,dt)
    %       % assumption: each element has only one bubble face
    %       if isempty(obj.Jmaster{:})
    %
    %         % loop over slave faces and:
    %         % 1) compute Aub, Abu and Abb local matrix from the neighbor cell
    %         % 2) compute local M, D and Db
    %         % 3) assemble static condensation blocks and Jmult
    %
    %         % differently from the base method of the mortar class, here global
    %         % dof indexing is used for the assembled matrices
    %
    %         % get number of index entries for sparse matrix
    %         nm = nnz(obj.mesh.elemConnectivity)*9*obj.mesh.nN(1);
    %         ns = obj.mesh.nEl(2)*9*obj.mesh.nN(2);
    %         % allocate condensation block for the inner slave domain
    %         poroSlave = getSolver(obj.solvers(2),obj.physic);
    %         cellsSlave = obj.mesh.f2c{2};
    %         mshSlave = obj.solvers(2).grid.topology;
    %         nuu = sum((mshSlave.nDim^2)*(mshSlave.cellNumVerts(cellsSlave)).^2);
    %         nul = 9*sum(mshSlave.cellNumVerts(cellsSlave));
    %         nll = 9*obj.mesh.nEl(2);
    %
    %         % get number of dofs for each block
    %         nDofMaster = obj.dofm(1).getNumDoF(obj.physic);
    %         nDofSlave = obj.dofm(2).getNumDoF(obj.physic);
    %         nDofMult = obj.mesh.nEl(2)*obj.dofm(2).getDoFperEnt(obj.physic);
    %
    %         masterFlag = false(obj.solvers(1).grid.topology.nSurfaces,1); % flags for inner master condensation
    %
    %         fld = [obj.dofm(1).getFieldId(obj.physic), ...
    %           obj.dofm(2).getFieldId(obj.physic)];
    %
    %         % define assemblers
    %
    %         for i = 1:obj.mesh.nEl(2)
    %           is = obj.mesh.getActiveCells(2,i);
    %           elemSlave = getElem(obj,2,is);
    %           masterElems = find(obj.mesh.elemConnectivity(:,is));
    %           if isempty(masterElems)
    %             continue
    %           end
    %           Lloc = zeros(3,3);
    %           %Compute Mortar Slave quantities
    % %           w = elemSlave.getDerBasisFAndDet(is);
    % %           posGP = getGPointsLocation(elemSlave,is);
    %           nodeSlave = obj.mesh.local2glob{2}(obj.mesh.msh(2).surfaces(is,:));
    % %           Nslave = getBasisFinGPoints(elemSlave); % Get slave basis functions
    % %           NbubbleSlave = getBubbleBasisFinGPoints(elemSlave);
    %           Dbloc = zeros(3,3);
    %           Dloc = zeros(3,3*size(Nslave,2));
    %           for im = masterElems'
    %             nodeMaster = obj.mesh.local2glob{1}(obj.mesh.msh(1).surfaces(im,:));
    %             % get building blocks to compute any basis function
    %             [Nslave,Nmaster,Nmult,Nbub] = ...
    %               getMortarBasisFunctions(obj.quadrature,is,im);
    %
    %             if isempty(Nmaster)
    %               continue
    %             end
    %
    %               % working with 3 dimensional dofs in contact mechanics
    %               Nm = Discretizer.reshapeBasisF(Nm,3);
    %               Ns = Discretizer.reshapeBasisF(Ns,3);
    %               Nbub = Discretizer.reshapeBasisF(Nbub,3);
    % %               Nbm = Discretizer.reshapeBasisF(Nbm,3);
    %               Nmult = Discretizer.reshapeBasisF(Nmult,3);
    %
    %               % compute slave mortar matrix
    %               Dtmp = pagemtimes(Nmult,'transpose',Ns,'none');
    %               Dtmp = Dtmp.*reshape(w(id),1,1,[]);
    %               Dtmp = sum(Dtmp,3);
    %               Dloc = Dloc+Dtmp;
    %
    %               % compute master mortar matrix
    %               Mtmp = pagemtimes(Nmult,'transpose',Nm,'none');
    %               Mtmp = Mtmp.*reshape(w(id),1,1,[]);
    %               Mloc = sum(Mtmp,3);
    %               nm = numel(Mloc);
    %               dofMaster = obj.dofm(1).getLocalDoF(nodeMaster,fld(1));
    %               [jMloc, iMloc] = meshgrid(dofMaster, dofId(i,3));
    %               iMvec(cm+1:cm+nm) = iMloc(:);
    %               jMvec(cm+1:cm+nm) = jMloc(:);
    %               Mvec(cm+1:cm+nm) = Mloc(:);
    %               cm = cm+nm;
    %               % [iMvec,jMvec,Mvec,cm] = Discretizer.computeLocalMatrix( ...
    %               % Mloc,iMvec,jMvec,Mvec,cm,w(id),i,nMaster);
    %
    %               % compute master Bubble matrix
    %               Mbloc = pagemtimes(Nmult,'transpose',Nbm,'none');
    %               Mbloc = Mbloc.*reshape(w(id),1,1,[]);
    %               Mbloc = sum(Mbloc,3);
    %
    %               % compute slave Bubble matrix
    %               Dbtmp = pagemtimes(Nmult,'transpose',Nbs,'none');
    %               Dbtmp = Dbtmp.*reshape(w(id),1,1,[]);
    %               Dbtmp = sum(Dbtmp,3);
    %               Dbloc = Dbloc+Dbtmp;
    %
    %               cellId = obj.mesh.f2c{1}(im);
    %               [dofRow,dofCol,Kub,Kbb] = computeLocalStiffBubble(poroMaster,cellId,dt);
    %               % extract only active bubble dofs
    %               faceId = dofId(obj.localFaceIndex{1}(im),3);
    %               Kub = Kub(:,faceId);
    %               Kbb = Kbb(faceId,faceId);
    %               invKbb = inv(Kbb);
    %
    %               if ~masterFlag(im) % if not yet computed
    %                 KuuLoc = -Kub*(invKbb)*Kub';
    %                 [jKMloc,iKMloc] = meshgrid(dofCol,dofRow);
    %                 nkm = numel(KuuLoc);
    %                 iKMvec(ckm+1:ckm+nkm) = iKMloc(:);   jKMvec(ckm+1:ckm+nkm) = jKMloc(:);
    %                 KMvec(ckm+1:ckm+nkm) = KuuLoc(:);
    %                 ckm = ckm + nkm;
    %                 masterFlag(im) = true;
    %               end
    %
    %               % condensation term coupling displacement with multipliers
    %               Bloc = + Mbloc*(invKbb)*Kub';
    %               [jBMloc,iBMloc] = meshgrid(dofRow,dofId(i,3));
    %               nbm = numel(Bloc);
    %               iBMvec(cbm+1:cbm+nbm) = iBMloc(:);   jBMvec(cbm+1:cbm+nbm) = jBMloc(:);
    %               BMvec(cbm+1:cbm+nbm) = Bloc(:);
    %               cbm = cbm+nbm;
    %
    %               % condensation term coupling mortar bubble with multiplier
    %               Lloc = Lloc - Mbloc*(invKbb)*Mbloc';
    %
    %               % sort out gauss points already used
    %               w = w(~id);
    %               posGP = posGP(~id,:);
    %               Nslave = Nslave(~id,:);
    %             else
    %               % pair of elements does not share support. update connectivity
    %               % matrix
    %               obj.mesh.elemConnectivity(im,is) = 0;
    %             end
    %             if all(id)
    %               break
    %             end
    %           end
    %
    %           % Dbloc is diagonal with entries that are all equal. we can cheaply
    %           % store it for future computations in case of non linear material
    %           obj.Dbubble(i) = Dbloc(1,1);
    %
    %           % slave inner condensation contribution
    %           cellId = obj.mesh.f2c{2}(is);
    %           [dofRow,dofCol,Kub,Kbb] = computeLocalStiffBubble(poroSlave,cellId,dt);
    %           [jKloc,iKMloc] = meshgrid(dofCol,dofRow);
    %           % extract only active bubble dofs
    %           faceId = dofId(obj.localFaceIndex{2}(is),3);
    %           Kub = Kub(:,faceId);
    %           Kbb = Kbb(faceId,faceId);
    %           invKbb = inv(Kbb);
    %           KuuLoc = -Kub*(invKbb)*Kub';
    %           nks = numel(KuuLoc);
    %           iKSvec(cks+1:cks+nks) = iKMloc(:);   jKSvec(cks+1:cks+nks) = jKloc(:);
    %           KSvec(cks+1:cks+nks) = KuuLoc(:);
    %
    %           % condensation term coupling slave displacement with multipliers
    %           Bloc = -Dbloc*(invKbb)*Kub';
    %           [jBMloc,iBMloc] = meshgrid(dofRow,dofId(i,3));
    %           nbs = numel(Bloc);
    %           iBSvec(cbs+1:cbs+nbs) = iBMloc(:);   jBSvec(cbs+1:cbs+nbs) = jBMloc(:);
    %           BSvec(cbs+1:cbs+nbs) = Bloc(:);
    %
    %           % Mortar Dslave matrix
    %           dofSlave = obj.dofm(2).getLocalDoF(nodeSlave,fld(2));
    %           [jDloc,iDloc] = meshgrid(dofSlave,dofId(i,3));
    %           nd = numel(Dloc);
    %           iDvec(cd+1:cd+nd) = iDloc(:);   jDvec(cd+1:cd+nd) = jDloc(:);
    %           Dvec(cd+1:cd+nd) = Dloc(:);
    %
    %           % multipliers condensation (includes also master contribution)
    %           Lloc = -Dbloc*(invKbb)*Dbloc'; %+ Lloc;
    %           [jLloc,iLloc] = meshgrid(dofId(i,3),dofId(i,3));
    %           nl = numel(Lloc);
    %           iLvec(cl+1:cl+nl) = iLloc(:);   jLvec(cl+1:cl+nl) = jLloc(:);
    %           Lvec(cl+1:cl+nl) = Lloc(:);
    %
    %           % update counters
    %           cd = cd+nd;
    %           cks = cks+nks;
    %           cl = cl+nl;
    %           cbs = cbs+nbs;
    %
    %           % track element not fully projected
    %           if ~all(id)
    %             error('%i GP not sorted for slave elem numb %i \n',sum(id),is);
    %           end
    %         end
    %
    %         % cut master index vectors for sparse matrix assembly
    %         iMvec = iMvec(1:cm); jMvec = jMvec(1:cm); Mvec = Mvec(1:cm);
    %         iBMvec = iBMvec(1:cbm); jBMvec = jBMvec(1:cbm); BMvec = BMvec(1:cbm);
    %
    %         % assemble mortar matrices in sparse format
    %         obj.Jmaster{1} = sparse(iMvec,jMvec,-Mvec,nDofMult,nDofMaster); %+ ...
    %            sparse(iBMvec,jBMvec,BMvec,nDofMult,nDofMaster);
    %         obj.Jslave{1} = sparse(iDvec,jDvec,Dvec,nDofMult,nDofSlave) + ...
    %           sparse(iBSvec,jBSvec,BSvec,nDofMult,nDofSlave);
    %         obj.Jmult{1} = sparse(iLvec,jLvec,Lvec,nDofMult,nDofMult);
    %         %       % apply condensation to slave inner block
    %         JcondSlave = sparse(iKSvec,jKSvec,KSvec,nDofSlave,nDofSlave);
    %         %JcondMaster = sparse(iKMvec,jKMvec,KMvec,nDofMaster,nDofMaster);
    %         poroSlave.J = poroSlave.J + JcondSlave;
    %         %poroMaster.J = poroMaster.J + JcondMaster;
    %       end
    %     end

    function computeMat(obj,dt)
      computeMortarMatricesTest(obj,dt);
    end

    %
    %     function out = isMatrixComputed(obj)
    %       out = all(cellfun(@(x) ~isempty(x), ...
    %         [obj.Jmaster(:); obj.Jslave(:); obj.Jmult(:)]));
    %     end

    function updateState(obj,dSol)
      % get Poromechanics object
      solv = obj.solvers(2).getSolver(obj.physic);

      % get increment of displacement
      du = solv.state.data.dispCurr - solv.state.data.dispConv;

      % update multipliers
      updateState@MeshGlue(obj,dSol);
      dmult = obj.multipliers(1).curr - obj.multipliers(1).prev;

      % Enhance strain with bubble contribution
      for i = 1:obj.mesh.msh(2).nSurfaces

        % get index of surface and neighbor cell
        el = obj.mesh.f2c{2}(i);

        l = solv.cell2stress(el);

        % get index of cell displacement dofs in du
        dof = getDoFID(solv.mesh,el);

        vtkId = solv.mesh.cellVTKType(el);
        elem = getElement(solv.elements,vtkId);
        nG = elem.GaussPts.nNode;

%         N = getDerBasisF(obj.elements.tetra,el);
%         B = zeros(6,4*obj.mesh.nDim);
%         B(obj.indBtetra(:,2)) = N(obj.indBtetra(:,1));
%         obj.state.data.curr.strain(l+1,:) = (B*dSol(dof))';
%         l = l + 1;

        % get local Db matrix
        Db = diag(repelem(obj.Dbubble(i),3));
        % the last input is the time increment. consider adding it as a
        % property of State.m
        [~,~,Kub,Kbb,Bb] = computeLocalStiffBubble(solv,el,0,'Bb');
        faceId = dofId(obj.localFaceIndex{2}(i),3);
        Kub = Kub(:,faceId);
        Kbb = Kbb(faceId,faceId);
        % get bubble dof
        duLoc = du(dof);
        dmultLoc = dmult(dofId(i,3));
        duBub = -Kbb\(Kub'*duLoc + Db'*dmultLoc);
        strainBubble = pagemtimes(Bb(:,faceId,:),duBub);
        solv.state.data.curr.strain((l+1):(l+nG),:) = ...
          solv.state.data.curr.strain((l+1):(l+nG),:) + ...
          reshape(strainBubble,6,nG)';
      end
    end
  end



  methods(Access = private)

    function setFaces(obj)
      % associate the slave element to the local face index where bubble is
      % located
      
      nElM = obj.solvers(1).grid.topology.nSurfaces;
      nElS = obj.solvers(2).grid.topology.nSurfaces;
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
        msh = obj.solvers(i).grid.topology;
        mapf2c = obj.mesh.f2c{i};
        for f = 1:obj.mesh.msh(i).nSurfaces
          % get global nodes of the cell
          nodeC = msh.cells(mapf2c(f),:);
          nodeF = obj.mesh.msh(i).surfaces(f,:);
          nodeF = obj.mesh.local2glob{i}(nodeF);
          cellNodes = sort(find(ismember(nodeC,nodeF)));
          %
          % Face ID lookup using row-wise comparison
          faceID = find(all(bsxfun(@eq, faceNodes, cellNodes), 2), 1);
          obj.localFaceIndex{i}(f) = faceID;
        end
      end
    end

    %     function computeLocStiffBubble
    %     end
    %
    %     function computeLoc
    %     end




  end


end

