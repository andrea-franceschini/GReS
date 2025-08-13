classdef ContactMortar < Mortar

  properties
    activeSet         % struct with stick, slip, open node id
    gap               % struct with global, normal and tangential gap
    Jmu               % master displacement block
    Jmt               % master multipliers block
    Jsu               % slave displacement block
    Jst               % slave multipliers block
    rhsUm             % master displacement rhs
    rhsUs             % slave displacement rhs
    rhsT              % traction
    traction          % global traction array
    iniTraction     
    physic = "Poromechanics"
    multiplierType = "P0"
  end
  
  methods
    function obj = ContactMortar(id,inputStruct,domains)

      obj@Mortar(inputStruct,domains);

      isFldMaster = isField(obj.dofm(1),obj.physic);
      isFldSlave = isField(obj.dofm(2),obj.physic);
      assert(isFldMaster,['%s not available for ' ...
        'master domain %i'],obj.physic,obj.idDomain(1));
      assert(isFldSlave,['%s not available for ' ...
        'slave domain %i'],obj.physic,obj.idDomain(2));

      %computing mortar matrices D and M and updating list of slave entitities
      computeMortarMatrices(obj);

      %remove inactive multipliers from D and M
      id = find(~any(obj.D,2));
      obj.D(id,:) = [];
      obj.M(id,:) = [];

      % initialize Jacobian and rhs for the interface (now that active
      % multipliers are known)
      initializeInterface(obj)

      finalizeInterface(obj.mesh,obj.solvers);
      
      % use the updated mesh object for the output of the interface
      if ~isempty(obj.outStruct)
        obj.outStruct.VTK = VTKOutput(obj.mesh.msh(2),obj.outStruct.name);
      end
    end

    function computeContactMatrices(obj)

      um = obj.solvers(1).state.data.dispCurr;
      us = obj.solvers(2).state.data.dispCurr;
      um_old = obj.solvers(1).state.data.dispConv;
      us_old = obj.solvers(2).state.data.dispConv;

      um_slip = um - um_old;
      us_slip = us - us_old;

      % Compute all contact related quantities
      % Everything is done in one single loop
      ncomp = obj.dofm(2).getDoFperEnt(obj.physic);

      % define matrix assemblers
      [asbMu,asbDu,absMt,absDt] = defineAssemblers(obj);

      for is = 1:obj.mesh.msh(2).nSurfaces
        masterElems = find(obj.mesh.elemConnectivity(:,is));
        if isempty(masterElems)
          continue
        end

        Rloc = getRotationMatrix(obj,is);

        elSlave = getElem(obj,2,is);
        nN = elSlave.nNode;
        switch obj.multiplierType
          case 'P0'
            Dloc = zeros(ncomp,ncomp*nN);
          otherwise
            Dloc = zeros(ncomp*nN,ncomp*nN);
        end

        for im = masterElems'

           % compute basis function matrices
          [Nslave,Nmaster,Nmult] = ...
            getMortarBasisFunctions(obj.quadrature,is,im);

          if isempty(Nmaster)
            % refine connectivity matrix
            %obj.mesh.elemConnectivity(im,is) = 0;
            continue
          end

          % reshape basis function matrices to match number of components
          [Nslave,Nmaster,Nmult] = ...
            obj.reshapeBasisFunctions(ncomp,Nslave,Nmaster,Nmult);

          % apply rotation matrix to multiplier bf matrices
          Nmult = pagemtimes(Rloc,Nmult);

          Mloc =  obj.quadrature.integrate(@(a,b) ...
            pagemtimes(a,'ctranspose',b,'none'), Nmaster,Nmult);
          asbMu.localAssembly(is,im,Mloc);

          Dloc = Dloc + ...
            obj.quadrature.integrate(@(a,b) pagemtimes(a,'ctranspose',b,'none'),...
            Nmult,Nslave);

          obj.rhsUm = Mloc*traction(get)
        end

        % if Dloc is empty, the current slave element is inactive remove
        % also corresponding master elements if it is not connected to any
        % other element
        if all(Dloc==0,"all")
          isInactiveSlave(is) = true;
        end

        asbDu.localAssembly(is,is,Dloc);
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
    
    function initializeInterface(obj)
      obj.traction = struct('prev',[],'curr',[]);
      obj.activeSet = struct('stick',[],'slip',[],'open',[],'tag',[]);
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
  end

  methods (Access = protected)

    function [dofRow,dofCol,mat] = computeLoc(obj,i1,i2,kernel,fl)

      if isnumeric(kernel)
        % local matrix is already provided
        mat = kernel;
      else
        % Compute local matrix using kernel as anonymous function
        mat = obj.quadrature.integrate(kernel);
      end

      dofMult = getMultiplierDoF(obj,i1);
      switch fl
        case 1  % Mu
          dofRow = obj.mesh.local2glob{1}(obj.mesh.msh(1).surfaces(i2,:));
          dofCol = dofMult;
        case 2  % Du
          dofRow = obj.mesh.local2glob{2}(obj.mesh.msh(2).surfaces(i2,:));
          dofCol = dofMult;
        case 3  % Mt
          dofRow = dofMult;
          dofCol = obj.mesh.local2glob{1}(obj.mesh.msh(1).surfaces(i2,:));
        case 4  % Dt
          dofRow = dofMult;
          dofCol = obj.mesh.local2glob{2}(obj.mesh.msh(2).surfaces(i2,:));
        case 5  % Dt
          dofRow = dofMult;
          dofCol = dofMult;
      end
    end

    function [asbMu,asbDu,asbMt,asbDt,asbA] = defineAssemblers(obj)

      % get number of index entries for sparse matrices
      % valid only for P0 multipliers!

      % number of master nodes attacched to slave elements
      nNmaster = obj.mesh.msh(1).surfaceNumVerts'*obj.mesh.elemConnectivity;

      N1 = sum(nNmaster);
      N2 = sum(obj.mesh.msh(2).surfaceNumVerts);

      nmu = (ncomp^2)*N1;
      nsu = ncomp^2*N2;

      nmt = nmu;
      nst = nsu;

      na = obj.mesh.mesh(2).nSurfaces*ncomp^2;

      nDofMaster = obj.dofm(1).getNumDoF(obj.physic);
      nDofSlave = obj.dofm(2).getNumDoF(obj.physic);
      nDofMult = getNumbMultipliers(obj);

      % define matrix assemblers
      locMu = @(iRow,iCol,kernel) ...
        computeLoc(obj,iRow,iCol,kernel,1);
      locDu = @(iRow,iCol,kernel) ...
        computeLoc(obj,iRow,iCol,kernel,2);
      locMt = @(iRow,iCol,kernel) ...
        computeLoc(obj,iRow,iCol,kernel,3);
      locDt = @(iRow,iCol,kernel) ...
        computeLoc(obj,iRow,iCol,kernel,4);
      locA = @(iRow,iCol,kernel) ...
        computeLoc(obj,iRow,iCol,kernel,5);

      asbMu = assembler(nmu,locMu,nDofMult,nDofMaster);
      asbDu = assembler(nsu,locDu,nDofMult,nDofSlave);
      asbMt = assembler(nmt,locMt,nDofMult,nDofMaster);
      asbDt = assembler(nst,locDt,nDofMult,nDofSlave);
      asbA = assembler(na,locA,nDofMult,nDofMult);
    end
  end
end

