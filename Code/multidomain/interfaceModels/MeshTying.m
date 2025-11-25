classdef MeshTying < InterfaceSolver

  % Standard mesh tying between non conforming interfaces
  
  properties

    coupledVariable
    state
    stateOld
    D
    M

  end

  methods
    function obj = MeshTying(domains,inputStruct)
      
      obj@InterfaceSolver(domains,inputStruct);

    end

    function setInterface(obj,input)

      % specify the variable to be coupled
      if isfield(input,"variable")
        sharedVars = getXMLData(input,[],"variable");
      else
        varMaster = getVariableNames(obj.domains(1).dofm);
        varSlave = getVariableNames(obj.domains(2).dofm);

        sharedVars = intersect(varSlave,varMaster);

      end

      if ~isscalar(sharedVars)
        error("Multiple variables are shared by the master and slave domain. " + ...
          "The 'variable' attribute is required to define which of the " + ...
          "available variable has to be coupled. If you wish to couple more than one variable field," + ...
          "define a MeshTying interface for each variable.")
      else
        obj.coupledVariable = sharedVars;
      end

      setInterface@InterfaceSolver(obj,input);

    end

    function initializeConstraint(obj)

      ncomp = obj.domains(2).dofm.getNumberOfComponents(obj.coupledVariable);
      obj.nMult = ncomp * getNumberOfEntities(obj.multiplierLocation,...
                                              obj.interfMesh.msh(2));

      obj.state.multipliers = zeros(obj.nMult,1);
      obj.state.iniMultipliers = zeros(obj.nMult,1);

      obj.stateOld = obj.state;

    end

    function advanceState(obj)
      % do nothing here
    end


    function updateState(obj,du)

      obj.multipliers.curr = obj.multipliers.curr + du(1:obj.nMult);

    end


    function assembleConstraint(obj)

      computeConstraintMatrices(obj);

      [rhsM,rhsS,rhsMult] = computeRhs(obj);

      fldMaster = obj.domains(1).getVariableId(obj.coupledVariable);
      fldSlave = obj.domains(2).getVariableId(obj.coupledVariable);

      % write jacobian to domain coupling blocks

      % master side
      obj.domain(1).Jum{obj.interfaceId(1)}{fldMaster} = obj.M';
      obj.domain(1).Jmu{obj.interfaceId(1)}{fldMaster} = obj.M;

      % slave side
      obj.domain(2).Jum{obj.interfaceId(2)}{fldSlave} = obj.D';
      obj.domain(2).Jmu{obj.interfaceId(2)}{fldSlave} = obj.D;

      rhsMaster = getRhs(obj.domain(1),fldMaster);
      rhsSlave = getRhs(obj.domain(2),fldSlave);

      obj.domain(1).rhs{obj.fldMaster} = rhsMaster + rhsM;
      obj.domain(2).rhs{obj.fldSlave} = rhsSlave + rhsS;
      obj.rhs = rhsMult;

    end


    function computeConstraintMatrices(obj)

      % This method computes the cross grid mortar matrices between
      % connected interfaces


      % number of components per dof of interpolated physics
      if ~isempty(obj.domains(2).dofm)
        ncomp = obj.domains(2).dofm.getNumberOfComponents(obj.coupledVariable);
      else
        ncomp = 1;
      end

      % get number of index entries for sparse matrices
      % overestimate number of sparse indices assuming all quadrilaterals
      nNmaster = obj.interfMesh.msh(1).surfaceNumVerts'*obj.interfMesh.elemConnectivity;

      switch obj.multiplierLocation
        case entityField.cell
          N1 = sum(nNmaster);
          N2 = sum(obj.interfMesh.msh(2).surfaceNumVerts);
        otherwise
          N1 = nNmaster*obj.interfMesh.msh(2).surfaceNumVerts;
          N2 = sum(obj.interfMesh.msh(2).surfaceNumVerts.^2);
      end

      nm = (ncomp^2)*N1;
      ns = ncomp^2*N2;

      nDofMaster = obj.dofm(1).getNumbDoF(obj.coupledVariable);
      nDofSlave = obj.dofm(2).getNumbDoF(obj.coupledVariable);
      nDofMult = obj.nMult;


      fldM = obj.domains(1).dofm.getVariableId(obj.coupledVariable);
      fldS = obj.domains(2).dofm.getVariableId(obj.coupledVariable);

      asbM = assembler(nm,nDofMult,nDofMaster);
      asbD = assembler(ns,nDofMult,nDofSlave);

      f = @(a,b) pagemtimes(a,'ctranspose',b,'none');

      % loop over pairs of connected master/slave elements
      for iPair = 1:obj.quadrature.numbMortarPairs

        is = obj.quadrature.mortarPairs(iPair,1);
        im = obj.quadrature.mortarPairs(iPair,2);

        % get dofs 
        nodeSlave = obj.interfMesh.local2glob{2}(obj.interfMesh.msh(2).surfaces(is,:));
        usDof = obj.dofm(2).getLocalDoF(fldS,nodeSlave);
        nodeMaster = obj.interfMesh.local2glob{1}(obj.interfMesh.msh(1).surfaces(im,:));
        umDof = obj.dofm(1).getLocalDoF(fldM,nodeMaster);
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
      nodes = obj.interfMesh.local2glob{side}(obj.interfMesh.msh(side).surfaces(iu,:));
      fld = obj.dofm(side).getFieldId(obj.physic);
      dofc = obj.dofm(side).getLocalDoF(nodes,fld);
      dofr = getMultiplierDoF(obj,imult);
    end


    function [rhsMaster,rhsSlave,rhsMult] = computeRhs(obj)

      actMult = getMultiplierDoF(obj);
      obj.rhsMaster = ...
        obj.Jmaster'*(obj.multipliers.curr(actMult)-obj.iniMultipliers(actMult));
      var = getState(obj.solvers(1).getSolver(obj.physic));
      ents = obj.dofm(1).getActiveEnts(obj.physic);
      obj.rhsMult = obj.rhsMult + obj.Jmaster*var(ents);
      
    end



    % function dofs = getMultiplierDoF(obj,is)
    %   % return multiplier dofs associated with slave element #is inside
    %   % obj.interfMesh.msh(2)
    % 
    %   nc = obj.dofm(2).getDoFperEnt(obj.physic);
    %   if nargin > 1
    %     assert(isscalar(is),'Input id must be a scalar integer \n')
    %     if strcmp(obj.multiplierType,'P0')
    %       dofs = dofId(is,nc);
    %     else
    %       nodes = obj.interfMesh.msh(2).surfaces(is,:);
    %       dofs = dofId(nodes,nc);
    %     end
    %   else
    %     if strcmp(obj.multiplierType,'P0')
    %       dofs = dofId((1:obj.interfMesh.msh(2).nSurfaces)',nc);
    %     else
    %       dofs = dofId((1:obj.interfMesh.msh(2).nNodes)',nc);
    %     end
    %   end
    % end


    % function nDof = getNumbMultipliers(obj)
    %   % nDoFAct -> number of currently active multipliers in the interface
    %   % nDoFTot -> total number of dof in the whole interface
    % 
    %   nc = obj.dofm(2).getDoFperEnt(obj.physic);
    %   switch obj.multiplierType
    %     case 'P0'
    %       nDof = nc*obj.interfMesh.msh(2).nSurfaces;
    %     case {'dual','standard'} % nodal multipliers
    %       % get number of nodes in active cells
    %       nDof = nc*obj.interfMesh.msh(2).nNodes;
    %   end
    % end


    function computeMortarInterpolation(obj)

      processMortarPairs(obj.quadrature); 

      inactiveMaster = ~ismember(1:obj.interfMesh.msh(1).nSurfaces,...
        obj.quadrature.mortarPairs(:,2));

      [~, ~, obj.quadrature.mortarPairs(:,2)] = ...
        unique(obj.quadrature.mortarPairs(:,2));

      inactiveSlave = ~ismember(1:obj.interfMesh.msh(2).nSurfaces,...
        obj.quadrature.mortarPairs(:,1));

      [~, ~, obj.quadrature.mortarPairs(:,1)] = ...
        unique(obj.quadrature.mortarPairs(:,1));

      % remove master elements
      obj.interfMesh.removeMortarSurface(1,inactiveMaster);

      % remove slave elements
      obj.interfMesh.removeMortarSurface(2,inactiveSlave);


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
        out.VTK = VTKOutput(obj.interfMesh.msh(2),out.name);
        out.tList = outState.timeList;
        obj.outStruct = out;
      end
    end

    function checkInterfaceDisjoint(obj)
      % check that the nodes of mortar and slave side are disjoint
      if obj.idDomain(1)==obj.idDomain(2)
        % interface defined within the same domain 
        out = setdiff(obj.interfMesh.local2glob{1},obj.interfMesh.local2glob{2});
        if ~all(ismember(obj.interfMesh.local2glob{1},out))
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

      obj.interfMesh = interfaceMesh(domains(1).grid.topology,domains(2).grid.topology,...
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

      obj.setQuadrature(quadType,nG,nInt);

      removeSlaveBCents(obj);

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
 

    function varargout = reshapeBasisFunctions(nc,varargin)
      assert(numel(varargin)==nargout);
      varargout = cell(1,nargout);
      for i = 1:numel(varargin)
        varargout{i} = reshapeBasisF(varargin{i},nc);
      end
    end


  end
end

