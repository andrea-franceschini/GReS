classdef MeshTyingTPFA < InterfaceSolver

  % Finite volume mesh tying between non conforming interfaces
  
  properties
    neigh                     % neighborship map between non-conforming cells [master,slave]
    faceArea                  % area of the intersection polygons
    faceCenter                % centroid of the intersection polygons
    halfTrans                 % Nface x 2 halfTransmissibilities at the intersecting faces [master,slave]
    flowModels = cell(2,1)    % handle to Flow models in master and slave side
    faceNormals               % normals of the intersection polygons
    elements

  end

  methods

    function obj = MeshTyingTPFA(id,domains,inputStruct)
      
      obj@InterfaceSolver(id,domains,inputStruct);

    end

    function registerInterface(obj,input)

      % check validity of the model
      validateInterface(obj)

      % number of interface variables
      obj.nMult = length(obj.faceArea);

      obj.state.interfacePressures = zeros(obj.nMult,1);

      obj.stateOld = obj.state;

    end


    function updateState(obj,du)

      obj.state.multipliers = obj.state.multipliers + du(1:obj.nMult);

    end


    function assembleConstraint(obj,varargin)

      % set up half face transmissibilities and compute rhs fluxes
      computeHalfTrans(obj);

      % reset the jacobian blocks
      obj.setJmu(MortarSide.slave, []);
      obj.setJmu(MortarSide.master, []);
      obj.setJum(MortarSide.slave, []);
      obj.setJum(MortarSide.master, []);

      assembleMatrices(obj); 
      [rhsMaster,rhsSlave,rhsInterf] = computeRhs(obj);

      addRhs(obj,MortarSide.master,rhsMaster);
      addRhs(obj,MortarSide.slave,rhsSlave);
      obj.rhsConstraint = rhsInterf;

    end


    function assembleMatrices(obj)

      % assemble half transmissibilities
      for side = [MortarSide.master,MortarSide.slave]

        s = getSide(side);

        dofm = obj.getDoFManager(side);

        idPressure = dofm.getVariableId(obj.coupledVariables);

        nDoF = dofm.getNumbDoF(obj.coupledVariables);


        Tfu = sparse((1:obj.nMult)',obj.neigh(:,s),obj.halfTrans(:,s),obj.nMult,nDoF);
        Tuu = sparse(obj.neigh(:,s),obj.neigh(:,s),obj.halfTrans(:,s),nDoF,nDoF);

        addJmu(obj,side,-Tfu);
        addJum(obj,side,-Tfu');
        obj.flowModels{s}.J{idPressure,idPressure} = ...
          obj.flowModels{s}.J{idPressure,idPressure} + Tuu;
        
      end

      obj.Jconstraint = spdiags(sum(obj.halfTrans,2),0,obj.nMult,obj.nMult);

    end

    function [rhsM,rhsS,rhsI] = computeRhs(obj)

      pM = obj.domains(1).state.data.(obj.coupledVariables);
      pS = obj.domains(2).state.data.(obj.coupledVariables);
      pI = obj.state.interfacePressure;
      nDoFM = getDoFManager(obj,MortarSide.master).getNumbDoF(obj.coupledVariables);
      nDoFS = getDoFManager(obj,MortarSide.slave).getNumbDoF(obj.coupledVariables);
      rhsM = zeros(nDoFM,1);
      rhsS = zeros(nDoFS,1);
      rhsM(obj.neigh(:,1)) = -obj.halfTrans(:,1).*pI;
      rhsS(obj.neigh(:,1)) = -obj.halfTrans(:,2).*pI;
      rhsI = - obj.halfTrans(:,1).*pM(obj.neigh(:,1))  ...
             - obj.halfTrans(:,2).*pS(obj.neigh(:,2)) ...
             + sum(obj.halfTrans,2).*pI;

    end


    function [surfaceStr,pointStr] = writeVTK(obj,fac,varargin)

      pCurr = obj.state.faacePressure;
      pOld = obj.stateOld.facePressure;
      p = fac*pCurr + (1-fac)*pOld;

      surfaceStr = [];
      pointStr = [];

      % append state variable to output structure
      surfaceStr(1).name = 'interfacePressure';
      surfaceStr(1).data = p;
    end


    function writeMatFile(obj,fac,tID)

      pCurr = obj.state.interfacePressure;
      pOld = obj.stateOld.multipliers;
      p = fac*pCurr + (1-fac)*pOld;

      obj.outstate.matFile(tID).interfacePressure = p;

    end


  end


  methods (Access = protected)

    function setMortarInterface(obj,varargin)

      interfaceInput = varargin{1};

      % setup the face neighborship at the interface
      masterSurf = getXMLData(interfaceInput.Master,[],"surfaceTag");
      slaveSurf = getXMLData(interfaceInput.Slave,[],"surfaceTag");

      % rough screen connectivity
      obj.interfMesh = InterfaceMesh(obj.domains(1).grid.topology,...
        obj.domains(2).grid.topology,...
        masterSurf,slaveSurf);

      checkInterfaceDisjoint(obj);

      % initialize the maps to store mortar quadrature infos
      nConnections = nnz(obj.interfMesh.elemConnectivity);

      obj.neigh = zeros(nConnections,2);
      obj.faceArea = zeros(nConnections,1);
      obj.faceCenter = zeros(nConnections,3);

      obj.elements = [Elements(obj.interfMesh.msh(1),3),...
                      Elements(obj.interfMesh.msh(2),3)];

      idx = 0;

      for is = 1:obj.interfMesh.msh(2).nSurfaces

        imList = find(obj.interfMesh.elemConnectivity(:,is));

        for j = 1:numel(imList)

          % add new pair to neighborship
          im = imList(j);
          
          idx = processFacePair(obj,is,im,idx);

        end
      end

      obj.neigh = obj.neigh(1:idx,:);
      obj.faceArea = obj.faceArea(1:idx);
      obj.faceCenter = obj.faceCenter(1:idx,:);


      % remove inactive interface elements
      mshMaster = getMesh(obj,MortarSide.master);
      inactiveMaster = ~ismember(1:mshMaster.nSurfaces,...
        obj.neigh(:,1));
      [~, ~, obj.neigh(:,1)] = unique(obj.neigh(:,1));
      mshSlave = getMesh(obj,MortarSide.slave);
      inactiveSlave = ~ismember(1:mshSlave.nSurfaces,...
        obj.neigh(:,2));
      [~, ~, obj.neigh(:,2)] = unique(obj.neigh(:,2));

      % remove master elements
      obj.interfMesh.removeMortarSurface(1,inactiveMaster);
      % remove slave elements
      obj.interfMesh.removeMortarSurface(2,inactiveSlave);

      % map interface surfaces to adjacent domain cells
      mesh3DMaster = obj.domains(1).grid.topology;
      mesh3DSlave = obj.domains(2).grid.topology;
      obj.interfMesh.buildFace2CellMap([mesh3DMaster,mesh3DSlave]);

      obj.neigh = [obj.interfMesh.f2c{1}(obj.neigh(:,1)),...
                   obj.interfMesh.f2c{2}(obj.neigh(:,2))];

    end



    function idx = processFacePair(obj,idS,idM,idx)

      mshS = getMesh(obj,MortarSide.slave);
      mshM = getMesh(obj,MortarSide.master);

      elemS = obj.elements(2).getElement(mshS.surfaceVTKType(idS));
      elemM = obj.elements(1).getElement(mshM.surfaceVTKType(idM));

      % compute auxiliairy plane on slave surface 
      P0 = mshS.surfaceCentroid(idS,:);
      nP = elemS.computeNormal(idS,elemS.centroid);

      coordS3D = FEM.getElementCoords(elemS,idS);
      coordS = pointToSurfaceProjection(P0,nP,coordS3D);

      coordM3D = FEM.getElementCoords(elemM,idM);
      coordM = pointToSurfaceProjection(P0,nP,coordM3D);

      [polyClip,isClipValid] = mxPolygonClip(coordS,coordM);

      if ~isClipValid || isempty(polyClip)
        return
      end

      polyClip = orderPointsCCW2D(polyClip);
      center2D = mxComputePolygonCentroid2D(polyClip);
      area = computePolygonArea2D(polyClip);

      % map the center back to 3D coords
      R = mxComputeRotationMat(nP);
      center3D = P0' + center2D(1)*R(:,2) + center2D(2)*R(:,3);

      % add new face to list
      idx = idx+1;
      obj.neigh(idx,:) = [idM,idS];
      obj.faceNormals(idx,:) = nP;
      obj.faceArea(idx) = area;
      obj.faceCenter(idx,:) = center3D;

    end

    function computeHalfTrans(obj)

      % half face transmissibility 
      % the property is a Nfaces x 2 matrix with the half face
      % transmissibility of each neigboring cell
      r = [1, 1, 1, 2, 2, 2, 3, 3, 3];
      c = [1, 2, 3, 1, 2, 3, 1, 2, 3];

      obj.halfTrans = zeros(length(obj.faceArea),2);

      
      % sign to ensure that the scalar product between c and normal is positive
      % remember: normal points outward the slave domain
      sgn = [-1,+1];

      for side = [MortarSide.master,MortarSide.slave]

        hT = zeros(size(obj.halfTrans,1),1);

        s = getSide(side);
        mesh = obj.domains(s).grid.topology;
        L = obj.faceCenter - mesh.cellCentroid(obj.neigh(:,s),:);
        L = sgn(s)*L;
        KMat = zeros(mesh.nCellTag,9);
        for i=1:mesh.nCellTag
          KMat(i,:) = obj.domains(s).materials.getMaterial(i).PorousRock.getPermVector();
        end

        for k=1:length(r)
          hT = hT + L(:,r(k)) .* KMat(mesh.cellTag(obj.neigh(:,s)),k) .* obj.faceNormals(:,c(k));
        end
        
        obj.halfTrans(:,s) = hT./sum(L.*L,2);

      end

    end 




  end

  methods (Access=private)

    function  validateInterface(obj)

      % validate the model and populate the obj.flowModels pointers

      var = obj.coupledVariables;
      assert(strcmp(obj.coupledVariables,"pressure"),...
        "MeshTyingTPFA can only couple pressure field with TPFA");

      for side = [MortarSide.master,MortarSide.slave]

        s = getSide(side);
        fieldLocation = getDoFManager(obj,side).getFieldLocation(var);
        assert(fieldLocation == entityField.cell,...
          "MeshTyingTPFA works only with cell pressure variables in master and slave domains")

        solvNames = obj.domains(s).solverNames;

        % not elegant but effective
        isFlow = contains(solvNames,"Flow");
        assert(sum(isFlow)==1,"Only one flow model can be coupled");

        obj.flowModels{s} = getPhysicsSolver(obj.domains(s),solvNames(isFlow)); 

      end

    end
  end

  methods (Static)

    function var = getCoupledVariables()
      var = [];
    end

  end
end



