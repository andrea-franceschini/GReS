classdef GrowningDomain < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here

  properties
    % solver       % Grid for each field
    grids        % Grid for each field
    nSolver      % number of the solver related to the grid
    nfields      % Number of fields of SPFLOW
    bound        % store the boundary conditions to remap
  end

  methods(Access = public)
    function obj = GrowningDomain(linSyst,bound)
      %ADAPGRID Construct an instance of this class
      obj.initLinSyst(linSyst);
      obj.initBoundary(bound);
    end

    function [state1,state2]=addCells(obj,linSyst,nsyst,cell,dir,state1,state2)
      %ADDCELL Use a surface as reference, and create cells above this
      % surface
      for i=1:length(cell)
        [state1,state2]=obj.addCell(linSyst,nsyst,cell(i),dir,state1,state2);
      end
      temp=linSyst.solver(1);
      temp.H=[];
      temp.P=[];
      temp.J=[];
    end

    function [state1,state2]=addCell(obj,linSyst,nsyst,cell,dir,state1,state2)
      %ADDCELL Add a cell in the grid and update the mesh

      % Find the position to where to grow the mesh.
      [n_cells,n_ijk] = obj.grids(nsyst).findNeighborhod( ...
        obj.grids(nsyst).ijk(cell,:));
      if obj.canGrow(dir,n_cells)
        % Find the information about the cell to be created.
        ijk_newCell=n_ijk(dir,:);
        ncell = obj.grids(nsyst).newCellData(...
          linSyst.solver(nsyst).mesh, linSyst.solver(nsyst).faces, ...
          ijk_newCell,dir,0.5,cell);

        % Update the mesh
        obj.updateCell(linSyst.solver(nsyst),ijk_newCell,ncell);

        % Update the states.
        [state1,state2]=obj.updateState(cell,state1,state2);

        % Update the IJK information.
        obj.grids(nsyst).updateIJK(ijk_newCell,ncell.gtCell+1, ...
          ncell.gtNode+ncell.nNodes,ncell.gtFace+ncell.nFaces);
      end
    end




    
    function grow(obj,nsyst,cell,dir,domain,state1,state2)
      %ADDCELL Use a surface as reference, and create cells above this
      % surface
      for i=1:length(cell)
        obj.addCellB(domain.solver(nsyst),nsyst,cell(i),dir,state1,state2);
      end
      map = domain.solver(nsyst);
      map.H=[];
      map.P=[];
      map.J=[];
    end

    function addCellB(obj,linSyst,nsyst,cell,dir,state1,state2)
      %ADDCELL Add a cell in the grid and update the mesh

      % Find the position to where to grow the mesh.
      [n_cells,n_ijk] = obj.grids(nsyst).findNeighborhod( ...
        obj.grids(nsyst).ijk(cell,:));
      if obj.canGrow(dir,n_cells)
        % Find the information about the cell to be created.
        ijk_newCell=n_ijk(dir,:);
        ncell = obj.grids(nsyst).newCellData(...
          linSyst.mesh, linSyst.faces, ...
          ijk_newCell,dir,0.5,cell);

        % Update the mesh
        obj.updateCell(linSyst,ijk_newCell,ncell);

        % Update the states.
        obj.updateState(cell,state1,state2);

        % Update the IJK information.
        obj.grids(nsyst).updateIJK(ijk_newCell,ncell.gtCell+1, ...
          ncell.gtNode+ncell.nNodes,ncell.gtFace+ncell.nFaces);
      end
    end




  end

  methods (Access = private)

    % Functions related to initialaze this class.
    function initLinSyst(obj,linSyst)
      obj.grids = containers.Map('KeyType','double','ValueType','any');
      numSol = linSyst.numSolvers;
      countfield = 0;
      for i=1:numSol
        if linSyst(i).dofm.fields.field == "SinglePhaseFlow" ...
            && linSyst(i).dofm.fields.scheme == "FV"
          countfield=countfield+1;
          tmpSolver = linSyst.solver(i);
          obj.nSolver(i) = i;
          obj.grids(countfield) = structuredGrid(tmpSolver.mesh,tmpSolver.faces);
        end
      end
      obj.nfields = countfield;
    end

    function initBoundary(obj,bord)
      bkeys = keys(bord.db);
      obj.bound = containers.Map('KeyType','char','ValueType','any');
      for i=1:length(bkeys)
        it = bord.db(bkeys{i});
        if it.physics == 'SinglePhaseFlow'
          obj.bound(bkeys{i})=GrowingBoundary(it.cond,it.data,'Constant');
        end
      end
    end




    % Functions related to update the physics.
    function updateCell(obj,linSyst,ijk_newCell,ncell)
      %UPDATECELL update the model.

      % Update Geometric information.
      obj.updateFaces(linSyst.faces,ijk_newCell,ncell);
      obj.updateElements(linSyst.elements,ncell);
      obj.updateMesh(linSyst.mesh,ncell,1);

      % Update Physics information.
      trans = linSyst.computeTransCell(ncell);
      obj.updateTrans(linSyst,ncell,trans);

      % Update DOFmanager.
      obj.updateDOFM(linSyst.dofm,ncell);
    end

    function updateFaces(obj,faces,refCell,ncell)
      %UPDATEFACES update the face information to add the new cell

      % Mask to some Variables.
      faces2Add = ncell.nFaces;

      % Find the cells around.
      [neighcells,~] = obj.grids(obj.nfields).findNeighborhod(refCell);

      % Some Important Variables.
      ref_pos = [6 5 1 2       % sth (1)
        3 4 8 7       % nth (2)
        5 8 4 1       % wst (3)
        2 3 7 6       % est (4)
        1 4 3 2       % bot (5)
        6 7 8 5       % top (6)
        ];

      faceCentroid(1:6,1:3)=0;
      faceNormal(1:6,1:3)=0;
      nod2Fac(1:6,1:4)=0;
      facNeig(1:6,1:2)=ncell.gtCell+1;
      for i=1:6
        pos = ref_pos(i,:);
        faceCentTemp=sum(ncell.coord(pos,:))/4;
        triagCent(1:4,1:3)=0;
        triagSurf(1:4)=0;
        for j=1:4
          triagCent(j,:)=sum([ faceCentTemp; ...
            ncell.coord(pos(j),:); ...
            ncell.coord(pos(mod(j,4)+1),:) ])/3;
          triagSurf(j)=norm(0.5*cross((ncell.coord(pos(j),:)-faceCentTemp),...
            (ncell.coord(pos(mod(j,4)+1),:)-faceCentTemp)));
        end
        faceSurf = sum(triagSurf);
        faceCentroid(i,:) = sum(triagSurf'.*triagCent)/faceSurf;

        [~,nodeC]=min(ncell.node(pos));
        nodeL = mod(nodeC - 2, 4) + 1;
        nodeR = mod(nodeC, 4) + 1;
        nodeA = mod(nodeC + 1, 4) + 1;
        if (ncell.node(pos(nodeL))<ncell.node(pos(nodeR)))
          nod2Fac(i,:)=[
            ncell.node(pos(nodeC)) ncell.node(pos(nodeL)) ...
            ncell.node(pos(nodeA)) ncell.node(pos(nodeR))];
          facNeig(i,1)=neighcells(i);   % <---(i,1) or (i,2) to check
        else
          nod2Fac(i,:)=[
            ncell.node(pos(nodeC)) ncell.node(pos(nodeR)) ...
            ncell.node(pos(nodeA)) ncell.node(pos(nodeL))];
          facNeig(i,2)=neighcells(i);   % <---(i,1) or (i,2) to check
        end
        v1 = ncell.coord(pos(nodeL),:)-ncell.coord(pos(nodeC),:);
        v2 = ncell.coord(pos(nodeA),:)-ncell.coord(pos(nodeC),:);
        faceNormal(i,:)=cross(v1,v2)/norm(cross(v1,v2)); % Normalize the vector
        faceNormal(i,:)=faceSurf*faceNormal(i,:); % storing the face surface in the vector
      end
      % Variables to add with simple construction.
      mapN2F(1:faces2Add)=max(faces.mapN2F)+4*(1:faces2Add);
      faces2Elements=[ncell.face' (1:6)'];
      mapF2E=max(faces.mapF2E)+6;

      % UPDATE THE FACE STRUCT.
      fNewFaces = neighcells==0;
      nodes2Faces = nod2Fac(fNewFaces,:)';
      faces.nodes2Faces = [faces.nodes2Faces; nodes2Faces(:)];
      faces.mapN2F = [faces.mapN2F; mapN2F'];
      faces.faces2Elements = [faces.faces2Elements; faces2Elements];
      faces.mapF2E = [faces.mapF2E; mapF2E];
      faces.faceNeighbors = [faces.faceNeighbors; facNeig(fNewFaces,:)];
      faces.faceCentroid = [faces.faceCentroid; faceCentroid(fNewFaces,:)];
      faces.faceNormal = [faces.faceNormal; faceNormal(fNewFaces,:)];
      faces.nFaces = faces.nFaces+faces2Add;

      facesNot2Add = 6 - faces2Add;
      existFaces = faces2Elements(~fNewFaces,1);
      for i=1:facesNot2Add
        loc=faces.faceNeighbors(existFaces(i),:)==0;
        faces.faceNeighbors(existFaces(i),loc)=ncell.gtCell+1;
      end
    end




  end


  methods(Static, Access = private)
    % Functions to check the mesh
    function flag = canGrow(dir,n_cells)
      %CANGROW Function to check if a cell can be add in the mesh, for
      %the direction specified.
      flag = n_cells(dir)==0;
    end

    % Functions related to update the mesh
    function updateMesh(mesh,ncell,cellTag)
      %UPDATEMESH update the face information to add the new cell
      coord(1:ncell.nNodes,1:3)=0;
      for i=1:ncell.nNodes
        nloc=ncell.node==ncell.gtNode+i;
        coord(i,:)=ncell.coord(nloc,:);
      end

      % Here is necessary to update the boundary.
      nodefaces = [6 5 1 2       % sth (1)
        3 4 8 7       % nth (2)
        5 8 4 1       % wst (3)
        2 3 7 6       % est (4)
        1 4 3 2       % bot (5)
        6 7 8 5       % top (6)
        ];
      refFaceLc = ncell.faceGrow+(-1).^(mod(ncell.faceGrow,2)+1);
      refFaceGL = ncell.refBound(refFaceLc);
      for i=1:6
        if (ncell.face(i)>ncell.gtFace)
          % Add boundary information if a face is created.
          mesh.nSurfaces = mesh.nSurfaces+1;
          mesh.surfaces(end+1,:)=ncell.node(nodefaces(i,:));
          mesh.surfaceTag(end+1,:)=mesh.surfaceTag(refFaceGL);
          mesh.surfaceNumVerts(end+1,:)=4;
          mesh.surfaceVTKType(end+1,:)=9;

          % % TODO -> acrescentar essas duas informacoes no computo da celula
          % mesh.surfaceArea(end+1,:) = ncell.surfaceArea(i);
          % mesh.surfaceCentroid(end+1,:) = ncell.surfaceCentroid(i,:);

          % Correcting in case a boundary in the different tag group.
          if (ncell.refBound(i)~=0)
            mesh.surfaceTag(end,:)=mesh.surfaceTag(ncell.refBound(i));
          end
        else
          % Delete boundary information if a face is already exist's.
          nodes = ncell.node(nodefaces(i,:));
          locs=sum(ismember(mesh.surfaces,nodes),2);
          loc=find(locs==4);
          mesh.nSurfaces = mesh.nSurfaces-1;
          mesh.surfaces(loc,:)=[];
          mesh.surfaceTag(loc,:)=[];
          mesh.surfaceNumVerts(loc,:)=[];
          mesh.surfaceVTKType(loc,:)=[];
        end
      end

      % UPDATE THE GRID STRUCT.
      mesh.nNodes = mesh.nNodes+ncell.nNodes;
      mesh.nCells = mesh.nCells+1;
      mesh.coordinates = [mesh.coordinates;coord];
      mesh.cells = [mesh.cells;ncell.node];
      mesh.cellTag = [mesh.cellTag;cellTag];
      mesh.cellNumVerts = [mesh.cellNumVerts;8];
      mesh.cellVTKType = [mesh.cellVTKType;12];
    end

    % Functions related to update the element
    function updateElements(grid,ncell)
      %UPDATEFACES update the elements information to add the new cell
      center = mean(ncell.coord);

      tets = [ 1 2 4 5;
        2 3 4 7;
        2 5 6 7;
        4 5 7 8;
        2 4 5 7 ];

      vol = 0;
      for i = 1:size(tets,1)
        A = ncell.coord(tets(i,1),:);
        B = ncell.coord(tets(i,2),:);
        C = ncell.coord(tets(i,3),:);
        D = ncell.coord(tets(i,4),:);

        % Compute volume of tetrahedron
        tet_vol = abs(dot(B - A, cross(C - A, D - A))) / 6;
        vol = vol + tet_vol;
      end

      % UPDATE THE GRID STRUCT.
      grid.nCellsByType(2) = grid.nCellsByType(2) + 1;
      grid.mesh.cellCentroid = [grid.mesh.cellCentroid;center];
      grid.mesh.cellVolume = [grid.mesh.cellVolume; vol];
    end

    % Functions related to update the transmissibility
    function updateTrans(linSys,cell,trans)
      cut = (cell.face>cell.gtFace);
      linSys.trans = [linSys.trans; trans(cut)];
      linSys.isIntFaces = [linSys.isIntFaces; false(cell.nFaces,1)];
      linSys.isIntFaces(cell.face(~cut)) = true;
    end

    % Functions related to update the degree of freedom's
    function updateDOFM(linSys,cell)
      linSys.numEntsField =linSys.numEntsField+1;
      linSys.totDoF = linSys.totDoF+1;
      linSys.numEntsSubdomain = linSys.numEntsSubdomain+1;
      linSys.fields.entCount = linSys.fields.entCount+1;
      linSys.fields.subdomainEnts{1}=[linSys.fields.subdomainEnts{1};cell.gtCell+1];
      linSys.fields.isEntActive=[linSys.fields.isEntActive;true];
      linSys.addCellTag(1);
    end




    % Functions related to the states
    function updateState(cell,state1,state2)
      labels = fieldnames(state1.data);
      for label=1:length(labels)
        data = getfield(state1.data,labels{label});
        state1.data=setfield(state1.data,labels{label},[data; data(cell)]);
        state2.data=setfield(state2.data,labels{label},[data; data(cell)]);
      end
    end




  end
end