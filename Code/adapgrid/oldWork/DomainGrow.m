classdef DomainGrow < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here

  properties(Access = public)
    simData      % Store a link to the simulation data
    sdomain      % grid for each field
  end

  properties(Access = private)
    nsolver      % number of the solver related to the flow.
  end

  methods(Access = public)
    function obj = DomainGrow(data)
      %ADAPGRID Construct an instance of this class

      % Store a pointer to the main domain
      obj.simData = data;

      % Create a skeleton structure to control the domain.
      numSol = data.numSolvers;
      for  i=1:numSol
        acess = data.solver(i);
        if acess.model.isFlow
          if acess.model.isFVTPFABased(acess.field)
            obj.sdomain = SDomain(acess);
            obj.nsolver = i;
            break;
          end
        end
      end
      
    end



    function grow(obj,boundary,state1,state2)
      % XXXX % Updates the domain mesh structure by adding cells outward
      % from the reference surface.

      [cell,direction] = obj.findBorderCell(boundary);

      for i=1:length(cell)
        obj.addCell(cell(i),direction(i),state1,state2);
      end
      lnk = obj.simData.solver(obj.nsolver);
      lnk.H=[];
      lnk.P=[];
      lnk.J=[];
    end

    function addCell(obj,cell,dir,state1,state2)
      %ADDCELL Add a cell in the grid and update the mesh

      scut = obj.simData.solver(obj.nsolver);

      % Find the position to where to grow the mesh.
      [n_cells,n_ijk] = obj.sdomain.findNeighborhodOfCell(cell);
      if obj.canGrow(dir,n_cells)
        % Find the information about the cell to be created.
        ijk_newCell=n_ijk(dir,:);
        ncell = obj.sdomain.newCellData(scut.mesh,scut.faces, ...
          ijk_newCell,dir,0.5,cell);

        % Update the mesh
        obj.updateCell(ncell);

        % Update the reference of the boundary.
        obj.sdomain.updateBorder(scut,ncell,cell);

        % Update the states.
        obj.updateState(cell,state1,state2);

        obj.updateBoundary();

        % Update the IJK information.
        % obj.sdomain.updateIJK(ijk_newCell,ncell.gtCell+1, ...
        %   ncell.gtNode+ncell.nNodes,ncell.gtFace+ncell.nFaces);
      end
    end

    function updateBoundary(obj)
      lnkDm = obj.simData.solver(obj.nsolver);
      keys = string(obj.sdomain.boundaries.keys);
      for i=keys
        lnkBc = obj.sdomain.boundaries(i);
        lnkBc.updateBorder(lnkDm.bcs.db(i));
      end
    end

    






    

    


  end



  methods (Access = private)

    function [cells, direction] = findBorderCell(obj,boundary)
      data  = obj.borderData(boundary);
      cells = data.cells;
      direction = data.direction;
    end


    % Functions related to update the physics.
    function updateCell(obj,ncell)
      %UPDATECELL update the model.

      % Definition of some pointers 
      lnk = obj.simData.solver(obj.nsolver);
      face = lnk.faces;
      elem = lnk.elements;
      mesh = lnk.mesh;
      dofm = lnk.dofm;

      % Update the information.
      % Find the cells around.
      [neighcells,~] = obj.sdomain.findNeighborhod(ncell.new.ijk);

      % -------------------------------------------------------------------
      % Update Geometric information.
      faceNew = obj.data2face(ncell,neighcells);

      fNewFaces = neighcells==0;
      nodes2Faces = faceNew.nod2Fac(fNewFaces,:)';
      face.nodes2Faces = [face.nodes2Faces; nodes2Faces(:)];
      face.mapN2F = [face.mapN2F; faceNew.mapN2F'];
      face.faces2Elements = [face.faces2Elements; faceNew.faces2Elements];
      face.mapF2E = [face.mapF2E; faceNew.mapF2E];
      face.faceNeighbors = [face.faceNeighbors; faceNew.facNeig(fNewFaces,:)];
      face.faceCentroid = [face.faceCentroid; ncell.surf.centroid(fNewFaces,:)];
      face.faceNormal = [face.faceNormal; ncell.surf.normal(fNewFaces,:).*ncell.surf.area(fNewFaces)'];
      face.nFaces = face.nFaces+ncell.new.nFaces;

      facesNot2Add = 6 - ncell.new.nFaces;
      existFaces = faceNew.faces2Elements(~fNewFaces,1);
      for i=1:facesNot2Add
        loc=face.faceNeighbors(existFaces(i),:)==0;
        face.faceNeighbors(existFaces(i),loc)=ncell.volu.gtCell+1;
      end

      elem.nCellsByType(2) = elem.nCellsByType(2) + 1;
      mesh.cellCentroid = [mesh.cellCentroid; ncell.volu.centroid];
      mesh.cellVolume = [mesh.cellVolume; ncell.volu.volume];

      nnodes = (ncell.point.node>ncell.point.gtNode);
      mesh.nNodes = mesh.nNodes+ncell.new.nNodes;
      mesh.nCells = mesh.nCells+1;
      mesh.coordinates = [mesh.coordinates;ncell.point.coord(nnodes',:)];
      mesh.cells = [mesh.cells;ncell.point.node];
      mesh.cellTag = [mesh.cellTag;ncell.new.cellTag];
      mesh.cellNumVerts = [mesh.cellNumVerts;8];
      mesh.cellVTKType = [mesh.cellVTKType;12];

      % -------------------------------------------------------------------
      % Update transmissibility information.
      trans = lnk.computeTransCell(ncell);
      lnk.trans = [lnk.trans; trans(fNewFaces)];
      lnk.isIntFaces = [lnk.isIntFaces; false(ncell.new.nFaces,1)];
      lnk.isIntFaces(ncell.surf.face(~fNewFaces)) = true;

      % -------------------------------------------------------------------
      % Update DOFmanager.
      dofm.numEntsField =dofm.numEntsField+1;
      dofm.totDoF = dofm.totDoF+1;
      dofm.numEntsSubdomain = dofm.numEntsSubdomain+1;
      dofm.fields.entCount = dofm.fields.entCount+1;
      dofm.fields.subdomainEnts{1}=[dofm.fields.subdomainEnts{1};ncell.volu.gtCell+1];
      dofm.fields.isEntActive=[dofm.fields.isEntActive;true];
      dofm.addCell(1,ncell.volu.gtCell+1);

      % Update the IJK  grid information.
      obj.sdomain.updateIJK(ncell.new.ijk,ncell.volu.gtCell+1, ...
        ncell.point.gtNode+ncell.new.nNodes,ncell.surf.gtFace+ncell.new.nFaces);
    end



    % % Functions related to update the physics.
    function data = data2face(obj,ncell,neighcells)
      shortcut = obj.simData.solver(obj.nsolver);
      % Update the information to Face.

      % bot = [1 2 3 4];
      % top = [5 6 7 8];
      % est = [2 6 7 3];
      % wst = [1 5 8 4];
      % sth = [1 2 6 5];
      % nth = [4 3 7 8];

      bot = [1 4 3 2];
      top = [6 7 8 5];
      est = [2 3 7 6];
      wst = [5 8 4 1];
      sth = [6 5 1 2];
      nth = [3 4 8 7];
      ref_pos = [wst; est; sth; nth; bot; top];

      data.nod2Fac(1:6,1:4)=0;
      data.facNeig(1:6,1:2)=ncell.volu.gtCell+1;
      for i=1:6
        pos = ref_pos(i,:);
        [~,nodeC]=min(ncell.point.node(pos));
        nodeL = mod(nodeC - 2, 4) + 1;
        nodeR = mod(nodeC, 4) + 1;
        nodeA = mod(nodeC + 1, 4) + 1;
        if (ncell.point.node(pos(nodeL))<ncell.point.node(pos(nodeR)))
          data.nod2Fac(i,:)=[
            ncell.point.node(pos(nodeC)) ncell.point.node(pos(nodeL)) ...
            ncell.point.node(pos(nodeA)) ncell.point.node(pos(nodeR))];
          data.facNeig(i,1)=neighcells(i);   % <---(i,1) or (i,2) to check
        else
          data.nod2Fac(i,:)=[
            ncell.point.node(pos(nodeC)) ncell.point.node(pos(nodeR)) ...
            ncell.point.node(pos(nodeA)) ncell.point.node(pos(nodeL))];
          data.facNeig(i,2)=neighcells(i);   % <---(i,1) or (i,2) to check
        end
      end
      data.mapN2F(1:ncell.new.nFaces)=max(shortcut.faces.mapN2F)+4*(1:ncell.new.nFaces);
      data.faces2Elements=[ncell.surf.face' (1:6)'];
      data.mapF2E=max(shortcut.faces.mapF2E)+6;
    end


    function data = borderData(obj,boundary)
      % Function to identify the grid information of a boundary
      lnkborder = obj.simData.bcs.db(boundary);
      lnkfaces = obj.simData.grid(obj.nsolver).faces;
      if strcmp(lnkborder.cond,'SurfBC')
        % assuming that the data is in the border!
        % TODO - still need to test if the neighbors is zeros,
        data.faces = lnkborder.data.entities;
        data.cells = sum(lnkfaces.faceNeighbors(data.faces,:),2)';
        loc = ismember(lnkfaces.faces2Elements(:,1),data.faces);
        data.direction = lnkfaces.faces2Elements(loc,2);        
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