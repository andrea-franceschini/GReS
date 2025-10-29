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

    function addCell(obj,refCell,dir,state1,state2)
      %ADDCELL Add a cell in the grid and update the mesh

      lnk = obj.simData.solver(obj.nsolver);

      % Find the position to where to grow the mesh.
      [n_cells,n_ijk] = obj.sdomain.findNeighborhodOfCell(refCell);
      if obj.canGrow(dir,n_cells)
        % Find the information about the cell to be created.
        [cellsNeigh,~] = obj.sdomain.findNeighborhod(n_ijk(dir,:));
        cell = HexaCell(lnk.mesh,lnk.faces,cellsNeigh,dir,0.5,lnk.mesh.cellTag(refCell));

        % Update the mesh
        obj.updateCell(cell,n_ijk(dir,:));

        % Update the reference of the boundary.
        obj.sdomain.updateBorder(lnk,cell,obj.sdomain.grid(refCell,:),n_ijk(dir,:));

        % Update the states.
        obj.updateState(cell,state1,state2);

        obj.updateBoundary();
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
    function updateCell(obj,ncell,ijk_newCell)
      %UPDATECELL update the model.

      % Definition of some pointers 
      lnk = obj.simData.solver(obj.nsolver);
      face = lnk.faces;
      % elem = lnk.elements;
      mesh = lnk.mesh;
      dofm = lnk.dofm;

      % Update the information.

      % -------------------------------------------------------------------
      % Update Geometric information.
      face.nodes2Faces = [face.nodes2Faces; ncell.newNodes2Faces];
      face.mapN2F = [face.mapN2F; face.mapN2F(end)+ncell.mapN2F'];
      face.faces2Elements = [face.faces2Elements; ncell.faces2Elements];
      face.mapF2E = [face.mapF2E; face.mapF2E(end)+ncell.mapF2E];
      face.faceCentroid = [face.faceCentroid; ncell.faceCentroid(ncell.facesNew,:)];
      face.faceNormal = [face.faceNormal;
        ncell.faceNormal(ncell.facesNew,:).*ncell.faceArea(ncell.facesNew)'];
      face.nFaces = face.nFaces+ncell.newFaces;
      face.faceNeighbors = ncell.faceNeighbors(face);

      mesh.cellCentroid = [mesh.cellCentroid; ncell.centroid];
      mesh.cellVolume = [mesh.cellVolume; ncell.volume];

      mesh.nNodes = mesh.nNodes+ncell.newNodes;
      mesh.nCells = mesh.nCells+1;
      mesh.coordinates = [mesh.coordinates;ncell.coordinates];
      mesh.cells = [mesh.cells;ncell.nodes];
      mesh.cellTag = [mesh.cellTag;ncell.cellTag];
      mesh.cellNumVerts = [mesh.cellNumVerts;ncell.cellNumVerts];
      mesh.cellVTKType = [mesh.cellVTKType;ncell.cellVTKType];

      % -------------------------------------------------------------------
      % Update transmissibility information.
      % trans = lnk.computeTransCell(ncell);
      trans = ncell.computeTrans(lnk.material.getMaterial(ncell.tag).PorousRock.getPermVector());      
      lnk.trans = [lnk.trans; trans(ncell.facesNew)];
      lnk.isIntFaces = [lnk.isIntFaces; false(ncell.newFaces,1)];
      lnk.isIntFaces(ncell.faces(~ncell.facesNew)) = true;

      % -------------------------------------------------------------------
      % Update DOFmanager.
      dofm.numEntsField =dofm.numEntsField+1;
      dofm.totDoF = dofm.totDoF+1;
      dofm.numEntsSubdomain = dofm.numEntsSubdomain+1;
      dofm.fields.entCount = dofm.fields.entCount+1;
      dofm.fields.subdomainEnts{1}=[dofm.fields.subdomainEnts{1};ncell.lastCell+1];
      dofm.fields.isEntActive=[dofm.fields.isEntActive;true];
      dofm.addCell(1,ncell.lastCell+1);

      % Update the IJK  grid information.
      obj.sdomain.updateIJK(ijk_newCell,ncell.lastCell+1, ...
        ncell.lastNode+ncell.newNodes,ncell.lastFace+ncell.newFaces);
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
        if size(data,1) == cell.lastFace
          loc=cell.newNodesMirrorNodes;        
        end
        if size(data,1) == cell.lastCell
          loc = cell.cellNeigh(cell.faceGrow);
        end
        state1.data=setfield(state1.data,labels{label},[data; data(loc,:)]);
        state2.data=setfield(state2.data,labels{label},[data; data(loc,:)]);
      end
    end

  end
end