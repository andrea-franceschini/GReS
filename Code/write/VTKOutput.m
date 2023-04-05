classdef VTKOutput < handle
  % VTKOUTPUT General VTK output class, can be subclassed for specialized types

  properties (SetAccess = public, GetAccess = public)

    % None

  end

  properties (Access = private)

    mesh;
    hasFaults = false;
    surfaceList;
    isFolderReady = false;
    folderName = 'vtkOutput';
    cellFileName = 'domain.vtu';
    surfaceFileName = 'fracture.vtu';
    timestep = 0;
    glo2loc;
    surfaceNumNodes;
    surfaceNumElems;
    surfaceNumVerts;
    surfaceTag;
    surfaceVTKType;
    surfaceCoord;
    surfaceElems;
    is2DMeshReady = false;
    NotANumber = 1.e100;
    nDataMax = 100;
    IDdata = 0;
    data;

  end

  methods (Access = public)

    function obj = VTKOutput(mesh, varargin)
      obj.mesh = mesh;
      if (nargin > 1)
        obj.folderName = varargin(1);
        if (nargin > 2)
          error('Too many inputs.');
        end
      end
      obj.data = repmat(struct('time', 0, 'vtm', []), obj.nDataMax, 1);
    end

    function delete(obj)
      mesh = [];
      surfaceList = [];
      glo2loc = [];
      surfaceCoord = [];
      surfaceElems = [];
      data = [];
    end

    function finalize(obj)
      obj.writePVDFile();
      obj.delete();
    end

    function writeVTKFile(obj, time, pointData3D, cellData3D, pointData2D, cellData2D)
      vtmFileName = writeVTUFile(obj, time, pointData3D, cellData3D, pointData2D, cellData2D);
      obj.IDdata = obj.IDdata + 1;
      if (obj.IDdata > obj.nDataMax)
        obj.nDataMax = 2*obj.nDataMax;
        obj.data(obj.nDataMax) = obj.data(1);
      end
      obj.data(obj.IDdata).time = time;
      obj.data(obj.IDdata).vtm = vtmFileName;
    end

    function setVTKFolder(obj, folderName)
      obj.folderName = folderName;
    end

    function setSurfaces(obj, surfaceList)
      obj.surfaceList = surfaceList;
      obj.hasFaults = ~isempty(surfaceList);
    end

  end

  methods (Access = private)

    function createVTKFolder(obj)
      if (isfolder(obj.folderName))
        rmdir(obj.folderName, 's');
      end

      status = mkdir(obj.folderName);
      if (status ~= 1)
        error('Unable to create folder for VTK output.');
      end
      obj.isFolderReady = true;
    end

    function writePVDFile(obj)
      docNode = com.mathworks.xml.XMLUtils.createDocument('VTKFile');
      toc = docNode.getDocumentElement;
      toc.setAttribute('type', 'Collection');
      toc.setAttribute('version', '1.0');
      blocks = docNode.createElement('Collection');

      for i = 1 : obj.IDdata
        block = docNode.createElement('DataSet');
        block.setAttribute('timestep', sprintf('%e', obj.data(i).time));
        block.setAttribute('file', obj.data(i).vtm);
        blocks.appendChild(block);
      end

      toc.appendChild(blocks);

      fileName = sprintf('%s.pvd', obj.folderName);
      xmlwrite(fileName, docNode);
    end

    function fileName = writeVTMFile(obj)
      outName = setOutputName(obj);

      docNode = com.mathworks.xml.XMLUtils.createDocument('VTKFile');
      toc = docNode.getDocumentElement;
      toc.setAttribute('type', 'vtkMultiBlockDataSet');
      toc.setAttribute('version', '1.0');
      blocks = docNode.createElement('vtkMultiBlockDataSet');

      block1 = docNode.createElement('Block');
      block1.setAttribute('name', 'CellElementRegion');
      dataset = docNode.createElement('DataSet');
      dataset.setAttribute('name', 'Domain');
      dataset.setAttribute('file', sprintf('%s/%s', outName, obj.cellFileName));
      block1.appendChild(dataset);
      blocks.appendChild(block1);

      if (obj.hasFaults)
        block2 = docNode.createElement('Block');
        block2.setAttribute('name', 'SurfaceElementRegion');
        dataset = docNode.createElement('DataSet');
        dataset.setAttribute('name', 'Fracture');
        dataset.setAttribute('file', sprintf('%s/%s', outName, obj.surfaceFileName));
        block2.appendChild(dataset);
        blocks.appendChild(block2);
      end

      toc.appendChild(blocks);

      fileName = sprintf('%s/%s.vtm', obj.folderName, outName);
      xmlwrite(fileName, docNode);
    end

    function create2DMesh(obj)
      nSurf = obj.mesh.nSurfaces;
      list = zeros(nSurf,1);
      for i = 1 : length(obj.surfaceList)
        list = list + obj.mesh.findSurfacesOfRegion(obj.surfaceList{i});
      end
      list = list == 1;
      srf = obj.mesh.surfaces(list,:);
      obj.surfaceNumVerts = obj.mesh.surfaceNumVerts(list);
      obj.surfaceTag = obj.mesh.surfaceTag(list);
      obj.surfaceVTKType = obj.mesh.surfaceVTKType(list);
      obj.glo2loc = zeros(obj.mesh.nNodes,1);
      for i = 1 : size(srf,1)
        for j = 1 : obj.surfaceNumVerts(i)
          if (obj.glo2loc(srf(i,j)) == 0)
            obj.glo2loc(srf(i,j)) = 1;
          end
        end
      end
      nLoc = sum(obj.glo2loc);
      obj.glo2loc(obj.glo2loc==1) = 1 : nLoc;
      obj.surfaceCoord = obj.mesh.coordinates(obj.glo2loc>0,:);
      obj.surfaceElems = srf;
      for i = 1 : size(srf,1)
        obj.surfaceElems(i,1:obj.surfaceNumVerts(i)) = obj.glo2loc(obj.surfaceElems(i,1:obj.surfaceNumVerts(i)));
      end
      obj.surfaceNumNodes = nLoc;
      obj.surfaceNumElems = size(srf,1);
      obj.is2DMeshReady = true;
    end

    function vtmFileName = writeVTUFile(obj, time, pointData3D, cellData3D, pointData2D, cellData2D)
      if (~obj.isFolderReady)
        createVTKFolder(obj);
      end

      outName = setOutputName(obj);

      if (obj.hasFaults)
        if (~obj.is2DMeshReady)
          create2DMesh(obj);
        end
      end

      vtmFileName = writeVTMFile(obj);
      status = mkdir(sprintf('%s/%s', obj.folderName, outName));
      if (status ~= 1)
        error('Unable to create folder for VTK output.');
      end

      % Add fake fields to all the meshes (2D and 3D) to be able to always visualize all of them
      if (isfield(pointData3D, 'name'))
        pointData3DNames = extractfield(pointData3D, 'name');
      else
        pointData3DNames = [];
      end
      if (isfield(cellData3D, 'name'))
        cellData3DNames = extractfield(cellData3D, 'name');
      else
        cellData3DNames = [];
      end
      if (isfield(pointData2D, 'name'))
        pointData2DNames = extractfield(pointData2D, 'name');
      else
        pointData2DNames = [];
      end
      if (isfield(cellData2D, 'name'))
        cellData2DNames = extractfield(cellData2D, 'name');
      else
        cellData2DNames = [];
      end
      pointDataNames = unique([pointData3DNames, pointData2DNames]);
      cellDataNames = unique([cellData3DNames, cellData2DNames]);

      for i = 1 : length(pointDataNames)
        if (~ismember(pointDataNames{i}, pointData3DNames))
          for j = 1 : length(pointData2D)
            if (strcmp(pointData2D(j).name, pointDataNames{i}))
              ndims = size(pointData2D(j).data,2);
            end
          end
          nfields = length(pointData3D);
          pointData3D = setfield(pointData3D, {nfields+1, 1}, 'name', pointDataNames{i});
          pointData3D = setfield(pointData3D, {nfields+1, 1}, 'data', ones(obj.mesh.nNodes,ndims)*obj.NotANumber);
        end
        if (~ismember(pointDataNames{i}, pointData2DNames))
          for j = 1 : length(pointData3D)
            if (strcmp(pointData3D(j).name, pointDataNames{i}))
              ndims = size(pointData3D(j).data,2);
            end
          end
          nfields = length(pointData2D);
          pointData2D = setfield(pointData2D, {nfields+1, 1}, 'name', pointDataNames{i});
          %pointData2D = setfield(pointData2D, {nfields+1, 1}, 'data', ones(obj.surfaceNumNodes,1)*obj.NotANumber);
          pointData2D = setfield(pointData2D, {nfields+1, 1}, 'data', ones(obj.mesh.nNodes,ndims)*obj.NotANumber);
        end
      end

      for i = 1 : length(cellDataNames)
        if (~ismember(cellDataNames{i}, cellData3DNames))
          for j = 1 : length(cellData2D)
            if (strcmp(cellData2D(j).name, cellDataNames{i}))
              ndims = size(cellData2D(j).data,2);
            end
          end
          nfields = length(cellData3D);
          cellData3D = setfield(cellData3D, {nfields+1, 1}, 'name', cellDataNames{i});
          cellData3D = setfield(cellData3D, {nfields+1, 1}, 'data', ones(obj.mesh.nCells,ndims)*obj.NotANumber);
        end
        if (~ismember(cellDataNames{i}, cellData2DNames))
          for j = 1 : length(cellData3D)
            if (strcmp(cellData3D(j).name, cellDataNames{i}))
              ndims = size(cellData3D(j).data,2);
            end
          end
          nfields = length(cellData2D);
          cellData2D = setfield(cellData2D, {nfields+1, 1}, 'name', cellDataNames{i});
          cellData2D = setfield(cellData2D, {nfields+1, 1}, 'data', ones(obj.surfaceNumElems,ndims)*obj.NotANumber);
        end
      end

      fname = sprintf('%s/%s/%s', obj.folderName, outName, obj.cellFileName);
      mxVTKWriter(fname, time, obj.mesh.coordinates, obj.mesh.cells, obj.mesh.cellVTKType, ...
                  obj.mesh.cellNumVerts, pointData3D, cellData3D);

      if (obj.hasFaults)

        for i = 1 : length(pointData2D)
          pointData2D(i).data = pointData2D(i).data(obj.glo2loc>0,:);
        end

        fname = sprintf('%s/%s/%s', obj.folderName, outName, obj.surfaceFileName);
        mxVTKWriter(fname, time, obj.surfaceCoord, obj.surfaceElems, obj.surfaceVTKType, ...
                    obj.surfaceNumVerts, pointData2D, cellData2D);
      end

      updateCounter(obj);
    end

    function outName = setOutputName(obj)
      outName = sprintf('output_%5.5i', obj.timestep);
    end

    function updateCounter(obj)
      obj.timestep = obj.timestep + 1;
    end

  end

end
