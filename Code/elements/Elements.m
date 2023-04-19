classdef Elements < handle
  % ELEMENT General element class
  
  properties (Access = public)
    cellCentroid
    vol
    tetra
    hexa
  end

  properties (Access = private)
    % Element type
%     elemType;
    % Mesh object
    mesh
    GaussPts = []
  end
   
  methods (Access = public)
    % Class constructor method   
    function obj = Elements(mesh,varargin)
      % Calling the function to set element data
      data = varargin;
      nIn = nargin;
      obj.setElementData(nIn,mesh,data);
      obj.computeCellProperties;
    end
    
%     function v = getVolume(obj,el)
%       v = obj.vol(el);
%     end
    
%     % Calling the specific element class based on elemType
%     function elemObj = getElement(obj,el)
%       if obj.mesh.cellVTKType(el) == 10
%         elemObj = obj.tetra;
%       elseif obj.mesh.cellVTKType(el) == 12
%         elemObj = obj.hexa;
%       end
%     end
    
%     function elemObj = getHexaObj(obj)
%       elemObj = obj.hexa;
%     end
    
%     function elemObj = getTetraObj(obj)
%       elemObj = obj.tetra;
%     end

  end
  
  methods (Access = private)
    % Assigns element data (properties of the 'Mesh' object)
    function setElementData(obj,nIn,mesh,data)
      obj.mesh = mesh;
      if nIn > 1 
        obj.GaussPts = data{1};
      end
      %
      if any(obj.mesh.cellVTKType == 10)
        obj.tetra = Tetrahedron(obj.mesh);
      end
      if any(obj.mesh.cellVTKType == 12)
        obj.hexa = Hexahedron(obj.mesh,obj.GaussPts);
      end
      obj.vol = zeros(obj.mesh.nCells,1);
      obj.cellCentroid = zeros(obj.mesh.nCells,3);
    end
    
    function computeCellProperties(obj)
      idCell = 1:obj.mesh.nCells;
      idTetra = idCell(obj.mesh.cellVTKType == 10);
      idHexa = idCell(obj.mesh.cellVTKType == 12);
      if ~isempty(idTetra)
        obj.vol(idTetra) = obj.tetra.findVolume(idTetra);
        obj.cellCentroid(idTetra,:) = computeCentroidGeneral(obj,idTetra);
      end
      if ~isempty(idHexa)
         [obj.vol(idHexa),obj.cellCentroid(idHexa,:)]= obj.hexa.findVolumeAndCentroid(idHexa);
%          obj.vol(idHexa) = tmpVol;
%          obj.cellCentroid(idHexa,:) = tmpCentroid;
      end
    end
    
    function tetraCtr = computeCentroidGeneral(obj,idTetra)
      tetraCtr = sparse(repelem(1:length(idTetra),obj.mesh.cellNumVerts(idTetra)), ...
          nonzeros((obj.mesh.cells(idTetra,:))'),repelem((obj.mesh.cellNumVerts(idTetra)).^(-1),obj.mesh.cellNumVerts(idTetra)),length(idTetra),obj.mesh.nNodes) ...
          * obj.mesh.coordinates;
%     end
    end
%     function computeTetraCentroid(obj)
%       obj.cellCentroid = sparse(repelem(1:obj.mesh.nCells,obj.mesh.cellNumVerts), ...
%           nonzeros((obj.mesh.cells)'),repelem((obj.mesh.cellNumVerts).^(-1),obj.mesh.cellNumVerts),obj.mesh.nCells,obj.mesh.nNodes) ...
%           * obj.mesh.coordinates;
%     end
  end
  
end