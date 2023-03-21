classdef Elements < handle
  % ELEMENT General element class
  
  properties (Access = public)
    cellCentroid
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
    end
    
    % Calling the specific element class based on elemType
    function el = getElement(obj,elemType)
       if elemType == 10
          el = Tetrahedron(obj.mesh);
       elseif elemType == 12
          el = Hexahedron(obj.mesh,obj.GaussPts);
       else 
          error('Element not available');
       end
    end

  end
  
  methods (Access = private)
    % Assigns element data (properties of the 'Mesh' object)
    function setElementData(obj,nIn,mesh,data)
      obj.mesh = mesh;
      if nIn > 1 
        obj.GaussPts = data{1};
      end
      computeCellCentroid(obj);
    end
    
    function computeCellCentroid(obj)
      obj.cellCentroid = sparse(repelem(1:obj.mesh.nCells,obj.mesh.cellNumVerts), ...
          nonzeros((obj.mesh.cells)'),repelem((obj.mesh.cellNumVerts).^(-1),obj.mesh.cellNumVerts),obj.mesh.nCells,obj.mesh.nNodes) ...
          * obj.mesh.coordinates;
    end
  end
  
end