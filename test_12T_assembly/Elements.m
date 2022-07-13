classdef Elements < handle
  % ELEMENT General element class

  properties (Access = private)
   % Element type
     elemType = [];
   % Element data
     data = [];
  end
   
  methods (Access = public)
      function obj = Elements(nNode,nodeCoords,nElem,elemTopol,surfTopol)
        obj.setElementData(nNode,nodeCoords,nElem,elemTopol,surfTopol)
      end
  
  
      function el = getElement(obj,elemType)
        if elemType == 10
            el = Tetrahedron(obj.data);
        elseif elemType == 12
            el = Hexahedron(obj.data);
        else 
           error('Element not available');
        end
     end
  end
  
  methods (Access = private)
      function obj = setElementData(obj,nNode,nodeCoords,nElem,elemTopol,surfTopol)
          obj.data = {nNode nodeCoords nElem elemTopol surfTopol};
      end
  end
  
end

