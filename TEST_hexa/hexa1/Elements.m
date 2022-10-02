classdef Elements < handle
  % ELEMENT General element class

  properties (Access = private)
    % Element type
    elemType = [];
    % Element data
    data = [];
  end
   
  methods (Access = public)
    % Class constructor method   
    function obj = Elements(varargin)
      % Calling the function to set element data    
      obj.setElementData(varargin)
    end
    
    % Calling the specific element class based on elemType
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
     % Assigns element data (properties of the 'Mesh' object)
     function obj = setElementData(obj,varargin)
        obj.data = varargin{1};
     end
  end
  
end

