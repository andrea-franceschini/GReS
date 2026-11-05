classdef MortarSide < uint8
  % Enum class to access master/slave quantities
  
  enumeration
    master (2)
    slave  (1)
  end
  
  methods

    function side = getSide(obj)
      side = uint8(obj);
    end

  end
end

