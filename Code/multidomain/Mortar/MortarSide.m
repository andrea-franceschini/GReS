classdef MortarSide
  % Enum class to access master/slave quantities
  
  enumeration
    master 
    slave
  end
  
  methods

    function side = getSide(obj)
      switch obj
        case MortarSide.master
          side = 1;
        case MortarSide.slave 
          side = 2;
      end
    end

  end
end

