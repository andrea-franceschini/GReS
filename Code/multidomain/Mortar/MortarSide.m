classdef MortarSide
  %MORTARSIDE Summary of this class goes here
  %   Detailed explanation goes here
  
  enumeration
    master (1)
    slave (2)
  end
  
  methods

    % additional logic for retrieving quantities from master and slave
    % mortar side

    function side = getSide(side)
      side = double(side);
    end

  end
end

