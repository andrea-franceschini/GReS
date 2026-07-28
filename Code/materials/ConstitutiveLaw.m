classdef ConstitutiveLaw < handle
  %CONSTITUTIVELAW Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    status = struct('curr',[],'conv',[])      % state variables
  end
  
  methods (Abstract)

    initializeStatus(obj)

    constitutiveUpdate(obj)

  end


  methods


    function commitStatus(obj)

      obj.status.conv = obj.status.curr;

    end

  end
end

