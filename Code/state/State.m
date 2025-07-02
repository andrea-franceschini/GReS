classdef State < handle & matlab.mixin.Copyable
  % Handle for state troughout simulation

  properties
    t = 0
    data
  end

  methods

  end

  methods (Access = protected)
    function cp = copyElement(obj)
      cp = State();
      cp.t = obj.t;
      cp.data = obj.data;
    end
  end
end

