classdef State < handle & matlab.mixin.Copyable
  % Handle for state troughout simulation

  properties 
    t = 0
    data = struct('curr',[],'init',[],'prev',[])
  end

  methods

  end

  methods (Access = protected)
    function cp = copyElement(obj,src)
      switch src
      cp = State();
      cp.t = obj.t;
      cp.data = obj.data;
      end

      function state = getState(obj,src)

        state = obj.data.(src);
        
      end


  end
end

