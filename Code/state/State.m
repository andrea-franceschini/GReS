classdef State < handle
  % State class 
  % Store solution fields (either primary or auxiliary) throughout the
  % simulation
  % It refers to a unique simulation time
  % Stores vartiables in a structure array with fields:
  % init: initial variable (balanced)
  % old: variables from last converged step
  % curr: variables at current iteration


  properties (Access = private)

    data = struct('init',[],'old',[],'curr',[])

  end

  properties 
    t = 0
  end


  methods

    function state = get(obj,src,var)
      if nargin == 2
        state = obj.data.(src);
      elseif nargin == 3
        if ~isfield(obj.data.(src),var)
          error("Variable '%s' does not exist in the State object",var)
        end
        state = obj.data.(src).(var);
      end
    end


    function set(obj,src,val,var)
      if nargin < 4
        obj.data.(src) = val;
      elseif nargin == 4
        if ~isfield(obj.data.(src),var)
          error("Variable '%s' does not exist in the State object",var)
        end
        obj.data.(src).(var) = val;
      end
    end


    function valInterp = interpolate(obj,fac,var)

      if fac < 0 || fac > 1
        error("Interpolation factor must be a real number from 0 to 1.")
      end

      % interpolate field in current and old state

      if nargin == 2
        fldNames = reshape(string(fieldnames(obj.data.curr)),1,[]);
        for f = fldNames
          valInterp.(f) = obj.interpolate(fac,f);
        end
      elseif nargin == 3
        vCurr = obj.get('curr',var);
        vOld = obj.get('old',var);
        valInterp = vCurr*fac + vOld*(1-fac);

      end

    end
  end

end

