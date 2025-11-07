classdef Logger < handle
  % Error logger in GReS
  
  properties (Access = private)
    verbosity = 0
  end

  methods 
    function obj = Logger()
    end

    function setVerbosity(obj,lev)
      obj.verbosity = lev;
    end

    function log(obj,loglevel,varargin)
      if obj.verbosity > loglevel
        fprinft(varargin{:});
      end
    end

    function error(obj,loglevel,varargin)
      if obj.verbosity > loglevel
        error(varargin{:});
      end
    end

    function warning(obj,loglevel,varargin)
      if obj.verbosity > loglevel
        error(varargin{:});
      end
    end

  end

  methods (Static)
    
    function welcomeMsg(obj)
    end

  end
  
end

