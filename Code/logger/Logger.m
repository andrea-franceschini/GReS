
classdef Logger < handle
  % Error logger in GReS
  
  properties (Access = private)
    verbosity = 1
  end

  methods 
    function obj = Logger()
    end

    function setVerbosity(obj,lev)
      obj.verbosity = lev;
    end

    function lev = getVerbosity(obj)
      lev = obj.verbosity;
    end

    function log(obj,varargin)
      if isnumeric(varargin{1})
        loglevel = varargin{1};
        a = 2;
      else
        loglevel = -1;
        a = 1;
      end
      if obj.verbosity >= loglevel
        fprintf(varargin{a:end});
        %fprintf("\n");
      end
    end

    function error(obj,varargin)
      if isnumeric(varargin{1})
        loglevel = varargin{1};
        a = 2;
      else
        loglevel = 0;
        a = 1;
      end
      if obj.verbosity > loglevel
        error(varargin{a:end});
      end
    end

    function warning(obj,varargin)
      if isnumeric(varargin{1})
        loglevel = varargin{1};
        a = 2;
      else
        loglevel = 0;
        a = 1;
      end
      if obj.verbosity > loglevel
        warning(varargin{a:end});
      end
    end

  end

  methods (Static)
    
    function welcomeMsg()
      disp(' ')
      disp('===============================================')
      disp('              Welcome to GReS!                 ')
      disp('===============================================')
      disp(' ')
      disp('To get started:')
      disp(' ')
      disp('Run  >> compileAll')
      disp('     to compile all existing MEX files in the project.')
      disp('     You only need to run it once!')
      disp(' ')
      disp('Check out the tutorials folder')
      disp('     for hands-on guides and usage examples.')
      disp('     <a href="matlab:open(strcat(gres_root,''/Tutorial/quickStart.mlx''))">Start the tutorial</a>')
      disp(' ')
      disp('Explore the Tests repository')
      disp('     to check available simulations in GReS.')
      disp(' ')
      disp(' ')
    end

  end

end

