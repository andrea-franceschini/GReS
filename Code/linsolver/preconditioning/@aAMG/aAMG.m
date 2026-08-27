classdef aAMG < preconditioner
%   aAMG adaptive Algebraic MultiGrid preconditioner.
%
%   This class implements an adaptive Algebraic MultiGrid (aAMG)
%   preconditioner derived from the abstract 'preconditioner' base class.
%   It is responsible for:
%     - Reading default AMG-related parameters from XML setup files based
%       on the physical problem (pressure, velocity, displacements).
%     - Allowing user overrides of those parameters via solver XML input.
%     - Constructing and storing the internal preconditioner operator(s),
%       along with any symmetry information and threading limits.
%     - Exposing a Compute method to build the preconditioner for a given
%       matrix and Apply_L / Apply_R handles to apply the preconditioner.
%
%   The constructor initializes defaults, enforces thread limits based on
%   system capabilities, and merges user-specified parameters. The class
%   keeps implementation details private and presents a simple interface
%   for preconditioner creation and application.

   properties (GetAccess = public,SetAccess = private)
       % Preconditioner
      Prec = []
      
      % Symmetry of the matrix on which the preconditioner has been
      % computed
      PrecSym = true

      % Physics
      phys
      
      % Params struct
      params
      
      % Max Threads
      maxThreads

      % Generalsolver
      generalsolver

   end

   methods (Access = public)

      % Function to compute the preconditioner
      Compute(obj,A,sym,varargin)

      % Getter for the function handle to apply the left preconditioner
      function x = ApplyLeft(obj,b,varargin)
         if nargin < 3
            error('Not enough arguments for AMG apply Left');
         end
         A = varargin{1};
         
         x = AMG_Vcycle(obj.Prec,A,b);
      end

      % Getter for the function handle to apply the right preconditioner
      function x = ApplyRight(obj,b,varargin)
         if nargin < 2
            error('Not enough arguments for AMG apply Right');
         end

         x = b;
      end

      % Constructor Function
      function obj = aAMG(debugflag,generalsolver,physname)

         % Call the constructor of the abstract class
         obj = obj@preconditioner();

         % Set the debugflag
         obj.DEBUGflag = debugflag;

         % Get General Solver
         obj.generalsolver = generalsolver;

         % Get default values
         if(contains(physname,'pressure') || contains(physname,'u'))
            obj.phys = 0;
            chronos_xml_default = fullfile(gres_root,'Code','linsolver','XML_setup','chronos_xml_setup_CFD.xml');
         elseif(contains(physname,'displacements')) 
            obj.phys = 1;
            chronos_xml_default = fullfile(gres_root,'Code','linsolver','XML_setup','chronos_xml_setup.xml');
         else
            disp(physname);
            error('Non supported Physics for preconditioner');
         end

         % Read Defaults
         data = readInput(chronos_xml_default);

         % Set maximum number of threads to use if the system provides less
         obj.maxThreads = maxNumCompThreads;

         % Get the different parameters 
         obj.params.amg      = data.amg;
         obj.params.smoother = data.smoother;
         obj.params.prolong  = data.prolong;
         obj.params.coarsen  = data.coarsen;
         obj.params.tspace   = data.tspace;
         obj.params.filter   = data.filter;
         obj.params.prolong.np = min(obj.params.prolong.np,obj.maxThreads);
         obj.params.filter.np = min(obj.params.filter.np,obj.maxThreads);

         % Get user prescribed values
         obj.params.smoother.nthread = min(obj.params.smoother.nthread,obj.maxThreads);
         obj.params = obj.getUserInput(obj.params,generalsolver.simparams.linSolverParams);
      end
   end

   methods (Access = private)

      % Function to get the user input parameters for the preconditioner
      [params] = getUserInput(obj,params,xml);
   end
end


