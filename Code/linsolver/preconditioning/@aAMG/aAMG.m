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

      % AMG info
      AMG_info

      % Hierarchy construction timing records
      buildTimes = struct('smoother', {}, ...
                          'testspace', {}, ...
                          'coarsening', {}, ...
                          'compatibleRelaxation', {}, ...
                          'prolongation', {}, ...
                          'filterProlongation', {}, ...
                          'filterCoarseLevel', {}, ...
                          'lowRank', {}, ...
                          'total', {})

   end

   methods (Access = public)

      % Function to compute the preconditioner
      A = Compute(obj,A,sym,varargin)

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

      % Print info
      function AMGinfo(obj)
         print_AMG_info(obj.AMG_info)
         if ~isempty(obj.buildTimes)
            obj.printBuildTimes();
         end
      end

      % Print hierarchy construction timing statistics
      function printBuildTimes(obj)
         n = numel(obj.buildTimes);
         avgSmoo = mean([obj.buildTimes.smoother]);
         avgTspa = mean([obj.buildTimes.testspace]);
         avgCoar = mean([obj.buildTimes.coarsening]);
         avgCR   = mean([obj.buildTimes.compatibleRelaxation]);
         avgProl = mean([obj.buildTimes.prolongation]);
         avgFilP = mean([obj.buildTimes.filterProlongation]);
         avgFilC = mean([obj.buildTimes.filterCoarseLevel]);
         avgLowR = mean([obj.buildTimes.lowRank]);
         avgTot  = mean([obj.buildTimes.total]);

         fprintf('\n');
         fprintf('-------------------------------------------------------------------\n');
         fprintf('             AMG HIERARCHY CONSTRUCTION TIMING (s)\n');
         fprintf('-------------------------------------------------------------------\n');
         fprintf(' Number of recorded constructions: %d\n\n', n);
         fprintf(' %-30s %15s %12s\n', 'Phase', 'Avg Time (s)', '% Total');
         fprintf(' %-30s %15.4e %11.2f%%\n', 'Smoother setup', avgSmoo, (avgSmoo/avgTot)*100);
         fprintf(' %-30s %15.4e %11.2f%%\n', 'Test space setup', avgTspa, (avgTspa/avgTot)*100);
         fprintf(' %-30s %15.4e %11.2f%%\n', 'Coarsening', avgCoar, (avgCoar/avgTot)*100);
         if avgCR > 0
            fprintf(' %-30s %15.4e %11.2f%%\n', 'Compatible relaxation', avgCR, (avgCR/avgTot)*100);
         end
         fprintf(' %-30s %15.4e %11.2f%%\n', 'Prolongation', avgProl, (avgProl/avgTot)*100);
         fprintf(' %-30s %15.4e %11.2f%%\n', 'Prolongation filtering', avgFilP, (avgFilP/avgTot)*100);
         fprintf(' %-30s %15.4e %11.2f%%\n', 'Coarse level filtering', avgFilC, (avgFilC/avgTot)*100);
         if avgLowR > 0
            fprintf(' %-30s %15.4e %11.2f%%\n', 'Low-rank updates', avgLowR, (avgLowR/avgTot)*100);
         end
         fprintf('-------------------------------------------------------------------\n');
         fprintf(' %-30s %15.4e %12s\n', 'Total construction time', avgTot, '100.00%');
         fprintf('-------------------------------------------------------------------\n\n');
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

      % Record hierarchy construction times from global timing variables
      function recordBuildTimes(obj, totalTime)
         amg_times;
         tRec = struct('smoother', T_smoo, ...
                       'testspace', T_tspa, ...
                       'coarsening', T_coar, ...
                       'compatibleRelaxation', T_CR, ...
                       'prolongation', T_prol, ...
                       'filterProlongation', T_FilProl, ...
                       'filterCoarseLevel', T_FilCLev, ...
                       'lowRank', T_LowRank, ...
                       'total', totalTime);
         if isempty(obj.buildTimes)
            obj.buildTimes = tRec;
         else
            obj.buildTimes(end+1) = tRec;
         end
      end
   end
end


