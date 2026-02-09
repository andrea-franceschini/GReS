classdef preconditioner < handle
   properties (Access = private)

      % Nonsymmetry tolerance
      nsyTol

      % Flag for debug
      DEBUGflag = false

      % General solver params
      generalsolver

      % Discretizer
      domain

      % Number of domains/interfaces
      nDom
      nInt

      % Flag to treat multiple domains as multiple domains
      multidomFlag = false

      % Flag to know if the problem has multiphysics
      multiPhysFlag = false

      % Max Threads
      maxThreads

      % Preconditioner Type
      PrecType

      % Preconditioner
      Prec = []

   end


   properties (GetAccess = public,SetAccess = private)

      % Physics
      phys

      % Params struct
      params

      % Test Space
      TV0

      % Preconditioner application
      Apply_L = []
      Apply_R = []

   end

   methods (Access = public)

      % Function to compute the preconditioner
      Compute(obj,A)

   end


   methods (Static,Access = public)

      % Constructor of the preconditioner object, specifies if it cannot be used as not supported
      function [obj, useChronos] = create(debugflag,nsyTol,generalsolver,usrInput)

         % Initialize an empty class
         obj = preconditioner.empty;
         useChronos = false;

         domainin = generalsolver.domains;

         % Check if the problem comes from multiphysics
         multiPhysFlag = false;
         if(domainin(1).dofm.getNumberOfVariables() > 1)
            multiPhysFlag = true;
            return
         end

         nDom = generalsolver.nDom;
         nInt = generalsolver.nInterf;

         % Check the number of interfaces and domains
         if nInt ~= 0
            interfacein = generalsolver.interfaces;
         end

         % Select the physics
         physname = domainin(1).dofm.getVariableNames();
         if ~multiPhysFlag
            % Supported Single Physics
            if(contains(physname,'pressure'))
               phys = 0;
            elseif(physname == 'displacements')
               phys = 1;
               % Check if there is contact 
               if any(cellfun(@(o) isa(o,'SolidMechanicsContact'),interfacein))
                  phys = 1.1;
               end
            else
               if debugflag
                  warning('Non supported Physics for preconditioner, falling back to matlab solver');
               end
               return
            end
         else
            % Supported MultiPhysics
            if(contains(physname,['pressure' 'displacements']))
               if domainin.dofm.getVariableNames(1) == 'pressure'
                  phys = 0;
               else
                  phys = 1;
               end
            else
               if debugflag
                  warning('Non supported Physics for linsolver, falling back to matlab solver');
               end
               return
            end
         end

         % Now the preconditioner can actually be built, the checks have been passed
         obj = preconditioner(debugflag,nsyTol,generalsolver,domainin,multiPhysFlag,phys,usrInput);
         useChronos = true;
      end

   end
   
   methods (Access = private)

      % Constructor Function
      function obj = preconditioner(debugflag,nsyTol,generalsolver,domainin,multiPhysFlag,phys,usrInput)

         % Use the debugflag set into the linearsolver
         obj.DEBUGflag = debugflag;
         obj.multiPhysFlag = multiPhysFlag;
         obj.generalsolver = generalsolver;
         obj.domain = generalsolver.domains;
         obj.nDom = generalsolver.nDom;
         obj.nInt = generalsolver.nInterf;
         obj.phys = phys;
         obj.nsyTol = nsyTol;

         % Get default values
         if obj.phys == 0
           chronos_xml_default = fullfile(gres_root,'Code','linsolver','XML_setup','chronos_xml_setup_CFD.xml');
         else
           chronos_xml_default = fullfile(gres_root,'Code','linsolver','XML_setup','chronos_xml_setup.xml');
         end

         % Read Defaults
         data = readstruct(chronos_xml_default,AttributeSuffix="");

         % Get the preconditioner type
         obj.PrecType = lower(data.preconditioner);

         % Get the different parameters according to the prectype
         switch obj.PrecType
            case 'amg'
               obj.params.amg      = data.amg;
               obj.params.smoother = data.smoother;
               obj.params.prolong  = data.prolong;
               obj.params.coarsen  = data.coarsen;
               obj.params.tspace   = data.tspace;
               obj.params.filter   = data.filter;
               obj.params.minIter  = 30;

            case 'fsai'
               obj.params.smoother = data.smoother;
               obj.params.minIter  = 300;
         end

         % Set maximum number of threads to use if the system provides less
         obj.maxThreads = maxNumCompThreads;
         obj.params.smoother.nthread = min(obj.params.smoother.nthread,obj.maxThreads);
         obj.params.prolong.np = min(obj.params.prolong.np,obj.maxThreads);
         obj.params.filter.np = min(obj.params.filter.np,obj.maxThreads);

         % Get user prescribed values
         obj.params = obj.getUserInput(obj.params,usrInput);

      end

      % Function to get the user input parameters for the preconditioner
      [params] = getUserInput(obj,params,xml);

      % Function to compute the preconditioner for the single block (single physics)
      computeSinglePhPrec(obj,A);

      % Function to compute the Reverse Agumented preconditioner for the lagrange multiplier case (single physics multi domain)
      computeRACP(obj,A)

      % Function to compute the MCP preconditioner for the multiphysics case
      computeMCP(obj,A)

      % Function to treat the Dirichlet boundary conditions
      treatDirBC(obj,A)
   end

end


