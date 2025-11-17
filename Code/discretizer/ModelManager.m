classdef ModelManager < handle
  % General model class
  % This class stores all the information related to the problem at hand
  % and store a map to existing domain solvers and 



  properties (GetAccess=private, SetAccess=private)
    domains        % db for physics solver in the model
    interfaces     % db for interface solvers in the model
    simparams      % simulation parameters
    t = 0          % simulation time
    tStep = 0      % simulation time step
  end


  methods (Access = public)
    function obj = ModelManager(varargin)

      obj.domains = containers.Map('KeyType','double','ValueType','any');
      obj.interfaces = containers.Map('KeyType','double','ValueType','any');

      if ~isempty(varargin)
        obj.simparams = varargin{1};
      end

    end


    function runProblem(obj,solutionScheme)

      % this method implement the simulation time loop
      % handle to a function implementing a specific type of solution
      % strategy for each time step

      while obj.t < obj.domainManager.simparams.tMax

        % Update the simulation time and time step ID
        obj.tStep = obj.tStep + 1;
        % new time update to fit the outTime list
        [obj.t, dt] = obj.updateTime(flConv, delta_t);

        isConverged = solutionScheme(obj.t,dt);

        delta_t = manageNextTimeStep(obj,dt,flConv);
        
      end
      %
    end

    % function printSolution()
    %   % replaces finalize state
    % end

    function createModel(obj,fileName)
      % utility to create an entire GReS model from a single xml file

      str = readstruct(fileName,AttributeSuffix="");

      obj.simparams = SimulationParameters(fileName);

      for i = 1:numel(str.Domain)
        createDomain(obj,str.Domain(i));
      end

      if isfield(str,"Interface")
        for iI = 1:numel(str.Interface)
          addInterfaceSolver(obj,str.Interface(iD));
        end
      else
        if numel(str.Domain) > 1
          warning("Running a problem with multiple domain " + ...
            "without any connecting interfaces");
        end
      end
    end

  
    function createDomain(obj,varargin)
      n = numel(keys(obj.domains));
      obj.domains(n+1) = Domain(varargin{:});
    end

    function addInterfaceSolver(obj,varargin)

      % TO DO 
    end

    function out = getDomain(obj,id)
      out = obj.domains(id);
    end

    function out = getInterface(obj,id)
      out = obj.interfaces(id);
    end


  end

  methods (Static)

    function solver = createSolverFromFile(fileName)

      inputStruct = readstruct(fileName,AttributeSuffix="");

      meshFile = getXMLData(inputStruct,[],"Geometry");
      mesh = Mesh();
      mesh.importMesh(meshFile);

      ng = getXMLData(inputStruct,1,"GaussPoints");
      
      elements = Elements(mesh,ng);
      faces = Faces(mesh)

      mat = Materials()



      
    end

  end
end
