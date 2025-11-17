classdef problemSolver < handle
  % General model class
  % This class stores all the information related to the problem at hand
  % and store a map to the existing domainSolver and interfaceSolver



  properties (GetAccess=private, SetAccess=private)
    physicsSolvers        % db for physics solver in the model
    interfaceSolvers      % db for interface solvers in the model
    simparams             % simulation parameters
    t = 0                 % simulation time
    tStep = 0             % simulation time step
  end


  methods (Access = public)
    function obj = problemSolver(varargin)

      obj.physicsSolvers = containers.Map('KeyType','double','ValueType','any');
      obj.interfaceSolvers = containers.Map('KeyType','double','ValueType','any');

      if ~isempty(varargin)
        obj.simparams = varargin{1};
      end

    end


    function runProblem(obj,solutionScheme)

      % this method implement the simulation time loop
      % handle to a function implementing a specific type of solution
      % strategy for each time step

      while obj.t < obj.domain.simparams.tMax

        % Update the simulation time and time step ID
        obj.tStep = obj.tStep + 1;
        % new time update to fit the outTime list
        [obj.t, dt] = obj.updateTime(flConv, delta_t);

        isConverged = solutionScheme(obj.t,dt);

        delta_t = manageNextTimeStep(obj,dt,flConv);
        
      end
      %
    end

    function printSolution()
      % replaces finalize state
    end

    function createModel(obj,fileName)
      % utility to create an entire GReS model from an xml file
      str = readstruct(fileName,AttributeSuffix="");

      obj.simparams = SimulationParameters(fileName);

      for iD = 1:numel(str.Domain)
        addPhysicsSolver(obj,fileName);
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

    function addPhysicsSolver(obj,varargin)

      % ways to manually add a physics solver 
      % 1) a scalar struct
      % 2) a key-value list of input, with additional solver specific input at the end

      if isscalar(varargin)
        % input is a xml file
        fileName = varargin{1};
        pSolv = obj.createSolverFromFile(fileName);
      else 
        % input is a list of key-value pair
        nIn = numel(varargin);
        assert(nIn >= 8,"Error in addPhysicsSolver(): " + ...
          "Not enough input parameters")
        fixedInput = struct(varargin{1:8});

        % read fixed input parameters
        solverName = fixedInput.solverName;
        grid = fixedInput.grid;
        material = fixedInput.materials;
        bcs = fixedInput.boundaryConditions;

        if nIn == 9
          % extra input is already a structure from a <Solver> xml field
          solverInput = varargin{9};
        else
          % extra input is a key-value pair (limited functionality)
          solverInput = struct(varargin{9:end});
        end
        pSolv = feval(solverName,obj,grid,material,bcs,solverInput);
      end

      % add solver to the database
      n = numel(keys(obj.physicsSolvers));
      obj.physicsSolvers(n+1) = pSolv;

    end

    function addInterfaceSolver(obj,varargin)

      % TO DO 
    end

    function out = getPhysicsSolver(obj,id)
      out = obj.physicsSolvers(id);
    end

    function out = getInterfaceSolver(obj,id)
      out = obj.interfaceSolvers(id);
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
