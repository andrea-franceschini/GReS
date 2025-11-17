classdef Domain < handle

  % General model class
  % This class stores all the information related to the problem at hand
  % and store a map to the existing domainSolver and interfaceSolver



  properties (GetAccess=private, SetAccess=private)
    physicsSolvers        % db for physics solver in the model
  end

  % domain specific properties
  properties (GetAccess=public, SetAccess=private)
    dofm
    bcs                 % boundary conditions
    outstate            % printing utilities
    materials           % materials
    % grid struct
    grid = struct('topology',[],'faces',[],'cells',[]) 
    elements            % db for finite elements
    faces               % db for finite volumes
    variables           % list of fields with the same order of the global system
  end

  properties (GetAccess=private, SetAccess=public)
    state               % class holding any state variable in the domain
    stateOld            % class holding any state variable in the domain at the last converged time step
  end


  methods (Access = public)
    function obj = Domain(varargin)

      obj.setDomain(varargin{:})

      % the database for physicsSolvers in the domain
      obj.physicsSolvers = containers.Map('KeyType','char','ValueType','any');

    end


    function setDomain(obj,varargin)

      % create a domain from a key-value set of parameters or an xml file

      if isscalar(varargin)
        % input is a struct coming from an xml file
        obj.createDomainFromFile(varargin{1});
      else 
        % input is a list of key-value pair
        nIn = numel(varargin);
        assert(nIn == 8,"Error in setDomain(): " + ...
          "Incorrect key-value inputs. Keys must be:" + ...
          "'grid', 'materials', 'boundaryConditions', 'output'");
        fixedInput = struct(varargin{1:6});

        % read fixed input parameters
        obj.grid = fixedInput.grid;
        obj.materials = fixedInput.materials;
        obj.bcs = fixedInput.boundaryConditions;
        obj.outstate = fixedInput.output;

        % create the DoFManager
        obj.dofm = DoFManagerNew(obj.grid);

        % create the State object
        obj.state = State();

      end

    end

    function addPhysicsSolver(obj,solverName,varargin)

      % add a physic solver to an existing domain

      % ways to manually add a physics solver 
      % 1) a scalar struct coming from an xml file
      % 2) a key-value list of input, with additional solver specific input at the end

      if isscalar(varargin)
        % input is a xml file
        solverInput = varargin{1};
      else 
        % input is a list of key-value pair of solver inputs
        solverInput = struct(varargin{:});
      end

      % create solver
      pSolv = feval(solverName,obj,solverInput);

      % add solver to the database
      n = numel(keys(obj.physicsSolvers));
      obj.physicsSolvers(n+1) = pSolv;

    end

    function addInterfaceSolver(obj,varargin)

      % TO DO
    end

    function out = getPhysicsSolver(obj,id)
      out = obj.domainSolvers(id);
    end

    function out = getInterfaceSolver(obj,id)
      out = obj.interfaceSolvers(id);
    end

    function out = getState(obj)
      out = obj.state;
    end

    function out = getOldState(obj)
      out = obj.stateOld;
    end


    function createDomainFromFile(obj,inputStruct)

      meshFile = getXMLData(inputStruct,[],"Geometry");
      mesh = Mesh();
      mesh.importMesh(meshFile);

      ng = getXMLData(inputStruct,1,"GaussPoints");

      obj.grid.topology = mesh;
      obj.grid.cells = Elements(mesh,ng);
      obj.grid.faces = Faces(mesh);

      obj.materials = Materials(inputStruct);
      obj.bcs = Boundaries(inputStruct,obj.grid);
      obj.outstate = OutState(mesh,inputStruct);

      % add physics solvers
      solvers = inputStruct.Solver;
      solvNames = fieldnames(solvers);

      % create the DoFManager
      obj.dofm = DoFManagerNew(obj.grid);

      % create the State object
      obj.state = State();

      for s = 1:numel(solvNames)
        % define solver (the xml field must match the name of the solver
        % class)
        obj.addPhysicsSolver(solvNames{s},solvers.(solvNames{s}));
      end

    end







  end

  methods (Static)

  end
end
