classdef testGrowning < matlab.unittest.TestCase

  properties
    model
    mesh
    elems
    faces
    material
    dofmanager
    linSyst
    boundary
    growning
    state
    % solver
    % simState
  end
  
  methods (TestClassSetup)
    % Shared setup for the entire test class
    function setup(testCase)
      testCase.model = ModelType("SinglePhaseFlow_FVTPFA");
      simParam = SimulationParameters('Inputs/simParam.dat',testCase.model);
      GaussPts = Gauss(12,2,3);

      % Create the Mesh related object
      testCase.mesh = Mesh();
      testCase.mesh.importGMSHmesh('Inputs/mesh.msh');      
      testCase.elems = Elements(testCase.mesh,GaussPts);
      testCase.faces = Faces(testCase.model,testCase.mesh);

      % Create an object of the Materials class
      testCase.material = Materials(testCase.model,'Inputs/materialsList.dat');

      % Degree of freedom manager
      testCase.dofmanager = DoFManager(testCase.mesh,testCase.model);

      % Create Discretizer
      grid = struct('topology',testCase.mesh,'cells',testCase.elems, ...
        'faces',testCase.faces);
      testCase.linSyst = Discretizer(testCase.model,simParam, ...
        testCase.dofmanager,grid,testCase.material,GaussPts);

      % Creating Initial and Boundaries conditions.
      testCase.state = testCase.linSyst.setState();
      testCase.state.pressure(:) = 1.e5;
      fileName = ["Inputs/BC_BoundA.dat","Inputs/BC_BoundB.dat"];
      testCase.boundary = Boundaries(fileName,testCase.model,grid);

      % % Create and set the print utility
      % printUtils = OutState(obj.model,obj.mesh, ...
      %   'Inputs/outTime.dat','folderName','Outputs');
      % 
      % obj.solver = FCSolver(obj.model,simParam, ...
      %   obj.dofmanager,grid,obj.material, ...
      %   obj.boundary,printUtils,obj.state,obj.linSyst, ...
      %   'GaussPts',GaussPts,'SaveRelError',true);
      % obj.simState = obj.solver.NonLinearLoop();
      % printUtils.finalize()

      % Creating the object to grow the mesh
      testCase.growning = GrowningDomain(testCase.linSyst,testCase.boundary);
    end
  end

  methods (TestMethodSetup)
    % Setup for each test
  end

  methods (Test)
    % Test methods

    function GrownMesh(testCase)
      statek=testCase.state;
      stateTmp=testCase.state;
      [statek,stateTmp] = testCase.growning.addCells(testCase.linSyst,1, ...
        [2 4 6 8],6,statek,stateTmp);

      verifyEqual(testCase,length(statek.pressure),12);
      verifyEqual(testCase,length(stateTmp.pressure),12);
    end
  end

end