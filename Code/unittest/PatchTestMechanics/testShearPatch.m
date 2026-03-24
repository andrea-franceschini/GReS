classdef testShearPatch < matlab.unittest.TestCase

  properties
    pathToFile
  end

  methods(TestClassSetup)
    % Shared setup for the entire test class
    function setupOnce(testCase)
      testCase.pathToFile = fullfile(gres_root,...
        "/Code/unittest/PatchTestMechanics/shearPatch.xml");
    end
  end

  methods(TestMethodSetup)
    % Setup for each test
  end

  methods(Test)
    % Test methods

    function testPatchXML(testCase)
      gresLog().setVerbosity(-2);
      input = readInput(testCase.pathToFile);
      simparams = SimulationParameters(input.SimulationParameters);
      mesh = structuredMesh(1,1,1,[0 1],[0 1],[0 1]);
      elems = Elements(mesh,2);
      faces = Faces(mesh);
      grid = struct('topology',mesh,'cells',elems,'faces',faces);
      mat = Materials(input.Materials);
      bc = Boundaries(grid,input.BoundaryConditions);
      printUtils = OutState('printTimes',1,'outputFile',"test");
      domain = Discretizer('Boundaries',bc,...
                           'Materials',mat,...
                           'Grid',grid);
      domain.addPhysicsSolvers(input.Solver);
      solver = NonLinearImplicit('simulationparameters',simparams,'domain',domain,'output',printUtils);
      solver.simulationLoop();
      gresLog().setVerbosity(-1);
      verifyEqual(testCase,domain.state.data.stress(:,5),1.5*ones(8,1),"AbsTol",1e-9)

      % validate the vtk output
      s = readstruct("test/output_00001/Domain_1.vtu","FileType","xml");
      v = s.UnstructuredGrid.Piece.PointData.DataArray.Text;
      v = str2num(v);
      vv = zeros(8,3);
      vv(5:end,1) = 1.0;
      verifyEqual(testCase,v,vv,"AbsTol",1e-9)

    end

    function testPatchKV(testCase)
      gresLog().setVerbosity(-2);
      simparams = SimulationParameters('Start',0.0,...
        'End',3e0,...
        'DtInit',1.0,...
        'DtMax',1.0,...
        'DtMin',1.0,...
        'RelativeTolerance',1.e-8,...
        'AbsoluteTolerance',1.e-9);

      mesh = structuredMesh(1,1,1,[0 1],[0 1],[0 1]);
      elems = Elements(mesh,2);
      faces = Faces(mesh);
      grid = struct('topology',mesh,'cells',elems,'faces',faces);
      mat = Materials();
      mat.addSolid('name',"solid",'cellTags',1);
      mat.addConstitutiveLaw("solid","Elastic",'youngModulus',1e0,'poissonRatio',0.0);

      bc = Boundaries(grid);
      bc.addBC('name',"fix_bottom",...
        'type',"Dirichlet",...
        'field',"surface",...
        'variable',"displacements",...
        'entityListType',"tags",...
        'entityList',1,...
        'components',[1,2,3]);

      bc.addBCEvent('fix_bottom','time',0,'value',0.0);

      bc.addBC('name','top_fix',...
        'type',"Dirichlet",...
        'field',"surface",...
        'variable',"displacements",...
        'entityListType',"tags",...
        'entityList',2,...
        'components',3);

      bc.addBCEvent('top_fix','time',0,'value',0.0);

      bc.addBC('name','top_disp',...
        'type',"Dirichlet",...
        'field',"surface",...
        'variable',"displacements",...
        'entityListType',"tags",...
        'entityList',2,...
        'components',1);

      bc.addBCEvent('top_disp','time',0,'value',0.0);
      bc.addBCEvent('top_disp','time',5.0,'value',5.0);


      domain = Discretizer('Boundaries',bc,...
        'Materials',mat,...
        'Grid',grid);
      domain.addPhysicsSolver('Poromechanics');
      solver = NonLinearImplicit('simulationparameters',simparams,'domain',domain);
      solver.simulationLoop();
      gresLog().setVerbosity(-1);
      verifyEqual(testCase,domain.state.data.stress(:,5),1.5*ones(8,1),"AbsTol",1e-9)


    end
  end

end