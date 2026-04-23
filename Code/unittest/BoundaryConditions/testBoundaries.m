classdef testBoundaries < matlab.unittest.TestCase

  properties
    bc
    pathToFile
  end

  methods(TestClassSetup)
    % Shared setup for the entire test class
    function setupOnce(testCase)
      testCase.pathToFile = fullfile(gres_root,...
        "/Code/unittest/BoundaryConditions/boundaryConditions.xml");
    end
  end

  methods(TestMethodSetup)
    % Setup for each test
  end

  methods(Test)
    % Test methods

    function readBCFile(testCase)
      % test bc added with  xml file
      grid = Grid();
      grid.importMesh("mesh.msh");
      testCase.bc = Boundaries(grid,testCase.pathToFile);

      validate(testCase)



    end

    function readBC(testCase)
      % test bc added with key-value input
      grid = Grid();
      grid.importMesh("mesh.msh");
   
      testCase.bc = Boundaries(grid);

      testCase.bc.addBC("name","bc1",...
        "type","Dirichlet",...
        "field","surface",...
        "variable","pressure",...
        "entityListType","bcList",...
        "entityList",[1,13,168]);
      testCase.bc.addBCEvent("bc1",'time',"2*t",'value',"10-z");


      testCase.bc.addBC('name',"bc2",...
        'type',"Neumann",...
        'field',"surface",...
        'entityListType',"tags",...
        'entityList',2,...
        'variable',"displacements",...
        'components',[1 3])
      testCase.bc.addBCEvent("bc2",'time',1,'value',-10);
      testCase.bc.addBCEvent("bc2",'time',3,'value',-5);
      testCase.bc.addBCEvent("bc2",'time',6,'value',"z");


      testCase.bc.addBC("name","bc3",...
        "type","Dirichlet",...
        "field","node",...
        "components",1,...
        "variable","displacements",...
        "entityListType","bcListFile",...
        "entityList",'entityList');
      testCase.bc.addBCEvent("bc3",'time',1,'value',5);


      testCase.bc.addBC("name","bc4",...
        "type","Dirichlet",...
        "field","surface",...
        "components","x",...
        "variable","displacements",...
        "entityListType","tags",...
        "entityList",4);
      testCase.bc.addBCEvent("bc4",'time',0,'value','bcvals_list.dat');


      testCase.bc.addBC("name","bc5",...
        "type","Dirichlet",...
        "field","surface",...
        "components",["x","y","z"],...
        "variable","displacements",...
        "entityListType","tags",...
        "entityList",1);
      testCase.bc.addBCEvent("bc5",'time',0,'value',0);
      testCase.bc.addBCEvent("bc5",'time',5,'value',4);
      testCase.bc.addBCEvent("bc5",'time',3,'value',2);


      testCase.bc.addBC("name","bc6",...
        "type","Source",...
        "field","cell",...
        "components",["x","y","z"],...
        "variable","displacements",...
        "entityListType","box",...
        "entityList",[0.2 0.8 0.2 0.8 3.0 6.0]);
      testCase.bc.addBCEvent("bc6",'time',0,'value',0);
      testCase.bc.addBCEvent("bc6",'time',1.0,'value',1.0);

      validate(testCase)

    end

  end


  methods

    function validate(testCase)

      bcs = testCase.bc;

      tol = 1e-6;

      bcs.initialize("bc1",'node')
      bcs.initialize("bc2",'node')
      bcs.initialize("bc3",'node')
      bcs.initialize("bc4",'node')
      bcs.initialize("bc5",'node')
      bcs.initialize("bc6",'node')

      e1 = bcs.getTargetEntities("bc1");
      e2 = bcs.getTargetEntities("bc2");
      e3 = bcs.getTargetEntities("bc3");
      e4 = bcs.getTargetEntities("bc4");
      e5 = bcs.getDofs("bc5");
      es6 = bcs.getSourceEntities("bc6");
      e6 = bcs.getTargetEntities("bc6");

      v1 = bcs.getVals("bc1",3);
      v21 = bcs.getVals("bc2",0.2);
      v22 = bcs.getVals("bc2",2);
      v23 = bcs.getVals("bc2",10);
      v3 = bcs.getVals("bc3",1.1);
      v4 = bcs.getVals("bc4",2);
      v5 = bcs.getVals("bc5",0.0);
      v6 = bcs.getVals("bc6",0.5);

      % verify bcs
      % bc1
      verifyEqual(testCase,e1([1; end]),[1;170],"AbsTol",tol)
      verifyEqual(testCase,v1([1; 5; end-1]),[60.0;0.0;34.5],"AbsTol",tol)
      % bc2
      verifyEqual(testCase,all(e2(1:9)==e2(10:18)),true,"AbsTol",tol)
      verifyEqual(testCase,v21([1; 5; end]),[-0.625; -1.25; -2.5],"AbsTol",tol)
      verifyEqual(testCase,mean(v22),-0.833333333333,"AbsTol",tol)
      verifyEqual(testCase,all(v23(1:9)==v23(10:18)),true,"AbsTol",tol)
      % bc3
      %verifyEqual(testCase,all(e3==load('entityList')),true,"AbsTol",tol)
      verifyEqual(testCase,all(v3==5),true,"AbsTol",tol)
      % bc5
      verifyEqual(testCase,e5([1;10;19]),[1;2;3],"AbsTol",tol)
      verifyEqual(testCase,max(v5),0,"AbsTol",tol)
      % bc6
      verifyEqual(testCase,mean(e6),1.016666666666667e+02,"AbsTol",tol)
      verifyEqual(testCase,mean(es6),  39.500000000000000,"AbsTol",tol)
      verifyEqual(testCase,mean(v6),0.023809523809524,"AbsTol",tol)
    end
  end

end
