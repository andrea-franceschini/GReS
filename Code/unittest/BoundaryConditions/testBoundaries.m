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
      mesh = Mesh();
      mesh.importMesh("mesh.msh");
      elems = Elements(mesh,2);
      faces = Faces(mesh);
      grid = struct("topology",mesh,"cells",elems,"faces",faces);
      testCase.bc = Boundaries(grid,testCase.pathToFile);

      validate(testCase)



    end

    function readBC(testCase)
      % test bc added with key-value input
      mesh = Mesh();
      mesh.importMesh("mesh.msh");
      elems = Elements(mesh,2);
      faces = Faces(mesh);
      grid = struct("topology",mesh,"cells",elems,"faces",faces);

      testCase.bc = Boundaries(grid);

      testCase.bc.addBC("name","kbc1",...
        "type","Dirichlet",...
        "field","surface",...
        "variable","pressure",...
        "entityListType","bcList",...
        "entityList",[1,13,168]);
      testCase.bc.addBCEvent("kbc1",'time',"2*t",'value',"10-z");


      testCase.bc.addBC("name","kbc2",...
        "type","Dirichlet",...
        "field","node",...
        "components",1,...
        "variable","displacements",...
        "entityListType","surfacetags",...
        "entityList",4);
      testCase.bc.addBCEvent("kbc2",'time',0,'value',"bcvals_list.dat");

      testCase.bc.addBC('name',"kbc3",...
        'type',"Neumann",...
        'field',"surface",...
        'entityListType',"surfaceTags",...
        'entityList',2,...
        'variable',"displacements",...
        'components',[1 3])
      testCase.bc.addBCEvent("kbc3",'time',1,'value',-10);
      testCase.bc.addBCEvent("kbc3",'time',3,'value',-5);
      testCase.bc.addBCEvent("kbc3",'time',6,'value',"z");

      testCase.bc.addBC("name","kbc4",...
        "type","Dirichlet",...
        "field","surface",...
        "components",["x","y","z"],...
        "variable","displacements",...
        "entityListType","surfaceTags",...
        "entityList",1);
      testCase.bc.addBCEvent("kbc4",'time',0,'value',0);
      testCase.bc.addBCEvent("kbc4",'time',5,'value',4);
      testCase.bc.addBCEvent("kbc4",'time',3,'value',4);

      testCase.bc.addBC("name","kbc4",...
        "type","Dirichlet",...
        "field","surface",...
        "components",["x","y","z"],...
        "variable","displacements",...
        "entityListType","surfaceTags",...
        "entityList",1);
      testCase.bc.addBCEvent("kbc4",'time',0,'value',0);
      testCase.bc.addBCEvent("kbc4",'time',5,'value',4);
      testCase.bc.addBCEvent("kbc4",'time',3,'value',4);

      validate(testCase)

    end

    function validate(testCase)

      tol = 1e-9;

      testCase.bc.computeTargetEntities("bc1",'node')
      testCase.bc.computeTargetEntities("bc2",'node')
      testCase.bc.computeTargetEntities("bc3",'node')
      testCase.bc.computeTargetEntities("bc4",'node')
      testCase.bc.computeTargetEntities("bc5",'node')

      e1 = testCase.bc.getTargetEntities("bc1");
      e2 = testCase.bc.getTargetEntities("bc2");
      e3 = testCase.bc.getTargetEntities("bc3");
      e4 = testCase.bc.getTargetEntities("bc4");
      e5 = testCase.bc.getTargetEntities("bc5");
      e6 = testCase.bc.getTargetEntities("bc6");

      v1 = testCase.bc.getVals("bc1",3);
      v21 = testCase.bc.getVals("bc2",0.2);
      v22 = testCase.bc.getVals("bc2",2);
      v23 = testCase.bc.getVals("bc2",10);
      v3 = testCase.bc.getVals("bc3",1.1);
      v4 = testCase.bc.getVals("bc4",2);
      v5 = testCase.bc.getVals("bc5",0.0);
      v6 = testCase.bc.getVals("bc6",0.5);

      verifyEqual(testCase,v1([1 5 end-1]),[60;0.0;34.5],"AbsTol",tol)
      verifyEqual(testCase,v21([1 5 end]),[-0.625 -1.25 -2.5],"AbsTol",tol)
      verifyEqual(testCase,v21([1 5 end]),[-0.625 -1.25 -2.5],"AbsTol",tol)
      verifyEqual(testCase,mean(v23),10.0,"AbsTol",1e-9)
      verifyEqual(testCase,length(v3),126,"AbsTol",1e-9)
      verifyEqual(testCase,mean(v3),5,"AbsTol",1e-9)
      verifyEqual(testCase,length(v4),126,"AbsTol",1e-9)
      verifyEqual(testCase,max(v5),0,"AbsTol",1e-9)
      verifyEqual(testCase,length(v5),12,"AbsTol",1e-9)

      verifyEqual(testCase,e1,[1;13;168],"AbsTol",1e-9)
      verifyEqual(testCase,e2,repmat((165:168)',2,1),"AbsTol",1e-9)
      verifyEqual(testCase,length(e3),126,"AbsTol",1e-9)
      verifyEqual(testCase,e5,repmat((1:4)',3,1),"AbsTol",1e-9)

      
      % verifyEqual(testCase,v1,[60;34.50;0.0],"AbsTol",1e-8)
      % verifyEqual(testCase,mean(v21),-10.0,"AbsTol",1e-9)
      % verifyEqual(testCase,mean(v22),-7.5,"AbsTol",1e-9)
      % verifyEqual(testCase,mean(v23),10.0,"AbsTol",1e-9)
      % verifyEqual(testCase,length(v3),126,"AbsTol",1e-9)
      % verifyEqual(testCase,mean(v3),5,"AbsTol",1e-9)
      % verifyEqual(testCase,length(v4),126,"AbsTol",1e-9)
      % verifyEqual(testCase,max(v5),0,"AbsTol",1e-9)
      % verifyEqual(testCase,length(v5),12,"AbsTol",1e-9)
      % 
      % verifyEqual(testCase,e1,[1;13;168],"AbsTol",1e-9)
      % verifyEqual(testCase,e2,repmat((165:168)',2,1),"AbsTol",1e-9)
      % verifyEqual(testCase,length(e3),126,"AbsTol",1e-9)
      % verifyEqual(testCase,e5,repmat((1:4)',3,1),"AbsTol",1e-9)
    end

  end

end
