classdef testDoFManager < matlab.unittest.TestCase
  
  properties
    model
    mesh
  end

  methods(TestClassSetup)
    % Shared setup for the entire test class
    function setup(testCase)
      testCase.model = ModelType(["Poromechanics_FEM","SinglePhaseFlow_FVTPFA"]);
      testCase.mesh = Mesh();
      testCase.mesh.importMesh('domain.msh');
    end
  end
   
  methods(Test)
    % Test methods
    function test1(testCase)
      dof = DoFManager(testCase.mesh,testCase.model); 
      verifyEqual(testCase,dof.getNumDoF('Poromechanics'),10527);
      verifyEqual(testCase,dof.getNumDoF('SinglePhaseFlow'),2800);
      verifyFalse(testCase,dof.isField('VariablySaturatedFlow'));
      verifyEqual(testCase,dof.getActiveSubdomain('Poromechanics'),1)
      verifyEqual(testCase,dof.getActiveSubdomain(["Poromechanics","SinglePhaseFlow"]),1)
      verifyEqual(testCase,dof.getDoF('Poromechanics',30),[88;89;90])
    end

    function test2(testCase)
      dof = DoFManager(testCase.mesh,testCase.model,'dof.dat');
      verifyEqual(testCase,dof.getNumDoF('Poromechanics'),10527);
      verifyEqual(testCase,dof.getNumDoF('SinglePhaseFlow'),2000);
      verifyFalse(testCase,dof.isField('VariablySaturatedFlow'));
      verifyEqual(testCase,dof.getActiveSubdomain('Poromechanics'),[1;2])
      verifyEqual(testCase,dof.getActiveSubdomain(["Poromechanics","SinglePhaseFlow"]),2)
      verifyEqual(testCase,dof.getDoF('Poromechanics',30),[88;89;90])
      verifyEqual(testCase,dof.getDoF('SinglePhaseFlow',1000),11127)
      verifyError(testCase,@() dof.getDoF('SinglePhaseFlow', 1),'dofError:inactiveEntity')
    end
  end
  
end