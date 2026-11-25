classdef gresUITest < matlab.unittest.TestCase
%GRESUITEST Unit tests for gresUI auto-scanner, tickless section-diff export, depth colors, layout order, consensus scheme, multi-XML I/O, and execution.

  methods (Test)

    function testScanCoreSchemas(testCase)
      schema = scanGReSParameters();

      spNames = [schema.simParams.name];
      testCase.verifyFalse(any(spNames == "LinearSolver"));
      testCase.verifyEqual(sort(spNames([schema.simParams.required])), ...
        sort(["End", "DtInit", "DtMin", "DtMax"]));

      solvers = schema.physicsSolvers;
      sNames = [solvers.name];
      testCase.verifyFalse(any(ismember(["PhysicsSolver", "SinglePhaseFlow"], sNames)));

      mpEntry = solvers(sNames == "MultiphasePoromechanics");
      testCase.verifyTrue(mpEntry.isCoupled);
      msParam = mpEntry.params([mpEntry.params.name] == "mechanicsSolver");
      testCase.verifyEqual(msParam.defaultValue, "Poromechanics");
      testCase.verifyTrue(all(ismember(["Poromechanics", "PreStressMechanics"], msParam.choices)));

      topNames = [schema.linearSolver.topLevel.name];
      testCase.verifyTrue(all(ismember( ...
        ["useMatlab", "useIterative", "scalingFlag", "useEW", "tol", "eta0", "saveMatType", "saveMatSolves", "cprType"], ...
        topNames)));

      flowSecs = schema.linearSolver.chronos.Flow;
      testCase.verifyEqual([flowSecs.section], ["amg", "smoother", "prolong", "coarsen", "tspace"]);
      smootherNames = [flowSecs([flowSecs.section] == "smoother").params.name];
      testCase.verifyFalse(any(ismember(["nupre", "nupost", "nthread", "method"], smootherNames)));
    end

    function testScanExtendedSchemas(testCase)
      schema = scanGReSParameters();

      schNames = [schema.solutionSchemes.name];
      testCase.verifyTrue(all(ismember(["NonLinearImplicit", "FixedStressSplit", "EvolvingGrid"], schNames)));
      fss = schema.solutionSchemes(schNames == "FixedStressSplit");
      testCase.verifyTrue(all(ismember(["maxiterations", "reltolerance"], [fss.params.name])));

      gNames = [schema.geometry.name];
      testCase.verifyEqual(gNames, ["StructuredMesh", "BlockStructuredMesh", "meshFile"]);
      sm = schema.geometry(gNames == "StructuredMesh");
      testCase.verifyEqual([sm.params.name], ["NX", "NY", "NZ", "LX", "LY", "LZ"]);
      testCase.verifyTrue(all([sm.params.required]));

      outNames = [schema.outputs.params.name];
      testCase.verifyTrue(all(ismember(["outputFile", "matFileName", "printTimes", "vtkFormat"], outNames)));
      testCase.verifyTrue(any([schema.outputs.restartParams.name] == "filename"));

      testCase.verifyTrue(all(ismember(["specificWeight", "compressibility", "dynamicViscosity"], [schema.materials.fluidParams.name])));
      lawNames = [schema.materials.constitutiveLaws.name];
      testCase.verifyTrue(all(ismember(["BiElastic", "DruckerPrager", "Elastic", "HypoPlastic", "SedimentMaterial"], lawNames)));
      testCase.verifyFalse(any(ismember(["TransvElastic", "SoftSoilCreep", "CoulombFriction"], lawNames)));
      testCase.verifyTrue(all(ismember(["preset", "tabular", "burdine", "mualem"], [schema.materials.vanGenuchten.type])));

      bcReq = [schema.boundaryConditions.bcParams.name];
      testCase.verifyTrue(all(ismember(["name", "variable", "field", "entityList", "entityListType", "type"], bcReq)));
      entParam = schema.boundaryConditions.bcParams(bcReq == "entityList");
      testCase.verifyTrue(contains(entParam.helpText, "1 = Bottom"));
      testCase.verifyTrue(contains(entParam.helpText, "2 = Top"));
      testCase.verifyEqual([schema.boundaryConditions.eventParams.name], ["time", "value"]);

      iModels = schema.interfaces.models;
      testCase.verifyTrue(all(ismember(["MeshTying", "MeshTyingTPFA", "SolidMechanicsContact"], [iModels.name])));
      testCase.verifyTrue(all(ismember(["masterDomain", "slaveDomain", "masterSurface", "slaveSurface", "multiplierType", "variable"], [schema.interfaces.commonParams.name])));
      testCase.verifyTrue(all(ismember(["SegmentBasedQuadrature", "ElementBasedQuadrature", "RBFquadrature"], [schema.interfaces.quadratures.name])));

      smc = iModels([iModels.name] == "SolidMechanicsContact");
      testCase.verifyEqual([smc.coulombParams.name], ["cohesion", "frictionAngle"]);
      testCase.verifyTrue(all([smc.coulombParams.required]));
      testCase.verifyEqual([smc.activeSetParams.name], ["resetActiveSet", "forceStickBoundary"]);
      testCase.verifyTrue(all(ismember(["sliding", "normalGap", "maxStateChange"], [smc.activeSetTolParams.name])));
      testCase.verifyEqual(smc.requiredField, "displacements");

      tpfa = iModels([iModels.name] == "MeshTyingTPFA");
      testCase.verifyFalse(tpfa.supportsQuadrature);
      testCase.verifyFalse(tpfa.supportsMultiplierType);
      testCase.verifyEqual(tpfa.requiredField, "pressure");
    end

    function testNoTickBoxesAndSectionDefaultDiffExport(testCase)
      app = gresUI('Visible', false);
      cleanupApp = onCleanup(@() delete(app));

      % Verify no blank per-parameter enable checkboxes exist in the figure
      allBoxes = findobj(app.Figure, 'Type', 'uicheckbox');
      boxTexts = strtrim(string({allBoxes.Text}));
      testCase.verifyFalse(any(boxTexts == ""));

      % Untouched BiotFullyCoupled & LinearSolver should emit empty params/topLevel
      app.selectPhysicsSolver("BiotFullyCoupled", 1);
      cfg0 = app.collectConfiguration();
      testCase.verifyEmpty(fieldnames(cfg0.domains(1).solver.params));
      testCase.verifyEmpty(fieldnames(cfg0.linearSolver.topLevel));

      % Modifying 1 parameter in BiotFullyCoupled (BulkModulus -> 1D) exports all visible section params
      bulkField = getWidgetForParam(app, 'BulkModulus');
      bulkField.Value = '1D';
      cfg1 = app.collectConfiguration();
      testCase.verifyEqual(cfg1.domains(1).solver.params.BulkModulus, "1D");
      testCase.verifyEqual(cfg1.domains(1).solver.params.mechanicsSolver, "Poromechanics");
      testCase.verifyEqual(cfg1.domains(1).solver.params.flowSolver, "SinglePhaseFlowFVTPFA");

      % Modifying 1 parameter in LinearSolver (useEW -> true) exports all visible LinearSolver topLevel params
      ewBox = getWidgetForParam(app, 'useEW');
      ewBox.Value = true;
      ewBox.ValueChangedFcn([], []);
      testCase.verifyEqual(getRowVis(app, 'eta0'), "on");
      testCase.verifyEqual(getRowVis(app, 'tol'), "off");
      cfg2 = app.collectConfiguration();
      testCase.verifyEqual(cfg2.linearSolver.topLevel.useEW, 1);
      testCase.verifyTrue(isfield(cfg2.linearSolver.topLevel, 'eta0'));
      testCase.verifyTrue(isfield(cfg2.linearSolver.topLevel, 'saveMatType'));
      testCase.verifyFalse(isfield(cfg2.linearSolver.topLevel, 'tol'));

      % TwoPhaseFlow reveals cprType; SinglePhaseFlowFVTPFA hides it; essential hidden on dirichlet BC
      app.selectPhysicsSolver("TwoPhaseFlow", 1);
      testCase.verifyEqual(getRowVis(app, 'cprType'), "on");
      app.selectPhysicsSolver("SinglePhaseFlowFVTPFA", 1);
      testCase.verifyEqual(getRowVis(app, 'cprType'), "off");

      app.addBoundaryCondition(struct('name', "bcD", 'type', "dirichlet"), 1);
      testCase.verifyEqual(getRowVis(app, 'essential'), "off");
      app.addBoundaryCondition(struct('name', "bcC", 'type', "customRobin", 'essential', 1), 1);
      essLbls = findobj(app.Figure, 'Type', 'uilabel', 'Text', 'essential');
      essVis = arrayfun(@(l) string(l.Parent.Visible), essLbls);
      testCase.verifyTrue(any(essVis == "on"));
      testCase.verifyTrue(any(essVis == "off"));
    end

    function testDepthColorsDomainOrderAndBottomScheme(testCase)
      app = gresUI('Visible', false);
      cleanupApp = onCleanup(@() delete(app));

      imgs = findobj(app.Figure, 'Type', 'uiimage');
      testCase.verifyNotEmpty(imgs);
      testCase.verifyTrue(endsWith(string(imgs(1).ImageSource), "logo_GReS.png"));

      % Add solid (with PorousRock), 2 BCs (with BCevent), Domain #2, and Interface
      app.selectPhysicsSolver("Poromechanics", 1);
      app.addSolid(struct('name', "s1", 'cellTags', 1, 'PorousRock', struct('porosity', 0.3)), 1);
      app.addBoundaryCondition(struct('name', "bc1", 'BCevent', struct('time', 0, 'value', 0)), 1);
      d2 = app.addDomain();
      app.selectPhysicsSolver("Poromechanics", d2);
      app.addInterface("MeshTying", struct('masterDomain', 2, 'slaveDomain', 1));

      % Verify all Remove buttons use uniform crimson red styling
      rmNames = ["Remove Domain", "Remove Solid", "Remove BC", "Remove Event", "Remove Interface"];
      for k = 1:numel(rmNames)
        rBtns = findobj(app.Figure, 'Type', 'uibutton', 'Text', char(rmNames(k)));
        testCase.verifyNotEmpty(rBtns);
        for b = 1:numel(rBtns)
          testCase.verifyEqual(rBtns(b).BackgroundColor, [0.80, 0.18, 0.18], 'AbsTol', 1e-6);
          testCase.verifyEqual(rBtns(b).FontColor, [1, 1, 1], 'AbsTol', 1e-6);
        end
      end

      % Verify top-level section order places Problem Solver at Row 6 (underneath everything)
      simBtn = findSectionButton(app.Figure, "Simulation Parameters");
      domBtn = findSectionButton(app.Figure, "Domains");
      probBtn = findSectionButton(app.Figure, "Problem Solver (SolutionScheme)");
      testCase.verifyEqual(simBtn.Parent.Layout.Row, 1);
      testCase.verifyEqual(domBtn.Parent.Layout.Row, 2);
      testCase.verifyEqual(probBtn.Parent.Layout.Row, 6);

      % Verify depth colors: Depth 0 Blue -> Depth 1 Green -> Depth 2 Yellow -> Depth 3 Lavender -> Depth 4 Rose
      bluePanel = findobj(simBtn.Parent, 'Type', 'uipanel', '-depth', 1);
      testCase.verifyEqual(bluePanel.BackgroundColor, [0.91, 0.95, 1.00], 'AbsTol', 1e-6);

      dom1Btn = findSectionButton(app.Figure, "Domain #1");
      dom2Btn = findSectionButton(app.Figure, "Domain #2");
      testCase.verifyEqual(string(dom1Btn.Text), "▼ Domain #1");
      testCase.verifyEqual(string(dom2Btn.Text), "▼ Domain #2");

      greenPanel = findobj(dom1Btn.Parent, 'Type', 'uipanel', '-depth', 1);
      testCase.verifyEqual(greenPanel.BackgroundColor, [0.91, 0.97, 0.92], 'AbsTol', 1e-6);

      dom1Body = greenPanel.Children(1);
      geomBtn = findSectionButton(dom1Body, "Geometry (Mesh)");
      bcSecBtn = findSectionButton(dom1Body, "Boundary Conditions");
      matBtn = findSectionButton(dom1Body, "Materials");
      physBtn = findSectionButton(dom1Body, "Physics Solver");
      testCase.verifyEqual(string(physBtn.Text), "▼ Physics Solver *");

      % Verify + Add buttons sit underneath their respective item lists
      addDomBtn = findobj(app.Figure, 'Type', 'uibutton', 'Text', '+ Add Domain');
      testCase.verifyEqual(addDomBtn(1).Parent.Layout.Row, 2);
      addIfBtn = findobj(app.Figure, 'Type', 'uibutton', 'Text', '+ Add Interface');
      testCase.verifyEqual(addIfBtn(1).Layout.Row, 2);
      addBcBtn = findobj(dom1Body, 'Type', 'uibutton', 'Text', '+ Add Boundary Condition');
      testCase.verifyEqual(addBcBtn(1).Layout.Row, 2);
      addSolidBtn = findobj(dom1Body, 'Type', 'uibutton', 'Text', '+ Add Solid Material');
      testCase.verifyEqual(addSolidBtn(1).Layout.Row, 4);
      addEvBtn = findobj(dom1Body, 'Type', 'uibutton', 'Text', '+ Add BCevent (time / value)');
      testCase.verifyEqual(addEvBtn(1).Layout.Row, 4);

      % Verify inside-domain order: Geometry (2) -> Boundary Conditions (3) -> Materials (4) -> Physics Solver (5)
      testCase.verifyEqual(geomBtn.Parent.Layout.Row, 2);
      testCase.verifyEqual(bcSecBtn.Parent.Layout.Row, 3);
      testCase.verifyEqual(matBtn.Parent.Layout.Row, 4);
      testCase.verifyEqual(physBtn.Parent.Layout.Row, 5);

      yellowPanel = findobj(geomBtn.Parent, 'Type', 'uipanel', '-depth', 1);
      testCase.verifyEqual(yellowPanel.BackgroundColor, [1.00, 0.97, 0.88], 'AbsTol', 1e-6);

      bc1Lbl = findobj(app.Figure, 'Type', 'uilabel', 'Text', 'Boundary Condition #1');
      lavPanel = bc1Lbl(1).Parent.Parent.Parent;
      testCase.verifyEqual(lavPanel.BackgroundColor, [0.95, 0.92, 0.99], 'AbsTol', 1e-6);

      entLbl = findobj(lavPanel, 'Type', 'uilabel', 'Text', 'entityList');
      helpBtns = findobj(entLbl(1).Parent, 'Type', 'uibutton', 'Tag', 'ParamHelpButton');
      testCase.verifyNumElements(helpBtns, 1);
      testCase.verifyEqual(string(helpBtns(1).Text), "?");
      testCase.verifyTrue(contains(string(helpBtns(1).Tooltip), "Bottom"));
      testCase.verifyTrue(contains(string(helpBtns(1).Tooltip), "Top"));
      testCase.verifyClass(helpBtns(1).ButtonPushedFcn, 'function_handle');

      evPanel = findobj(app.Figure, 'Type', 'uipanel', 'Title', 'BCevent #1');
      testCase.verifyEqual(evPanel(1).BackgroundColor, [1.00, 0.93, 0.93], 'AbsTol', 1e-6);

      % Verify footer buttons include multi-XML load, script generation, and RUN
      testCase.verifyNotEmpty(findobj(app.Figure, 'Type', 'uibutton', 'Text', 'RUN'));
      testCase.verifyNotEmpty(findobj(app.Figure, 'Type', 'uibutton', 'Text', 'Load XML(s)...'));
      testCase.verifyNotEmpty(findobj(app.Figure, 'Type', 'uibutton', 'Text', 'Generate XML & Script...'));
    end

    function testConsensusAndCompatibilityWarnings(testCase)
      app = gresUI('Visible', false);
      cleanupApp = onCleanup(@() delete(app));

      % Single domain BiotFixedStressSplit -> auto-selects FixedStressSplit
      app.selectPhysicsSolver("BiotFixedStressSplit", 1);
      testCase.verifyEqual(app.collectConfiguration().solutionScheme.name, "FixedStressSplit");
      [sev, ~, ~] = app.checkSchemeCompatibility();
      testCase.verifyEqual(sev, "ok");

      % Adding a second domain auto-switches consensus to NonLinearImplicit, and flags BiotFixedStressSplit conflict
      d2 = app.addDomain();
      testCase.verifyEqual(app.collectConfiguration().solutionScheme.name, "NonLinearImplicit");
      [sev, ~, ~] = app.checkSchemeCompatibility();
      testCase.verifyEqual(sev, "error");

      % Switching Domain 1 to BiotFullyCoupled & Domain 2 to Poromechanics clears conflict
      app.selectPhysicsSolver("BiotFullyCoupled", 1);
      app.selectPhysicsSolver("Poromechanics", d2);
      [sev, ~, ~] = app.checkSchemeCompatibility();
      testCase.verifyEqual(sev, "ok");

      % Connecting Poromechanics (displacements) and SinglePhaseFlowFVTPFA (pressure) via MeshTying warns
      app.selectPhysicsSolver("Poromechanics", 1);
      app.selectPhysicsSolver("SinglePhaseFlowFVTPFA", d2);
      app.addInterface("MeshTying", struct('masterDomain', 2, 'slaveDomain', 1));
      [sev, ~, ~] = app.checkSchemeCompatibility();
      testCase.verifyEqual(sev, "warning");

      % Remove Domain 2 & Interface, check Sedimentation -> auto-selects EvolvingGrid
      app2 = gresUI('Visible', false);
      cleanupApp2 = onCleanup(@() delete(app2));
      app2.selectPhysicsSolver("Sedimentation", 1);
      testCase.verifyEqual(app2.collectConfiguration().solutionScheme.name, "EvolvingGrid");
      [sev2, ~, ~] = app2.checkSchemeCompatibility();
      testCase.verifyEqual(sev2, "ok");

      % Manual override to NonLinearImplicit warns & suggests EvolvingGrid
      app2.selectSolutionScheme("NonLinearImplicit");
      [sev3, ~, sug3] = app2.checkSchemeCompatibility();
      testCase.verifyEqual(sev3, "warning");
      testCase.verifyEqual(sug3, "EvolvingGrid");
    end

    function testLoadUserXmlAndMultiDomainRoundTrip(testCase)
      userXml = fullfile(gres_root(), 'Tests', 'testUI', 'problem.xml');
      testCase.verifyTrue(isfile(userXml));
      appUser = gresUI('Visible', false);
      cleanupUser = onCleanup(@() delete(appUser));
      appUser.loadXmlFile(userXml);
      cfgUser = appUser.collectConfiguration();
      testCase.verifyEqual(cfgUser.domains(1).solver.name, "Poromechanics");
      testCase.verifyEqual(cfgUser.domains(1).boundaryConditions(1).params.components, "x,y,z");

      tmpDir = tempname;
      mkdir(tmpDir);
      cleanupDir = onCleanup(@() rmdir(tmpDir, 's'));

      app = buildSampleTwoDomainContactApp();
      cleanupApp = onCleanup(@() delete(app));

      uniDir = fullfile(tmpDir, 'unified');
      mkdir(uniDir);
      [uFiles, uScript] = app.exportToPath(fullfile(uniDir, 'problem.xml'), 'unified');
      testCase.verifyTrue(isfile(uFiles(1)));
      testCase.verifyTrue(isfile(uScript));
      scriptTxt = fileread(char(uScript));
      testCase.verifyTrue(contains(scriptTxt, '[domains, interfaces] = buildModel(params);'));
      testCase.verifyTrue(contains(scriptTxt, "'interface', interfaces"));

      parsed = readInput(char(uFiles(1)));
      testCase.verifyEqual(numel(parsed.Domain), 2);
      testCase.verifyTrue(isfield(parsed, 'Interface'));
      testCase.verifyTrue(isfield(parsed.Interface, 'SolidMechanicsContact'));

      modDir = fullfile(tmpDir, 'modular');
      mkdir(modDir);
      [mFiles, mScript] = app.exportToPath(modDir, 'modular');
      testCase.verifyTrue(isfile(mScript));
      testCase.verifyTrue(any(endsWith(mFiles, "simparams.xml")));
      testCase.verifyTrue(any(endsWith(mFiles, "domain1.xml")));
      testCase.verifyTrue(any(endsWith(mFiles, "domain2.xml")));
      testCase.verifyTrue(any(endsWith(mFiles, "interfaces.xml")));
      testCase.verifyTrue(any(endsWith(mFiles, "problem.xml")));

      appReload = gresUI('Visible', false);
      cleanupReload = onCleanup(@() delete(appReload));
      appReload.loadXmlFile(fullfile(modDir, 'problem.xml'));
      cfgReload = appReload.collectConfiguration();
      testCase.verifyEqual(numel(cfgReload.domains), 2);
      testCase.verifyEqual(numel(cfgReload.interfaces), 1);
      testCase.verifyEqual(cfgReload.interfaces(1).type, "SolidMechanicsContact");
      testCase.verifyEqual(cfgReload.interfaces(1).coulomb.frictionAngle, 30, 'AbsTol', 1e-12);
      testCase.verifyEqual(cfgReload.interfaces(1).activeSet.tolerances.sliding, 1e-4, 'AbsTol', 1e-12);
    end

    function testLoadMultipleXmlFiles(testCase)
      % Load standalone split XML files together from Tests/MandelBiot/Input
      mbDir = fullfile(gres_root(), 'Tests', 'MandelBiot', 'Input');
      appM = gresUI('Visible', false);
      cleanupM = onCleanup(@() delete(appM));
      appM.loadXmlFiles([ ...
        string(fullfile(mbDir, 'boundaries.xml')), ...
        string(fullfile(mbDir, 'materials.xml')), ...
        string(fullfile(mbDir, 'simparam.xml'))]);
      cfgM = appM.collectConfiguration();
      testCase.verifyEqual(cfgM.simParams.End, 5.01, 'AbsTol', 1e-12);
      testCase.verifyEqual(cfgM.domains(1).materials.Fluid.compressibility, 4.4e-7, 'AbsTol', 1e-12);
      testCase.verifyEqual(cfgM.domains(1).materials.Solids(1).constitutive.name, "Elastic");
      testCase.verifyEqual(cfgM.domains(1).materials.Solids(1).constitutive.params.youngModulus, 1e5, 'AbsTol', 1e-12);
      testCase.verifyEqual(cfgM.domains(1).materials.Solids(1).porousRock.params.porosity, 0.375, 'AbsTol', 1e-12);
      testCase.verifyNumElements(cfgM.domains(1).boundaryConditions, 4);
      testCase.verifyEqual(cfgM.domains(1).boundaryConditions(1).params.type, "dirichlet");
      testCase.verifyEqual(cfgM.domains(1).boundaryConditions(1).params.entityListType, "tags");

      % Export a 2-domain modular folder and reload in two multi-select ways
      tmpDir = tempname;
      mkdir(tmpDir);
      cleanupDir = onCleanup(@() rmdir(tmpDir, 's'));
      appSrc = buildSampleTwoDomainContactApp();
      cleanupSrc = onCleanup(@() delete(appSrc));
      mFiles = appSrc.exportToPath(tmpDir, 'modular');

      % Reload modular components WITHOUT problem.xml
      noProbFiles = mFiles(~endsWith(mFiles, "problem.xml"));
      appNoProb = gresUI('Visible', false);
      cleanupNoProb = onCleanup(@() delete(appNoProb));
      appNoProb.loadXmlFiles(noProbFiles);
      cfgNoProb = appNoProb.collectConfiguration();
      testCase.verifyNumElements(cfgNoProb.domains, 2);
      testCase.verifyNumElements(cfgNoProb.interfaces, 1);
      testCase.verifyEqual(cfgNoProb.interfaces(1).type, "SolidMechanicsContact");
      testCase.verifyEqual(cfgNoProb.output.params.matFileName, "Output/contact");

      % Reload ALL generated files (Cmd+A: problem.xml + sub-files) without duplicates
      appAll = gresUI('Visible', false);
      cleanupAll = onCleanup(@() delete(appAll));
      appAll.loadXmlFiles(mFiles);
      cfgAll = appAll.collectConfiguration();
      testCase.verifyNumElements(cfgAll.domains, 2);
      testCase.verifyNumElements(cfgAll.interfaces, 1);
      testCase.verifyNumElements(cfgAll.domains(1).boundaryConditions, 1);
      testCase.verifyNumElements(cfgAll.domains(2).boundaryConditions, 1);
      testCase.verifyEqual(cfgAll.output.params.matFileName, "Output/contact");

      % Combine two standalone <BoundaryConditions> XML files in a single loadXmlFiles call
      bcFile1 = fullfile(tmpDir, 'extra_bc1.xml');
      bcFile2 = fullfile(tmpDir, 'extra_bc2.xml');
      fid1 = fopen(bcFile1, 'w');
      fprintf(fid1, '<BoundaryConditions><BC name="b1" type="dirichlet" field="surface" variable="displacements" entityListType="tags" entityList="1"><BCevent time="0" value="0"/></BC></BoundaryConditions>');
      fclose(fid1);
      fid2 = fopen(bcFile2, 'w');
      fprintf(fid2, '<BoundaryConditions><BC name="b2" type="neumann" field="surface" variable="displacements" entityListType="tags" entityList="2"><BCevent time="0" value="-5"/></BC></BoundaryConditions>');
      fclose(fid2);
      appBcs = gresUI('Visible', false);
      cleanupBcs = onCleanup(@() delete(appBcs));
      appBcs.loadXmlFiles([string(bcFile1), string(bcFile2)]);
      cfgBcs = appBcs.collectConfiguration();
      testCase.verifyNumElements(cfgBcs.domains(1).boundaryConditions, 2);
      testCase.verifyEqual(cfgBcs.domains(1).boundaryConditions(1).params.name, "b1");
      testCase.verifyEqual(cfgBcs.domains(1).boundaryConditions(2).params.name, "b2");
    end

    function testEndToEndSimulationRun(testCase)
      tmpDir = tempname;
      mkdir(tmpDir);
      cleanupDir = onCleanup(@() rmdir(tmpDir, 's'));

      app = buildSampleFlowApp();
      cleanupApp = onCleanup(@() delete(app));

      solver = app.runSimulation(tmpDir);
      testCase.verifyClass(solver, ?NonLinearImplicit);
      testCase.verifyTrue(isfile(fullfile(tmpDir, 'problem.xml')));
      testCase.verifyTrue(isfile(fullfile(tmpDir, 'runSimulation.m')));
    end

    function testMultiDomainMeshTyingSimulationRun(testCase)
      tmpDir = tempname;
      mkdir(tmpDir);
      cleanupDir = onCleanup(@() rmdir(tmpDir, 's'));

      app = gresUI('Visible', false);
      cleanupApp = onCleanup(@() delete(app));

      app.selectPhysicsSolver("SinglePhaseFlowFEM", 1);
      app.setGeometry('StructuredMesh', struct('NX', 2, 'NY', 2, 'NZ', 2, 'LX', [0 1], 'LY', [0 1], 'LZ', [0 1]), 1);
      app.setFluid(struct('specificWeight', 0, 'compressibility', 4.4e-7, 'dynamicViscosity', 1e-3), 1);
      app.addSolid(struct('name', "rock1", 'cellTags', 1, 'PorousRock', struct('porosity', 0.3, 'permeability', 1e-12)), 1);
      app.addBoundaryCondition(struct('name', "leftP", 'type', "dirichlet", 'field', "surface", ...
        'variable', "pressure", 'entityListType', "tags", 'entityList', "5", ...
        'BCevent', struct('time', 0.0, 'value', 0.0)), 1);

      d2 = app.addDomain();
      app.selectPhysicsSolver("SinglePhaseFlowFEM", d2);
      app.setGeometry('StructuredMesh', struct('NX', 2, 'NY', 2, 'NZ', 2, 'LX', [1 2], 'LY', [0 1], 'LZ', [0 1]), d2);
      app.setFluid(struct('specificWeight', 0, 'compressibility', 4.4e-7, 'dynamicViscosity', 1e-3), d2);
      app.addSolid(struct('name', "rock2", 'cellTags', 1, 'PorousRock', struct('porosity', 0.3, 'permeability', 1e-12)), d2);
      app.addBoundaryCondition(struct('name', "rightP", 'type', "dirichlet", 'field', "surface", ...
        'variable', "pressure", 'entityListType', "tags", 'entityList', "6", ...
        'BCevent', struct('time', 0.0, 'value', 10.0)), d2);

      app.addInterface("MeshTying", struct( ...
        'masterDomain', 2, 'slaveDomain', 1, ...
        'masterSurface', 5, 'slaveSurface', 6, ...
        'multiplierType', "P0", 'variable', "pressure"));

      tmpSeed = [tempname '.xml'];
      fid = fopen(tmpSeed, 'w');
      fprintf(fid, '<SimulationParameters Start="0" End="1.0" DtInit="0.5" DtMin="0.1" DtMax="0.5"/>');
      fclose(fid);
      cleanupSeed = onCleanup(@() delete(tmpSeed));
      app.loadXmlFile(tmpSeed);

      solver = app.runSimulation(tmpDir);
      testCase.verifyClass(solver, ?NonLinearImplicit);
      testCase.verifyEqual(solver.nDom, 2);
      testCase.verifyEqual(solver.nInterf, 1);
    end

  end

end


function app = buildSampleFlowApp()
app = gresUI('Visible', false);
app.selectPhysicsSolver("SinglePhaseFlowFVTPFA");
app.selectSolutionScheme("NonLinearImplicit");

app.setGeometry('StructuredMesh', struct('NX', 2, 'NY', 2, 'NZ', 2, 'LX', [0 1], 'LY', [0 1], 'LZ', [0 1]));

app.setFluid(struct('specificWeight', 9.81, 'compressibility', 4.4e-7, 'dynamicViscosity', 1e-3));
app.addSolid(struct( ...
  'name', "sand", ...
  'cellTags', 1, ...
  'PorousRock', struct('porosity', 0.375, 'permeability', 1e-12)));

app.addBoundaryCondition(struct( ...
  'name', "botP", 'type', "dirichlet", 'field', "surface", ...
  'variable', "pressure", 'entityListType', "tags", 'entityList', "1", ...
  'BCevent', struct('time', 0.0, 'value', 0.0)));
app.addBoundaryCondition(struct( ...
  'name', "topP", 'type', "dirichlet", 'field', "surface", ...
  'variable', "pressure", 'entityListType', "tags", 'entityList', "2", ...
  'BCevent', struct('time', 0.0, 'value', 10.0)));

app.setOutput(struct('printTimes', 1.0, 'matFileName', "Output/results"));

tmpSeed = [tempname '.xml'];
fid = fopen(tmpSeed, 'w');
fprintf(fid, '<SimulationParameters Start="0" End="1.0" DtInit="0.5" DtMin="0.1" DtMax="0.5"/>');
fclose(fid);
cleanupSeed = onCleanup(@() delete(tmpSeed));
app.loadXmlFile(tmpSeed);
end


function app = buildSampleTwoDomainContactApp()
app = gresUI('Visible', false);
app.selectPhysicsSolver("Poromechanics", 1);
app.setGeometry('StructuredMesh', struct('NX', 2, 'NY', 2, 'NZ', 2, 'LX', [0 1], 'LY', [0 1], 'LZ', [0 1]), 1);
app.addSolid(struct('name', "dom1Mat", 'cellTags', 1, 'Constitutive', struct('Elastic', struct('youngModulus', 1e4, 'poissonRatio', 0.25))), 1);
app.addBoundaryCondition(struct('name', "fix1", 'variable', "displacements", 'field', "surface", ...
  'entityListType', "tags", 'entityList', "1", 'type', "dirichlet", 'components', "x,y,z", ...
  'BCevent', struct('time', 0, 'value', 0)), 1);

d2 = app.addDomain();
app.selectPhysicsSolver("Poromechanics", d2);
app.setGeometry('StructuredMesh', struct('NX', 2, 'NY', 2, 'NZ', 2, 'LX', [0 1], 'LY', [0 1], 'LZ', [1 2]), d2);
app.addSolid(struct('name', "dom2Mat", 'cellTags', 1, 'Constitutive', struct('Elastic', struct('youngModulus', 1e4, 'poissonRatio', 0.25))), d2);
app.addBoundaryCondition(struct('name', "load2", 'variable', "displacements", 'field', "surface", ...
  'entityListType', "tags", 'entityList', "2", 'type', "neumann", 'components', "z", ...
  'BCevent', struct('time', 0, 'value', -1)), d2);

app.addInterface("SolidMechanicsContact", struct( ...
  'masterDomain', 2, 'slaveDomain', 1, 'masterSurface', 1, 'slaveSurface', 2, ...
  'multiplierType', "P0", 'augmentationParameter', 1e3, 'isSmooth', 1, ...
  'Quadrature', struct('type', "SegmentBasedQuadrature", 'gaussOrder', 4), ...
  'Coulomb', struct('cohesion', 0.1, 'frictionAngle', 30), ...
  'ActiveSet', struct('resetActiveSet', 0, 'Tolerances', struct('sliding', 1e-4))));

app.setOutput(struct('printTimes', [0.5, 1.0], 'matFileName', "Output/contact"));

tmpSeed = [tempname '.xml'];
fid = fopen(tmpSeed, 'w');
fprintf(fid, '<SimulationParameters Start="0" End="1.0" DtInit="0.5" DtMin="0.1" DtMax="0.5"/>');
fclose(fid);
cleanupSeed = onCleanup(@() delete(tmpSeed));
app.loadXmlFile(tmpSeed);
end


function vis = getRowVis(app, fieldName)
lbls = findobj(app.Figure, 'Type', 'uilabel', 'Text', fieldName);
grid = lbls(1).Parent;
vis = string(grid.Visible);
end


function w = getWidgetForParam(app, fieldName)
lbls = findobj(app.Figure, 'Type', 'uilabel', 'Text', fieldName);
grid = lbls(1).Parent;
kids = grid.Children;
w = kids(arrayfun(@(c) ~isa(c, 'matlab.ui.control.Label') && ~isa(c, 'matlab.ui.control.Button'), kids));
w = w(1);
end


function btn = findSectionButton(root, prefixText)
btns = findobj(root, 'Type', 'uibutton');
for i = 1:numel(btns)
  if contains(string(btns(i).Text), prefixText)
    btn = btns(i);
    return;
  end
end
error('Button containing "%s" not found.', prefixText);
end
