classdef gresUI < handle
%GRESUI Auto-maintaining graphical interface for configuring and running GReS simulations.
%
%   app = gresUI()
%   app = gresUI('Visible', false)
%
%   Automatically discovers SimulationParameters, Geometry, PhysicsSolvers,
%   SolutionSchemes, Materials, BoundaryConditions, Interfaces, OutState, and
%   permitted LinearSolver/Chronos options directly from the GReS codebase.
%   Sections without compulsory (*) fields are exported automatically whenever
%   at least one parameter differs from its default or is explicitly loaded.

  properties (SetAccess = private, GetAccess = public)
    Figure
    Schema
  end

  properties (Access = private)
    RootGrid
    ScrollGrid

    LogoImage
    ProblemSolverSection
    SchemeDropdown
    SchemeConsensusLabel
    SchemeSubSection
    WarningPanel
    WarningLabel
    ApplySuggestedBtn

    ExportModeDropdown
    RunButton
    StatusLabel

    SimSection
    DomainsSection
    DomainsListGrid
    InterfacesSection
    InterfacesListGrid
    OutputSection
    LinearSection
    LinearTopSection
    LinearFlowSection
    LinearMechSection

    % Global parameter control maps
    SimControls = struct()
    SchemeControls = struct()
    OutputControls = struct()
    RestartControls = struct()
    LinearTopControls = struct()
    ChronosFlowControls = struct()
    ChronosMechControls = struct()

    % Output dynamic containers
    OutputEnableBox
    OutputRowsGrid
    RestartEnableBox
    RestartSection
    RestartRowsGrid

    % Dynamic domain & interface cards
    DomainCards = {}
    InterfaceCards = {}

    CurrentSchemeName (1,1) string = ""
    SuggestedScheme (1,1) string = ""
    UserOverrodeScheme (1,1) logical = false
    NextCardId (1,1) double = 1
  end

  methods (Access = public)

    function obj = gresUI(opts)
      arguments
        opts.Visible (1,1) logical = true
      end

      obj.Schema = scanGReSParameters();
      obj.buildUI(opts.Visible);
    end

    function delete(obj)
      if ~isempty(obj.Figure) && isvalid(obj.Figure)
        delete(obj.Figure);
      end
    end

    function domIdx = addDomain(obj)
      obj.appendDomainCard();
      domIdx = numel(obj.DomainCards);
    end

    function removeDomain(obj, domIdx)
      if domIdx < 1 || domIdx > numel(obj.DomainCards) || isscalar(obj.DomainCards)
        return;
      end
      obj.removeDomainById(obj.DomainCards{domIdx}.Id);
    end

    function selectPhysicsSolver(obj, solverName, domIdx)
      if nargin < 3
        domIdx = 1;
      end
      domCard = obj.DomainCards{domIdx};
      domCard.SolverDropdown.Value = char(string(solverName));
      obj.onDomainSolverChanged(domCard.Id);
    end

    function selectSolutionScheme(obj, schemeName)
      obj.SchemeDropdown.Value = char(string(schemeName));
      obj.onUserSchemeChanged();
    end

    function setGeometry(obj, geomType, paramsStruct, domIdx)
      if nargin < 4
        domIdx = 1;
      end
      domId = obj.DomainCards{domIdx}.Id;
      obj.DomainCards{domIdx}.GeomEnableBox.Value = true;
      obj.onDomainGeomEnableChanged(domId);
      obj.DomainCards{domIdx}.GeomModeDropdown.Value = char(string(geomType));
      obj.onDomainGeomModeChanged(domId);
      obj.DomainCards{domIdx}.GeomControls = applyStructToGroup(obj.DomainCards{domIdx}.GeomControls, paramsStruct);
    end

    function setFluid(obj, fluidStruct, domIdx)
      if nargin < 3
        domIdx = 1;
      end
      domId = obj.DomainCards{domIdx}.Id;
      obj.DomainCards{domIdx}.FluidEnableBox.Value = true;
      obj.onDomainFluidEnableChanged(domId);
      obj.DomainCards{domIdx}.FluidControls = applyStructToGroup(obj.DomainCards{domIdx}.FluidControls, fluidStruct);
    end

    function addSolid(obj, solidStruct, domIdx)
      if nargin < 3
        domIdx = 1;
      end
      if nargin < 2
        solidStruct = struct();
      end
      domId = obj.DomainCards{domIdx}.Id;
      sId = obj.appendSolidCard(domId);
      [~, sIdx] = obj.findSolidCard(domId, sId);
      sCard = obj.DomainCards{domIdx}.SolidCards{sIdx};
      sCard.HeaderControls = applyStructToGroup(sCard.HeaderControls, solidStruct);
      obj.DomainCards{domIdx}.SolidCards{sIdx} = sCard;

      if isfield(solidStruct, 'Constitutive') && isstruct(solidStruct.Constitutive)
        laws = fieldnames(solidStruct.Constitutive);
        if ~isempty(laws)
          sCard.LawDropdown.Value = laws{1};
          obj.onSolidLawChanged(domId, sId);
          sCard = obj.DomainCards{domIdx}.SolidCards{sIdx};
          sCard.LawControls = applyStructToGroup(sCard.LawControls, solidStruct.Constitutive.(laws{1}));
          obj.DomainCards{domIdx}.SolidCards{sIdx} = sCard;
        end
      end

      if isfield(solidStruct, 'PorousRock') && isstruct(solidStruct.PorousRock)
        sCard = obj.DomainCards{domIdx}.SolidCards{sIdx};
        sCard.PorousEnableBox.Value = true;
        obj.onSolidPorousChanged(domId, sId);
        sCard = obj.DomainCards{domIdx}.SolidCards{sIdx};
        sCard.PorousControls = applyStructToGroup(sCard.PorousControls, solidStruct.PorousRock);
        obj.DomainCards{domIdx}.SolidCards{sIdx} = sCard;
        if isfield(solidStruct.PorousRock, 'Curves') && isstruct(solidStruct.PorousRock.Curves)
          sCard.CurveEnableBox.Value = true;
          obj.onSolidCurveEnableChanged(domId, sId);
          if isfield(solidStruct.PorousRock.Curves, 'type')
            sCard.CurveTypeDropdown.Value = char(string(solidStruct.PorousRock.Curves.type));
            obj.onSolidCurveTypeChanged(domId, sId);
          end
          sCard = obj.DomainCards{domIdx}.SolidCards{sIdx};
          sCard.CurveControls = applyStructToGroup(sCard.CurveControls, solidStruct.PorousRock.Curves);
          obj.DomainCards{domIdx}.SolidCards{sIdx} = sCard;
        end
      end
    end

    function addBoundaryCondition(obj, bcStruct, domIdx)
      if nargin < 3
        domIdx = 1;
      end
      if nargin < 2
        bcStruct = struct();
      end
      domId = obj.DomainCards{domIdx}.Id;
      bcId = obj.appendBCCard(domId);
      [~, bIdx] = obj.findBCCard(domId, bcId);
      obj.DomainCards{domIdx}.BCCards{bIdx}.Controls = applyStructToGroup( ...
        obj.DomainCards{domIdx}.BCCards{bIdx}.Controls, bcStruct);
      obj.onBCTypeChanged(domId, bcId);

      if isfield(bcStruct, 'BCevent') && ~isempty(bcStruct.BCevent)
        evs = bcStruct.BCevent;
        for i = 1:numel(evs)
          obj.appendBCEventRow(domId, bcId, evs(i));
        end
      end
    end

    function ifaceIdx = addInterface(obj, interfType, interfStruct)
      if nargin < 2 || strlength(string(interfType)) == 0
        interfType = obj.Schema.interfaces.models(1).name;
      end
      if nargin < 3
        interfStruct = struct();
      end
      ifaceId = obj.appendInterfaceCard();
      idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
      obj.InterfaceCards{idx}.ModelDropdown.Value = char(string(interfType));
      obj.onInterfaceModelChanged(ifaceId);

      idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
      card = obj.InterfaceCards{idx};
      card.CommonControls = applyStructToGroup(card.CommonControls, interfStruct);
      card.ModelControls = applyStructToGroup(card.ModelControls, interfStruct);
      obj.InterfaceCards{idx} = card;

      if isfield(interfStruct, 'Quadrature') && isstruct(interfStruct.Quadrature) && card.ModelEntry.supportsQuadrature
        card.QuadEnableBox.Value = true;
        obj.onInterfaceQuadEnableChanged(ifaceId);
        if isfield(interfStruct.Quadrature, 'type')
          card = obj.InterfaceCards{obj.findCardIndex(obj.InterfaceCards, ifaceId)};
          card.QuadTypeDropdown.Value = char(string(interfStruct.Quadrature.type));
          obj.onInterfaceQuadTypeChanged(ifaceId);
        end
        idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
        obj.InterfaceCards{idx}.QuadControls = applyStructToGroup(obj.InterfaceCards{idx}.QuadControls, interfStruct.Quadrature);
      end

      if isfield(interfStruct, 'Coulomb') && isstruct(interfStruct.Coulomb)
        idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
        obj.InterfaceCards{idx}.CoulombControls = applyStructToGroup(obj.InterfaceCards{idx}.CoulombControls, interfStruct.Coulomb);
      end

      if isfield(interfStruct, 'ActiveSet') && isstruct(interfStruct.ActiveSet) && ~isempty(card.ModelEntry.activeSetParams)
        card = obj.InterfaceCards{obj.findCardIndex(obj.InterfaceCards, ifaceId)};
        card.ActiveSetEnableBox.Value = true;
        obj.onInterfaceActiveSetChanged(ifaceId);
        idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
        obj.InterfaceCards{idx}.ActiveSetControls = applyStructToGroup(obj.InterfaceCards{idx}.ActiveSetControls, interfStruct.ActiveSet);
        if isfield(interfStruct.ActiveSet, 'Tolerances') && isstruct(interfStruct.ActiveSet.Tolerances)
          card = obj.InterfaceCards{idx};
          card.TolEnableBox.Value = true;
          obj.onInterfaceTolChanged(ifaceId);
          idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
          obj.InterfaceCards{idx}.TolControls = applyStructToGroup(obj.InterfaceCards{idx}.TolControls, interfStruct.ActiveSet.Tolerances);
        end
      end

      obj.onInterfaceDomainsOrModelUpdated(ifaceId);
      ifaceIdx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
    end

    function setOutput(obj, outStruct)
      obj.OutputEnableBox.Value = true;
      obj.onOutputEnableChanged();
      obj.OutputControls = applyStructToGroup(obj.OutputControls, outStruct);
      if isfield(outStruct, 'restart') && isstruct(outStruct.restart)
        obj.RestartEnableBox.Value = true;
        obj.onRestartEnableChanged();
        obj.RestartControls = applyStructToGroup(obj.RestartControls, outStruct.restart);
      end
    end

    function consensus = computeConsensusScheme(obj)
      % Determine the Problem Solver (SolutionScheme) agreed upon by all domains & interfaces.
      consensus = "NonLinearImplicit";
      if isscalar(obj.DomainCards) && isempty(obj.InterfaceCards)
        sName = obj.DomainCards{1}.CurrentSolverName;
        if sName == "BiotFixedStressSplit" && any([obj.Schema.solutionSchemes.name] == "FixedStressSplit")
          consensus = "FixedStressSplit";
        elseif sName == "Sedimentation" && any([obj.Schema.solutionSchemes.name] == "EvolvingGrid")
          consensus = "EvolvingGrid";
        end
      end
    end

    function [severity, message, suggestedScheme] = checkSchemeCompatibility(obj)
      % Check compatibility among SolutionScheme, domain PhysicsSolvers, and Interfaces.
      severity = "ok";
      message = "";
      suggestedScheme = "";

      nDom = numel(obj.DomainCards);
      nInterf = numel(obj.InterfaceCards);
      scheme = obj.CurrentSchemeName;
      if nDom == 0 || strlength(scheme) == 0
        return;
      end

      solverNames = strings(1, nDom);
      for d = 1:nDom
        solverNames(d) = obj.DomainCards{d}.CurrentSolverName;
      end

      if scheme == "FixedStressSplit"
        if nDom > 1 || nInterf > 0
          severity = "error";
          suggestedScheme = "NonLinearImplicit";
          message = "FixedStressSplit supports only 1 domain and 0 interfaces. Switch Problem Solver to NonLinearImplicit.";
          return;
        elseif solverNames(1) ~= "BiotFixedStressSplit"
          severity = "error";
          suggestedScheme = obj.computeConsensusScheme();
          message = sprintf("FixedStressSplit requires BiotFixedStressSplit physics solver (currently '%s').", solverNames(1));
          return;
        end
      elseif any(solverNames == "BiotFixedStressSplit")
        severity = "error";
        if nDom == 1 && nInterf == 0
          suggestedScheme = "FixedStressSplit";
          message = "BiotFixedStressSplit requires FixedStressSplit Problem Solver (or switch domain solver to BiotFullyCoupled).";
        else
          message = "BiotFixedStressSplit only works with FixedStressSplit on 1 domain and 0 interfaces (switch solver to BiotFullyCoupled).";
        end
        return;
      elseif scheme == "EvolvingGrid"
        if nDom > 1 || nInterf > 0
          severity = "error";
          suggestedScheme = "NonLinearImplicit";
          message = "EvolvingGrid supports only 1 domain and 0 interfaces. Switch Problem Solver to NonLinearImplicit.";
          return;
        elseif solverNames(1) ~= "Sedimentation"
          severity = "warning";
          suggestedScheme = "NonLinearImplicit";
          message = sprintf("EvolvingGrid is designed for Sedimentation (currently '%s'). Recommended: NonLinearImplicit.", solverNames(1));
          return;
        end
      elseif any(solverNames == "Sedimentation")
        if nDom > 1 || nInterf > 0
          severity = "error";
          message = "Sedimentation only supports a single domain without interfaces.";
        else
          severity = "warning";
          suggestedScheme = "EvolvingGrid";
          message = "Sedimentation is recommended with EvolvingGrid Problem Solver.";
        end
        return;
      end

      for i = 1:nInterf
        ic = obj.InterfaceCards{i};
        mDom = str2double(string(ic.CommonControls.masterDomain.widget.Value));
        sDom = str2double(string(ic.CommonControls.slaveDomain.widget.Value));
        reqFld = ic.ModelEntry.requiredField;

        if strlength(reqFld) > 0
          for dIdx = unique([mDom, sDom])
            if dIdx >= 1 && dIdx <= nDom
              flds = obj.getDomainActiveFields(obj.DomainCards{dIdx});
              if ~any(flds == reqFld)
                severity = "warning";
                message = sprintf("Interface #%d (%s) requires '%s', which is missing in Domain #%d (%s).", ...
                  i, ic.ModelEntry.name, reqFld, dIdx, obj.DomainCards{dIdx}.CurrentSolverName);
                return;
              end
            end
          end
        elseif mDom >= 1 && mDom <= nDom && sDom >= 1 && sDom <= nDom
          fM = obj.getDomainActiveFields(obj.DomainCards{mDom});
          fS = obj.getDomainActiveFields(obj.DomainCards{sDom});
          if isempty(intersect(fM, fS))
            severity = "warning";
            message = sprintf("Interface #%d (%s) connects Domain #%d and Domain #%d with no shared field.", ...
              i, ic.ModelEntry.name, mDom, sDom);
            return;
          end
        end
      end
    end

    function [isValid, missingList] = validateInputs(obj, opts)
      arguments
        obj
        opts.RequireGeometry (1,1) logical = false
      end

      missingList = string.empty;

      [okSim, missSim] = validateControlGroup(obj.SimControls);
      if ~okSim
        missingList = [missingList, "SimulationParameters." + missSim];
        obj.setSectionCollapsed(obj.SimSection, false);
      end

      nDom = numel(obj.DomainCards);
      for d = 1:nDom
        dc = obj.DomainCards{d};
        domPrefix = sprintf('Domain#%d.', d);
        expandDomain = false;

        if opts.RequireGeometry && ~dc.GeomEnableBox.Value
          missingList(end+1) = string([domPrefix 'Geometry (required to run simulation)']); %#ok<AGROW>
          obj.setSectionCollapsed(dc.GeomSection, false);
          expandDomain = true;
        elseif dc.GeomEnableBox.Value
          [okG, missG] = validateControlGroup(dc.GeomControls);
          if ~okG
            missingList = [missingList, domPrefix + "Geometry." + missG]; %#ok<AGROW>
            obj.setSectionCollapsed(dc.GeomSection, false);
            expandDomain = true;
          end
        end

        for b = 1:numel(dc.BCCards)
          bc = dc.BCCards{b};
          [okB, missB] = validateControlGroup(bc.Controls);
          if ~okB
            missingList = [missingList, sprintf('%sBC#%d.%s', domPrefix, b, strjoin(missB, ','))]; %#ok<AGROW>
            obj.setSectionCollapsed(dc.BoundariesSection, false);
            expandDomain = true;
          end
          for e = 1:numel(bc.Events)
            [okE, missE] = validateControlGroup(bc.Events{e}.Controls);
            if ~okE
              missingList = [missingList, sprintf('%sBC#%d.Event#%d.%s', domPrefix, b, e, strjoin(missE, ','))]; %#ok<AGROW>
              obj.setSectionCollapsed(dc.BoundariesSection, false);
              expandDomain = true;
            end
          end
        end

        for s = 1:numel(dc.SolidCards)
          c = dc.SolidCards{s};
          [okS, missS] = validateControlGroup(c.HeaderControls);
          if ~okS
            missingList = [missingList, sprintf('%sSolid#%d.%s', domPrefix, s, strjoin(missS, ','))]; %#ok<AGROW>
            obj.setSectionCollapsed(dc.MaterialsSection, false);
            expandDomain = true;
          end
          [okL, missL] = validateControlGroup(c.LawControls);
          if ~okL
            missingList = [missingList, sprintf('%sSolid#%d.%s.%s', domPrefix, s, c.LawDropdown.Value, strjoin(missL, ','))]; %#ok<AGROW>
            obj.setSectionCollapsed(dc.MaterialsSection, false);
            expandDomain = true;
          end
          if c.PorousEnableBox.Value
            [okP, missP] = validateControlGroup(c.PorousControls);
            if ~okP
              missingList = [missingList, sprintf('%sSolid#%d.PorousRock.%s', domPrefix, s, strjoin(missP, ','))]; %#ok<AGROW>
              obj.setSectionCollapsed(dc.MaterialsSection, false);
              expandDomain = true;
            end
            if c.CurveEnableBox.Value
              [okC, missC] = validateControlGroup(c.CurveControls);
              if ~okC
                missingList = [missingList, sprintf('%sSolid#%d.Curves.%s', domPrefix, s, strjoin(missC, ','))]; %#ok<AGROW>
                obj.setSectionCollapsed(dc.MaterialsSection, false);
                expandDomain = true;
              end
            end
          end
        end

        [okSol, missSol] = validateControlGroup(dc.SolverControls);
        if ~okSol
          missingList = [missingList, domPrefix + dc.CurrentSolverName + "." + missSol]; %#ok<AGROW>
          obj.setSectionCollapsed(dc.PhysicsSection, false);
          expandDomain = true;
        end

        if ~isempty(dc.CoupledMechSection)
          [okM, missM] = validateControlGroup(dc.CoupledMechControls);
          if ~okM
            missingList = [missingList, domPrefix + dc.CurrentMechSolver + "." + missM]; %#ok<AGROW>
            obj.setSectionCollapsed(dc.PhysicsSection, false);
            obj.setSectionCollapsed(dc.CoupledMechSection, false);
            expandDomain = true;
          end
        end

        if ~isempty(dc.CoupledFlowSection)
          [okF, missF] = validateControlGroup(dc.CoupledFlowControls);
          if ~okF
            missingList = [missingList, domPrefix + dc.CurrentFlowSolver + "." + missF]; %#ok<AGROW>
            obj.setSectionCollapsed(dc.PhysicsSection, false);
            obj.setSectionCollapsed(dc.CoupledFlowSection, false);
            expandDomain = true;
          end
        end

        if expandDomain
          obj.setSectionCollapsed(obj.DomainsSection, false);
          obj.setSectionCollapsed(dc.CardSection, false);
        end
      end

      for i = 1:numel(obj.InterfaceCards)
        ic = obj.InterfaceCards{i};
        [okIc, missIc] = validateControlGroup(ic.CommonControls);
        [okIm, missIm] = validateControlGroup(ic.ModelControls);
        [okCo, missCo] = validateControlGroup(ic.CoulombControls);
        missAll = [missIc, missIm, missCo];
        if ic.QuadEnableBox.Value
          [~, missQ] = validateControlGroup(ic.QuadControls);
          missAll = [missAll, missQ]; %#ok<AGROW>
        end
        if ~isempty(ic.ActiveSetEnableBox) && ic.ActiveSetEnableBox.Value
          [~, missA] = validateControlGroup(ic.ActiveSetControls);
          missAll = [missAll, missA]; %#ok<AGROW>
          if ~isempty(ic.TolEnableBox) && ic.TolEnableBox.Value
            [~, missT] = validateControlGroup(ic.TolControls);
            missAll = [missAll, missT]; %#ok<AGROW>
          end
        end
        if ~okIc || ~okIm || ~okCo || ~isempty(missAll)
          missingList = [missingList, sprintf('Interface#%d.%s', i, strjoin(missAll, ','))]; %#ok<AGROW>
          obj.setSectionCollapsed(obj.InterfacesSection, false);
        end
      end

      if obj.OutputEnableBox.Value && obj.RestartEnableBox.Value
        [okR, missR] = validateControlGroup(obj.RestartControls);
        if ~okR
          missingList = [missingList, "Output.restart." + missR];
          obj.setSectionCollapsed(obj.OutputSection, false);
          obj.setSectionCollapsed(obj.RestartSection, false);
        end
      end

      [okSch, missSch] = validateControlGroup(obj.SchemeControls);
      if ~okSch
        missingList = [missingList, obj.CurrentSchemeName + "." + missSch];
        obj.setSectionCollapsed(obj.ProblemSolverSection, false);
        obj.setSectionCollapsed(obj.SchemeSubSection, false);
      end

      isValid = isempty(missingList);
      [compSev, compMsg, ~] = obj.checkSchemeCompatibility();
      if ~isValid
        obj.setStatus("Missing/invalid compulsory (*) fields: " + strjoin(missingList, ", "), [0.75 0.15 0.15]);
      elseif compSev ~= "ok"
        obj.setStatus("Fields valid. Compatibility warning: " + compMsg, [0.75 0.45 0.05]);
      else
        obj.setStatus("Validation passed. Ready to export or run.", [0.1 0.5 0.2]);
      end
    end

    function cfg = collectConfiguration(obj)
      cfg = struct();
      cfg.simParams = collectGroupValues(obj.SimControls, false);

      cfg.solutionScheme = struct( ...
        'name', obj.CurrentSchemeName, ...
        'params', collectGroupValues(obj.SchemeControls, false));

      domainsArr = struct('solver', {}, 'geometry', {}, 'materials', {}, 'boundaryConditions', {});
      for d = 1:numel(obj.DomainCards)
        dc = obj.DomainCards{d};
        dEntry = struct('solver', struct(), 'geometry', [], 'materials', [], 'boundaryConditions', []);

        sol = struct();
        sol.name = dc.CurrentSolverName;
        sol.params = collectGroupValues(dc.SolverControls, false);
        if obj.isDomainCoupled(dc)
          sol.Mechanics = collectGroupValues(dc.CoupledMechControls, false);
          sol.Flow = collectGroupValues(dc.CoupledFlowControls, false);
        end
        dEntry.solver = sol;

        if dc.GeomEnableBox.Value
          dEntry.geometry = struct( ...
            'type', string(dc.GeomModeDropdown.Value), ...
            'params', collectGroupValues(dc.GeomControls, true));
        end

        mat = struct();
        hasMat = false;
        if dc.FluidEnableBox.Value
          mat.Fluid = collectGroupValues(dc.FluidControls, true);
          hasMat = true;
        end
        if ~isempty(dc.SolidCards)
          solidsArr = struct('params', {}, 'constitutive', {}, 'porousRock', {});
          for s = 1:numel(dc.SolidCards)
            c = dc.SolidCards{s};
            sEntry = struct('params', collectGroupValues(c.HeaderControls, true), ...
              'constitutive', [], 'porousRock', []);
            lawName = string(c.LawDropdown.Value);
            if lawName ~= "None"
              sEntry.constitutive = struct('name', lawName, 'params', collectGroupValues(c.LawControls, true));
            end
            if c.PorousEnableBox.Value
              pr = struct('params', collectGroupValues(c.PorousControls, true), 'curves', []);
              if c.CurveEnableBox.Value
                cv = collectGroupValues(c.CurveControls, true);
                cv.type = string(c.CurveTypeDropdown.Value);
                pr.curves = cv;
              end
              sEntry.porousRock = pr;
            end
            solidsArr(end+1) = sEntry; %#ok<AGROW>
          end
          mat.Solids = solidsArr;
          hasMat = true;
        end
        if hasMat
          dEntry.materials = mat;
        end

        if ~isempty(dc.BCCards)
          bcsArr = struct('params', {}, 'events', {});
          for b = 1:numel(dc.BCCards)
            bc = dc.BCCards{b};
            evArr = struct([]);
            for e = 1:numel(bc.Events)
              evVals = collectGroupValues(bc.Events{e}.Controls, true);
              if isempty(evArr)
                evArr = evVals;
              else
                evArr(end+1) = evVals; %#ok<AGROW>
              end
            end
            bcsArr(end+1) = struct('params', collectGroupValues(bc.Controls, true), 'events', evArr); %#ok<AGROW>
          end
          dEntry.boundaryConditions = bcsArr;
        end

        domainsArr(end+1) = dEntry; %#ok<AGROW>
      end
      cfg.domains = domainsArr;

      if ~isempty(obj.InterfaceCards)
        ifacesArr = struct('type', {}, 'params', {}, 'quadrature', {}, 'coulomb', {}, 'activeSet', {});
        for i = 1:numel(obj.InterfaceCards)
          ic = obj.InterfaceCards{i};
          iEntry = struct( ...
            'type', ic.ModelEntry.name, ...
            'params', mergeStructFields( ...
              collectGroupValues(ic.CommonControls, false), ...
              collectGroupValues(ic.ModelControls, false)), ...
            'quadrature', [], ...
            'coulomb', [], ...
            'activeSet', []);

          if ic.ModelEntry.supportsQuadrature && ic.QuadEnableBox.Value
            iEntry.quadrature = struct( ...
              'type', string(ic.QuadTypeDropdown.Value), ...
              'params', collectGroupValues(ic.QuadControls, true));
          end

          if ~isempty(ic.ModelEntry.coulombParams)
            iEntry.coulomb = collectGroupValues(ic.CoulombControls, true);
          end

          if ~isempty(ic.ActiveSetEnableBox) && ic.ActiveSetEnableBox.Value
            asCfg = struct('params', collectGroupValues(ic.ActiveSetControls, true), 'tolerances', []);
            if ~isempty(ic.TolEnableBox) && ic.TolEnableBox.Value
              asCfg.tolerances = collectGroupValues(ic.TolControls, true);
            end
            iEntry.activeSet = asCfg;
          end

          ifacesArr(end+1) = iEntry; %#ok<AGROW>
        end
        cfg.interfaces = ifacesArr;
      end

      if obj.OutputEnableBox.Value
        outCfg = struct('params', collectGroupValues(obj.OutputControls, true));
        if obj.RestartEnableBox.Value
          outCfg.restart = collectGroupValues(obj.RestartControls, true);
        end
        cfg.output = outCfg;
      end

      lin = struct();
      lin.topLevel = collectGroupValues(obj.LinearTopControls, false);
      if obj.isLinearPhysicsVisible(obj.LinearFlowSection)
        lin.Flow = collectChronosValues(obj.ChronosFlowControls);
      end
      if obj.isLinearPhysicsVisible(obj.LinearMechSection)
        lin.Mechanics = collectChronosValues(obj.ChronosMechControls);
      end
      cfg.linearSolver = lin;
    end

    function [writtenFiles, scriptPath] = exportToPath(obj, outputPath, mode)
      [isValid, missingList] = obj.validateInputs('RequireGeometry', false);
      if ~isValid
        error('gresUI:ValidationFailed', ...
          'Cannot export XML. Missing/invalid compulsory fields: %s', strjoin(missingList, ', '));
      end

      cfg = obj.collectConfiguration();
      writtenFiles = exportGReSXML(cfg, outputPath, mode);

      if strcmpi(string(mode), "unified")
        [outDir, nm, ext] = fileparts(char(outputPath));
        if isempty(outDir)
          outDir = pwd;
        end
        probName = string([nm ext]);
      else
        if isfolder(outputPath)
          outDir = char(outputPath);
        else
          outDir = fileparts(char(outputPath));
        end
        probName = "problem.xml";
      end

      scriptPath = exportGReSScript(cfg, outDir, probName);
      [compSev, compMsg, ~] = obj.checkSchemeCompatibility();
      if compSev ~= "ok"
        obj.setStatus("Generated XML & runSimulation.m (Warning: " + compMsg + ")", [0.75 0.45 0.05]);
      else
        obj.setStatus("Generated XML & runSimulation.m in: " + string(outDir), [0.1 0.5 0.2]);
      end
    end

    function solver = runSimulation(obj, targetDir)
      [isValid, missingList] = obj.validateInputs('RequireGeometry', true);
      if ~isValid
        error('gresUI:ValidationFailed', ...
          'Cannot run simulation. Missing/invalid compulsory fields: %s', strjoin(missingList, ', '));
      end

      if ~isfolder(targetDir)
        mkdir(targetDir);
      end

      mode = string(obj.ExportModeDropdown.Value);
      if mode == "unified"
        xmlTarget = fullfile(char(targetDir), 'problem.xml');
      else
        xmlTarget = char(targetDir);
      end
      [~, ~] = obj.exportToPath(xmlTarget, mode);

      obj.RunButton.Enable = 'off';
      obj.setStatus("Running simulation in " + string(targetDir) + "...", [0.1 0.35 0.65]);
      drawnow;

      oldDir = pwd;
      cleanupDir = onCleanup(@() obj.finishRunCleanup(oldDir));
      cd(char(targetDir));

      tStart = tic;
      params = readInput('problem.xml');
      simparams = SimulationParameters(params.SimulationParameters);
      hasInterf = isfield(params, 'Interface') && isstruct(params.Interface) && ~isempty(fieldnames(params.Interface));
      if hasInterf
        [domains, interfaces] = buildModel(params);
        schemeArgs = {'simulationparameters', simparams, 'domains', domains, 'interface', interfaces};
      else
        domains = buildModel(params);
        schemeArgs = {'simulationparameters', simparams, 'domains', domains};
      end

      if isfield(params, 'Output')
        ensureOutputParentDir(params.Output);
        out = OutState(params.Output);
        schemeArgs = [schemeArgs, {'output', out}];
      end

      cfg = obj.collectConfiguration();
      if isfield(cfg.solutionScheme, 'params') && isstruct(cfg.solutionScheme.params)
        fn = fieldnames(cfg.solutionScheme.params);
        for i = 1:numel(fn)
          schemeArgs = [schemeArgs, {fn{i}, cfg.solutionScheme.params.(fn{i})}]; %#ok<AGROW>
        end
      end

      solver = feval(char(cfg.solutionScheme.name), schemeArgs{:});
      solver.simulationLoop();

      elapsed = toc(tStart);
      obj.setStatus(sprintf('Simulation completed in %.2f s.', elapsed), [0.1 0.5 0.2]);
      clear cleanupDir;
    end

    function loadXmlFile(obj, xmlPath)
      obj.loadXmlFiles(xmlPath);
    end

    function loadXmlFiles(obj, xmlPaths)
      if ischar(xmlPaths)
        paths = string(xmlPaths);
      else
        paths = reshape(string(xmlPaths), 1, []);
      end
      paths = strtrim(paths);
      paths = paths(strlength(paths) > 0);
      if isempty(paths)
        return;
      end

      nF = numel(paths);
      canonPaths = strings(1, nF);
      rootTags = strings(1, nF);
      priorities = zeros(1, nF);
      for i = 1:nF
        canonPaths(i) = canonicalFilePath(paths(i));
        rootTags(i) = readXmlRootTag(paths(i));
        priorities(i) = getXmlLoadPriority(rootTags(i));
      end

      [~, order] = sort(priorities);
      paths = paths(order);
      canonPaths = canonPaths(order);
      rootTags = rootTags(order);

      visited = string.empty;
      clearedMatDom1 = false;
      clearedBcDom1 = false;
      clearedInterfaces = false;

      % Load <Problem> files first (recording any resolved fileName references as visited)
      for i = 1:nF
        if rootTags(i) ~= "Problem" || any(visited == canonPaths(i))
          continue;
        end
        visited(end+1) = canonPaths(i); %#ok<AGROW>
        [baseDir, ~, ~] = fileparts(char(paths(i)));
        data = readInput(char(paths(i)));

        if isfield(data, 'SimulationParameters')
          [spNode, ~, visited] = resolveXmlNodeWithDir(data.SimulationParameters, baseDir, visited);
          obj.applySimParamsStruct(spNode);
        end

        if isfield(data, 'Domain') && isstruct(data.Domain)
          obj.UserOverrodeScheme = false;
          rawDoms = data.Domain;
          nD = numel(rawDoms);
          obj.resetDomainsToCount(nD);
          clearedMatDom1 = true;
          clearedBcDom1 = true;
          for d = 1:nD
            [domNode, domDir, visited] = resolveXmlNodeWithDir(rawDoms(d), baseDir, visited);
            visited = obj.applyDomainStruct(domNode, domDir, d, visited);
          end
        else
          if isfield(data, 'Solver')
            [solNode, ~, visited] = resolveXmlNodeWithDir(data.Solver, baseDir, visited);
            obj.applySolverStruct(solNode, 1);
          end
          if isfield(data, 'Geometry')
            [gNode, ~, visited] = resolveXmlNodeWithDir(data.Geometry, baseDir, visited);
            obj.applyGeometryStruct(gNode, 1);
          end
          if isfield(data, 'Materials')
            if ~clearedMatDom1
              obj.clearDomainMaterials(1);
              clearedMatDom1 = true;
            end
            [mNode, ~, visited] = resolveXmlNodeWithDir(data.Materials, baseDir, visited);
            obj.applyMaterialsStruct(mNode, 1);
          end
          if isfield(data, 'BoundaryConditions')
            if ~clearedBcDom1
              obj.clearDomainBoundaries(1);
              clearedBcDom1 = true;
            end
            [bNode, ~, visited] = resolveXmlNodeWithDir(data.BoundaryConditions, baseDir, visited);
            obj.applyBoundariesStruct(bNode, 1);
          end
        end

        if isfield(data, 'Interface')
          if ~clearedInterfaces
            obj.clearAllInterfaces();
            clearedInterfaces = true;
          end
          [iNode, ~, visited] = resolveXmlNodeWithDir(data.Interface, baseDir, visited);
          obj.applyInterfacesStruct(iNode);
        end

        if isfield(data, 'Output')
          [oNode, ~, visited] = resolveXmlNodeWithDir(data.Output, baseDir, visited);
          obj.setOutput(oNode);
        end
      end

      % Load any unvisited standalone <Domain> files across Domains 1..nD
      domIndices = find(rootTags == "Domain");
      unvisitedDom = domIndices(arrayfun(@(k) ~any(visited == canonPaths(k)), domIndices));
      if ~isempty(unvisitedDom)
        obj.UserOverrodeScheme = false;
        nD = numel(unvisitedDom);
        obj.resetDomainsToCount(nD);
        clearedMatDom1 = true;
        clearedBcDom1 = true;
        for d = 1:nD
          k = unvisitedDom(d);
          visited(end+1) = canonPaths(k); %#ok<AGROW>
          [domDir, ~, ~] = fileparts(char(paths(k)));
          domNode = readInput(char(paths(k)));
          visited = obj.applyDomainStruct(domNode, domDir, d, visited);
        end
      end

      % Load remaining standalone component XML files in priority order
      for i = 1:nF
        if rootTags(i) == "Problem" || rootTags(i) == "Domain" || any(visited == canonPaths(i))
          continue;
        end
        visited(end+1) = canonPaths(i); %#ok<AGROW>
        [baseDir, ~, ~] = fileparts(char(paths(i)));
        data = readInput(char(paths(i)));
        rTag = rootTags(i);

        switch rTag
          case "SimulationParameters"
            spNode = resolveXmlNodeWithDir(data, baseDir);
            obj.applySimParamsStruct(spNode);
          case "LinearSolver"
            lsNode = resolveXmlNodeWithDir(data, baseDir);
            obj.applyLinearSolverStruct(lsNode);
          case "Solver"
            solNode = resolveXmlNodeWithDir(data, baseDir);
            obj.applySolverStruct(solNode, 1);
          case "Geometry"
            gNode = resolveXmlNodeWithDir(data, baseDir);
            obj.applyGeometryStruct(gNode, 1);
          case "Materials"
            if ~clearedMatDom1
              obj.clearDomainMaterials(1);
              clearedMatDom1 = true;
            end
            mNode = resolveXmlNodeWithDir(data, baseDir);
            obj.applyMaterialsStruct(mNode, 1);
          case "BoundaryConditions"
            if ~clearedBcDom1
              obj.clearDomainBoundaries(1);
              clearedBcDom1 = true;
            end
            bNode = resolveXmlNodeWithDir(data, baseDir);
            obj.applyBoundariesStruct(bNode, 1);
          case "Interface"
            if ~clearedInterfaces
              obj.clearAllInterfaces();
              clearedInterfaces = true;
            end
            iNode = resolveXmlNodeWithDir(data, baseDir);
            obj.applyInterfacesStruct(iNode);
          case "Output"
            oNode = resolveXmlNodeWithDir(data, baseDir);
            obj.setOutput(oNode);
          otherwise
            if isfield(data, 'End') || isfield(data, 'DtInit') || isfield(data, 'Start')
              obj.applySimParamsStruct(data);
            end
            if isfield(data, 'LinearSolver')
              obj.applyLinearSolverStruct(data.LinearSolver);
            end
            knownSolvers = [obj.Schema.physicsSolvers.name];
            topFields = string(fieldnames(data));
            if any(ismember(topFields, knownSolvers))
              obj.applySolverStruct(data, 1);
            end
            if isfield(data, 'StructuredMesh') || isfield(data, 'BlockStructuredMesh') || isfield(data, 'meshFile')
              obj.applyGeometryStruct(data, 1);
            end
            if isfield(data, 'Fluid') || isfield(data, 'Solid')
              if ~clearedMatDom1
                obj.clearDomainMaterials(1);
                clearedMatDom1 = true;
              end
              obj.applyMaterialsStruct(data, 1);
            end
            if isfield(data, 'BC')
              if ~clearedBcDom1
                obj.clearDomainBoundaries(1);
                clearedBcDom1 = true;
              end
              obj.applyBoundariesStruct(data, 1);
            end
            knownIfaces = [obj.Schema.interfaces.models.name];
            if any(ismember(topFields, knownIfaces))
              if ~clearedInterfaces
                obj.clearAllInterfaces();
                clearedInterfaces = true;
              end
              obj.applyInterfacesStruct(data);
            end
            if isfield(data, 'outputFile') || isfield(data, 'matFileName') || isfield(data, 'printTimes')
              obj.setOutput(data);
            end
        end
      end

      obj.updateGlobalLinearSolverVisibility();
      obj.refreshInterfaceDomainDropdowns();
      obj.syncConsensusScheme();
      if isscalar(paths)
        obj.setStatus("Loaded XML: " + paths(1), [0.1 0.35 0.65]);
      else
        obj.setStatus(sprintf('Loaded %d XML files.', numel(paths)), [0.1 0.35 0.65]);
      end
    end

  end

  methods (Access = private)

    function finishRunCleanup(obj, oldDir)
      cd(oldDir);
      if ~isempty(obj.RunButton) && isvalid(obj.RunButton)
        obj.RunButton.Enable = 'on';
      end
    end

    function buildUI(obj, isVis)
      visStr = 'off';
      if isVis
        visStr = 'on';
      end

      obj.Figure = uifigure( ...
        'Name', 'GReS UI — Simulation Configurator & Runner', ...
        'Position', [90, 50, 1040, 840], ...
        'Visible', visStr);

      obj.RootGrid = uigridlayout(obj.Figure, [3, 1]);
      obj.RootGrid.RowHeight = {'fit', '1x', 'fit'};
      obj.RootGrid.Padding = [12, 12, 12, 12];
      obj.RootGrid.RowSpacing = 10;

      obj.buildTopBar();

      scrollPanel = uipanel(obj.RootGrid, 'BorderType', 'none');
      scrollPanel.Layout.Row = 2;
      obj.ScrollGrid = uigridlayout(scrollPanel, [6, 1]);
      obj.ScrollGrid.RowHeight = repmat({'fit'}, 1, 6);
      obj.ScrollGrid.RowSpacing = 8;
      obj.ScrollGrid.Padding = [2, 2, 10, 2];
      obj.ScrollGrid.Scrollable = 'on';

      obj.buildSimParamsSection(1);
      obj.buildDomainsSection(2);
      obj.buildInterfacesSection(3);
      obj.buildOutputSection(4);
      obj.buildLinearSolverSection(5);
      obj.buildProblemSolverSection(6);

      obj.buildFooterBar();
      obj.appendDomainCard();
      obj.syncConsensusScheme();
    end

    function buildTopBar(obj)
      topPanel = uipanel(obj.RootGrid, 'BorderType', 'line', 'BackgroundColor', [0.97, 0.98, 1.0]);
      topPanel.Layout.Row = 1;

      tGrid = uigridlayout(topPanel, [1, 2]);
      tGrid.RowHeight = {56};
      tGrid.Padding = [10, 8, 10, 8];
      tGrid.ColumnSpacing = 14;
      tGrid.BackgroundColor = [0.97, 0.98, 1.0];

      logoPath = fullfile(gres_root(), 'logo_GReS.png');
      if isfile(logoPath)
        tGrid.ColumnWidth = {140, '1x'};
        obj.LogoImage = uiimage(tGrid, 'ImageSource', logoPath, 'ScaleMethod', 'fit');
        obj.LogoImage.Layout.Column = 1;
      else
        tGrid.ColumnWidth = {0, '1x'};
      end

      txtGrid = uigridlayout(tGrid, [2, 1]);
      txtGrid.Layout.Column = 2;
      txtGrid.RowHeight = {'fit', 'fit'};
      txtGrid.Padding = [0, 4, 0, 2];
      txtGrid.RowSpacing = 3;
      txtGrid.BackgroundColor = [0.97, 0.98, 1.0];

      titleLbl = uilabel(txtGrid, 'Text', 'GReS Interactive Simulation Configurator', ...
        'FontSize', 16, 'FontWeight', 'bold', 'FontColor', [0.15, 0.28, 0.48]);
      titleLbl.Layout.Row = 1;

      noteLbl = uilabel(txtGrid, ...
        'Text', '* = compulsory parameter. Edit any field to activate its section; untouched optional sections stay omitted from XML.', ...
        'FontSize', 11, 'FontColor', [0.35 0.40 0.48], 'WordWrap', 'on');
      noteLbl.Layout.Row = 2;
    end

    function buildSimParamsSection(obj, rowIdx)
      hasReq = any([obj.Schema.simParams.required]);
      titleStr = "Simulation Parameters";
      if hasReq
        titleStr = titleStr + " *";
      end
      obj.SimSection = obj.createCollapsibleSection(obj.ScrollGrid, rowIdx, titleStr, ~hasReq, 0);
      obj.SimControls = obj.populateParamRows(obj.SimSection.BodyGrid, obj.Schema.simParams, @() [], obj.SimSection.Tint);
    end

    function buildDomainsSection(obj, rowIdx)
      obj.DomainsSection = obj.createCollapsibleSection(obj.ScrollGrid, rowIdx, "Domains *", false, 0);
      dGrid = obj.DomainsSection.BodyGrid;
      dGrid.RowHeight = {'fit', 'fit'};
      dGrid.RowSpacing = 8;

      obj.DomainsListGrid = uigridlayout(dGrid, [1, 1]);
      obj.DomainsListGrid.Layout.Row = 1;
      obj.DomainsListGrid.Padding = [0, 0, 0, 0];
      obj.DomainsListGrid.RowSpacing = 8;
      obj.DomainsListGrid.BackgroundColor = obj.DomainsSection.Tint;

      bar = uigridlayout(dGrid, [1, 2]);
      bar.Layout.Row = 2;
      bar.ColumnWidth = {'fit', '1x'};
      bar.Padding = [0, 0, 0, 0];
      bar.BackgroundColor = obj.DomainsSection.Tint;

      addDomBtn = uibutton(bar, 'Text', '+ Add Domain', 'FontWeight', 'bold', ...
        'ButtonPushedFcn', @(~,~) obj.appendDomainCard());
      addDomBtn.Layout.Column = 1;

      infoLbl = uilabel(bar, 'Text', 'Each domain configures Geometry, Boundary Conditions, Materials, and Physics Solver.', ...
        'FontSize', 11, 'FontColor', [0.3 0.35 0.42]);
      infoLbl.Layout.Column = 2;
    end

    function buildInterfacesSection(obj, rowIdx)
      obj.InterfacesSection = obj.createCollapsibleSection(obj.ScrollGrid, rowIdx, "Interfaces (Mortar / Contact / Mesh Tying)", true, 0);
      iGrid = obj.InterfacesSection.BodyGrid;
      iGrid.RowHeight = {'fit', 'fit'};
      iGrid.RowSpacing = 6;

      obj.InterfacesListGrid = uigridlayout(iGrid, [1, 1]);
      obj.InterfacesListGrid.Layout.Row = 1;
      obj.InterfacesListGrid.RowHeight = {0};
      obj.InterfacesListGrid.Padding = [0, 0, 0, 0];
      obj.InterfacesListGrid.RowSpacing = 8;
      obj.InterfacesListGrid.BackgroundColor = obj.InterfacesSection.Tint;

      addIfBtn = uibutton(iGrid, 'Text', '+ Add Interface', ...
        'ButtonPushedFcn', @(~,~) obj.appendInterfaceCard());
      addIfBtn.Layout.Row = 2;
    end

    function buildOutputSection(obj, rowIdx)
      obj.OutputSection = obj.createCollapsibleSection(obj.ScrollGrid, rowIdx, "Output", true, 0);
      oGrid = obj.OutputSection.BodyGrid;
      oGrid.RowHeight = {'fit', 'fit', 'fit'};
      oGrid.RowSpacing = 6;

      obj.OutputEnableBox = uicheckbox(oGrid, 'Text', 'Enable Output', 'Value', false, ...
        'FontWeight', 'bold', 'ValueChangedFcn', @(~,~) obj.onOutputEnableChanged());
      obj.OutputEnableBox.Layout.Row = 1;

      obj.OutputRowsGrid = uigridlayout(oGrid, [1, 1]);
      obj.OutputRowsGrid.Layout.Row = 2;
      obj.OutputRowsGrid.Padding = [0, 0, 0, 0];
      obj.OutputRowsGrid.BackgroundColor = obj.OutputSection.Tint;
      obj.OutputControls = obj.populateParamRows(obj.OutputRowsGrid, obj.Schema.outputs.params, @() [], obj.OutputSection.Tint);

      obj.RestartSection = obj.createCollapsibleSection(oGrid, 3, "Restart Settings", true, 1);
      rGrid = obj.RestartSection.BodyGrid;
      rGrid.RowHeight = {'fit', 'fit'};
      obj.RestartEnableBox = uicheckbox(rGrid, 'Text', 'Enable Restart', 'Value', false, ...
        'FontWeight', 'bold', 'ValueChangedFcn', @(~,~) obj.onRestartEnableChanged());
      obj.RestartEnableBox.Layout.Row = 1;

      obj.RestartRowsGrid = uigridlayout(rGrid, [1, 1]);
      obj.RestartRowsGrid.Layout.Row = 2;
      obj.RestartRowsGrid.Padding = [0, 0, 0, 0];
      obj.RestartRowsGrid.BackgroundColor = obj.RestartSection.Tint;
      obj.RestartControls = obj.populateParamRows(obj.RestartRowsGrid, obj.Schema.outputs.restartParams, @() [], obj.RestartSection.Tint);

      obj.onOutputEnableChanged();
      obj.onRestartEnableChanged();
    end

    function buildLinearSolverSection(obj, rowIdx)
      obj.LinearSection = obj.createCollapsibleSection(obj.ScrollGrid, rowIdx, "Linear Solver & Chronos Preconditioner", true, 0);
      linGrid = obj.LinearSection.BodyGrid;
      linGrid.RowHeight = {'fit', 'fit', 'fit'};
      linGrid.RowSpacing = 6;

      obj.LinearTopSection = obj.createCollapsibleSection(linGrid, 1, "General Options", true, 1);
      obj.LinearTopControls = obj.populateParamRows(obj.LinearTopSection.BodyGrid, ...
        obj.Schema.linearSolver.topLevel, @() obj.onLinearOptionChanged(), obj.LinearTopSection.Tint);

      obj.LinearFlowSection = obj.createCollapsibleSection(linGrid, 2, "Flow Preconditioner (Chronos AMG)", true, 1);
      obj.ChronosFlowControls = obj.populateChronosSection(obj.LinearFlowSection.BodyGrid, ...
        obj.Schema.linearSolver.chronos.Flow, 2);

      obj.LinearMechSection = obj.createCollapsibleSection(linGrid, 3, "Mechanics Preconditioner (Chronos AMG)", true, 1);
      obj.ChronosMechControls = obj.populateChronosSection(obj.LinearMechSection.BodyGrid, ...
        obj.Schema.linearSolver.chronos.Mechanics, 2);
    end

    function buildProblemSolverSection(obj, rowIdx)
      obj.ProblemSolverSection = obj.createCollapsibleSection(obj.ScrollGrid, rowIdx, ...
        "Problem Solver (SolutionScheme) *", false, 0);
      pGrid = obj.ProblemSolverSection.BodyGrid;
      pGrid.RowHeight = {'fit', 0, 0};
      pGrid.RowSpacing = 6;

      hdr = uigridlayout(pGrid, [1, 3]);
      hdr.Layout.Row = 1;
      hdr.ColumnWidth = {'fit', '1x', 'fit'};
      hdr.Padding = [0, 0, 0, 0];
      hdr.BackgroundColor = obj.ProblemSolverSection.Tint;

      lbl = uilabel(hdr, 'Text', 'Problem Solver (SolutionScheme) *:', 'FontWeight', 'bold');
      lbl.Layout.Column = 1;

      schemeNames = cellstr([obj.Schema.solutionSchemes.name]);
      obj.SchemeDropdown = uidropdown(hdr, ...
        'Items', schemeNames, 'Value', 'NonLinearImplicit', ...
        'ValueChangedFcn', @(~,~) obj.onUserSchemeChanged());
      obj.SchemeDropdown.Layout.Column = 2;

      obj.SchemeConsensusLabel = uilabel(hdr, 'Text', 'Auto-selected from domain physics solvers', ...
        'FontSize', 11, 'FontAngle', 'italic', 'FontColor', [0.2 0.38 0.58]);
      obj.SchemeConsensusLabel.Layout.Column = 3;

      % Compatibility warning banner (hidden when compatible)
      obj.WarningPanel = uipanel(pGrid, 'BorderType', 'line', ...
        'BackgroundColor', [1.0, 0.95, 0.84], 'HighlightColor', [0.85, 0.55, 0.10], ...
        'Visible', 'off');
      obj.WarningPanel.Layout.Row = 2;

      wGrid = uigridlayout(obj.WarningPanel, [1, 2]);
      wGrid.ColumnWidth = {'1x', 'fit'};
      wGrid.Padding = [8, 4, 8, 4];
      wGrid.BackgroundColor = [1.0, 0.95, 0.84];

      obj.WarningLabel = uilabel(wGrid, 'Text', '', 'FontWeight', 'bold', ...
        'FontColor', [0.55, 0.25, 0.0], 'WordWrap', 'on');
      obj.WarningLabel.Layout.Column = 1;

      obj.ApplySuggestedBtn = uibutton(wGrid, 'Text', 'Use Suggested', ...
        'FontWeight', 'bold', 'BackgroundColor', [0.90, 0.65, 0.15], ...
        'FontColor', [1, 1, 1], 'ButtonPushedFcn', @(~,~) obj.onApplySuggestedClicked());
      obj.ApplySuggestedBtn.Layout.Column = 2;

      obj.SchemeSubSection = obj.createCollapsibleSection(pGrid, 3, "Solution Scheme Options", true, 1);
    end

    function buildFooterBar(obj)
      fGrid = uigridlayout(obj.RootGrid, [2, 6]);
      fGrid.Layout.Row = 3;
      fGrid.RowHeight = {'fit', 'fit'};
      fGrid.ColumnWidth = {'fit', 'fit', '1x', 'fit', 'fit', 'fit'};
      fGrid.Padding = [2, 4, 2, 0];

      loadBtn = uibutton(fGrid, 'Text', 'Load XML(s)...', ...
        'ButtonPushedFcn', @(~,~) obj.onLoadClicked());
      loadBtn.Layout.Row = 1;
      loadBtn.Layout.Column = 1;

      valBtn = uibutton(fGrid, 'Text', 'Validate', ...
        'ButtonPushedFcn', @(~,~) obj.validateInputs());
      valBtn.Layout.Row = 1;
      valBtn.Layout.Column = 2;

      obj.ExportModeDropdown = uidropdown(fGrid, ...
        'Items', {'Unified XML + Script', 'Modular XMLs + Script'}, ...
        'ItemsData', {'unified', 'modular'}, ...
        'Value', 'unified');
      obj.ExportModeDropdown.Layout.Row = 1;
      obj.ExportModeDropdown.Layout.Column = 4;

      expBtn = uibutton(fGrid, 'Text', 'Generate XML & Script...', ...
        'ButtonPushedFcn', @(~,~) obj.onExportClicked());
      expBtn.Layout.Row = 1;
      expBtn.Layout.Column = 5;

      obj.RunButton = uibutton(fGrid, 'Text', 'RUN', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.15, 0.55, 0.25], ...
        'FontColor', [1, 1, 1], ...
        'ButtonPushedFcn', @(~,~) obj.onRunClicked());
      obj.RunButton.Layout.Row = 1;
      obj.RunButton.Layout.Column = 6;

      obj.StatusLabel = uilabel(fGrid, 'Text', 'Ready.', 'FontSize', 11);
      obj.StatusLabel.Layout.Row = 2;
      obj.StatusLabel.Layout.Column = [1, 6];
    end

    function domId = appendDomainCard(obj)
      domId = obj.NextCardId;
      obj.NextCardId = obj.NextCardId + 1;
      nD = numel(obj.DomainCards) + 1;

      obj.DomainsListGrid.RowHeight = repmat({'fit'}, 1, nD);
      cardSec = obj.createCollapsibleSection(obj.DomainsListGrid, nD, sprintf('Domain #%d', nD), false, 1);
      bodyGrid = cardSec.BodyGrid;
      bodyGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit'};
      bodyGrid.RowSpacing = 6;

      % Row 1: Domain header bar (summary + red Remove Domain button)
      bar = uigridlayout(bodyGrid, [1, 2]);
      bar.Layout.Row = 1;
      bar.ColumnWidth = {'1x', 'fit'};
      bar.Padding = [0, 0, 0, 0];
      bar.BackgroundColor = cardSec.Tint;

      domNote = uilabel(bar, 'Text', 'Order: Geometry → Boundary Conditions → Materials → Physics Solver', ...
        'FontSize', 11, 'FontAngle', 'italic', 'FontColor', [0.18, 0.42, 0.26]);
      domNote.Layout.Column = 1;

      rmBtn = createDeleteButton(bar, 'Remove Domain', @(~,~) obj.removeDomainById(domId));
      rmBtn.Layout.Column = 2;

      % Row 2: Geometry (Mesh) — Depth 2 Yellow
      geomSec = obj.createCollapsibleSection(bodyGrid, 2, "Geometry (Mesh)", true, 2);
      gGrid = geomSec.BodyGrid;
      gGrid.RowHeight = {'fit', 'fit'};
      gGrid.RowSpacing = 6;

      gHdr = uigridlayout(gGrid, [1, 3]);
      gHdr.Layout.Row = 1;
      gHdr.ColumnWidth = {'fit', 'fit', '1x'};
      gHdr.Padding = [0, 0, 0, 0];
      gHdr.BackgroundColor = geomSec.Tint;

      gEnBox = uicheckbox(gHdr, 'Text', 'Include Geometry', 'Value', false, ...
        'FontWeight', 'bold', 'ValueChangedFcn', @(~,~) obj.onDomainGeomEnableChanged(domId));
      gEnBox.Layout.Column = 1;

      gModeLbl = uilabel(gHdr, 'Text', 'Mesh Type:');
      gModeLbl.Layout.Column = 2;

      geomModes = cellstr([obj.Schema.geometry.name]);
      gModeDrop = uidropdown(gHdr, 'Items', geomModes, 'Value', geomModes{1}, ...
        'ValueChangedFcn', @(~,~) obj.onDomainGeomModeChanged(domId));
      gModeDrop.Layout.Column = 3;

      gRows = uigridlayout(gGrid, [1, 1]);
      gRows.Layout.Row = 2;
      gRows.Padding = [0, 0, 0, 0];
      gRows.BackgroundColor = geomSec.Tint;

      % Row 3: Boundary Conditions — Depth 2 Yellow
      bcSec = obj.createCollapsibleSection(bodyGrid, 3, "Boundary Conditions", true, 2);
      bGrid = bcSec.BodyGrid;
      bGrid.RowHeight = {'fit', 'fit'};
      bGrid.RowSpacing = 6;

      bcListGrid = uigridlayout(bGrid, [1, 1]);
      bcListGrid.Layout.Row = 1;
      bcListGrid.RowHeight = {0};
      bcListGrid.Padding = [0, 0, 0, 0];
      bcListGrid.RowSpacing = 6;
      bcListGrid.BackgroundColor = bcSec.Tint;

      addBcBtn = uibutton(bGrid, 'Text', '+ Add Boundary Condition', ...
        'ButtonPushedFcn', @(~,~) obj.appendBCCard(domId));
      addBcBtn.Layout.Row = 2;

      % Row 4: Materials — Depth 2 Yellow
      matSec = obj.createCollapsibleSection(bodyGrid, 4, "Materials", true, 2);
      mGrid = matSec.BodyGrid;
      mGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
      mGrid.RowSpacing = 6;

      flBox = uicheckbox(mGrid, 'Text', 'Enable Fluid Phase', 'Value', false, ...
        'FontWeight', 'bold', 'ValueChangedFcn', @(~,~) obj.onDomainFluidEnableChanged(domId));
      flBox.Layout.Row = 1;

      flRows = uigridlayout(mGrid, [1, 1]);
      flRows.Layout.Row = 2;
      flRows.Padding = [12, 0, 0, 0];
      flRows.BackgroundColor = matSec.Tint;
      flCtrls = obj.populateParamRows(flRows, obj.Schema.materials.fluidParams, @() [], matSec.Tint);

      solidsGrid = uigridlayout(mGrid, [1, 1]);
      solidsGrid.Layout.Row = 3;
      solidsGrid.RowHeight = {0};
      solidsGrid.Padding = [0, 0, 0, 0];
      solidsGrid.RowSpacing = 6;
      solidsGrid.BackgroundColor = matSec.Tint;

      addSolidBtn = uibutton(mGrid, 'Text', '+ Add Solid Material', ...
        'ButtonPushedFcn', @(~,~) obj.appendSolidCard(domId));
      addSolidBtn.Layout.Row = 4;

      % Row 5: Physics Solver * — Depth 2 Yellow (expanded so compulsory solver choice is visible)
      physSec = obj.createCollapsibleSection(bodyGrid, 5, "Physics Solver *", false, 2);
      pBody = physSec.BodyGrid;
      pBody.RowHeight = {'fit', 'fit'};
      pBody.RowSpacing = 6;

      pBar = uigridlayout(pBody, [1, 3]);
      pBar.Layout.Row = 1;
      pBar.ColumnWidth = {'fit', '1x', 'fit'};
      pBar.Padding = [0, 0, 0, 0];
      pBar.BackgroundColor = physSec.Tint;

      sLbl = uilabel(pBar, 'Text', 'Physics Solver *:', 'FontWeight', 'bold');
      sLbl.Layout.Column = 1;

      solverNames = cellstr([obj.Schema.physicsSolvers.name]);
      sDrop = uidropdown(pBar, 'Items', solverNames, 'Value', solverNames{1}, ...
        'ValueChangedFcn', @(~,~) obj.onDomainSolverChanged(domId));
      sDrop.Layout.Column = 2;

      fldLbl = uilabel(pBar, 'Text', '', 'FontColor', [0.35 0.25 0.05]);
      fldLbl.Layout.Column = 3;

      solverBodyGrid = uigridlayout(pBody, [1, 1]);
      solverBodyGrid.Layout.Row = 2;
      solverBodyGrid.Padding = [0, 0, 0, 0];
      solverBodyGrid.RowSpacing = 6;
      solverBodyGrid.BackgroundColor = physSec.Tint;

      dc = struct( ...
        'Id', domId, ...
        'CardSection', cardSec, ...
        'RemoveBtn', rmBtn, ...
        'GeomSection', geomSec, ...
        'GeomEnableBox', gEnBox, ...
        'GeomModeDropdown', gModeDrop, ...
        'GeomRowsGrid', gRows, ...
        'GeomControls', struct(), ...
        'BoundariesSection', bcSec, ...
        'BCListGrid', bcListGrid, ...
        'BCCards', {{}}, ...
        'MaterialsSection', matSec, ...
        'FluidEnableBox', flBox, ...
        'FluidRowsGrid', flRows, ...
        'FluidControls', flCtrls, ...
        'SolidsListGrid', solidsGrid, ...
        'SolidCards', {{}}, ...
        'PhysicsSection', physSec, ...
        'SolverDropdown', sDrop, ...
        'FieldsLabel', fldLbl, ...
        'SolverBodyGrid', solverBodyGrid, ...
        'SolverControls', struct(), ...
        'CoupledMechSection', [], ...
        'CoupledFlowSection', [], ...
        'CoupledMechControls', struct(), ...
        'CoupledFlowControls', struct(), ...
        'CurrentSolverName', string(solverNames{1}), ...
        'CurrentMechSolver', "", ...
        'CurrentFlowSolver', "");

      obj.DomainCards{end+1} = dc;
      obj.onDomainGeomModeChanged(domId);
      obj.onDomainGeomEnableChanged(domId);
      obj.onDomainFluidEnableChanged(domId);
      obj.onDomainSolverChanged(domId);
    end

    function removeDomainById(obj, domId)
      if isscalar(obj.DomainCards)
        return;
      end
      idx = obj.findCardIndex(obj.DomainCards, domId);
      if isempty(idx)
        return;
      end
      delete(obj.DomainCards{idx}.CardSection.Container);
      obj.DomainCards(idx) = [];
      obj.relayoutDomainCards();
      obj.updateGlobalLinearSolverVisibility();
      obj.refreshInterfaceDomainDropdowns();
      obj.syncConsensusScheme();
    end

    function relayoutDomainCards(obj)
      nD = numel(obj.DomainCards);
      obj.DomainsListGrid.RowHeight = repmat({'fit'}, 1, nD);
      for i = 1:nD
        dc = obj.DomainCards{i};
        dc.CardSection.Container.Layout.Row = i;
        dc.CardSection.BaseTitle = sprintf('Domain #%d', i);
        isCol = startsWith(string(dc.CardSection.HeaderButton.Text), "▶");
        dc.CardSection.Collapsed = isCol;
        obj.applySectionState(dc.CardSection);
        dc.RemoveBtn.Enable = boolOnOff(nD > 1);
        obj.DomainCards{i} = dc;
      end
    end

    function onDomainSolverChanged(obj, domId)
      idx = obj.findCardIndex(obj.DomainCards, domId);
      dc = obj.DomainCards{idx};
      sName = string(dc.SolverDropdown.Value);
      dc.CurrentSolverName = sName;
      sEntry = obj.getSolverEntry(sName);
      dc.FieldsLabel.Text = char("Fields: [" + strjoin(sEntry.fields, ", ") + "]");

      dc = obj.rebuildDomainPhysicsPanel(dc, sEntry);
      obj.DomainCards{idx} = dc;
      obj.refreshDomainBCVariableChoices(domId);
      obj.relayoutDomainCards();
      obj.updateGlobalLinearSolverVisibility();
      obj.refreshInterfaceDomainDropdowns();
      obj.syncConsensusScheme();
    end

    function dc = rebuildDomainPhysicsPanel(obj, dc, sEntry)
      bodyGrid = dc.SolverBodyGrid;
      delete(bodyGrid.Children);
      dc.SolverControls = struct();
      dc.CoupledMechControls = struct();
      dc.CoupledFlowControls = struct();
      dc.CoupledMechSection = [];
      dc.CoupledFlowSection = [];
      dc.CurrentMechSolver = "";
      dc.CurrentFlowSolver = "";

      if sEntry.hasCustomSchema
        bodyGrid.RowHeight = {'fit'};
        msg = uilabel(bodyGrid, ...
          'Text', sprintf('%s uses a custom domain/table XML schema (left to expert/manual XML setup).', sEntry.name), ...
          'FontAngle', 'italic', 'WordWrap', 'on');
        msg.Layout.Row = 1;
        hasReq = false;
      elseif sEntry.isCoupled
        bodyGrid.RowHeight = {'fit', 'fit', 'fit'};
        topParamGrid = uigridlayout(bodyGrid, [1, 1]);
        topParamGrid.Layout.Row = 1;
        topParamGrid.Padding = [0, 0, 0, 0];
        topParamGrid.BackgroundColor = dc.PhysicsSection.Tint;

        domId = dc.Id;
        dc.SolverControls = obj.populateParamRows(topParamGrid, sEntry.params, ...
          @() obj.onDomainCoupledSubSolverChanged(domId), dc.PhysicsSection.Tint);

        dc.CoupledMechSection = obj.createCollapsibleSection(bodyGrid, 2, "Coupled Mechanics Sub-Solver", true, 3);
        dc.CoupledFlowSection = obj.createCollapsibleSection(bodyGrid, 3, "Coupled Flow Sub-Solver", true, 3);
        dc = obj.populateDomainCoupledSubSolvers(dc);
        hasReq = ~isempty(sEntry.params) && any([sEntry.params.required]);
      else
        bodyGrid.RowHeight = {'fit'};
        paramGrid = uigridlayout(bodyGrid, [1, 1]);
        paramGrid.Layout.Row = 1;
        paramGrid.Padding = [0, 0, 0, 0];
        paramGrid.BackgroundColor = dc.PhysicsSection.Tint;
        dc.SolverControls = obj.populateParamRows(paramGrid, sEntry.params, @() [], dc.PhysicsSection.Tint);
        hasReq = ~isempty(sEntry.params) && any([sEntry.params.required]);
      end

      dc.PhysicsSection.BaseTitle = "Physics Solver *";
      isCol = startsWith(string(dc.PhysicsSection.HeaderButton.Text), "▶");
      dc.PhysicsSection.Collapsed = isCol && ~hasReq;
      obj.applySectionState(dc.PhysicsSection);
    end

    function onDomainCoupledSubSolverChanged(obj, domId)
      idx = obj.findCardIndex(obj.DomainCards, domId);
      dc = obj.populateDomainCoupledSubSolvers(obj.DomainCards{idx});
      obj.DomainCards{idx} = dc;
      obj.refreshDomainBCVariableChoices(domId);
      obj.updateGlobalLinearSolverVisibility();
      obj.refreshInterfaceDomainDropdowns();
      obj.syncConsensusScheme();
    end

    function dc = populateDomainCoupledSubSolvers(obj, dc)
      if ~obj.isDomainCoupled(dc)
        return;
      end
      mechName = string(dc.SolverControls.mechanicsSolver.widget.Value);
      flowName = string(dc.SolverControls.flowSolver.widget.Value);
      dc.CurrentMechSolver = mechName;
      dc.CurrentFlowSolver = flowName;

      mechEntry = obj.getSolverEntry(mechName);
      flowEntry = obj.getSolverEntry(flowName);

      dc.CoupledMechSection.BaseTitle = "Mechanics Sub-Solver (" + mechName + ")";
      dc.CoupledMechControls = obj.populateParamRows(dc.CoupledMechSection.BodyGrid, ...
        mechEntry.params, @() [], dc.CoupledMechSection.Tint);
      mechReq = ~isempty(mechEntry.params) && any([mechEntry.params.required]);
      obj.setSectionCollapsed(dc.CoupledMechSection, ~mechReq);

      dc.CoupledFlowSection.BaseTitle = "Flow Sub-Solver (" + flowName + ")";
      dc.CoupledFlowControls = obj.populateParamRows(dc.CoupledFlowSection.BodyGrid, ...
        flowEntry.params, @() [], dc.CoupledFlowSection.Tint);
      flowReq = ~isempty(flowEntry.params) && any([flowEntry.params.required]);
      obj.setSectionCollapsed(dc.CoupledFlowSection, ~flowReq);

      allFields = unique([mechEntry.fields, flowEntry.fields]);
      dc.FieldsLabel.Text = char("Fields: [" + strjoin(allFields, ", ") + "]");
    end

    function refreshDomainBCVariableChoices(obj, domId)
      idx = obj.findCardIndex(obj.DomainCards, domId);
      dc = obj.DomainCards{idx};
      activeFields = obj.getDomainActiveFields(dc);
      if isempty(activeFields)
        return;
      end
      baseItems = cellstr(activeFields);
      for b = 1:numel(dc.BCCards)
        bc = dc.BCCards{b};
        if isfield(bc.Controls, 'variable')
          vw = bc.Controls.variable.widget;
          curr = string(vw.Value);
          if strlength(curr) > 0 && ~any(strcmp(baseItems, char(curr)))
            vw.Items = [{char(curr)}, baseItems];
          else
            vw.Items = baseItems;
          end
          if strlength(curr) == 0 || ~any(strcmp(vw.Items, char(curr)))
            vw.Value = vw.Items{1};
          end
        end
      end
    end

    function onDomainGeomEnableChanged(obj, domId)
      idx = obj.findCardIndex(obj.DomainCards, domId);
      dc = obj.DomainCards{idx};
      isOn = dc.GeomEnableBox.Value;
      dc.GeomModeDropdown.Enable = boolOnOff(isOn);
      if isOn
        dc.GeomRowsGrid.Visible = 'on';
        dc.GeomSection.BodyGrid.RowHeight{2} = 'fit';
      else
        dc.GeomRowsGrid.Visible = 'off';
        dc.GeomSection.BodyGrid.RowHeight{2} = 0;
      end
    end

    function onDomainGeomModeChanged(obj, domId)
      idx = obj.findCardIndex(obj.DomainCards, domId);
      dc = obj.DomainCards{idx};
      gMode = string(dc.GeomModeDropdown.Value);
      gEntry = obj.Schema.geometry([obj.Schema.geometry.name] == gMode);
      dc.GeomControls = obj.populateParamRows(dc.GeomRowsGrid, gEntry.params, @() [], dc.GeomSection.Tint);
      obj.DomainCards{idx} = dc;
    end

    function onDomainFluidEnableChanged(obj, domId)
      idx = obj.findCardIndex(obj.DomainCards, domId);
      dc = obj.DomainCards{idx};
      if dc.FluidEnableBox.Value
        dc.FluidRowsGrid.Visible = 'on';
        dc.MaterialsSection.BodyGrid.RowHeight{2} = 'fit';
      else
        dc.FluidRowsGrid.Visible = 'off';
        dc.MaterialsSection.BodyGrid.RowHeight{2} = 0;
      end
    end

    function sId = appendSolidCard(obj, domId)
      idx = obj.findCardIndex(obj.DomainCards, domId);
      dc = obj.DomainCards{idx};
      sId = obj.NextCardId;
      obj.NextCardId = obj.NextCardId + 1;
      nCards = numel(dc.SolidCards) + 1;

      [tint, accent] = getDepthColorPair(3);
      [roseTint, roseAccent] = getDepthColorPair(4);
      dc.SolidsListGrid.RowHeight = repmat({'fit'}, 1, nCards);

      outerPanel = uipanel(dc.SolidsListGrid, 'Title', '', ...
        'BorderType', 'line', 'BorderWidth', 2, ...
        'HighlightColor', accent, 'BackgroundColor', tint);
      outerPanel.Layout.Row = nCards;

      cGrid = uigridlayout(outerPanel, [6, 1]);
      cGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
      cGrid.RowSpacing = 4;
      cGrid.Padding = [8, 8, 8, 8];
      cGrid.BackgroundColor = tint;

      hdrGrid = uigridlayout(cGrid, [1, 2]);
      hdrGrid.Layout.Row = 1;
      hdrGrid.ColumnWidth = {'1x', 'fit'};
      hdrGrid.Padding = [8, 4, 6, 4];
      hdrGrid.BackgroundColor = accent;

      titleLbl = uilabel(hdrGrid, 'Text', sprintf('Solid Material #%d', nCards), ...
        'FontWeight', 'bold', 'FontColor', [1, 1, 1]);
      titleLbl.Layout.Column = 1;

      rmBtn = createDeleteButton(hdrGrid, 'Remove Solid', @(~,~) obj.removeSolidCard(domId, sId));
      rmBtn.Layout.Column = 2;

      bar = uigridlayout(cGrid, [1, 2]);
      bar.Layout.Row = 2;
      bar.ColumnWidth = {'fit', '1x'};
      bar.Padding = [0, 0, 0, 0];
      bar.BackgroundColor = tint;

      lawLbl = uilabel(bar, 'Text', 'Constitutive Law:', 'FontWeight', 'bold');
      lawLbl.Layout.Column = 1;

      lawItems = [{'None'}, cellstr([obj.Schema.materials.constitutiveLaws.name])];
      lawDrop = uidropdown(bar, 'Items', lawItems, 'Value', 'None', ...
        'ValueChangedFcn', @(~,~) obj.onSolidLawChanged(domId, sId));
      lawDrop.Layout.Column = 2;

      hdrRows = uigridlayout(cGrid, [1, 1]);
      hdrRows.Layout.Row = 3;
      hdrRows.Padding = [0, 0, 0, 0];
      hdrRows.BackgroundColor = tint;
      hdrCtrls = obj.populateParamRows(hdrRows, obj.Schema.materials.solidParams, @() [], tint);

      lawRows = uigridlayout(cGrid, [1, 1]);
      lawRows.Layout.Row = 4;
      lawRows.Padding = [0, 0, 0, 0];
      lawRows.BackgroundColor = tint;

      prBox = uicheckbox(cGrid, 'Text', 'Enable PorousRock', 'Value', false, ...
        'FontWeight', 'bold', 'ValueChangedFcn', @(~,~) obj.onSolidPorousChanged(domId, sId));
      prBox.Layout.Row = 5;

      prPanel = uipanel(cGrid, 'Title', 'PorousRock & Capillary Curves', ...
        'BorderType', 'line', 'BorderWidth', 1.5, ...
        'HighlightColor', roseAccent, 'BackgroundColor', roseTint);
      prPanel.Layout.Row = 6;

      prWrap = uigridlayout(prPanel, [3, 1]);
      prWrap.RowHeight = {'fit', 'fit', 'fit'};
      prWrap.Padding = [8, 6, 8, 6];
      prWrap.BackgroundColor = roseTint;

      prRows = uigridlayout(prWrap, [1, 1]);
      prRows.Layout.Row = 1;
      prRows.Padding = [0, 0, 0, 0];
      prRows.BackgroundColor = roseTint;
      prCtrls = obj.populateParamRows(prRows, obj.Schema.materials.porousRockParams, @() [], roseTint);

      cvBar = uigridlayout(prWrap, [1, 2]);
      cvBar.Layout.Row = 2;
      cvBar.ColumnWidth = {'fit', '1x'};
      cvBar.Padding = [0, 0, 0, 0];
      cvBar.BackgroundColor = roseTint;

      cvBox = uicheckbox(cvBar, 'Text', 'Enable Capillary Curves (VanGenuchten):', 'Value', false, ...
        'ValueChangedFcn', @(~,~) obj.onSolidCurveEnableChanged(domId, sId));
      cvBox.Layout.Column = 1;

      cvTypes = cellstr([obj.Schema.materials.vanGenuchten.type]);
      cvDrop = uidropdown(cvBar, 'Items', cvTypes, 'Value', cvTypes{end}, ...
        'ValueChangedFcn', @(~,~) obj.onSolidCurveTypeChanged(domId, sId));
      cvDrop.Layout.Column = 2;

      cvRows = uigridlayout(prWrap, [1, 1]);
      cvRows.Layout.Row = 3;
      cvRows.Padding = [12, 0, 0, 0];
      cvRows.BackgroundColor = roseTint;

      sCard = struct( ...
        'Id', sId, ...
        'Panel', outerPanel, ...
        'TitleLabel', titleLbl, ...
        'Grid', cGrid, ...
        'HeaderControls', hdrCtrls, ...
        'LawDropdown', lawDrop, ...
        'LawRowsGrid', lawRows, ...
        'LawControls', struct(), ...
        'PorousEnableBox', prBox, ...
        'PorousPanel', prPanel, ...
        'PorousWrap', prWrap, ...
        'PorousControls', prCtrls, ...
        'CurveEnableBox', cvBox, ...
        'CurveTypeDropdown', cvDrop, ...
        'CurveRowsGrid', cvRows, ...
        'CurveControls', struct(), ...
        'TintColor', tint, ...
        'RoseTintColor', roseTint);

      dc.SolidCards{end+1} = sCard;
      obj.DomainCards{idx} = dc;
      obj.onSolidLawChanged(domId, sId);
      obj.onSolidPorousChanged(domId, sId);
      obj.onSolidCurveTypeChanged(domId, sId);
      obj.onSolidCurveEnableChanged(domId, sId);
    end

    function removeSolidCard(obj, domId, sId)
      [dIdx, sIdx] = obj.findSolidCard(domId, sId);
      if isempty(sIdx)
        return;
      end
      dc = obj.DomainCards{dIdx};
      delete(dc.SolidCards{sIdx}.Panel);
      dc.SolidCards(sIdx) = [];
      relayoutAccentCards(dc.SolidsListGrid, dc.SolidCards, 'Solid Material');
      obj.DomainCards{dIdx} = dc;
    end

    function onSolidLawChanged(obj, domId, sId)
      [dIdx, sIdx] = obj.findSolidCard(domId, sId);
      dc = obj.DomainCards{dIdx};
      sCard = dc.SolidCards{sIdx};
      lawName = string(sCard.LawDropdown.Value);
      if lawName == "None"
        delete(sCard.LawRowsGrid.Children);
        sCard.LawControls = struct();
        sCard.Grid.RowHeight{4} = 0;
        sCard.LawRowsGrid.Visible = 'off';
      else
        laws = obj.Schema.materials.constitutiveLaws;
        lEntry = laws([laws.name] == lawName);
        sCard.LawControls = obj.populateParamRows(sCard.LawRowsGrid, lEntry.params, @() [], sCard.TintColor);
        sCard.LawRowsGrid.Visible = 'on';
        sCard.Grid.RowHeight{4} = 'fit';
      end
      dc.SolidCards{sIdx} = sCard;
      obj.DomainCards{dIdx} = dc;
    end

    function onSolidPorousChanged(obj, domId, sId)
      [dIdx, sIdx] = obj.findSolidCard(domId, sId);
      dc = obj.DomainCards{dIdx};
      sCard = dc.SolidCards{sIdx};
      if sCard.PorousEnableBox.Value
        sCard.PorousPanel.Visible = 'on';
        sCard.Grid.RowHeight{6} = 'fit';
        if ~dc.FluidEnableBox.Value
          dc.FluidEnableBox.Value = true;
          obj.onDomainFluidEnableChanged(domId);
        end
      else
        sCard.PorousPanel.Visible = 'off';
        sCard.Grid.RowHeight{6} = 0;
      end
    end

    function onSolidCurveEnableChanged(obj, domId, sId)
      [dIdx, sIdx] = obj.findSolidCard(domId, sId);
      sCard = obj.DomainCards{dIdx}.SolidCards{sIdx};
      isOn = sCard.CurveEnableBox.Value;
      sCard.CurveTypeDropdown.Enable = boolOnOff(isOn);
      if isOn
        sCard.CurveRowsGrid.Visible = 'on';
        sCard.PorousWrap.RowHeight{3} = 'fit';
      else
        sCard.CurveRowsGrid.Visible = 'off';
        sCard.PorousWrap.RowHeight{3} = 0;
      end
    end

    function onSolidCurveTypeChanged(obj, domId, sId)
      [dIdx, sIdx] = obj.findSolidCard(domId, sId);
      dc = obj.DomainCards{dIdx};
      sCard = dc.SolidCards{sIdx};
      cType = string(sCard.CurveTypeDropdown.Value);
      vg = obj.Schema.materials.vanGenuchten;
      vEntry = vg([vg.type] == cType);
      sCard.CurveControls = obj.populateParamRows(sCard.CurveRowsGrid, vEntry.params, @() [], sCard.RoseTintColor);
      dc.SolidCards{sIdx} = sCard;
      obj.DomainCards{dIdx} = dc;
    end

    function bcId = appendBCCard(obj, domId)
      idx = obj.findCardIndex(obj.DomainCards, domId);
      dc = obj.DomainCards{idx};
      bcId = obj.NextCardId;
      obj.NextCardId = obj.NextCardId + 1;
      nCards = numel(dc.BCCards) + 1;

      [tint, accent] = getDepthColorPair(3);
      dc.BCListGrid.RowHeight = repmat({'fit'}, 1, nCards);

      outerPanel = uipanel(dc.BCListGrid, 'Title', '', ...
        'BorderType', 'line', 'BorderWidth', 2, ...
        'HighlightColor', accent, 'BackgroundColor', tint);
      outerPanel.Layout.Row = nCards;

      bGrid = uigridlayout(outerPanel, [4, 1]);
      bGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
      bGrid.RowSpacing = 5;
      bGrid.Padding = [8, 8, 8, 8];
      bGrid.BackgroundColor = tint;

      hdrGrid = uigridlayout(bGrid, [1, 2]);
      hdrGrid.Layout.Row = 1;
      hdrGrid.ColumnWidth = {'1x', 'fit'};
      hdrGrid.Padding = [8, 4, 6, 4];
      hdrGrid.BackgroundColor = accent;

      titleLbl = uilabel(hdrGrid, 'Text', sprintf('Boundary Condition #%d', nCards), ...
        'FontWeight', 'bold', 'FontColor', [1, 1, 1]);
      titleLbl.Layout.Column = 1;

      rmBtn = createDeleteButton(hdrGrid, 'Remove BC', @(~,~) obj.removeBCCard(domId, bcId));
      rmBtn.Layout.Column = 2;

      bcRows = uigridlayout(bGrid, [1, 1]);
      bcRows.Layout.Row = 2;
      bcRows.Padding = [0, 0, 0, 0];
      bcRows.BackgroundColor = tint;

      bcParams = obj.Schema.boundaryConditions.bcParams;
      activeFields = obj.getDomainActiveFields(dc);
      varIdx = find([bcParams.name] == "variable", 1);
      if ~isempty(varIdx) && ~isempty(activeFields)
        bcParams(varIdx).type = "editableEnum";
        bcParams(varIdx).choices = activeFields;
        bcParams(varIdx).defaultValue = activeFields(1);
        bcParams(varIdx).allowedText = "active fields [" + strjoin(activeFields, ", ") + "] or string";
      end
      bcCtrls = obj.populateParamRows(bcRows, bcParams, @() obj.onBCTypeChanged(domId, bcId), tint);

      evListGrid = uigridlayout(bGrid, [1, 1]);
      evListGrid.Layout.Row = 3;
      evListGrid.RowHeight = {0};
      evListGrid.Padding = [12, 0, 0, 0];
      evListGrid.RowSpacing = 4;
      evListGrid.BackgroundColor = tint;

      addEvBtn = uibutton(bGrid, 'Text', '+ Add BCevent (time / value)', ...
        'ButtonPushedFcn', @(~,~) obj.appendBCEventRow(domId, bcId, []));
      addEvBtn.Layout.Row = 4;

      bcCard = struct( ...
        'Id', bcId, ...
        'Panel', outerPanel, ...
        'Grid', bGrid, ...
        'HeaderGrid', hdrGrid, ...
        'TitleLabel', titleLbl, ...
        'RowsGrid', bcRows, ...
        'Controls', bcCtrls, ...
        'EventsGrid', evListGrid, ...
        'AccentColor', accent, ...
        'Events', {{}});

      dc.BCCards{end+1} = bcCard;
      obj.DomainCards{idx} = dc;
      obj.onBCTypeChanged(domId, bcId);
    end

    function onBCTypeChanged(obj, domId, bcId)
      [dIdx, bIdx] = obj.findBCCard(domId, bcId);
      if isempty(bIdx)
        return;
      end
      bc = obj.DomainCards{dIdx}.BCCards{bIdx};
      if isfield(bc.Controls, 'type') && isfield(bc.Controls, 'essential')
        bcType = lower(strtrim(string(bc.Controls.type.widget.Value)));
        isCustom = ~ismember(bcType, ["dirichlet", "neumann"]);
        setParamRowVisible(bc.Controls.essential, isCustom);
      end
    end

    function removeBCCard(obj, domId, bcId)
      [dIdx, bIdx] = obj.findBCCard(domId, bcId);
      if isempty(bIdx)
        return;
      end
      dc = obj.DomainCards{dIdx};
      delete(dc.BCCards{bIdx}.Panel);
      dc.BCCards(bIdx) = [];
      relayoutAccentCards(dc.BCListGrid, dc.BCCards, 'Boundary Condition');
      obj.DomainCards{dIdx} = dc;
    end

    function appendBCEventRow(obj, domId, bcId, seedEvent)
      [dIdx, bIdx] = obj.findBCCard(domId, bcId);
      dc = obj.DomainCards{dIdx};
      bcCard = dc.BCCards{bIdx};
      evId = obj.NextCardId;
      obj.NextCardId = obj.NextCardId + 1;
      nEv = numel(bcCard.Events) + 1;

      [evBg, evAccent] = getDepthColorPair(4);
      bcCard.EventsGrid.RowHeight = repmat({'fit'}, 1, nEv);
      evPanel = uipanel(bcCard.EventsGrid, 'Title', sprintf('BCevent #%d', nEv), ...
        'BorderType', 'line', 'BorderWidth', 1.5, ...
        'HighlightColor', evAccent, 'BackgroundColor', evBg);
      evPanel.Layout.Row = nEv;

      evGrid = uigridlayout(evPanel, [2, 1]);
      evGrid.RowHeight = {'fit', 'fit'};
      evGrid.Padding = [6, 4, 6, 4];
      evGrid.BackgroundColor = evBg;

      rmBar = uigridlayout(evGrid, [1, 2]);
      rmBar.Layout.Row = 1;
      rmBar.ColumnWidth = {'1x', 'fit'};
      rmBar.Padding = [0, 0, 0, 0];
      rmBar.BackgroundColor = evBg;

      rmBtn = createDeleteButton(rmBar, 'Remove Event', @(~,~) obj.removeBCEventRow(domId, bcId, evId));
      rmBtn.Layout.Column = 2;

      rGrid = uigridlayout(evGrid, [1, 1]);
      rGrid.Layout.Row = 2;
      rGrid.Padding = [0, 0, 0, 0];
      rGrid.BackgroundColor = evBg;
      evCtrls = obj.populateParamRows(rGrid, obj.Schema.boundaryConditions.eventParams, @() [], evBg);

      if ~isempty(seedEvent) && isstruct(seedEvent)
        evCtrls = applyStructToGroup(evCtrls, seedEvent);
      end

      bcCard.Events{end+1} = struct('Id', evId, 'Panel', evPanel, 'Controls', evCtrls);
      dc.BCCards{bIdx} = bcCard;
      obj.DomainCards{dIdx} = dc;
    end

    function removeBCEventRow(obj, domId, bcId, evId)
      [dIdx, bIdx] = obj.findBCCard(domId, bcId);
      dc = obj.DomainCards{dIdx};
      bcCard = dc.BCCards{bIdx};
      evIdx = obj.findCardIndex(bcCard.Events, evId);
      if isempty(evIdx)
        return;
      end
      delete(bcCard.Events{evIdx}.Panel);
      bcCard.Events(evIdx) = [];
      relayoutPanelList(bcCard.EventsGrid, bcCard.Events, 'BCevent');
      dc.BCCards{bIdx} = bcCard;
      obj.DomainCards{dIdx} = dc;
    end

    function ifaceId = appendInterfaceCard(obj)
      ifaceId = obj.NextCardId;
      obj.NextCardId = obj.NextCardId + 1;
      nIf = numel(obj.InterfaceCards) + 1;

      [tint, accent] = getDepthColorPair(1);
      [tolTint, tolAccent] = getDepthColorPair(3);
      obj.InterfacesListGrid.RowHeight = repmat({'fit'}, 1, nIf);

      outerPanel = uipanel(obj.InterfacesListGrid, 'Title', '', ...
        'BorderType', 'line', 'BorderWidth', 2, ...
        'HighlightColor', accent, 'BackgroundColor', tint);
      outerPanel.Layout.Row = nIf;

      cGrid = uigridlayout(outerPanel, [7, 1]);
      cGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
      cGrid.RowSpacing = 5;
      cGrid.Padding = [8, 8, 8, 8];
      cGrid.BackgroundColor = tint;

      hdrBar = uigridlayout(cGrid, [1, 2]);
      hdrBar.Layout.Row = 1;
      hdrBar.ColumnWidth = {'1x', 'fit'};
      hdrBar.Padding = [8, 4, 6, 4];
      hdrBar.BackgroundColor = accent;

      titleLbl = uilabel(hdrBar, 'Text', sprintf('Interface #%d', nIf), ...
        'FontWeight', 'bold', 'FontColor', [1, 1, 1]);
      titleLbl.Layout.Column = 1;

      rmBtn = createDeleteButton(hdrBar, 'Remove Interface', @(~,~) obj.removeInterfaceById(ifaceId));
      rmBtn.Layout.Column = 2;

      mBar = uigridlayout(cGrid, [1, 2]);
      mBar.Layout.Row = 2;
      mBar.ColumnWidth = {'fit', '1x'};
      mBar.Padding = [0, 0, 0, 0];
      mBar.BackgroundColor = tint;

      mLbl = uilabel(mBar, 'Text', 'Interface Model *:', 'FontWeight', 'bold');
      mLbl.Layout.Column = 1;

      modelNames = cellstr([obj.Schema.interfaces.models.name]);
      mDrop = uidropdown(mBar, 'Items', modelNames, 'Value', modelNames{1}, ...
        'ValueChangedFcn', @(~,~) obj.onInterfaceModelChanged(ifaceId));
      mDrop.Layout.Column = 2;

      commonRows = uigridlayout(cGrid, [1, 1]);
      commonRows.Layout.Row = 3;
      commonRows.Padding = [0, 0, 0, 0];
      commonRows.BackgroundColor = tint;
      commonCtrls = obj.populateParamRows(commonRows, obj.Schema.interfaces.commonParams, ...
        @() obj.onInterfaceDomainsOrModelUpdated(ifaceId), tint);

      modelRows = uigridlayout(cGrid, [1, 1]);
      modelRows.Layout.Row = 4;
      modelRows.Padding = [0, 0, 0, 0];
      modelRows.BackgroundColor = tint;

      % Depth 2 Yellow subsections
      quadSec = obj.createCollapsibleSection(cGrid, 5, "Mortar Quadrature (optional)", true, 2);
      qGrid = quadSec.BodyGrid;
      qGrid.RowHeight = {'fit', 'fit'};
      qHdr = uigridlayout(qGrid, [1, 3]);
      qHdr.Layout.Row = 1;
      qHdr.ColumnWidth = {'fit', 'fit', '1x'};
      qHdr.Padding = [0, 0, 0, 0];
      qHdr.BackgroundColor = quadSec.Tint;

      qEnBox = uicheckbox(qHdr, 'Text', 'Customize Quadrature', 'Value', false, ...
        'FontWeight', 'bold', 'ValueChangedFcn', @(~,~) obj.onInterfaceQuadEnableChanged(ifaceId));
      qEnBox.Layout.Column = 1;
      qTLbl = uilabel(qHdr, 'Text', 'Quadrature Type:');
      qTLbl.Layout.Column = 2;
      quadNames = cellstr([obj.Schema.interfaces.quadratures.name]);
      qDrop = uidropdown(qHdr, 'Items', quadNames, 'Value', quadNames{1}, ...
        'ValueChangedFcn', @(~,~) obj.onInterfaceQuadTypeChanged(ifaceId));
      qDrop.Layout.Column = 3;

      qRows = uigridlayout(qGrid, [1, 1]);
      qRows.Layout.Row = 2;
      qRows.Padding = [0, 0, 0, 0];
      qRows.BackgroundColor = quadSec.Tint;

      coulombSec = obj.createCollapsibleSection(cGrid, 6, "Coulomb Friction *", false, 2);

      asSec = obj.createCollapsibleSection(cGrid, 7, "ActiveSet & Tolerances (optional)", true, 2);
      asGrid = asSec.BodyGrid;
      asGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
      asEnBox = uicheckbox(asGrid, 'Text', 'Customize ActiveSet', 'Value', false, ...
        'FontWeight', 'bold', 'ValueChangedFcn', @(~,~) obj.onInterfaceActiveSetChanged(ifaceId));
      asEnBox.Layout.Row = 1;
      asRows = uigridlayout(asGrid, [1, 1]);
      asRows.Layout.Row = 2;
      asRows.Padding = [0, 0, 0, 0];
      asRows.BackgroundColor = asSec.Tint;

      tolEnBox = uicheckbox(asGrid, 'Text', 'Customize Tolerances', 'Value', false, ...
        'FontWeight', 'bold', 'ValueChangedFcn', @(~,~) obj.onInterfaceTolChanged(ifaceId));
      tolEnBox.Layout.Row = 3;

      tolPanel = uipanel(asGrid, 'Title', 'Tolerances', ...
        'BorderType', 'line', 'BorderWidth', 1.5, ...
        'HighlightColor', tolAccent, 'BackgroundColor', tolTint);
      tolPanel.Layout.Row = 4;

      tolRows = uigridlayout(tolPanel, [1, 1]);
      tolRows.Padding = [8, 6, 8, 6];
      tolRows.BackgroundColor = tolTint;

      ic = struct( ...
        'Id', ifaceId, ...
        'Panel', outerPanel, ...
        'TitleLabel', titleLbl, ...
        'Grid', cGrid, ...
        'TintColor', tint, ...
        'TolTintColor', tolTint, ...
        'ModelDropdown', mDrop, ...
        'ModelEntry', obj.Schema.interfaces.models(1), ...
        'CommonControls', commonCtrls, ...
        'ModelRowsGrid', modelRows, ...
        'ModelControls', struct(), ...
        'QuadSection', quadSec, ...
        'QuadEnableBox', qEnBox, ...
        'QuadTypeDropdown', qDrop, ...
        'QuadRowsGrid', qRows, ...
        'QuadControls', struct(), ...
        'CoulombSection', coulombSec, ...
        'CoulombControls', struct(), ...
        'ActiveSetSection', asSec, ...
        'ActiveSetEnableBox', asEnBox, ...
        'ActiveSetRowsGrid', asRows, ...
        'ActiveSetControls', struct(), ...
        'TolEnableBox', tolEnBox, ...
        'TolPanel', tolPanel, ...
        'TolRowsGrid', tolRows, ...
        'TolControls', struct());

      obj.InterfaceCards{end+1} = ic;
      obj.refreshInterfaceDomainDropdowns();
      obj.onInterfaceModelChanged(ifaceId);
      obj.onInterfaceQuadTypeChanged(ifaceId);
      obj.onInterfaceQuadEnableChanged(ifaceId);
      obj.syncConsensusScheme();
    end

    function removeInterfaceById(obj, ifaceId)
      idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
      if isempty(idx)
        return;
      end
      delete(obj.InterfaceCards{idx}.Panel);
      obj.InterfaceCards(idx) = [];
      relayoutAccentCards(obj.InterfacesListGrid, obj.InterfaceCards, 'Interface');
      obj.syncConsensusScheme();
    end

    function clearAllInterfaces(obj)
      for i = 1:numel(obj.InterfaceCards)
        delete(obj.InterfaceCards{i}.Panel);
      end
      obj.InterfaceCards = {};
      obj.InterfacesListGrid.RowHeight = {0};
    end

    function onInterfaceModelChanged(obj, ifaceId)
      idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
      ic = obj.InterfaceCards{idx};
      mName = string(ic.ModelDropdown.Value);
      models = obj.Schema.interfaces.models;
      mEntry = models([models.name] == mName);
      ic.ModelEntry = mEntry;

      ic.ModelControls = obj.populateParamRows(ic.ModelRowsGrid, mEntry.params, @() [], ic.TintColor);

      if isfield(ic.CommonControls, 'multiplierType')
        mtCtrl = ic.CommonControls.multiplierType;
        setParamRowVisible(mtCtrl, mEntry.supportsMultiplierType);
        if mEntry.supportsMultiplierType
          items = cellstr(mEntry.multiplierChoices);
          mtCtrl.widget.Items = items;
          if ~any(strcmp(items, mtCtrl.widget.Value))
            mtCtrl.widget.Value = items{1};
          end
        end
      end

      setGridRowVisible(ic.Grid, ic.QuadSection, 5, mEntry.supportsQuadrature);

      hasCoulomb = ~isempty(mEntry.coulombParams);
      setGridRowVisible(ic.Grid, ic.CoulombSection, 6, hasCoulomb);
      if hasCoulomb
        ic.CoulombControls = obj.populateParamRows(ic.CoulombSection.BodyGrid, ...
          mEntry.coulombParams, @() [], ic.CoulombSection.Tint);
      else
        delete(ic.CoulombSection.BodyGrid.Children);
        ic.CoulombControls = struct();
      end

      hasAS = ~isempty(mEntry.activeSetParams);
      setGridRowVisible(ic.Grid, ic.ActiveSetSection, 7, hasAS);
      if hasAS
        ic.ActiveSetControls = obj.populateParamRows(ic.ActiveSetRowsGrid, ...
          mEntry.activeSetParams, @() [], ic.ActiveSetSection.Tint);
        ic.TolControls = obj.populateParamRows(ic.TolRowsGrid, ...
          mEntry.activeSetTolParams, @() [], ic.TolTintColor);
      else
        delete(ic.ActiveSetRowsGrid.Children);
        delete(ic.TolRowsGrid.Children);
        ic.ActiveSetControls = struct();
        ic.TolControls = struct();
      end

      obj.InterfaceCards{idx} = ic;
      if hasAS
        obj.onInterfaceActiveSetChanged(ifaceId);
      end
      obj.onInterfaceDomainsOrModelUpdated(ifaceId);
    end

    function onInterfaceDomainsOrModelUpdated(obj, ifaceId)
      idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
      if isempty(idx)
        return;
      end
      ic = obj.InterfaceCards{idx};
      if isfield(ic.CommonControls, 'variable')
        nDom = numel(obj.DomainCards);
        mDom = str2double(string(ic.CommonControls.masterDomain.widget.Value));
        sDom = str2double(string(ic.CommonControls.slaveDomain.widget.Value));
        choices = string.empty;
        if mDom >= 1 && mDom <= nDom && sDom >= 1 && sDom <= nDom
          fM = obj.getDomainActiveFields(obj.DomainCards{mDom});
          fS = obj.getDomainActiveFields(obj.DomainCards{sDom});
          choices = intersect(fM, fS, 'stable');
        end
        if isempty(choices)
          choices = ["displacements", "pressure"];
        end
        varW = ic.CommonControls.variable.widget;
        curr = string(varW.Value);
        items = [{''}, cellstr(choices)];
        if strlength(curr) > 0 && ~any(strcmp(items, char(curr)))
          items = [{char(curr)}, items];
        end
        varW.Items = items;
      end
      obj.updateCompatibilityBanner();
    end

    function onInterfaceQuadEnableChanged(obj, ifaceId)
      idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
      ic = obj.InterfaceCards{idx};
      isOn = ic.QuadEnableBox.Value;
      ic.QuadTypeDropdown.Enable = boolOnOff(isOn);
      if isOn
        ic.QuadRowsGrid.Visible = 'on';
        ic.QuadSection.BodyGrid.RowHeight{2} = 'fit';
      else
        ic.QuadRowsGrid.Visible = 'off';
        ic.QuadSection.BodyGrid.RowHeight{2} = 0;
      end
    end

    function onInterfaceQuadTypeChanged(obj, ifaceId)
      idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
      ic = obj.InterfaceCards{idx};
      qType = string(ic.QuadTypeDropdown.Value);
      quads = obj.Schema.interfaces.quadratures;
      qEntry = quads([quads.name] == qType);
      ic.QuadControls = obj.populateParamRows(ic.QuadRowsGrid, qEntry.params, @() [], ic.QuadSection.Tint);
      obj.InterfaceCards{idx} = ic;
    end

    function onInterfaceActiveSetChanged(obj, ifaceId)
      idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
      ic = obj.InterfaceCards{idx};
      isOn = ic.ActiveSetEnableBox.Value;
      asGrid = ic.ActiveSetSection.BodyGrid;
      if isOn
        ic.ActiveSetRowsGrid.Visible = 'on';
        ic.TolEnableBox.Visible = 'on';
        asGrid.RowHeight{2} = 'fit';
        asGrid.RowHeight{3} = 'fit';
      else
        ic.ActiveSetRowsGrid.Visible = 'off';
        ic.TolEnableBox.Visible = 'off';
        asGrid.RowHeight{2} = 0;
        asGrid.RowHeight{3} = 0;
      end
      obj.onInterfaceTolChanged(ifaceId);
    end

    function onInterfaceTolChanged(obj, ifaceId)
      idx = obj.findCardIndex(obj.InterfaceCards, ifaceId);
      ic = obj.InterfaceCards{idx};
      showTol = ic.ActiveSetEnableBox.Value && ic.TolEnableBox.Value;
      if showTol
        ic.TolPanel.Visible = 'on';
        ic.ActiveSetSection.BodyGrid.RowHeight{4} = 'fit';
      else
        ic.TolPanel.Visible = 'off';
        ic.ActiveSetSection.BodyGrid.RowHeight{4} = 0;
      end
    end

    function refreshInterfaceDomainDropdowns(obj)
      nDom = max(1, numel(obj.DomainCards));
      domItems = arrayfun(@(k) sprintf('%d', k), 1:nDom, 'UniformOutput', false);
      for i = 1:numel(obj.InterfaceCards)
        ic = obj.InterfaceCards{i};
        for key = ["masterDomain", "slaveDomain"]
          if isfield(ic.CommonControls, key)
            dd = ic.CommonControls.(key).widget;
            prev = dd.Value;
            dd.Items = domItems;
            if any(strcmp(domItems, prev))
              dd.Value = prev;
            elseif strcmp(key, 'masterDomain') && nDom >= 2
              dd.Value = '2';
            else
              dd.Value = '1';
            end
          end
        end
        obj.onInterfaceDomainsOrModelUpdated(ic.Id);
      end
    end

    function syncConsensusScheme(obj)
      consensus = obj.computeConsensusScheme();
      if ~obj.UserOverrodeScheme
        obj.SchemeDropdown.Value = char(consensus);
      end
      obj.rebuildSchemeOptionsPanel();
      if string(obj.SchemeDropdown.Value) == consensus
        obj.UserOverrodeScheme = false;
        obj.SchemeConsensusLabel.Text = 'Auto-selected from domain physics solvers';
      else
        obj.SchemeConsensusLabel.Text = char("Manual override (suggested: " + consensus + ")");
      end
      obj.updateCompatibilityBanner();
    end

    function onUserSchemeChanged(obj)
      consensus = obj.computeConsensusScheme();
      chosen = string(obj.SchemeDropdown.Value);
      obj.UserOverrodeScheme = (chosen ~= consensus);
      obj.rebuildSchemeOptionsPanel();
      if ~obj.UserOverrodeScheme
        obj.SchemeConsensusLabel.Text = 'Auto-selected from domain physics solvers';
      else
        obj.SchemeConsensusLabel.Text = char("Manual override (suggested: " + consensus + ")");
      end
      obj.updateCompatibilityBanner();
    end

    function rebuildSchemeOptionsPanel(obj)
      obj.CurrentSchemeName = string(obj.SchemeDropdown.Value);
      schEntry = obj.Schema.solutionSchemes([obj.Schema.solutionSchemes.name] == obj.CurrentSchemeName);

      pGrid = obj.ProblemSolverSection.BodyGrid;
      bodyGrid = obj.SchemeSubSection.BodyGrid;
      delete(bodyGrid.Children);
      obj.SchemeControls = struct();

      if isempty(schEntry.params)
        setGridRowVisible(pGrid, obj.SchemeSubSection, 3, false);
      else
        setGridRowVisible(pGrid, obj.SchemeSubSection, 3, true);
        schHasReq = any([schEntry.params.required]);
        schTitle = "Solution Scheme Options (" + schEntry.name + ")";
        if schHasReq
          schTitle = schTitle + " *";
        end
        obj.SchemeSubSection.BaseTitle = schTitle;
        obj.SchemeControls = obj.populateParamRows(bodyGrid, schEntry.params, @() [], obj.SchemeSubSection.Tint);
        obj.setSectionCollapsed(obj.SchemeSubSection, ~schHasReq);
      end
    end

    function onApplySuggestedClicked(obj)
      if strlength(obj.SuggestedScheme) > 0
        obj.UserOverrodeScheme = false;
        obj.selectSolutionScheme(obj.SuggestedScheme);
      end
    end

    function updateCompatibilityBanner(obj)
      [severity, msg, suggested] = obj.checkSchemeCompatibility();
      obj.SuggestedScheme = suggested;
      pGrid = obj.ProblemSolverSection.BodyGrid;

      if severity == "ok"
        obj.WarningPanel.Visible = 'off';
        pGrid.RowHeight{2} = 0;
        return;
      end

      if severity == "error"
        bg = [1.0, 0.90, 0.90];
        bd = [0.80, 0.18, 0.18];
        fg = [0.55, 0.05, 0.05];
      else
        bg = [1.0, 0.95, 0.84];
        bd = [0.85, 0.55, 0.10];
        fg = [0.55, 0.25, 0.0];
      end

      obj.WarningPanel.BackgroundColor = bg;
      obj.WarningPanel.HighlightColor = bd;
      wGrid = obj.WarningPanel.Children(1);
      wGrid.BackgroundColor = bg;
      obj.WarningLabel.Text = char(msg);
      obj.WarningLabel.FontColor = fg;

      if strlength(suggested) > 0
        obj.ApplySuggestedBtn.Visible = 'on';
        obj.ApplySuggestedBtn.Text = char("Use " + suggested);
        wGrid.ColumnWidth = {'1x', 'fit'};
      else
        obj.ApplySuggestedBtn.Visible = 'off';
        wGrid.ColumnWidth = {'1x', 0};
      end

      obj.WarningPanel.Visible = 'on';
      pGrid.RowHeight{2} = 'fit';
    end

    function onOutputEnableChanged(obj)
      isOn = obj.OutputEnableBox.Value;
      oGrid = obj.OutputSection.BodyGrid;
      if isOn
        obj.OutputRowsGrid.Visible = 'on';
        obj.RestartSection.Container.Visible = 'on';
        oGrid.RowHeight{2} = 'fit';
        oGrid.RowHeight{3} = 'fit';
      else
        obj.OutputRowsGrid.Visible = 'off';
        obj.RestartSection.Container.Visible = 'off';
        oGrid.RowHeight{2} = 0;
        oGrid.RowHeight{3} = 0;
      end
    end

    function onRestartEnableChanged(obj)
      isOn = obj.RestartEnableBox.Value;
      rGrid = obj.RestartSection.BodyGrid;
      if isOn
        obj.RestartRowsGrid.Visible = 'on';
        rGrid.RowHeight{2} = 'fit';
      else
        obj.RestartRowsGrid.Visible = 'off';
        rGrid.RowHeight{2} = 0;
      end
    end

    function onLinearOptionChanged(obj)
      useMatlab = isfield(obj.LinearTopControls, 'useMatlab') && ...
        logical(obj.LinearTopControls.useMatlab.widget.Value);

      fields = obj.getAllActiveFields();
      hasFlow = ~useMatlab && any(ismember(fields, ["pressure", "u"]));
      hasMech = ~useMatlab && any(fields == "displacements");
      setGridRowVisible(obj.LinearSection.BodyGrid, obj.LinearFlowSection, 2, hasFlow);
      setGridRowVisible(obj.LinearSection.BodyGrid, obj.LinearMechSection, 3, hasMech);

      if isfield(obj.LinearTopControls, 'saveMatType') && isfield(obj.LinearTopControls, 'saveMatSolves')
        smt = obj.LinearTopControls.saveMatType;
        sms = obj.LinearTopControls.saveMatSolves;
        isSel = strcmp(string(smt.widget.Value), "selected");
        setParamRowVisible(sms, isSel);
      end

      if isfield(obj.LinearTopControls, 'useEW')
        ewActive = logical(obj.LinearTopControls.useEW.widget.Value);
        if isfield(obj.LinearTopControls, 'eta0')
          setParamRowVisible(obj.LinearTopControls.eta0, ewActive);
        end
        if isfield(obj.LinearTopControls, 'tol')
          setParamRowVisible(obj.LinearTopControls.tol, ~ewActive);
        end
      end
    end

    function updateGlobalLinearSolverVisibility(obj)
      fields = obj.getAllActiveFields();
      hasSat = any(fields == "saturation");
      if isfield(obj.LinearTopControls, 'cprType')
        setParamRowVisible(obj.LinearTopControls.cprType, hasSat);
      end
      obj.onLinearOptionChanged();
    end

    function onExportClicked(obj)
      [isValid, ~] = obj.validateInputs('RequireGeometry', false);
      if ~isValid
        return;
      end

      mode = string(obj.ExportModeDropdown.Value);
      if mode == "unified"
        [file, path] = uiputfile('*.xml', 'Export GReS Unified Problem XML', 'problem.xml');
        if isequal(file, 0)
          return;
        end
        target = fullfile(path, file);
      else
        selDir = uigetdir(pwd, 'Select Folder for Modular XMLs + runSimulation.m');
        if isequal(selDir, 0)
          return;
        end
        target = selDir;
      end

      obj.exportToPath(target, mode);
    end

    function onRunClicked(obj)
      [isValid, ~] = obj.validateInputs('RequireGeometry', true);
      if ~isValid
        return;
      end

      [compSev, compMsg, ~] = obj.checkSchemeCompatibility();
      if compSev ~= "ok" && strcmp(obj.Figure.Visible, 'on')
        choice = uiconfirm(obj.Figure, char(compMsg + " Continue anyway?"), ...
          'Problem Solver Compatibility Warning', ...
          'Options', {'Run Anyway', 'Cancel'}, 'DefaultOption', 2, 'Icon', 'warning');
        if strcmp(choice, 'Cancel')
          return;
        end
      end

      selDir = uigetdir(pwd, 'Select Working Folder to Save & Run Simulation');
      if isequal(selDir, 0)
        return;
      end
      obj.runSimulation(selDir);
    end

    function onLoadClicked(obj)
      [file, path] = uigetfile('*.xml', 'Load GReS XML File(s)', 'MultiSelect', 'on');
      if isequal(file, 0)
        return;
      end
      paths = string(fullfile(path, file));
      obj.loadXmlFiles(paths);
    end

    function resetDomainsToCount(obj, nTarget)
      while numel(obj.DomainCards) > 1
        obj.removeDomainById(obj.DomainCards{end}.Id);
      end
      obj.resetDomainContent(1);
      for d = 2:nTarget
        obj.appendDomainCard();
      end
    end

    function resetDomainContent(obj, domIdx)
      dc = obj.DomainCards{domIdx};
      dc.GeomEnableBox.Value = false;
      obj.onDomainGeomEnableChanged(dc.Id);
      obj.DomainCards{domIdx} = dc;
      obj.clearDomainMaterials(domIdx);
      obj.clearDomainBoundaries(domIdx);
    end

    function clearDomainMaterials(obj, domIdx)
      dc = obj.DomainCards{domIdx};
      dc.FluidEnableBox.Value = false;
      obj.onDomainFluidEnableChanged(dc.Id);
      for i = 1:numel(dc.SolidCards)
        delete(dc.SolidCards{i}.Panel);
      end
      dc.SolidCards = {};
      dc.SolidsListGrid.RowHeight = {0};
      obj.DomainCards{domIdx} = dc;
    end

    function clearDomainBoundaries(obj, domIdx)
      dc = obj.DomainCards{domIdx};
      for i = 1:numel(dc.BCCards)
        delete(dc.BCCards{i}.Panel);
      end
      dc.BCCards = {};
      dc.BCListGrid.RowHeight = {0};
      obj.DomainCards{domIdx} = dc;
    end

    function visited = applyDomainStruct(obj, domNode, domDir, domIdx, visited)
      if nargin < 5
        visited = string.empty;
      end
      if isfield(domNode, 'Solver')
        [solNode, ~, visited] = resolveXmlNodeWithDir(domNode.Solver, domDir, visited);
        obj.applySolverStruct(solNode, domIdx);
      end
      if isfield(domNode, 'Geometry')
        [gNode, ~, visited] = resolveXmlNodeWithDir(domNode.Geometry, domDir, visited);
        obj.applyGeometryStruct(gNode, domIdx);
      end
      if isfield(domNode, 'Materials')
        [mNode, ~, visited] = resolveXmlNodeWithDir(domNode.Materials, domDir, visited);
        obj.applyMaterialsStruct(mNode, domIdx);
      end
      if isfield(domNode, 'BoundaryConditions')
        [bNode, ~, visited] = resolveXmlNodeWithDir(domNode.BoundaryConditions, domDir, visited);
        obj.applyBoundariesStruct(bNode, domIdx);
      end
    end

    function applySimParamsStruct(obj, spStruct)
      obj.SimControls = applyStructToGroup(obj.SimControls, spStruct);
      if isfield(spStruct, 'LinearSolver') && isstruct(spStruct.LinearSolver)
        obj.applyLinearSolverStruct(spStruct.LinearSolver);
      end
    end

    function applyLinearSolverStruct(obj, ls)
      if ~isstruct(ls)
        return;
      end
      obj.LinearTopControls = applyStructToGroup(obj.LinearTopControls, ls);
      if isfield(ls, 'Flow') && isstruct(ls.Flow)
        obj.ChronosFlowControls = applyChronosStruct(obj.ChronosFlowControls, ls.Flow);
      end
      if isfield(ls, 'Mechanics') && isstruct(ls.Mechanics)
        obj.ChronosMechControls = applyChronosStruct(obj.ChronosMechControls, ls.Mechanics);
      end
      obj.onLinearOptionChanged();
    end

    function applyGeometryStruct(obj, gStruct, domIdx)
      if isfield(gStruct, 'meshFile')
        obj.setGeometry('meshFile', struct('meshFile', gStruct.meshFile), domIdx);
      elseif isfield(gStruct, 'StructuredMesh') && isstruct(gStruct.StructuredMesh)
        obj.setGeometry('StructuredMesh', gStruct.StructuredMesh, domIdx);
      elseif isfield(gStruct, 'BlockStructuredMesh') && isstruct(gStruct.BlockStructuredMesh)
        obj.setGeometry('BlockStructuredMesh', gStruct.BlockStructuredMesh, domIdx);
      end
    end

    function applyMaterialsStruct(obj, matStruct, domIdx)
      if isfield(matStruct, 'Fluid') && isstruct(matStruct.Fluid)
        obj.setFluid(matStruct.Fluid, domIdx);
      end
      if isfield(matStruct, 'Solid') && isstruct(matStruct.Solid)
        for i = 1:numel(matStruct.Solid)
          obj.addSolid(matStruct.Solid(i), domIdx);
        end
      end
    end

    function applyBoundariesStruct(obj, bcStruct, domIdx)
      if isfield(bcStruct, 'BC') && isstruct(bcStruct.BC)
        for i = 1:numel(bcStruct.BC)
          obj.addBoundaryCondition(bcStruct.BC(i), domIdx);
        end
      end
    end

    function applySolverStruct(obj, solStruct, domIdx)
      known = [obj.Schema.physicsSolvers.name];
      fn = string(fieldnames(solStruct));
      match = fn(ismember(fn, known));
      if isempty(match)
        return;
      end
      sName = match(1);
      obj.selectPhysicsSolver(sName, domIdx);
      sBody = solStruct.(char(sName));
      if ~isstruct(sBody)
        return;
      end
      dc = obj.DomainCards{domIdx};
      dc.SolverControls = applyStructToGroup(dc.SolverControls, sBody);
      obj.DomainCards{domIdx} = dc;
      if obj.isDomainCoupled(dc)
        obj.onDomainCoupledSubSolverChanged(dc.Id);
        dc = obj.DomainCards{domIdx};
        if isfield(sBody, 'Mechanics') && isstruct(sBody.Mechanics)
          dc.CoupledMechControls = applyStructToGroup(dc.CoupledMechControls, sBody.Mechanics);
        end
        if isfield(sBody, 'Flow') && isstruct(sBody.Flow)
          dc.CoupledFlowControls = applyStructToGroup(dc.CoupledFlowControls, sBody.Flow);
        end
        obj.DomainCards{domIdx} = dc;
      end
    end

    function applyInterfacesStruct(obj, interfStruct)
      known = [obj.Schema.interfaces.models.name];
      fn = string(fieldnames(interfStruct));
      for k = 1:numel(fn)
        fName = fn(k);
        if ~any(known == fName)
          continue;
        end
        items = interfStruct.(char(fName));
        if ~isstruct(items)
          continue;
        end
        for m = 1:numel(items)
          obj.addInterface(fName, items(m));
        end
      end
    end

    function fields = getDomainActiveFields(obj, dc)
      sEntry = obj.getSolverEntry(dc.CurrentSolverName);
      if obj.isDomainCoupled(dc) && strlength(dc.CurrentMechSolver) > 0 && strlength(dc.CurrentFlowSolver) > 0
        mEntry = obj.getSolverEntry(dc.CurrentMechSolver);
        fEntry = obj.getSolverEntry(dc.CurrentFlowSolver);
        fields = unique([mEntry.fields, fEntry.fields]);
      else
        fields = sEntry.fields;
      end
    end

    function allFields = getAllActiveFields(obj)
      allFields = string.empty;
      for d = 1:numel(obj.DomainCards)
        allFields = [allFields, obj.getDomainActiveFields(obj.DomainCards{d})]; %#ok<AGROW>
      end
      allFields = unique(allFields);
    end

    function sec = createCollapsibleSection(obj, parentGrid, rowIdx, titleText, startCollapsed, depth)
      if nargin < 6
        depth = 0;
      end
      [tint, accent] = getDepthColorPair(depth);

      container = uigridlayout(parentGrid, [2, 1]);
      container.Layout.Row = rowIdx;
      container.RowSpacing = 2;
      container.Padding = [0, 0, 0, 0];

      headerBtn = uibutton(container, 'HorizontalAlignment', 'left', 'FontWeight', 'bold', ...
        'BackgroundColor', accent, 'FontColor', [1, 1, 1]);
      headerBtn.Layout.Row = 1;

      bodyPanel = uipanel(container, 'BorderType', 'line', 'BorderWidth', 1.5, ...
        'HighlightColor', accent, 'BackgroundColor', tint);
      bodyPanel.Layout.Row = 2;

      bodyGrid = uigridlayout(bodyPanel, [1, 1]);
      bodyGrid.RowHeight = {'fit'};
      bodyGrid.ColumnWidth = {'1x'};
      bodyGrid.Padding = [10, 8, 10, 8];
      bodyGrid.RowSpacing = 6;
      bodyGrid.BackgroundColor = tint;

      sec = struct( ...
        'Container', container, ...
        'HeaderButton', headerBtn, ...
        'BodyPanel', bodyPanel, ...
        'BodyGrid', bodyGrid, ...
        'BaseTitle', string(titleText), ...
        'Collapsed', startCollapsed, ...
        'Depth', depth, ...
        'Tint', tint, ...
        'Accent', accent);

      headerBtn.ButtonPushedFcn = @(~,~) obj.toggleSection(sec);
      obj.applySectionState(sec);
    end

    function toggleSection(obj, sec)
      currTxt = string(sec.HeaderButton.Text);
      currentlyCollapsed = startsWith(currTxt, "▶");
      sec.BaseTitle = regexprep(currTxt, '^[▶▼]\s*', '');
      obj.setSectionCollapsed(sec, ~currentlyCollapsed);
    end

    function setSectionCollapsed(obj, sec, collapsed)
      sec.Collapsed = collapsed;
      obj.applySectionState(sec);
    end

    function applySectionState(~, sec)
      if sec.Collapsed
        sec.HeaderButton.Text = char("▶ " + sec.BaseTitle);
        sec.BodyPanel.Visible = 'off';
        sec.Container.RowHeight = {'fit', 0};
      else
        sec.HeaderButton.Text = char("▼ " + sec.BaseTitle);
        sec.BodyPanel.Visible = 'on';
        sec.Container.RowHeight = {'fit', 'fit'};
      end
    end

    function ctrlMap = populateChronosSection(obj, parentGrid, chronosSections, subDepth)
      ctrlMap = struct();
      nSec = numel(chronosSections);
      if nSec == 0
        return;
      end
      parentGrid.RowHeight = repmat({'fit'}, 1, nSec);
      for i = 1:nSec
        secName = char(chronosSections(i).section);
        sub = obj.createCollapsibleSection(parentGrid, i, string(secName), true, subDepth);
        ctrlMap.(secName) = obj.populateParamRows(sub.BodyGrid, chronosSections(i).params, @() [], sub.Tint);
      end
    end

    function ctrlMap = populateParamRows(obj, parentGrid, paramList, onChangeCallback, bgColor)
      delete(parentGrid.Children);
      ctrlMap = struct();

      nP = numel(paramList);
      if nP == 0
        parentGrid.RowHeight = {'fit'};
        lbl = uilabel(parentGrid, 'Text', 'No extra parameters for this selection.', ...
          'FontAngle', 'italic', 'FontColor', [0.45 0.45 0.45]);
        lbl.Layout.Row = 1;
        return;
      end

      parentGrid.RowHeight = repmat({'fit'}, 1, nP);
      parentGrid.ColumnWidth = {'1x'};
      parentGrid.RowSpacing = 4;

      for i = 1:nP
        p = paramList(i);
        hasHelp = strlength(p.helpText) > 0;
        if hasHelp
          rowGrid = uigridlayout(parentGrid, [2, 4]);
          rowGrid.ColumnWidth = {18, '1x', 24, 240};
          widgetCol = 4;
        else
          rowGrid = uigridlayout(parentGrid, [2, 3]);
          rowGrid.ColumnWidth = {18, '1x', 240};
          widgetCol = 3;
        end
        rowGrid.Layout.Row = i;
        rowGrid.RowHeight = {'fit', 'fit'};
        rowGrid.Padding = [2, 2, 2, 4];
        rowGrid.RowSpacing = 2;
        if nargin >= 5 && ~isempty(bgColor)
          rowGrid.BackgroundColor = bgColor;
        end

        if p.required
          starLbl = uilabel(rowGrid, 'Text', '*', 'FontWeight', 'bold', ...
            'FontSize', 14, 'FontColor', [0.85 0.1 0.1], 'HorizontalAlignment', 'center');
          starLbl.Layout.Row = 1;
          starLbl.Layout.Column = 1;
        end

        nameLbl = uilabel(rowGrid, 'Text', char(p.name), 'FontWeight', 'bold');
        nameLbl.Layout.Row = 1;
        nameLbl.Layout.Column = 2;

        if hasHelp
          helpMsg = char(p.helpText);
          helpTitle = sprintf('%s — Suggestions & Guide', p.name);
          helpBtn = uibutton(rowGrid, 'Text', '?', 'Tag', 'ParamHelpButton', ...
            'FontWeight', 'bold', 'BackgroundColor', [0.22, 0.46, 0.76], ...
            'FontColor', [1, 1, 1], 'Tooltip', helpMsg, ...
            'ButtonPushedFcn', @(~,~) uialert(obj.Figure, helpMsg, helpTitle, 'Icon', 'info'));
          helpBtn.Layout.Row = 1;
          helpBtn.Layout.Column = 3;
        end

        infoText = sprintf('%s  |  Allowed: %s (Default: %s)', ...
          p.description, p.allowedText, p.defaultText);
        infoLbl = uilabel(rowGrid, 'Text', infoText, 'FontSize', 11, ...
          'FontColor', [0.35 0.35 0.35], 'WordWrap', 'on');
        infoLbl.Layout.Row = 2;
        infoLbl.Layout.Column = [2, widgetCol];

        widget = createParamWidget(rowGrid, p, onChangeCallback);
        widget.Layout.Row = 1;
        widget.Layout.Column = widgetCol;

        initVal = readWidgetValue(widget, p);
        ctrlMap.(char(p.name)) = struct( ...
          'param', p, ...
          'widget', widget, ...
          'nameLabel', nameLbl, ...
          'rowGrid', rowGrid, ...
          'rowIdx', i, ...
          'parentGrid', parentGrid, ...
          'initialValue', {initVal}, ...
          'explicitlyLoaded', false);
      end
    end

    function entry = getSolverEntry(obj, solverName)
      idx = find([obj.Schema.physicsSolvers.name] == string(solverName), 1);
      if isempty(idx)
        error('gresUI:UnknownSolver', 'Solver "%s" not found in schema.', solverName);
      end
      entry = obj.Schema.physicsSolvers(idx);
    end

    function tf = isDomainCoupled(~, dc)
      tf = ~isempty(dc.CoupledMechSection) && ~isempty(dc.CoupledFlowSection);
    end

    function tf = isLinearPhysicsVisible(~, sec)
      tf = strcmp(sec.Container.Visible, 'on');
    end

    function setStatus(obj, msg, colorRGB)
      obj.StatusLabel.Text = char(msg);
      obj.StatusLabel.FontColor = colorRGB;
    end

    function idx = findCardIndex(~, cardsCell, cardId)
      idx = [];
      for i = 1:numel(cardsCell)
        if cardsCell{i}.Id == cardId
          idx = i;
          return;
        end
      end
    end

    function [dIdx, sIdx] = findSolidCard(obj, domId, sId)
      dIdx = obj.findCardIndex(obj.DomainCards, domId);
      sIdx = obj.findCardIndex(obj.DomainCards{dIdx}.SolidCards, sId);
    end

    function [dIdx, bIdx] = findBCCard(obj, domId, bcId)
      dIdx = obj.findCardIndex(obj.DomainCards, domId);
      bIdx = obj.findCardIndex(obj.DomainCards{dIdx}.BCCards, bcId);
    end

  end

end


function btn = createDeleteButton(parent, labelText, callback)
btn = uibutton(parent, 'Text', char(labelText), ...
  'BackgroundColor', [0.80, 0.18, 0.18], ...
  'FontColor', [1, 1, 1], ...
  'FontWeight', 'bold', ...
  'ButtonPushedFcn', callback);
end


function [tint, accent] = getDepthColorPair(depth)
% Structured nesting palette: Blue (outer) -> Green -> Yellow -> Lavender -> Rose.
tints = [ ...
  0.91, 0.95, 1.00; ... % Depth 0: Blue
  0.91, 0.97, 0.92; ... % Depth 1: Green
  1.00, 0.97, 0.88; ... % Depth 2: Yellow
  0.95, 0.92, 0.99; ... % Depth 3: Lavender
  1.00, 0.93, 0.93];    % Depth 4: Rose
accents = [ ...
  0.22, 0.46, 0.76; ...
  0.20, 0.56, 0.34; ...
  0.78, 0.54, 0.12; ...
  0.48, 0.32, 0.72; ...
  0.76, 0.30, 0.38];
cIdx = mod(depth, size(tints, 1)) + 1;
tint = tints(cIdx, :);
accent = accents(cIdx, :);
end


function relayoutAccentCards(listGrid, cardsCell, titlePrefix)
nC = numel(cardsCell);
if nC == 0
  listGrid.RowHeight = {0};
  return;
end
listGrid.RowHeight = repmat({'fit'}, 1, nC);
for i = 1:nC
  cardsCell{i}.Panel.Layout.Row = i;
  cardsCell{i}.TitleLabel.Text = sprintf('%s #%d', titlePrefix, i);
end
end


function relayoutPanelList(listGrid, cardsCell, titlePrefix)
nC = numel(cardsCell);
if nC == 0
  listGrid.RowHeight = {0};
  return;
end
listGrid.RowHeight = repmat({'fit'}, 1, nC);
for i = 1:nC
  cardsCell{i}.Panel.Layout.Row = i;
  cardsCell{i}.Panel.Title = sprintf('%s #%d', titlePrefix, i);
end
end


function widget = createParamWidget(parentGrid, p, extraCallback)
switch p.type
  case "domainRef"
    defVal = char(string(p.defaultValue));
    widget = uidropdown(parentGrid, 'Items', {'1'}, 'Value', '1', ...
      'ValueChangedFcn', @(~,~) extraCallback());
    if strcmp(defVal, '2')
      widget.Items = {'1', '2'};
      widget.Value = '2';
    end

  case "enum"
    items = cellstr(p.choices);
    defVal = char(p.defaultValue);
    if ~any(strcmp(items, defVal))
      defVal = items{1};
    end
    widget = uidropdown(parentGrid, 'Items', items, 'Value', defVal, ...
      'ValueChangedFcn', @(~,~) extraCallback());

  case "editableEnum"
    items = cellstr(p.choices);
    defVal = char(p.defaultValue);
    if isempty(defVal)
      items = [{''}, items];
      defVal = '';
    elseif ~any(strcmp(items, defVal))
      items = [{defVal}, items];
    end
    widget = uidropdown(parentGrid, 'Items', items, 'Value', defVal, ...
      'Editable', 'on', 'ValueChangedFcn', @(~,~) extraCallback());

  case "flag01"
    defBool = false;
    if ~isempty(p.defaultValue)
      defBool = logical(p.defaultValue);
    end
    widget = uicheckbox(parentGrid, 'Text', 'Enabled (1)', 'Value', defBool, ...
      'ValueChangedFcn', @(~,~) extraCallback());

  case "double"
    if p.required
      widget = uieditfield(parentGrid, 'numeric', 'AllowEmpty', 'on', 'Value', [], ...
        'Placeholder', 'Required', 'ValueChangedFcn', @(~,~) extraCallback());
    elseif isempty(p.defaultValue)
      widget = uieditfield(parentGrid, 'numeric', 'AllowEmpty', 'on', 'Value', [], ...
        'Placeholder', 'Optional', 'ValueChangedFcn', @(~,~) extraCallback());
    else
      widget = uieditfield(parentGrid, 'numeric', 'AllowEmpty', 'on', ...
        'Value', double(p.defaultValue), 'ValueChangedFcn', @(~,~) extraCallback());
    end

  case "vector"
    defTxt = formatValueForWidgetText(p.defaultValue);
    widget = uieditfield(parentGrid, 'text', 'Value', defTxt, ...
      'ValueChangedFcn', @(~,~) extraCallback());

  case "file"
    subGrid = uigridlayout(parentGrid, [1, 2]);
    subGrid.ColumnWidth = {'1x', 'fit'};
    subGrid.Padding = [0, 0, 0, 0];
    subGrid.ColumnSpacing = 4;
    txtField = uieditfield(subGrid, 'text', 'Value', '', ...
      'ValueChangedFcn', @(~,~) extraCallback());
    txtField.Layout.Column = 1;
    brBtn = uibutton(subGrid, 'Text', 'Browse...', ...
      'ButtonPushedFcn', @(~,~) browseIntoField(txtField, extraCallback));
    brBtn.Layout.Column = 2;
    widget = subGrid;

  otherwise
    defTxt = formatValueForWidgetText(p.defaultValue);
    widget = uieditfield(parentGrid, 'text', 'Value', defTxt, ...
      'ValueChangedFcn', @(~,~) extraCallback());
end
end


function browseIntoField(txtField, onEdit)
[file, path] = uigetfile({'*.msh;*.vtk;*', 'Mesh Files'}, 'Select Mesh File');
if isequal(file, 0)
  return;
end
txtField.Value = fullfile(path, file);
onEdit();
end


function tf = isControlModified(c)
if c.explicitlyLoaded
  tf = true;
  return;
end
val = readWidgetValue(c.widget, c.param);
def = c.initialValue;
if isempty(val) && isempty(def)
  tf = false;
elseif isempty(val) ~= isempty(def)
  tf = true;
elseif (isnumeric(val) || islogical(val)) && (isnumeric(def) || islogical(def))
  tf = ~isequal(size(val), size(def)) || ...
    any(abs(double(val(:)) - double(def(:))) > 1e-12 * max(1, abs(double(def(:)))));
else
  tf = ~strcmp(string(val), string(def));
end
end


function [isValid, missingNames] = validateControlGroup(ctrlMap)
isValid = true;
missingNames = string.empty;
fNames = fieldnames(ctrlMap);

for i = 1:numel(fNames)
  c = ctrlMap.(fNames{i});
  p = c.param;
  if strcmp(c.rowGrid.Visible, 'off')
    continue;
  end

  val = readWidgetValue(c.widget, p);
  bad = false;
  if p.required && isempty(val)
    bad = true;
  elseif p.type == "vector" && isempty(val) && ...
      strlength(strtrim(string(c.widget.Value))) > 0 && ~strcmpi(strtrim(c.widget.Value), 'all')
    bad = true;
  end

  if bad
    isValid = false;
    missingNames(end+1) = p.name; %#ok<AGROW>
    c.nameLabel.FontColor = [0.85 0.1 0.1];
  else
    c.nameLabel.FontColor = [0 0 0];
  end
end
end


function out = collectGroupValues(ctrlMap, forceInclude)
if nargin < 2
  forceInclude = false;
end
out = struct();
fNames = fieldnames(ctrlMap);
if isempty(fNames)
  return;
end

shouldExport = forceInclude;
if ~shouldExport
  for i = 1:numel(fNames)
    c = ctrlMap.(fNames{i});
    if strcmp(c.rowGrid.Visible, 'off')
      continue;
    end
    if c.param.required || isControlModified(c)
      shouldExport = true;
      break;
    end
  end
end

if ~shouldExport
  return;
end

for i = 1:numel(fNames)
  fn = fNames{i};
  c = ctrlMap.(fn);
  if strcmp(c.rowGrid.Visible, 'off')
    continue;
  end
  val = readWidgetValue(c.widget, c.param);
  if ~isempty(val)
    out.(fn) = val;
  end
end
end


function out = collectChronosValues(chronosGroup)
out = struct();
secNames = fieldnames(chronosGroup);
for i = 1:numel(secNames)
  sec = secNames{i};
  secVals = collectGroupValues(chronosGroup.(sec), false);
  if ~isempty(fieldnames(secVals))
    out.(sec) = secVals;
  end
end
end


function val = readWidgetValue(widget, p)
switch p.type
  case "domainRef"
    val = str2double(string(widget.Value));
  case "enum"
    val = string(widget.Value);
  case "editableEnum"
    raw = strtrim(string(widget.Value));
    if raw == ""
      val = [];
    else
      val = raw;
    end
  case "flag01"
    val = double(widget.Value);
  case "double"
    val = widget.Value;
  case "vector"
    raw = strtrim(string(widget.Value));
    if raw == "" || strcmpi(raw, "all")
      val = [];
    else
      num = str2num(char(strrep(raw, ';', ','))); %#ok<ST2NM>
      if isnumeric(num) && ~isempty(num)
        val = num;
      else
        val = [];
      end
    end
  case "file"
    txtChild = findobj(widget, 'Type', 'uieditfield');
    raw = strtrim(string(txtChild.Value));
    if raw == ""
      val = [];
    else
      val = raw;
    end
  otherwise
    raw = strtrim(string(widget.Value));
    if raw == ""
      val = [];
    else
      val = raw;
    end
end
end


function ctrlMap = applyStructToGroup(ctrlMap, dataStruct)
if ~isstruct(dataStruct)
  return;
end
fNames = fieldnames(ctrlMap);
for i = 1:numel(fNames)
  fn = fNames{i};
  if ~isfield(dataStruct, fn)
    continue;
  end
  val = dataStruct.(fn);
  if isempty(val) || any(ismissing(val), 'all')
    continue;
  end

  c = ctrlMap.(fn);
  p = c.param;
  c.explicitlyLoaded = true;

  switch p.type
    case "domainRef"
      sVal = formatValueForWidgetText(val);
      if ~any(strcmp(c.widget.Items, sVal))
        c.widget.Items = [c.widget.Items, {sVal}];
      end
      c.widget.Value = sVal;
    case "enum"
      sVal = formatValueForWidgetText(val);
      if strcmp(fn, 'entityListType') && strcmpi(sVal, 'tag')
        sVal = 'tags';
      end
      matchIdx = find(strcmpi(c.widget.Items, sVal), 1);
      if ~isempty(matchIdx)
        c.widget.Value = c.widget.Items{matchIdx};
      end
    case "editableEnum"
      sVal = formatValueForWidgetText(val);
      matchIdx = find(strcmpi(c.widget.Items, sVal), 1);
      if ~isempty(matchIdx)
        c.widget.Value = c.widget.Items{matchIdx};
      else
        c.widget.Items = [c.widget.Items, {sVal}];
        c.widget.Value = sVal;
      end
    case "flag01"
      if islogical(val) || isnumeric(val)
        c.widget.Value = logical(val(1));
      else
        sv = lower(strtrim(string(val)));
        c.widget.Value = (sv == "true" || sv == "1");
      end
    case "double"
      if isnumeric(val) && isscalar(val)
        c.widget.Value = double(val);
      elseif isstring(val) || ischar(val)
        nVal = str2double(string(val));
        if ~isnan(nVal)
          c.widget.Value = nVal;
        end
      end
    case "file"
      txtChild = findobj(c.widget, 'Type', 'uieditfield');
      txtChild.Value = formatValueForWidgetText(val);
    otherwise
      c.widget.Value = formatValueForWidgetText(val);
  end
  ctrlMap.(fn) = c;
end
end


function txt = formatValueForWidgetText(val)
if isempty(val)
  txt = '';
elseif islogical(val)
  parts = arrayfun(@(x) sprintf('%d', double(x)), val(:)', 'UniformOutput', false);
  txt = strjoin(parts, ',');
elseif isnumeric(val)
  parts = arrayfun(@(x) sprintf('%.15g', x), val(:)', 'UniformOutput', false);
  txt = strjoin(parts, ',');
elseif isstring(val) || ischar(val) || iscellstr(val)
  v = string(val);
  txt = char(strjoin(v(:), ","));
else
  txt = char(string(val));
end
end


function chronosGroup = applyChronosStruct(chronosGroup, physStruct)
secNames = fieldnames(chronosGroup);
for i = 1:numel(secNames)
  sec = secNames{i};
  if isfield(physStruct, sec) && isstruct(physStruct.(sec))
    chronosGroup.(sec) = applyStructToGroup(chronosGroup.(sec), physStruct.(sec));
  end
end
end


function [node, nodeDir, visited] = resolveXmlNodeWithDir(node, baseDir, visited)
if nargin < 3
  visited = string.empty;
end
nodeDir = baseDir;
if isstruct(node) && isscalar(node) && isfield(node, 'fileName')
  refRaw = char(string(node.fileName));
  if ~startsWith(refRaw, '/') && ~startsWith(refRaw, '\') && ~contains(refRaw, ':')
    cand = fullfile(baseDir, refRaw);
    if isfile(cand)
      refPath = cand;
    elseif ~isempty(baseDir) && isfile(fullfile(fileparts(baseDir), refRaw))
      refPath = fullfile(fileparts(baseDir), refRaw);
    else
      refPath = cand;
    end
  else
    refPath = refRaw;
  end
  canonRef = canonicalFilePath(refPath);
  if strlength(canonRef) > 0 && ~any(visited == canonRef)
    visited(end+1) = canonRef;
  end
  node = readInput(refPath);
  [nodeDir, ~, ~] = fileparts(refPath);
end
end


function cp = canonicalFilePath(p)
d = dir(char(p));
if ~isempty(d) && ~d(1).isdir
  cp = string(fullfile(d(1).folder, d(1).name));
else
  cp = string(p);
end
end


function tag = readXmlRootTag(xmlPath)
txt = fileread(char(xmlPath));
txt = regexprep(txt, '<\?[\s\S]*?\?>|<!--[\s\S]*?-->|<!DOCTYPE[\s\S]*?>', '');
tok = regexp(txt, '<\s*([A-Za-z_][\w.:-]*)', 'tokens', 'once');
if isempty(tok)
  tag = "";
else
  tag = string(tok{1});
end
end


function pr = getXmlLoadPriority(tag)
switch tag
  case "Problem"
    pr = 1;
  case "Domain"
    pr = 2;
  case {"SimulationParameters", "LinearSolver"}
    pr = 3;
  case "Solver"
    pr = 4;
  case {"Geometry", "Materials", "BoundaryConditions"}
    pr = 5;
  case {"Interface", "Output"}
    pr = 6;
  otherwise
    pr = 7;
end
end


function out = mergeStructFields(a, b)
out = a;
if ~isstruct(b)
  return;
end
fn = fieldnames(b);
for i = 1:numel(fn)
  out.(fn{i}) = b.(fn{i});
end
end


function setParamRowVisible(ctrl, isVis)
if isVis
  ctrl.rowGrid.Visible = 'on';
  ctrl.parentGrid.RowHeight{ctrl.rowIdx} = 'fit';
else
  ctrl.rowGrid.Visible = 'off';
  ctrl.parentGrid.RowHeight{ctrl.rowIdx} = 0;
end
end


function setGridRowVisible(parentGrid, sec, rowIdx, isVis)
if isVis
  sec.Container.Visible = 'on';
  parentGrid.RowHeight{rowIdx} = 'fit';
else
  sec.Container.Visible = 'off';
  parentGrid.RowHeight{rowIdx} = 0;
end
end


function str = boolOnOff(tf)
if tf
  str = 'on';
else
  str = 'off';
end
end


function ensureOutputParentDir(outStruct)
if ~isstruct(outStruct) || ~isfield(outStruct, 'matFileName')
  return;
end
matDir = fileparts(char(string(outStruct.matFileName)));
if ~isempty(matDir) && ~isfolder(matDir)
  mkdir(matDir);
end
end
