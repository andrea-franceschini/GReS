classdef VariablySaturatedFlowFEM < SinglePhaseFlowFEM
  % VARIABLYSATURATEDFLOWFEM Variably Saturated flow solver using FEM
  % Implements Richards equations for unsaturated flow using Finite Elements

  properties
    NLscheme      % newton or picard
  end

  methods (Access = public)
    function obj = VariablySaturatedFlowFEM(domain)
      obj@SinglePhaseFlowFEM(domain);
    end

    function registerSolver(obj,varargin)
      registerSolver@SinglePhaseFlowFEM(obj,varargin{:});
      state = getState(obj);
      state.saturation = zeros(obj.grid.cells.num,1);
      setState(obj,state);

      input = readInput(struct('NLscheme',"newton"),varargin{:});
      obj.NLscheme = input.NLscheme;
    end

    function assembleSystem(obj,dt)
      dofm = obj.domain.dofm;
      materials = obj.domain.materials;
      coordinates = obj.grid.coordinates;
      cells = obj.grid.cells;
      subCells = dofm.getFieldCells(obj.fieldId);
      Ndof = dofm.getNumbDoF(obj.fieldId);
      ents = dofm.getActiveEntities(obj.fieldId);
      
      p = obj.getState(obj.getField());
      pOld = [];
      if ~obj.steadyState
        pOld = obj.getStateOld(obj.getField());
      end
      
      gamma = materials.getFluid().getSpecificWeight();
      mu = materials.getFluid().getDynViscosity();
      beta = materials.getFluid().getFluidCompressibility();
      coupledRegions = dofm.getTargetRegions([obj.getField(),"displacements"]);
      
      nEntries = sum(cells.numVerts(subCells).^2);
      asbH = assembler(nEntries, Ndof, Ndof);
      asbP = assembler(nEntries, Ndof, Ndof);
      asbJh = assembler(nEntries, Ndof, Ndof);
      asbJp = assembler(nEntries, Ndof, Ndof);
      
      for vtkId = cells.vtkTypes
        cellList = obj.grid.getCellsByVTKId(vtkId, subCells);
        elem = FiniteElementType.create(vtkId, obj.grid, obj.gaussOrder);
        N = getBasisFinGPoints(elem);
        
        for tag = 1:cells.nTag
          cellListTag = cellList(cells.tag(cellList)==tag);
          if isempty(cellListTag), continue; end
          
          mat = materials.getMaterial(tag);
          permMat = mat.PorousRock.getPermMatrix() / mu;
          poro = mat.PorousRock.getPorosity();
          alpha = getRockCompressibility(obj, tag, coupledRegions);
          Sws = mat.PorousRock.getMaxSaturation();
          Swr = mat.PorousRock.getResidualSaturation();
          
          topol = obj.grid.getCellNodes(cellListTag);
          
          for i = 1:numel(cellListTag)
            nodes = topol(i,:);
            coords = coordinates(nodes,:);
            [gradN, dJWeighed] = getDerBasisFAndDet(elem, coords);
            
            pLoc = p(nodes);
            
            HLoc = zeros(numel(nodes));
            PLoc = zeros(numel(nodes));
            JhLoc = zeros(numel(nodes));
            JpLoc = zeros(numel(nodes));
            
            for gp = 1:numel(dJWeighed)
              N_gp = N(gp,:);
              gradN_gp = gradN(:,:,gp);
              dV = dJWeighed(gp);
              
              p_gp = N_gp * pLoc;
              
              [Sw_raw, dSw_raw, d2Sw_raw] = mat.PorousRock.Curves.computeSwAnddSw(p_gp);
              Sw = Swr + (Sws - Swr) * Sw_raw;
              dSw = -(Sws - Swr) * dSw_raw; 
              d2Sw = (Sws - Swr) * d2Sw_raw;
              
              [kr, dkr_raw] = mat.PorousRock.Curves.computeRelativePermeability(p_gp);
              dkr = -dkr_raw;
              
              H_gp = gradN_gp' * (kr * permMat) * gradN_gp * dV;
              HLoc = HLoc + H_gp;
              
              if ~obj.steadyState
                C_gp = alpha * Sw + poro * beta * Sw + poro * dSw;
                P_gp = C_gp * diag(N_gp) * dV; 
                PLoc = PLoc + P_gp;
              end
              
              if obj.isNewtonNLSolver()
                gradP_gp = gradN_gp * pLoc;
                if gamma > 0
                  gradP_gp = gradP_gp + [0; 0; gamma];
                end
                Jh_gp = gradN_gp' * (dkr * permMat) * gradP_gp * N_gp * dV;
                JhLoc = JhLoc + Jh_gp;
                
                if ~obj.steadyState
                  pOldLoc = pOld(nodes);
                  p_gp_old = N_gp * pOldLoc;
                  pRate_gp = (p_gp - p_gp_old) / dt;
                  dC_dp = beta * alpha * Sw + (2 * alpha + beta * poro) * dSw + poro * d2Sw;
                  Jp_gp = dC_dp * pRate_gp * diag(N_gp) * dV;
                  JpLoc = JpLoc + Jp_gp;
                end
              end
            end
            
            dof = dofm.getLocalDoF(obj.fieldId, nodes);
            asbH.localAssembly(dof, dof, HLoc);
            if ~obj.steadyState
              asbP.localAssembly(dof, dof, PLoc);
            end
            if obj.isNewtonNLSolver()
              asbJh.localAssembly(dof, dof, JhLoc);
              if ~obj.steadyState
                asbJp.localAssembly(dof, dof, JpLoc);
              end
            end
          end
        end
      end
      
      obj.H = asbH.sparseAssembly();
      if ~obj.steadyState
        obj.P = asbP.sparseAssembly();
      else
        obj.P = sparse(Ndof, Ndof);
      end
      
      entType = dofm.getFieldLocation(obj.fieldId);
      z = getLocation(entType, obj.grid, ents);
      z = z(:,3);
      
      if obj.steadyState
        rhs = obj.H * (p(ents) + gamma * z);
        J = obj.H;
      else
        rhs = obj.H * (p(ents) + gamma * z) + (obj.P / dt) * (p(ents) - pOld(ents));
        J = obj.H + obj.P / dt;
      end
      
      if obj.isNewtonNLSolver()
        Jh = asbJh.sparseAssembly();
        J = J + Jh;
        if ~obj.steadyState
          Jp = asbJp.sparseAssembly();
          J = J + Jp;
        end
      end
      
      obj.domain.J{obj.fieldId, obj.fieldId} = J;
      obj.domain.rhs{obj.fieldId} = rhs;
    end

    function flux = computeFlux(obj,p,varargin)
      dofm = obj.domain.dofm;
      materials = obj.domain.materials;
      coordinates = obj.grid.coordinates;
      cells = obj.grid.cells;
      subCells = dofm.getFieldCells(obj.fieldId);
      gamma = materials.getFluid().getSpecificWeight();
      mu = materials.getFluid().getDynViscosity();

      flux = zeros(3*cells.num,1);

      for vtkId = cells.vtkTypes
        cellList = obj.grid.getCellsByVTKId(vtkId,subCells);
        elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);
        N = getBasisFinGPoints(elem);

        for tag = 1:cells.nTag
          cellListTag = cellList(cells.tag(cellList)==tag);
          if isempty(cellListTag), continue; end
          mat = materials.getMaterial(tag);
          permMat = mat.PorousRock.getPermMatrix() / mu;

          topol = obj.grid.getCellNodes(cellListTag);

          for i = 1:numel(cellListTag)
            nodes = topol(i,:);
            coords = coordinates(nodes,:);
            el = cellListTag(i);

            [gradN,dJWeighed] = getDerBasisFAndDet(elem,coords);
            pLoc = p(nodes);
            
            fluxAvg = zeros(3,1);
            totalVol = sum(dJWeighed);
            
            for gp = 1:numel(dJWeighed)
              N_gp = N(gp,:);
              gradN_gp = gradN(:,:,gp);
              dV = dJWeighed(gp);
              
              p_gp = N_gp * pLoc;
              [kr, ~] = mat.PorousRock.Curves.computeRelativePermeability(p_gp);
              
              gradP = gradN_gp * pLoc;
              if gamma > 0
                gradP = gradP + [0; 0; gamma];
              end
              
              fluxGP = (kr * permMat) * gradP;
              fluxAvg = fluxAvg + fluxGP * dV;
            end
            
            flux(3*el-2:3*el,1) = fluxAvg / totalVol;
          end
        end
      end
    end

    function Sw = computeSaturation(obj,p)
      cells = obj.grid.cells;
      Sw = zeros(cells.num,1);
      subCells = obj.domain.dofm.getFieldCells(obj.getField());
      
      for el = subCells'
        tag = cells.tag(el);
        mat = obj.domain.materials.getMaterial(tag);
        Sws = mat.PorousRock.getMaxSaturation();
        Swr = mat.PorousRock.getResidualSaturation();
        
        nodes = obj.grid.getCellNodes(el);
        pCell = mean(p(nodes));
        
        [Sw_raw, ~, ~] = mat.PorousRock.Curves.computeSwAnddSw(pCell);
        Sw(el) = Swr + (Sws - Swr) * Sw_raw;
      end
    end

    function states = finalizeState(obj,p,t)
      states = finalizeState@SinglePhaseFlowFEM(obj,p,t);
      states.saturation = computeSaturation(obj,p);
    end

    function updateState(obj,dSol)
      dofm = obj.domain.dofm;
      if nargin > 1
        ents = dofm.getActiveEntities(obj.getField());
        state = getState(obj);
        p = state.pressure;
        state.pressure(ents) = p(ents) + dSol(dofm.getDoF(obj.fieldId));
        state.saturation = computeSaturation(obj,state.pressure); 
        setState(obj,state);
      end
    end

    function writeSolution(obj,fac,tID)
      state = obj.domain.state.interpolate(fac);
      obj.domain.outstate.results(tID).pressure = state.pressure;
      obj.domain.outstate.results(tID).saturation = state.saturation;
    end

    function [cellStr,pointStr] = buildPrintStruct(obj,state)
      [cellStr, pointStr] = buildPrintStruct@SinglePhaseFlowFEM(obj,state);
      nCells = numel(cellStr);
      cellStr(nCells+1).name = 'saturation';
      cellStr(nCells+1).data = state.saturation;
    end

    function out = isFEM(obj)
      out = true;
    end

    function out = isTPFA(obj)
      out = false;
    end

    function str = typeDiscretization(obj)
      str = "FEM";
    end

    function out = isLinear(obj)
      out = false;
    end

    function out = isNewtonNLSolver(obj)
      out = strcmp("newton",obj.NLscheme);
    end
  end
end