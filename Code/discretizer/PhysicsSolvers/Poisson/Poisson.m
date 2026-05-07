classdef Poisson < PhysicsSolver
  % FEM solver for Poisson equation
  
  properties
    F
    analyticalSolution
    fieldId
    gaussOrder
  end


  methods
    function obj = Poisson(domain)
      obj@PhysicsSolver(domain)
    end

    function registerSolver(obj,varargin)

      nTags = obj.grid.cells.nTag;

      default = struct('targetRegions',1:nTags,'gaussOrder',0);

      params = readInput(default,varargin{:});

      obj.gaussOrder = params.gaussOrder;
      dofm = obj.domain.dofm;

      % register scalar poisson variable on target regions
      dofm.registerVariable(obj.getField(),entityField.node,1,params.targetRegions);

      % store the id of the field in the degree of freedom manager
      obj.fieldId = dofm.getVariableId(obj.getField());

      % initialize the state object
      n = getNumbDoF(dofm,obj.fieldId);

      s = getState(obj);
      s.u = zeros(n,1);
      s.err = zeros(n,1);
      setState(obj,s);

    end

    function assembleSystem(obj,varargin)

      obj.domain.J{obj.fieldId,obj.fieldId} = computeMat(obj);
      obj.domain.rhs{obj.fieldId} = computeRhs(obj);

    end

    function J = computeMat(obj,varargin)

      subCells = obj.domain.dofm.getFieldCells(obj.fieldId);
      cells = obj.grid.cells;
      
      n = sum(cells.numVerts(subCells).^2);
      Ndof = obj.domain.dofm.getNumbDoF(obj.getField());
      asbJ = assembler(n,Ndof,Ndof);
      coordinates = obj.grid.coordinates;
      cells = obj.grid.cells;

      % loop over cells
      for vtkId = cells.vtkTypes

        tmp = obj.grid.getCellsByVTKId(vtkId);
        subCellsLoc = intersect(subCells,tmp,'sorted');
        elem = FiniteElementType.create(vtkId,obj.grid,'gaussOrder',obj.gaussOrder);

        % get node topology for given vtk type
        topol = obj.grid.getCellNodes(subCellsLoc);

        for i = 1:numel(subCellsLoc)
          el = subCellsLoc(i);
          nodes = topol(el,:);
          coords = coordinates(nodes,:);
          % get dof id and local matrix
          [gradN,dJW] = getDerBasisFAndDet(elem,coords);
          Jloc = pagemtimes(gradN,'ctranspose',gradN,'none');
          Jloc = Jloc.*reshape(dJW,1,1,[]);
          Jloc = sum(Jloc,3);
          % get global DoF
          dof = obj.domain.dofm.getLocalDoF(obj.fieldId,nodes);
          asbJ.localAssembly(dof,dof,Jloc);
        end
      end
      % populate stiffness matrix
      J = asbJ.sparseAssembly();
    end

    function rhs = computeRhs(obj,varargin)
      ents = obj.domain.dofm.getActiveEntities(obj.getField());
      J = getJacobian(obj);
      f = computeForcingTerm(obj);
      u = getState(obj,"u");
      rhs = J*u(ents) + f(ents);
    end


    function updateState(obj,dSol)
      % Update state structure with last solution increment
      ents = obj.domain.dofm.getActiveEntities(obj.fieldId);
      if nargin > 1
        % update solution
        u = getState(obj,"u");
        u(ents) = u(ents) + dSol(getDoF(obj.domain.dofm,obj.fieldId));
        setState(obj,u,"u");
      end

      % compute error if analytical solution is available
      if ~isempty(obj.analyticalSolution)
        analSol = computeAnal(obj,ents,'u',0);
        err = u(ents) - analSol;
        setState(obj,err,"err");
      end
    end

    % function applyBC(obj,bcId,t)
    % 
    % 
    %   if ~strcmp(bc.getField(bcId),"node")
    %     error('BC entitiy %s is not available for %s field',cond,obj.getField());
    %   end
    % 
    %   applyBC@PhysicsSolver(obj,bcId,t);
    % 
    %   % % get bcDofs and bcVals
    %   % [bcDofs,bcVals] = getBC(obj,bcId,t);
    %   % 
    %   % bcType = obj.domain.bcs.getType(bcId);
    %   % 
    %   % switch bcType
    %   %   case 'Dirichlet'
    %   %     applyDirBC(obj,bcId,bcDofs);
    %   %   case 'Neumann'
    %   %     applyNeuBC(obj,bcId,bcDofs,bcVals);
    %   %   otherwise
    %   %     error("Error in %s: Boundary condition type '%s' is not " + ...
    %   %       "available in %s",class(obj),bcType);
    %   % end
    % 
    % end

    function advanceState(obj)
      % do nothing
    end

    % 
    % function [dof,vals] = getBC(obj,bcId,t)
    % 
    %   bc = obj.domain.bcs;
    %   dof = bc.getBCentities(bcId);
    %   vals = bc.getVals(bcId,t);
    % 
    % end

    % function applyDirVal(obj,bcId,t)
    % 
    %   [bcDofs,bcVals] = getBC(obj,bcId,t);
    % 
    %   obj.domain.state.data.u(bcDofs) = bcVals;
    % 
    % end

    function [cellData,pointData] = writeVTK(obj,fac,varargin)

      sOld = getStateOld(obj);
      sNew = getState(obj);

      u = sNew.data.u*fac+sOld.data.u*(1-fac);

      [cellData,pointData] = buildPrintStruct(obj,u);

    end

    function writeSolution(obj,fac,tID)

      uOld = getStateOld(obj,obj.getField());
      uCurr = getState(obj,obj.getField());

      obj.domain.outstate.results(tID).(obj.getField()) = uCurr*fac+uOld*(1-fac);
    
    end



    function [cellStr,pointStr] = buildPrintStruct(obj,var)
      nPointData = 1;
      if ~isempty(obj.analyticalSolution)
        nPointData = nPointData + 1;
      end
      pointStr = repmat(struct('name', 1, 'data', 1), nPointData, 1);
      cellStr = [];
      % Displacement
      pointStr(1).name = 'u';
      pointStr(1).data = var;
      if ~isempty(obj.analyticalSolution)
        pointStr(2).name = 'abs_error';
        pointStr(2).data = abs(var-computeAnal(obj,1:obj.grid.nNodes,'u',0));
      end
    end

    function setAnalSolution(obj,u,f,dux,duy,duz)
      %
      obj.analyticalSolution.u = u;
      obj.analyticalSolution.f = f;
      %
      if nargin > 3
        assert(nargin == 6,'Incorrect number of input.');
        obj.analyticalSolution.dux = dux;
        obj.analyticalSolution.duy = duy;
        obj.analyticalSolution.duz = duz;
      end
    end

    function anal = computeAnal(obj,list,var,mode)
      if mode == 0
        x = obj.grid.coordinates(list,:);
      elseif mode == 1
        x = list;
      else 
        error('Mode input must be 0 (nodes list) or 1 (coordinate input list)');
      end

      [x,y,z] = deal(x(:,1),x(:,2),x(:,3));
      switch var
        case 'u'
          f = obj.analyticalSolution.u;
        case 'f'
          f = obj.analyticalSolution.f;
        case 'grad_x'
          f = obj.analyticalSolution.dux;
        case 'grad_y'
          f = obj.analyticalSolution.duy;
        case 'grad_z'
          f = obj.analyticalSolution.duz;
        otherwise
          error('Incorrect input for flag')
      end
      anal = f(x,y,z);
      anal = reshape(anal,[],1);      % make column array
    end

    function F = computeForcingTerm(obj)
      % classical
      % general sparse assembly loop over elements for Poromechanics
      subCells = obj.domain.dofm.getFieldCells(obj.getField());
      cells = obj.grid.cells;
      coordinates = obj.grid.coordinates;
      F = zeros(obj.domain.dofm.getNumbDoF(obj.getField()),1);

      for vtkId = cells.vtkTypes

        subCellsLoc = obj.grid.getCellsByVTKId(vtkId,subCells);
        elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);

        % get node topology for given vtk type
        topol = obj.grid.getCellNodes(subCellsLoc);

        for i = 1:numel(subCellsLoc)
          el = subCellsLoc(i);
          nodes = topol(el,:);
          coords = coordinates(nodes,:);
          c_gp = getGPointsLocation(elem,el);
          f_gp = computeAnal(obj,c_gp,'f',1);
          [~,dJW] = getDerBasisFAndDet(elem,coords);
          N = getBasisFinGPoints(elem);
          floc = N'*(f_gp.*dJW');
          % proper integration
          F(nodes) = F(nodes) + floc;
        end
      end

    end



    function [L2err,H1err] = computeError(obj)
      assert(~isempty(obj.analyticalSolution),['Missing analytical solution for ' ...
        'Poisson model \n'])
      L2err = 0;
      H1err = 0;

      u = getState(obj,"u");

      coordinates = obj.grid.coordinates;
      cells = obj.grid.cells;

      for vtkId = cells.vtkTypes

        cellList = obj.grid.getCellsByVTKId(vtkId);
        elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);

        % get node topology for given vtk type
        topol = obj.grid.getCellNodes(cellList);

        for i = 1:numel(cellList)

          el = cellList(i);
          nodes = topol(el,:);
          coords = coordinates(nodes,:);
          N = getBasisFinGPoints(elem);
          [gradN,dJW] = getDerBasisFAndDet(elem,coords);
          c_gp = getGPointsLocation(elem,el);
          u_gp = N*u(nodes);
          err = u_gp - computeAnal(obj,c_gp,'u',1);
          err_2 = err.^2;        % squared value of error
          L2errLoc = sum(err_2.*reshape(dJW,[],1));
          grad_uh = pagemtimes(gradN,u(nodes));
          grad_uh = squeeze(permute(grad_uh,[3 1 2]));

          grad_u = [computeAnal(obj,c_gp,'grad_x',1),...
                    computeAnal(obj,c_gp,'grad_y',1),...
                    computeAnal(obj,c_gp,'grad_z',1)];

          grad_err = grad_uh - grad_u;
          grad_err2 = sum(grad_err.^2,2);
          semiH1errLoc = sum(grad_err2.*reshape(dJW,[],1));

          L2err = L2err + L2errLoc;
          H1err = H1err + L2errLoc + semiH1errLoc;
        end

        L2err = sqrt(L2err);
        H1err = sqrt(H1err);
      end
    end

      function out = isSymmetric(obj)
        out = true;
      end
    end

    methods (Static)
      function out = getField()
      out = "u";
    end
  end
end

