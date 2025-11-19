classdef Poisson < PhysicsSolver
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    F
    analyticalSolution
    fieldId
  end


  methods
    function obj = Poisson(domain)
      obj@PhysicsSolver(domain)
    end

    function registerSolver(obj,solverInput)

      nTags = obj.mesh.nCellTag;

      if ~isempty(solverInput)
        targetRegions = getXMLData(solverInput,1:nTags,"targetRegions");
      else
        targetRegions = 1:nTags;
      end

      % register nodal displacements on target regions
      obj.dofm.registerVariable(obj.getField(),entityField.node,1,targetRegions);

      % store the id of the field in the degree of freedom manager
      obj.fieldId = obj.dofm.getVariableId(obj.variableName);

      % initialize the state object
      n = getNumbDoF(obj.dofm,obj.fieldId);

      state = getState(obj);
      state.data.u = zeros(n,1);
      state.data.err = zeros(n,1);


    end

    function assembleSystem(obj,varargin)

      obj.J{obj.fldId,obj.fldId} = computeMat(obj);
      obj.rhs{obj.fldId} = computeRhs(obj);

    end

    function J = computeMat(obj,varargin)

      subCells = obj.dofm.getFieldCells(obj.fieldId);
      n = sum(obj.mesh.cellNumVerts(subCells).^2);
      Ndof = obj.dofm.getNumDoF(obj.field);
      asbJ = assembler(n,Ndof,Ndof, @(el) computeLocalMatrix(obj,el));
      % loop over cells
      for el = subCells'
        % get dof id and local matrix
        asbJ.localAssembly(el);
      end
      % populate stiffness matrix
      J = asbJ.sparseAssembly();
    end


    function [dofr,dofc,matLoc] = computeLocalMatrix(obj,elID)
      elem = getElementByID(obj.elements,elID);
      [gradN,dJW] = getDerBasisFAndDet(elem,elID,1);
      matLoc = pagemtimes(gradN,'ctranspose',gradN,'none');
      matLoc = matLoc.*reshape(dJW,1,1,[]);
      matLoc = sum(matLoc,3);
      % get global DoF
      nodes = obj.mesh.cells(elID,1:obj.mesh.cellNumVerts(elID));
      dof = obj.dofm.getLocalDoF(nodes,obj.fldId);
      dofr = dof; dofc = dof;
    end

    function rhs = computeRhs(obj,varargin)
      ents = obj.dofm.getActiveEnts(obj.getField());
      f = computeForcingTerm(obj);
      rhs = obj.J{obj.fldId,obj.fldId}*obj.state.data.u(ents) + f(ents);
    end

    function initState(obj)
      % add poromechanics fields to state structure
      obj.state.data.u = zeros(obj.mesh.nNodes,1);
      obj.state.data.err = zeros(obj.mesh.nNodes,1);
    end

    function updateState(obj,dSol)
      % Update state structure with last solution increment
      ents = obj.dofm.getActiveEntities(obj.fieldId);
      if nargin > 1
        % update solution
        obj.state.data.u(ents) = obj.state.data.u(ents) + dSol(getDoF(obj.dofm,obj.fieldId));
      end

      % compute error if analytical solution is available
      if ~isempty(obj.analyticalSolution)
        analSol = computeAnal(obj,ents,'u',0);
        obj.state.data.err(ents) = obj.state.data.u(ents) - analSol;
      end
    end


    function [dof,vals] = getBC(obj,bcId,t)

      if ~strcmp(obj.bcs.getCond(bcId),"NodeBC")
          error('BC entitiy %s is not available for %s field',cond,obj.field);
      end

      dof = obj.bcs.getBCentities(bcId);
      vals = obj.bcs.getVals(bcId,t);

    end

    function applyDirVal(obj,dof,vals)
      obj.state.data.u(dof) = vals;
    end

    function [cellData,pointData] = printState(obj,sOld,sNew,t)
      % append state variable to output structure
      switch nargin
        case 2
          var = sOld.data.u;
        case 4
          % linearly interpolate state variables containing print time
          fac = (t - sOld.t)/(sNew.t - sOld.t);
          var = sNew.data.u*fac+sOld.data.u*(1-fac);
        otherwise
          error('Wrong number of input arguments');
      end
      [cellData,pointData] = buildPrintStruct(obj,var);
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
        pointStr(2).data = abs(var-computeAnal(obj,1:obj.mesh.nNodes,'u',0));
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
        x = obj.mesh.coordinates(list,:);
      elseif mode == 1
        x = list;
      else 
        error('Mode input must be 0 (nodes list) or 1 (coordinate input list)');
      end
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
      anal = arrayfun(@(i) f(x(i,:)),1:size(x,1));
      anal = reshape(anal,[],1);      % make column array
    end

    function F = computeForcingTerm(obj)
      % classical
      % general sparse assembly loop over elements for Poromechanics
      subCells = obj.dofm.getFieldCells(obj.field);
      F = zeros(obj.dofm.getNumDoF('Poisson'),1);
      % loop over cells
      for el = subCells'
        % get dof id and local matrix
        id = obj.mesh.cells(el,:);
        elem = getElementByID(obj.elements,el);
        c_gp = getGPointsLocation(elem,el);
        f_gp = computeAnal(obj,c_gp,'f',1);
        [~,dJW] = getDerBasisFAndDet(elem,el,1);
        N = getBasisFinGPoints(elem);
        floc = N'*(f_gp.*dJW');
        % proper integration
        F(id) = F(id) + floc;
      end
    end

    function [L2err,H1err] = computeError(obj)
      assert(~isempty(obj.analyticalSolution),['Missing analytical solution for ' ...
        'Poisson model \n'])
      L2err = 0;
      H1err = 0;
      for el = 1:obj.mesh.nCells
        vtkId = obj.mesh.cellVTKType(el);
        elem = getElement(obj.elements,vtkId);
        N = getBasisFinGPoints(elem);
        [gradN,dJW] = getDerBasisFAndDet(elem,el,1);
        dofId = obj.mesh.cells(el,:);
        locErr = (obj.state.data.err(dofId));
        norm_e = N*locErr.^2;        % squared value of error
        L2errLoc = sum(norm_e.*reshape(dJW,[],1));
        grad_e = pagemtimes(gradN,locErr);
        ngp = size(grad_e,3);
        norm_grad_e = zeros(ngp,1);
        for i = 1:ngp
          norm_grad_e(i) = norm(grad_e(:,:,i),2);
        end
        norm_e_grad_e = norm_e + norm_grad_e;
        H1errLoc = sum(norm_e_grad_e.*reshape(dJW,[],1));
        L2err = L2err + L2errLoc;
        H1err = H1err + H1errLoc;
      end
      L2err = sqrt(L2err);
      H1err = sqrt(H1err);
    end


  function [L2err,H1err] = computeError_v2(obj)
    assert(~isempty(obj.analyticalSolution),['Missing analytical solution for ' ...
      'Poisson model \n'])
    L2err = 0;
    H1err = 0;
    for el = 1:obj.mesh.nCells
      vtkId = obj.mesh.cellVTKType(el);
      elem = getElement(obj.elements,vtkId);
      N = getBasisFinGPoints(elem);
      [gradN,dJW] = getDerBasisFAndDet(elem,el,1);
      c_gp = getGPointsLocation(elem,el);
      dofId = obj.mesh.cells(el,:);
      u_gp = N*obj.state.data.u(dofId);
      err = u_gp - computeAnal(obj,c_gp,'u',1);
      err_2 = err.^2;        % squared value of error
      L2errLoc = sum(err_2.*reshape(dJW,[],1));

      grad_uh = pagemtimes(gradN,obj.state.data.u(dofId));
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

  methods (Static)
    function out = getField()
      out = "poissonField";
    end
  end
end

