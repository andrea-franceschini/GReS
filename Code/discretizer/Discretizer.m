classdef Discretizer < handle
   % General discretizer class
   properties (GetAccess=public, SetAccess=private)
     solver      % database for physics solvers in the model
     dofm        % dofManager
     numSolvers  % number of solvers discretized
     mod
     fields
     grid
   end

   properties (GetAccess=public, SetAccess=public)
     interfaceList = []; 
     interfaces = []
     interfaceSurf            % surfaceTag in each surface of interfaceSurf
     % empty - single domain simulation
     % not empty - call to mesh glue instances 
     state
   end

   properties (Access = private)
     solverMap
   end

   methods (Access = public)
      function obj = Discretizer(symmod,simParams,dofManager,grid,mat)
         %UNTITLED Construct an instance of this class
         %   Detailed explanation goes here
         obj.mod = symmod;
         obj.dofm = dofManager;
         obj.grid = grid;
         obj.solver = containers.Map('KeyType','double','ValueType','any');
         obj.setDiscretizer(symmod,simParams,dofManager,grid,mat);
         obj.checkTimeDependence(symmod,mat,simParams);
      end
      
      function applyBC(obj,bound,t,idDom)
        % Apply boundary condition to blocks of physical solver
        % ents: id of constrained entity
        % vals: value to assign to each entity
        bcList = bound.db.keys;
        % get entities and values of boundary condition
        for bc = string(bcList)
          field = bound.getPhysics(bc);
          % get id of constrained entities and corresponding BC values
          [bcEnts,bcVals] = getBC(getSolver(obj,field),bound,bc,t);

          %removeInterfaceBC(obj,bcEnts,bcVals)
          % apply Boundary conditions to each Jacobian/rhs block
          for f = obj.fields
            if ~isCoupled(obj,field,f)
              continue
              % skip pair of uncoupled physics
            end
            switch bound.getType(bc)
              case 'Dir'
                if nargin > 3
                  assert(~isempty(obj.interfaceList),['Too many input arguments: ' ...
                    'invalid domain id input for single domain BC imposition']);
                  for i = 1:length(obj.interfaceList)
                    [bcEnts,bcVals] = obj.interfaces(i).removeSlaveBCdofs(field,[bcEnts,bcVals],idDom);
                  end
                end
                applyDirBC(obj.getSolver({field,f}),field,bcEnts,bcVals);
              case {'Neu','VolumeForce'}
                applyNeuBC(obj.getSolver({field,f}),bcEnts,bcVals);
            end
          end
        end
      end

      function applyDirVal(obj,bound,t,idDom)
         % Apply boundary condition to blocks of physical solver
         % ents: id of constrained entity
         % vals: value to assign to each entity
         bcList = bound.db.keys;
         % get entities and values of boundary condition
         for bc = string(bcList)
            if ~strcmp(bound.getType(bc),'Dir')
               continue
            end
            field = bound.getPhysics(bc);
            [bcEnts,bcVals] = getBC(getSolver(obj,field),bound,bc,t);
            if nargin > 3
              assert(~isempty(obj.interfaceList),['Too many input arguments: ' ...
                'invalid domain id input for single domain BC imposition']);
              for i = 1:length(obj.interfaceList)
                [bcEnts,bcVals] = obj.interfaces(i).removeSlaveBCdofs(field,[bcEnts,bcVals],idDom);
              end
            end
            getSolver(obj,field).applyDirVal(bcEnts,bcVals);
         end
      end

      function J = assembleJacobian(obj)
         % put together jacobian blocks of SinglePhysicsSolver and
         % CoupledSolver in the model
         % Use the ordering specified in DoF manager class   
         switch obj.dofm.ordering
            case 'field'
               nFld = numel(obj.fields);
               J = cell(nFld,nFld);
               for i = 1:nFld
                  for j = 1:nFld
                     J{i,j} = getJacobian(obj.getSolver({obj.fields(i),obj.fields(j)}),obj.fields(i));
                  end
               end
            otherwise
               error('Invalid DoF manager ordering')
         end
      end

      function rhs = assembleRhs(obj)
         % put together rhs blocks of SinglePhysicsSolver and
         % CoupledSolver in the model
         nFld = numel(obj.fields);
         rhs = cell(nFld,1);
         for i = 1:nFld
            rhs{i} = zeros(getNumDoF(obj.dofm,obj.fields(i)),1);
            for j = 1:nFld
               rhs{i} = rhs{i} + ...
                  getRhs(getSolver(obj,{obj.fields(i),obj.fields(j)}),obj.fields(i));
            end  
         end
      end

      % function dSol = solve(obj,J,rhs)
      %    % assemble and solve whole linear system
      %    J = assembleJacobian(obj);
      %    rhs = assembleRhs(obj);
      %    dSol = J\-rhs;
      % end


      function out = getSolver(obj,fldList)
         fldList = string(fldList);
         % map single field or pair of field to db position
         if isscalar(obj.fields) % singlePhysic model
            v = 0;
         else
            nF = numel(obj.fields); % multiPhysic model
            cs = cumsum(nF:-1:1);
            v = [0 cs(1:end-1)];
         end
         fldList = unique(fldList);          
         % get solver from database
         if isscalar(fldList) % query to single physic solver
            fldId = obj.dofm.getFieldId(fldList);
            id = 1+v(fldId);
         else                 % query to coupled solver
            fldId = obj.dofm.getFieldId(fldList);
            fldId = sort(fldId);
            id = v(fldId(1))+fldId(2);
         end
         out = obj.solver(id);
      end


      function addInterface(obj,interfId,interf)
        if ~ismember(interfId,obj.interfaceList)
          obj.interfaceList = sort([obj.interfaceList interfId]);
          obj.interfaces = [obj.interfaces interf];
        end
      end


      function computeMatricesAndRhs(obj,stateOld,dt)
         % loop trough solver database and compute non-costant jacobian
         % blocks and rhs block
         for i = 1:obj.numSolvers
           computeMat(obj.solver(i),stateOld,dt);
           computeRhs(obj.solver(i),stateOld,dt);
         end
      end

      function out = isCoupled(obj,field1,field2)
         % check if input fields are coupled, i.e exist cells having both
         % fields activated
         sub1 = getActiveSubdomain(obj.dofm,field1);
         sub2 = getActiveSubdomain(obj.dofm,field2);
         out = any(intersect(sub1,sub2));
      end

      function initState(obj)
         % loop trough active single physics solver and update the state class
         % accordingly
         for i = 1:numel(obj.fields)
            % loop trough active fields and update the state structure
            initState(obj.getSolver(obj.fields(i)));
         end
      end

      function updateState(obj,du)
         % update current state
         for i = 1:numel(obj.fields)
            obj.getSolver(obj.fields(i)).updateState(du);
         end
      end
   end

   methods(Access = private)
      function setDiscretizer(obj,symmod,params,dofManager,grid,mat)
        obj.setSolverMap();
        flds = getFieldList(obj.dofm);
        nF = numel(flds);
        % loop over all fields and define corresponding models
        k = 0;
        % create the handle to state object that will be shared across all physical
        % modules
        stat = State();
        for i = 1:nF
          for j = i:nF
            k = k+1;
            addPhysics(obj,k,flds(i),flds(j),symmod,params,dofManager,grid,mat,stat);
          end
        end
        obj.state = stat;
        obj.fields = flds;
        obj.numSolvers = k;
      end

      function checkTimeDependence(obj,mod,mat,parm)
        % check if there is any time dependence in the input model
        % no time dependence in absence of flow and
        % incompressible single phase flow model.
        if ~isSinglePhaseFlow(mod)
          % Biot model is time dependent
          setTimeDependence(parm,false);
          return
        else
          % check if fluid is incompressible
          beta = getFluidCompressibility(mat.getFluid());
          if beta < eps
            setTimeDependence(parm,false);
          end
        end
      end

      function addPhysics(obj,id,f1,f2,mod,parm,dof,grid,mat,state)
        % Add new key to solver database
        % Prepare input fields for solver definition
        if ~isCoupled(obj,f1,f2)
          return
        end
        f = join(unique(sort({char(f1),char(f2)})));

        if ~obj.solverMap.isKey(string(f{:}))
          error('A physical module coupling %s with %s is not available.',f1,f2)
        else
          solv = obj.solverMap(f{:});
        end

        obj.solver(id) = solv(mod,parm,dof,grid,mat,state);
      end

      function setSolverMap(obj)

        obj.solverMap = containers.Map('KeyType','char','ValueType','any');

        subClasses = [findSubClasses('SinglePhysics','SinglePhysics'), ...
          findSubClasses('CouplingPhysics','CouplingPhysics')];

        for i = 1:numel(subClasses)
          obj.solverMap = feval([subClasses{i} '.registerSolver'],...
            obj.solverMap,subClasses{i});
        end
      end

   end

   methods (Static)
%      function [row,col,val,c] = computeLocalMatrix(mat,row,col,val,c,w,dofRow,dofCol)
%        % shortcut for assemblying local matrix contributions in sparse format
%        mat = mat.*reshape(w,1,1,[]);
%        mat = sum(mat,3);
%        n = numel(mat);
%        [J, I] = meshgrid(1:size(mat,2), 1:size(mat,1));
%        row(c+1:c+n) = dofRow(I);
%        col(c+1:c+n) = dofCol(J);
%        val(c+1:c+n) = mat(:);
%        c = c+n;
%      end

     function mat_new = expandMat(mat,n)
       % Get the size of the original matrix
       [s1, s2] = size(mat);

       % Initialize the sparse matrix: row indices, column indices, and values
       rows = [];
       cols = [];
       values = [];

       % Loop through the original matrix and populate the sparse matrix
       for s = n-1:-1:0
         % Get the row and column indices for the block
         r1 = n*(1:s1) - s;
         r2 = n*(1:s2) - s;
         [colIdx,rowIdx]  = meshgrid(r2,r1);
         % Add the values from the original matrix to the sparse matrix
         rows = [rows; rowIdx(:)];
         cols = [cols; colIdx(:)];
         values = [values; mat(:)];
       end

       % Create the sparse matrix directly from the row, column indices, and values
       mat_new = sparse(rows, cols, values, s1 * n, s2 * n);
     end
   end
end