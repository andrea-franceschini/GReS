classdef Discretizer < handle
   % General discretizer class
   properties (GetAccess=public, SetAccess=private)
      solver      % database for physics solvers in the model
      dofm        % dofManager 
      numSolvers  % number of solvers discretized
      mod
      fields
   end

   methods (Access = public)
      function obj = Discretizer(symmod,simParams,dofManager,grid,mat,varargin)
         %UNTITLED Construct an instance of this class
         %   Detailed explanation goes here
         obj.mod = symmod;
         obj.dofm = dofManager;
         obj.solver = containers.Map('KeyType','double','ValueType','any');
         obj.setDiscretizer(symmod,simParams,dofManager,grid,mat,varargin);
         obj.checkTimeDependence(symmod,mat,simParams);
      end
      
      function applyBC(obj,bound,t,state)
         % Apply boundary condition to blocks of physical solver
         % ents: id of constrained entity
         % vals: value to assign to each entity
         bcList = bound.db.keys;
         % get entities and values of boundary condition
         for bc = string(bcList)
            field = translatePhysic(bound.getPhysics(bc),obj.mod);
            % get id of constrained entities and corresponding BC values
            [bcEnts,bcVals] = getBC(getSolver(obj,field),bound,bc,t,state);
            % apply Boundary conditions to each Jacobian/rhs block
            for f = obj.fields
               if ~isCoupled(obj,field,f)
                  continue
                  % skip pair of uncoupled physics
               end
               switch bound.getType(bc)
                  case {'Dir','Spg'}
                     applyDirBC(obj.getSolver({field,f}),field,bcEnts,bcVals);
                  case {'Neu','VolumeForce'}
                     applyNeuBC(obj.getSolver({field,f}),bcEnts,bcVals);
               end
            end
         end
      end

      function state = applyDirVal(obj,bound,t,state)
         % Apply boundary condition to blocks of physical solver
         % ents: id of constrained entity
         % vals: value to assign to each entity
         bcList = bound.db.keys;
         % get entities and values of boundary condition
         for bc = string(bcList)
            if ~strcmp(bound.getType(bc),'Dir')
               continue
            end
            field = translatePhysic(bound.getPhysics(bc),obj.mod);
            state = getSolver(obj,field).applyDirVal(bound,bc,t,state);
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
         % convert cell to matrix
         J = cell2mat(J);
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
         rhs = cell2mat(rhs);
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


      function stateTmp = computeMatricesAndRhs(obj,stateTmp,statek,dt)
         % loop trough solver database and compute non-costant jacobian
         % blocks and rhs block
         for i = 1:obj.numSolvers
            stateTmp = computeMat(obj.solver(i),stateTmp,statek,dt);
            stateTmp = computeRhs(obj.solver(i),stateTmp,statek,dt);
         end
      end

      function out = isCoupled(obj,field1,field2)
         % check if input fields are coupled, i.e exist cells having both
         % fields activated
         sub1 = getActiveSubdomain(obj.dofm,field1);
         sub2 = getActiveSubdomain(obj.dofm,field2);
         out = any(intersect(sub1,sub2));
      end

      function state = setState(obj)
         % loop trough active single physics solver and update the state class
         % accordingly
         state = struct();
         state.t = 0;
         for i = 1:numel(obj.fields)
            % loop trough active fields and update the state structure
            state = setState(obj.getSolver(obj.fields(i)),state);
         end
      end

      function state = updateState(obj,state,du)
         % update current state
         for i = 1:numel(obj.fields)
            state = obj.getSolver(obj.fields(i)).updateState(state,du);
         end
      end
   end

   methods(Access = private)
      function setDiscretizer(obj,symmod,params,dofManager,grid,mat,data)
         flds = getFieldList(obj.dofm); 
         nF = numel(flds);
         % loop over all fields and define corresponding models
         k = 0;
         for i = 1:nF
            for j = i:nF
               k = k+1;
               addPhysics(obj,k,flds(i),flds(j),symmod,params,dofManager,grid,mat,data);
            end
         end
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

      function addPhysics(obj,id,f1,f2,mod,parm,dof,grid,mat,data)
         % Add new key to solver database
         % Prepare input fields for solver definition
         if ~isCoupled(obj,f1,f2)
            return
         end
         f = sort({char(f1),char(f2)});
         f = join(f,'_');
         switch f{:}
            case 'SPFlow_SPFlow'
               obj.solver(id) = SPFlow(mod,parm,dof,grid,mat,data);
            case 'Poromechanics_Poromechanics'
               obj.solver(id) = Poromechanics(mod,parm,dof,grid,mat,data);
            case 'Poromechanics_SPFlow'
               assert(isSinglePhaseFlow(mod),['Coupling between' ...
                  'poromechanics and unsaturated flow is not yet implemented']);
               obj.solver(id) = Biot(mod,parm,dof,grid,mat,data);
            case 'VSFlow_VSFlow'
               obj.solver(id) = VSFlow(mod,parm,dof,grid,mat,data);
            otherwise
               error('A physical module coupling %s with %s is not available!',f1,f2)
         end
      end

   end
end