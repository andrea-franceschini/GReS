classdef Discretizer < handle
   % General discretizer class
   properties (Access = public)
      solver      % database for physics solvers in the model
      dofm        % dofManager 
      numSolvers  % number of solvers discretized
   end

   methods (Access = public)
      function obj = Discretizer(symmod,simParams,dofManager,grid,mat,varargin)
         %UNTITLED Construct an instance of this class
         %   Detailed explanation goes here
         obj.dofm = dofManager;
         obj.solver = containers.Map('KeyType','double','ValueType','any');
         obj.setDiscretizer(symmod,simParams,dofManager,grid,mat,varargin);
         obj.checkTimeDependence(symmod,mat,simParams);
      end

      function J = assembleJacobian(obj)
         % put together jacobian blocks of SinglePhysicsSolver and
         % CoupledSolver in the model
         % According to the ordering specified in DoF manager class         
      end

      function rhs = assembleRhs(obj)
         % put together rhs blocks of SinglePhysicsSolver and
         % CoupledSolver in the model
      end


      function dSol = solve(obj)
         % assemble and solve linear system
         J = assembleJacobian(obj);
         rhs = assembleRhs(obj);
         dSol = J\-rhs;
      end


      function out = getSolver(obj,fldList)
         % map single field or pair of field to db position
         if isscalar(getFieldList(obj.dofm)) % singlePhysic model
            v = 0;
         else
            nF = numel(getFieldList(obj.dofm)); % multiPhysic model
            cs = cumsum(nF:-1:1);
            v = [0 cs(1:end-1)];
         end
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

      function computeLinearMatrices(obj,stateTmp,statek,dt)
         % loop trough solver database and compute costant jacobian blocks
         for i = 1:obj.numSolvers
            if isLinear(obj.solver(i))
               obj.solver(i).computeMat(stateTmp,statek,dt);
            end
         end
      end

      function computeNLMatricesAndRhs(obj,stateTmp,statek,dt)
         % loop trough solver database and compute non-costant jacobian
         % blocks and rhs block
         for i = 1:obj.numSolvers
            if ~isLinear(obj.solver(i))
               obj.solver(id).computeMat(stateTmp,statek,dt);
            end
            obj.solver(id).computeRhs(stateTmp,statek,dt);
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
         fldList = getFieldList(obj.dofm);
         for i = 1:numel(fldList)
            % loop trough active fields and update the state structure
            setState(obj.getSolver(fldList(i)),state);
         end
      end       
   end

   methods(Access = private)
      function setDiscretizer(obj,symmod,params,dofManager,grid,mat,data)
         fields = getFieldList(obj.dofm); 
         nF = numel(fields);
         % loop over all fields and define corresponding models
         k = 0;
         for i = 1:nF
            for j = i:nF
               k = k+1;
               addSolver(obj,k,fields(i),fields(j),symmod,params,dofManager,grid,mat,data);
            end
         end
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

      function addSolver(obj,id,f1,f2,mod,parm,dof,grid,mat,data)
         % Add new key to solver database
         % Prepare input fields for solver definition
         if ~isCoupled(obj,f1,f2)
            return
         end
         f = sort({char(f1),char(f2)});
         f = join(f,'_');
         switch f{:}
            case 'SPFlow_SPFlow'
               if isSinglePhaseFlow(mod)
                  obj.solver(id) = SPFlow(mod,parm,dof,grid,mat,data);
               elseif isVariabSatFlow(mod)
                  obj.solver(id) = VSFlow(mod,parm,dof,grid,mat,data);
               end
            case 'Poromechanics_Poromechanics'
               obj.solver(id) = Poromechanics(mod,parm,dof,grid,mat,data);
            case 'Poromechanics_SPFlow'
               assert(isSinglePhaseFlow(mod),['Coupling between' ...
                  'poromechanics and unsaturated flow is not yet implemented']);
               obj.solver(id) = Biot(mod,parm,dof,grid,mat,data);
            otherwise
               error('Solver coupling %s with %s does not exist!',f1,f2)
         end
      end

   end
end