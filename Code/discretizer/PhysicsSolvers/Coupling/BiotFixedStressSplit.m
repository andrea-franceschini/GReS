classdef BiotFixedStressSplit < BiotFullyCoupled
    % Biot model subclass of Coupled Poromechanics with SinglePhaseFlow 
    % This solver is meant for application of the fixed stress split scheme
    % Reference: Castelletto et al. (2015). Accuracy and convergence
    % properties of the fixed‐stress iterative solution of two‐way coupled
    % poromechanics solved in time with backward implicit euler

    properties
        conv           % struct with converged variables at last fixed stress iteration
    end


    methods (Access = public)
      function obj = BiotFixedStressSplit(domain)

        obj@BiotFullyCoupled(domain);

      end

      function registerSolver(obj,input)

        registerSolver@BiotFullyCoupled(obj,input);

        default = struct('BulkModulus',"3D");
        params = readInput(default,input);
        blkMod = char(params.BulkModulus);

        obj.KdrType = str2double(blkMod(1));

      end


      function assembleSystem(obj,dt,var)

        % compute coupling matrix obj.Q
        computeMat(obj);

        % get Jacobian and rhs from single physics solvers
        if strcmpi(var,"pressure")
          obj.flowSolver.assembleSystem(dt);
          computeRelaxationMatrix(obj);
          rhsFlow = computeRhsFlow(obj,dt);
          % add the relaxation matrix 
          obj.domain.J{obj.fldFlow,obj.fldFlow} = ...
            obj.domain.J{obj.fldFlow,obj.fldFlow} + obj.R/dt;
          obj.domain.rhs{obj.fldFlow} = obj.domain.rhs{obj.fldFlow} + rhsFlow;
        elseif strcmpi(var,"displacements")
          obj.mechSolver.assembleSystem(dt);
          rhsMech = computeRhsMech(obj);
          obj.domain.rhs{obj.fldMech} = obj.domain.rhs{obj.fldMech} + rhsMech;
        else
          error("Unknown variable for BiotFixedStressSplit solver." + ...
            "\nAvailable variables are: '%s' and '%s'",...
            obj.mechSolver.getField(),obj.flowSolver.getField())
        end

      end

      function computeMat(obj)

        % call method according to the discretization technique chosen
        if isempty(obj.Q) || ~isLinear(obj)
          computeMatBiot(obj)
        end
      end

      function rhsMech = computeRhsMech(obj,dt)

        pCurr = getState(obj,"pressure");
        p0 = getStateInit(obj,"pressure");

        entsFlow = obj.domain.dofm.getActiveEntities(obj.fldFlow,1);

        % remember the minus sign (look at biot)
        rhsMech = - obj.Q * (pCurr(entsFlow)-p0(entsFlow));

      end

      function rhsFlow = computeRhsFlow(obj,dt)

        % retrieve State variables
          pCurr = getState(obj,"pressure");
          uOld = getStateOld(obj,"displacements");
          pConv = obj.conv.pressure;
          uConv = obj.conv.displacements;

          % select active coefficients of solution vectors
          entsPoro = obj.domain.dofm.getActiveEntities(obj.fldMech,1);
          entsFlow = obj.domain.dofm.getActiveEntities(obj.fldFlow,1);


          rhsFlow = (obj.Q'/dt) * (uConv(entsPoro) - uOld(entsPoro));
          % new contribution of fixed stress split algorithm
          rhsFlow = rhsFlow  + (obj.R/dt) * (pCurr(entsFlow) - pConv(entsFlow));
        end

    

        function updateState(obj,solution,var)
          v = getState(obj,var);
          ents = obj.domain.dofm.getActiveEntities(var,1);
          v(ents) = v(ents) + solution;
          setState(obj,v,var);
        end


        function advanceState(obj,varargin)
          if nargin == 1
            % advance state after fixed stress algorithm converged
            advanceState@BiotFullyCoupled(obj)
          else
            % advance individual physics when newton convergence is
            % achiedeved
            varName = varargin{1};
            % save converged variable
            obj.conv.(varName) = getState(obj,varName);
          end
        end

    end
    
    % methods (Static)
    % 
    %     function out = getField()
    %       out = [Poromechanics.getField(), SinglePhaseFlow.getField()];
    %     end
    % 
    %     function out = isSymmetric()
    %       out = false;
    %     end
    % 
    % end

end


