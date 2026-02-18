classdef NonLinearImplicit < SolutionScheme
  % Class for solving non linear implicit problems with changing
  % configuration

  % The time loop is implemented in the base class SolutionScheme

  properties (Access = protected)
    iterNL = 0          % nonlinear iteration number
    iterConfig = 0      % configuration iteration number
    targetVariables     % variables currently solved for
  end



  methods (Access = public)

    function converged = solveStep(obj,varargin)

      % set target variable for each domain and interface
      obj.targetVariables = [varargin{:}];
      if ~isempty(varargin)
        error("solveStep with target variables is not yet implemented")
      end

      hasConfigurationChanged = true;
      absTol = obj.simparams.absTol;
      obj.iterConfig = 0;

      while (hasConfigurationChanged) && (obj.iterConfig < obj.simparams.itMaxConfig)

        gresLog().log(0,'\nConfiguration iteration n. %i \n', obj.iterConfig);
        obj.iterConfig = obj.iterConfig + 1;

        for i = 1:obj.nDom
          obj.domains(i).applyDirVal(obj.t);
        end

        gresLog().log(1,'Iter     ||rhs||     ||rhs||/||rhs_0||\n');

        for i = 1:obj.nDom
          obj.domains(i).assembleSystem(obj.dt);
        end

        for i = 1:obj.nInterf
          obj.interfaces{i}.assembleConstraint();
        end

        for i = 1:obj.nDom
          obj.domains(i).applyBC(obj.t);
        end

        rhs = assembleRhs(obj);
        rhsNorm = norm(cell2mat(rhs),2);
        rhsNormIt0 = rhsNorm;

        tolWeigh = obj.simparams.relTol*rhsNorm;

        gresLog().log(1,'0     %e     %e\n',rhsNorm,rhsNorm/rhsNormIt0);

        newtonConv = false;

        % reset non linear iteration counter
        obj.iterNL = 0;

        %%% NEWTON LOOP %%%
        while (~newtonConv) && (obj.iterNL < obj.simparams.itMaxNR)

          obj.iterNL = obj.iterNL + 1;

          J = assembleJacobian(obj);

          % solve linear system
          du = solve(obj,J,rhs);

          c = 0;

          % update simulation state with linear system solution
          for i = 1:obj.nDom
            if obj.nDom == 1 && obj.nInterf == 0
              sol = du;
            else
              nDof = obj.domains(i).getNumbDoF();
              sol = du(c+1:c+nDof);
              c = c + nDof;
            end
            obj.domains(i).updateState(sol);
          end

          for i = 1:obj.nInterf
            nDof = obj.interfaces{i}.getNumbDoF();
            sol = du(c+1:c+nDof);
            obj.interfaces{i}.updateState(sol);
            c = c + nDof;
          end

          % reassemble system
          for i = 1:obj.nDom
            obj.domains(i).assembleSystem(obj.dt);
          end

          for i = 1:obj.nInterf
            obj.interfaces{i}.assembleConstraint();
          end

          for i = 1:obj.nDom
            obj.domains(i).applyBC(obj.t);
          end

          rhs = assembleRhs(obj);
          rhsNorm = norm(cell2mat(rhs),2);
          gresLog().log(1,'%d     %e     %e\n',obj.iterNL,rhsNorm,rhsNorm/rhsNormIt0);

          % Check for convergence
          newtonConv = (rhsNorm < tolWeigh || rhsNorm < absTol);

        end % end newton loop

        if newtonConv % Newton Convergence

          hasConfigurationChanged = false;

          % update the active set
          for i = 1:obj.nDom
            hasConfigurationChanged = any([hasConfigurationChanged; ...
              obj.domains(i).updateConfiguration()]);
          end

          for i = 1:obj.nInterf
            hasConfigurationChanged = any([hasConfigurationChanged; ...
              obj.interfaces{i}.updateConfiguration()]);
          end

        else

          break

        end

      end

      converged = newtonConv && ~hasConfigurationChanged;

    end
  end




  methods (Access = protected)


    function setLinearSolver(obj)
      % Check if there is manual input from the user, if not use defaults
      start_dir = pwd;
      chronos_xml = fullfile(start_dir,'linsolver.xml');
      if(isfile(chronos_xml))
        obj.linsolver = linearSolver(obj.domains,obj.interfaces,chronos_xml);
      else
        if gresLog().getVerbosity > 2
          fprintf('Using default values for linsolver\n');
        end
        obj.linsolver = linearSolver(obj.domains,obj.interfaces);
      end
    end



    function out = computeRhsNorm(obj)

      %Return maximum norm of the entire domain
      rhsNorm = zeros(obj.nDom,1);
      for i = 1:obj.nDom
        nRhs = length(obj.domains(i).dofm.subList);
        rhsNorm_loc = zeros(nRhs,1);
        for j = 1:nRhs
          rhsNorm_loc(j) = norm(obj.domains(i).rhs{j}, obj.simparams.pNorm);
        end
        rhsNorm(i) = sqrt(sum(rhsNorm_loc.^2));
      end
      out = norm(rhsNorm);
    end



    function J = assembleJacobian(obj)

      J = cell(obj.nVars + obj.nInterf);

      k = 0;

      for iD = 1:obj.nDom

        dom = obj.domains(iD);
        nV = dom.dofm.getNumberOfVariables;

        % inner domain blocks
        J(k+1:k+nV,k+1:k+nV) = getJacobian(dom);

        for iI = 1:numel(dom.interfaceList)

          q = dom.interfaceList(iI);

          % domain coupling blocks
          [J(k+1:k+nV,obj.nVars+q), J(obj.nVars+q,k+1:k+nV)] = ...
            getInterfaceJacobian(dom,iI);

        end

        k = k + nV;

      end

      for iI = 1:obj.nInterf
        interf = obj.interfaces{iI};
        % constraint blocks
        J{obj.nVars+iI,obj.nVars+iI} = getJacobian(interf);
      end

    end


    function rhs = assembleRhs(obj)
      % assemble blocks of rhs for multidomain system

      % each variable field of each domain represents a cell row
      rhs = cell(obj.nVars + obj.nInterf,1);

      k = 0;

      for iD = 1:obj.nDom
        nV = obj.domains(iD).dofm.getNumberOfVariables;
        rhs(k+1:k+nV) = getRhs(obj.domains(iD));
        k = k+nV;
      end

      for iI = 1:obj.nInterf

        rhs{k+1} = getRhs(obj.interfaces{iI});

        % each interface has one only multiplier field!
        k = k+1;
      end
    end


  end

end

