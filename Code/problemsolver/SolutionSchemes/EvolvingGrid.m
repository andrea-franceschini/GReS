classdef EvolvingGrid < SolutionScheme
  % Class for solving non linear implicit problems with changing
  % configuration

  % The time loop is implemented in the base class SolutionScheme

  properties (SetAccess = protected,GetAccess=public)
    iterNL = 0          % nonlinear iteration number
    targetVariables     % variables currently solved for
    field
    physics
  end

  properties (Access = private)
    printGrow
    printCount = 0
  end

  methods (Access = public)
    function obj = EvolvingGrid(varargin)
      obj@SolutionScheme(varargin{:});
    end

    function converged = solveStep(obj,varargin)
      absTol = obj.simparams.absTol;
      converged = false;

      gresLog().log(1,'Iter     ||rhs||     ||rhs||/||rhs_0||\n');

      obj.physics.applyDirVal(obj.t);
      obj.physics.assembleSystem(obj.dt);
      obj.physics.applyBC(obj.field,obj.t);

      rhs = obj.domains.rhs;
      rhsNorm = norm(cell2mat(rhs),2);
      rhsNormIt0 = rhsNorm;

      tolWeigh = obj.simparams.relTol*rhsNorm;

      gresLog().log(1,'0     %e     %e\n',rhsNorm,rhsNorm/rhsNormIt0);

      % reset non linear iteration counter
      obj.iterNL = 0;

      %%% NEWTON LOOP %%%
      while (~converged) && (obj.iterNL < obj.simparams.itMaxNR)
        obj.iterNL = obj.iterNL + 1;

        J = obj.domains.J;

        % solve linear system
        du = solve(obj,J,rhs);

        % update simulation state with linear system solution
        obj.physics.updateState(du);

        % reassemble system
        obj.physics.assembleSystem(obj.dt);

        % apply boundary condition
        obj.physics.applyBC(obj.t);

        rhs = obj.domains.rhs;
        rhsNorm = norm(cell2mat(rhs),2);
        gresLog().log(1,'%d     %e     %e\n',obj.iterNL,rhsNorm,rhsNorm/rhsNormIt0);

        % Check for convergence
        converged = (rhsNorm < tolWeigh || rhsNorm < absTol);
      end
    end

  end




  methods (Access = protected)

    function setSolutionScheme(obj,varargin)
      % Check that we have an even number of inputs
      if mod(length(varargin), 2) ~= 0
        error('Arguments must come in key-value pairs.');
      end

      % Loop through the key-value pairs
      for k = 1:2:length(varargin)
        key = varargin{k};
        value = varargin{k+1};

        if isempty(value)
          continue
        end

        if ~ischar(key) && ~isstring(key)
          error('Keys must be strings');
        end

        switch lower(key)
          % case 'simulationparameters'
          %   assert(isa(value, 'SimulationParameters')|| isempty(value),msg)
          %   obj.simparams = value;
          case 'simulationparameters'
            obj.simparams = value;
          case 'output'
            obj.output = value;
          case {'domain','domains'}
            obj.domains = value;
          case 'growprint'
            obj.printGrow = value;
          otherwise
            error('Unknown input key %s for SolutionScheme \n', key);
        end
      end

      obj.nDom = numel(obj.domains);
      obj.nInterf = numel(obj.interfaces);

      assert(~isempty(obj.simparams),"Input 'simulationParameters'" + ...
        " is required for SolutionScheme")
      assert(obj.nDom == 1,"Only one domain is admitted when using EvolvingGrid");
      assert(obj.nInterf == 0,"EvolvingGrid do not alound interfaces to another domains");

      obj.nVars = 0;
      obj.domains(1).domainId = 1;
      obj.domains(1).simparams = obj.simparams;
      obj.domains(1).outstate = obj.output;
      if isempty(obj.domains(1).stateOld)
        obj.domains(1).stateOld = copy(obj.domains(1).getState());
      end
      obj.nVars = obj.nVars + obj.domains(1).dofm.getNumberOfVariables();

      obj.attemptedReset = false;

      solv = obj.domains.solverNames;
      obj.physics = obj.domains.getPhysicsSolver(solv);
      obj.field = obj.physics.getField;
    end

    % Sets the linear solver and checks for eventual user input parameters
    function setLinearSolver(obj,varargin)
      if isempty(varargin)
        physname = [];
      else
        % check if the user provided the physics
        physname = varargin{1};
      end
      obj.linsolver = linearSolver(obj,physname);
    end

    function printState(obj)

      if isempty(obj.output)
        return
      end

      print = false;
      newcells = obj.physics.domain.state.data.newcells~=0;
      printByGrow = and(obj.printGrow,newcells);

      if obj.output.timeID <= length(obj.output.timeList)

        outTime = obj.output.timeList(obj.output.timeID);

        % loop over print times contained in the current time step
        while outTime <= obj.t

          assert(outTime >= obj.tOld, 'Print time %f out of range (%f - %f)',...
            outTime, obj.tOld, obj.t);
          assert(obj.t - obj.tOld > eps('double'),...
            'Time step is too small for printing purposes');

          % compute factor to interpolate current and old state variables
          fac = (outTime - obj.tOld)/(obj.t - obj.tOld);

          if isnan(fac) || isinf(fac)
            fac = 1;
          end

          % print vtk
          % create new vtm file
          obj.printVTK(fac,outTime);

          % write results to MAT-file
          obj.printMAT(fac);

          % move to next print time
          obj.output.timeID = obj.output.timeID + 1;

          if obj.output.timeID > length(obj.output.timeList)
            print = true;
            break
          else
            outTime = obj.output.timeList(obj.output.timeID);
          end
        end

        if and(~print,printByGrow)
          fac = 0;  % 0=stateOld, 1=stateNew 

          % write Vtk results
          obj.printVTK(fac,outTime);

          % write results to MAT-file
          obj.printMAT(fac);

          % update the print
          obj.printCount = obj.printCount+1;
        end

      end

    end







  end

  methods (Access = private)

    function printVTK(obj,fac,outTime)
      
      if obj.output.writeVtk
        % set folders
        tID = obj.output.timeID + obj.printCount;
        obj.output.prepareOutputFolders(tID);

        obj.output.vtkFile = com.mathworks.xml.XMLUtils.createDocument('VTKFile');
        toc = obj.output.vtkFile.getDocumentElement;
        toc.setAttribute('type', 'vtkMultiBlockDataSet');
        toc.setAttribute('version', '1.0');
        blocks = obj.output.vtkFile.createElement('vtkMultiBlockDataSet');

        % append blocks looping into domains and interfaces
        for i = 1:obj.nDom
          vtmBlock = obj.domains(i).writeVTK(fac,outTime);
          if ~isempty(vtmBlock)
            blocks.appendChild(vtmBlock);
          end
        end
        %
        for i = 1:obj.nInterf
          vtmBlock = obj.interfaces{i}.writeVTKfile(fac,outTime);
          if ~isempty(vtmBlock)
            blocks.appendChild(vtmBlock);
          end
        end

        toc.appendChild(blocks);
        obj.output.writeVTMFile();
      end
    end

    function printMAT(obj,fac)
      if obj.output.writeSolution
        for i = 1:obj.nDom
          obj.domains(i).writeSolution(fac,obj.printCount+obj.output.timeID);
        end

        for i = 1:obj.nInterf
          obj.interfaces{i}.writeSolution(fac,obj.printCount+obj.output.timeID);
        end
      end
    end
  end


end