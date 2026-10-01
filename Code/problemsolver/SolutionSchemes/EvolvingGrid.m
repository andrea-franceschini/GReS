classdef EvolvingGrid < SolutionScheme
  % Class for solving non linear implicit problems with changing
  % configuration

  % The time loop is implemented in the base class SolutionScheme

  properties (SetAccess = protected,GetAccess=public)
    iterNL = 0          % nonlinear iteration number
    targetVariables     % variables currently solved for
    field
    physics
    isBackStep = false
    freezeCount = 0
    freezeAt
  end

  properties (Access = private)
    printGrow
    printInterv
    printCount = 0
    printList = []
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
        % converged = rhsNorm < absTol;
      end
    end

  end




  methods (Access = protected)
    function setSolutionScheme(obj,varargin)
      default = struct('simulationparameters',SimulationParameters.empty,...
                       'output',missing,...
                       'domains',Discretizer.empty,...
                       'growprint',0,...
                       'intervalprint',missing,...
                       'freezeAt',missing);
      params = readInput(default,varargin{:});

      obj.simparams = params.simulationparameters;
      obj.domains = params.domains;
      if ~ismissing(params.output)
        obj.output = params.output;
      end

      obj.printGrow = logical(params.growprint);
      obj.printInterv = params.intervalprint;
      if ismissing(obj.printInterv)
        % Define an interval outside of the simulation
        tf = params.simulationparameters.tMax+params.simulationparameters.dtMax;
        obj.printInterv = [tf, 2*tf];
      end

      if ismissing(params.freezeAt)
        obj.freezeAt = params.simulationparameters.dtMin + ...
          (params.simulationparameters.dtMax-params.simulationparameters.dtMin)/10;
      else
        obj.freezeAt = params.freezeAt;
      end

      obj.nDom = numel(obj.domains);
      obj.nInterf = numel(obj.interfaces);      
      obj.nVars = 1;

      assert(~isempty(obj.simparams),"Input 'simulationParameters'" + ...
        " is required for SolutionScheme")
      assert(obj.nDom == 1,"Only one domain is admitted when using EvolvingGrid");
    end
    
    function setLinearSolver(obj,varargin)
      % Sets the linear solver and checks for eventual user input parameters
      if isempty(varargin)
        physname = [];
      else
        % check if the user provided the physics
        physname = varargin{1};
      end

      solv = obj.domains.solverNames;
      obj.physics = obj.domains.getPhysicsSolver(solv);
      obj.field = obj.physics.getField;

      obj.linsolver = linearSolver(obj,physname);
    end

    function manageNextTimeStep(obj,flConv)
      if ~flConv
        % BACKSTEP
        % newton did not converge or configuration changed too many times
        obj.t = obj.tOld;
        obj.tStep = obj.tStep - 1;

        dt = obj.dt/obj.simparams.divFac;  % Time increment chop
        if (dt<obj.freezeAt) && ~obj.isBackStep
          dt = obj.dt;  % Time increment chop
          obj.physics.cbFreeze = true;
          obj.freezeCount = obj.freezeCount +1;
        end
        % obj.dt = obj.dt/obj.simparams.divFac;  % Time increment chop

        obj.dt = dt;  % Time increment chop
        obj.dtSave = obj.dt;
        goBackState(obj);
        obj.t = obj.t + obj.dt;
        if obj.dt < obj.simparams.dtMin
          error('Minimum time step reached')
        else
          gresLog().log(0,'\n %s \n','BACKSTEP')
        end
        obj.isBackStep = true;
        % return
      else
        % TIME STEP CONVERGED - advance to the next time step
        printState(obj);
        advanceState(obj);

        % go to next time step
        obj.dt = getNextDt(obj);
        obj.tOld = obj.t;
        obj.t = obj.t + obj.dt;

        % allow new survival attempts on new time steps
        if obj.simparams.attemptSimplestConfiguration
          obj.attemptedReset = false;
        end
        obj.isBackStep = false;
      end
    end

    function printState(obj)
      if isempty(obj.output)
        return
      end

      tf = obj.t;
      t0 = tf-obj.dt;
      flagPrintExtra = false;
      if and(tf>=obj.printInterv(1),t0<=obj.printInterv(2))
        flagPrintExtra = true;        
      end

      if getState(obj.physics.domain,'newcells')~=0
        flagPrintExtra = or(flagPrintExtra,obj.printGrow);
      end

      if flagPrintExtra
        % if printByGrow
        fac = 0;  % 0=stateOld, 1=stateNew - Because the mesh update
        % happens after the print, to plot after the grow, i plot the
        % oldstate in the next time step.

        pos = obj.output.timeID + obj.printCount;

        % write Vtk results
        obj.printVTK(fac,obj.tOld,pos);

        % write results to MAT-file
        obj.printMAT(fac,pos);

        % update the print
        obj.printList(pos)=obj.tOld;
        obj.printCount = obj.printCount+1;
      end

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
          pos = obj.output.timeID + obj.printCount;
          obj.printVTK(fac,outTime,pos);
          obj.printList(pos)=outTime;
          print = true;

          % write results to MAT-file
          obj.printMAT(fac,pos);

          % move to next print time
          obj.output.timeID = obj.output.timeID + 1;

          if obj.output.timeID > length(obj.output.timeList)
            break
          else
            outTime = obj.output.timeList(obj.output.timeID);
          end
        end
      end
    end

    function printVTK(obj,fac,outTime,tID)
      if obj.output.writeVtk
        % set folders
        obj.output.prepareOutputFolders(tID);

        obj.output.vtkFile = com.mathworks.xml.XMLUtils.createDocument('VTKFile');
        toc = obj.output.vtkFile.getDocumentElement;
        toc.setAttribute('type', 'vtkMultiBlockDataSet');
        toc.setAttribute('version', '1.0');
        blocks = obj.output.vtkFile.createElement('vtkMultiBlockDataSet');

        vtmBlock = obj.output.vtkFile.createElement('Block');
        [cellData,pointData] = obj.physics.writeVTK(fac,outTime);

        cellData = OutState.printMeshData(obj.physics.grid,cellData);

        % write dataset to vtmBlock
        obj.output.writeVTKfile(vtmBlock,'Sedimentation',obj.physics.grid,...
          outTime, pointData, cellData, [], [],tID);

        blocks.appendChild(vtmBlock);

        toc.appendChild(blocks);
        obj.output.writeVTMFile(tID);
      end
    end

    function printMAT(obj,fac,timeID)
      if obj.output.writeSolution
        obj.physics.writeSolution(fac,timeID);
      end
    end

    function finalize(obj)
      if ~isempty(obj.output)
        obj.output.savePvd(obj.printList);

        obj.output.saveHistory;
      end
    end

  end



end