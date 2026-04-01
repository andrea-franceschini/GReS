classdef PreStressMechanics < Poromechanics

  % This model allow computing initial stress distribution to balance
  % gravitational forces and other possible external sources

  % The call stress = PreStressMechanics.getInitialStress() solves a
  % stand-alone one-step problem and return the initial stress on each cell
  % in the model, without going trough a full simulationLoop().

  properties
    specWeight         % stores specific weight in each cell of the model
    gravity = false
  end

  properties (Access = private)
    fieldId = 1
  end

  methods (Access = public)

    function obj = PreStressMechanics(domain)

      % call physicsSolver constructor
      obj@Poromechanics(domain);

    end


    function registerSolver(obj,varargin)

      registerSolver@Poromechanics(obj,varargin)

      setGravity(obj,varargin)

    end



    function assembleSystem(obj,dt)
      % compute the displacements matrices and rhs in the domain
      assembleSystem@Poromechanics(obj,dt);
      obj.domain.rhs{obj.fieldId} = ...
        obj.domain.rhs{obj.fieldId}  + computeGravityRhs(obj);

    end




    function rhsGrav = computeGravityRhs(obj)

      % add gravity contribution to the rhs
      nDoF = obj.domain.dofm.getNumbDoF(obj.fieldId);
      rhsGrav = zeros(nDoF,1);

      for el = 1:obj.mesh.nCells
        nodes = obj.mesh.cells(el,:);
        dof = dofId(nodes,3);
        % only affects z contribution
        dof = dof(3:3:end);
        vtkId = obj.mesh.cellVTKType(el);
        elem = getElement(obj.elements,vtkId);
        N = getBasisFinGPoints(elem);
        dJWeighed = getDerBasisFAndDet(elem,el,3);
        b = obj.specWeight(el)*(N'*dJWeighed');
        rhsGrav(dof) = rhsGrav(dof) + b;
      end

    end





    function initialStress = getInitialStress(obj,varargin)


      % setup a self-contained one-step GReS simulation to compute initial stress

      % default parameters for the simulation
      parm = struct('Start',0.0,'End',1.0,'DtInit',1.0,'DtMax',1.0,'DtMin',0.01);
      default = struct('simulationparameters',SimulationParameters(parm),...
        'domains',obj.domain);

      params = readInput(default,varargin{:});

      solv = NonLinearImplicit(params);

      assert(solv.nDom == 1 && solv.nInterf == 0, ...
        "PreStressMechanics can be used only as a stand-alone PhysicsSolver within a single domain");

      solv.simulationLoop();

      [initialStress,~] = obj.finalizeState(obj.getState());

    end


  end

  methods (Access=private)

    function setGravity(obj,input)

      % compute the specific weight in each cell of the model

      if any(ismissing(input))
        % no gravity in the model
        return
      else
        obj.gravity = true;
        % read input for gravity
        default = struct('hydrostaticPressure',0,...
                         'waterDepth',0.0,...
                         'obg',missing);
        grav = readInput(default,input);
        grav.hydrostaticPressure = logical(grav.hydrostaticPressure);
      end

      obj.specWeight = getSpecificWeight(obj,grav);


    end


    function gamma = getSpecificWeight(obj,grav)

      mat = obj.domain.materials;

      gamma = zeros(obj.mesh.nCells,1);
      poro = zeros(obj.mesh.nCells,1);

      if ~any(ismissing(grav.obg))
        %  specific weight computed from overburden gradient
        try
          obg = str2func(['@(x,y,z)', char(grav.obg)]);
        catch
          error("The overburden gradient obg field must be a valid function handle!")
        end

        C = obj.mesh.cellCentroid;
        gamma = obj(C(:,1),C(:,2),C(:,3));

      else
        % specific weight computed from material user-input

        for i = 1:obj.mesh.nCellTag

          cellsID = obj.mesh.cellTag == i;
          gamma(cellsID) = mat.getSpecificWeight(i);

        end

      end

      % submerge specific weight
      if grav.hydrostaticPressure

        fluid = mat.getFluid();
        if numel(fluid) == 0
          error("hydrostatic contribution to gravitational forces requires definition of fluid and porous rock properties");
        end

        gammaW = fluid.getSpecificWeight();

        for i = 1:obj.mesh.nCellTag

          cellsID = obj.mesh.cellTag == i;
          rock = mat.getPorousRock(i);
          poro(cellsID) = rock.getPorosity();

        end

        maxZ = max(obj.mesh.cellCentroid(:,3));
        depth = abs(obj.mesh.cellCentroid - maxZ);
        sID = depth > grav.waterDepth;

        % compute submerged specific weight (for cells below waterlevel)
        gamma(sID) = (1 - poro(sID)) * ( gamma(sID) - gammaW );
       
      end

    end

  end
end

