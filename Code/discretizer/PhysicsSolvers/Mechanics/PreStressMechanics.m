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

      setGravity(obj,varargin{:})

    end



    function assembleSystem(obj,dt)
      % compute the displacements matrices and rhs in the domain
      assembleSystem@Poromechanics(obj,dt);

      if ~isempty(obj.specWeight)
        % gravitational forces
        obj.domain.rhs{obj.fieldId} = ...
          obj.domain.rhs{obj.fieldId}  + computeGravityRhs(obj);
      end

    end




    function rhsGrav = computeGravityRhs(obj)


      % add gravity contribution to the rhs
      nDoF = obj.domain.dofm.getNumbDoF(obj.fieldId);
      rhsGrav = zeros(nDoF,1);
      coordinates = obj.grid.coordinates;
      cells = obj.grid.cells;
      dofm = obj.domain.dofm;
      subCells = dofm.getFieldCells(obj.fieldId);

      for vtkId = cells.vtkTypes

        subCellsLoc = obj.grid.getCellsByVTKId(vtkId,subCells);
        elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);
        N = elem.getBasisFinGPoints();

        % get node topology for given vtk type
        topol = obj.grid.getCellNodes(subCellsLoc);

        for i = 1:numel(subCellsLoc)

          % assembly loop for homogeneous element type
          el = subCellsLoc(i);
          nodes = topol(i,:);
          dof = dofm.getLocalDoF(obj.fieldId,nodes);
          coord = coordinates(nodes,:);
          dof = dof(3:3:end);
          [~,dJWeighed] = getDerBasisFAndDet(elem,coord);
          b = obj.specWeight(el)*(N'*dJWeighed');
          rhsGrav(dof) = rhsGrav(dof) + b;

        end

      end
    end





    function initialStress = getInitialStress(obj,varargin)
      %
      % getInitialStress  Compute initial stress via a one-step GReS simulation.
      % DEFAULT SETUP
      % initialStress = getInitialStress(obj)   % default setup
      % SETUP WITH USER DEFINED PARAMETERS
      % initialStress = getInitialStress(obj, 'simulationparameters', sp, 'output', out)
      %
      % Accept only a single domain with no interfaces.
      %
      % Output:
      %   initialStress - [nGP x 6] matrix with stress tensor components

      % default parameters for the simulation
      parm = struct('Start',0.0,'End',1.0,'DtInit',1.0,'DtMax',1.0,'DtMin',0.01);
      default = struct('simulationparameters',SimulationParameters(parm),...
        'domains',obj.domain);

      params = readInput(default,varargin{:});

      solv = NonLinearImplicit(params);

      assert(solv.nDom == 1 && solv.nInterf == 0, ...
        "getInitialStress can be used only as a stand-alone PhysicsSolver within a single domain");

      solv.simulationLoop();

      initialStress = obj.domain.state.data.stress;

    end


  end

  methods (Access=private)

    function setGravity(obj,varargin)

      % compute the specific weight in each cell of the model

      if any(ismissing(varargin))
        % no gravity in the model
        return
      else
        obj.gravity = true;
        % read input for gravity
        default = struct('hydrostaticPressure',0,...
                         'waterDepth',0.0,...
                         'obg',missing);
        grav = readInput(default,varargin{:});
        grav.hydrostaticPressure = logical(grav.hydrostaticPressure);
      end

      obj.specWeight = getSpecificWeight(obj,grav);


    end


    function gamma = getSpecificWeight(obj,grav)

      tol = 1e-6;

      mat = obj.domain.materials;
      cells = obj.domain.grid.cells;

      gamma = zeros(cells.num,1);
      poro = zeros(cells.num,1);

      if ~any(ismissing(grav.obg))
        %  specific weight computed from overburden gradient
        try
          obg = str2func(['@(x,y,z)', char(grav.obg)]);
        catch
          error("The overburden gradient obg field must be a valid function handle!")
        end

        C = cells.center;
        gamma = obg(C(:,1),C(:,2),C(:,3));

      else
        % specific weight computed from material user-input

        for i = 1:cells.nTag

          cellsID = cells.tag == i;
          gamma(cellsID) = mat.getSpecificWeight(i);

        end

      end

      % submerge specific weight
      if grav.hydrostaticPressure

        fluid = mat.getFluid();
        if numel(fluid) == 0
          error("Hydrostatic contribution to gravitational forces requires definition of fluid and porous rock properties.");
        end

        gammaW = fluid.getSpecificWeight();

        for i = 1:cells.nTag
          cellsID = cells.tag == i;
          rock = mat.getPorousRock(i);
          poro(cellsID) = rock.getPorosity();

        end

        maxZ = max(cells.center(:,3));
        depth = abs(cells.center(:,3) - maxZ);
        sID = depth > grav.waterDepth - tol;

        % compute submerged specific weight (for cells below waterlevel)
        gamma(sID) = (1 - poro(sID)) .* ( gamma(sID) - gammaW );
       
      end

    end

  end
end

