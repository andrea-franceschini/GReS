classdef Sedimentation < PhysicsSolver
  % Sedimentation model subclass

  properties
     grid = gridForSedimentation()
     mapSediments = SedimentsMap()
     bcsTypes
  end

  properties (Access = protected)
    fieldId
  end

  methods (Access = public)
    function obj = Sedimentation(domain)
      obj@PhysicsSolver(domain);
    end

    function registerSolver(obj,input)
      % setup the solver with custom input

      if ~(isfield(input,'domain') || isfield(input,'Domain'))
        error("Domain for the simulation not defined!");
      end
      obj.grid = gridForSedimentation(input.domain);
      
      % Initialize the sedimentation control
      obj.mapSediments = SedimentsMap(obj.grid.ncells(1:2),input.sediment_map);

      % Initialize the BCs
      obj.prepareBC(input.boundary);


      % Initialize the states


    end

    function assembleSystem(obj,dt)
      obj.Solver.assembleSystem(dt);
    end

    function applyBC(obj,bcId,t)
      obj.Solver.applyBC(bcId,t);
    end

    function applyDirVal(obj,bcId,t)
      bcVar = obj.bcs.getVariable(bcId);
      if ~strcmp(bcVar,obj.getField()) 
        return 
      end
      [bcDofs,bcVals] = getBC(obj,bcId,t);
      if size(bcVals,2)==2
        % skip BC assigned to external surfaces
        return
      end
      state = getState(obj);
      state.data.pressure(bcDofs) = bcVals;
    end

    function updateState(obj,solution)
      obj.Solver.updateState(solution);
    end

    function [cellData,pointData] = writeVTK(obj,fac,t)
      [cellData,pointData] = obj.Solver.writeVTK(fac,t);
    end

    function writeMatFile(obj,fac,tID)
      obj.Solver.writeMatFile(fac,tID);
    end

    function out = isLinear(obj)
      out = obj.Solver.isLinear();
    end
  end

  methods  (Access = private)
    function prepareBC(obj,data)

      if ~isfield(data,"BC")
        error("No boundary was defined for the simulation");
      end
      nbcs=length(data.BC);
      bc = struct('data',[],'cond','SurfBC',...
        'type',[],'variable',[],'surface',[]);
      % obj.bcs = Boundaries();

      % bcs = repmat(, nbcs, 1);
      for i=1:nbcs
        % bcsNames{i} = data.BC(i).name;
        bc.type = data.BC.type;
        bc.variable = data.BC.variable;
        bc.surface = data.BC.surface;
        obj.bcs.db(data.BC(i).name) = bc;
      end

      % bcsNames = [data.BC.name];
      % 
      % 
      % obj.bcs = Boundaries("sedimentation",bcsNames,bcs);
      % 
      % case 3
      %     if strcmp(varargin{1},"sedimentation")
      %       ptNames = varargin{2};
      %       ptSt = varargin{3};
      %       for i=1:length(ptSt)
      %         obj.db(ptNames(i)) = ptSt(i);
      %       end
      %     end

      
    end
  end



  methods (Static)

    function out = getField()
      out = "pressure";
    end

    function out = isSymmetric()
      out = obj.Solver.isSymmetric();
    end

  end

end


