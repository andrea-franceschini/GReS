classdef (Abstract) Constraint < handle
  %CONSTRAINT Summary of this class goes here
  %   Detailed explanation goes here

  properties

    interface
    coupledVariables

    quadrature                % quadrature owned by this constraint



  end


  methods (Abstract)

    registerConstraint(obj,varargin)

    assembleConstraint(obj,varargin)

    updateState(obj,solution)

    % update the output structures for printing purposes
    [cellData,pointData] = writeVTK(obj,t);

    % write history to MAT-file
    writeSolution(obj,t,tID);

  end



  methods
    function obj = Constraint(interface)
      if nargin > 0
        obj.interface = interface;
      end
    end


    function quadrature = setMortarQuadrature(obj,multiplierType,varargin)

      default = struct('type',"SegmentBasedQuadrature");
      options = readInput(default,varargin{:});

      quadType = options.type;
      grid = InterfaceGrid(obj.interface.grid);
      quadrature = feval(quadType,multiplierType,grid,options);
      quadrature.processMortarTriplets();
      obj.quadrature = quadrature;

    end

  end
end

