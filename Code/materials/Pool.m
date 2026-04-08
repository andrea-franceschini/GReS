classdef Pool < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here

  properties
    data
    mat
    prop
  end

  methods (Access = public)
    function obj = Pool()
      obj.data = cell([]);
      obj.mat = cell([]);
      obj.prop = cell([]);
    end

    function registerProp(obj,matType,propType,data)
    end

    function data = getProp(obj,matType,propType)
      data =[];
    end
  end

  methods (Access = private)
    function addMaterials(obj,input)

      default = struct('Fluid',struct(),...
                       'Solid',struct());
      input = readInput(default,input);

      % order matters: some solid PorousRock properties depend on the fluid
      if numel(fieldnames(input.Fluid)) > 0
        addFluid(obj,input.Fluid);
      end

      if  numel(fieldnames(input.Solid)) > 0
        for i = 1:numel(input.Solid)
          addSolid(obj,input.Solid(i));
        end
      end

    end
  end
end