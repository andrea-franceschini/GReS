classdef ModelType < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = private)
    isPoro = false
    isSPFlow = false
  end
  
  methods (Access = public)
    function obj = ModelType(str)
      if strcmp(str,'Poromechanics')
        obj.isPoro = true;
      elseif strcmp(str,'SinglePhaseFlow')
        obj.isSPFlow = true;
      elseif strcmp(str,'CoupFlowPoro')
        obj.isPoro = true;
        obj.isSPFlow = true;
      end
    end
    
    function out = isPoromechanics(obj)
      out = obj.isPoro;
    end
    
    function out = isSinglePhaseFlow(obj)
      out = obj.isSPFlow;
    end
    
    function out = isCoupFlowPoro(obj)
      out = obj.isPoro && obj.isSPFlow;
    end
  end
end