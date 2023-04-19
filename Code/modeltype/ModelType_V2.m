classdef ModelType < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = private)
    isPoro = false
    isSPFlow = false
    isFEM = false;
    isFVTPFA = false;
  end
  
  methods (Access = public)
    function obj = ModelType(str)
      [str1,str2] = strtok(str,'_');
      if strcmpi(str1,'Poromechanics')
        obj.isPoro = true;
      elseif strcmpi(str1,'SinglePhaseFlow')
        obj.isSPFlow = true;
        if strcmpi(str2(2:end),'FEM')
          obj.isFEM = true;
        elseif strcmpi(str2(2:end),'FVTPFA')
          obj.isFVTPFA = true;
        else
          error('Discretization scheme for the single phase flow not defined or invalid');
        end
      elseif strcmpi(str,'CoupFlowPoro')
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
    
    function out = isFEMBased(obj)
      out = obj.isFEM;
    end
    
    function out = isFVTPFABased(obj)
      out = obj.isFVTPFA;
    end
  end
end