classdef ModelType < uint32
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  enumeration
    Poromechanics       (1)
    SinglePhaseFlow     (2)
%     SteadyState         (4)
%     Transient           (8)
  end
  
  methods (Static)
    function res = isPoromechanics(obj1)
%       res = boolean(bitand(obj1,ModelType.Poromechanics));
%       res = bitand(obj1,ModelType.Poromechanics);
%       res = bitand(obj1,1);
      res = ismember(obj1,[1 3]);
    end
    %
    function res = isSinglePhaseFlow(obj1)
%       res = boolean(bitand(obj1,ModelType.SinglePhaseFlow));
%       res = bitand(obj1,ModelType.SinglePhaseFlow);
%         res = bitand(obj1,2);
        res = ismember(obj1,[2 3]);
    end
    %
    function res = isCoupFlowPoro(obj1)
%       resFlow = boolean(bitand(obj1,ModelType.SinglePhaseFlow));
%       resPoro = boolean(bitand(obj1,ModelType.Poromechanics));
%       resFlow = bitand(obj1,ModelType.SinglePhaseFlow);
%       resPoro = bitand(obj1,ModelType.Poromechanics);
%       res = bitand(obj1,3);
%       resPoro = bitand(obj1,ModelType.Poromechanics);
%       res = all([resFlow,resPoro]);
      res = obj1 == 3;
    end
    %
%     function res = isSteadyState(obj1)
% %       res = boolean(bitand(obj1,ModelType.SteadyState));
%       res = bitand(obj1,ModelType.SteadyState);
%     end
%     %
%     function res = isTransient(obj1)
% %       res = boolean(bitand(obj1,ModelType.Transient));
%       res = bitand(obj1,ModelType.Transient);
%     end
  end
end