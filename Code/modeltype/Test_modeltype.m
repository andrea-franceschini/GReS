close all
clear all
% models = ["SinglePhaseFlow_FVTPFA";
%           "VariabSatFlow_FEM"];
models = "VariabSatFlow_FVTPFA";
m = ModelType(models);
out = isPoromechanics(m)
out = isSinglePhaseFlow(m)
out = isVariabSatFlow(m)
out = isFEMBased(m,'Poro')
out = isFEMBased(m,'Flow')
out = isFVTPFABased(m,'Poro')
out = isFVTPFABased(m,'Flow')
% tic
% m1 = ModelType.Poromechanics;
% m2 = ModelType.SinglePhaseFlow;
% m3 = bitor(ModelType.SinglePhaseFlow,ModelType.Poromechanics);
% for i=1:5000
%   t1 = ModelType.isPoromechanics(m1);
%   t2 = ModelType.isSinglePhaseFlow(m1);
%   t3 = ModelType.isPoromechanics(m3);
%   t4 = ModelType.isSinglePhaseFlow(m3);
%   t5 = ModelType.isCoupFlowPoro(m3);
% end
% totTime1 = toc

% tic
% m1a = ModelType('Poromechanics');
% m2a = ModelType('SinglePhaseFlow');
% m3a = ModelType('CoupFlowPoro');
% for i=1:5000
%   t1 = isPoromechanics(m1a);
%   t2 = isSinglePhaseFlow(m1a);
%   t3 = isPoromechanics(m3a);
%   t4 = isSinglePhaseFlow(m3a);
%   t5 = isCoupFlowPoro(m3a);
% end
% totTime2 = toc

% switch m1
%   case ModelType.SinglePhaseFlow
%     a = 0;
%   case ModelType.Poromechanics
%     a = 10;
% end
% if bitand(m3,ModelType.Poromechanics)
%   b = 1
% end
% if ModelType.isPoromechanics(m4)
%   c = 1
% end