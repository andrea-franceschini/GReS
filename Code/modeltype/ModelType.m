classdef ModelType < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = private)
    ModSettings = zeros(3,1); % [Poromechanics; Flow; Poisson]
    % ModSettings(1) -> Poromechanics
    %    0   -> The model is inactive
    %  10-19 -> Continuous mechanics
    %
    % ModSettings(2) -> Flow
    %    0   -> The model is inactive
    %  10-19 -> Single-phase flow
    %  20-29 -> Variably saturated flow (Richards eq.)

    % ModSettings(3) -> Poisson
    %    0   -> The model is inactive
    %  10-19 -> Poisson problem with non homogeneous bcs
    %   ...
    %
    % Discretization schemes:
    %    0   -> FEM      (Finite Element Method)
    %    1   -> FV-TPFA  (Finite Volume with Two-Point Flux Approximation)
  end
  
  methods (Access = public)
    function obj = ModelType(str)
      setModelType(obj,str);
    end
    
    function out = isPoromechanics(obj)
      out = false;
      if obj.ModSettings(1) > 0
        out = true;
      end
    end
    
    function out = isFlow(obj)
      out = false;
      if obj.ModSettings(2) > 0
        out = true;
      end
    end

    function out = isPoisson(obj)
      out = false;
      if obj.ModSettings(3) > 0
        out = true;
      end
    end
    
    function out = isSinglePhaseFlow(obj)
      out = false;
      if floor(obj.ModSettings(2)/10) == 1
        out = true;
      end
    end
    
    function out = isVariabSatFlow(obj)
      out = false;
      if floor(obj.ModSettings(2)/10) == 2
        out = true;
      end
    end
    
%     function out = isCoupFlowPoro(obj)
%       out = false;
%       if all(obj.ModSettings(1:2))
%         out = true;
%       end
%     end
    
    function out = isFEMBased(obj,modStr)
      out = false;
      r = ModelType.findIDPhysics(modStr);
      if obj.ModSettings(r) > 0 && rem(obj.ModSettings(r),10) == 0
        out = true;
      end
    end
    
    function out = isFVTPFABased(obj,modStr)
      out = false;
      r = ModelType.findIDPhysics(modStr);
      if obj.ModSettings(r) > 0 && rem(obj.ModSettings(r),10) == 1
        out = true;
      end
    end

    function physics = getAvailPhysics(obj)
       % using cellfun to vectorize the checks on active physical modules
       physicsList = {
          "Poromechanics", @isPoromechanics
          "SinglePhaseFlow", @isSinglePhaseFlow;
          "VariablySaturatedFlow", @isVariabSatFlow;
          "Poisson", @isPoisson;
          };
       physics = physicsList(cellfun(@(f) f(obj), physicsList(:, 2)), 1);
    end

  function checkAvailPhysics(obj,physic)
    % check if input physic is available in the model or correctly typed
    phList = string(getAvailPhysics(obj));
    check = ismember(physic,phList);
    if check
      return
    else
      error(['Input physic %s is not available in the model.\n' ...
        'Available physics are:' ...
        ' %s'],physic,strjoin(phList,','));
    end
  end
end
  
  methods (Access = private)
    function setModelType(obj,str)
      for i=1:length(str)
        spltStr = split(str(i),'_');
        %
        switch spltStr(1)
          case 'Poromechanics'
            r = 1;
            v = 10;
          case 'SinglePhaseFlow'
            r = 2;
            v = 10;
          case 'VariabSatFlow'
            r = 2;
            v = 20;
          case 'Poisson'
            r = 3;
            v = 10;
          otherwise
            error(['%s model is invalid\n', ...
              'Accepted physics are: Poromechanics, SinglePhaseFlow,\n', ...
              'VariabSatFlow, Poisson'],spltStr(1));
        end
        if obj.ModSettings(r) ~= 0
          s = ModelType.findPhysicsFromID(r);
          error('Multiple %s models have been declared',s);
        else
          obj.ModSettings(r) = v;
        end
        %
        switch spltStr(2)
          case 'FEM'
            % Nothing to do
          case 'FVTPFA'
            obj.ModSettings(r) = obj.ModSettings(r) + 1;
          otherwise
            error(['%s spatial discretization scheme is invalid\n', ...
              'Accepted schemes are:\n', ...
              '  FEM    -> Finite Element Model\n', ...
              '  FVTPFA -> Finite Volume Method with Two-Point Flux Approximation'],spltStr(2));
        end
        %
        checkCompatibleSchemes(obj);
      end
    end
    
    function checkCompatibleSchemes(obj)
      if rem(obj.ModSettings(1),10) == 1
        error('FVTPFA is not implemented for Poromechanics')
      end
    end
  end
  
  methods (Static = true)
    function r = findIDPhysics(str)
      switch str
        case 'Poromechanics'
          r = 1;
        case {'Flow','SinglePhaseFlow','VariablySaturatedFlow'}
          r = 2;
        case 'Poisson'
          r = 3;
        otherwise
          error(['%s model is invalid\n', ...
            'Accepted keywords are: Poro, Flow and Thermal,\n'],str);
      end
    end
    
    function s = findPhysicsFromID(id)
      switch id
        case 1
          s = 'Poromechanics';
        case 2
          s = 'Flow';
        case 3
          s = 'Poisson';
      end
    end
  end

end