classdef TabularCurve < handle

  properties (Access = private)
    tabW
    derivW
    derivW2
    nPoints
  end

  methods (Access = public)
    % Class constructor method
    function obj = TabularCurve(fID, matFileName)
      % Calling the function to set the object properties 
      obj.readTabularCurve(fID, matFileName);
    end
    
    function varargout = interpTable(obj,x)
      varargout{1} = zeros(length(x),1);
      %
      x = -x;
      xNeg = x <= 0;
      xOutL = x > 0 & x < obj.tabW(1,1);
      xIn = x >= obj.tabW(1,1) & x <= obj.tabW(end,1);
      xOutR = x > obj.tabW(end,1);
      %
      varargout{1}(xNeg) = 1;
      %
      dyL = (obj.tabW(1,2) - 1)/obj.tabW(1,1);
      varargout{1}(xOutL) = 1 + dyL*x(xOutL);
      %
      varargout{1}(xIn) = interp1(obj.tabW(:,1),obj.tabW(:,2),x(xIn));
      %
      dyR = (obj.tabW(end,2) - obj.tabW(end-1,2))/(obj.tabW(end,1) - obj.tabW(end-1,1));
      varargout{1}(xOutR) = obj.tabW(end-1,2) + dyR*(x(xOutR) - obj.tabW(end-1,1));
      idOutR = find(xOutR);
      yNeg = varargout{1}(xOutR) < 0;
      varargout{1}(idOutR(yNeg)) = 0;
      if nargout > 1 % compute first derivative
        varargout{2} = zeros(length(x),1);
        %
        varargout{2}(xNeg) = 0;
        %
        varargout{2}(xOutL) = dyL;
        %
        [~,~,binID] = histcounts(x(xIn),obj.tabW(:,1));
        varargout{2}(xIn) = obj.derivW(binID);
        %
        varargout{2}(idOutR) = dyR;
        varargout{2}(idOutR(yNeg)) = 0;
      end
      if nargout > 2 %compute also the second derivative
          varargout{3} = zeros(length(x),1);
          %
          varargout{3}(xNeg) = 0;
          %
          varargout{3}(xOutL) = 0;
          %
          % compute histogram columns for second derivative
          nPts = length(obj.tabW);
          tabWinterp = (interp1(linspace(0,1,nPts),obj.tabW(:,1),linspace(0,1,nPts-1)))';
          [~,~,binID] = histcounts(x(xIn),tabWinterp);
          varargout{3}(xIn) = obj.derivW2(binID);
          %
          varargout{3}(idOutR) = 0;
          varargout{3}(idOutR(yNeg)) = 0;
      end
    end
%     function [y, dy] = interpTable(obj, x)
%       x = -x;
%       if x < 0
%         y = 1;
%         dy = 0;
%       elseif (x >= obj.tabW(1,1)) & (x <= obj.tabW(end,1))  % Enter and table look-up
%         p1 = 1;
%         p2 = obj.nPoints;
%         while (p2 - p1) > 1
%           pos = floor((p2 + p1)/2);
%           if obj.tabW(pos,1) >= x
%             p2 = pos;
%           else
%             p1 = pos;
%           end
%         end
%         dy = (obj.tabW(p2,2) - obj.tabW(p1,2))/(obj.tabW(p2,1) - obj.tabW(p1,1));
%         y = obj.tabW(p1,2) + dy*(x - obj.tabW(p1,1));
%       elseif x < obj.tabW(1,1) % Out of range to the left
%         dy = (obj.tabW(1,2) - 1)/obj.tabW(1,1);
%         y = 1 + dy*x;
%       elseif x > obj.tabW(end,1) % Out of range to the right
%         dy = (obj.tabW(end,2) - obj.tabW(end-1,2))/(obj.tabW(end,1) - obj.tabW(end-1,1));
%         y = obj.tabW(end-1,2) + dy*(x - obj.tabW(end-1,1));
%         if y < 0
%           y = 0;
%           dy = 0;
%         end
%       end
%     end
  end

  methods (Access = private)
    function readTabularCurve(obj,fID, matFileName)
      curveFname = readToken(fID, matFileName);
      obj.tabW = load(curveFname);
      obj.nPoints = size(obj.tabW,1);
      obj.derivW = diff(obj.tabW(:,2))./diff(obj.tabW(:,1));
      obj.derivW2 = ((obj.tabW(3:end,2)-obj.tabW(2:end-1,2))./(obj.tabW(3:end,1)-obj.tabW(2:end-1,1)) - ...
          (obj.tabW(2:end-1,2)-obj.tabW(1:end-2,2))./(obj.tabW(2:end-1,1)-obj.tabW(1:end-2,1)))./(0.5*...
          (obj.tabW(3:end,1)-obj.tabW(1:end-2,1)));
    end
  end
end