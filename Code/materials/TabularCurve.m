classdef TabularCurve < handle
  % ELASTIC ISOTROPIC material class

  properties (Access = private)
    tabW
    derivW
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
      dy = (obj.tabW(1,2) - 1)/obj.tabW(1,1);
      varargout{1}(xOutL) = 1 + dy*x(xOutL);
      %
      varargout{1}(xIn) = interp1(obj.tabW(:,1),obj.tabW(:,2),x(xIn));
      %
      dy = (obj.tabW(end,2) - obj.tabW(end-1,2))/(obj.tabW(end,1) - obj.tabW(end-1,1));
      varargout{1}(xOutR) = obj.tabW(end-1,2) + dy*(x(xOutR) - obj.tabW(end-1,1));
      idOutR = find(xOutR);
      yNeg = varargout{1}(xOutR) < 0;
      varargout{1}(idOutR(yNeg)) = 0;
      if nargout == 2
        varargout{2} = zeros(length(x),1);
        %
        varargout{2}(xNeg) = 0;
        %
        varargout{2}(xOutL) = dy;
        %
        [~,~,binID] = histcounts(x(xIn),obj.tabW(:,1));
        varargout{2}(xIn) = obj.derivW(binID);
        %
        varargout{2}(idOutR) = dy;
        varargout{2}(idOutR(yNeg)) = 0;
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
    end
  end
end